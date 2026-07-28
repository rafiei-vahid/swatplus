      module command_wave_module
!!    swatplus_perf OpenMP wavefront support.
!!    The SWAT+ daily step walks objects in a single topological order (cmd_next).
!!    Objects that share a topological "level" (longest path from a headwater leaf)
!!    are mutually independent: none is upstream of another, so they can run
!!    concurrently. This module computes that level for every object and groups the
!!    HRU objects by level, so the land phase can be driven wave-by-wave under an
!!    !$omp parallel do (HRU->HRU landscape routing is respected because a receiving
!!    HRU lands on a strictly higher level than every HRU that feeds it).
!!    Building this is read-only w.r.t. simulation state, so it does not change output.
      implicit none

      integer :: hru_nwave = 0                                   !number of HRU levels (max HRU cmd_order)
      integer, dimension(:),   allocatable :: hru_wave_cnt       !(level) number of HRU objects at that level
      integer, dimension(:,:), allocatable :: hru_wave_obj       !(level,k) -> object (icmd) index of k-th HRU
      logical :: hru_wave_ready = .false.

      !! ISOLATION EXPERIMENT support (SWATPLUS_HRU_SERIAL=1): per-level flag, true when
      !! the level contains at least one HRU object. Used by command.f90 to run those
      !! levels on the master thread while non-HRU levels stay parallel, isolating
      !! "HRU parallel" from "routing parallel". All .false. unless the env var is set,
      !! so the shipped behaviour is untouched.
      logical, dimension(:), allocatable :: obj_wave_has_hru

      !! Phase C full-DAG wave: ALL command objects (hru/ru/channel/res/aqu/...) bucketed
      !! by cmd_order level. Same-level objects are mutually independent -> run concurrently.
      integer :: obj_nwave = 0                                   !number of object levels (max cmd_order)
      integer, dimension(:),   allocatable :: obj_wave_cnt       !(level) number of command objects at that level
      integer, dimension(:,:), allocatable :: obj_wave_obj       !(level,k) -> command object index of k-th object

      contains

      subroutine command_wave_build
!!    Compute ob(:)%cmd_order = longest path from a headwater leaf (fixpoint over the
!!    DAG), then bucket HRU-type objects by level. Called once after connectivity is
!!    final; safe to call again (idempotent) - it reallocates.
      use hydrograph_module, only : ob, sp_ob, sp_ob1, ru_def, ru_elem
      implicit none
      integer :: ic, in, iob, newlev, lev, k, maxcnt, npass
      integer :: iru_l, ise_l, ie
      integer :: env_st
      character(len=8) :: hru_ser_env
      logical :: changed

      if (sp_ob%objs <= 0) return

      !! longest-path levels: leaves = 1, else 1 + max(level of receiving objects)
      do ic = 1, sp_ob%objs
        ob(ic)%cmd_order = 0
      end do
      changed = .true.
      npass = 0
      do while (changed)
        changed = .false.
        npass = npass + 1
        ic = sp_ob1%objs
        do while (ic /= 0)
          newlev = 1
          do in = 1, ob(ic)%rcv_tot
            iob = ob(ic)%obj_in(in)
            if (iob >= 1 .and. iob <= sp_ob%objs) then
              if (ob(iob)%cmd_order + 1 > newlev) newlev = ob(iob)%cmd_order + 1
            end if
          end do

          !! ------------------------------------------------------------------
          !! A routing unit consumes its element HRUs' daily hydrographs
          !! (ru_control: iob = ru_elem(ru_def(iru)%num(ielem))%obj, then
          !! ht1 = ob(iob)%hd(...)), but that coupling lives in ru_def/ru_elem --
          !! a MEMBERSHIP table, not a .con record -- so it never appears in
          !! obj_in and rcv_tot is 0. Ordering by flow edges alone therefore put
          !! every element-only ru on level 1 beside the very HRUs it reads.
          !!
          !! Stock SWAT+ guards this in hyd_connect (".. subbasin has to be in
          !! parallel order after elements in the subbasin"); rewriting the level
          !! assignment as a longest-path fixpoint dropped the guard, and the
          !! wavefront then read one-day-stale hydrographs -- ~95% median daily
          !! error on the event-driven constituents (orgn/sedp/nh3) while flow and
          !! no3, being baseflow-dominated and autocorrelated, moved only ~2%.
          !!
          !! Take the max over the ELEMENT objects too. This is stronger than the
          !! stock "iorder = 1" floor, which only forces level >= 2 and would break
          !! if an element were itself above level 1.
          if (ob(ic)%typ == "ru") then
            iru_l = ob(ic)%num
            if (iru_l >= 1 .and. allocated(ru_def)) then
              if (iru_l <= size(ru_def)) then
                if (allocated(ru_def(iru_l)%num)) then
                  do ie = 1, ru_def(iru_l)%num_tot
                    ise_l = ru_def(iru_l)%num(ie)
                    if (ise_l >= 1 .and. allocated(ru_elem)) then
                      if (ise_l <= size(ru_elem)) then
                        iob = ru_elem(ise_l)%obj
                        if (iob >= 1 .and. iob <= sp_ob%objs) then
                          if (ob(iob)%cmd_order + 1 > newlev)                     &
                            newlev = ob(iob)%cmd_order + 1
                        end if
                      end if
                    end if
                  end do
                end if
              end if
            end if
          end if
          if (newlev /= ob(ic)%cmd_order) then
            ob(ic)%cmd_order = newlev
            changed = .true.
          end if
          ic = ob(ic)%cmd_next
        end do
        if (npass > sp_ob%objs + 2) exit   !safety: cannot exceed object count for a DAG
      end do

      !! HRU wave buckets
      hru_nwave = 0
      do ic = 1, sp_ob%objs
        if (ob(ic)%typ == "hru" .and. ob(ic)%cmd_order > hru_nwave) hru_nwave = ob(ic)%cmd_order
      end do
      if (hru_nwave <= 0) then
        hru_wave_ready = .true.
        return
      end if

      if (allocated(hru_wave_cnt)) deallocate (hru_wave_cnt)
      allocate (hru_wave_cnt(hru_nwave))
      hru_wave_cnt = 0
      do ic = 1, sp_ob%objs
        if (ob(ic)%typ == "hru") then
          lev = ob(ic)%cmd_order
          hru_wave_cnt(lev) = hru_wave_cnt(lev) + 1
        end if
      end do

      maxcnt = 0
      do lev = 1, hru_nwave
        if (hru_wave_cnt(lev) > maxcnt) maxcnt = hru_wave_cnt(lev)
      end do

      if (allocated(hru_wave_obj)) deallocate (hru_wave_obj)
      allocate (hru_wave_obj(hru_nwave, maxcnt))
      hru_wave_obj = 0
      hru_wave_cnt = 0
      do ic = 1, sp_ob%objs
        if (ob(ic)%typ == "hru") then
          lev = ob(ic)%cmd_order
          k = hru_wave_cnt(lev) + 1
          hru_wave_cnt(lev) = k
          hru_wave_obj(lev, k) = ic
        end if
      end do

      !! Phase C: bucket ALL command objects by level (walk cmd_next = exact command set).
      obj_nwave = 0
      ic = sp_ob1%objs
      do while (ic /= 0)
        if (ob(ic)%cmd_order > obj_nwave) obj_nwave = ob(ic)%cmd_order
        ic = ob(ic)%cmd_next
      end do
      if (obj_nwave > 0) then
        if (allocated(obj_wave_cnt)) deallocate (obj_wave_cnt)
        allocate (obj_wave_cnt(obj_nwave)); obj_wave_cnt = 0
        ic = sp_ob1%objs
        do while (ic /= 0)
          lev = ob(ic)%cmd_order
          if (lev >= 1) obj_wave_cnt(lev) = obj_wave_cnt(lev) + 1
          ic = ob(ic)%cmd_next
        end do
        maxcnt = 0
        do lev = 1, obj_nwave
          if (obj_wave_cnt(lev) > maxcnt) maxcnt = obj_wave_cnt(lev)
        end do
        if (allocated(obj_wave_obj)) deallocate (obj_wave_obj)
        allocate (obj_wave_obj(obj_nwave, maxcnt)); obj_wave_obj = 0

        !! ISOLATION EXPERIMENT: mark levels containing at least one HRU. Populated ONLY
        !! when SWATPLUS_HRU_SERIAL=1, so the flag stays all-.false. and the shipped
        !! scheduling is bit-for-bit unchanged in normal operation.
        if (allocated(obj_wave_has_hru)) deallocate (obj_wave_has_hru)
        allocate (obj_wave_has_hru(obj_nwave)); obj_wave_has_hru = .false.
        hru_ser_env = " "
        call get_environment_variable ("SWATPLUS_HRU_SERIAL", hru_ser_env, status=env_st)
        if (env_st == 0 .and. (hru_ser_env(1:1) == "1" .or. hru_ser_env(1:1) == "y")) then
          ic = sp_ob1%objs
          do while (ic /= 0)
            lev = ob(ic)%cmd_order
            if (lev >= 1 .and. ob(ic)%typ == "hru") obj_wave_has_hru(lev) = .true.
            ic = ob(ic)%cmd_next
          end do
        end if
        obj_wave_cnt = 0
        ic = sp_ob1%objs
        do while (ic /= 0)
          lev = ob(ic)%cmd_order
          if (lev >= 1) then
            k = obj_wave_cnt(lev) + 1
            obj_wave_cnt(lev) = k
            obj_wave_obj(lev, k) = ic
          end if
          ic = ob(ic)%cmd_next
        end do
      end if

      hru_wave_ready = .true.

      !! swatplus_perf diagnostic: dump the wave histogram so we can see HRU-phase
      !! parallelism width (objects per level). Separate file - not a simulation output.
      open (9123, file = "openmp_waves.out", status = "replace")
      write (9123, '(a)') "swatplus_perf HRU wavefront (cmd_order levels)"
      write (9123, '(a,i0)') "total objects: ", sp_ob%objs
      write (9123, '(a,i0)') "total HRUs:    ", sp_ob%hru
      write (9123, '(a,i0)') "HRU waves:     ", hru_nwave
      write (9123, '(a)') "level   n_hru"
      do lev = 1, hru_nwave
        write (9123, '(i5,3x,i6)') lev, hru_wave_cnt(lev)
      end do
      write (9123, '(a,i0)') "full-DAG object waves: ", obj_nwave
      write (9123, '(a)') "level   n_obj"
      do lev = 1, obj_nwave
        write (9123, '(i5,3x,i6)') lev, obj_wave_cnt(lev)
      end do
      close (9123)

      return
      end subroutine command_wave_build

      logical function lev_is_hru (lev)
      implicit none
      integer, intent(in) :: lev
      lev_is_hru = .false.
      if (.not. allocated(obj_wave_has_hru)) return
      if (lev < 1 .or. lev > size(obj_wave_has_hru)) return
      lev_is_hru = obj_wave_has_hru(lev)
      return
      end function lev_is_hru

      end module command_wave_module
