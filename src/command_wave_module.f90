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
      use hydrograph_module, only : ob, sp_ob, sp_ob1
      implicit none
      integer :: ic, in, iob, newlev, lev, k, maxcnt, npass
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
          if (ob(ic)%rcv_tot == 0) then
            newlev = 1
          else
            newlev = 1
            do in = 1, ob(ic)%rcv_tot
              iob = ob(ic)%obj_in(in)
              if (iob >= 1 .and. iob <= sp_ob%objs) then
                if (ob(iob)%cmd_order + 1 > newlev) newlev = ob(iob)%cmd_order + 1
              end if
            end do
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

      end module command_wave_module
