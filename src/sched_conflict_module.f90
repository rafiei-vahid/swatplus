      module sched_conflict_module
!!    ------------------------------------------------------------------------------
!!    STAGE 0 AUDIT — conflict census for augmented-dependency scheduling.
!!
!!    The OpenMP wavefront (command_wave_module) builds levels from HYDROLOGIC FLOW
!!    edges only (ob(ic)%obj_in, from the .con files). That is not sufficient: several
!!    routing objects write shared state indexed by ANOTHER object's id, so two objects
!!    the flow graph believes are independent can collide on the same cell. Measured
!!    consequence on calibrated basins (organic-N median daily error vs serial, full
!!    wavefront): 28% at 137 wetlands, 57% at 1315, 97% at 1838 — the error scales with
!!    lateral-coupling density while flow stays at ~1-2%.
!!
!!    This module does NOT change any behaviour. It counts, before a single simulation
!!    day is run, how much conflict a model actually contains and what enforcing it
!!    would cost in parallelism — the go/no-go for building the real scheduler.
!!
!!    Every coupling below is derivable at BUILD time: the membership tables are filled
!!    during init and never mutated during simulation. Deliberately implemented in
!!    Fortran rather than re-derived in Python from chan-surf.lin / ls_unit.ele, because
!!    reproducing the ru->HRU expansion (sd_channel_surf_link) elsewhere is exactly where
!!    a duplicated-logic bug would hide and silently understate conflict density.
!!
!!    Domains:
!!      FPHRU  channels writing wet(iihru)/wet_in_d(iihru) at a floodplain HRU
!!             (sd_channel_sediment3:158-160), membership sd_ch(ich)%fp%hru
!!      LSU    channels writing lsu_wb_d(ilsu) (ch_temp:201/212/216), membership via
!!             ob%obtyp_in == "ru"
!!      CHSTOR recall objects writing ch_stor(ichan) (recall_nut/salt/cs)
!!      WST    objects writing wst(iwst)%weat -- the write-back has since been REMOVED
!!             from sd_channel_control3/res_control, so this domain should now report
!!             zero multi-writer cells. Kept as a regression sentinel.
!!      AQCH   aquifer->channel precedence via sd_ch(ich)%aqu_link (from aqu-cha.lin,
!!             NOT a .con edge, so the flow graph cannot enforce it)
!!      RUELEM routing-unit -> element-HRU precedence: ru_control reads each element's
!!             ob(iob)%hd(...) through ru_def/ru_elem, a MEMBERSHIP table rather than a
!!             .con record, so rcv_tot is 0 and flow-only ordering put every ru on
!!             level 1 beside the HRUs it consumes. This is the domain that actually
!!             broke the full wavefront (~95% median daily error on orgn/sedp/nh3);
!!             the audit was blind to it until 2026-07-28 precisely because it had no
!!             RUELEM domain -- see OVERLAP_DIAGNOSIS_2026-07-28.md.
!!
!!    Enabled by SWATPLUS_SCHED_AUDIT=1; writes conflict_audit.out and stops.
!!    ------------------------------------------------------------------------------
      implicit none

      integer, parameter :: NDOM = 6
      integer :: sched_pred_dropped = 0
      integer, parameter :: DOM_FPHRU = 1, DOM_LSU = 2, DOM_CHSTOR = 3,             &
                            DOM_WST = 4, DOM_AQCH = 5, DOM_RUELEM = 6
      character(len=6), dimension(NDOM), parameter :: DOM_NAME =                    &
        (/ "FPHRU ", "LSU   ", "CHSTOR", "WST   ", "AQCH  ", "RUELEM" /)

      contains

      logical function sched_audit_requested ()
      implicit none
      character(len=8) :: env
      integer :: st
      env = " "
      call get_environment_variable ("SWATPLUS_SCHED_AUDIT", env, status=st)
      sched_audit_requested = (st == 0 .and. (env(1:1) == "1" .or. env(1:1) == "y"))
      return
      end function sched_audit_requested


!!    Serial rank: position in the cmd_next walk. This is THE reference order — any
!!    orientation of a conflict must follow it or bit-identity is unreachable.
      subroutine sched_serial_rank (rank, nranked)
      use hydrograph_module, only : ob, sp_ob, sp_ob1
      implicit none
      integer, dimension(:), intent(out) :: rank
      integer, intent(out) :: nranked
      integer :: ic, r
      rank = 0
      r = 0
      ic = sp_ob1%objs
      do while (ic /= 0)
        r = r + 1
        rank(ic) = r
        ic = ob(ic)%cmd_next
      end do
      nranked = r
      return
      end subroutine sched_serial_rank


!!    Level widths for a given level assignment, plus the realizable time on P threads
!!    (sum over levels of ceil(width/P)) — the parallelism metric that matters.
      subroutine sched_level_stats (lev, nobj, nlev, maxw, nnarrow, t1, t8)
      implicit none
      integer, dimension(:), intent(in) :: lev
      integer, intent(in) :: nobj
      integer, intent(out) :: nlev, maxw, nnarrow, t1, t8
      integer, allocatable :: w(:)
      integer :: i, L
      nlev = 0
      do i = 1, nobj
        if (lev(i) > nlev) nlev = lev(i)
      end do
      if (nlev <= 0) then
        maxw = 0; nnarrow = 0; t1 = 0; t8 = 0
        return
      end if
      allocate (w(nlev))
      w = 0
      do i = 1, nobj
        if (lev(i) >= 1) w(lev(i)) = w(lev(i)) + 1
      end do
      maxw = 0; nnarrow = 0; t1 = 0; t8 = 0
      do L = 1, nlev
        if (w(L) > maxw) maxw = w(L)
        if (w(L) > 0 .and. w(L) <= 16) nnarrow = nnarrow + 1
        t1 = t1 + w(L)
        t8 = t8 + (w(L) + 7) / 8
      end do
      deallocate (w)
      return
      end subroutine sched_level_stats


!!    ------------------------------------------------------------------------------
!!    The audit. Builds the writer lists per (domain, cell), reports multiplicity and
!!    co-level conflict pairs, then recomputes levels with the conflict edges added
!!    (chain of k-1 edges per cell, ordered by serial rank -- the transitive reduction,
!!    O(k) not O(k^2)) and reports what that costs.
!!    ------------------------------------------------------------------------------
      subroutine sched_conflict_report
      use hydrograph_module, only : ob, sp_ob, sp_ob1, ru_def, ru_elem
      use sd_channel_module, only : sd_ch
      implicit none

      integer :: iru_a, ise_a, ie_a
      integer, allocatable :: rank(:), lev_flow(:), lev_ord(:)
      integer, allocatable :: wr_dom(:), wr_cell(:), wr_obj(:)
      integer, allocatable :: khist(:)
      integer :: nobj, nranked, nwr, i, j, ic, in, ihru, iihru, ilsu, iaq, ich
      integer :: kmax_d, npair_d, ncolev_d, ncell_d, nmulti_d, d
      integer :: nlev_f, maxw_f, nnar_f, t1_f, t8_f
      integer :: nlev_o, maxw_o, nnar_o, t1_o, t8_o
      integer :: newlev, iob, npass, ci, cj, k, kk, tmp
      integer, allocatable :: cellw(:)
      logical :: changed
      integer, parameter :: IU = 9124

      nobj = sp_ob%objs
      if (nobj <= 0) return

      allocate (rank(nobj), lev_flow(nobj), lev_ord(nobj))
      call sched_serial_rank (rank, nranked)

      !! ---- collect (domain, cell, writer) triples -------------------------------
      !! Two passes: count, then fill — same idiom as command_wave_module's bucketing.
      nwr = 0
      do i = 1, 2
        if (i == 2) then
          allocate (wr_dom(nwr), wr_cell(nwr), wr_obj(nwr))
          nwr = 0
        end if
        ic = sp_ob1%objs
        do while (ic /= 0)
          !! FPHRU: a channel writes wet()/wet_in_d() at each of its floodplain HRUs
          if (ob(ic)%typ == "chandeg") then
            ich = ob(ic)%props
            if (ich >= 1 .and. allocated(sd_ch)) then
              if (ich <= size(sd_ch)) then
                if (allocated(sd_ch(ich)%fp%hru)) then
                  do ihru = 1, sd_ch(ich)%fp%hru_tot
                    iihru = sd_ch(ich)%fp%hru(ihru)
                    if (iihru >= 1) then
                      nwr = nwr + 1
                      if (i == 2) then
                        wr_dom(nwr) = DOM_FPHRU; wr_cell(nwr) = iihru; wr_obj(nwr) = ic
                      end if
                    end if
                  end do
                end if
                !! AQCH: aquifer -> channel precedence, absent from the flow graph
                iaq = sd_ch(ich)%aqu_link
                if (iaq >= 1) then
                  nwr = nwr + 1
                  if (i == 2) then
                    wr_dom(nwr) = DOM_AQCH; wr_cell(nwr) = iaq; wr_obj(nwr) = ic
                  end if
                end if
              end if
            end if
            !! LSU: a channel writes lsu_wb_d() for each routing unit feeding it
            if (allocated(ob(ic)%obtyp_in)) then
              do in = 1, ob(ic)%rcv_tot
                if (ob(ic)%obtyp_in(in) == "ru") then
                  ilsu = ob(ic)%obtypno_in(in)
                  if (ilsu >= 1) then
                    nwr = nwr + 1
                    if (i == 2) then
                      wr_dom(nwr) = DOM_LSU; wr_cell(nwr) = ilsu; wr_obj(nwr) = ic
                    end if
                  end if
                end if
              end do
            end if
          end if
          !! CHSTOR: a recall writes ch_stor() of its outflow channel
          if (ob(ic)%typ == "recall") then
            if (allocated(ob(ic)%obtypno_out)) then
              if (size(ob(ic)%obtypno_out) >= 1) then
                if (ob(ic)%obtypno_out(1) >= 1) then
                  nwr = nwr + 1
                  if (i == 2) then
                    wr_dom(nwr) = DOM_CHSTOR
                    wr_cell(nwr) = ob(ic)%obtypno_out(1); wr_obj(nwr) = ic
                  end if
                end if
              end if
            end if
          end if
          !! WST: NO LONGER POPULATED. This domain censused "wst(iwst)%weat = w", the
          !! write-back of the shared weather record from sd_channel_control3 and
          !! res_control. Both write-backs have been deleted, so there is no write left to
          !! order and recording the objects that merely SHARE a station overstated the
          !! conflict graph badly -- 45 multi-writer cells and 22,239 pairs on basin
          !! 02297310, which alone accounted for most of the gap between S8_flow and
          !! S8_ord. The domain id is kept so the report keeps a stable column layout and
          !! so the block can be restored verbatim if the write-back ever comes back.
          !! RUELEM: a routing unit reads each element HRU's daily hydrograph. Encoded
          !! like AQCH -- cell = the PRODUCER object (the element), writer = the
          !! CONSUMER (the ru) -- because it is a precedence constraint, not mutual
          !! exclusion. Expanded from the live ru_def/ru_elem tables, never re-derived.
          if (ob(ic)%typ == "ru") then
            iru_a = ob(ic)%num
            if (iru_a >= 1 .and. allocated(ru_def)) then
              if (iru_a <= size(ru_def)) then
                if (allocated(ru_def(iru_a)%num)) then
                  do ie_a = 1, ru_def(iru_a)%num_tot
                    ise_a = ru_def(iru_a)%num(ie_a)
                    if (ise_a >= 1 .and. allocated(ru_elem)) then
                      if (ise_a <= size(ru_elem)) then
                        if (ru_elem(ise_a)%obj >= 1) then
                          nwr = nwr + 1
                          if (i == 2) then
                            wr_dom(nwr) = DOM_RUELEM
                            wr_cell(nwr) = ru_elem(ise_a)%obj
                            wr_obj(nwr) = ic
                          end if
                        end if
                      end if
                    end if
                  end do
                end if
              end if
            end if
          end if
          ic = ob(ic)%cmd_next
        end do
      end do

      !! ---- flow-only levels (recompute so the audit is self-contained) ----------
      call sched_flow_levels (lev_flow, nobj)

      open (IU, file="conflict_audit.out", status="replace")
      write (IU,'(a)') "# SWATGenX scheduling-conflict audit (Stage 0, no simulation)"
      write (IU,'(a,i0)') "objects ", nobj
      write (IU,'(a,i0)') "writer_records ", nwr

      !! ---- per-domain census -----------------------------------------------------
      allocate (cellw(nobj))
      write (IU,'(a)') "# domain n_cells n_multi kmax n_pairs n_colevel_pairs"
      do d = 1, NDOM
        call sched_domain_census (wr_dom, wr_cell, wr_obj, nwr, d, lev_flow,        &
                                  ncell_d, nmulti_d, kmax_d, npair_d, ncolev_d)
        write (IU,'(a,1x,5(i0,1x))') trim(DOM_NAME(d)), ncell_d, nmulti_d,          &
                                     kmax_d, npair_d, ncolev_d
      end do

      !! ---- augmented levels ------------------------------------------------------
      call sched_ord_levels (wr_dom, wr_cell, wr_obj, nwr, rank, lev_ord, nobj)

      call sched_level_stats (lev_flow, nobj, nlev_f, maxw_f, nnar_f, t1_f, t8_f)
      call sched_level_stats (lev_ord,  nobj, nlev_o, maxw_o, nnar_o, t1_o, t8_o)

      write (IU,'(a)') "# schedule nlev maxwidth n_narrow work T8"
      write (IU,'(a,1x,5(i0,1x))') "G_flow", nlev_f, maxw_f, nnar_f, t1_f, t8_f
      write (IU,'(a,1x,5(i0,1x))') "G_ord ", nlev_o, maxw_o, nnar_o, t1_o, t8_o
      if (t8_f > 0 .and. t8_o > 0) then
        write (IU,'(a,f8.4)') "S8_flow ", real(t1_f) / real(t8_f)
        write (IU,'(a,f8.4)') "S8_ord  ", real(t1_o) / real(t8_o)
        write (IU,'(a,f8.4)') "S8_ratio", (real(t1_o)/real(t8_o)) /                 &
                                          (real(t1_f)/real(t8_f))
      end if
      write (IU,'(a,i0)') "pred_edges_dropped ", sched_pred_dropped
      close (IU)

      deallocate (rank, lev_flow, lev_ord, cellw)
      if (allocated(wr_dom)) deallocate (wr_dom, wr_cell, wr_obj)
      return
      end subroutine sched_conflict_report


!!    Longest-path levels from flow edges only — mirrors command_wave_module:38-66.
      subroutine sched_flow_levels (lev, nobj)
      use hydrograph_module, only : ob, sp_ob, sp_ob1
      implicit none
      integer, dimension(:), intent(out) :: lev
      integer, intent(in) :: nobj
      integer :: ic, in, iob, newlev, npass
      logical :: changed
      lev = 0
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
            if (iob >= 1 .and. iob <= nobj) then
              if (lev(iob) + 1 > newlev) newlev = lev(iob) + 1
            end if
          end do
          if (newlev /= lev(ic)) then
            lev(ic) = newlev
            changed = .true.
          end if
          ic = ob(ic)%cmd_next
        end do
        if (npass > nobj + 2) exit
      end do
      return
      end subroutine sched_flow_levels


!!    Levels with conflict edges added. For each (domain,cell) the writers are ordered
!!    by serial rank and chained w(i) -> w(i+1): the transitive reduction, which gives
!!    the identical constraint at O(k) instead of O(k^2) edges. Orientation by serial
!!    rank makes G_flow U C_ord acyclic by construction, since the rank is a linear
!!    extension of the flow graph.
      subroutine sched_ord_levels (wr_dom, wr_cell, wr_obj, nwr, rank, lev, nobj)
      use hydrograph_module, only : ob, sp_ob1
      implicit none
      integer, dimension(:), intent(in) :: wr_dom, wr_cell, wr_obj, rank
      integer, intent(in) :: nwr, nobj
      integer, dimension(:), intent(out) :: lev
      integer, allocatable :: pred(:,:), npred(:)
      !! A routing unit gets one predecessor per element HRU, so this has to comfortably
      !! exceed the largest ru element count -- at 64 the edges were silently dropped and the
      !! reported cost of ordering came out too low. Overflows are counted and reported.
      integer, parameter :: MAXPRED = 1024
      integer :: i, j, ic, in, iob, newlev, npass, p
      logical :: changed

      allocate (pred(nobj, MAXPRED), npred(nobj))
      pred = 0; npred = 0
      call sched_chain_edges (wr_dom, wr_cell, wr_obj, nwr, rank, pred, npred,      &
                              nobj, MAXPRED)

      lev = 0
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
            if (iob >= 1 .and. iob <= nobj) then
              if (lev(iob) + 1 > newlev) newlev = lev(iob) + 1
            end if
          end do
          do p = 1, npred(ic)
            iob = pred(ic, p)
            if (iob >= 1 .and. iob <= nobj) then
              if (lev(iob) + 1 > newlev) newlev = lev(iob) + 1
            end if
          end do
          if (newlev /= lev(ic)) then
            lev(ic) = newlev
            changed = .true.
          end if
          ic = ob(ic)%cmd_next
        end do
        if (npass > nobj + 2) exit
      end do

      deallocate (pred, npred)
      return
      end subroutine sched_ord_levels


!!    Build the chain edges.
!!
!!    Mutual-conflict domains get the transitive reduction: writers of a cell are ordered by
!!    serial rank and chained w(i) -> w(i+1), giving the identical constraint at O(k) instead
!!    of O(k^2) edges. Orientation by serial rank keeps G_flow U C_ord acyclic, since the rank
!!    is a linear extension of the flow graph.
!!
!!    AQCH and RUELEM are PRECEDENCE domains, not mutual-conflict ones: the cell is the
!!    producer (an aquifer, an element HRU) and the writer field is the consumer (a channel, a
!!    routing unit). They get one real edge each, producer -> consumer. RUELEM is the coupling
!!    that command_wave_module now enforces in its own level fixpoint, so emitting it here is
!!    what makes G_ord the schedule the engine actually runs rather than a hypothetical one.
!!    Previously both were skipped, which silently left the reported cost of ordering
!!    understated for AQCH and undefined for RUELEM.
      subroutine sched_chain_edges (wr_dom, wr_cell, wr_obj, nwr, rank, pred, npred, &
                                    nobj, maxpred)
      implicit none
      integer, dimension(:), intent(in) :: wr_dom, wr_cell, wr_obj, rank
      integer, intent(in) :: nwr, nobj, maxpred
      integer, dimension(:,:), intent(inout) :: pred
      integer, dimension(:), intent(inout) :: npred
      integer, allocatable :: grp(:)
      integer :: d, i, j, n, a, b, tmp

      allocate (grp(nwr))

      !! precedence domains: one edge per record, producer(cell) -> consumer(obj)
      do i = 1, nwr
        if (wr_dom(i) /= DOM_AQCH .and. wr_dom(i) /= DOM_RUELEM) cycle
        a = wr_cell(i); b = wr_obj(i)
        if (a < 1 .or. a > nobj .or. b < 1 .or. b > nobj) cycle
        if (npred(b) < maxpred) then
          npred(b) = npred(b) + 1
          pred(b, npred(b)) = a
        else
          sched_pred_dropped = sched_pred_dropped + 1
        end if
      end do

      do d = 1, NDOM
        if (d == DOM_AQCH .or. d == DOM_RUELEM) cycle
        !! gather writers of each cell of this domain, cell by cell
        do i = 1, nwr
          if (wr_dom(i) /= d) cycle
          if (wr_cell(i) < 0) cycle
          n = 0
          do j = i, nwr
            if (wr_dom(j) == d .and. wr_cell(j) == wr_cell(i)) then
              n = n + 1
              grp(n) = wr_obj(j)
            end if
          end do
          if (n < 2) cycle
          !! insertion sort by serial rank (n is small; deterministic)
          do a = 2, n
            tmp = grp(a)
            b = a - 1
            do while (b >= 1)
              if (rank(grp(b)) <= rank(tmp)) exit
              grp(b+1) = grp(b)
              b = b - 1
            end do
            grp(b+1) = tmp
          end do
          do a = 1, n - 1
            b = grp(a+1)
            if (npred(b) < maxpred) then
              npred(b) = npred(b) + 1
              pred(b, npred(b)) = grp(a)
            end if
          end do
        end do
      end do
      deallocate (grp)
      return
      end subroutine sched_chain_edges


!!    Per-domain census: cells, multi-writer cells, k_max, conflict pairs, and the
!!    number of those pairs already CO-LEVEL in the flow schedule (the only pairs that
!!    actually cost anything to order).
      subroutine sched_domain_census (wr_dom, wr_cell, wr_obj, nwr, d, lev_flow,    &
                                      ncell, nmulti, kmax, npair, ncolev)
      implicit none
      integer, dimension(:), intent(in) :: wr_dom, wr_cell, wr_obj, lev_flow
      integer, intent(in) :: nwr, d
      integer, intent(out) :: ncell, nmulti, kmax, npair, ncolev
      integer :: i, j, a, b, n
      logical, allocatable :: seen(:)

      allocate (seen(nwr))
      seen = .false.
      ncell = 0; nmulti = 0; kmax = 0; npair = 0; ncolev = 0
      do i = 1, nwr
        if (wr_dom(i) /= d .or. seen(i)) cycle
        n = 0
        do j = i, nwr
          if (wr_dom(j) == d .and. wr_cell(j) == wr_cell(i)) then
            seen(j) = .true.
            n = n + 1
          end if
        end do
        ncell = ncell + 1
        if (n > kmax) kmax = n
        if (n > 1) then
          nmulti = nmulti + 1
          npair = npair + n * (n - 1) / 2
          !! count co-level pairs
          do a = i, nwr
            if (wr_dom(a) /= d .or. wr_cell(a) /= wr_cell(i)) cycle
            do b = a + 1, nwr
              if (wr_dom(b) /= d .or. wr_cell(b) /= wr_cell(i)) cycle
              if (lev_flow(wr_obj(a)) == lev_flow(wr_obj(b))) ncolev = ncolev + 1
            end do
          end do
        end if
      end do
      deallocate (seen)
      return
      end subroutine sched_domain_census

      end module sched_conflict_module
