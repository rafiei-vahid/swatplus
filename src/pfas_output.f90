      subroutine pfas_output

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    End-of-run PFAS land-phase output + mass-balance verification for the
!!    surface-water-only PFAS implementation.  Writes:
!!      * pfas_hru_aa.txt  : per-HRU, per-PFAS run-cumulative land losses (kg/ha)
!!                           surface runoff / lateral / leach / sediment + initial
!!                           and final soil pool, with the per-HRU closure residual.
!!      * pfas_balance.out : basin-total (area-weighted, kg) initial vs final soil
!!                           mass and losses by pathway, with the closure residual
!!                           and percent — the land-phase mass-balance proof.
!!
!!    Closure check: the soil pool changes ONLY via the four loss pathways (no
!!    land-phase source in this configuration), so for every HRU and PFAS
!!        init_pool - final_pool  ==  surq + latq + perc + sed   (to FP round-off).
!!    A non-zero residual would mean the decrement/accumulation bookkeeping is
!!    inconsistent.  Called once from main after time_control.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      use pfas_module
      use pfas_output_module
      use hydrograph_module, only : sp_ob
      use soil_module, only : soil
      use hru_module, only : hru
      use pfas_cha_module, only : pfdiag_in, pfdiag_out, pfdiag_bury,        &
     &                            pfdiag_active
      use constituent_mass_module, only : ch_pfas_water, ch_pfas_benthic

      implicit none

      integer :: j, ly, k, nly, nhru, jworst
      real :: finalm, lossm, resid, rel, maxres, maxrel, area
      real(8), dimension(:), allocatable :: tinit, tfinal, tsurq, tlatq, tperc, tsed

      if (npfas <= 0) return

      !! in-stream PFAS routing mass-balance summary (no-op if pfas_cha never ran)
      block
        integer :: ic, ipf
        real(8) :: stored
        stored = 0.d0
        if (allocated(ch_pfas_water)) then
          do ic = 1, ubound(ch_pfas_water, 1)
            do ipf = 1, npfas
              stored = stored + ch_pfas_water(ic)%pfas(ipf)                  &
     &                        + ch_pfas_benthic(ic)%pfas(ipf)
            end do
          end do
          open (7710, file="pfas_cha_balance.out")
          write (7710,'(a)') "PFAS in-stream (channel) routing mass balance "  &
     &      // "(run-cumulative, kg)"
          write (7710,'(a,i12)')   " reach-days with PFAS inflow : ", pfdiag_active
          write (7710,'(a,es15.7)')" cumulative reach inflow     : ", pfdiag_in
          write (7710,'(a,es15.7)')" cumulative reach outflow    : ", pfdiag_out
          write (7710,'(a,es15.7)')" cumulative burial (sink)    : ", pfdiag_bury
          write (7710,'(a,es15.7)')" final channel storage       : ", real(stored)
          write (7710,'(a)') " note: in - out == bury + final storage (no "    &
     &      // "decay/volat); reach in/out are pass-through sums across reaches"
          close (7710)
        end if
      end block

      if (.not. allocated(pfas_init_hru)) return
      nhru = sp_ob%hru
      if (nhru <= 0) return

      allocate (tinit(npfas), tfinal(npfas), tsurq(npfas), tlatq(npfas),    &
     &          tperc(npfas), tsed(npfas))
      tinit = 0.d0; tfinal = 0.d0; tsurq = 0.d0; tlatq = 0.d0
      tperc = 0.d0; tsed = 0.d0
      maxres = 0.; maxrel = 0.; jworst = 0

      open (7701, file="pfas_hru_aa.txt")
      write (7701,'(a)') "hru   pfas   init_kgha   final_kgha   surq_kgha   "  &
     &  // "latq_kgha   perc_kgha   sed_kgha   resid_kgha"

      do j = 1, nhru
        if (.not. allocated(pfas_soil_hru(j)%ly)) cycle
        nly = soil(j)%nly
        if (nly < 1) nly = 1
        area = hru(j)%area_ha
        do k = 1, npfas
          finalm = 0.
          do ly = 1, nly
            finalm = finalm + pfas_soil_hru(j)%ly(ly)%sol_pfas(k)
          end do
          lossm = hpfasb_d(j)%surq_a(k) + hpfasb_d(j)%latq_a(k)            &
     &          + hpfasb_d(j)%perc_a(k) + hpfasb_d(j)%sed_a(k)
          resid = (pfas_init_hru(j,k) - finalm) - lossm

          tinit(k)  = tinit(k)  + real(pfas_init_hru(j,k),8) * area
          tfinal(k) = tfinal(k) + real(finalm,8) * area
          tsurq(k)  = tsurq(k)  + real(hpfasb_d(j)%surq_a(k),8) * area
          tlatq(k)  = tlatq(k)  + real(hpfasb_d(j)%latq_a(k),8) * area
          tperc(k)  = tperc(k)  + real(hpfasb_d(j)%perc_a(k),8) * area
          tsed(k)   = tsed(k)   + real(hpfasb_d(j)%sed_a(k),8) * area

          if (abs(resid) > maxres) then
            maxres = abs(resid); jworst = j
          end if
          if (pfas_init_hru(j,k) > 0.) then
            rel = abs(resid) / pfas_init_hru(j,k)
            if (rel > maxrel) maxrel = rel
            write (7701,'(i8,i5,7es13.5)') j, k, pfas_init_hru(j,k),       &
     &        finalm, hpfasb_d(j)%surq_a(k), hpfasb_d(j)%latq_a(k),        &
     &        hpfasb_d(j)%perc_a(k), hpfasb_d(j)%sed_a(k), resid
          end if
        end do
      end do
      close (7701)

      open (7702, file="pfas_balance.out")
      write (7702,'(a)') "PFAS surface-water land-phase mass balance "      &
     &  // "(area-weighted, kg over the simulation)"
      write (7702,'(a)') "pfas  name              init_kg       "          &
     &  // "final_kg      surq_kg       latq_kg       perc_kg       "       &
     &  // "sed_kg        loss_kg       resid_kg     resid_pct"
      do k = 1, npfas
        write (7702,'(i4,2x,a16,8es14.5,f10.4)') k, pfasdb(k)%name,        &
     &    tinit(k), tfinal(k), tsurq(k), tlatq(k), tperc(k), tsed(k),      &
     &    tsurq(k)+tlatq(k)+tperc(k)+tsed(k),                              &
     &    tinit(k)-tfinal(k)-(tsurq(k)+tlatq(k)+tperc(k)+tsed(k)),         &
     &    real(100.d0*(tinit(k)-tfinal(k)-(tsurq(k)+tlatq(k)+tperc(k)      &
     &      +tsed(k))) / max(tinit(k),1.d-30))
      end do
      write (7702,'(/,a,es13.5)') "max per-HRU abs closure residual (kg/ha): ", maxres
      write (7702,'(a,es13.5)')   "max per-HRU rel closure residual (-):     ", maxrel
      write (7702,'(a,i9)')       "worst HRU index:                          ", jworst
      close (7702)

      deallocate (tinit, tfinal, tsurq, tlatq, tperc, tsed)

      return
      end subroutine pfas_output
