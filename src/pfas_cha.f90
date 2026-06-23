      subroutine pfas_cha

!!     ~ ~ ~ PURPOSE ~ ~ ~
!!     In-stream + reservoir PFAS fate-and-transport for modern (free-form)
!!     SWAT+ -- SERIAL channel phase.  Computes the daily reach PFAS balance
!!     (soluble + sorbed) using a LINEAR-Koc partition and benthic exchange.
!!
!!     Ported from the SWAT2012 PFAS implementation (Vahid Rafiei):
!!     rtpfas.f (channels) / lakeqpfas.f (reservoirs).  Those legacy routines
!!     use a LINEAR Koc downstream (no Freundlich / Langmuir in-stream) -- i.e.
!!     they are physically identical to the SWAT+ pesticide in-stream model
!!     ch_rtpest.f90: settling, resuspension, sediment-water diffusion, and
!!     burial driven by frsol = 1/(1+kd*sedcon).  PFAS therefore drops the
!!     three pesticide-only loss terms (no chemical/biological reaction, no
!!     volatilization, no metabolite daughters) and keeps everything else.
!!
!!     ~ ~ ~ REUSE NOTE (read before integrating) ~ ~ ~
!!     The CLEANEST integration is NOT to ship this file at all but to route
!!     PFAS THROUGH ch_rtpest as a pesticide-type constituent, because the
!!     SWAT+ pesticide path already carries the full benthic machinery, the
!!     per-reach storage (ch_water / ch_benthic), the inflow hydrograph
!!     (hcs1/hcs2%pest), and the daily/monthly/yearly output (chpst_*).
!!     The exact constituent -> pesticide-database mapping that makes a PFAS
!!     compound behave like rtpfas.f is:
!!
!!       legacy rtpfas variable        SWAT+ pesticide-db / channel field
!!       --------------------------    ----------------------------------
!!       chpfas_koc(jrch) (m3/g)       kd = pestdb(jpst)%koc * carbon/100
!!       chpfas_stl(jrch) (m/day)      pestdb(jpst)%aq_settle
!!       chpfas_rsp(jrch) (m/day)      pestdb(jpst)%aq_resus
!!       chpfas_mix(jrch) (m/day)      sd_ch(jrch)%aq_mix(ipest)
!!       sedpfas_bry(jrch) (m/day)     pestdb(jpst)%ben_bury
!!       sedpfas_act(jrch) (m)         pestdb(jpst)%ben_act_dep
!!       --- ZERO THESE for PFAS (no analogue in rtpfas.f) ---
!!       (none)                        pestcp(jpst)%decay_a   = 1.
!!       (none)                        pestcp(jpst)%decay_b   = 1.
!!       (none)                        pestdb(jpst)%aq_volat  = 0.
!!       (none)                        pestcp(jpst)%num_metab = 0
!!       varoute(34,:)/varoute(35,:)   hcs1%pest(ipest)  (sol+sorbed combined,
!!                                     repartitioned each day by frsol/frsrb)
!!       chpfas_conc*rchwtr            ch_water(jrch)%pest(ipest)   (kg)
!!       sedpfas_conc*bedvol           ch_benthic(jrch)%pest(ipest) (kg)
!!       solpfaso / sorpfaso           frsol*hcs2%pest / frsrb*hcs2%pest
!!
!!     With those settings ch_rtpest reproduces rtpfas.f aside from the koc
!!     carbon-scaling and the (PFAS-absent) reaction/volat terms, and
!!     lakeqpfas.f maps the same way onto res_pest.f90.  If a SEPARATE PFAS
!!     reach output is wanted (distinct columns, no carbon-scaling on koc,
!!     PFAS named from pfasdb rather than pestdb), use THIS routine: it is
!!     the linear-Koc subset of ch_rtpest wired to the pfas_module /
!!     pfas_cha_module containers, run per PFAS in the SERIAL channel phase.
!!
!!     ~ ~ ~ REENTRANCY ~ ~ ~
!!     This is the SERIAL channel phase.  It writes shared per-reach state
!!     (ch_pfas_water/benthic, chpfas_d) and is NOT called from the parallel
!!     land phase, so no OpenMP guard is needed here.  The parallel HRU/soil
!!     PFAS equilibrium lives in pfas_module (pfas_partition).
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use channel_data_module
      use channel_module
      use sd_channel_module
      use pfas_cha_module
      use pfas_module, only : npfas, pfas_num
      use hydrograph_module, only : jrch, ht1, ht2, ch_stor
      use constituent_mass_module, only : ch_pfas_water, ch_pfas_benthic, hcs1, hcs2

      implicit none

      integer :: ipf = 0        !none          |PFAS counter - sequential
      integer :: jpf = 0        !none          |PFAS counter from data base
      real :: pfin = 0.         !kg            |total PFAS transported into reach during time step
      real :: kd = 0.           !(mg/kg)/(mg/L)|koc * carbon
      real :: depth = 0.        !m             |depth of water in reach
      real :: chpfmass = 0.     !kg            |mass of PFAS in reach water column
      real :: sedpfmass = 0.    !kg            |mass of PFAS in bed sediment
      real :: fd2 = 0.          !none          |benthic sorbed/total partition factor
      real :: solmax = 0.       !kg            |max soluble PFAS at solubility limit
      real :: sedcon = 0.       !g/m^3         |sediment concentration
      real :: tday = 0.         !none          |flow duration (fraction of 24 hr)
      real :: por = 0.          !none          |porosity of bottom sediments
      real :: rto_out = 0.      !none          |ratio of outflow to (outflow + storage)
      !! wtrin is a module global (channel_module) -- shared with ch_rtpest

      !! zero daily outputs for this reach
      chpfas_d(jrch) = chpfasz

      !! initialize depth of water for PFAS calculations
      depth = rcurv%dep
      if (depth < 0.01) then
        depth = .01
      endif

      do ipf = 1, npfas
        jpf = pfas_num(ipf)

        !! volume of water entering reach and stored in reach
        wtrin = ht1%flo + ch_stor(jrch)%flo

        !! PFAS transported into reach during day (kg; sol+sorbed combined)
        pfin = hcs1%pfas(ipf)

        !! calculate mass of PFAS in reach water column
        chpfmass = pfin + ch_pfas_water(jrch)%pfas(ipf)

        !! calculate mass of PFAS in bed sediment
        sedpfmass = ch_pfas_benthic(jrch)%pfas(ipf)

        if (chpfmass + sedpfmass < 1.e-12) then
          ch_pfas_water(jrch)%pfas(ipf) = 0.
          ch_pfas_benthic(jrch)%pfas(ipf) = 0.
        end if
        if (chpfmass + sedpfmass < 1.e-12) cycle

        !!in-stream processes
        if (wtrin / 86400. > 1.e-9) then
          !! calculate sediment concentration (g/m^3)
          sedcon = ht1%sed / wtrin * 1.e6

          !! set kd (linear Koc * organic carbon fraction)
          kd = pfas_chadb(jpf)%koc * sd_ch(jrch)%carbon / 100.

          !! calculate fraction of soluble and sorbed PFAS
          if (kd > 0.) then
            frsol = 1. / (1. + kd * sedcon)
          else
            frsol = 1.
          end if
          frsrb = 1. - frsol

          !! ASSUME DENSITY=2.65E6; KD2=KD1 (benthic partition)
          por = 1. - sd_ch(jrch)%ch_bd / 2.65
          fd2 = 1. / (por + kd)

          !! calculate flow duration
          tday = rttime / 24.0
          if (tday > 1.0) tday = 1.0

          !! -----------------------------------------------------------
          !! NOTE: pesticide reaction (decay_a) and volatilization
          !! (aq_volat) terms of ch_rtpest are intentionally OMITTED --
          !! PFAS are non-volatile and non-degradable in this model.
          !! -----------------------------------------------------------

          !! calculate amount of PFAS removed from reach by settling
          chpfas_d(jrch)%pfas(ipf)%settle = pfas_chadb(jpf)%aq_settle *  &
     &        frsrb * chpfmass * tday / depth
          if (chpfas_d(jrch)%pfas(ipf)%settle > frsrb * chpfmass) then
            chpfas_d(jrch)%pfas(ipf)%settle = frsrb * chpfmass
            chpfmass = chpfmass - chpfas_d(jrch)%pfas(ipf)%settle
          else
            chpfmass = chpfmass - chpfas_d(jrch)%pfas(ipf)%settle
          end if
          sedpfmass = sedpfmass + chpfas_d(jrch)%pfas(ipf)%settle

          !! calculate resuspension of PFAS in reach
          chpfas_d(jrch)%pfas(ipf)%resus = pfas_chadb(jpf)%aq_resus *    &
     &        sedpfmass * tday / depth
          if (chpfas_d(jrch)%pfas(ipf)%resus > sedpfmass) then
            chpfas_d(jrch)%pfas(ipf)%resus = sedpfmass
            sedpfmass = 0.
          else
            sedpfmass = sedpfmass - chpfas_d(jrch)%pfas(ipf)%resus
          end if
          chpfmass = chpfmass + chpfas_d(jrch)%pfas(ipf)%resus

          !! calculate diffusion of PFAS between reach water and sediment
          !! (aq_mix_pfas is PFAS-mol_wt-based + dimensioned by npfas; aq_mix is
          !!  pesticide-dimensioned -- see sd_channel_module / pfas_cha_read)
          chpfas_d(jrch)%pfas(ipf)%difus = sd_ch(jrch)%aq_mix_pfas(ipf) *  &
     &        (fd2 * sedpfmass - frsol * chpfmass) * tday / depth
          if (chpfas_d(jrch)%pfas(ipf)%difus > 0.) then
            if (chpfas_d(jrch)%pfas(ipf)%difus > sedpfmass) then
              chpfas_d(jrch)%pfas(ipf)%difus = sedpfmass
              sedpfmass = 0.
            else
              sedpfmass = sedpfmass - Abs(chpfas_d(jrch)%pfas(ipf)%difus)
            end if
            chpfmass = chpfmass + Abs(chpfas_d(jrch)%pfas(ipf)%difus)
          else
            if (Abs(chpfas_d(jrch)%pfas(ipf)%difus) > chpfmass) then
              chpfas_d(jrch)%pfas(ipf)%difus = -chpfmass
              chpfmass = 0.
            else
              chpfmass = chpfmass - Abs(chpfas_d(jrch)%pfas(ipf)%difus)
            end if
            sedpfmass = sedpfmass + Abs(chpfas_d(jrch)%pfas(ipf)%difus)
          end if

          !! calculate removal of PFAS from active sediment layer by burial
          chpfas_d(jrch)%pfas(ipf)%bury = pfas_chadb(jpf)%ben_bury *     &
     &        sedpfmass / pfas_chadb(jpf)%ben_act_dep
          if (chpfas_d(jrch)%pfas(ipf)%bury > sedpfmass) then
            chpfas_d(jrch)%pfas(ipf)%bury = sedpfmass
            sedpfmass = 0.
          else
            sedpfmass = sedpfmass - chpfas_d(jrch)%pfas(ipf)%bury
          end if

          !! verify that water concentration is at or below solubility
          solmax = pfasdb_sol(jpf) * wtrin
          if (solmax < chpfmass * frsol) then
            sedpfmass = sedpfmass + (chpfmass * frsol - solmax)
            chpfmass = chpfmass - (chpfmass * frsol - solmax)
          end if

        else
          !!insignificant flow -- all PFAS settles to bed
          sedpfmass = sedpfmass + chpfmass
          chpfmass = 0.
        end if

        !! benthic reaction term OMITTED for PFAS (non-degradable)

        !! set new water-column mass (in + store) after processes
        !! (matches ch_rtpest exactly for cross-validation parity: the dead
        !!  write to hcs1 and the no-zero dry branch mirror ch_rtpest L200-204)
        if (wtrin > 1.e-6) then
          hcs1%pfas(ipf) = chpfmass
        else
          sedpfmass = sedpfmass + chpfmass
        end if
        ch_pfas_benthic(jrch)%pfas(ipf) = sedpfmass

        !! calculate outflow and storage in water column
        rto_out = ht2%flo / (1.e-6 + ht2%flo + ch_stor(jrch)%flo)
        rto_out = Min (1., rto_out)
        hcs2%pfas(ipf) = rto_out * chpfmass
        ch_pfas_water(jrch)%pfas(ipf) = (1. - rto_out) * chpfmass

        !! -----------------------------------------------------------
        !! daily reach PFAS output (kg); soluble/sorbed split of outflow
        !! -----------------------------------------------------------
        chpfas_d(jrch)%pfas(ipf)%tot_in  = pfin
        chpfas_d(jrch)%pfas(ipf)%sol_out = frsol * hcs2%pfas(ipf)
        chpfas_d(jrch)%pfas(ipf)%sor_out = frsrb * hcs2%pfas(ipf)
        chpfas_d(jrch)%pfas(ipf)%water   = ch_pfas_water(jrch)%pfas(ipf)
        chpfas_d(jrch)%pfas(ipf)%benthic = ch_pfas_benthic(jrch)%pfas(ipf)

        !! total water-column concentration of reach outflow (ng/L); kg/m3*1e9
        !! = ng/L. The in-stream calibration target vs grab samples. Guarded at
        !! 1 m3/day outflow so essentially-dry headwater reaches (tiny denominator
        !! -> unphysical spikes) report 0 rather than blowing up; flowing reaches,
        !! incl. all gauged/monitored mainstem reaches, are unaffected.
        if (ht2%flo > 1.0) then
          chpfas_d(jrch)%pfas(ipf)%conc = hcs2%pfas(ipf) / ht2%flo * 1.e9
        else
          chpfas_d(jrch)%pfas(ipf)%conc = 0.
        end if

        !! in-stream routing run-cumulative mass balance (-> pfas_cha_balance.out)
        pfdiag_in    = pfdiag_in   + chpfas_d(jrch)%pfas(ipf)%tot_in
        pfdiag_out   = pfdiag_out  + chpfas_d(jrch)%pfas(ipf)%sol_out          &
     &                             + chpfas_d(jrch)%pfas(ipf)%sor_out
        pfdiag_bury  = pfdiag_bury + chpfas_d(jrch)%pfas(ipf)%bury
        if (chpfas_d(jrch)%pfas(ipf)%tot_in > 1.e-12)                          &
     &      pfdiag_active = pfdiag_active + 1

      end do

      return

      contains

!!    solubility lookup (mg/L from pfas_module%pfasdb) without pulling the
!!    whole soil database type into the channel-phase namespace
      real function pfasdb_sol(idb)
        use pfas_module, only : pfasdb
        integer, intent(in) :: idb
        pfasdb_sol = pfasdb(idb)%sol
      end function pfasdb_sol

      end subroutine pfas_cha
