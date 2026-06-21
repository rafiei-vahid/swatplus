      subroutine pfas_xval_selftest

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Self-contained numerical cross-validation of the dedicated in-stream PFAS
!!    routing (pfas_cha) against the trusted SWAT+ pesticide channel routine
!!    (ch_rtpest).  Sets up ONE synthetic reach with identical geometry, the
!!    same per-day inflow load, and matched routing parameters, then drives a
!!    sweep of channel-days through BOTH routines and writes the per-day
!!    comparison of all 9 shared mass terms to pfas_xval.csv.
!!
!!    A matched pesticide (decay_a=1, decay_b=1, aq_volat=0, num_metab=0) reduces
!!    ch_rtpest EXACTLY to pfas_cha's linear-Koc settle/resus/diffuse/bury
!!    physics, so the two must agree to floating-point tolerance.
!!
!!    Runs BEFORE any model init (called from main, gated by sentinel file
!!    "pfas_xval.run") and STOPs -- so no model state is required or touched.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      use hydrograph_module, only : jrch, ht1, ht2, ch_stor
      use channel_module,    only : rttime
      use sd_channel_module, only : sd_ch, rcurv
      use constituent_mass_module
      use pesticide_data_module
      use ch_pesticide_module, only : chpst, chpstz, chpst_d,             &
     &                                frsol_p => frsol, frsrb_p => frsrb
      use pfas_module
      use pfas_cha_module

      implicit none

      external :: ch_rtpest, pfas_cha

      integer :: id = 0
      real    :: amix = 0.            !! shared aquatic mixing velocity
      real    :: load = 0.            !! daily PFAS/pesticide load into reach (kg)
      real    :: pv(9) = 0., fv(9) = 0., dmax = 0., gmax = 0., grel = 0.
      real    :: mw_g = 500.13        !! g/mol (PFOS); pfasdb%mw = mw_g/1000

      !! ---- matched routing parameters (pesticide units == pfas units) ----
      real, parameter :: p_koc      = 100.0
      real, parameter :: p_solub    = 680.0
      real, parameter :: p_settle   = 1.0
      real, parameter :: p_resus    = 0.05
      real, parameter :: p_bury     = 0.001
      real, parameter :: p_actdep   = 0.1
      real, parameter :: p_carbon   = 2.0      !! % bed carbon
      real, parameter :: p_chbd     = 1.25     !! t/m3

!!    ~ ~ ~ allocate + populate the minimal globals both routines read ~ ~ ~

      jrch = 1

      !! pesticide database (1 compound), neutralised to PFAS-equivalent physics
      cs_db%num_pests = 1
      allocate (cs_db%pest_num(0:1), source = 0)
      cs_db%pest_num(1) = 1
      allocate (pestdb(1))
      pestdb(1)%koc        = p_koc
      pestdb(1)%solub      = p_solub
      pestdb(1)%aq_settle  = p_settle
      pestdb(1)%aq_resus   = p_resus
      pestdb(1)%ben_bury   = p_bury
      pestdb(1)%ben_act_dep= p_actdep
      pestdb(1)%aq_volat   = 0.
      pestdb(1)%mol_wt     = mw_g
      allocate (pestcp(1))
      pestcp(1)%decay_a    = 1.0        !! no aqueous decay
      pestcp(1)%decay_b    = 1.0        !! no benthic decay
      pestcp(1)%num_metab  = 0          !! no daughters

      !! PFAS database (1 compound) matched to the pesticide
      npfas = 1
      allocate (pfas_num(1));  pfas_num(1) = 1
      allocate (pfasdb(1))
      pfasdb(1)%mw  = mw_g / 1000.      !! kg/mol
      pfasdb(1)%sol = p_solub
      allocate (pfas_chadb(1))
      pfas_chadb(1)%koc        = p_koc
      pfas_chadb(1)%aq_settle  = p_settle
      pfas_chadb(1)%aq_resus   = p_resus
      pfas_chadb(1)%ben_bury   = p_bury
      pfas_chadb(1)%ben_act_dep= p_actdep

      !! channel geometry + shared aquatic mixing velocity (identical for both)
      amix = mw_g ** (-.6666) * (1. - p_chbd / 2.65) * (69.35 / 365)
      allocate (sd_ch(0:1))
      sd_ch(1)%carbon = p_carbon
      sd_ch(1)%ch_bd  = p_chbd
      sd_ch(1)%chw    = 5.0
      sd_ch(1)%chl    = 2.0
      allocate (sd_ch(1)%aq_mix(1));       sd_ch(1)%aq_mix(1)      = amix
      allocate (sd_ch(1)%aq_mix_pfas(1));  sd_ch(1)%aq_mix_pfas(1) = amix

      !! channel storage + flow scratch
      allocate (ch_stor(0:1))
      ch_stor(1)%flo = 0.

      !! channel water/benthic pools (start clean, evolve in step)
      allocate (ch_water(0:1), ch_benthic(0:1))
      allocate (ch_pfas_water(0:1), ch_pfas_benthic(0:1))
      allocate (ch_water(1)%pest(1),   source = 0.)
      allocate (ch_benthic(1)%pest(1), source = 0.)
      allocate (ch_pfas_water(1)%pfas(1),   source = 0.)
      allocate (ch_pfas_benthic(1)%pfas(1), source = 0.)

      !! inflow/outflow constituent hydrographs
      allocate (hcs1%pest(1), hcs2%pest(1), source = 0.)
      allocate (hcs1%pfas(1), hcs2%pfas(1), source = 0.)

      !! daily reach output accumulators
      allocate (chpstz%pest(1), chpst%pest(1))
      allocate (chpst_d(0:1))
      allocate (chpst_d(1)%pest(1))
      allocate (chpfasz%pfas(1))
      allocate (chpfas_d(0:1))
      allocate (chpfas_d(1)%pfas(1))

      open (9990, file='pfas_xval.csv')
      write (9990,'(a)') 'day,load_kg,tot_in_pest,tot_in_pfas,'//          &
     &  'settle_p,settle_f,resus_p,resus_f,difus_p,difus_f,'//             &
     &  'bury_p,bury_f,water_p,water_f,benthic_p,benthic_f,maxabsdiff'

!!    ~ ~ ~ drive a sweep of synthetic channel-days through both routines ~ ~ ~
      do id = 1, 60
        !! varying hydrology (a dry day at id=30 to exercise no-flow branches)
        if (id == 30) then
          ht1%flo = 0.;  ht1%sed = 0.;  ht2%flo = 0.;  ch_stor(1)%flo = 0.
        else
          ht1%flo = 1.0e4 * real(1 + mod(id, 12))        !! m3 inflow
          ht1%sed = 5.0  * real(mod(id, 8))              !! tons sediment
          ht2%flo = 0.75 * ht1%flo                       !! m3 outflow
          ch_stor(1)%flo = 3.0e3 * real(1 + mod(id, 4))  !! m3 storage
        end if
        rcurv%dep = 0.4 + 0.15 * real(mod(id, 6))        !! m depth
        rttime    = 6.0 + real(mod(id, 18))              !! hr travel time

        !! identical daily inflow load into BOTH routines
        load = 0.02 * real(1 + mod(id, 7))
        hcs1%pest(1) = load
        hcs1%pfas(1) = load

        !! zero the pesticide working accumulator each day (ch_rtpest only zeros
        !! chpst_d, not chpst, so on a fully-dry reach it leaks stale process
        !! terms; pfas_cha zeros chpfas_d every day -> mirror that for fairness)
        chpst%pest(1) = chpstz%pest(1)
        call ch_rtpest
        !! replicate the caller-side (sd_channel_control3) pesticide output
        !! assembly: ch_rtpest leaves process terms in chpst + outflow in
        !! hcs2%pest / ch_water / ch_benthic; the caller builds the rest.
        pv = (/ load,                                                      &
     &          chpst%pest(1)%settle, chpst%pest(1)%resus,                 &
     &          chpst%pest(1)%difus,  chpst%pest(1)%bury,                  &
     &          ch_water(1)%pest(1),  ch_benthic(1)%pest(1),               &
     &          frsol_p * hcs2%pest(1), frsrb_p * hcs2%pest(1) /)

        call pfas_cha
        !! pfas_cha assembles chpfas_d fully on its own
        fv = (/ chpfas_d(1)%pfas(1)%tot_in,  chpfas_d(1)%pfas(1)%settle,   &
     &          chpfas_d(1)%pfas(1)%resus,   chpfas_d(1)%pfas(1)%difus,    &
     &          chpfas_d(1)%pfas(1)%bury,    chpfas_d(1)%pfas(1)%water,    &
     &          chpfas_d(1)%pfas(1)%benthic, chpfas_d(1)%pfas(1)%sol_out,  &
     &          chpfas_d(1)%pfas(1)%sor_out /)
        dmax = maxval (abs (pv - fv))
        if (dmax > gmax) gmax = dmax
        if (maxval(abs(pv)) > 1.e-9) grel = max (grel, dmax / maxval(abs(pv)))

        !! columns: pest vs pfas for tot_in,settle,resus,difus,bury,water,benthic
        write (9990,'(i4,16(",",es15.7),",",es15.7)') id, load,           &
     &    pv(1), fv(1), pv(2), fv(2), pv(3), fv(3), pv(4), fv(4),          &
     &    pv(5), fv(5), pv(6), fv(6), pv(7), fv(7), dmax
      end do
      close (9990)

      write (*,*) '======================================================'
      write (*,*) ' PFAS in-stream cross-validation vs ch_rtpest (60 days)'
      write (*,*) '   max abs diff (kg) : ', gmax
      write (*,*) '   max rel diff      : ', grel
      if (gmax < 1.e-6) then
        write (*,*) '   VERDICT           :  PASS (pfas_cha == ch_rtpest)'
      else
        write (*,*) '   VERDICT           :  FAIL'
      end if
      write (*,*) '======================================================'

      stop
      end subroutine pfas_xval_selftest
