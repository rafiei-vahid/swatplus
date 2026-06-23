      subroutine cha_pfas_output(jrch)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Daily per-reach in-stream PFAS output for the surface-water PFAS module.
!!    Writes channel_pfas_day.txt (+ .csv if csvout) with, per channel and PFAS
!!    compound: the daily reach mass balance (kg) and the outflow water
!!    concentration (ng/L) -- the in-stream calibration target vs grab samples.
!!
!!    Files are opened lazily on the first call (units 7110 txt / 7114 csv) so
!!    no header_pest / print.prt change is needed.  Daily only for now; monthly/
!!    yearly aggregation can be added later via chpfas_m/y/a.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      use pfas_cha_module, only : chpfas_d
      use pfas_module, only : npfas, pfasdb, pfas_num
      use time_module
      use basin_module
      use hydrograph_module, only : sp_ob1, ob

      implicit none

      integer, intent (in) :: jrch
      integer :: ipf = 0
      integer :: jpf = 0
      integer :: iob = 0
      logical, save :: opened = .false.

      if (npfas <= 0) return

      !! lazy-open the output files + headers (once)
      if (.not. opened) then
        open (7110, file="channel_pfas_day.txt", recl = 800)
        write (7110,1001) "jday","mon","day","yr","unit","gis_id","name",   &
     &    "pfas","tot_in_kg","sol_out_kg","sor_out_kg","settle_kg",         &
     &    "resus_kg","diffuse_kg","bury_kg","water_kg","benthic_kg",        &
     &    "conc_ngL"
        if (pco%csvout == "y") then
          open (7114, file="channel_pfas_day.csv", recl = 800)
          write (7114,'(*(G0.6,:","))') "jday","mon","day","yr","unit",     &
     &      "gis_id","name","pfas","tot_in_kg","sol_out_kg","sor_out_kg",   &
     &      "settle_kg","resus_kg","diffuse_kg","bury_kg","water_kg",       &
     &      "benthic_kg","conc_ngL"
        end if
        opened = .true.
      end if

      !! daily print (respect the daily print interval)
      if (pco%day_print == "y" .and. pco%int_day_cur == pco%int_day) then
        iob = sp_ob1%chandeg + jrch - 1
        do ipf = 1, npfas
          jpf = pfas_num(ipf)
          write (7110,100) time%day, time%mo, time%day_mo, time%yrc, jrch,  &
     &      ob(iob)%gis_id, ob(iob)%name, pfasdb(jpf)%name,                 &
     &      chpfas_d(jrch)%pfas(ipf)%tot_in,  chpfas_d(jrch)%pfas(ipf)%sol_out,  &
     &      chpfas_d(jrch)%pfas(ipf)%sor_out, chpfas_d(jrch)%pfas(ipf)%settle,   &
     &      chpfas_d(jrch)%pfas(ipf)%resus,   chpfas_d(jrch)%pfas(ipf)%difus,    &
     &      chpfas_d(jrch)%pfas(ipf)%bury,    chpfas_d(jrch)%pfas(ipf)%water,    &
     &      chpfas_d(jrch)%pfas(ipf)%benthic, chpfas_d(jrch)%pfas(ipf)%conc
          if (pco%csvout == "y") then
            write (7114,'(*(G0.6,:","))') time%day, time%mo, time%day_mo,   &
     &        time%yrc, jrch, ob(iob)%gis_id, ob(iob)%name, pfasdb(jpf)%name, &
     &        chpfas_d(jrch)%pfas(ipf)%tot_in,  chpfas_d(jrch)%pfas(ipf)%sol_out, &
     &        chpfas_d(jrch)%pfas(ipf)%sor_out, chpfas_d(jrch)%pfas(ipf)%settle,  &
     &        chpfas_d(jrch)%pfas(ipf)%resus,   chpfas_d(jrch)%pfas(ipf)%difus,   &
     &        chpfas_d(jrch)%pfas(ipf)%bury,    chpfas_d(jrch)%pfas(ipf)%water,   &
     &        chpfas_d(jrch)%pfas(ipf)%benthic, chpfas_d(jrch)%pfas(ipf)%conc
          end if
        end do
      end if

100   format (4i6,2i9,2x,a16,2x,a16,10es15.6)
1001  format (4a6,2a9,2x,a16,2x,a16,10a15)

      return
      end subroutine cha_pfas_output
