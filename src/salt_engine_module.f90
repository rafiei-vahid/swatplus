      module salt_engine_module

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    optional salt-engine switches, read from salt_engine.ctl in the run folder ONLY when that file exists. With the file
!!    absent nothing in this module runs and no output changes. One key per change; an unknown key or value stops the run.
!!
!!    road_timing hru_gauge (E3): each simulated year's yearly road-salt load is spread over the days of the HRU's OWN weather
!!    station gauges (its pgage and tgage), not over the deposition station's index used as a gauge index, and every simulated
!!    year receives its whole load:
!!      Q  days with precipitation > 0 and (tmax+tmin)/2 < 0, weighted by precipitation;
!!      C  if no Q day: days with (tmax+tmin)/2 < 0, equal weights;
!!      W  if no C day: 1 January - 31 March and 1 - 31 December of that year, equal weights.
!!    A day with any of its three values at -99 is missing and never qualifies for Q or C. Only the simulated days of a year are
!!    used (the engine holds no record before the first simulated day), so a partial year receives its whole load inside its
!!    simulated days (a convention). The road-salt year is taken from the calendar, not from atmodep_cont%ts.
!!
!!    salt_chemistry off (E6): salt_chem_hru and salt_chem_aqu keep only their negative-mass clamp and the concentration refresh;
!!    the precipitation-dissolution package, the cation exchange and the write-back are skipped, so every ion is moved only by
!!    the fluxes, the minerals keep their fractions and no dissolution is booked. `salt_chemistry on` equals an absent key.
!!    salt_chem_soil_single (the plant salt-stress estimate) is not affected.
!!
!!    road_load hru_file (E1, needs road_timing hru_gauge): each HRU's yearly road salt per ion comes from salt_road_hru (kg/ha;
!!    line 2 'yr_init num', line 3 'hru year' + the ion names in constituents.cs order, then one row per HRU and year) instead
!!    of its deposition station's salt_road load; timing stays E3's. An HRU without rows gets zero and is counted.
!!
!!    road_washoff impervious (E2, needs road_timing hru_gauge, daily runs): on an urban HRU the share frac_dc_imp
!!    (urban.urb, urbdb%fcimp) of each day's road salt goes to an impervious store instead of soil layer 1; on a day
!!    hru_urban washes urban solids off, the same fraction 1 - Exp(-xx) of the store washes into urbqsalt (build-up/wash-off),
!!    or the whole store on a USGS-regression runoff day. The store and the washed mass are printed as extra columns of the
!!    HRU and basin salt files only with the key.
!!
!!    road_washoff_scale m (E2.4, needs road_washoff impervious): the store's share becomes min(1, m x frac_dc_imp), m a real in
!!    [0, 2]; without the key m = 1, which is E2 itself.
!!
!!    aquifer_salt hold_initial (E5): each aquifer's salt concentration stays at its initial value C0. Recharge salt is booked out
!!    (rchout) instead of added; the mass is set to C0 x V before discharge and seepage and again after them, the two settings
!!    printed as a signed 'hold'. Columns rchout_<ion>, hold_<ion> (aquifer files) and rchout_kg, hold_kg (basin) only with the key.

      implicit none

      logical :: road_timing_hru_gauge = .false.        !! E3 key set in salt_engine.ctl
      logical :: salt_chemistry_off = .false.           !! E6 key: salt_chemistry off (no equilibrium chemistry, no exchange)
      logical :: road_load_hru_file = .false.           !! E1 key: road_load hru_file (each HRU's own yearly load, salt_road_hru)

      real, dimension(:,:,:), allocatable :: rdhru      !! E1 (ion, hru, series year) road salt, kg/ha
      logical, dimension(:), allocatable :: rdhru_has   !! E1 (hru) the file gives this HRU a load
      integer :: rdhru_yr0 = 0                          !! E1 first calendar year of the file's series
      integer :: rdhru_num = 0                          !! E1 number of years in the file's series

      logical :: road_washoff_imp = .false.             !! E2 key: road_washoff impervious
      integer :: road_washoff_lines = 0                 !! E2 road_washoff lines in salt_engine.ctl (more than one stops)
      real :: road_washoff_m = 1.                       !! E2.4 road_washoff_scale m (1 without the key)
      integer :: road_washoff_scale_lines = 0           !! E2.4 road_washoff_scale lines (more than one stops)
      real, dimension(:,:), allocatable :: rdimp        !! E2 (ion, hru) road salt on directly connected impervious area, kg/ha
      integer :: n_e2 = 0                               !! E2 ions printed per extra HRU column group (0 without the key)
      integer :: n_e2b = 0                              !! E2 extra basin columns (0 without the key, else 2)
      character(len=16), dimension(16) :: e2_hdr_hru = ""  !! E2 HRU header names: imps_<ion> then impw_<ion>
      character(len=16), dimension(2) :: e2_hdr_bsn = (/"        impst_kg", "        impwo_kg"/)
      real, dimension(2) :: impb = 0.                   !! E2 basin, day: store (kg), washed (kg)
      real, dimension(2) :: impb_mo = 0.
      real, dimension(2) :: impb_yr = 0.
      real, dimension(2) :: impb_aa = 0.

      logical :: aqu_hold_c0 = .false.                  !! E5 key: aquifer_salt hold_initial
      integer :: aqu_hold_lines = 0                     !! E5 aquifer_salt lines (more than one stops)
      real, dimension(:,:), allocatable :: aqu_c0       !! E5 (ion, aquifer) initial concentration, g/m3
      logical, dimension(:), allocatable :: aqu_c0_set  !! E5 (aquifer) C0 recorded from its initial salt record
      integer :: n_e5 = 0                               !! E5 ions per extra aquifer column group (0 without the key)
      integer :: n_e5b = 0                              !! E5 extra basin columns (0 without the key, else 2)
      character(len=16), dimension(16) :: e5_hdr_aqu = ""  !! E5 aquifer header names: rchout_<ion> then hold_<ion>
      character(len=16), dimension(2) :: e5_hdr_bsn = (/"       rchout_kg", "         hold_kg"/)
      real, dimension(2) :: hldb = 0.                   !! E5 basin, day: rchout (kg), hold (kg)
      real, dimension(2) :: hldb_mo = 0.
      real, dimension(2) :: hldb_yr = 0.
      real, dimension(2) :: hldb_aa = 0.

      type road_timing_weights
        real, dimension(:,:), allocatable :: w                  !! (day, road-salt year) fraction of the yearly load
        character(len=1), dimension(:), allocatable :: method   !! (road-salt year) Q, C or W; - for a year not simulated
      end type road_timing_weights
      type (road_timing_weights), dimension(:), allocatable :: rdw   !! per weather station
      logical, dimension(:), allocatable :: rdw_nogauge              !! station without a measured pcp or tmp gauge

      contains

      subroutine salt_engine_ctl_read
!!    reads salt_engine.ctl if it exists: one `key value` per line; blank lines and lines starting with ! or # are skipped
      logical :: i_exist
      integer :: iu
      integer :: eof
      integer :: iline
      character(len=256) :: line
      character(len=64) :: key
      character(len=64) :: val

      inquire (file="salt_engine.ctl", exist=i_exist)
      if (.not. i_exist) return
      open (newunit=iu, file="salt_engine.ctl", status="old", action="read")
      iline = 0
      do
        read (iu, "(a)", iostat=eof) line
        if (eof /= 0) exit
        iline = iline + 1
        line = adjustl(line)
        if (len_trim(line) == 0) cycle
        if (line(1:1) == "!" .or. line(1:1) == "#") cycle
        key = ""
        val = ""
        read (line, *, iostat=eof) key, val
        if (eof /= 0 .or. len_trim(val) == 0) then
          write (*, "(a,i0,2a)") "salt_engine.ctl line ", iline, ": expected 'key value', read: ", trim(line)
          error stop "salt_engine.ctl: unreadable line"
        end if
        select case (trim(key))
        case ("road_timing")
          if (trim(val) == "hru_gauge") then
            road_timing_hru_gauge = .true.
          else
            write (*, "(a,i0,2a)") "salt_engine.ctl line ", iline, ": unknown value for road_timing (known: hru_gauge): ", &
                                   trim(line)
            error stop "salt_engine.ctl: unknown value"
          end if
        case ("salt_chemistry")
          if (trim(val) == "off") then
            salt_chemistry_off = .true.
          else if (trim(val) /= "on") then
            write (*, "(a,i0,2a)") "salt_engine.ctl line ", iline, ": unknown value for salt_chemistry (known: off, on): ", &
                                   trim(line)
            error stop "salt_engine.ctl: unknown value"
          end if
        case ("road_load")
          if (trim(val) == "hru_file") then
            road_load_hru_file = .true.
          else
            write (*, "(a,i0,2a)") "salt_engine.ctl line ", iline, ": unknown value for road_load (known: hru_file): ", &
                                   trim(line)
            error stop "salt_engine.ctl: unknown value"
          end if
        case ("road_washoff")
          road_washoff_lines = road_washoff_lines + 1
          if (road_washoff_lines > 1) then
            write (*, "(a,i0,2a)") "salt_engine.ctl line ", iline, ": road_washoff given twice: ", trim(line)
            error stop "salt_engine.ctl: road_washoff given twice"
          end if
          if (trim(val) == "impervious") then
            road_washoff_imp = .true.
          else
            write (*, "(a,i0,2a)") "salt_engine.ctl line ", iline, ": unknown value for road_washoff (known: impervious): ", &
                                   trim(line)
            error stop "salt_engine.ctl: unknown value"
          end if
        case ("road_washoff_scale")
          road_washoff_scale_lines = road_washoff_scale_lines + 1
          if (road_washoff_scale_lines > 1) then
            write (*, "(a,i0,2a)") "salt_engine.ctl line ", iline, ": road_washoff_scale given twice: ", trim(line)
            error stop "salt_engine.ctl: road_washoff_scale given twice"
          end if
          read (val, *, iostat=eof) road_washoff_m
          if (eof /= 0) then
            write (*, "(a,i0,2a)") "salt_engine.ctl line ", iline, ": road_washoff_scale is not a number: ", trim(line)
            error stop "salt_engine.ctl: road_washoff_scale is not a number"
          end if
          if (road_washoff_m < 0. .or. road_washoff_m > 2.) then
            write (*, "(a,i0,2a)") "salt_engine.ctl line ", iline, ": road_washoff_scale outside [0, 2]: ", trim(line)
            error stop "salt_engine.ctl: road_washoff_scale outside [0, 2]"
          end if
        case ("aquifer_salt")
          aqu_hold_lines = aqu_hold_lines + 1
          if (aqu_hold_lines > 1) then
            write (*, "(a,i0,2a)") "salt_engine.ctl line ", iline, ": aquifer_salt given twice: ", trim(line)
            error stop "salt_engine.ctl: aquifer_salt given twice"
          end if
          if (trim(val) == "hold_initial") then
            aqu_hold_c0 = .true.
          else
            write (*, "(a,i0,2a)") "salt_engine.ctl line ", iline, ": unknown value for aquifer_salt (known: hold_initial): ", &
                                   trim(line)
            error stop "salt_engine.ctl: unknown value"
          end if
        case default
          write (*, "(a,i0,2a)") "salt_engine.ctl line ", iline,                                                        &
                                 ": unknown key (known: road_timing, salt_chemistry, road_load, road_washoff, road_washoff_scale, " // &
                                 "aquifer_salt): ", &
                                 trim(line)
          error stop "salt_engine.ctl: unknown key"
        end select
      end do
      close (iu)
      end subroutine salt_engine_ctl_read

      subroutine salt_imp_init
!!    E2 init (serial, from salt_hru_init): the impervious store, zero at the start; every urban HRU's frac_dc_imp must lie
!!    in [0, 1]
      use hru_module, only : hru
      use hydrograph_module, only : sp_ob
      use urban_data_module, only : urbdb
      use constituent_mass_module, only : cs_db

      integer :: j
      integer :: ulu

      allocate (rdimp(cs_db%num_salts, sp_ob%hru), source = 0.)
      do j = 1, sp_ob%hru
        ulu = hru(j)%luse%urb_lu
        if (ulu > 0) then
          if (urbdb(ulu)%fcimp < 0. .or. urbdb(ulu)%fcimp > 1.) then
            write (*, "(a,i0,a,f0.6)") "road_washoff impervious: HRU ", j, " has frac_dc_imp outside [0, 1]: ", urbdb(ulu)%fcimp
            error stop "road_washoff impervious: frac_dc_imp outside [0, 1]"
          end if
        end if
      end do
      end subroutine salt_imp_init

      subroutine salt_road_timing_build
!!    E3 read phase (serial): the weights of every weather station with measured gauges, for every simulated year, plus
!!    salt_road_timing.txt. Called from salt_roadsalt_read for the yearly form only, before the base's own partition.
      use climate_module, only : wst, pcp, tmp, pcp_n, tmp_n, atmodep_cont
      use maximum_data_module, only : db_mx
      use time_module, only : time

      integer :: s
      integer :: gp
      integer :: gt
      integer :: iyr
      integer :: iroad
      integer :: kp
      integer :: kt
      integer :: iday
      integer :: day_lo
      integer :: day_hi
      integer :: ndays_yr
      integer :: leap
      integer :: iu
      real :: p
      real :: tx
      real :: tn
      real :: psum
      logical, dimension(366) :: is_q
      logical, dimension(366) :: is_c
      logical, dimension(366) :: is_w

      !! every simulated year needs a road-salt value (no silent zero year)
      do iyr = time%yrc_start, time%yrc_end
        iroad = iyr - atmodep_cont%yr_init + 1
        if (iroad < 1 .or. iroad > atmodep_cont%num) then
          write (*, "(a,i0,a,i0,a,i0)") "road_timing hru_gauge: salt_road has no value for simulated year ", iyr, &
                                        "; its yearly series covers ", atmodep_cont%yr_init, " - ",                 &
                                        atmodep_cont%yr_init + atmodep_cont%num - 1
          error stop "road_timing hru_gauge: a simulated year has no road-salt value"
        end if
      end do

      allocate (rdw(db_mx%wst))
      allocate (rdw_nogauge(db_mx%wst), source = .false.)
      open (newunit=iu, file="salt_road_timing.txt", status="replace", action="write")
      write (iu, "(a)") "salt_road_timing.txt: road_timing hru_gauge (salt_engine.ctl); one row per weather station and" // &
                        " simulated year; method Q snow days, C cold days, W winter window"
      write (iu, "(a)") "station pcp_gage tmp_gage year method days w_sum"

      do s = 1, db_mx%wst
        gp = wst(s)%wco%pgage
        gt = wst(s)%wco%tgage
        if (gp == 0 .or. gt == 0) then
          !! generated weather: the stop is raised in salt_road_timing_check, after the HRUs are tied to their stations
          rdw_nogauge(s) = .true.
          cycle
        end if
        if (pcp(gp)%tstep > 0) then
          write (*, "(4a)") "road_timing hru_gauge: station ", trim(wst(s)%name), " has a sub-daily precipitation gauge ", &
                            trim(pcp_n(gp))
          error stop "road_timing hru_gauge needs daily precipitation gauges"
        end if
        allocate (rdw(s)%w(366, atmodep_cont%num), source = 0.)
        allocate (rdw(s)%method(atmodep_cont%num), source = "-")

        do iyr = time%yrc_start, time%yrc_end
          !! no read past either gauge's record (F4)
          if (iyr < pcp(gp)%start_yr .or. iyr > pcp(gp)%start_yr + pcp(gp)%nbyr - 1) then
            write (*, "(4a,i0)") "road_timing hru_gauge: precipitation gauge ", trim(pcp_n(gp)), " (station ",  &
                                 trim(wst(s)%name) // ") has no record for simulated year ", iyr
            error stop "road_timing hru_gauge: a gauge record does not cover the simulation"
          end if
          if (iyr < tmp(gt)%start_yr .or. iyr > tmp(gt)%start_yr + tmp(gt)%nbyr - 1) then
            write (*, "(4a,i0)") "road_timing hru_gauge: temperature gauge ", trim(tmp_n(gt)), " (station ",   &
                                 trim(wst(s)%name) // ") has no record for simulated year ", iyr
            error stop "road_timing hru_gauge: a gauge record does not cover the simulation"
          end if
          !! the engine's own column for this year (climate_control.f90:80-81: time%yrs - yrs_start)
          kp = (iyr - time%yrc_start + 1) - pcp(gp)%yrs_start
          kt = (iyr - time%yrc_start + 1) - tmp(gt)%yrs_start
          iroad = iyr - atmodep_cont%yr_init + 1

          leap = 0
          if ((mod(iyr, 4) == 0 .and. mod(iyr, 100) /= 0) .or. mod(iyr, 400) == 0) leap = 1
          ndays_yr = 365 + leap
          day_lo = 1
          if (iyr == time%yrc_start) day_lo = time%day_start
          day_hi = ndays_yr
          if (iyr == time%yrc_end .and. time%day_end > 0) day_hi = min(time%day_end, ndays_yr)

          is_q = .false.
          is_c = .false.
          is_w = .false.
          psum = 0.
          do iday = day_lo, day_hi
            is_w(iday) = iday <= 90 + leap .or. iday >= 335 + leap
            p = pcp(gp)%ts(iday, kp)
            tx = tmp(gt)%ts(iday, kt)
            tn = tmp(gt)%ts2(iday, kt)
            if (p <= -99. .or. tx <= -99. .or. tn <= -99.) cycle
            if ((tx + tn) / 2. < 0.) then
              is_c(iday) = .true.
              if (p > 0.) then
                is_q(iday) = .true.
                psum = psum + p
              end if
            end if
          end do

          if (count(is_q) > 0) then
            rdw(s)%method(iroad) = "Q"
            do iday = day_lo, day_hi
              if (is_q(iday)) rdw(s)%w(iday, iroad) = pcp(gp)%ts(iday, kp) / psum
            end do
          else if (count(is_c) > 0) then
            rdw(s)%method(iroad) = "C"
            where (is_c) rdw(s)%w(:, iroad) = 1. / real(count(is_c))
          else if (count(is_w) > 0) then
            rdw(s)%method(iroad) = "W"
            where (is_w) rdw(s)%w(:, iroad) = 1. / real(count(is_w))
          else
            write (*, "(3a,i0,a,i0,a,i0,a)") "road_timing hru_gauge: station ", trim(wst(s)%name),           &
                                             " has no snow, cold or winter-window day in simulated year ", iyr, &
                                             " (simulated days ", day_lo, " - ", day_hi, ")"
            error stop "road_timing hru_gauge: a simulated year has no day to receive its road salt"
          end if
          write (iu, "(a,1x,a,1x,a,1x,i0,1x,a1,1x,i0,1x,f11.8)") trim(wst(s)%name), trim(pcp_n(gp)), trim(tmp_n(gt)), &
                 iyr, rdw(s)%method(iroad), count(rdw(s)%w(:, iroad) > 0.), sum(rdw(s)%w(:, iroad))
        end do
      end do
      close (iu)
      end subroutine salt_road_timing_build

      subroutine salt_road_timing_check
!!    E3, after hyd_connect and hrudb_init: a weather station without measured gauges stops the run, naming every such station
!!    with the number of HRUs that use it (0 = the station can be dropped; otherwise give it gauges).
      use climate_module, only : wst
      use hydrograph_module, only : ob, sp_ob
      use hru_module, only : hru

      integer :: j
      integer :: s
      integer, dimension(:), allocatable :: nhru

      if (.not. allocated(rdw_nogauge)) return
      if (.not. any(rdw_nogauge)) return
      allocate (nhru(size(rdw_nogauge)), source = 0)
      do j = 1, sp_ob%hru
        s = ob(hru(j)%obj_no)%wst
        if (s > 0 .and. s <= size(rdw_nogauge)) then
          if (rdw_nogauge(s)) nhru(s) = nhru(s) + 1
        end if
      end do
      do s = 1, size(rdw_nogauge)
        if (rdw_nogauge(s)) write (*, "(7a,i0,a)") "road_timing hru_gauge: station ", trim(wst(s)%name),            &
            " has no measured gauge (pgage ", trim(wst(s)%wco_c%pgage), ", tgage ", trim(wst(s)%wco_c%tgage),     &
            "), used by ", nhru(s), " HRUs"
      end do
      error stop "road_timing hru_gauge needs measured precipitation and temperature gauges at every weather station"
      end subroutine salt_road_timing_check

      subroutine salt_road_hru_read
!!    E1, in salt_hru_init after hyd_connect and hrudb_init (serial): reads salt_road_hru, checks it loudly, and reports the HRUs
!!    without a load and the file / station ratio per station and simulated year to simulation.out and salt_road_hru_log.txt.
      use constituent_mass_module, only : cs_db
      use hydrograph_module, only : ob, sp_ob
      use hru_module, only : hru
      use climate_module, only : wst, rdapp_salt, atmodep_cont
      use time_module, only : time

      integer :: iu
      integer :: eof
      integer :: iline
      integer :: nion
      integer :: i
      integer :: j
      integer :: k
      integer :: ih
      integer :: iyr
      integer :: icl
      integer :: iadep
      integer :: iroad
      integer :: nmiss
      integer :: nhas
      integer :: nzero
      integer :: lg
      logical :: i_exist
      character(len=512) :: line
      character(len=16), dimension(:), allocatable :: tok
      real, dimension(:), allocatable :: vals
      logical, dimension(:,:), allocatable :: seen
      real :: a_miss
      real :: a_all
      real :: a_zero
      real :: f_tot
      real :: s_tot

      inquire (file="salt_road_hru", exist=i_exist)
      if (.not. i_exist) then
        write (*,"(a)") "salt_engine.ctl road_load hru_file: there is no salt_road_hru file in the run folder"
        error stop "road_load hru_file needs salt_road_hru"
      end if
      nion = cs_db%num_salts
      allocate (tok(2 + nion))
      allocate (vals(nion))
      open (newunit=iu, file="salt_road_hru", status="old", action="read")
      read (iu, "(a)", iostat=eof) line
      read (iu, *, iostat=eof) rdhru_yr0, rdhru_num
      if (eof /= 0 .or. rdhru_num < 1) then
        write (*,"(a)") "salt_road_hru line 2: expected 'yr_init num' with num >= 1"
        error stop "salt_road_hru: unreadable series line"
      end if
      read (iu, "(a)", iostat=eof) line
      tok = ""
      read (line, *, iostat=eof) tok
      if (n_tokens(line) /= 2 + nion) then
        write (*,"(a,i0,a,i0,2a)") "salt_road_hru line 3: ", n_tokens(line), " names, expected ", 2 + nion, &
                                    " ('hru year' + the ion names): ", trim(line)
        error stop "salt_road_hru: header has the wrong number of names"
      end if
      if (eof /= 0 .or. trim(tok(1)) /= "hru" .or. trim(tok(2)) /= "year") then
        write (*,"(2a)") "salt_road_hru line 3: expected 'hru year' then the ion names, read: ", trim(line)
        error stop "salt_road_hru: bad header"
      end if
      do i = 1, nion
        if (trim(tok(2 + i)) /= trim(cs_db%salts(i))) then
          write (*,"(a,i0,5a)") "salt_road_hru line 3: ion ", i, " is '", trim(tok(2 + i)), "', constituents.cs has '", &
                                trim(cs_db%salts(i)), "'"
          error stop "salt_road_hru: ion names differ from constituents.cs"
        end if
      end do
      do iyr = time%yrc_start, time%yrc_end
        if (iyr < rdhru_yr0 .or. iyr > rdhru_yr0 + rdhru_num - 1) then
          write (*,"(a,i0,a,i0,a,i0)") "salt_road_hru has no year ", iyr, " (simulated); its series covers ", rdhru_yr0, &
                                       " - ", rdhru_yr0 + rdhru_num - 1
          error stop "salt_road_hru: a simulated year is outside its series"
        end if
      end do

      allocate (rdhru(nion, sp_ob%hru, rdhru_num), source = 0.)
      allocate (rdhru_has(sp_ob%hru), source = .false.)
      allocate (seen(sp_ob%hru, rdhru_num), source = .false.)
      iline = 3
      do
        read (iu, "(a)", iostat=eof) line
        if (eof /= 0) exit
        iline = iline + 1
        if (len_trim(line) == 0) cycle
        !! a list-directed read ignores trailing extras, so a row whose columns drifted at the right edge would pass (Lane S)
        if (n_tokens(line) /= 2 + nion) then
          write (*,"(a,i0,a,i0,a,i0,2a)") "salt_road_hru line ", iline, ": ", n_tokens(line), " values, expected ", 2 + nion, &
                                          ": ", trim(line)
          error stop "salt_road_hru: a row has the wrong number of values"
        end if
        read (line, *, iostat=eof) ih, iyr, vals
        if (eof /= 0) then
          write (*,"(a,i0,2a)") "salt_road_hru line ", iline, ": unreadable row: ", trim(line)
          error stop "salt_road_hru: unreadable row"
        end if
        if (ih < 1 .or. ih > sp_ob%hru) then
          write (*,"(a,i0,a,i0,a,i0)") "salt_road_hru line ", iline, ": hru ", ih, " is outside 1 - ", sp_ob%hru
          error stop "salt_road_hru: hru id out of range"
        end if
        if (iyr < rdhru_yr0 .or. iyr > rdhru_yr0 + rdhru_num - 1) then
          write (*,"(a,i0,a,i0,a)") "salt_road_hru line ", iline, ": year ", iyr, " is outside the series on line 2"
          error stop "salt_road_hru: year out of range"
        end if
        k = iyr - rdhru_yr0 + 1
        if (seen(ih, k)) then
          write (*,"(a,i0,a,i0,a,i0)") "salt_road_hru line ", iline, ": hru ", ih, " year ", iyr, " appears twice"
          error stop "salt_road_hru: duplicated row"
        end if
        if (any(vals < 0.)) then
          write (*,"(a,i0,a)") "salt_road_hru line ", iline, ": a negative load"
          error stop "salt_road_hru: negative load"
        end if
        seen(ih, k) = .true.
        rdhru_has(ih) = .true.
        rdhru(:, ih, k) = vals
      end do
      close (iu)
      do j = 1, sp_ob%hru
        if (rdhru_has(j) .and. .not. all(seen(j, :))) then
          k = findloc(seen(j, :), .false., dim=1)
          write (*,"(a,i0,a,i0)") "salt_road_hru: hru ", j, " has rows but none for year ", rdhru_yr0 + k - 1
          error stop "salt_road_hru: an HRU lacks a year of the series"
        end if
      end do

      !! report: the HRUs without a load, and the file / station Cl ratio per station and simulated year
      nhas = count(rdhru_has)
      nmiss = sp_ob%hru - nhas
      a_all = 0.
      a_miss = 0.
      a_zero = 0.
      nzero = 0
      do j = 1, sp_ob%hru
        a_all = a_all + hru(j)%area_ha
        if (.not. rdhru_has(j)) a_miss = a_miss + hru(j)%area_ha
        !! a builder writes a row for every HRU, zeros included, so rows alone cannot show a class written as zero (Lane SD)
        if (rdhru_has(j) .and. all(rdhru(:, j, :) == 0.)) then
          nzero = nzero + 1
          a_zero = a_zero + hru(j)%area_ha
        end if
      end do
      write (9003,"(a,i0,a,i0,a,i0,a,f10.8,a,i0,a,f10.8)") "salt_road_hru: HRUs with rows ", nhas, " of ", sp_ob%hru,     &
            "; without rows ", nmiss, ", area share ", a_miss / max(a_all, 1.e-12), "; with rows but zero load ", nzero,  &
            ", area share ", a_zero / max(a_all, 1.e-12)
      open (newunit=lg, file="salt_road_hru_log.txt", status="replace", action="write")
      write (lg,"(a,i0,a,i0,a,i0,a,f10.8,a,i0,a,f10.8)") "salt_road_hru: HRUs with rows ", nhas, " of ", sp_ob%hru,       &
            "; without rows ", nmiss, ", area share ", a_miss / max(a_all, 1.e-12), "; with rows but zero load ", nzero,  &
            ", area share ", a_zero / max(a_all, 1.e-12)
      k = 0
      do j = 1, sp_ob%hru
        if (.not. rdhru_has(j) .and. k < 50) then
          write (lg,"(a,i0)") "without rows: hru ", j
          k = k + 1
        end if
      end do
      k = 0
      do j = 1, sp_ob%hru
        if (rdhru_has(j) .and. all(rdhru(:, j, :) == 0.) .and. k < 50) then
          write (lg,"(a,i0)") "zero load: hru ", j
          k = k + 1
        end if
      end do
      icl = 0
      do i = 1, nion
        if (trim(cs_db%salts(i)) == "cl") icl = i
      end do
      if (icl == 0) then
        write (9003,"(a)") "salt_road_hru: no 'cl' ion in constituents.cs, so no file / station ratio is reported"
      else
        do iadep = 1, atmodep_cont%num_sta
          do iyr = time%yrc_start, time%yrc_end
            iroad = iyr - atmodep_cont%yr_init + 1
            k = iyr - rdhru_yr0 + 1
            f_tot = 0.
            s_tot = 0.
            do j = 1, sp_ob%hru
              if (wst(ob(hru(j)%obj_no)%wst)%wco%atmodep /= iadep) cycle
              f_tot = f_tot + rdhru(icl, j, k) * hru(j)%area_ha
              s_tot = s_tot + rdapp_salt(iadep)%salt(icl)%roadyr(iroad) * hru(j)%area_ha
            end do
            if (s_tot > 0.) then
              write (9003,"(a,i0,1x,i0,1x,es14.7)") "salt_road_hru ratio ", iadep, iyr, f_tot / s_tot
              write (lg,"(a,i0,1x,i0,1x,es14.7)") "salt_road_hru ratio ", iadep, iyr, f_tot / s_tot
            end if
          end do
        end do
      end if
      close (lg)
      end subroutine salt_road_hru_read

      integer function n_tokens(s)
!!    number of blank- or tab-separated tokens in s
      character(len=*), intent(in) :: s
      integer :: i
      logical :: in_tok
      n_tokens = 0
      in_tok = .false.
      do i = 1, len_trim(s)
        if (s(i:i) == " " .or. s(i:i) == char(9)) then
          in_tok = .false.
        else if (.not. in_tok) then
          in_tok = .true.
          n_tokens = n_tokens + 1
        end if
      end do
      end function n_tokens

      end module salt_engine_module
