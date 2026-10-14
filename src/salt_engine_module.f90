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

      implicit none

      logical :: road_timing_hru_gauge = .false.        !! E3 key set in salt_engine.ctl

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
        case default
          write (*, "(a,i0,2a)") "salt_engine.ctl line ", iline, ": unknown key (known: road_timing): ", trim(line)
          error stop "salt_engine.ctl: unknown key"
        end select
      end do
      close (iu)
      end subroutine salt_engine_ctl_read

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

      end module salt_engine_module
