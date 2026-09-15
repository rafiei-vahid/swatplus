      subroutine pfas_dep_read

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    RESEARCH BUILD (Q645 year-by-year deposition deck, 2026-09-15; base 46af1a8).
!!    Reads the OPTIONAL atmospheric-deposition schedule pfas_dep.dat once, at the
!!    start of the run (called from pfas_read after pfas_calib.dat).
!!
!!    ABSENT file -> pfas_dep_on stays 0, nothing is allocated, and no other code
!!    path of the PFAS land phase changes: pfas_dep returns at its first statement
!!    and pfas_output takes its original branches (the byte-identity control
!!    against the base engine).
!!
!!    File layout (list-directed, one record per line):
!!      line 1   title (free text, one token or quoted)
!!      line 2   mode daily        or   mode year_start
!!      line 3   header (free text, e.g. year pfasid flux_ug_ha_yr)
!!      then     <calendar year> <pfasid> <flux: micrograms PFAS per hectare per year>
!!    Years not listed deposit nothing. The flux is basin-uniform (every HRU) and
!!    enters the TOP soil layer. It is NOT multiplied by pfas_soil_scale: that
!!    multiplier scales only the initial pools read from pfas_hru.ini, so a flux
!!    here is an absolute input (the DEP-UB lesson: an ini increment was scaled).
!!    Any malformed line stops the run with a message (research input, fail loud).
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      use pfas_module
      use pfas_output_module
      use hydrograph_module, only : sp_ob

      implicit none

      character (len=80) :: titldum = ""   !!          |title line
      character (len=80) :: header = ""    !!          |header line
      character (len=20) :: ctok = ""      !!          |"mode" keyword
      character (len=20) :: cmode = ""     !!          |daily | year_start
      logical :: i_exist = .false.         !!          |file-existence flag
      integer :: eof = 0                   !!          |read status
      integer :: n = 0                     !! none     |number of schedule records
      integer :: i = 0                     !! none     |record counter
      integer :: nhru = 0                  !! none     |number of HRUs
      integer :: yr = 0                    !! yr       |scratch year
      integer :: pid = 0                   !! none     |scratch pfasid
      real(8) :: fl = 0.d0                 !! ug/ha/yr |scratch flux
      real(8) :: total = 0.d0              !! ug/ha    |schedule sum (printed)
      integer, dimension(:), allocatable :: ryr, rpid
      real(8), dimension(:), allocatable :: rfl

      pfas_dep_on = 0
      if (npfas <= 0) return
      inquire (file="pfas_dep.dat", exist=i_exist)
      if (.not. i_exist) return

      open (173, file="pfas_dep.dat")
      read (173,*,iostat=eof) titldum
      if (eof /= 0) call dep_fail("missing title line")
      read (173,*,iostat=eof) ctok, cmode
      if (eof /= 0 .or. trim(ctok) /= "mode") call dep_fail("line 2 must read: mode daily|year_start")
      select case (trim(cmode))
      case ("daily")
        pfas_dep_mode = 1
      case ("year_start")
        pfas_dep_mode = 2
      case default
        call dep_fail("unknown mode (expected daily or year_start)")
      end select
      read (173,*,iostat=eof) header
      if (eof /= 0) call dep_fail("missing header line")

      !! pass 1: count records
      n = 0
      do
        read (173,*,iostat=eof) yr, pid, fl
        if (eof < 0) exit
        if (eof > 0) call dep_fail("unreadable record (expected: year pfasid flux)")
        n = n + 1
      end do
      if (n == 0) call dep_fail("no deposition records")

      !! pass 2: store and validate
      allocate (ryr(n), rpid(n), rfl(n))
      rewind (173)
      read (173,*) titldum
      read (173,*) ctok, cmode
      read (173,*) header
      do i = 1, n
        read (173,*) ryr(i), rpid(i), rfl(i)
        if (rpid(i) < 1 .or. rpid(i) > npfas) call dep_fail("pfasid outside 1..npfas")
        if (.not. (rfl(i) >= 0.d0)) call dep_fail("flux must be >= 0 (and a number)")
      end do
      close (173)

      pfas_dep_yr0 = minval(ryr)
      pfas_dep_yr1 = maxval(ryr)
      allocate (pfas_dep_flux(pfas_dep_yr0:pfas_dep_yr1, npfas), source=-1.d0)
      do i = 1, n
        if (pfas_dep_flux(ryr(i), rpid(i)) >= 0.d0) call dep_fail("duplicate year/pfasid record")
        pfas_dep_flux(ryr(i), rpid(i)) = rfl(i)
      end do
      where (pfas_dep_flux < 0.d0) pfas_dep_flux = 0.d0
      total = sum(rfl)

      nhru = max(sp_ob%hru, 1)
      allocate (pfas_dep_lastyr(nhru), source=0)
      allocate (pfas_dep_target(nhru, npfas), source=0.d0)
      allocate (pfas_dep_applied(nhru, npfas), source=0.d0)
      pfas_dep_on = 1

      write (*,'(a,a,a,i0,a,i0,a,i0,a,es14.6,a)') " PFAS deposition ON (pfas_dep.dat): mode ", &
        trim(cmode), ", years ", pfas_dep_yr0, "-", pfas_dep_yr1, ", records ", n,           &
        ", schedule sum ", total, " ug/ha, not scaled by pfas_soil_scale"
      deallocate (ryr, rpid, rfl)
      return

      contains

        subroutine dep_fail(msg)
          character (len=*), intent(in) :: msg
          write (*,'(a,a)') " pfas_dep.dat ERROR: ", msg
          stop 1
        end subroutine dep_fail

      end subroutine pfas_dep_read


      subroutine pfas_dep

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    RESEARCH BUILD (Q645). Atmospheric deposition to the top-layer PFAS pool of
!!    the current HRU (ihru). Called once per HRU per day from hru_control, just
!!    before pfas_lch, so the day's deposit takes part in that day's partitioning,
!!    runoff, erosion and leaching. No-op unless pfas_dep.dat was read.
!!
!!      mode daily      : flux(yr) / days-in-year (365 or 366, the engine's own
!!                        calendar, ndays(13)) on every simulated day of year yr
!!      mode year_start : flux(yr) on the first simulated day of year yr
!!
!!    Compensated summation: the double-precision target accumulates the exact
!!    schedule; the single-precision pool receives (target - applied), and
!!    "applied" records what actually entered the pool after rounding. A daily
!!    increment on a large surface pool is only ~10 units in the last place of a
!!    4-byte real, so plain repeated addition could bias the deposited mass by
!!    per cent; with the feedback the applied total stays within one rounding
!!    unit of the target. pfas_output closes the mass balance on "applied" and
!!    reports the target beside it.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      use hru_module, only : ihru
      use time_module
      use pfas_module
      use pfas_output_module

      implicit none

      integer :: j = 0          !! none   |HRU number
      integer :: k = 0          !! none   |sequential PFAS counter
      real :: old = 0.          !! kg/ha  |top-layer pool before the deposit
      real(8) :: add = 0.d0     !! kg/ha  |today's scheduled deposit

      if (pfas_dep_on == 0) return
      if (npfas == 0) return
      j = ihru
      if (time%yrc < pfas_dep_yr0 .or. time%yrc > pfas_dep_yr1) return
      if (.not. allocated(pfas_soil_hru(j)%ly)) return
      if (pfas_dep_mode == 2) then
        if (pfas_dep_lastyr(j) == time%yrc) return
        pfas_dep_lastyr(j) = time%yrc
      end if

      do k = 1, npfas
        if (pfas_dep_mode == 1) then
          add = pfas_dep_flux(time%yrc, k) * 1.d-9 / dble(ndays(13))
        else
          add = pfas_dep_flux(time%yrc, k) * 1.d-9
        end if
        if (add <= 0.d0) cycle
        pfas_dep_target(j,k) = pfas_dep_target(j,k) + add
        old = pfas_soil_hru(j)%ly(1)%sol_pfas(k)
        pfas_soil_hru(j)%ly(1)%sol_pfas(k) =                                   &
          real(dble(old) + (pfas_dep_target(j,k) - pfas_dep_applied(j,k)))
        pfas_dep_applied(j,k) = pfas_dep_applied(j,k)                          &
          + (dble(pfas_soil_hru(j)%ly(1)%sol_pfas(k)) - dble(old))
        pfas_flag(j) = 1
      end do

      return
      end subroutine pfas_dep
