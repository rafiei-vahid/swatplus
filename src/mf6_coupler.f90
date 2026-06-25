      module mf6_coupler
      !! ------------------------------------------------------------------
      !!  SWAT+ <-> MODFLOW 6 daily two-way coupler (BMI / libmf6.so).
      !!
      !!  Design (see vadose-zone-coupling-gap memory + paper plan):
      !!    - GWF (flow) is stepped DAILY in lockstep with the SWAT+ day loop
      !!      (water table / recharge / baseflow are fast and matter daily).
      !!    - GWT (PFAS transport) is advanced at a coarser (monthly) cadence
      !!      because the plume moves slowly -- the classic MODFLOW->MT3D/RT3D
      !!      flow-then-transport split (Bailey's SWAT-MODFLOW-NWT-RT3D).
      !!
      !!  Activation: presence of "mf6.con" in the SWAT+ TxtInOut.  This keeps
      !!  the engine backward compatible (bsn_cc is read as one record, so a
      !!  new flag there would break every existing codes.bsn).
      !!
      !!  mf6.con format (free, '!' comments ignored):
      !!    line 1: workspace        relative dir holding mfsim.nam (e.g. ./mf6)
      !!    line 2: gwf_cadence      flow step interval, days     (default 1)
      !!    line 3: gwt_cadence      transport step interval, days(default 30)
      !!    line 4: lib path         libmf6.so (default below)    (optional)
      !!
      !!  IMPORTANT -- why dlopen() and not link-time -lmf6:
      !!    libmf6.so statically embeds the Intel Fortran runtime, exporting
      !!    for__* symbols (for__file_info_hash_table, for__aio_global_mutex,
      !!    ...).  Link-time loading puts those in the GLOBAL symbol scope where
      !!    they INTERPOSE the engine's own Fortran I/O runtime -> SWAT+ file
      !!    buffers get freed by the wrong runtime copy -> munmap_chunk crash
      !!    during climate-file reading.  We therefore dlopen libmf6.so with
      !!    RTLD_LOCAL so its runtime stays private and never interposes.
      !! ------------------------------------------------------------------
      use iso_c_binding
      implicit none

      private
      public :: mf6_coupler_init, mf6_coupler_step, mf6_coupler_finalize
      public :: mf6_active, mf6_baseflow_active, mf6_channel_baseflow

      logical :: mf6_active = .false.    !! set .true. once MF6 is initialized
      character(len=256) :: mf6_ws = "." !! MF6 workspace (relative to TxtInOut)
      character(len=512) :: mf6_lib = "/data/SWATGenXApp/codes/bin/libmf6.so"
      integer :: gwf_cadence = 1         !! flow step interval (days)
      integer :: gwt_cadence = 30        !! transport step interval (days)
      integer :: mf6_daycount = 0        !! engine days since coupling start
      real(c_double) :: mf6_tend = 0.0_c_double  !! MF6 end time (days)

      !! ---- dlopen/dlsym (resolve libmf6.so BMI entry points privately) ----
      integer(c_int), parameter :: RTLD_NOW_LOCAL = 2   !! RTLD_NOW(2)|RTLD_LOCAL(0)
      type(c_ptr) :: mf6_handle = c_null_ptr

      interface
        function c_dlopen(filename, mode) bind(c, name="dlopen") result(h)
          import :: c_char, c_int, c_ptr
          character(kind=c_char), intent(in) :: filename(*)
          integer(c_int), value :: mode
          type(c_ptr) :: h
        end function
        function c_dlsym(handle, symbol) bind(c, name="dlsym") result(p)
          import :: c_char, c_ptr, c_funptr
          type(c_ptr), value :: handle
          character(kind=c_char), intent(in) :: symbol(*)
          type(c_funptr) :: p
        end function
        function c_dlclose(handle) bind(c, name="dlclose") result(s)
          import :: c_ptr, c_int
          type(c_ptr), value :: handle
          integer(c_int) :: s
        end function
      end interface

      !! ---- BMI function signatures (C, return int) ------------------------
      abstract interface
        integer(c_int) function bmi_noarg() bind(c)
          import :: c_int
        end function
        integer(c_int) function bmi_dbl(t) bind(c)
          import :: c_int, c_double
          real(c_double), intent(inout) :: t
        end function
        integer(c_int) function bmi_dt(dt) bind(c)   !! prepare_time_step(dt)
          import :: c_int, c_double
          real(c_double), intent(in) :: dt
        end function
        integer(c_int) function bmi_getptr(addr, ptr) bind(c)
          import :: c_int, c_char, c_ptr
          character(kind=c_char), intent(in) :: addr(*)
          type(c_ptr), intent(inout) :: ptr
        end function
      end interface

      procedure(bmi_noarg), pointer :: p_initialize => null()
      procedure(bmi_noarg), pointer :: p_update => null()
      procedure(bmi_noarg), pointer :: p_finalize => null()
      procedure(bmi_dbl),   pointer :: p_get_current_time => null()
      procedure(bmi_dbl),   pointer :: p_get_end_time => null()
      !! XMI granular stepping (lets us inject recharge between prepare & solve)
      procedure(bmi_dt),    pointer :: p_prepare_time_step => null()
      procedure(bmi_noarg), pointer :: p_do_time_step => null()
      procedure(bmi_noarg), pointer :: p_finalize_time_step => null()
      procedure(bmi_getptr),pointer :: p_get_value_ptr_double => null()

      !! ---- recharge down-coupling (M2): area-weighted HRU -> MF6 cell map --
      integer :: n_map = 0          !! number of (cell,hru) overlap entries
      integer :: n_rch = 0          !! size of MF6 RECHARGE array (NROW*NCOL)
      integer, allocatable :: map_idx(:)   !! 0-based RECHARGE index per entry
      integer, allocatable :: map_hru(:)   !! SWAT+ HRU id per entry
      real,    allocatable :: map_w(:)     !! overlap_area / cell_area weight
      real(c_double), pointer :: rch_arr(:) => null()  !! BMI ptr to RECHARGE
      logical :: rch_ready = .false.
      character(len=256) :: rch_addr = "MODFLOW_SFR/RCHA_0/RECHARGE"
      real(c_double) :: rch_depth_cum = 0.0_c_double   !! diagnostic (sum m/day)

      !! ---- baseflow up-coupling (M3): MF6 SFR GWFLOW -> SWAT+ channels ----
      integer :: n_bf = 0           !! number of reach->channel links
      integer :: n_sfr = 0          !! number of SFR reaches (GWFLOW size)
      integer :: max_gis = 0        !! largest SWAT+ channel gis id in the map
      integer, allocatable :: bf_reach(:)  !! 0-based SFR reach index per link
      integer, allocatable :: bf_gis(:)    !! SWAT+ channel gis id per link
      real(c_double), pointer :: gwf_arr(:) => null()   !! BMI ptr to SFR GWFLOW
      real(c_double), allocatable :: bf_chan(:)  !! baseflow per gis id (m3/day)
      logical :: bf_ready = .false.
      character(len=256) :: gwf_addr = "MODFLOW_SFR/SFR_0/GWFLOW"
      real(c_double) :: bf_cum_gain = 0.0_c_double, bf_cum_loss = 0.0_c_double

      contains

      !! ================================================================
      subroutine mf6_coupler_init
      !! Detect mf6.con, dlopen libmf6.so (RTLD_LOCAL), resolve BMI entry
      !! points, initialize MF6 inside its workspace.
      use ifport, only : chdir
      implicit none
      integer :: iu, istat, ios
      logical :: have_con
      logical :: halt_save(5)
      character(len=256) :: cwd0
      character(len=600) :: tmp

      inquire (file="mf6.con", exist=have_con)
      if (.not. have_con) return        !! coupling not requested -> no-op
      call trace("init: mf6.con found")

      !! --- read control file ---
      open (newunit=iu, file="mf6.con", status="old", action="read")
      call read_setting(iu, mf6_ws)
      call read_int_setting(iu, gwf_cadence)
      call read_int_setting(iu, gwt_cadence)
      tmp = ""
      call read_setting(iu, tmp)
      if (len_trim(tmp) > 0) mf6_lib = tmp
      close (iu)
      if (gwf_cadence < 1) gwf_cadence = 1
      if (gwt_cadence < 1) gwt_cadence = gwf_cadence
      call trace("init: read mf6.con ws="//trim(mf6_ws)//" lib="//trim(mf6_lib))

      !! --- dlopen libmf6.so privately + resolve symbols ---
      mf6_handle = c_dlopen(trim(mf6_lib)//c_null_char, RTLD_NOW_LOCAL)
      if (.not. c_associated(mf6_handle)) then
        write (*,*) "MF6 COUPLER: dlopen failed for ", trim(mf6_lib)
        return
      end if
      call trace("init: dlopen OK")
      if (.not. bind_bmi()) then
        write (*,*) "MF6 COUPLER: could not resolve BMI symbols in libmf6.so"
        return
      end if
      call trace("init: BMI symbols resolved")

      !! --- initialize MF6 from inside its workspace (reads mfsim.nam in CWD) ---
      istat = getcwd_safe(cwd0)
      call trace("init: cwd0="//trim(cwd0))
      istat = chdir(trim(mf6_ws))
      call trace("init: chdir to workspace done")
      if (istat /= 0) then
        write (*,*) "MF6 COUPLER: cannot chdir to workspace ", trim(mf6_ws)
        return
      end if
      call trace("init: calling MF6 BMI initialize()...")
      call fpe_off(halt_save)        !! engine runs -fpe0; MF6 handles Inf/NaN itself
      ios = p_initialize()
      if (ios == 0) ios = ios + p_get_end_time(mf6_tend)
      call fpe_restore(halt_save)
      call trace("init: MF6 initialize returned")
      istat = chdir(trim(cwd0))

      if (ios /= 0) then
        write (*,*) "MF6 COUPLER: BMI initialize FAILED, status=", ios
        mf6_active = .false.
        return
      end if
      mf6_active = .true.
      mf6_daycount = 0
      write (*,'(a)')      " MF6 COUPLER: active (daily SWAT+ <-> MODFLOW 6)"
      write (*,'(a,a)')    "   workspace   : ", trim(mf6_ws)
      write (*,'(a,i0,a,i0,a)') "   cadence     : flow=", gwf_cadence, &
                                " d, transport=", gwt_cadence, " d"
      write (*,'(a,f10.1,a)')   "   MF6 end time: ", mf6_tend, " days"

      !! --- M2: set up the HRU -> MF6 recharge down-coupling (optional) ---
      call mf6_recharge_setup
      !! --- M3: set up the MF6 SFR baseflow -> SWAT+ channel up-coupling ---
      call mf6_baseflow_setup
      end subroutine mf6_coupler_init

      !! ================================================================
      subroutine mf6_recharge_setup
      !! Read mf6_recharge.map (area-weighted HRU->cell) and bind a pointer to
      !! the MF6 RECHARGE array, so each day we can overwrite recharge with the
      !! SWAT+ soil-profile percolation.  No-op (rch_ready stays .false.) if the
      !! map file is absent or the RECHARGE pointer can't be obtained -- in that
      !! case the coupler still steps MF6 on the model's own recharge.
      implicit none
      integer :: iu, ios, e, idx, hru
      integer :: nmap, nrch, ncolx, nhrux
      real :: w
      logical :: have_map
      type(c_ptr) :: cptr

      inquire (file="mf6_recharge.map", exist=have_map)
      if (.not. have_map) then
        call trace("recharge: no mf6_recharge.map -> MF6 uses its own recharge")
        return
      end if
      open (newunit=iu, file="mf6_recharge.map", status="old", action="read")
      read (iu,*,iostat=ios) nmap, nrch, ncolx, nhrux
      if (ios /= 0 .or. nmap <= 0) then; close(iu); return; end if
      n_map = nmap; n_rch = nrch
      allocate (map_idx(n_map), map_hru(n_map), map_w(n_map))
      do e = 1, n_map
        read (iu,*,iostat=ios) idx, hru, w
        if (ios /= 0) then
          call trace("recharge: map truncated -> disabled"); close(iu)
          deallocate (map_idx, map_hru, map_w); return
        end if
        map_idx(e) = idx; map_hru(e) = hru; map_w(e) = w
      end do
      close (iu)

      !! bind a Fortran pointer to MF6's RECHARGE array (rate, m/day)
      cptr = c_null_ptr
      ios = p_get_value_ptr_double(trim(rch_addr)//c_null_char, cptr)
      if (ios /= 0 .or. .not. c_associated(cptr)) then
        call trace("recharge: get_value_ptr_double failed for "//trim(rch_addr))
        deallocate (map_idx, map_hru, map_w); return
      end if
      call c_f_pointer(cptr, rch_arr, [n_rch])
      rch_ready = .true.
      write (*,'(a,i0,a,i0,a)') "   recharge map: ", n_map, &
        " HRU-cell links onto ", n_rch, " MF6 cells (SWAT+ drives recharge)"
      end subroutine mf6_recharge_setup

      !! ================================================================
      subroutine mf6_coupler_step
      !! Advance MF6 one coupling step in lockstep with the SWAT+ day.
      !! M1: step flow when due and confirm MF6 time advances.
      !! M2 will push recharge before the flow step; M3 will pull
      !! baseflow + PFAS after it and trigger transport at gwt_cadence.
      use ifport, only : chdir
      implicit none
      integer :: istat, ios
      logical :: halt_save(5)
      real(c_double) :: tnow, dt0
      character(len=256) :: cwd0

      if (.not. mf6_active) return
      mf6_daycount = mf6_daycount + 1
      if (mod(mf6_daycount, gwf_cadence) /= 0) return

      istat = getcwd_safe(cwd0)
      istat = chdir(trim(mf6_ws))
      call fpe_off(halt_save)
      ios = p_get_current_time(tnow)
      if (tnow < mf6_tend) then
        !! granular step so we can inject recharge between prepare and solve
        dt0 = 0.0_c_double
        ios = p_prepare_time_step(dt0)        !! rp: MF6 (re)loads its recharge
        if (rch_ready) call push_recharge     !! overwrite with SWAT+ percolation
        ios = p_do_time_step()
        if (ios /= 0) write (*,*) "MF6 COUPLER: do_time_step failed day ", &
                                  mf6_daycount, " status ", ios
        ios = p_finalize_time_step()
        if (bf_ready) call pull_baseflow      !! aggregate MF6 SFR baseflow -> channels
        ios = p_get_current_time(tnow)
      end if
      call fpe_restore(halt_save)
      istat = chdir(trim(cwd0))

      if (mod(mf6_daycount, 365) == 0) then
        write (*,'(a,i0,a,f9.1)') " MF6 COUPLER: engine day ", mf6_daycount, &
          "  MF6 t(d)=", tnow
        if (rch_ready) write (*,'(a,es11.3)') &
          "   cum recharge depth (m) = ", rch_depth_cum
        if (bf_ready) write (*,'(a,es11.3,a,es11.3,a)') &
          "   cum SFR exchange: gaining ", bf_cum_gain, "  losing ", bf_cum_loss, " (m3)"
      end if
      end subroutine mf6_coupler_step

      !! ================================================================
      subroutine push_recharge
      !! Overwrite MF6 RECHARGE (rate, m/day) with area-weighted SWAT+ soil-
      !! profile percolation sepbtm (mm/day):
      !!   RECHARGE[idx] = sum_HRU( sepbtm(hru)/1000 * weight )
      !! Called after prepare_time_step (which reloaded the model's own
      !! recharge) and before do_time_step, so SWAT+ fully drives recharge.
      use hru_module, only : sepbtm
      use hydrograph_module, only : sp_ob
      implicit none
      integer :: e, h, nhru
      real(c_double) :: dsum
      nhru = sp_ob%hru
      rch_arr(:) = 0.0_c_double
      do e = 1, n_map
        h = map_hru(e)
        if (h >= 1 .and. h <= nhru) &
          rch_arr(map_idx(e)+1) = rch_arr(map_idx(e)+1) &
            + real(sepbtm(h), c_double) / 1000.0_c_double * real(map_w(e), c_double)
      end do
      dsum = 0.0_c_double
      do e = 1, n_rch
        dsum = dsum + rch_arr(e)
      end do
      rch_depth_cum = rch_depth_cum + dsum
      end subroutine push_recharge

      !! ================================================================
      subroutine mf6_baseflow_setup
      !! Read mf6_baseflow.map (SFR reach -> SWAT+ channel gis id) and bind a
      !! pointer to the SFR GWFLOW array (per-reach groundwater<->stream
      !! exchange, m3/day).  No-op if the map is absent or the pointer fails.
      implicit none
      integer :: iu, ios, e, rea, gis, nbf, nsfr
      logical :: have_map
      type(c_ptr) :: cptr

      inquire (file="mf6_baseflow.map", exist=have_map)
      if (.not. have_map) then
        call trace("baseflow: no mf6_baseflow.map -> up-coupling off")
        return
      end if
      open (newunit=iu, file="mf6_baseflow.map", status="old", action="read")
      read (iu,*,iostat=ios) nbf, nsfr
      if (ios /= 0 .or. nbf <= 0) then; close(iu); return; end if
      n_bf = nbf; n_sfr = nsfr
      allocate (bf_reach(n_bf), bf_gis(n_bf))
      max_gis = 0
      do e = 1, n_bf
        read (iu,*,iostat=ios) rea, gis
        if (ios /= 0) then
          call trace("baseflow: map truncated -> disabled"); close(iu)
          deallocate (bf_reach, bf_gis); return
        end if
        bf_reach(e) = rea; bf_gis(e) = gis
        if (gis > max_gis) max_gis = gis
      end do
      close (iu)

      cptr = c_null_ptr
      ios = p_get_value_ptr_double(trim(gwf_addr)//c_null_char, cptr)
      if (ios /= 0 .or. .not. c_associated(cptr)) then
        call trace("baseflow: get_value_ptr_double failed for "//trim(gwf_addr))
        deallocate (bf_reach, bf_gis); return
      end if
      call c_f_pointer(cptr, gwf_arr, [n_sfr])
      allocate (bf_chan(0:max_gis)); bf_chan = 0.0_c_double
      bf_ready = .true.
      write (*,'(a,i0,a,i0,a)') "   baseflow map: ", n_bf, &
        " SFR reaches -> SWAT+ channels (gis 1..", max_gis, ")"
      end subroutine mf6_baseflow_setup

      !! ================================================================
      subroutine pull_baseflow
      !! Aggregate per-reach SFR<->aquifer exchange (GWFLOW, m3/day) onto
      !! SWAT+ channels by gis id.  bf_chan(gis) holds the day's groundwater
      !! discharge INTO that channel (m3/day, +ve = baseflow gain), ready for
      !! injection (M3b).  SIGN (verified vs the GWF budget): MF6 SFR GWFLOW>0 =
      !! reach LOSES to aquifer; GWFLOW<0 = aquifer feeds reach (baseflow).  So
      !! baseflow-into-channel = -GWFLOW.
      implicit none
      integer :: e, r
      real(c_double) :: bf
      bf_chan = 0.0_c_double
      do e = 1, n_bf
        r = bf_reach(e)
        if (r >= 0 .and. r < n_sfr) then
          bf = -gwf_arr(r+1)                  !! aquifer -> stream, +ve = gain
          bf_chan(bf_gis(e)) = bf_chan(bf_gis(e)) + bf
          if (bf > 0.0_c_double) then
            bf_cum_gain = bf_cum_gain + bf    !! baseflow (aquifer -> stream)
          else
            bf_cum_loss = bf_cum_loss + bf    !! seepage (stream -> aquifer)
          end if
        end if
      end do
      end subroutine pull_baseflow

      !! ================================================================
      subroutine mf6_coupler_finalize
      use ifport, only : chdir
      implicit none
      integer :: istat, ios
      logical :: halt_save(5)
      character(len=256) :: cwd0
      if (.not. mf6_active) return
      istat = getcwd_safe(cwd0)
      istat = chdir(trim(mf6_ws))
      call fpe_off(halt_save)
      ios = p_finalize()
      call fpe_restore(halt_save)
      istat = chdir(trim(cwd0))
      if (c_associated(mf6_handle)) istat = c_dlclose(mf6_handle)
      write (*,'(a,i0)') " MF6 COUPLER: finalized, status=", ios
      mf6_active = .false.
      end subroutine mf6_coupler_finalize

      !! ================================================================
      logical function mf6_baseflow_active()
      !! true when MF6 supplies the channel baseflow (so SWAT+ should skip its
      !! own aquifer->channel return flow to avoid double counting)
      mf6_baseflow_active = mf6_active .and. bf_ready
      end function mf6_baseflow_active

      real(c_double) function mf6_channel_baseflow(gis)
      !! groundwater discharge into the SWAT+ channel with this gis id
      !! (m3/day; +ve = aquifer feeds stream, -ve = stream seeps to aquifer).
      !! From the previous day's MF6 solve (one-day explicit lag).
      integer, intent(in) :: gis
      mf6_channel_baseflow = 0.0_c_double
      if (bf_ready .and. gis >= 0 .and. gis <= max_gis) &
        mf6_channel_baseflow = bf_chan(gis)
      end function mf6_channel_baseflow

      !! ---- helpers --------------------------------------------------
      subroutine fpe_off(saved)
      !! The engine is built with -fpe0 (trap on invalid/divzero/overflow),
      !! a PROCESS-WIDE FPU setting.  MODFLOW 6 legitimately produces and
      !! handles Inf/NaN (dry cells, etc.); under -fpe0 those would SIGFPE.
      !! Disable FP halting around MF6 BMI calls, saving the engine's modes.
      use ieee_exceptions, only : ieee_get_halting_mode, ieee_set_halting_mode, ieee_all
      logical, intent(out) :: saved(:)
      call ieee_get_halting_mode(ieee_all, saved)
      call ieee_set_halting_mode(ieee_all, .false.)
      end subroutine fpe_off

      subroutine fpe_restore(saved)
      use ieee_exceptions, only : ieee_set_halting_mode, ieee_all
      logical, intent(in) :: saved(:)
      call ieee_set_halting_mode(ieee_all, saved)
      end subroutine fpe_restore

      subroutine trace(msg)
      !! unbuffered diagnostic to stderr (survives a hard crash)
      use iso_fortran_env, only : error_unit
      character(len=*), intent(in) :: msg
      write (error_unit,'(a)') " [mf6_coupler] "//msg
      flush (error_unit)
      end subroutine trace

      logical function bind_bmi()
      !! resolve the five BMI entry points; .false. if any is missing
      bind_bmi = .false.
      if (.not. resolve("initialize",       p_initialize))       return
      if (.not. resolve("update",           p_update))           return
      if (.not. resolve("finalize",         p_finalize))         return
      if (.not. resolve_dbl("get_current_time", p_get_current_time)) return
      if (.not. resolve_dbl("get_end_time",     p_get_end_time))     return
      if (.not. resolve("do_time_step",       p_do_time_step))       return
      if (.not. resolve("finalize_time_step", p_finalize_time_step)) return
      if (.not. resolve_dt("prepare_time_step", p_prepare_time_step)) return
      if (.not. resolve_ptr("get_value_ptr_double", p_get_value_ptr_double)) return
      bind_bmi = .true.
      end function bind_bmi

      logical function resolve_dt(name, ptr)
      character(len=*), intent(in) :: name
      procedure(bmi_dt), pointer, intent(out) :: ptr
      type(c_funptr) :: fp
      fp = c_dlsym(mf6_handle, trim(name)//c_null_char)
      resolve_dt = c_associated(fp)
      if (resolve_dt) call c_f_procpointer(fp, ptr)
      end function resolve_dt

      logical function resolve_ptr(name, ptr)
      character(len=*), intent(in) :: name
      procedure(bmi_getptr), pointer, intent(out) :: ptr
      type(c_funptr) :: fp
      fp = c_dlsym(mf6_handle, trim(name)//c_null_char)
      resolve_ptr = c_associated(fp)
      if (resolve_ptr) call c_f_procpointer(fp, ptr)
      end function resolve_ptr

      logical function resolve(name, ptr)
      character(len=*), intent(in) :: name
      procedure(bmi_noarg), pointer, intent(out) :: ptr
      type(c_funptr) :: fp
      fp = c_dlsym(mf6_handle, trim(name)//c_null_char)
      resolve = c_associated(fp)
      if (resolve) call c_f_procpointer(fp, ptr)
      end function resolve

      logical function resolve_dbl(name, ptr)
      character(len=*), intent(in) :: name
      procedure(bmi_dbl), pointer, intent(out) :: ptr
      type(c_funptr) :: fp
      fp = c_dlsym(mf6_handle, trim(name)//c_null_char)
      resolve_dbl = c_associated(fp)
      if (resolve_dbl) call c_f_procpointer(fp, ptr)
      end function resolve_dbl

      integer function getcwd_safe(path)
      use ifport, only : getcwd
      character(len=*), intent(out) :: path
      getcwd_safe = getcwd(path)
      end function getcwd_safe

      subroutine read_setting(iu, val)
      !! read next non-comment, non-blank token into val.  NB: we extract the
      !! first whitespace-delimited token MANUALLY -- a list-directed read
      !! (read(line,*) val) truncates a path at the first "/" because "/"
      !! terminates list-directed input, which silently turned "./mf6" -> ".".
      integer, intent(in) :: iu
      character(len=*), intent(out) :: val
      character(len=600) :: line
      integer :: ios, ix
      val = ""
      do
        read (iu,'(a)',iostat=ios) line
        if (ios /= 0) return
        line = adjustl(line)
        if (len_trim(line) == 0) cycle
        if (line(1:1) == "!") cycle
        ix = index(line, "!")               ! strip inline comment
        if (ix > 0) line = line(:ix-1)
        ix = scan(trim(line), " "//char(9))  ! first space or tab
        if (ix > 0) then
          val = line(:ix-1)
        else
          val = trim(line)
        end if
        return
      end do
      end subroutine read_setting

      subroutine read_int_setting(iu, ival)
      integer, intent(in) :: iu
      integer, intent(inout) :: ival
      character(len=256) :: tok
      integer :: ios, tmp
      tok = ""
      call read_setting(iu, tok)
      if (len_trim(tok) == 0) return
      read (tok,*,iostat=ios) tmp
      if (ios == 0) ival = tmp
      end subroutine read_int_setting

      end module mf6_coupler
