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
      public :: mf6_pfas_disch_active, mf6_channel_pfas, mf6_pfas_compound
      public :: mf6_withhold_aquifer, mf6_note_unmet

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
      procedure(bmi_getptr),pointer :: p_get_value_ptr_int => null()

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
      real :: rch_mult = 1.0        !! recharge multiplier (mf6.con line 4; OAT knob)

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
      !! Q1130 M1: SWAT+ aquifer return NOT added to channels while coupled (MF6
      !! carries that groundwater), accumulated per receiving object (m3) so each
      !! channel writes only its own element when channels run in parallel.
      real(c_double), allocatable :: aq_withheld_ob(:)
      !! Q1130 M2: the channel floor (sd_channel_control3) discards the part of an MF6
      !! loss that the SWAT+ channel could not supply; MF6 had already moved that water
      !! into the aquifer, so without a return it is created. It is recorded per channel
      !! gis and taken back from MF6 RECHARGE at the next MF6 step, in the columns of that
      !! channel's losing reaches, in proportion to their loss.
      real(c_double), allocatable :: unmet_chan(:)      !! (0:max_gis) unmet m3 since last step
      real(c_double), allocatable :: gwf_prev(:)        !! GWFLOW of the last pulled step
      real(c_double) :: cb_carried = 0.0_c_double, cb_unplaced = 0.0_c_double
      integer :: cb_negrch = 0                          !! column-steps with RECHARGE < 0
      integer(c_int), pointer :: cb_igwf(:) => null(), cb_nodeuser(:) => null()
      real(c_double), pointer :: cb_delr(:) => null(), cb_delc(:) => null()
      integer :: cb_nrow = 0, cb_ncol = 0
      logical :: cb_reduced = .false., cb_ready = .false.
      integer :: cb_iu = 0                              !! mf6_floor.out (day, gis, m3); NEWUNIT gives a NEGATIVE unit, 0 = not open
      !! Q1144: one coupled day may take several MF6 time steps (an ATS cut, NSTP > 1). The day's RECHARGE and SRC rates
      !! are kept from its first sub-step and re-applied to every later one; mf6_clock.out records each coupled day.
      real(c_double), allocatable :: rch_day(:), src_day(:)
      integer :: clk_iu = 0                             !! mf6_clock.out (day, t_before, t_after, nsub); NEWUNIT is NEGATIVE, 0 = not open
      integer :: clk_days = 0, clk_subdays = 0, clk_maxsub = 0

      !! ---- Q1137: SWAT+ drives SFR's inflows, so SFR and SWAT+ route the same water ----
      !! On exactly when mf6_drive.map is present. The deck is then cut at SWAT+ channel boundaries (Q1157), and the map names, for
      !! every connection the cut removed into a mapped reach: the entry reach, its source channel, a static length weight and the
      !! source channel's exit reach. Each coupled day every mapped reach takes its water from SWAT+ only: INFLOW at a channel's entry
      !! reaches = SWAT+'s inflow from each upstream channel; RUNOFF on all of its reaches = its other (lateral) inflow, split by reach
      !! LENGTH. SFR adds the day's groundwater gains itself (ob%hin_d excludes the channel's own MF6 exchange). Stock MF6 books each
      !! cut channel outlet as EXT-OUTFLOW, so SFR's budget and SFT stay closed. Design: q1137/DESIGN.txt Addenda 1-3.
      logical :: drv_ready = .false.
      real(c_double), pointer :: drv_inflow(:) => null(), drv_runoff(:) => null()
      real(c_double), pointer :: drv_length(:) => null(), drv_rain(:) => null()
      integer(c_int), pointer :: drv_ia(:) => null(), drv_ja(:) => null(), drv_idir(:) => null()
      integer, allocatable :: drv_key(:)                 !! (1:n_sfr) the reach's channel gis; -1 = unmapped (its own network)
      integer, allocatable :: drv_exit(:)                !! (1:n_sfr) the reach by which its water leaves its network (1-based)
      integer, allocatable :: drv_ptr(:), drv_reach(:)   !! each gis's reaches (1-based), CSR over gis 0..max_gis
      real(c_double), allocatable :: drv_wlen(:)         !! those reaches' LENGTH shares of their channel
      integer, allocatable :: gis_ob(:)                  !! (0:max_gis) SWAT+ object of the channel with this gis, 0 = none
      integer, allocatable :: drv_order(:)               !! driven channels, each two-exit source before its receiver where possible
      integer :: n_order = 0, n_ent = 0
      integer, allocatable :: ent_ptr(:)                 !! map entries grouped by receiving gis, CSR over gis 0..max_gis
      integer, allocatable :: ent_reach(:), ent_src(:), ent_exit(:)  !! entry reach (1-based), source gis, source exit reach (1-based)
      real(c_double), allocatable :: ent_w(:)            !! the map's static length weight within its (source, receiver) pair
      logical, allocatable :: ent_dyn(:)                 !! a two-exit entry whose source is ordered first: split by the day's water
      real(c_double), allocatable :: drv_s(:)            !! (1:n_sfr) the day's injection that drains to each exit reach
      real(c_double), allocatable :: inflow_day(:), runoff_day(:)
      integer :: drv_iu = 0                              !! mf6_drive.out (day, gis, volumes); NEWUNIT is NEGATIVE, 0 = not open
      real(c_double) :: drv_cum_swat = 0.0_c_double, drv_cum_in = 0.0_c_double, drv_cum_ro = 0.0_c_double
      real(c_double) :: drv_cum_tolat = 0.0_c_double, drv_cum_aqu = 0.0_c_double, drv_cum_moved = 0.0_c_double
      integer :: drv_raindays = 0, drv_fallback = 0
      logical :: drv_have_map = .false.
      character(len=200) :: drv_mode = ""               !! the run's mode line: printed once, and first in mf6_drive.out

      !! ---- PFAS down-coupling (M4b): SWAT+ leaching -> GWT SRC mass source --
      integer :: n_leach = 0        !! number of HRU->SRC-cell leaching links
      integer :: n_src = 0          !! number of SRC cells (BOUND size)
      integer :: pfas_k = 1         !! which simulated PFAS compound (PFOS) feeds GWT
      integer, allocatable :: leach_src(:)  !! 0-based SRC cell index per link
      integer, allocatable :: leach_hru(:)  !! SWAT+ HRU id per link
      real,    allocatable :: leach_ov(:)   !! overlap area (ha) per link
      real(c_double), pointer :: src_arr(:) => null()   !! BMI ptr to SRC BOUND
      logical :: src_ready = .false.
      character(len=256) :: src_addr = "PFAS/SRC_PFAS/BOUND"
      real(c_double) :: pfas_src_cum = 0.0_c_double   !! cumulative ug loaded

      !! ---- PFAS up-coupling (M4c): GW-discharged PFAS -> SWAT+ channels ----
      real(c_double), pointer :: gwt_conc(:) => null()  !! BMI ptr GWT concentration (ng/L)
      integer(c_int), pointer :: igwf(:) => null()      !! BMI ptr SFR reach -> GWF node
      real(c_double), allocatable :: pfas_chan(:)       !! GW PFAS to channel (kg/day) by gis
      logical :: pdis_ready = .false.
      character(len=256) :: gwtx_addr = "PFAS/X"
      character(len=256) :: igwf_addr = "MODFLOW_SFR/SFR_0/IGWFNODE"
      real(c_double) :: pfas_disch_cum = 0.0_c_double   !! cumulative kg discharged

      contains

      !! ================================================================
      subroutine mf6_coupler_init
      !! Detect mf6.con, dlopen libmf6.so (RTLD_LOCAL), resolve BMI entry
      !! points, initialize MF6 inside its workspace.
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
      use ifport, only : chdir          !! gfortran provides CHDIR as a GNU intrinsic
#endif
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
      !! line 4 (optional): recharge multiplier (calibration / OAT sensitivity knob)
      tmp = ""
      call read_setting(iu, tmp)
      if (len_trim(tmp) > 0) then
        read (tmp,*,iostat=ios) rch_mult
        if (ios /= 0) rch_mult = 1.0
      end if
      !! line 5 (optional): libmf6.so path
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
        call coupler_fatal("dlopen failed for "//trim(mf6_lib))
      end if
      call trace("init: dlopen OK")
      if (.not. bind_bmi()) then
        call coupler_fatal("could not resolve BMI symbols in "//trim(mf6_lib))
      end if
      call trace("init: BMI symbols resolved")

      !! --- initialize MF6 from inside its workspace (reads mfsim.nam in CWD) ---
      istat = getcwd_safe(cwd0)
      call trace("init: cwd0="//trim(cwd0))
      istat = chdir(trim(mf6_ws))
      call trace("init: chdir to workspace done")
      if (istat /= 0) then
        call coupler_fatal("cannot chdir to workspace "//trim(mf6_ws))
      end if
      call trace("init: calling MF6 BMI initialize()...")
      call fpe_off(halt_save)        !! engine runs -fpe0; MF6 handles Inf/NaN itself
      ios = p_initialize()
      if (ios == 0) ios = ios + p_get_end_time(mf6_tend)
      call fpe_restore(halt_save)
      call trace("init: MF6 initialize returned")
      istat = chdir(trim(cwd0))

      if (ios /= 0) then
        mf6_active = .false.
        call coupler_fatal("BMI initialize failed (see the MF6 listing in the workspace)")
      end if
      mf6_active = .true.
      mf6_daycount = 0
      write (*,'(a)')      " MF6 COUPLER: active (daily SWAT+ <-> MODFLOW 6)"
      !! Q1137: the run states its SFR drive mode itself, as the coupler's first line (the switch is a file's presence)
      call drive_mode_line
      write (*,'(a,a)')    "   workspace   : ", trim(mf6_ws)
      write (*,'(a,i0,a,i0,a)') "   cadence     : flow=", gwf_cadence, &
                                " d, transport=", gwt_cadence, " d"
      write (*,'(a,f10.1,a)')   "   MF6 end time: ", mf6_tend, " days"

      !! --- M2: set up the HRU -> MF6 recharge down-coupling (optional) ---
      call mf6_recharge_setup
      !! --- M3: set up the MF6 SFR baseflow -> SWAT+ channel up-coupling ---
      call mf6_baseflow_setup
      !! --- M4b: set up the SWAT+ PFAS leaching -> GWT SRC mass source ---
      call mf6_pfas_src_setup
      !! --- M4c: set up the GW PFAS discharge -> SWAT+ channel up-coupling ---
      call mf6_pfas_disch_setup
      !! --- Q1137: SWAT+ drives SFR's inflows (mf6_drive.map), or check that the deck was not cut without its map ---
      call mf6_drive_setup
      !! --- Q1144: the per-day MF6 clock record ---
      open (newunit=clk_iu, file="mf6_clock.out", status="replace", action="write")
      write (clk_iu,'(a)') "day t_before t_after nsub   (Q1144: MF6 time at the start and end of each coupled day; MF6 time steps it took)"
      flush (clk_iu)
      end subroutine mf6_coupler_init

      !! ================================================================
      subroutine mf6_pfas_disch_setup
      !! Bind pointers to the GWT concentration (PFAS/X) and the SFR reach->GWF
      !! node map (IGWFNODE).  Reuses the M3 baseflow reach->channel map, so it
      !! needs bf_ready.  Each gaining reach discharges PFAS to its channel at
      !! rate baseflow * GW concentration.  No-op if any pointer is unavailable.
      implicit none
      integer :: ios
      type(c_ptr) :: cx, ci
      if (.not. bf_ready) then
        call trace("pfas-disch: baseflow map not ready -> PFAS up off"); return
      end if
      cx = c_null_ptr; ci = c_null_ptr
      ios = p_get_value_ptr_double(trim(gwtx_addr)//c_null_char, cx)
      if (ios /= 0 .or. .not. c_associated(cx)) then
        call trace("pfas-disch: no GWT concentration ("//trim(gwtx_addr)//")"); return
      end if
      ios = p_get_value_ptr_int(trim(igwf_addr)//c_null_char, ci)
      if (ios /= 0 .or. .not. c_associated(ci)) then
        call trace("pfas-disch: no IGWFNODE ("//trim(igwf_addr)//")"); return
      end if
      call c_f_pointer(ci, igwf, [n_sfr])
      !! GWT concentration array length is unknown a priori; use the largest
      !! node referenced by the reaches as a safe lower bound and bind generously.
      call c_f_pointer(cx, gwt_conc, [maxval(igwf)])
      allocate (pfas_chan(0:max_gis)); pfas_chan = 0.0_c_double
      pdis_ready = .true.
      write (*,'(a)') "   pfas discharge: GW concentration -> channels (PFAS up-coupling on)"
      end subroutine mf6_pfas_disch_setup

      !! ================================================================
      subroutine mf6_pfas_src_setup
      !! Read pfas_leach.map (HRU -> SRC cell, overlap ha) and bind a pointer to
      !! the GWT SRC mass-loading array (PFAS/SRC_PFAS/BOUND).  No-op if absent.
      implicit none
      integer :: iu, ios, e, sidx, hru, nleach, nsrc
      real :: ov
      logical :: have_map
      type(c_ptr) :: cptr

      inquire (file="pfas_leach.map", exist=have_map)
      if (.not. have_map) then
        call trace("pfas-src: no pfas_leach.map -> PFAS down-coupling off")
        return
      end if
      open (newunit=iu, file="pfas_leach.map", status="old", action="read")
      read (iu,*,iostat=ios) nleach, nsrc
      if (ios /= 0 .or. nleach <= 0) call coupler_fatal("pfas_leach.map header unreadable or count <= 0")
      n_leach = nleach; n_src = nsrc
      allocate (leach_src(n_leach), leach_hru(n_leach), leach_ov(n_leach))
      do e = 1, n_leach
        read (iu,*,iostat=ios) sidx, hru, ov
        if (ios /= 0) then
          call coupler_fatal("pfas_leach.map truncated (fewer rows than its header count)")
        end if
        leach_src(e) = sidx; leach_hru(e) = hru; leach_ov(e) = ov
      end do
      close (iu)

      cptr = c_null_ptr
      ios = p_get_value_ptr_double(trim(src_addr)//c_null_char, cptr)
      if (ios /= 0 .or. .not. c_associated(cptr)) then
        call coupler_fatal("cannot bind MF6 variable "//trim(src_addr))
      end if
      call c_f_pointer(cptr, src_arr, [n_src])
      src_ready = .true.
      write (*,'(a,i0,a,i0,a)') "   pfas leach map: ", n_leach, &
        " HRU-cell links onto ", n_src, " GWT source cells (SWAT+ drives PFAS)"
      end subroutine mf6_pfas_src_setup

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
      if (ios /= 0 .or. nmap <= 0) call coupler_fatal("mf6_recharge.map header unreadable or count <= 0")
      n_map = nmap; n_rch = nrch
      allocate (map_idx(n_map), map_hru(n_map), map_w(n_map))
      do e = 1, n_map
        read (iu,*,iostat=ios) idx, hru, w
        if (ios /= 0) then
          call coupler_fatal("mf6_recharge.map truncated (fewer rows than its header count)")
        end if
        map_idx(e) = idx; map_hru(e) = hru; map_w(e) = w
      end do
      close (iu)

      !! bind a Fortran pointer to MF6's RECHARGE array (rate, m/day)
      cptr = c_null_ptr
      ios = p_get_value_ptr_double(trim(rch_addr)//c_null_char, cptr)
      if (ios /= 0 .or. .not. c_associated(cptr)) then
        call coupler_fatal("cannot bind MF6 variable "//trim(rch_addr))
      end if
      call c_f_pointer(cptr, rch_arr, [n_rch])
      rch_ready = .true.
      write (*,'(a,i0,a,i0,a)') "   recharge map: ", n_map, &
        " HRU-cell links onto ", n_rch, " MF6 cells (SWAT+ drives recharge)"
      end subroutine mf6_recharge_setup

      !! ================================================================
      subroutine mf6_coupler_step
      !! Advance MF6 through one coupled SWAT+ day, in as many MF6 time steps as MF6's own clock takes to reach
      !! the day end (Q1144): push recharge (M2) + PFAS sources before the first step, hold them for the day,
      !! and pull baseflow + PFAS (M3) after each step, weighted by its share of the day.
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
      use ifport, only : chdir          !! gfortran provides CHDIR as a GNU intrinsic
#endif
      implicit none
      integer :: istat, ios, nsub, e
      logical :: halt_save(5)
      real(c_double) :: tnow, dt0, tprev, tday, tsub, w
      character(len=256) :: cwd0

      if (.not. mf6_active) return
      mf6_daycount = mf6_daycount + 1
      if (mod(mf6_daycount, gwf_cadence) /= 0) return

      istat = getcwd_safe(cwd0)
      istat = chdir(trim(mf6_ws))
      call fpe_off(halt_save)
      ios = p_get_current_time(tnow)
      if (tnow < mf6_tend) then
        tprev = tnow                          !! MF6 time at the start of this coupled day
        !! Q1144: the coupled day ends here (or where MF6 ends, if sooner). MF6 is stepped until its OWN clock reaches it:
        !! an ATS cut or NSTP > 1 gives the day several MF6 time steps. One step per day left MF6 up to 26 periods behind
        !! SWAT+ with rc 0 (Q1040(b) OAT, Q1130 R7b-1); Q1040(c) stopped such a run at its first cut instead.
        tday = min(tprev + real(gwf_cadence, c_double), mf6_tend)
        nsub = 0
        do
          nsub = nsub + 1
          tsub = tnow
          !! granular step so we can inject recharge between prepare and solve
          dt0 = 0.0_c_double
          ios = p_prepare_time_step(dt0)        !! rp: MF6 (re)loads recharge + sources
          if (nsub == 1) then
            if (rch_ready) call push_recharge     !! overwrite recharge with SWAT+ percolation
            if (bf_ready) call carry_back_unmet   !! Q1130 M2: take back what the channel floor discarded
            if (src_ready) call push_pfas_src     !! overwrite GWT source with SWAT+ PFAS leaching
            if (drv_ready) call push_drive(.true.)  !! Q1137: SWAT+'s inflows of this day become SFR's INFLOW and RUNOFF
            !! Q1144: keep the day's rates, so each later sub-step of this day carries the same ones
            if (rch_ready) then
              if (.not. allocated(rch_day)) allocate (rch_day(n_rch))
              rch_day(:) = rch_arr(1:n_rch)
            end if
            if (src_ready) then
              if (.not. allocated(src_day)) allocate (src_day(n_src))
              src_day(:) = src_arr(1:n_src)
            end if
          else
            !! Q1144: MF6 may have reloaded its own recharge/sources at a period start; the day's rates go back in
            if (rch_ready) rch_arr(1:n_rch) = rch_day(:)
            if (src_ready) src_arr(1:n_src) = src_day(:)
            if (drv_ready) call push_drive(.false.)
          end if
          ios = p_do_time_step()
          if (ios /= 0) then
            write (*,*) "MF6 COUPLER: do_time_step failed day ", mf6_daycount, " status ", ios
            call coupler_fatal("MF6 do_time_step failed on a coupled day")
          end if
          ios = p_finalize_time_step()
          if (ios /= 0) then
            !! Q1040(c): MF6 reports non-convergence here, not at do_time_step. Discarding it let a diverged
            !! run print a full summary with rc 0 (Q1040(b) OAT: 6 of 13 gfortran cases).
            write (*,*) "MF6 COUPLER: finalize_time_step failed day ", mf6_daycount, " status ", ios
            call coupler_fatal("MF6 finalize_time_step failed on a coupled day (MF6 did not converge; see mfsim.lst)")
          end if
          ios = p_get_current_time(tnow)
          if (.not. (tnow > tsub)) call clock_stop(tprev, tnow, tday, nsub, &
            "MF6 clock did not advance within a coupled day")
          if (tnow > tday + 1.0e-6_c_double) call clock_stop(tprev, tnow, tday, nsub, &
            "MF6 clock overshot the SWAT+ day end (an MF6 time step is longer than the coupling step)")
          if (nsub == 1 .and. tnow >= tday - 1.0e-6_c_double) then
            !! the whole day in ONE MF6 step: the pulls are today's, unchanged
            if (bf_ready) call pull_baseflow      !! aggregate MF6 SFR baseflow -> channels
            if (pdis_ready) call pull_pfas_discharge  !! GW PFAS discharge -> channels
            exit
          end if
          !! Q1144: a sub-step; its pulls count with its share of the day
          w = (tnow - tsub) / (tday - tprev)
          if (bf_ready) call pull_baseflow_part(w, nsub == 1)
          if (pdis_ready) call pull_pfas_discharge_part(w, nsub == 1)
          if (tnow >= tday - 1.0e-6_c_double) then
            if (pdis_ready) then
              do e = 0, max_gis
                pfas_disch_cum = pfas_disch_cum + pfas_chan(e)
              end do
            end if
            exit
          end if
          if (nsub >= 10000) call clock_stop(tprev, tnow, tday, nsub, &
            "MF6 clock cannot reach the SWAT+ day end (more than 10000 sub-steps)")
        end do
        clk_days = clk_days + 1
        if (nsub > 1) clk_subdays = clk_subdays + 1
        clk_maxsub = max(clk_maxsub, nsub)
        if (clk_iu /= 0) then
          write (clk_iu,'(i0,1x,f18.8,1x,f18.8,1x,i0)') mf6_daycount, tprev, tnow, nsub
          flush (clk_iu)
        end if
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
        if (src_ready) write (*,'(a,es11.3,a)') &
          "   cum PFAS loaded to aquifer = ", pfas_src_cum, " (ug)"
        if (pdis_ready) write (*,'(a,es11.3,a)') &
          "   cum PFAS discharged to streams = ", pfas_disch_cum, " (kg)"
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
            + real(sepbtm(h), c_double) / 1000.0_c_double * real(map_w(e), c_double) &
              * real(rch_mult, c_double)
      end do
      dsum = 0.0_c_double
      do e = 1, n_rch
        dsum = dsum + rch_arr(e)
      end do
      rch_depth_cum = rch_depth_cum + dsum
      end subroutine push_recharge

      !! ================================================================
      subroutine push_pfas_src
      !! Overwrite the GWT SRC mass-loading (ug/day) with SWAT+ soil-profile
      !! PFAS leaching: SRC[cell] = sum_HRU( hpfasb_d(hru)%perc(pfas_k)[kg/ha]
      !! * overlap_ha ) * 1e9.  (GWT conc = ng/L, length = m => GWT mass unit =
      !! ng/L*m^3 = ug; 1 kg = 1e9 ug.)  Called between prepare & do_time_step.
      use pfas_output_module, only : hpfasb_d
      use pfas_module, only : npfas
      use hydrograph_module, only : sp_ob
      implicit none
      integer :: e, h, nhru
      real(c_double) :: dsum, mass
      nhru = sp_ob%hru
      src_arr(:) = 0.0_c_double
      if (pfas_k < 1 .or. pfas_k > npfas) return
      do e = 1, n_leach
        h = leach_hru(e)
        if (h >= 1 .and. h <= nhru) then
          if (allocated(hpfasb_d(h)%perc)) then
            mass = real(hpfasb_d(h)%perc(pfas_k), c_double) &   ! kg/ha
                 * real(leach_ov(e), c_double) * 1.0e9_c_double  ! * ha * 1e9 = ug/day
            src_arr(leach_src(e)+1) = src_arr(leach_src(e)+1) + mass
          end if
        end if
      end do
      dsum = 0.0_c_double
      do e = 1, n_src
        dsum = dsum + src_arr(e)
      end do
      pfas_src_cum = pfas_src_cum + dsum
      end subroutine push_pfas_src

      !! ================================================================
      subroutine pull_pfas_discharge
      !! GW PFAS mass discharged to each SWAT+ channel (kg/day):
      !!   pfas_chan(gis) = sum_reach( baseflow[m3/d] * GWconc[ng/L] ) * 1e-9
      !! over gaining reaches (baseflow = -GWFLOW > 0).  ng/L * m3 * 1000 L/m3 *
      !! 1e-12 kg/ng = 1e-9.  GW concentration read at the reach's GWF node.
      implicit none
      integer :: e, r, node
      real(c_double) :: bf, conc
      pfas_chan = 0.0_c_double
      do e = 1, n_bf
        r = bf_reach(e)
        if (r >= 0 .and. r < n_sfr) then
          bf = -gwf_arr(r+1)                       ! baseflow, +ve = aquifer->stream
          if (bf > 0.0_c_double) then
            node = igwf(r+1)                        ! reach's GWF/GWT node (reduced)
            if (node >= 1 .and. node <= size(gwt_conc)) then
              conc = gwt_conc(node)                 ! ng/L
              if (conc > 0.0_c_double) then
                pfas_chan(bf_gis(e)) = pfas_chan(bf_gis(e)) &
                  + bf * conc * 1.0e-9_c_double      ! kg/day
              end if
            end if
          end if
        end if
      end do
      do e = 0, max_gis
        pfas_disch_cum = pfas_disch_cum + pfas_chan(e)
      end do
      end subroutine pull_pfas_discharge

      !! ================================================================
      subroutine pull_pfas_discharge_part(w, first)
      !! Q1144: pull_pfas_discharge for ONE MF6 sub-step of a coupled day that took several. Each sub-step adds its
      !! discharge with its share w of the day, so pfas_chan is the day's average (kg/day). The day's total goes into
      !! pfas_disch_cum once, when the day completes (mf6_coupler_step). A one-step day uses pull_pfas_discharge.
      implicit none
      real(c_double), intent(in) :: w
      logical, intent(in) :: first
      integer :: e, r, node
      real(c_double) :: bf, conc
      if (first) pfas_chan = 0.0_c_double
      do e = 1, n_bf
        r = bf_reach(e)
        if (r >= 0 .and. r < n_sfr) then
          bf = -gwf_arr(r+1)                       ! baseflow, +ve = aquifer->stream
          if (bf > 0.0_c_double) then
            node = igwf(r+1)                        ! reach's GWF/GWT node (reduced)
            if (node >= 1 .and. node <= size(gwt_conc)) then
              conc = gwt_conc(node)                 ! ng/L
              if (conc > 0.0_c_double) then
                pfas_chan(bf_gis(e)) = pfas_chan(bf_gis(e)) &
                  + bf * conc * 1.0e-9_c_double * w  ! kg/day, this sub-step's share
              end if
            end if
          end if
        end if
      end do
      end subroutine pull_pfas_discharge_part

      !! ================================================================
      subroutine mf6_baseflow_setup
      !! Read mf6_baseflow.map (SFR reach -> SWAT+ channel gis id) and bind a
      !! pointer to the SFR GWFLOW array (per-reach groundwater<->stream
      !! exchange, m3/day).  No-op if the map is absent or the pointer fails.
      use hydrograph_module, only : sp_ob
      implicit none
      integer :: iu, ios, e, rea, gis, nbf, nsfr, nhru_cov
      logical :: have_map
      logical, allocatable :: hru_seen(:)
      type(c_ptr) :: cptr

      inquire (file="mf6_baseflow.map", exist=have_map)
      if (.not. have_map) then
        call trace("baseflow: no mf6_baseflow.map -> up-coupling off")
        return
      end if
      open (newunit=iu, file="mf6_baseflow.map", status="old", action="read")
      read (iu,*,iostat=ios) nbf, nsfr
      if (ios /= 0 .or. nbf <= 0) call coupler_fatal("mf6_baseflow.map header unreadable or count <= 0")
      n_bf = nbf; n_sfr = nsfr
      allocate (bf_reach(n_bf), bf_gis(n_bf))
      max_gis = 0
      do e = 1, n_bf
        read (iu,*,iostat=ios) rea, gis
        if (ios /= 0) then
          call coupler_fatal("mf6_baseflow.map truncated (fewer rows than its header count)")
        end if
        bf_reach(e) = rea; bf_gis(e) = gis
        if (gis > max_gis) max_gis = gis
      end do
      close (iu)

      cptr = c_null_ptr
      ios = p_get_value_ptr_double(trim(gwf_addr)//c_null_char, cptr)
      if (ios /= 0 .or. .not. c_associated(cptr)) then
        call coupler_fatal("cannot bind MF6 variable "//trim(gwf_addr))
      end if
      call c_f_pointer(cptr, gwf_arr, [n_sfr])
      allocate (bf_chan(0:max_gis)); bf_chan = 0.0_c_double
      allocate (aq_withheld_ob(sp_ob%objs)); aq_withheld_ob = 0.0_c_double
      allocate (unmet_chan(0:max_gis)); unmet_chan = 0.0_c_double
      allocate (gwf_prev(n_sfr)); gwf_prev = 0.0_c_double
      bf_ready = .true.
      call mf6_carryback_setup
      write (*,'(a,i0,a,i0,a)') "   baseflow map: ", n_bf, &
        " SFR reaches -> SWAT+ channels (gis 1..", max_gis, ")"
      !! Q1130 M1 replaces the SWAT+ aquifer return for the WHOLE model, which is right
      !! only when MF6 receives the recharge of every HRU; print the coverage so a
      !! partial MF6 domain is visible.
      if (rch_ready) then
        nhru_cov = 0
        allocate (hru_seen(sp_ob%hru)); hru_seen = .false.
        do e = 1, n_map
          if (map_hru(e) >= 1 .and. map_hru(e) <= sp_ob%hru) hru_seen(map_hru(e)) = .true.
        end do
        nhru_cov = count(hru_seen)
        deallocate (hru_seen)
        write (*,'(a,i0,a,i0,a)') "   aquifer return: replaced by MF6 at channels (recharge map covers ", &
          nhru_cov, " of ", sp_ob%hru, " HRUs)"
        if (nhru_cov < sp_ob%hru) then
          !! Q1130: an uncovered HRU recharges only the SWAT+ aquifer, whose return is no longer
          !! added to channels, so coupling a partial domain would silently drop that baseflow.
          write (*,'(a,i0,a,i0,a)') " MF6 COUPLER: ERROR recharge map covers ", nhru_cov, " of ", &
            sp_ob%hru, " HRUs; the aquifer return of the rest would be dropped. Refusing to couple."
          error stop 3
        end if
      else
        write (*,'(a)') "   aquifer return: replaced by MF6 at channels (WARNING: no recharge map, MF6 recharge is its own)"
      end if
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
      gwf_prev(:) = gwf_arr(1:n_sfr)          !! Q1130 M2: the losses the next carry-back is weighted by
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
      subroutine pull_baseflow_part(w, first)
      !! Q1144: pull_baseflow for ONE MF6 sub-step of a coupled day that took several. Each sub-step adds its exchange
      !! with its share w of the day, so bf_chan is the day's time average (m3/day), the cumulative gain/loss counts
      !! each sub-step's volume once, and gwf_prev (Q1130 M2's carry-back weights) is the day's average GWFLOW.
      !! A one-step day uses pull_baseflow, unchanged.
      implicit none
      real(c_double), intent(in) :: w
      logical, intent(in) :: first
      integer :: e, r
      real(c_double) :: bf
      if (first) then
        bf_chan = 0.0_c_double
        gwf_prev(:) = 0.0_c_double
      end if
      gwf_prev(:) = gwf_prev(:) + w * gwf_arr(1:n_sfr)
      do e = 1, n_bf
        r = bf_reach(e)
        if (r >= 0 .and. r < n_sfr) then
          bf = -gwf_arr(r+1) * w              !! aquifer -> stream, +ve = gain; this sub-step's share
          bf_chan(bf_gis(e)) = bf_chan(bf_gis(e)) + bf
          if (bf > 0.0_c_double) then
            bf_cum_gain = bf_cum_gain + bf
          else
            bf_cum_loss = bf_cum_loss + bf
          end if
        end if
      end do
      end subroutine pull_baseflow_part

      !! ================================================================
      subroutine mf6_coupler_finalize
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
      use ifport, only : chdir          !! gfortran provides CHDIR as a GNU intrinsic
#endif
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
      write (*,'(a,i0,a)') " MF6 COUPLER summary over ", mf6_daycount, " coupled days:"
      if (rch_ready) write (*,'(a,es12.4,a)') &
        "   recharge (sum of daily basin depths) = ", rch_depth_cum, " m"
      if (bf_ready) write (*,'(a,es12.4,a,es12.4,a)') &
        "   SFR exchange: gaining ", bf_cum_gain, "  losing ", bf_cum_loss, " m3"
      !! Q1130: es16.8 (9 significant digits) so each summary number checks against the per-day file and MF6's listing
      if (bf_ready) write (*,'(a,es16.8,a)') &
        "   SWAT+ aquifer return withheld from channels (replaced by MF6) = ", sum(aq_withheld_ob), " m3"
      if (bf_ready) then
        write (*,'(a,es16.8,a,es16.8,a,es16.8,a)') "   channel floor: unmet MF6 loss carried back to MF6 ", &
          cb_carried, " m3, unplaced ", cb_unplaced, " m3, left at end ", sum(unmet_chan), " m3"
        write (*,'(a,i0)') "   channel floor: column-steps with negative RECHARGE after carry-back = ", cb_negrch
        if (cb_iu /= 0) close (cb_iu)
      end if
      if (clk_iu /= 0) then
        write (*,'(a,i0,a,i0,a,i0)') "   MF6 clock: ", clk_days, " coupled days, ", clk_subdays, &
          " with sub-steps, max ", clk_maxsub
        close (clk_iu)
      end if
      if (drv_ready) then
        !! Q1137: es16.8, so each total checks against the per-day file
        write (*,'(a,es16.8,a,es16.8,a,es16.8,a)') "   SFR drive: SWAT+ channel inflow ", drv_cum_swat, &
          " m3; injected INFLOW ", drv_cum_in, " m3, RUNOFF ", drv_cum_ro, " m3"
        write (*,'(a,es16.8,a,es16.8,a)') "   SFR drive: channel inflow with no SFR entry (sent as RUNOFF) ", drv_cum_tolat, &
          " m3; aquifer inflow seen (withheld by M1, so 0) ", drv_cum_aqu, " m3"
        write (*,'(a,es16.8,a,i0,a,i0)') "   SFR drive: two-exit splits moved ", drv_cum_moved, &
          " m3 away from the length split; length fallback pair-days ", drv_fallback, "; days with SFR RAIN /= 0 ", drv_raindays
        if (drv_iu /= 0) close (drv_iu)
      end if
      if (src_ready) write (*,'(a,es12.4,a)') &
        "   PFAS loaded to aquifer    = ", pfas_src_cum, " ug"
      if (pdis_ready) write (*,'(a,es12.4,a)') &
        "   PFAS discharged to stream = ", pfas_disch_cum, " kg"
      mf6_active = .false.
      end subroutine mf6_coupler_finalize

      !! ================================================================
      subroutine mf6_carryback_setup
      !! Q1130 M2: bind what the carry-back needs to place a reach's volume in its RECHARGE
      !! column: SFR reach -> GWF node (IGWFNODE), reduced -> user node (NODEUSER), grid
      !! size and cell widths. If any is missing the carry-back is off and every unmet
      !! volume is reported as UNPLACED (created water), never silently dropped.
      implicit none
      integer :: ios, nodes, nodesuser
      type(c_ptr) :: cp
      integer(c_int), pointer :: sc(:)
      open (newunit=cb_iu, file="mf6_floor.out", status="replace", action="write")
      write (cb_iu,'(a)') "day gis unmet_m3   (Q1130: MF6 loss the SWAT+ channel could not supply)"
      if (.not. rch_ready) then
        write (*,'(a)') "   channel floor: WARNING no RECHARGE binding, unmet loss is reported, not returned"
        return
      end if
      nodes = int_scalar("MODFLOW_SFR/DIS/NODES", ios);         if (ios /= 0) return
      nodesuser = int_scalar("MODFLOW_SFR/DIS/NODESUSER", ios); if (ios /= 0) return
      cb_nrow = int_scalar("MODFLOW_SFR/DIS/NROW", ios);        if (ios /= 0) return
      cb_ncol = int_scalar("MODFLOW_SFR/DIS/NCOL", ios);        if (ios /= 0) return
      if (cb_nrow * cb_ncol /= n_rch) then
        write (*,'(a)') "   channel floor: WARNING NROW*NCOL /= RECHARGE size, carry-back off"
        return
      end if
      cp = c_null_ptr; ios = p_get_value_ptr_int("MODFLOW_SFR/SFR_0/IGWFNODE"//c_null_char, cp)
      if (ios /= 0 .or. .not. c_associated(cp)) return
      call c_f_pointer(cp, cb_igwf, [n_sfr])
      cb_reduced = (nodes < nodesuser)
      if (cb_reduced) then
        cp = c_null_ptr; ios = p_get_value_ptr_int("MODFLOW_SFR/DIS/NODEUSER"//c_null_char, cp)
        if (ios /= 0 .or. .not. c_associated(cp)) return
        call c_f_pointer(cp, cb_nodeuser, [nodes])
      end if
      cp = c_null_ptr; ios = p_get_value_ptr_double("MODFLOW_SFR/DIS/DELR"//c_null_char, cp)
      if (ios /= 0 .or. .not. c_associated(cp)) return
      call c_f_pointer(cp, cb_delr, [cb_ncol])
      cp = c_null_ptr; ios = p_get_value_ptr_double("MODFLOW_SFR/DIS/DELC"//c_null_char, cp)
      if (ios /= 0 .or. .not. c_associated(cp)) return
      call c_f_pointer(cp, cb_delc, [cb_nrow])
      cb_ready = .true.
      write (*,'(a)') "   channel floor: unmet MF6 loss is carried back to MF6 RECHARGE at the next step"
      contains
        integer function int_scalar(addr, ist)
          character(len=*), intent(in) :: addr
          integer, intent(out) :: ist
          type(c_ptr) :: p
          p = c_null_ptr
          ist = p_get_value_ptr_int(trim(addr)//c_null_char, p)
          int_scalar = 0
          if (ist /= 0 .or. .not. c_associated(p)) then
            ist = 1; call trace("carry-back: no "//trim(addr)); return
          end if
          call c_f_pointer(p, sc, [1])
          int_scalar = int(sc(1))
        end function int_scalar
      end subroutine mf6_carryback_setup

      !! ================================================================
      subroutine carry_back_unmet
      !! Q1130 M2: before MF6 solves this step, subtract each channel's unmet volume from
      !! the RECHARGE of the columns under that channel's LOSING reaches of the previous
      !! step (GWFLOW > 0), weighted by their loss. Writes mf6_floor.out.
      implicit none
      integer :: g, e, r, node, un, rc, i, j
      real(c_double) :: tot, v
      do g = 0, max_gis
        if (unmet_chan(g) <= 0.0_c_double) cycle
        if (cb_iu /= 0) write (cb_iu,'(i0,1x,i0,1x,es16.8)') mf6_daycount, g, unmet_chan(g)
        tot = 0.0_c_double
        do e = 1, n_bf
          r = bf_reach(e)
          if (bf_gis(e) == g .and. r >= 0 .and. r < n_sfr) then
            if (gwf_prev(r+1) > 0.0_c_double) tot = tot + gwf_prev(r+1)
          end if
        end do
        if (.not. cb_ready .or. tot <= 0.0_c_double) then
          cb_unplaced = cb_unplaced + unmet_chan(g)
        else
          do e = 1, n_bf
            r = bf_reach(e)
            if (bf_gis(e) /= g .or. r < 0 .or. r >= n_sfr) cycle
            if (gwf_prev(r+1) <= 0.0_c_double) cycle
            v = unmet_chan(g) * gwf_prev(r+1) / tot
            node = cb_igwf(r+1)
            if (node < 1) then
              cb_unplaced = cb_unplaced + v; cycle
            end if
            un = node
            if (cb_reduced) un = cb_nodeuser(node)
            rc = mod(un - 1, cb_nrow * cb_ncol)            !! 0-based cell within its layer
            i = rc / cb_ncol + 1; j = mod(rc, cb_ncol) + 1
            rch_arr(rc+1) = rch_arr(rc+1) - v / (cb_delr(j) * cb_delc(i))
            cb_carried = cb_carried + v
          end do
        end if
        unmet_chan(g) = 0.0_c_double
      end do
      cb_negrch = cb_negrch + count(rch_arr(1:n_rch) < 0.0_c_double)
      end subroutine carry_back_unmet

      !! ================================================================
      subroutine mf6_note_unmet(gis, vol)
      !! Q1130 M2: the channel with this gis floored its flow, discarding vol m3 of MF6
      !! loss. One element per channel gis, so parallel channel routing never shares one.
      implicit none
      integer, intent(in) :: gis
      real, intent(in) :: vol
      if (allocated(unmet_chan)) then
        if (gis >= 0 .and. gis <= max_gis) unmet_chan(gis) = unmet_chan(gis) + real(vol, c_double)
      end if
      end subroutine mf6_note_unmet

      !! ================================================================
      subroutine mf6_withhold_aquifer(iobj, vol)
      !! Q1130 M1: record an aquifer hydrograph (m3) that command_object did NOT add to
      !! channel object iobj because MF6 supplies that channel's groundwater exchange.
      !! One element per receiving object, so parallel channel routing never shares one.
      implicit none
      integer, intent(in) :: iobj
      real, intent(in) :: vol
      if (allocated(aq_withheld_ob)) then
        if (iobj >= 1 .and. iobj <= size(aq_withheld_ob)) &
          aq_withheld_ob(iobj) = aq_withheld_ob(iobj) + real(vol, c_double)
      end if
      end subroutine mf6_withhold_aquifer

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

      logical function mf6_pfas_disch_active()
      !! true when MF6 supplies the channel groundwater-PFAS discharge
      mf6_pfas_disch_active = mf6_active .and. pdis_ready
      end function mf6_pfas_disch_active

      integer function mf6_pfas_compound()
      !! which simulated PFAS compound (index into hcs1%pfas) the GWT carries
      mf6_pfas_compound = pfas_k
      end function mf6_pfas_compound

      real(c_double) function mf6_channel_pfas(gis)
      !! groundwater PFAS mass discharged into the channel (kg/day) for the
      !! coupled compound, from the previous day's MF6 solve.
      integer, intent(in) :: gis
      mf6_channel_pfas = 0.0_c_double
      if (pdis_ready .and. gis >= 0 .and. gis <= max_gis) &
        mf6_channel_pfas = pfas_chan(gis)
      end function mf6_channel_pfas

      !! ================================================================
      subroutine drive_mode_line
      !! Q1137: print the run's SFR drive mode once, as the coupler's first line. The switch is the presence of mf6_drive.map,
      !! which a run's record does not show unless the run says which side it ran on (Lane S): ON names the map's sha256, computed
      !! from the bytes read here, and its header's entry count. A run without this line is INDETERMINATE for every Q1137 row.
      use iso_fortran_env, only : int64
      implicit none
      character(len=64) :: hx
      integer(int64) :: nb
      integer :: iu, ios, nent
      inquire (file="mf6_drive.map", exist=drv_have_map)
      if (.not. drv_have_map) then
        drv_mode = " MF6 COUPLER: SFR drive OFF (no mf6_drive.map)"
      else
        call file_sha256("mf6_drive.map", hx, nb)
        if (len_trim(hx) == 0) call coupler_fatal("mf6_drive.map is present but cannot be read")
        nent = -1
        open (newunit=iu, file="mf6_drive.map", status="old", action="read")
        read (iu,*,iostat=ios) nent
        close (iu)
        if (ios /= 0 .or. nent < 0) call coupler_fatal("mf6_drive.map header unreadable")
        write (drv_mode,'(a,a,a,i0,a)') " MF6 COUPLER: SFR drive ON, mf6_drive.map sha256 ", hx, ", ", nent, " entries"
      end if
      write (*,'(a)') trim(drv_mode)
      end subroutine drive_mode_line

      !! ================================================================
      subroutine mf6_drive_setup
      !! Q1137. With mf6_drive.map: bind SFR's INFLOW, RUNOFF, LENGTH, RAIN and connectivity, check that the deck is cut at channel
      !! boundaries, read the map, and order the channels upstream first. Without it: check that the deck was not cut and its map
      !! lost (Lane S). Every failure is a named stop, because a half-applied drive routes water twice or not at all.
      use hydrograph_module, only : ob
      implicit none
      integer :: iu, ios, e, g, h, n, m, i, k, it, nja, nhdr, nsfr_hdr, r0, u0, ncross, nmulti, nnoob, nrain, nlinks
      integer :: npair, ntwo, nmiss_swat, nmiss_map, ne, in, nforced, ndyn, nstat
      integer, allocatable :: pos(:)
      logical, allocatable :: two(:)
      integer, allocatable :: cnt(:), t_reach(:), t_src(:), t_exit(:), t_g(:)
      real(c_double), allocatable :: t_w(:)
      logical, allocatable :: done(:)
      logical :: progress, ready, found
      real(c_double) :: lsum, w, wsum
      character(len=240) :: msg

      if (.not. bf_ready) then
        if (drv_have_map) call coupler_fatal( &
          "mf6_drive.map is present but mf6_baseflow.map is not: the SFR drive needs the reach -> channel map")
        return
      end if
      !! connectivity, and each reach's key: its channel's gis if mapped, -1 if not (an unmapped reach is its own network)
      call bind_sfr_int("IA", drv_ia, n_sfr + 1)
      nja = drv_ia(n_sfr + 1) - 1
      call bind_sfr_int("JA", drv_ja, nja)
      call bind_sfr_int("IDIR", drv_idir, nja)
      allocate (drv_key(n_sfr)); drv_key = -1
      do e = 1, n_bf
        n = bf_reach(e) + 1
        if (n < 1 .or. n > n_sfr) call coupler_fatal("mf6_baseflow.map names a reach outside SFR")
        if (drv_key(n) /= -1) call coupler_fatal("mf6_baseflow.map maps one SFR reach to two channels")
        drv_key(n) = bf_gis(e)
      end do
      allocate (cnt(0:max_gis)); cnt = 0
      do n = 1, n_sfr
        if (drv_key(n) >= 0) cnt(drv_key(n)) = cnt(drv_key(n)) + 1
      end do
      !! connections between reaches of different keys (each is listed by both of its reaches), and reaches with > 1 downstream
      ncross = 0; nmulti = 0
      do n = 1, n_sfr
        k = 0
        do i = drv_ia(n) + 1, drv_ia(n + 1) - 1
          m = drv_ja(i)
          if (m < 1 .or. m > n_sfr) call coupler_fatal("MF6's SFR JA names a reach outside SFR")
          if (drv_key(m) /= drv_key(n)) ncross = ncross + 1
          if (drv_idir(i) < 0) k = k + 1
        end do
        if (k > 1) nmulti = nmulti + 1
      end do
      ncross = ncross / 2

      if (.not. drv_have_map) then
        !! Lane S: a cut deck whose map was lost. SWAT+ routes water between channels that both have SFR reaches, yet SFR connects
        !! no reach across channels: its inter-channel routing is gone and nothing would replace it.
        nlinks = 0
        do i = 1, size(ob)
          if (ob(i)%typ /= "chandeg" .and. ob(i)%typ /= "cha") cycle
          g = int(ob(i)%gis_id)
          if (g < 0 .or. g > max_gis) cycle
          if (cnt(g) == 0) cycle
          do in = 1, ob(i)%rcv_tot
            if (ob(i)%obtyp_in(in) /= "chandeg" .and. ob(i)%obtyp_in(in) /= "cha") cycle
            h = int(ob(ob(i)%obj_in(in))%gis_id)
            if (h < 0 .or. h > max_gis .or. h == g) cycle
            if (cnt(h) > 0) nlinks = nlinks + 1
          end do
        end do
        if (nlinks > 0 .and. ncross == 0) then
          write (msg,'(a,i0,a)') "SFR connects no reach across channels while SWAT+ routes ", nlinks, &
            " links between channels with SFR reaches: a cut deck without its mf6_drive.map (Q1157 writes both)"
          call coupler_fatal(trim(msg))
        end if
        deallocate (drv_key, cnt)
        return
      end if

      if (gwf_cadence /= 1) call coupler_fatal( &
        "the SFR drive sets one SWAT+ day's inflow per MF6 day; mf6.con gwf_cadence must be 1")
      if (ncross > 0) then
        write (msg,'(a,i0,a)') "mf6_drive.map is present but the deck still connects ", ncross, &
          " SFR reach pairs across channels; a driven deck must be cut at channel boundaries (Q1157)"
        call coupler_fatal(trim(msg))
      end if
      if (nmulti > 0) then
        write (msg,'(a,i0,a)') "the SFR drive needs at most one downstream reach per reach; ", nmulti, " reaches have more"
        call coupler_fatal(trim(msg))
      end if
      call bind_sfr_dbl("INFLOW", drv_inflow)
      call bind_sfr_dbl("RUNOFF", drv_runoff)
      call bind_sfr_dbl("LENGTH", drv_length)
      call bind_sfr_dbl("RAIN", drv_rain)

      !! each reach's exit: follow its single downstream link to the reach that has none
      allocate (drv_exit(n_sfr))
      do n = 1, n_sfr
        m = n; k = 0
        do it = 1, n_sfr
          k = 0
          do i = drv_ia(m) + 1, drv_ia(m + 1) - 1
            if (drv_idir(i) < 0) k = drv_ja(i)
          end do
          if (k == 0) exit
          m = k
        end do
        if (k /= 0) call coupler_fatal("SFR downstream links form a loop; the SFR drive cannot find a reach's exit")
        drv_exit(n) = m
      end do

      !! each channel's reaches and their LENGTH shares (RUNOFF is split by length, DESIGN.txt 3.4)
      allocate (drv_ptr(0:max_gis + 1)); drv_ptr(0) = 1
      do g = 0, max_gis
        drv_ptr(g + 1) = drv_ptr(g) + cnt(g)
      end do
      allocate (drv_reach(drv_ptr(max_gis + 1) - 1), drv_wlen(drv_ptr(max_gis + 1) - 1)); cnt = 0
      do n = 1, n_sfr
        g = drv_key(n)
        if (g < 0) cycle
        drv_reach(drv_ptr(g) + cnt(g)) = n; cnt(g) = cnt(g) + 1
      end do
      do g = 0, max_gis
        if (drv_ptr(g + 1) == drv_ptr(g)) cycle
        lsum = 0.0_c_double
        do k = drv_ptr(g), drv_ptr(g + 1) - 1
          lsum = lsum + drv_length(drv_reach(k))
        end do
        if (.not. (lsum > 0.0_c_double)) call coupler_fatal( &
          "a channel's SFR reaches have no LENGTH; the SFR drive splits a channel's lateral inflow by reach length")
        do k = drv_ptr(g), drv_ptr(g + 1) - 1
          drv_wlen(k) = drv_length(drv_reach(k)) / lsum
        end do
      end do

      !! the SWAT+ channel object of each driven gis
      allocate (gis_ob(0:max_gis)); gis_ob = 0
      do i = 1, size(ob)
        if (ob(i)%typ /= "chandeg" .and. ob(i)%typ /= "cha") cycle
        g = int(ob(i)%gis_id)
        if (g < 0 .or. g > max_gis) cycle
        if (drv_ptr(g + 1) == drv_ptr(g)) cycle
        if (gis_ob(g) /= 0) call coupler_fatal("two SWAT+ channels carry the gis id of one driven SFR channel")
        gis_ob(g) = i
      end do
      nnoob = 0
      do g = 0, max_gis
        if (drv_ptr(g + 1) > drv_ptr(g) .and. gis_ob(g) == 0) nnoob = nnoob + 1
      end do

      !! the map: 'n_entries n_sfr', then per removed connection 'r0 h_gis w u0' (0-based reaches; DESIGN.txt Addenda 2, 3b)
      open (newunit=iu, file="mf6_drive.map", status="old", action="read")
      read (iu,*,iostat=ios) nhdr, nsfr_hdr
      if (ios /= 0 .or. nhdr < 0) call coupler_fatal("mf6_drive.map header unreadable")
      if (nsfr_hdr /= n_sfr) call coupler_fatal("mf6_drive.map was written for another SFR (its reach count differs)")
      n_ent = nhdr
      allocate (t_reach(n_ent), t_src(n_ent), t_exit(n_ent), t_g(n_ent), t_w(n_ent))
      do e = 1, n_ent
        read (iu,*,iostat=ios) r0, h, w, u0
        if (ios /= 0) call coupler_fatal("mf6_drive.map truncated or unreadable (fewer rows than its header count)")
        if (r0 < 0 .or. r0 >= n_sfr .or. u0 < 0 .or. u0 >= n_sfr) call coupler_fatal("mf6_drive.map names a reach outside SFR")
        g = drv_key(r0 + 1)
        if (g < 0) call coupler_fatal("mf6_drive.map: an entry reach is not mapped to a SWAT+ channel")
        if (h < 0 .or. drv_key(u0 + 1) /= h) call coupler_fatal("mf6_drive.map: an exit reach is not in its source channel")
        if (h == g) call coupler_fatal("mf6_drive.map: an entry's source channel is its own channel")
        if (drv_exit(u0 + 1) /= u0 + 1) call coupler_fatal("mf6_drive.map: an exit reach still has a downstream reach")
        if (.not. (w > 0.0_c_double .and. w <= 1.0_c_double)) call coupler_fatal("mf6_drive.map: a weight is outside (0, 1]")
        t_reach(e) = r0 + 1; t_src(e) = h; t_w(e) = w; t_exit(e) = u0 + 1; t_g(e) = g
      end do
      read (iu,*,iostat=ios) r0
      if (ios == 0) call coupler_fatal("mf6_drive.map has more rows than its header count")
      close (iu)
      do e = 1, n_ent
        do k = e + 1, n_ent
          if (t_exit(k) == t_exit(e)) call coupler_fatal("mf6_drive.map: one exit reach feeds two entries")
        end do
      end do
      !! group the entries by receiving channel
      allocate (ent_ptr(0:max_gis + 1)); cnt = 0
      do e = 1, n_ent
        cnt(t_g(e)) = cnt(t_g(e)) + 1
      end do
      ent_ptr(0) = 1
      do g = 0, max_gis
        ent_ptr(g + 1) = ent_ptr(g) + cnt(g)
      end do
      allocate (ent_reach(n_ent), ent_src(n_ent), ent_exit(n_ent), ent_w(n_ent)); cnt = 0
      do e = 1, n_ent
        g = t_g(e); k = ent_ptr(g) + cnt(g); cnt(g) = cnt(g) + 1
        ent_reach(k) = t_reach(e); ent_src(k) = t_src(e); ent_exit(k) = t_exit(e); ent_w(k) = t_w(e)
      end do
      !! each (source, receiver) pair's static weights sum to 1; count pairs, two-exit pairs and where the networks disagree
      npair = 0; ntwo = 0; nmiss_map = 0
      allocate (two(n_ent)); two = .false.
      do g = 0, max_gis
        do e = ent_ptr(g), ent_ptr(g + 1) - 1
          found = .false.
          do k = ent_ptr(g), e - 1
            if (ent_src(k) == ent_src(e)) found = .true.
          end do
          if (found) cycle
          npair = npair + 1; wsum = 0.0_c_double; ne = 0
          do k = e, ent_ptr(g + 1) - 1
            if (ent_src(k) /= ent_src(e)) cycle
            wsum = wsum + ent_w(k); ne = ne + 1
          end do
          if (ne > 1) then
            ntwo = ntwo + 1
            do k = e, ent_ptr(g + 1) - 1
              if (ent_src(k) == ent_src(e)) two(k) = .true.
            end do
          end if
          if (abs(wsum - 1.0_c_double) > 1.0e-9_c_double) call coupler_fatal( &
            "mf6_drive.map: the weights of a source -> receiver pair do not sum to 1")
          found = .false.
          if (gis_ob(g) > 0) then
            i = gis_ob(g)
            do in = 1, ob(i)%rcv_tot
              if (ob(i)%obtyp_in(in) /= "chandeg" .and. ob(i)%obtyp_in(in) /= "cha") cycle
              if (int(ob(ob(i)%obj_in(in))%gis_id) == ent_src(e)) found = .true.
            end do
          end if
          if (.not. found) nmiss_map = nmiss_map + 1
        end do
      end do
      nmiss_swat = 0
      do g = 0, max_gis
        if (gis_ob(g) == 0) cycle
        i = gis_ob(g)
        do in = 1, ob(i)%rcv_tot
          if (ob(i)%obtyp_in(in) /= "chandeg" .and. ob(i)%obtyp_in(in) /= "cha") cycle
          h = int(ob(ob(i)%obj_in(in))%gis_id)
          found = .false.
          do e = ent_ptr(g), ent_ptr(g + 1) - 1
            if (ent_src(e) == h) found = .true.
          end do
          if (.not. found) nmiss_swat = nmiss_swat + 1
        end do
      end do
      !! the order: a two-exit pair splits by the water injected into its source that day, so that source must be complete (every
      !! inflow of its own set) before the receiver is. Only the two-exit pairs impose order; single-exit pairs need nothing. SFR's
      !! reaches can leave a SWAT+ channel and re-enter it (the rogue deck: 234 channels in cycles of its h -> g graph), so a cycle
      !! can remain among these dependencies too. It is broken deterministically (the smallest gis waiting goes next), and a pair
      !! whose source is not yet complete keeps the map's static length weight on every day: counted and reported (Addendum 4).
      allocate (done(0:max_gis), drv_order(max_gis + 1), pos(0:max_gis)); n_order = 0; nforced = 0; pos = 0
      do g = 0, max_gis
        done(g) = (drv_ptr(g + 1) == drv_ptr(g))
      end do
      do
        if (all(done)) exit
        progress = .false.
        do g = 0, max_gis
          if (done(g)) cycle
          ready = .true.
          do e = ent_ptr(g), ent_ptr(g + 1) - 1
            if (two(e) .and. .not. done(ent_src(e))) ready = .false.
          end do
          if (.not. ready) cycle
          done(g) = .true.; n_order = n_order + 1; drv_order(n_order) = g; pos(g) = n_order; progress = .true.
        end do
        if (.not. progress) then
          g = findloc(done, .false., dim=1) - 1         !! done is (0:max_gis); findloc counts from 1
          done(g) = .true.; n_order = n_order + 1; drv_order(n_order) = g; pos(g) = n_order; nforced = nforced + 1
        end if
      end do
      allocate (ent_dyn(n_ent)); ndyn = 0; nstat = 0
      do g = 0, max_gis
        do e = ent_ptr(g), ent_ptr(g + 1) - 1
          ent_dyn(e) = two(e) .and. pos(ent_src(e)) < pos(g)
          if (two(e) .and. ent_dyn(e)) ndyn = ndyn + 1
          if (two(e) .and. .not. ent_dyn(e)) nstat = nstat + 1
        end do
      end do
      allocate (drv_s(n_sfr), inflow_day(n_sfr), runoff_day(n_sfr))
      nrain = count(drv_rain(1:n_sfr) /= 0.0_c_double)
      open (newunit=drv_iu, file="mf6_drive.out", status="replace", action="write")
      write (drv_iu,'(a)') trim(adjustl(drv_mode))
      write (drv_iu,'(a,i0,a)') "day gis swat_in_m3 inflow_m3 runoff_m3 chan_to_runoff_m3 aqu_m3   (Q1137, per MF6 day and " // &
        "driven channel: SWAT+'s inflow of that day; inflow + runoff = swat_in - aqu; SFR reaches with RAIN /= 0 at setup: ", &
        nrain, ")"
      flush (drv_iu)
      write (*,'(a,i0,a,i0,a,i0,a)') "   SFR drive: SWAT+ sets INFLOW/RUNOFF on ", size(drv_reach), " reaches of ", n_order, &
        " channels (", nnoob, " of them have no SWAT+ channel object; their reaches get none)"
      write (*,'(a,i0,a,i0,a,i0,a)') "   SFR drive: ", n_ent, " map entries in ", npair, " source -> receiver pairs, ", ntwo, &
        " split over two or more exits"
      write (*,'(a,i0,a,i0,a,i0,a)') "   SFR drive: two-exit entries split by the day's water ", ndyn, &
        ", by the static length weight ", nstat, " (their source waits on them; ", nforced, " dependency cycles broken)"
      write (*,'(a,i0,a,i0,a,i0)') "   SFR drive: the networks disagree on ", nmiss_swat, &
        " SWAT+ channel links (sent as RUNOFF) and ", nmiss_map, " map pairs (no SWAT+ link); SFR reaches with RAIN /= 0: ", nrain
      deallocate (cnt, t_reach, t_src, t_exit, t_g, t_w, done, pos, two)
      drv_ready = .true.
      end subroutine mf6_drive_setup

      !! ================================================================
      subroutine push_drive(first)
      !! Q1137: set SFR INFLOW and RUNOFF from SWAT+'s inflows of this day, after prepare_time_step and before the solve. The day's
      !! first MF6 step computes them; its later sub-steps re-apply them (Q1144). Each two-exit source is processed before its
      !! receiver where the dependencies allow, so its split follows the water injected into the source that day (DESIGN.txt
      !! Addendum 3b); the map's length weight is used for a pair whose source waits on it (a cycle, Addendum 4) and, counted, on a
      !! day the source received none.
      use hydrograph_module, only : ob
      implicit none
      logical, intent(in) :: first
      integer :: k, kk, g, iob, in, e, h, ne, r
      real(c_double) :: q, swat, lat, tolat, aqu, ssum, wd, inj, ro
      logical :: dyn
      if (.not. first) then
        drv_inflow(1:n_sfr) = inflow_day(:)
        drv_runoff(1:n_sfr) = runoff_day(:)
        return
      end if
      if (any(drv_rain(1:n_sfr) /= 0.0_c_double)) drv_raindays = drv_raindays + 1
      drv_s(:) = 0.0_c_double
      do k = 1, size(drv_reach)
        drv_inflow(drv_reach(k)) = 0.0_c_double
        drv_runoff(drv_reach(k)) = 0.0_c_double
      end do
      do k = 1, n_order
        g = drv_order(k); iob = gis_ob(g)
        swat = 0.0_c_double; lat = 0.0_c_double; tolat = 0.0_c_double; aqu = 0.0_c_double
        if (iob > 0) then
          if (allocated(ob(iob)%hin_d)) then
            do in = 1, ob(iob)%rcv_tot
              q = real(ob(iob)%hin_d(in)%flo, c_double)
              swat = swat + q
              if (ob(iob)%obtyp_in(in) == "aqu") then
                aqu = aqu + q                          !! Q1130 M1 withholds it (hin_d = 0): injected nowhere
                cycle
              end if
              if (ob(iob)%obtyp_in(in) /= "chandeg" .and. ob(iob)%obtyp_in(in) /= "cha") then
                lat = lat + q                          !! hru, ru, res, recall, ...: lateral
                cycle
              end if
              h = int(ob(ob(iob)%obj_in(in))%gis_id)
              ne = 0; ssum = 0.0_c_double; dyn = .false.
              do e = ent_ptr(g), ent_ptr(g + 1) - 1
                if (ent_src(e) /= h) cycle
                ne = ne + 1; ssum = ssum + drv_s(ent_exit(e)); dyn = ent_dyn(e)   !! one pair: one flag
              end do
              if (ne == 0) then
                lat = lat + q; tolat = tolat + q       !! SWAT+ links h -> g, SFR has no entry for it: lateral, and counted
                cycle
              end if
              if (ne > 1 .and. dyn .and. .not. (ssum > 0.0_c_double)) drv_fallback = drv_fallback + 1
              do e = ent_ptr(g), ent_ptr(g + 1) - 1
                if (ent_src(e) /= h) cycle
                if (ne == 1) then
                  wd = 1.0_c_double
                else if (ent_dyn(e) .and. ssum > 0.0_c_double) then
                  wd = drv_s(ent_exit(e)) / ssum
                  drv_cum_moved = drv_cum_moved + 0.5_c_double * abs(wd - ent_w(e)) * q
                else
                  wd = ent_w(e)
                end if
                drv_inflow(ent_reach(e)) = drv_inflow(ent_reach(e)) + wd * q
              end do
            end do
          end if
        end if
        inj = 0.0_c_double; ro = 0.0_c_double
        do kk = drv_ptr(g), drv_ptr(g + 1) - 1
          r = drv_reach(kk)
          drv_runoff(r) = lat * drv_wlen(kk)
          inj = inj + drv_inflow(r); ro = ro + drv_runoff(r)
          drv_s(drv_exit(r)) = drv_s(drv_exit(r)) + drv_inflow(r) + drv_runoff(r)
        end do
        if (iob > 0) then
          if (drv_iu /= 0) write (drv_iu,'(i0,1x,i0,5(1x,es16.8))') mf6_daycount, g, swat, inj, ro, tolat, aqu
          drv_cum_swat = drv_cum_swat + swat; drv_cum_in = drv_cum_in + inj; drv_cum_ro = drv_cum_ro + ro
          drv_cum_tolat = drv_cum_tolat + tolat; drv_cum_aqu = drv_cum_aqu + aqu
        end if
      end do
      if (drv_iu /= 0) flush (drv_iu)
      inflow_day(:) = drv_inflow(1:n_sfr)
      runoff_day(:) = drv_runoff(1:n_sfr)
      end subroutine push_drive

      !! ================================================================
      subroutine bind_sfr_dbl(name, p)
      !! Q1137: bind a per-reach double array of the SFR package GWFLOW belongs to; a named stop if MF6 does not expose it.
      character(len=*), intent(in) :: name
      real(c_double), pointer, intent(inout) :: p(:)
      type(c_ptr) :: cp
      integer :: ios
      character(len=256) :: addr
      addr = gwf_addr(1:index(gwf_addr, "/", back=.true.))//name
      cp = c_null_ptr
      ios = p_get_value_ptr_double(trim(addr)//c_null_char, cp)
      if (ios /= 0 .or. .not. c_associated(cp)) call coupler_fatal("SFR drive: cannot bind MF6 variable "//trim(addr))
      call c_f_pointer(cp, p, [n_sfr])
      end subroutine bind_sfr_dbl

      subroutine bind_sfr_int(name, p, n)
      !! Q1137: bind an integer array of that SFR package (IA, JA, IDIR) with n elements; a named stop if MF6 does not expose it.
      character(len=*), intent(in) :: name
      integer(c_int), pointer, intent(inout) :: p(:)
      integer, intent(in) :: n
      type(c_ptr) :: cp
      integer :: ios
      character(len=256) :: addr
      addr = gwf_addr(1:index(gwf_addr, "/", back=.true.))//name
      cp = c_null_ptr
      ios = p_get_value_ptr_int(trim(addr)//c_null_char, cp)
      if (ios /= 0 .or. .not. c_associated(cp)) call coupler_fatal("SFR drive: cannot bind MF6 variable "//trim(addr))
      if (n < 1) call coupler_fatal("SFR drive: MF6 variable "//trim(addr)//" is empty")
      call c_f_pointer(cp, p, [n])
      end subroutine bind_sfr_int

      !! ================================================================
      subroutine file_sha256(path, hex, nbytes)
      !! Q1137: SHA-256 (FIPS 180-4) of a file's bytes, as 64 lowercase hex characters, so the run log names the exact map the
      !! coupler read. 32-bit words are held in int64 and masked, so no signed overflow can occur. hex = "" if the file cannot be
      !! read. Checked against sha256sum on files of 0-100000 bytes, every padding boundary included (q1137 scratch test).
      use iso_fortran_env, only : int8, int64
      character(len=*), intent(in) :: path
      character(len=64), intent(out) :: hex
      integer(int64), intent(out) :: nbytes
      integer(int64), parameter :: M32 = 4294967295_int64
      integer(int64), parameter :: K(64) = [ &
        int(z'428a2f98',int64), int(z'71374491',int64), int(z'b5c0fbcf',int64), int(z'e9b5dba5',int64), &
        int(z'3956c25b',int64), int(z'59f111f1',int64), int(z'923f82a4',int64), int(z'ab1c5ed5',int64), &
        int(z'd807aa98',int64), int(z'12835b01',int64), int(z'243185be',int64), int(z'550c7dc3',int64), &
        int(z'72be5d74',int64), int(z'80deb1fe',int64), int(z'9bdc06a7',int64), int(z'c19bf174',int64), &
        int(z'e49b69c1',int64), int(z'efbe4786',int64), int(z'0fc19dc6',int64), int(z'240ca1cc',int64), &
        int(z'2de92c6f',int64), int(z'4a7484aa',int64), int(z'5cb0a9dc',int64), int(z'76f988da',int64), &
        int(z'983e5152',int64), int(z'a831c66d',int64), int(z'b00327c8',int64), int(z'bf597fc7',int64), &
        int(z'c6e00bf3',int64), int(z'd5a79147',int64), int(z'06ca6351',int64), int(z'14292967',int64), &
        int(z'27b70a85',int64), int(z'2e1b2138',int64), int(z'4d2c6dfc',int64), int(z'53380d13',int64), &
        int(z'650a7354',int64), int(z'766a0abb',int64), int(z'81c2c92e',int64), int(z'92722c85',int64), &
        int(z'a2bfe8a1',int64), int(z'a81a664b',int64), int(z'c24b8b70',int64), int(z'c76c51a3',int64), &
        int(z'd192e819',int64), int(z'd6990624',int64), int(z'f40e3585',int64), int(z'106aa070',int64), &
        int(z'19a4c116',int64), int(z'1e376c08',int64), int(z'2748774c',int64), int(z'34b0bcb5',int64), &
        int(z'391c0cb3',int64), int(z'4ed8aa4a',int64), int(z'5b9cca4f',int64), int(z'682e6ff3',int64), &
        int(z'748f82ee',int64), int(z'78a5636f',int64), int(z'84c87814',int64), int(z'8cc70208',int64), &
        int(z'90befffa',int64), int(z'a4506ceb',int64), int(z'bef9a3f7',int64), int(z'c67178f2',int64) ]
      character(len=16), parameter :: HD = "0123456789abcdef"
      integer(int64) :: H(8), W(64), a, b, c, d, e, f, g, hh, t1, t2, fsize, total, bitlen
      integer(int8), allocatable :: buf(:)
      integer(int64), allocatable :: m(:)
      integer :: iu, ios, i, j, blk, o
      hex = ""; nbytes = -1
      open (newunit=iu, file=path, access="stream", form="unformatted", action="read", status="old", iostat=ios)
      if (ios /= 0) return
      inquire (unit=iu, size=fsize)
      allocate (buf(max(fsize, 1_int64)))
      if (fsize > 0) read (iu, iostat=ios) buf(1:fsize)
      close (iu)
      if (ios /= 0 .or. fsize < 0) return
      nbytes = fsize
      total = ((fsize + 8) / 64 + 1) * 64                    !! message + 0x80 + zeros + 8-byte length, a multiple of 64
      allocate (m(total)); m = 0
      do i = 1, int(fsize)
        m(i) = iand(int(buf(i), int64), 255_int64)
      end do
      m(fsize + 1) = 128
      bitlen = fsize * 8
      do j = 0, 7
        m(total - j) = iand(ishft(bitlen, -8 * j), 255_int64)
      end do
      H = [ int(z'6a09e667',int64), int(z'bb67ae85',int64), int(z'3c6ef372',int64), int(z'a54ff53a',int64), &
            int(z'510e527f',int64), int(z'9b05688c',int64), int(z'1f83d9ab',int64), int(z'5be0cd19',int64) ]
      do blk = 0, int(total / 64) - 1
        o = blk * 64
        do i = 1, 16
          W(i) = ior(ior(ishft(m(o+4*i-3), 24), ishft(m(o+4*i-2), 16)), ior(ishft(m(o+4*i-1), 8), m(o+4*i)))
        end do
        do i = 17, 64
          W(i) = iand(sig1(W(i-2)) + W(i-7) + sig0(W(i-15)) + W(i-16), M32)
        end do
        a = H(1); b = H(2); c = H(3); d = H(4); e = H(5); f = H(6); g = H(7); hh = H(8)
        do i = 1, 64
          t1 = iand(hh + bsig1(e) + ieor(iand(e, f), iand(iand(not(e), M32), g)) + K(i) + W(i), M32)
          t2 = iand(bsig0(a) + ieor(ieor(iand(a, b), iand(a, c)), iand(b, c)), M32)
          hh = g; g = f; f = e; e = iand(d + t1, M32); d = c; c = b; b = a; a = iand(t1 + t2, M32)
        end do
        H(1) = iand(H(1) + a, M32); H(2) = iand(H(2) + b, M32); H(3) = iand(H(3) + c, M32); H(4) = iand(H(4) + d, M32)
        H(5) = iand(H(5) + e, M32); H(6) = iand(H(6) + f, M32); H(7) = iand(H(7) + g, M32); H(8) = iand(H(8) + hh, M32)
      end do
      do i = 1, 8
        do j = 1, 8
          o = int(iand(ishft(H(i), -4 * (8 - j)), 15_int64)) + 1
          hex(8*(i-1)+j:8*(i-1)+j) = HD(o:o)
        end do
      end do
      contains
        integer(int64) function rotr(x, n)
          integer(int64), intent(in) :: x
          integer, intent(in) :: n
          rotr = iand(ior(ishft(x, -n), ishft(x, 32 - n)), M32)
        end function rotr
        integer(int64) function bsig0(x)
          integer(int64), intent(in) :: x
          bsig0 = ieor(ieor(rotr(x, 2), rotr(x, 13)), rotr(x, 22))
        end function bsig0
        integer(int64) function bsig1(x)
          integer(int64), intent(in) :: x
          bsig1 = ieor(ieor(rotr(x, 6), rotr(x, 11)), rotr(x, 25))
        end function bsig1
        integer(int64) function sig0(x)
          integer(int64), intent(in) :: x
          sig0 = ieor(ieor(rotr(x, 7), rotr(x, 18)), ishft(x, -3))
        end function sig0
        integer(int64) function sig1(x)
          integer(int64), intent(in) :: x
          sig1 = ieor(ieor(rotr(x, 17), rotr(x, 19)), ishft(x, -10))
        end function sig1
      end subroutine file_sha256

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

      subroutine clock_stop(t0, tnow, tday, nsub, msg)
      !! Q1144: this coupled day cannot be completed on MF6's own clock. Print where the clock stands, keep the
      !! day rows already written, and stop with a named error.
      real(c_double), intent(in) :: t0, tnow, tday
      integer, intent(in) :: nsub
      character(len=*), intent(in) :: msg
      write (*,'(a,i0,a,f16.6,a,f16.6,a,f16.6,a,i0)') " MF6 COUPLER: clock day ", mf6_daycount, &
        "  MF6 t before ", t0, "  now ", tnow, "  day end ", tday, "  sub-steps ", nsub
      if (clk_iu /= 0) flush (clk_iu)
      call coupler_fatal(msg)
      end subroutine clock_stop

      subroutine coupler_fatal(msg)
      !! Q1040(a): mf6.con requested MODFLOW 6 coupling and it cannot be delivered. Stop with a named
      !! error and a non-zero status instead of continuing UNCOUPLED with a success status.
      use iso_fortran_env, only : error_unit
      character(len=*), intent(in) :: msg
      write (*,'(a)') " MF6 COUPLER FATAL: "//msg
      write (error_unit,'(a)') " MF6 COUPLER FATAL: "//msg
      flush (error_unit)
      error stop 2
      end subroutine coupler_fatal

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
      if (.not. resolve_ptr("get_value_ptr_int", p_get_value_ptr_int)) return
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
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
      use ifport, only : getcwd          !! gfortran provides GETCWD as a GNU intrinsic
#endif
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
