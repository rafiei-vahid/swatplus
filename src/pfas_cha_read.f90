      subroutine pfas_cha_read

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Read the per-PFAS in-stream (channel) transport parameters and the
!!    optional initial reach water/benthic PFAS concentrations, then ALLOCATE
!!    and INITIALIZE every container the serial PFAS channel phase (pfas_cha)
!!    needs:
!!      (1) pfas_chadb(:)            per-PFAS routing params (koc, settle,
!!                                   resus, bury, active bed depth)
!!      (2) ch_pfas_water/benthic    per-channel water + bed PFAS pools (kg),
!!          (+ _init copies)         %pfas dimensioned by npfas
!!      (3) chpfas_d/m/y/a(:)        daily/monthly/yearly/avg reach balance,
!!          chpfas / chpfasz         working + zero accumulators
!!      (4) the %pfas slot on every constituent hydrograph that the channel
!!          command carries: hcs1/hcs2/hcs3, hin_csz, and obcs(:)%hin/%hd/
!!          %hin_sur/%hin_lat/%hin_til  (so hcs1%pfas can be loaded from the
!!          HRU and obcs(icmd)%hd(1)%pfas can receive the routed mass)
!!
!!    This is the PFAS analogue of pest_cha_res_read.f90 + the pesticide
!!    channel-pool allocation in sd_channel_read.f90 and the pesticide
!!    hydrograph allocation in hyd_read_connect.f90 / hyd_connect.f90.  It is
!!    called ONCE during initialization, AFTER pfas_read (npfas + pfasdb set)
!!    and AFTER sd_channel_read (sp_ob%chandeg + sd_ch geometry set) and AFTER
!!    hyd_read_connect (obcs allocated).
!!
!!    SWAT+ reader conventions followed:
!!      * inquire(file=...) guard so a missing file is a clean no-op that
!!        still routes PFAS with built-in PFOS-like defaults,
!!      * title + header line skips, list-directed (free) reads,
!!      * allocation off sp_ob%chandeg / npfas with db_mx bookkeeping.
!!
!!    ~ ~ ~ INPUT FILE : "pfas_cha.dat"  (per-PFAS in-stream parameters) ~ ~ ~
!!    Analogous to pesticide.pes (aquatic block) re-cast per PFAS.  Two header
!!    lines, then one record per PFAS compound.  List-directed columns:
!!
!!      id  name        koc        settle    resus     bury     act_dep    &
!!      water_ppt  benthic_ppt
!!      (-) (a16)     m^3/g       m/day     m/day    m/day      m         &
!!        ng/L       ng/g
!!
!!      id          sequential PFAS index (1..npfas); a record with id<=0
!!                  ends the list (mirrors the readpfas.f sentinel).  id maps
!!                  to the database compound via pfas_num (identity crosswalk).
!!      name        compound name (informational; copied to pfas_chadb%name)
!!      koc         linear water-sediment partition (m^3/g)  -> %koc
!!      settle      settling velocity of sorbed PFAS (m/day) -> %aq_settle
!!      resus       resuspension velocity (m/day)            -> %aq_resus
!!      bury        burial velocity in bed sediment (m/day)  -> %ben_bury
!!      act_dep     active bed-sediment layer depth (m)      -> %ben_act_dep
!!      water_ppt   initial reach water-column conc (ng/L = ug/m^3); converted
!!                  to kg/m^3 by 1.e-12 then multiplied by reach water volume
!!      benthic_ppt initial bed-sediment conc (ng/g = ug/kg); converted to
!!                  kg/kg by 1.e-12 then multiplied by bed sediment mass
!!
!!    A record may omit the trailing water_ppt / benthic_ppt columns; the
!!    list-directed read leaves them at their default (0.) and the reach
!!    starts PFAS-free.  If the whole file is absent, ALL compounds get the
!!    PFOS-like defaults below and the reach starts PFAS-free, so a model with
!!    no pfas_cha.dat still routes the HRU loads through pfas_cha correctly.
!!
!!    Default (PFOS-like) in-stream parameters, used for any compound not
!!    found in the file (and for every compound when the file is absent).
!!    These mirror the rtpfas.f / pesticide aquatic defaults: a small
!!    settling/resuspension pair, modest burial, 0.1 m active layer, and a
!!    Koc so that with ~1% channel carbon and a few g/m^3 of suspended
!!    sediment most PFAS stays soluble (frsol ~ 1):
!!      koc=1.e-5 m^3/g, settle=1.0, resus=0.05, bury=0.001 m/day,
!!      act_dep=0.1 m.  Initial water/benthic conc default to zero.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      use pfas_module, only : npfas, pfas_num, pfasdb, pfas_koc_scale
      use pfas_cha_module
      use constituent_mass_module
      use sd_channel_module, only : sd_ch
      use hydrograph_module, only : sp_ob, ch_stor

      implicit none

      character (len=80) :: titldum = ""   !!        |title line of file
      character (len=80) :: header = ""    !!        |header line of file
      integer :: eof = 0                   !!        |end-of-file flag
      integer :: imax = 0                  !! none   |count of compounds overridden from file
      integer :: id = 0                    !! none   |PFAS id read from a record
      integer :: ipf = 0                   !! none   |sequential PFAS counter
      integer :: jpf = 0                   !! none   |database PFAS index
      integer :: ich = 0                   !! none   |channel counter
      integer :: iob = 0                   !! none   |spatial-object counter
      integer :: ihyd = 0                  !! none   |hydrograph counter
      integer :: nhyds = 0                 !! none   |number of hyds on an object
      logical :: i_exist = .false.         !! none   |file-existence flag
      real :: water_ppt = 0.               !! ng/L   |initial water conc (scratch)
      real :: benthic_ppt = 0.             !! ng/g   |initial benthic conc (scratch)
      real :: bedmass = 0.                 !! kg     |active-bed sediment mass
      real :: por = 0.                     !! none   |porosity of bottom sediments
      type (pfas_cha_db) :: dbrec          !!        |scratch param record

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    nothing to do if no PFAS compounds are simulated
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      if (npfas <= 0) return

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    (1) ALLOCATE per-PFAS routing database + seed PFOS-like defaults
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      if (.not. allocated(pfas_chadb)) allocate (pfas_chadb(npfas))
      do ipf = 1, npfas
        pfas_chadb(ipf)%name        = ""
        pfas_chadb(ipf)%koc         = 1.e-5
        pfas_chadb(ipf)%aq_settle   = 1.0
        pfas_chadb(ipf)%aq_resus    = 0.05
        pfas_chadb(ipf)%ben_bury    = 0.001
        pfas_chadb(ipf)%ben_act_dep = 0.1
      end do

      !! per-PFAS initial reach concentrations (default reach starts clean)
      if (.not. allocated(pfas_water_ini)) allocate (pfas_water_ini(npfas))
      do ipf = 1, npfas
        if (.not. allocated(pfas_water_ini(ipf)%water))                   &
     &    allocate (pfas_water_ini(ipf)%water(1), source = 0.)
        if (.not. allocated(pfas_water_ini(ipf)%benthic))                 &
     &    allocate (pfas_water_ini(ipf)%benthic(1), source = 0.)
        pfas_water_ini(ipf)%name = ""
        pfas_water_ini(ipf)%water(1) = 0.
        pfas_water_ini(ipf)%benthic(1) = 0.
      end do

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    (2) READ pfas_cha.dat (overrides the defaults for listed compounds)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      eof = 0
      imax = 0
      inquire (file="pfas_cha.dat", exist=i_exist)
      if (i_exist) then
        do
          open (172, file="pfas_cha.dat")
          read (172,*,iostat=eof) titldum
          if (eof < 0) exit
          read (172,*,iostat=eof) header
          if (eof < 0) exit
          do
            !! defaults so omitted trailing columns stay sane
            dbrec%name        = ""
            dbrec%koc         = 1.e-5
            dbrec%aq_settle   = 1.0
            dbrec%aq_resus    = 0.05
            dbrec%ben_bury    = 0.001
            dbrec%ben_act_dep = 0.1
            water_ppt   = 0.
            benthic_ppt = 0.
            read (172,*,iostat=eof) id, dbrec%name, dbrec%koc,            &
     &        dbrec%aq_settle, dbrec%aq_resus, dbrec%ben_bury,            &
     &        dbrec%ben_act_dep, water_ppt, benthic_ppt
            if (eof < 0) exit
            if (id <= 0) exit
            if (id > npfas) cycle
            pfas_chadb(id)%name        = dbrec%name
            pfas_chadb(id)%koc         = dbrec%koc
            pfas_chadb(id)%aq_settle   = dbrec%aq_settle
            pfas_chadb(id)%aq_resus    = dbrec%aq_resus
            pfas_chadb(id)%ben_bury    = dbrec%ben_bury
            if (dbrec%ben_act_dep > 1.e-6) then
              pfas_chadb(id)%ben_act_dep = dbrec%ben_act_dep
            end if
            pfas_water_ini(id)%name       = dbrec%name
            pfas_water_ini(id)%water(1)   = water_ppt    !! ng/L
            pfas_water_ini(id)%benthic(1) = benthic_ppt  !! ng/g
            imax = imax + 1
          end do
          exit
        end do
        close (172)
      end if
      !! imax (count of compounds overridden from file) is informational only;
      !! no db_mx slot is required since pfas_chadb is sized by npfas, not by a
      !! separate file-record max.  Add a db_mx%pfas_cha field if a count is
      !! wanted in the input-summary output.

      !! calibration: scale the in-stream koc (partition) by pfas_koc_scale
      do ipf = 1, npfas
        pfas_chadb(ipf)%koc = pfas_chadb(ipf)%koc * pfas_koc_scale
      end do

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    (3) ALLOCATE per-channel water + benthic PFAS pools and the daily/
!!        monthly/yearly/avg reach-balance accumulators.  Channel index 0 is
!!        kept (0:chandeg) to mirror ch_water / chpst_d allocation.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      if (.not. allocated(ch_pfas_water))                                 &
     &  allocate (ch_pfas_water(0:sp_ob%chandeg))
      if (.not. allocated(ch_pfas_benthic))                              &
     &  allocate (ch_pfas_benthic(0:sp_ob%chandeg))
      if (.not. allocated(ch_pfas_water_init))                            &
     &  allocate (ch_pfas_water_init(0:sp_ob%chandeg))
      if (.not. allocated(ch_pfas_benthic_init))                          &
     &  allocate (ch_pfas_benthic_init(0:sp_ob%chandeg))

      if (.not. allocated(chpfas_d)) allocate (chpfas_d(0:sp_ob%chandeg))
      if (.not. allocated(chpfas_m)) allocate (chpfas_m(0:sp_ob%chandeg))
      if (.not. allocated(chpfas_y)) allocate (chpfas_y(0:sp_ob%chandeg))
      if (.not. allocated(chpfas_a)) allocate (chpfas_a(0:sp_ob%chandeg))

      !! working + zero accumulators (scalar pfas_cha_output)
      if (.not. allocated(chpfas%pfas))  allocate (chpfas%pfas(npfas))
      if (.not. allocated(chpfasz%pfas)) allocate (chpfasz%pfas(npfas))
      do ipf = 1, npfas
        chpfas%pfas(ipf)  = ch_pfasbz
        chpfasz%pfas(ipf) = ch_pfasbz
      end do

      do ich = 0, sp_ob%chandeg
        allocate (ch_pfas_water(ich)%pfas(npfas),         source = 0.)
        allocate (ch_pfas_benthic(ich)%pfas(npfas),       source = 0.)
        allocate (ch_pfas_water_init(ich)%pfas(npfas),    source = 0.)
        allocate (ch_pfas_benthic_init(ich)%pfas(npfas),  source = 0.)
        allocate (chpfas_d(ich)%pfas(npfas))
        allocate (chpfas_m(ich)%pfas(npfas))
        allocate (chpfas_y(ich)%pfas(npfas))
        allocate (chpfas_a(ich)%pfas(npfas))
        do ipf = 1, npfas
          chpfas_d(ich)%pfas(ipf) = ch_pfasbz
          chpfas_m(ich)%pfas(ipf) = ch_pfasbz
          chpfas_y(ich)%pfas(ipf) = ch_pfasbz
          chpfas_a(ich)%pfas(ipf) = ch_pfasbz
        end do
      end do

      !! initialize the reach water / bed pools from the read concentrations.
      !!   water  pool (kg) = (ng/L)*1.e-12 (kg/L -> kg/m^3 since 1 L=1.e-3 m^3?
      !!     ng/L = 1.e-9 g/L = 1.e-12 kg/L = 1.e-9 kg/m^3) * water volume(m^3)
      !!   benthic pool(kg) = (ng/g)*1.e-12 (kg/kg) * active-bed sediment mass(kg)
      !! active-bed sediment mass = chw(m) * chl(km)*1000(m) * act_dep(m)
      !!                            * ch_bd(t/m^3)*1000(kg/t) * (1-porosity)
      do ich = 1, sp_ob%chandeg
        por = 1. - sd_ch(ich)%ch_bd / 2.65
        if (por < 0.) por = 0.
        if (por > 1.) por = 1.
        !! PFAS aquatic mixing velocity (diffusion), dimensioned by npfas.
        !! Same formula as sd_hydsed_init for pesticides (aq_mix), but with PFAS
        !! molar mass converted kg/mol -> g/mol (mw*1000) so it matches a
        !! pesticide of identical mol_wt exactly (ch_rtpest cross-validation).
        if (.not. allocated(sd_ch(ich)%aq_mix_pfas))                       &
     &    allocate (sd_ch(ich)%aq_mix_pfas(npfas), source = 0.)
        do ipf = 1, npfas
          jpf = pfas_num(ipf)
          sd_ch(ich)%aq_mix_pfas(ipf) = (pfasdb(jpf)%mw * 1000.)           &
     &        ** (-.6666) * (1. - sd_ch(ich)%ch_bd / 2.65) * (69.35 / 365)
          !! water column (kg): ng/L -> kg/m^3 is 1.e-9
          ch_pfas_water(ich)%pfas(ipf) =                                  &
     &        pfas_water_ini(jpf)%water(1) * 1.e-9 * ch_stor(ich)%flo
          !! bed sediment mass in active layer (kg)
          bedmass = sd_ch(ich)%chw * (sd_ch(ich)%chl * 1000.)             &
     &              * pfas_chadb(ipf)%ben_act_dep                          &
     &              * (sd_ch(ich)%ch_bd * 1000.) * (1. - por)
          !! benthic (kg): ng/g -> kg/kg is 1.e-12
          ch_pfas_benthic(ich)%pfas(ipf) =                                &
     &        pfas_water_ini(jpf)%benthic(1) * 1.e-12 * bedmass
          ch_pfas_water_init(ich)%pfas(ipf) = ch_pfas_water(ich)%pfas(ipf)
          ch_pfas_benthic_init(ich)%pfas(ipf) =                           &
     &        ch_pfas_benthic(ich)%pfas(ipf)
        end do
      end do

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    (3b) PFAS POINT SOURCES (optional pfas_source.dat): per-channel constant
!!         daily load (kg/day) -- WWTP effluent, contaminated-site leachate, AFFF.
!!         Format: title, header, then "channel pfas_id load_kgday [name]" rows;
!!         a row with channel<=0 (or EOF) ends the list. Missing file -> no sources.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      if (.not. allocated(pfas_src_load))                                   &
     &  allocate (pfas_src_load(sp_ob%chandeg, npfas), source = 0.)
      inquire (file="pfas_source.dat", exist=i_exist)
      if (i_exist) then
        open (173, file="pfas_source.dat")
        read (173,*,iostat=eof) titldum
        read (173,*,iostat=eof) header
        do
          read (173,*,iostat=eof) ich, id, water_ppt     !ich=channel, id=pfas, water_ppt reused as load kg/day
          if (eof < 0) exit
          if (ich <= 0) exit
          if (ich <= sp_ob%chandeg .and. id >= 1 .and. id <= npfas)         &
     &      pfas_src_load(ich, id) = pfas_src_load(ich, id) + water_ppt
        end do
        close (173)
      end if

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    (4) ALLOCATE the %pfas slot on the constituent hydrographs that the
!!        channel command carries.  Mirrors the %pest allocation in
!!        hyd_connect.f90 (hcs*, hin_csz) and hyd_read_connect.f90 (obcs).
!!        PFAS allocation is INDEPENDENT of cs_db%num_tot, so it loops over
!!        all objects here using obcs_alloc as the "this object has obcs"
!!        guard (set in hyd_read_connect when cs_db%num_tot > 0).  If no
!!        generic constituents exist, obcs may be unallocated; the PFAS HRU
!!        load path uses obcs, so obcs must exist whenever npfas>0 -- this is
!!        ensured by the cs_db%num_tot bump done when PFAS is active (see
!!        integration note).  We guard every access with allocated().
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      if (.not. allocated(hcs1%pfas)) allocate (hcs1%pfas(npfas), source=0.)
      if (.not. allocated(hcs2%pfas)) allocate (hcs2%pfas(npfas), source=0.)
      if (.not. allocated(hcs3%pfas)) allocate (hcs3%pfas(npfas), source=0.)
      if (.not. allocated(hin_csz%pfas))                                  &
     &  allocate (hin_csz%pfas(npfas), source=0.)
      hcs1%pfas = 0.
      hcs2%pfas = 0.
      hcs3%pfas = 0.
      hin_csz%pfas = 0.

      if (allocated(obcs)) then
        do iob = 1, sp_ob%objs
          if (.not. allocated(obcs_alloc)) exit
          if (obcs_alloc(iob) /= 1) cycle

          if (allocated(obcs(iob)%hin)) then
            if (.not. allocated(obcs(iob)%hin(1)%pfas))                   &
     &        allocate (obcs(iob)%hin(1)%pfas(npfas), source = 0.)
          end if
          if (allocated(obcs(iob)%hin_sur)) then
            if (.not. allocated(obcs(iob)%hin_sur(1)%pfas))               &
     &        allocate (obcs(iob)%hin_sur(1)%pfas(npfas), source = 0.)
          end if
          if (allocated(obcs(iob)%hin_lat)) then
            if (.not. allocated(obcs(iob)%hin_lat(1)%pfas))               &
     &        allocate (obcs(iob)%hin_lat(1)%pfas(npfas), source = 0.)
          end if
          if (allocated(obcs(iob)%hin_til)) then
            if (.not. allocated(obcs(iob)%hin_til(1)%pfas))               &
     &        allocate (obcs(iob)%hin_til(1)%pfas(npfas), source = 0.)
          end if
          if (allocated(obcs(iob)%hin_aqu)) then
            if (.not. allocated(obcs(iob)%hin_aqu(1)%pfas))               &
     &        allocate (obcs(iob)%hin_aqu(1)%pfas(npfas), source = 0.)
          end if

          if (allocated(obcs(iob)%hd)) then
            nhyds = size(obcs(iob)%hd)
            do ihyd = 1, nhyds
              if (.not. allocated(obcs(iob)%hd(ihyd)%pfas))               &
     &          allocate (obcs(iob)%hd(ihyd)%pfas(npfas), source = 0.)
            end do
          end if
        end do
      end if

      return
      end subroutine pfas_cha_read
