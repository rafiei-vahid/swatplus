      subroutine pfas_read

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Read PFAS input for surface-water-only PFAS fate-and-transport in
!!    modern (free-form) SWAT+ and allocate every module container in
!!    pfas_module.  This is the SWAT+ port of the SWAT2012 readpfas.f
!!    (global pfas.dat database) plus the *.chm PFAS block of readchm.f
!!    (per-HRU per-layer soil pools).  Scope is surface water only: no
!!    groundwater / aquifer PFAS state is read here.
!!
!!    It performs three jobs, in order:
!!      (1) read the global per-PFAS database  -> pfasdb(:) + npfas
!!      (2) allocate the per-HRU soil column   -> pfas_soil_hru(:)%ly(:)
!!      (3) read per-HRU per-layer initial soil pools, Freundlich kf/n,
!!          enrichment, contaminated-site count and grain diameter.
!!
!!    SWAT+ reader conventions followed (see pest_parm_read.f90,
!!    cs_hru_read.f90, pesticide_init.f90):
!!      * inquire(file=...) guard so a missing file is a clean no-op,
!!      * two-pass read (count -> allocate -> rewind -> fill),
!!      * title + header line skips, list-directed (free) reads,
!!      * allocation off sp_ob%hru and soil(ihru)%nly,
!!      * db_mx bookkeeping for the database count.
!!
!!    ~ ~ ~ INPUT FILE 1 : "pfas.dat"  (global PFAS database) ~ ~ ~
!!    Analogous to readpfas.f / pesticide.pes.  One title line, one header
!!    line, then one record per PFAS compound.  List-directed columns:
!!
!!      id   name        mw        sol       kl         lm        percop
!!      (-)  (a16)     kg/mol     mg/L      L/nmol    nmol/m^2     0-1
!!
!!      id      sequential PFAS index (1..npfas); a record with id<=0 ends
!!              the list (mirrors readpfas.f "if (ip == 0) exit").
!!      name    PFAS compound name                       -> pfasdb%name
!!      mw      molecular weight, kg/mol  (solver "a")   -> pfasdb%mw
!!      sol     max aqueous solubility, mg/L             -> pfasdb%sol
!!      kl      Langmuir K_L, L/nmol      (solver "n")   -> pfasdb%kl
!!      lm      Langmuir Gamma_max, nmol/m^2 (solver "h")-> pfasdb%lm
!!      percop  percolation/runoff partition coeff, 0-1  -> pfasdb%percop
!!
!!    ~ ~ ~ INPUT FILE 2 : "pfas_hru.ini"  (per-HRU soil pools) ~ ~ ~
!!    Analogous to the *.chm PFAS block of readchm.f, re-cast in the
!!    consolidated single-file style of pest_hru.ini / cs_hru.ini.  One
!!    record set per HRU, in HRU order (record k -> HRU k).  Layout:
!!
!!      line: title
!!      line: header
!!      then, repeated once per HRU (ihru = 1 .. nhru):
!!        * one header/name line  : "hru  <id>   nly  <nly>"  (skipped;
!!                                   nly is taken from soil(ihru)%nly)
!!        * num_pconta line       : integer, contaminated-site count
!!                                   (0 -> coerced to 1, as in readchm.f)
!!        * sol_d50 line          : real, median grain diameter (mm),
!!                                   one value per layer (1..nly)
!!        * then, repeated once per simulated PFAS (ip = 1 .. npfas):
!!            - pfasid + sol_pfas(1..nly)   ! kg/ha   layer PFAS mass
!!            - pfasid + kf(1..nly)         ! (nmol/kg)/(nM)^n Freundlich c
!!            - pfasid + nf(1..nly)         ! none    Freundlich exponent m
!!            - pfasid + enr                ! none    enrichment ratio
!!          where pfasid (1..npfas) selects the database compound; a
!!          record with pfasid<=0 marks "no PFAS in this HRU" and the
!!          remaining PFAS records for that HRU are skipped.
!!
!!    Soil-pool unit note (byte-faithful to readchm.f): the file carries
!!    initial soil PFAS as micrograms/ha; readchm.f scaled by 1.e-9 to
!!    kg/ha (sol_pfas = 1.e-9 * XX).  That same 1.e-9 conversion is applied
!!    here so the stored pfas_soil pools are kg/ha, matching the units the
!!    pfas_partition solver expects for "totmass".
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      use pfas_module
      use pfas_output_module
      use hydrograph_module, only : sp_ob
      use soil_module, only : soil
      use maximum_data_module
      use input_file_module

      implicit none

      character (len=80) :: titldum = ""   !!          |title line of file
      character (len=80) :: header = ""    !!          |header line of file
      integer :: eof = 0                   !!          |end-of-file flag
      integer :: imax = 0                  !! none     |number of records in database file
      integer :: ip = 0                    !! none     |database record counter
      integer :: id = 0                    !! none     |PFAS id read from a record
      integer :: ihru = 0                  !! none     |HRU counter
      integer :: nhru = 0                  !! none     |number of HRUs
      integer :: ly = 0                    !! none     |soil-layer counter
      integer :: nly = 0                   !! none     |number of layers in current HRU
      integer :: k = 0                     !! none     |PFAS slot counter
      integer :: pfasid = 0                !! none     |PFAS id read from a soil record
      integer :: nconta = 0                !! none     |contaminated-site count read from file
      integer :: rdly = 0                  !! none     |layer count declared in the file header line
      integer :: idd = 0                   !! none     |scratch HRU id from header line
      integer :: mly = 0                   !! none     |min(rdly, nly) layers actually stored
      character (len=20) :: c1 = ""        !! none     |scratch header token ("hru")
      character (len=20) :: c2 = ""        !! none     |scratch header token ("nly")
      logical :: i_exist = .false.         !! none     |file-existence flag
      real, dimension(:), allocatable :: xx  !! varies |per-layer scratch read buffer
      real :: enrtmp = 0.                  !! none     |enrichment scratch
      type (pfas_db) :: dbrec              !!          |scratch database record

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    (1) GLOBAL PFAS DATABASE  ->  pfasdb(:), npfas
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      eof = 0
      imax = 0
      npfas = 0

      inquire (file="pfas.dat", exist=i_exist)
      if (.not. i_exist) then
        !! no database -> PFAS feature inactive; nothing to allocate
        allocate (pfasdb(0:0))
        npfas = 0
        db_mx%pestparm = db_mx%pestparm   !! (no-op: keep pest count intact)
        return
      end if

      do
        open (170, file="pfas.dat")
        read (170,*,iostat=eof) titldum
        if (eof < 0) exit
        read (170,*,iostat=eof) header
        if (eof < 0) exit

        !! pass 1: count compound records (stop at id<=0 or EOF)
        do
          read (170,*,iostat=eof) id
          if (eof < 0) exit
          if (id <= 0) exit
          imax = imax + 1
        end do

        npfas = imax
        allocate (pfasdb(0:imax))

        !! sequential PFAS -> pfasdb index crosswalk (identity here, but
        !! the indirection mirrors readpfas.f nopfase/npfasno and lets the
        !! soil reader address compounds by id)
        if (allocated(pfas_num)) deallocate (pfas_num)
        allocate (pfas_num(0:imax))
        do ip = 0, imax
          pfas_num(ip) = ip
        end do

        !! pass 2: fill the database
        rewind (170)
        read (170,*,iostat=eof) titldum
        if (eof < 0) exit
        read (170,*,iostat=eof) header
        if (eof < 0) exit
        do ip = 1, imax
          read (170,*,iostat=eof) id, dbrec%name, dbrec%mw, dbrec%sol,     &
     &                            dbrec%kl, dbrec%lm, dbrec%percop
          if (eof < 0) exit
          if (id <= 0) exit
          pfasdb(id)%name   = dbrec%name
          pfasdb(id)%mw     = dbrec%mw
          pfasdb(id)%sol    = dbrec%sol
          pfasdb(id)%kl     = dbrec%kl
          pfasdb(id)%lm     = dbrec%lm
          pfasdb(id)%percop = dbrec%percop
        end do
        exit
      end do
      close (170)

      if (.not. allocated(pfasdb)) then
        allocate (pfasdb(0:0))
        npfas = 0
      end if

      !! nothing more to do if no PFAS compounds were defined
      if (npfas <= 0) return

      !! optional calibration multipliers: pfas_calib.dat = one line
      !! "soil_scale koc_scale" (global). Missing -> both stay 1.0 (no scaling).
      inquire (file="pfas_calib.dat", exist=i_exist)
      if (i_exist) then
        open (171, file="pfas_calib.dat")
        read (171,*,iostat=eof) pfas_soil_scale, pfas_koc_scale
        if (eof /= 0) then
          pfas_soil_scale = 1.0
          pfas_koc_scale  = 1.0
        end if
        close (171)
        if (pfas_soil_scale < 0.) pfas_soil_scale = 1.0
        if (pfas_koc_scale  < 0.) pfas_koc_scale  = 1.0
      end if

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    (2) ALLOCATE PER-HRU SOIL COLUMN  ->  pfas_soil_hru(:)%ly(:)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      nhru = sp_ob%hru
      if (nhru <= 0) return

      if (.not. allocated(pfas_soil_hru)) allocate (pfas_soil_hru(nhru))
      if (.not. allocated(pfas_flag))     allocate (pfas_flag(nhru))
      pfas_flag = 0

      !! per-HRU daily land-loss accumulators (written by pfas_lch / pfas_sed)
      if (.not. allocated(hpfasb_d)) allocate (hpfasb_d(nhru))

      do ihru = 1, nhru
        nly = soil(ihru)%nly
        if (nly < 1) nly = 1
        pfas_soil_hru(ihru)%num_pconta = 1
        allocate (pfas_soil_hru(ihru)%ly(nly))
        do ly = 1, nly
          allocate (pfas_soil_hru(ihru)%ly(ly)%sol_pfas(npfas), source=0.)
          allocate (pfas_soil_hru(ihru)%ly(ly)%kf(npfas),       source=0.)
          allocate (pfas_soil_hru(ihru)%ly(ly)%nf(npfas),       source=0.)
          allocate (pfas_soil_hru(ihru)%ly(ly)%enr(npfas),      source=0.)
          allocate (pfas_soil_hru(ihru)%ly(ly)%cw(npfas),       source=0.)
          pfas_soil_hru(ihru)%ly(ly)%sol_d50 = 0.
        end do
        allocate (hpfasb_d(ihru)%surq(npfas), source=0.)
        allocate (hpfasb_d(ihru)%latq(npfas), source=0.)
        allocate (hpfasb_d(ihru)%perc(npfas), source=0.)
        allocate (hpfasb_d(ihru)%sed(npfas),  source=0.)
        allocate (hpfasb_d(ihru)%surq_a(npfas), source=0.)
        allocate (hpfasb_d(ihru)%latq_a(npfas), source=0.)
        allocate (hpfasb_d(ihru)%perc_a(npfas), source=0.)
        allocate (hpfasb_d(ihru)%sed_a(npfas),  source=0.)
      end do

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    (3) READ PER-HRU SOIL POOLS  ->  pfas_soil_hru(:)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      eof = 0
      inquire (file="pfas_hru.ini", exist=i_exist)
      if (.not. i_exist) then
        !! pools stay at the zero initial condition allocated above
        return
      end if

      do
        open (171, file="pfas_hru.ini")
        read (171,*,iostat=eof) titldum            !! title
        if (eof < 0) exit
        read (171,*,iostat=eof) header             !! header
        if (eof < 0) exit

        do ihru = 1, nhru
          nly = soil(ihru)%nly
          if (nly < 1) nly = 1

          !! per-HRU header: "hru <id> nly <rdly>". We read the file's OWN declared
          !! layer count (rdly) and parse exactly rdly values per subsequent line, so
          !! the list-directed reads consume exactly what the writer emitted and can
          !! never drift across line boundaries even if rdly /= soil(ihru)%nly. Values
          !! are stored into the first mly = min(rdly, nly) pool layers.
          read (171,*,iostat=eof) c1, idd, c2, rdly
          if (eof < 0) exit
          if (rdly < 1) rdly = 1
          mly = min(rdly, nly)
          if (allocated(xx)) deallocate (xx)
          allocate (xx(max(rdly, nly)), source = 0.)

          !! contaminated-site count (0 -> 1, as in readchm.f)
          read (171,*,iostat=eof) nconta
          if (eof < 0) exit
          if (nconta <= 0) nconta = 1
          pfas_soil_hru(ihru)%num_pconta = nconta

          !! median grain diameter d50 (mm), one value per file layer
          read (171,*,iostat=eof) (xx(ly), ly = 1, rdly)
          if (eof < 0) exit
          do ly = 1, mly
            pfas_soil_hru(ihru)%ly(ly)%sol_d50 = xx(ly)
          end do

          !! per-PFAS blocks: sol_pfas, kf, nf, enr (in that order)
          do k = 1, npfas

            !! -- initial soil PFAS mass per layer (microg/ha -> kg/ha) --
            read (171,*,iostat=eof) pfasid, (xx(ly), ly = 1, rdly)
            if (eof < 0) exit
            if (pfasid <= 0) then
              !! "no PFAS in this HRU" sentinel: leave pools at zero and
              !! skip the rest of this HRU's PFAS records
              exit
            end if
            pfas_flag(ihru) = 1
            do ly = 1, mly
              pfas_soil_hru(ihru)%ly(ly)%sol_pfas(pfasid) =               &
     &                                  pfas_soil_scale * 1.e-9 * xx(ly)
            end do

            !! -- Freundlich coefficient kf (solver "c") per layer --
            read (171,*,iostat=eof) id, (xx(ly), ly = 1, rdly)
            if (eof < 0) exit
            do ly = 1, mly
              pfas_soil_hru(ihru)%ly(ly)%kf(pfasid) = xx(ly)
            end do

            !! -- Freundlich exponent nf (solver "m") per layer --
            read (171,*,iostat=eof) id, (xx(ly), ly = 1, rdly)
            if (eof < 0) exit
            do ly = 1, mly
              pfas_soil_hru(ihru)%ly(ly)%nf(pfasid) = xx(ly)
            end do

            !! -- enrichment ratio (single value, broadcast to layers) --
            read (171,*,iostat=eof) id, enrtmp
            if (eof < 0) exit
            do ly = 1, nly
              pfas_soil_hru(ihru)%ly(ly)%enr(pfasid) = enrtmp
            end do

          end do    !! PFAS block

          if (eof < 0) exit
        end do      !! HRU loop

        exit
      end do
      close (171)

      if (allocated(xx)) deallocate (xx)

      !! capture initial per-HRU soil PFAS mass (kg/ha, summed over layers) as the
      !! reference for the end-of-run mass-balance check (pfas_output)
      if (.not. allocated(pfas_init_hru)) allocate (pfas_init_hru(nhru, npfas), source=0.)
      do ihru = 1, nhru
        nly = soil(ihru)%nly
        if (nly < 1) nly = 1
        do k = 1, npfas
          do ly = 1, nly
            pfas_init_hru(ihru,k) = pfas_init_hru(ihru,k)                  &
     &                            + pfas_soil_hru(ihru)%ly(ly)%sol_pfas(k)
          end do
        end do
      end do

      return
      end subroutine pfas_read
