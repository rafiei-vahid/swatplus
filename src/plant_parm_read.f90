      subroutine plant_parm_read
      
      use input_file_module
      use maximum_data_module
      use plant_data_module
      use basin_module
      
      implicit none 
      
      external :: search
      integer :: ic !none       |plant counter
      character (len=80) :: titldum  !           |title of file
      character (len=80) :: header  !           |header of file
      integer :: eof !           |end of file
      integer :: imax !none       |determine max number for array (imax) and total number in file
      integer :: mpl !           | 
      logical :: i_exist              !none       |check to determine if file exists
      character (len=4000) :: plrec   !           |one plants.plt record, parsed on its own (Q1072)
      character (len=200) :: plmsg    !           |iomsg of a failed parse (Q1072)
      character (len=40) :: cdum(3)   !           |text columns skipped when reading pl_class (Q1072)
      real :: rdum(50)                !           |numeric columns skipped when reading pl_class (Q1072)
      real :: days_mat_r              !           |days_mat as the editor writes it, a real (Q1072)
      real :: mat_yrs_r               !           |yrs_mat as the editor writes it, a real (Q1072)
      integer :: ios                  !           |iostat of the record parse (Q1072)
      character (len=4000) :: plhdr   !           |the header record; its column count decides the layout (Q1072)
      integer :: plcols               !           |header columns: 54 (53 values + description) or 57 (+ 3 lignin) (Q1072)
      integer :: k                    !           |character / column counter (Q1072)
      real :: rlig(3)                 !           |lignin columns skipped when reading pl_class on the 57 layout (Q1072)

      
      eof = 0
      imax = 0
      mpl = 0

      inquire (file=in_parmdb%plants_plt, exist=i_exist)
      if (.not. i_exist .or. in_parmdb%plants_plt == " null") then
        allocate (pldb(0:0))
        allocate (plcp(0:0))
        allocate (pl_class(0:0))
        if (bsn_cc%cswat == 2) allocate (cswat_1_part_fracs(0:0))
      else
      do
        open (104,file=in_parmdb%plants_plt)
        read (104,*,iostat=eof) titldum
        if (eof < 0) exit
        read (104,*,iostat=eof) header
        if (eof < 0) exit
          do while (eof == 0)
            read (104,*,iostat=eof) titldum
            if (eof < 0) exit
            imax = imax + 1
          end do
        allocate (pldb(0:imax))
        allocate (plcp(0:imax))
        allocate (pl_class(0:imax))
        if (bsn_cc%cswat == 2) allocate (cswat_1_part_fracs(0:imax))
        
        rewind (104)
        read (104,*,iostat=eof) titldum
        if (eof < 0) exit
        read (104,'(a)',iostat=eof,iomsg=plmsg) plhdr
        if (eof < 0) exit
        if (eof > 0) then
          write (*,'(a,i6,2a)') ' plants.plt header: record read iostat ', eof, ' -- ', trim(plmsg)
          error stop 'plants.plt: the header could not be read (Q1072)'
        end if
        if (plhdr(len(plhdr):len(plhdr)) /= ' ') error stop 'plants.plt: the header is longer than the read buffer (Q1072)'
        !! Q1072: the header decides the layout. 54 columns = 53 values + description (the editor); 57 = 53 values +
        !! avg_lig_frac, ab_lig_frac, bg_lig_frac + description (upstream's refdata/Ames_sub1, carbon = 2), whose three
        !! lignin columns the old whole-type read took into res_part_fracs. Any other count is a layout this reader does
        !! not know, so it stops rather than guess.
        plcols = 0
        do k = 1, len_trim(plhdr)
          if (plhdr(k:k) /= ' ' .and. plhdr(k:k) /= char(9)) then
            if (k == 1) then
              plcols = plcols + 1
            else if (plhdr(k-1:k-1) == ' ' .or. plhdr(k-1:k-1) == char(9)) then
              plcols = plcols + 1
            end if
          end if
        end do
        if (plcols /= 54 .and. plcols /= 57) then
          write (*,'(a,i5,a)') ' plants.plt header has ', plcols, &
            ' columns; this reader knows 54 (53 values + description) and 57 (+ avg/ab/bg_lig_frac) (Q1072)'
          error stop 'plants.plt: unknown column layout (Q1072)'
        end if
        
        do ic = 1, imax
          !! Q1072: read the row as ONE record, then parse its 53 values from that string. The old list-directed read of
          !! the whole plant_db stopped at item 5 under gfortran, because the editor writes days_mat and yrs_mat as reals
          !! ("110.00000"), and it silently kept every later item at its type default (iostat 5010, only eof < 0 was
          !! checked). It could never fit the file on any compiler either: the type's 54th item is res_part_fracs, while the
          !! file's 54th column is the description. Parsing a single record also stops a short or malformed row from
          !! running on into the next one. Any failure now STOPS the run and names the row. On the 57-column layout the three
          !! lignin columns (54-56) are read into res_part_fracs as the old read did; on 54 they keep the type defaults.
          read (104,'(a)',iostat=eof,iomsg=plmsg) plrec
          if (eof < 0) exit
          if (eof > 0) then     !! an I/O error on the record itself must stop too, never reuse the previous record (Lane S)
            write (*,'(a,i5,a,i6,2a)') ' plants.plt row ', ic, ': record read iostat ', eof, ' -- ', trim(plmsg)
            error stop 'plants.plt: a record could not be read (Q1072)'
          end if
          if (plrec(len(plrec):len(plrec)) /= ' ') then
            write (*,'(a,i5,a,i6,a)') ' plants.plt row ', ic, ' is longer than ', len(plrec), ' characters (Q1072)'
            error stop 'plants.plt: a record is longer than the read buffer'
          end if
          read (plrec,*,iostat=ios,iomsg=plmsg) pldb(ic)%plantnm, pldb(ic)%typ, pldb(ic)%trig, pldb(ic)%nfix_co,      &
            days_mat_r, pldb(ic)%bio_e, pldb(ic)%hvsti, pldb(ic)%blai, pldb(ic)%frgrw1, pldb(ic)%laimx1,                  &
            pldb(ic)%frgrw2, pldb(ic)%laimx2, pldb(ic)%dlai, pldb(ic)%dlai_rate, pldb(ic)%chtmx, pldb(ic)%rdmx,           &
            pldb(ic)%t_opt, pldb(ic)%t_base, pldb(ic)%cnyld, pldb(ic)%cpyld, pldb(ic)%pltnfr1, pldb(ic)%pltnfr2,          &
            pldb(ic)%pltnfr3, pldb(ic)%pltpfr1, pldb(ic)%pltpfr2, pldb(ic)%pltpfr3, pldb(ic)%wsyf, pldb(ic)%usle_c,       &
            pldb(ic)%gsi, pldb(ic)%vpdfr, pldb(ic)%gmaxfr, pldb(ic)%wavp, pldb(ic)%co2hi, pldb(ic)%bioehi,                &
            pldb(ic)%rsdco_pl, pldb(ic)%alai_min, pldb(ic)%laixco_tree, mat_yrs_r, pldb(ic)%bmx_peren,                    &
            pldb(ic)%ext_coef, pldb(ic)%leaf_tov_min, pldb(ic)%leaf_tov_max, pldb(ic)%bm_dieoff, pldb(ic)%rsr1,           &
            pldb(ic)%rsr2, pldb(ic)%pop1, pldb(ic)%frlai1, pldb(ic)%pop2, pldb(ic)%frlai2, pldb(ic)%frsw_gro,             &
            pldb(ic)%aeration, pldb(ic)%rsd_pctcov, pldb(ic)%rsd_covfac
          if (ios /= 0) then
            write (*,'(a,i5,a,i6,2a)') ' plants.plt row ', ic, ': iostat ', ios, ' -- ', trim(plmsg)
            write (*,'(2a)') ' record: ', trim(plrec)
            error stop 'plants.plt: a row could not be read (Q1072)'
          end if
          pldb(ic)%days_mat = nint(days_mat_r)
          pldb(ic)%mat_yrs = nint(mat_yrs_r)
          if (plcols == 57) then
            read (plrec,*,iostat=ios,iomsg=plmsg) cdum, rdum, pldb(ic)%res_part_fracs%meta_frac,                       &
              pldb(ic)%res_part_fracs%str_frac, pldb(ic)%res_part_fracs%lig_frac
            if (ios /= 0) then
              write (*,'(a,i5,a,i6,2a)') ' plants.plt row ', ic, ' lignin columns 54-56: iostat ', ios, ' -- ', trim(plmsg)
              write (*,'(2a)') ' record: ', trim(plrec)
              error stop 'plants.plt: the lignin columns (54-56) could not be read (Q1072)'
            end if
          end if
          if (bsn_cc%nam1 /= 0) then
            !! the plant class is the description column, 54 or 57: the token the old read gave pl_class on 57 (Q1072)
            read (plrec,*,iostat=ios,iomsg=plmsg) cdum, rdum, (rlig(k), k = 1, plcols - 54), pl_class(ic)
            if (ios /= 0) then
              write (*,'(a,i5,a,i6,2a)') ' plants.plt row ', ic, ' plant class: iostat ', ios, ' -- ', trim(plmsg)
              error stop 'plants.plt: the plant class column could not be read (Q1072)'
            end if
          end if
          pldb(ic)%mat_yrs = Max (1, pldb(ic)%mat_yrs)
          if (bsn_cc%cswat == 2) then
            cswat_1_part_fracs(ic)%lig_frac_blg = pldb(ic)%res_part_fracs%lig_frac
            cswat_1_part_fracs(ic)%lig_frac_abg = pldb(ic)%res_part_fracs%str_frac
            cswat_1_part_fracs(ic)%str_frac_blg = cswat_1_part_fracs(ic)%lig_frac_blg / .80 
            cswat_1_part_fracs(ic)%str_frac_abg = cswat_1_part_fracs(ic)%lig_frac_abg / .80 
            cswat_1_part_fracs(ic)%meta_frac_blg = 1.0 - cswat_1_part_fracs(ic)%str_frac_blg 
            cswat_1_part_fracs(ic)%meta_frac_abg = 1.0 - cswat_1_part_fracs(ic)%str_frac_abg  

          else
            pldb(ic)%res_part_fracs%meta_frac = 0.85
            pldb(ic)%res_part_fracs%str_frac = 0.15 
            pldb(ic)%res_part_fracs%lig_frac = 0.12
          endif
            
              
        end do
        
        exit
      enddo
      endif

      db_mx%plantparm = imax
      
      close (104)
      return

      end subroutine plant_parm_read
