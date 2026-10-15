      subroutine salt_roadsalt_read
      !read in road salt loadings (kg/ha)
  
      use basin_module
      use input_file_module
      use climate_module
      use time_module
      use maximum_data_module
      use constituent_mass_module
      use salt_engine_module
!$    use omp_lib
      
      implicit none
      
      character (len=4) :: salt_ion  !           |
      character (len=15) :: station_name  !       |
      integer :: eof !           |end of file
      integer :: iadep !           |counter
      integer :: imo !           |counter
      integer :: iyr !           |counter
      logical :: i_exist              !none       |check to determine if file exists
      integer :: isalt !           |salt ion counter
      integer :: year_index
      integer :: yr_weat
      integer :: year_days
      integer :: iday
      integer :: ie2
      integer :: day_flag(366)
      real :: day_precip
      real :: day_temp
      real :: year_precip
      real :: day_fraction
      integer :: n_omp                !none       |OpenMP threads this run will use (1 without OpenMP)
           
      eof = 0

      !! optional salt-engine switches: salt_engine.ctl is read only when it exists (salt_engine_module)
      call salt_engine_ctl_read
      if (road_timing_hru_gauge) then
        inquire (file='salt_road',exist=i_exist)
        if (cs_db%num_salts == 0 .or. .not. i_exist) then
          write (*,"(a)") "salt_engine.ctl road_timing hru_gauge: the deck has no salt ions or no salt_road file"
          error stop "road_timing hru_gauge needs salt ions and a salt_road file"
        end if
        if (atmodep_cont%timestep /= "yr") then
          write (*,"(3a)") "salt_engine.ctl road_timing hru_gauge: the deposition timestep is '",                &
                           trim(atmodep_cont%timestep), "'; it applies to the yearly form 'yr' only"
          error stop "road_timing hru_gauge needs the yearly road-salt form"
        end if
      end if
      if (road_load_hru_file .and. .not. road_timing_hru_gauge) then
        write (*,"(a)") "salt_engine.ctl road_load hru_file: it takes E3's timing, so it needs road_timing hru_gauge too"
        error stop "road_load hru_file needs road_timing hru_gauge"
      end if
      if (road_washoff_scale_lines > 0 .and. .not. road_washoff_imp) then
        write (*,"(a)") "salt_engine.ctl road_washoff_scale: it scales road_washoff impervious, so it needs that key too"
        error stop "road_washoff_scale needs road_washoff impervious"
      end if
      if (road_washoff_imp) then
        if (.not. road_timing_hru_gauge) then
          write (*,"(a)") "salt_engine.ctl road_washoff impervious: it splits E3's daily load, so it needs road_timing hru_gauge too"
          error stop "road_washoff impervious needs road_timing hru_gauge"
        end if
        if (time%step > 1) then
          write (*,"(a)") "salt_engine.ctl road_washoff impervious: daily runs only (hru_urbanhr has no impervious store)"
          error stop "road_washoff impervious needs a daily time step"
        end if
        if (cs_db%num_salts < 1 .or. cs_db%num_salts > 8) then
          write (*,"(a,i0)") "salt_engine.ctl road_washoff impervious: salt ions must number 1 to 8, the deck has ", cs_db%num_salts
          error stop "road_washoff impervious needs 1 to 8 salt ions"
        end if
        n_e2 = cs_db%num_salts
        n_e2b = 2
        do ie2 = 1, n_e2
          e2_hdr_hru(ie2) = "  imps_" // trim(cs_db%salts(ie2))
          e2_hdr_hru(n_e2 + ie2) = "  impw_" // trim(cs_db%salts(ie2))
        end do
      end if
      if (aqu_hold_c0) then
        if (bsn_cc%gwflow == 1) then
          write (*,"(a)") "salt_engine.ctl aquifer_salt hold_initial: the gwflow aquifer is a different module, not covered"
          error stop "aquifer_salt hold_initial: gwflow is not covered"
        end if
        if (cs_db%num_salts < 1 .or. cs_db%num_salts > 8) then
          write (*,"(a,i0)") "salt_engine.ctl aquifer_salt hold_initial: salt ions must number 1 to 8, the deck has ", cs_db%num_salts
          error stop "aquifer_salt hold_initial needs 1 to 8 salt ions"
        end if
        n_e5 = cs_db%num_salts
        n_e5b = 2
        do ie2 = 1, n_e5
          e5_hdr_aqu(ie2) = "  rchout_" // trim(cs_db%salts(ie2))
          e5_hdr_aqu(n_e5 + ie2) = "  hold_" // trim(cs_db%salts(ie2))
        end do
      end if
      if (salt_chemistry_off .and. cs_db%num_salts == 0) then
        write (*,"(a)") "salt_engine.ctl salt_chemistry off: the deck has no salt ions"
        error stop "salt_chemistry off needs salt ions"
      end if
      !! E-TS: the salt chemistry keeps its scratch arrays and index counters in salt_data_module (and the saved Sol_*_p arrays
      !! of salt_chem_hru), shared by every OpenMP thread, so a run with salt ions, the chemistry on and more than one thread stops
      !! here instead of racing (it crashed on 3 of 3 runs at 4 threads). Without OpenMP the sentinel line compiles out: 1 thread.
      n_omp = 1
!$    n_omp = omp_get_max_threads()
      if (cs_db%num_salts > 0 .and. .not. salt_chemistry_off .and. n_omp > 1) then
        write (*,"(a,i0,a)") "salt chemistry is on with ", n_omp, " OpenMP threads, and it is not thread-safe: " // &
          "run with OMP_NUM_THREADS=1, or set 'salt_chemistry off' in salt_engine.ctl"
        error stop "salt chemistry needs OMP_NUM_THREADS=1"
      end if

      !only proceed if there are salt ions in the simulation
      if(cs_db%num_salts > 0) then
      
      inquire (file='salt_road',exist=i_exist)
      if(i_exist) then
        
        !open the file; skip first two lines (commentary)
        open(5051,file='salt_road')
        read(5051,*)
        read(5051,*)
        read(5051,*)
        read(5051,*)
      
        !allocate arrays
        allocate (rdapp_salt(0:atmodep_cont%num_sta))

        !! road_timing hru_gauge: the daily weights of each weather station's own gauges, built before the stations' loads
        if (road_timing_hru_gauge) call salt_road_timing_build

        !loop through the stations (num_sta is set in cli_read_atmodep subroutine)
        do iadep = 1, atmodep_cont%num_sta
          
          !allocate arrays
          allocate (rdapp_salt(iadep)%salt(cs_db%num_salts))
          
          !average annual values
          if (atmodep_cont%timestep == "aa") then
            read(5051,*) station_name !station name --> already read in cli_read_atmodep
            do isalt=1,cs_db%num_salts
              read(5051,*) salt_ion,rdapp_salt(iadep)%salt(isalt)%road
            enddo  
          endif
          
          !monthly values
          if (atmodep_cont%timestep == "mo") then
            read(5051,*) station_name !station name
            do isalt=1,cs_db%num_salts
              allocate (rdapp_salt(iadep)%salt(isalt)%roadmo(atmodep_cont%num), source = 0.)
              read(5051,*) salt_ion,(rdapp_salt(iadep)%salt(isalt)%roadmo(imo),imo=1,atmodep_cont%num)
            enddo 
          end if
                              
          !yearly values
          if (atmodep_cont%timestep == "yr") then
            read(5051,*) station_name !station name
            do isalt=1,cs_db%num_salts
              allocate (rdapp_salt(iadep)%salt(isalt)%roadyr(atmodep_cont%num), source = 0.)
              read(5051,*) salt_ion,(rdapp_salt(iadep)%salt(isalt)%roadyr(iyr),iyr=1,atmodep_cont%num)
            enddo 
            !loop through the years, partitioning annual salt load to daily salt loadings
            !(based on precipitation and temperature)
            do isalt=1,cs_db%num_salts
              allocate (rdapp_salt(iadep)%salt(isalt)%roadday(366,atmodep_cont%num), source = 0.)
              rdapp_salt(iadep)%salt(isalt)%roadday = 0.
                                    enddo
            !! with road_timing hru_gauge the applier reads salt_engine_module's weights, never roadday, so this partition
            !! (the deposition index read as a gauge index) is not computed
            if (.not. road_timing_hru_gauge) then
            year_index = atmodep_cont%yr_init
            yr_weat = 1 !weather year index
            do iyr=1,atmodep_cont%num
              !number of days in the year
              if(mod(year_index,4) == 0) then
                year_days = 366
              else
                year_days = 365
              endif
              if(year_index .ge. time%yrc_start) then
                year_precip = 0.
                day_flag = 0
                !first, determine days of snow, and total precipitation on days of snow
                do iday=1,year_days
                  !determine if there is snow (likely) 
                  day_precip = pcp(iadep)%ts(iday,yr_weat) !mm
                  day_temp = (tmp(iadep)%ts(iday,yr_weat) + tmp(iadep)%ts2(iday,yr_weat)) / 2 !average daily temp (oC)
                  if(day_temp.ne.-99) then
                    if(day_precip > 0 .and. day_temp < 0) then
                      day_flag(iday) = 1
                      year_precip = year_precip + day_precip !mm
                    endif
                  endif
                enddo
                !second, calculate daily salt loading (using  fractions (fraction of yearly snow that falls on each day)
                do iday=1,year_days
                  if(day_flag(iday) == 1) then
                    day_fraction = pcp(iadep)%ts(iday,yr_weat) / year_precip
                    do isalt=1,cs_db%num_salts
                      rdapp_salt(iadep)%salt(isalt)%roadday(iday,iyr) = rdapp_salt(iadep)%salt(isalt)%roadyr(iyr) * day_fraction !kg/ha
                    enddo
                  endif
                enddo
                yr_weat = yr_weat + 1
              else !no precipitation data --> set daily values to 0
                do iday=1,year_days
                  do isalt=1,cs_db%num_salts
                    rdapp_salt(iadep)%salt(isalt)%roadday(iday,iyr) = 0.
                  enddo
                enddo
              endif
              year_index = year_index + 1
            enddo 
            end if
          end if
          
        end do !go to next station
      
      endif
      endif
           
      return
      end subroutine salt_roadsalt_read
