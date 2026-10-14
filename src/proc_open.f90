      subroutine proc_open

      use basin_module
      use netcdf_output_module, only: netcdf_output_init
      use output_filter_module, only: output_filter_read
      implicit none
      
      external :: header_aquifer, header_channel, header_const, header_hyd, header_lu_change, header_mgt, &
                  header_path, header_pest, header_reservoir, header_salt, header_sd_channel, header_snutc, &
                  header_water_allocation, header_wetland, header_write, header_yield, &
                  output_landscape_init, search

      !! swatplus_perf: optional print_filter.prt (channel_sd channel list)
      call output_filter_read()

      !! swatplus_perf: output_filter
      call output_filter_read()
      
      !! write headers in output files (swatplus_perf: netcdf)
      if (pco%cdfout == "y") then
        call netcdf_output_init
        !! fort-leak-fix: suppress text-only outputs with no NetCDF backend (avoid fort.<N>)
        pco%mgtout = "n"
        pco%crop_yld = "n"
        pco%hydcon = "n"
        pco%fdcout = "n"
        !! basin and landscape-unit nb / ls / pw and the basin aquifer are NOT suppressed here: they have NetCDF streams, opened by
        !! netcdf_output_init just above, and every text write in basin_output, lsu_output and basin_aquifer_output sits in the
        !! else of an `if (pco%cdfout == "y")`, so no fort.* unit is written. Suppressing them after the streams were opened left
        !! those files full of zeros (the streams flushed every day, nothing was ever staged).
        pco%sd_chan_bsn%d = "n"; pco%sd_chan_bsn%m = "n"; pco%sd_chan_bsn%y = "n"; pco%sd_chan_bsn%a = "n"
        pco%res%d = "n"; pco%res%m = "n"; pco%res%y = "n"; pco%res%a = "n"
        pco%hyd%d = "n"; pco%hyd%m = "n"; pco%hyd%y = "n"; pco%hyd%a = "n"
        pco%ru%d = "n"; pco%ru%m = "n"; pco%ru%y = "n"; pco%ru%a = "n"
      else
      call output_landscape_init
      call header_channel
      call header_aquifer
      call header_sd_channel
      call header_mgt
      call header_lu_change
      call header_yield
      call header_hyd
      call header_reservoir
      call header_wetland
      ! call header_snutc
      call header_water_allocation
      
      call header_pest
      call header_path
      call header_salt !rtb salt
      call header_const !rtb cs

      call header_write
      end if
           
      return
      
      end subroutine proc_open