      subroutine command
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    for every day of simulation, this subroutine steps through the command
!!    lines in the watershed configuration (.fig) file. Depending on the 
!!    command code on the .fig file line, a command loop is accessed
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    SWAT: subbasin, route, routres, transfer, recmon
!!    SWAT: recepic, save, recday, recyear

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use time_module
      use hydrograph_module
      use ru_module
      use channel_module
      use hru_lte_module
      use aquifer_module
      use sd_channel_module
      use reservoir_module
      use organic_mineral_mass_module
      use constituent_mass_module
      use hru_module, only : ihru, hru
      use basin_module
      use netcdf_output_module
      use maximum_data_module
      use gwflow_module
      use soil_module
      use recall_module
      use water_allocation_module
      use command_wave_module
      implicit none
      
      external :: aqu_1d_control, aqu_cs_output, aqu_pesticide_output, aqu_salt_output, aquifer_output, &
                  ch_cs_output, ch_salt_output, cha_pesticide_output, channel_output, constit_hyd_mult, &
                  cs_str_output, flow_dur_curve, gwflow_simulate, hru_carbon_output, hru_control, &
                  hru_cs_output, hru_lte_control, hru_lte_output, hru_output, hru_pathogen_output, &
                  hru_pesticide_output, hru_salt_output, hydin_output, hydout_output, manure_demand_output, &
                  manure_source_output, obj_output, recall_nut, recall_output, res_control, res_cs_output, &
                  res_pesticide_output, res_salt_output, reservoir_output, ru_control, ru_cs_output, &
                  ru_output, ru_salt_output, sd_chanbud_output, sd_chanmorph_output, sd_channel_control3, &
                  sd_channel_output, wallo_allo_output, wallo_treat_output, wallo_trn_output, &
                  wallo_use_output, wet_cs_output, wet_salt_output, wetland_output, basin_aqu_pest_output, &
                  basin_aquifer_output, basin_ch_pest_output, basin_chanbud_output, basin_chanmorph_output, &
                  basin_channel_output, basin_ls_pest_output, basin_output, basin_recall_output, &
                  basin_res_pest_output, basin_reservoir_output, basin_sdchannel_output, cs_balance, &
                  lsu_output, salt_balance, hyddep_output, recall_salt, recall_cs, soil_nutcarb_write, &
                  soil_carbvar_write

      real, dimension(time%step) :: hyd_flo     !flow hydrograph
      integer :: in = 0               !              | 
      integer :: iob = 0              !              |
      integer :: iday = 0             !              |
      integer :: isd = 0              !none          |counter
      integer :: ires = 0             !none          |reservoir number
      integer :: irec = 0             !              |
      integer :: iout = 0             !none          |counter
      integer :: ihtyp = 0            !              |
      integer :: iaq = 0              !none          |counter
      integer :: j = 0                !none          |counter
      integer :: ihyd = 0             !              |
      integer :: idr = 0              !              |
      integer :: iwro = 0             !              |
      real :: conv = 0.               !              |
      real :: frac_in = 0.            !              |
      integer :: ts1 = 0
      integer :: ts2 = 0
      integer :: iw = 0               !              |counter for water allocation object
      integer :: iwallo = 0           !              |variable to pass to wallo_control
      integer :: i_count = 0          !rtb gwflow
      integer :: i_mfl = 0            !rtb gwflow    |counter
      integer :: i_chan = 0           !rtb gwflow    |counter
      integer :: iob_chan = 0        !rtb gwflow    |ob index for channel
      real :: sumflo = 0.
      integer :: ic_walk = 0          !swatplus_perf: serial-walk index (module icmd is set inside command_object)
      integer :: lev = 0              !swatplus_perf: wavefront level
      integer :: k = 0                !swatplus_perf: index within a wave
      logical :: use_wave = .false.   !swatplus_perf: drive HRU land phase wave-by-wave
      external :: command_object

      icmd = sp_ob1%objs
      wallo(:)%trn_cur = 1
      if (allocated(res_ob)) res_ob(:)%wallo_call = 0

      !! swatplus_perf: per-object body moved to reentrant command_object(ic).
      !! Walk drives it serially here; the parallel HRU wavefront (added next)
      !! calls the same routine. ob(icmd)%cmd_next uses the module icmd, which
      !! command_object may advance for gwflow sub-objects (preserves old behavior).
      !! swatplus_perf: build the HRU wavefront index once (topological levels).
      if (.not. hru_wave_ready) call command_wave_build

      !! swatplus_perf wavefront: run the HRU land phase wave-by-wave (objects at the
      !! same cmd_order level are mutually independent), then walk the remaining
      !! (routing/channel/reservoir) objects serially, skipping the already-done HRUs.
      !! Fall back to a pure serial walk when water/manure allocation is active, since
      !! allocation is interleaved global state that the wave pre-pass would reorder.
      use_wave = (hru_wave_ready .and. hru_nwave > 0 .and.                  &
                  db_mx%wallo_db == 0 .and. db_mx%mallo_db == 0)

      if (use_wave) then
        do lev = 1, hru_nwave
          do k = 1, hru_wave_cnt(lev)
            call command_object (hru_wave_obj(lev, k))
          end do
        end do
      end if

      ic_walk = sp_ob1%objs
      do while (ic_walk /= 0)
        if (use_wave .and. ob(ic_walk)%typ == "hru") then
          ic_walk = ob(ic_walk)%cmd_next          ! HRU already done in the wave pre-pass
        else
          call command_object (ic_walk)
          ic_walk = ob(icmd)%cmd_next             ! module icmd may be gwflow-advanced
        end if
      end do

      !! write object output for entire simulation (fort-leak-fix: no NetCDF backend)
      if (pco%cdfout /= "y") call obj_output
      
      !! print all output files
      if (time%yrs > pco%nyskip) then
      
        !! print water allocation output
        do iwro =1, db_mx%wallo_db
          call wallo_allo_output (iwro)
          call wallo_trn_output (iwro)
          call wallo_treat_output (iwro)
          call wallo_use_output (iwro)
          !call wallo_osrc_output (iwro)
          !call wallo_odmd_output (iwro)
        end do
        
        !! print manure allocation output
        do iwro =1, db_mx%mallo_db
          call manure_source_output (iwro)
          call manure_demand_output (iwro)
        end do
        
        do isd = 1, sp_ob%hru_lte
          call hru_lte_output (isd)
        end do
        
        do ihru = 1, sp_ob%hru
          call hru_output (ihru)
          call hru_carbon_output (ihru)
          if (hru(ihru)%dbs%surf_stor > 0) then
            call wetland_output(ihru)
            if (cs_db%num_salts > 0) then !rtb salt
              call wet_salt_output(ihru)
            endif
            if (cs_db%num_cs > 0) then !rtb cs
              call wet_cs_output(ihru)
            endif
          end if
          if (cs_db%num_tot > 0) then 
            call hru_pesticide_output (ihru)
            call hru_pathogen_output (ihru)
          end if
          if (cs_db%num_salts > 0) then !rtb salt
            call hru_salt_output(ihru)
          endif
          if (cs_db%num_cs > 0) then !rtb cs
            call hru_cs_output(ihru)
          endif
          !sum annual for SWIFT input
          if (bsn_cc%swift_out == 1) then
            icmd = hru(ihru)%obj_no
            do ihyd = 1, 5
              ob(icmd)%hd_aa(ihyd) = ob(icmd)%hd_aa(ihyd) + ob(icmd)%hd(ihyd)
            end do
          end if
                         
          ! Call soil_nutcarb_write for specified output for hru_cb in print.prt
          if (pco%cb_hru%d == "y") call soil_nutcarb_write(" d")
          if (pco%cb_hru%d == "l") call soil_nutcarb_write("dl")
          if (pco%cb_hru%m == "y" .and. time%end_mo == 1) call soil_nutcarb_write(" m")
          if (pco%cb_hru%m == "l" .and. time%end_mo == 1) call soil_nutcarb_write("ml")
          if (pco%cb_hru%y == "y" .and. time%end_yr == 1) call soil_nutcarb_write(" y") 
          if (pco%cb_hru%y == "l" .and. time%end_yr == 1) call soil_nutcarb_write("yl") 

          ! Call soil_carbvar_write for specified output for hru_cb_vars in print.prt
          if (bsn_cc%cswat == 1) then
            if (pco%cb_vars_hru%d == "y") call soil_carbvar_write(" d")
            if (pco%cb_vars_hru%d == "l") call soil_carbvar_write("dl")
            if (pco%cb_vars_hru%m == "y" .and. time%end_mo == 1) call soil_carbvar_write(" m")
            if (pco%cb_vars_hru%m == "l" .and. time%end_mo == 1) call soil_carbvar_write("ml")
            if (pco%cb_vars_hru%y == "y" .and. time%end_yr == 1) call soil_carbvar_write(" y")
            if (pco%cb_vars_hru%y == "l" .and. time%end_yr == 1) call soil_carbvar_write("yl")
          endif
        
        end do      ! hru loop  
        !! swatplus_perf: netcdf (hru flush)
        if (pco%cdfout == "y") then
          call nc_flush_daily_hru()
          if (time%end_mo == 1) call nc_flush_monthly_hru()
          if (time%end_yr == 1) call nc_flush_yearly_hru()
          if (time%end_sim == 1) call nc_flush_aa_hru()
        end if
        
        do iaq = 1, sp_ob%aqu
          call aquifer_output (iaq)
          if (cs_db%num_salts > 0) then !rtb salt
            call aqu_salt_output (iaq)
          endif
          if (cs_db%num_cs > 0) then !rtb cs
            call aqu_cs_output(iaq)
          endif  
          if (cs_db%num_tot > 0) then 
            call aqu_pesticide_output (iaq)
          end if       
        end do
        
        do jrch = 1, sp_ob%chan
          call channel_output (jrch)
        end do
                
        do jrch = 1, sp_ob%chandeg
          call sd_chanmorph_output (jrch)
          call sd_chanbud_output (jrch)
          call sd_channel_output (jrch)
          if (cs_db%num_tot > 0) then 
            call cha_pesticide_output (jrch)   
            !call ch_pathogen_output (jrch)
          end if   
          if (cs_db%num_salts > 0) then !rtb salt
            call ch_salt_output (jrch)
          endif
          if (cs_db%num_cs > 0) then
            call ch_cs_output (jrch) !rtb cs
          endif
        end do
        if(cs_db%num_cs > 0) then
          call cs_str_output !rtb cs
        endif
        

        do j = 1, sp_ob%res
          call reservoir_output(j)
         if (cs_db%num_tot > 0) then 
            call res_pesticide_output (j)
            if (cs_db%num_salts > 0) then !rtb salt
              call res_salt_output (j)
            endif
            if (cs_db%num_cs > 0) then !rtb cs
              call res_cs_output (j)
            endif
            !call res_pathogen_output (j)
          end if       
        end do 
        
        do j = 1, sp_ob%ru
          call ru_output(j)
          if(cs_db%num_salts > 0) then !rtb salt
            call ru_salt_output(j)
          endif
          if(cs_db%num_cs > 0) then !rtb cs
            call ru_cs_output(j)
          endif
        end do
        
        do j = 1, sp_ob%recall
          call recall_output (j)
        end do

        call hydin_output   !if all output is no, then don"t call
        !call hcsin_output  gives allocate error
        if (sp_ob%chandeg > 0 .and. cs_db%num_pests > 0) call basin_ch_pest_output  
        if (sp_ob%res > 0 .and. cs_db%num_pests > 0) call basin_res_pest_output     
        if (sp_ob%hru > 0 .and. cs_db%num_pests > 0) call basin_ls_pest_output
        if (sp_ob%aqu > 0 .and. cs_db%num_pests > 0) call basin_aqu_pest_output
        if (db_mx%lsu_elem > 0) call basin_output
        if (db_mx%lsu_out > 0) call lsu_output
        if (db_mx%aqu_elem > 0) call basin_aquifer_output
        !if (sp_ob%aqu > 0) call basin_aquifer_output !rtb - otherwise, aquifer output is not called
        if (sp_ob%res > 0) call basin_reservoir_output
        if (sp_ob%chan > 0) call basin_channel_output
        if (sp_ob%chandeg > 0) call basin_chanmorph_output
        if (sp_ob%chandeg > 0) call basin_chanbud_output
        if (sp_ob%chandeg > 0) call basin_sdchannel_output
        if (sp_ob%recall > 0) call basin_recall_output
        !call soil_nutcarb_output
        !call lsreg_output
        !call region_aquifer_output
        !call region_reservoir_output
        !call region_channel_output
        !call region_recall_output
        
        if(cs_db%num_salts > 0) call salt_balance !rtb salt
        if(cs_db%num_cs > 0) call cs_balance !rtb cs
        
        !! swatplus_perf: netcdf (non-hru flush)
        if (pco%cdfout == "y") then
          call nc_flush_daily_basin()
          call nc_flush_daily_lsu()
          call nc_flush_daily_aqu()
          call nc_flush_daily_sd()
          call nc_flush_daily_channel()
          if (time%end_mo == 1) then
            call nc_flush_monthly_basin()
            call nc_flush_monthly_lsu()
            call nc_flush_monthly_aqu()
            call nc_flush_monthly_sd()
            call nc_flush_monthly_channel()
          end if
          if (time%end_yr == 1) then
            call nc_flush_yearly_basin()
            call nc_flush_yearly_lsu()
            call nc_flush_yearly_aqu()
            call nc_flush_yearly_sd()
            call nc_flush_yearly_channel()
          end if
          if (time%end_sim == 1) then
            call nc_flush_aa_basin()
            call nc_flush_aa_lsu()
            call nc_flush_aa_aqu()
            call nc_flush_aa_sd()
            call nc_flush_aa_channel()
          end if
        end if
      end if

      gw_daycount = gw_daycount + 1
      
      !rtb hydrograph separation
      !write out hydrograph components for all channels
      if (bsn_cc%gwflow == 1) then
      do i_chan=1,sp_ob%chandeg
        if(hydsep_flag(i_chan) == 1) then
          iob_chan = sp_ob1%chandeg + i_chan - 1
          write(out_hyd_sep,8102) time%day,time%mo,time%day_mo,time%yrc, &
            i_chan,ob(iob_chan)%gis_id,ob(iob_chan)%name, &
            (hyd_sep_array(i_chan,i_count),i_count=1,7)
        endif
      enddo
      endif
      !zero out arrays for next day
      icmd = sp_ob1%objs
      do while (icmd /= 0)
        ob(icmd)%hdsep%flo_surq = 0.
        ob(icmd)%hdsep%flo_latq = 0.
        ob(icmd)%hdsep%flo_gwsw = 0.
        ob(icmd)%hdsep%flo_swgw = 0.
        ob(icmd)%hdsep%flo_satex = 0.
        ob(icmd)%hdsep%flo_satexsw = 0.
        ob(icmd)%hdsep%flo_tile = 0.
        ob(icmd)%hdsep_in%flo_surq = 0.
        ob(icmd)%hdsep_in%flo_latq = 0.
        ob(icmd)%hdsep_in%flo_gwsw = 0.
        ob(icmd)%hdsep_in%flo_swgw = 0.
        ob(icmd)%hdsep_in%flo_satex = 0.
        ob(icmd)%hdsep_in%flo_satexsw = 0.
        ob(icmd)%hdsep_in%flo_tile = 0.  
        icmd = ob(icmd)%cmd_next
      enddo
      
102   format(i6,11x,i3,8x,i5,5x,1000(f16.4))
103   format(4i6,2i8,2x,a,35f12.3)
8102  format(4i6,2i8,a18,7e13.4)      

      return
      end subroutine command