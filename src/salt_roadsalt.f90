      subroutine salt_roadsalt
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine adds salt from applied road salt to the soil profile

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      use basin_module
      use organic_mineral_mass_module
      use hydrograph_module
      use hru_module, only : hru, ihru
      use climate_module
      use output_landscape_module
      use salt_module
      use constituent_mass_module
      use time_module, only : time
      use salt_engine_module, only : road_timing_hru_gauge, rdw
      
      implicit none

      integer :: iadep !            |
      integer :: j !none        |counter
      integer :: iob !            |
      integer :: ist !            |
      integer :: isalt !            |salt ion counter
      integer :: iyr_road !       |road-salt year of this calendar year (road_timing hru_gauge)
      real :: const !            |constant used for rate, days, etc

      j = ihru
      
      if(cs_db%num_salts > 0) then
      
      iob = hru(j)%obj_no
      iwst = ob(iob)%wst
      iadep = wst(iwst)%wco%atmodep
      ist = atmodep_cont%ts
        
      !! salt_engine.ctl road_timing hru_gauge (E3): this calendar year's load at today's weight of the HRU's own gauges;
      !! the year comes from the calendar, not from atmodep_cont%ts (every simulated year was checked at read)
      if (road_timing_hru_gauge) then
        iyr_road = time%yrc - atmodep_cont%yr_init + 1
        do isalt = 1, cs_db%num_salts
          hsaltb_d(j)%salt(isalt)%road = rdapp_salt(iadep)%salt(isalt)%roadyr(iyr_road) *                     &
                                         rdw(iwst)%w(time%day, iyr_road)  !kg/ha
          cs_soil(j)%ly(1)%salt(isalt) = cs_soil(j)%ly(1)%salt(isalt) + hsaltb_d(j)%salt(isalt)%road
        enddo
        return
      end if

      !add applied road salt to the soil profile
      if (ist > 0 .and. ist <= atmodep_cont%num) then
      
        !monthly values
        if (atmodep_cont%timestep == "mo") then
          const = float (ndays(time%mo + 1) - ndays(time%mo))
          !loop through the salt ions
          do isalt=1,cs_db%num_salts 
            !add to salt balance arrays
            hsaltb_d(j)%salt(isalt)%road = rdapp_salt(iadep)%salt(isalt)%roadmo(ist) / const  !kg/ha
            !add both to soil profile  
            cs_soil(j)%ly(1)%salt(isalt) = cs_soil(j)%ly(1)%salt(isalt) + hsaltb_d(j)%salt(isalt)%road
          enddo
        endif 
        
        !yearly values
        if (atmodep_cont%timestep == "yr") then
          !loop through the salt ions
          do isalt=1,cs_db%num_salts 
            !add to salt balance arrays
            hsaltb_d(j)%salt(isalt)%road = rdapp_salt(iadep)%salt(isalt)%roadday(time%day,ist)  !kg/ha
            !add both to soil profile  
            cs_soil(j)%ly(1)%salt(isalt) = cs_soil(j)%ly(1)%salt(isalt) + hsaltb_d(j)%salt(isalt)%road
          enddo
        endif
        
      endif
      
      !annual average values
      if (atmodep_cont%timestep == "aa") then
        !loop through the salt ions
        do isalt=1,cs_db%num_salts 
          !add to salt balance arrays
          hsaltb_d(j)%salt(isalt)%road = rdapp_salt(iadep)%salt(isalt)%road / 365.  !kg/ha
          !add both to soil profile  
          cs_soil(j)%ly(1)%salt(isalt) = cs_soil(j)%ly(1)%salt(isalt) + hsaltb_d(j)%salt(isalt)%road
        enddo
      endif
       
      endif
      
      return
      end subroutine salt_roadsalt