       subroutine co2_read
      
       use input_file_module
       use basin_module
       use time_module
       use climate_module
       use output_path_module
       
       implicit none
       
       character (len=80) :: titldum  !             |title of file
       character (len=80) :: header  !             |header
       character (len=1) :: done
       integer :: eof !             |end of file
       logical :: i_exist               !             |check to determine if file exists
       integer :: iyr_start
       integer :: iyrc
       integer :: itot
       integer :: iyr
       integer :: iyr_co2
       real :: co2_end
              
      type co2
        integer :: iyr
        real :: co2
      end type co2
    
      type co2_annual
        !! Q870: DEFAULT-INITIALISED, and the default is the whole fix.
        !!
        !! When co2_yr.dat is absent -- the normal case; no SWATGenX deck ships one -- the
        !! reading branch below never executes and `yrs` was never assigned. Line ~79 then
        !! evaluates `co2_inc%co2_yr(co2_inc%yrs)%iyr`, indexing the UNALLOCATED co2_yr with
        !! whatever the stack held.
        !!
        !! The code already handles yrs == 0 correctly (it falls back to bsn_prm%co2); that
        !! branch was simply unreachable on a compiler that does not zero the component.
        !! gfortran's -finit-local-zero does NOT cover derived-type components: measured
        !! 2026-09-21 on box 2, SIGSEGV rc 139 at co2_read -> proc_bsn -> MAIN__, on a modified
        !! deck, an unmodified deck and a debug build alike. ifx evidently zeroes it, so the
        !! shipped engine takes the correct branch BY COMPILER BEHAVIOUR rather than by
        !! construction. This removes that dependence; it does not change any result where
        !! co2_yr.dat exists, which is the acceptance's byte-identity row.
        integer :: yrs = 0
        type (co2), dimension(:), allocatable :: co2_yr
      end type co2_annual
      type (co2_annual) :: co2_inc
      
      !! Q1281: co2_end is assigned only while the record index is inside co2_yr.dat but read once it runs past
      !! the end, which a file with gapped years reaches first. The OpenMP strip removed its implicit-SAVE '= 0.'.
      co2_end = 0.

      !! output annual CO2 
      call open_output_file(2222, "co2.out")
      write (2222,*) "         YR    CO2(ppm)"
         
       eof = 0
      
    !! read CO2 yearly values     
      inquire (file="co2_yr.dat", exist=i_exist)
        if (.not. i_exist .or. "co2_yr.dat" == " null") then
      else 
       do
          open (107,file="co2_yr.dat")
!! open output file
          read (107,*,iostat=eof) titldum
          if (eof < 0) exit
          read (107,*,iostat=eof) co2_inc%yrs
          if (eof < 0) exit
          read (107,*,iostat=eof) header
          if (eof < 0) exit
                
          allocate (co2_inc%co2_yr(co2_inc%yrs))
                   
          do itot = 1, co2_inc%yrs
            read (107,*,iostat=eof) co2_inc%co2_yr(itot)
            if (eof < 0) exit
          end do

          close (107)
          exit
        enddo
        endif
        
        allocate (co2y(time%nbyr), source = 0.)
        
        done = "n"
        !! no co2 file - default to basin parm
        if (co2_inc%yrs == 0) then
          do iyr = 1, time%nbyr
            co2y(iyr) = bsn_prm%co2
          end do
          done = "y"
        end if
        
        if (done == "n") then
          !! entire co2 period is before start of simulation
          if (time%yrc >= co2_inc%co2_yr(co2_inc%yrs)%iyr) then
            do iyr = 1, time%nbyr
              co2y(iyr) = co2_inc%co2_yr(co2_inc%yrs)%co2
            end do
            done = "y"
          end if
        end if
          
        if (done == "n") then
          iyr_start = time%yrc_start - co2_inc%co2_yr(1)%iyr + 1
          iyr_co2 = iyr_start
          do iyr = 1, time%nbyr
            !! co2 record hasn't started so set to basin parm co2
            if (iyr_co2 <= 0) then
              co2y(iyr) = co2_inc%co2_yr(1)%co2
            else
              if (iyr_co2 <= co2_inc%yrs) then
                co2y(iyr) = co2_inc%co2_yr(iyr_co2)%co2
                co2_end = co2y(iyr)
              else
                co2y(iyr) = co2_end
              end if
            end if
            iyr_co2 = iyr_co2 + 1
          end do
        end if
            
        !! write to co2.out
        iyrc = time%yrc_start
        do iyr = 1, time%nbyr
          write (2222,*) iyrc, co2y(iyr)
          iyrc = iyrc + 1
        end do
            
       return
      end subroutine co2_read
