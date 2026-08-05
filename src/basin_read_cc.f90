       subroutine basin_read_cc
      
       use input_file_module
       use basin_module
      
       implicit none
       
       character (len=80) :: titldum  !             |title of file
       character (len=80) :: header  !             |header
       integer :: eof !             |end of file
       logical :: i_exist               !             |check to determine if file exists
       
       eof = 0
      
       !! read basin
       inquire (file=in_basin%codes_bas, exist=i_exist)
       if (i_exist .or. in_basin%codes_bas /= "null") then      
       do 
         open (107,file=in_basin%codes_bas)
         read (107,*,iostat=eof) titldum
         if (eof < 0) exit
         read (107,*,iostat=eof) header
         if (eof < 0) exit
         read (107,*,iostat=eof) bsn_cc
         if (eof < 0) exit
         exit
       enddo
       endif
       
       if (bsn_cc%pet == 3) then 
        open (140,file = 'pet.cli')
       do
        read (140,*,iostat=eof) titldum
        if (eof < 0) exit
        read (140,*,iostat=eof) header
        if (eof < 0) exit
        read (140,*,iostat = eof) titldum
        exit
       end do
       end if

       close(107)

       !! cfac override, RELOCATED FROM ero_cfactor.f90:51 (2026-08-03).
       !! It was executed there on every HRU on every day, from inside the HRU-parallel
       !! region, so every worker thread wrote this one address concurrently -- ThreadSanitizer
       !! reported it 14 times with both accesses at that line. Every thread wrote the SAME
       !! constant, so the race could not alter results; moving it here removes the race with
       !! no change of value.
       !! Behaviour is unchanged: bsn_cc%cfac is referenced NOWHERE else in the engine except
       !! the `if (bsn_cc%cfac == 0)` test two lines below the old assignment, which the
       !! assignment made unreachable. The configured value was already dead.
       bsn_cc%cfac = 1

       return
      end subroutine basin_read_cc
