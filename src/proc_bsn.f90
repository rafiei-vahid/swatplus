      subroutine proc_bsn

      use time_module
      use output_path_module
      
      implicit none
      
      external :: basin_print_codes_read, basin_prm_default, basin_read_cc, basin_read_objs, &
                  basin_read_prm, carbon_bsn_read, carbon_layers_read, &
                  co2_read, readcio_read, time_read
      
      call readcio_read

!!!  open file to print all output files that are written
      call open_output_file(9000, "files_out.out")
      write (9000,*) "files_out.out - OUTPUT FILES WRITTEN"      
!!!  open diagnostics.out file to print problems with various files
     call open_output_file(9001, "diagnostics.out", 8000)  !!ext to 8000 recl per email 8/2/21 - Kai-Uwe
     write (9001,*) "DIAGNOSTICS.OUT FILE" 
!!!  open drainage areas output file
     !! Q1049 fork fix J4 (Lane P ruling J4-a): hyd_connect writes one list-directed record per object, and an object with
     !! ~1,400+ inflows needs ~110-130 KB. At RECL 80000 gfortran STOPS ("End of record") and ifx wraps silently (37 of 223
     !! served models). 2**24 is > 100x the longest measured record and keeps a bounded record buffer (a larger RECL, or no
     !! RECL, which makes ifx wrap list-directed output at its short default, are both avoided on purpose).
     call open_output_file(9004, "area_calc.out", 16777216)
                
      call basin_read_cc
      call basin_read_objs
      call time_read
      
      !if (time%step > 0) then
        time%dtm = 1440. / time%step
      !end if
             
      call basin_read_prm
      call basin_prm_default
      call basin_print_codes_read
      call co2_read
      call carbon_bsn_read
      call carbon_layers_read

      return
      
      end subroutine proc_bsn
