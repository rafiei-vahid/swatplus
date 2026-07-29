       subroutine readcio_read 
    
       use input_file_module
       use output_path_module

       implicit none
           
       character (len=80) :: titldum
       character (len=15) :: name
       character (len=256) :: out_path_value
       character (len=512) :: line_buffer
       integer :: eof
       integer :: idx
       logical :: i_exist              !none       |check to determine if file exists
       integer :: i
       
       eof = 0
       
       !! read file.cio
       inquire (file="file.cio", exist=i_exist)
       if (i_exist ) then
         open (107,file="file.cio")
         read (107,*) titldum
      do i = 1, 31
         read (107,*,iostat=eof) name, in_sim  
         if (eof < 0) exit
         !! The basin line is read as TEXT, not straight into in_basin, so that a file.cio
         !! written before carbon.bsn existed still loads. in_basin gained a third field
         !! (carbon_bsn) upstream; a list-directed read of the whole derived type demands
         !! three tokens, so a two-token line silently consumes the FIRST TOKEN OF THE NEXT
         !! LINE. Every following filename then shifts by one and no database is read -- the
         !! symptom is an unrelated "array bound 0" crash much later (hru_db, urbdb, ru_elem).
         !! Every model written by SWAT+ Editor 3.0.8, which is what SWATGenX pins, has the
         !! two-token form. Parsing tokens keeps both vintages working and needs no editor
         !! upgrade; carbon_bsn keeps its default when the token is absent.
         read (107,'(a)',iostat=eof) line_buffer
         if (eof < 0) exit
         call cio_basin_line (line_buffer)
         read (107,*,iostat=eof) name, in_cli
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_con
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_cha
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_res
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_ru
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_hru
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_exco
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_rec
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_delr
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_aqu
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_herd
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_watrts
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_link
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_hyd
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_str
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_parmdb
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_ops
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_lum
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_chg
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_init
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_sol
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_cond
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_regs
         if (eof < 0) exit
!!!!!weather path code
         read (107,*,iostat=eof) name, in_path_pcp
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_path_tmp
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_path_slr
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_path_hmd
         if (eof < 0) exit
         read (107,*,iostat=eof) name, in_path_wnd
         if (eof < 0) exit
!!!!!weather path code
!!!!!output path code
         ! Read the whole line to handle spaces and avoid '/' termination in list-directed input
         read (107,'(A)',iostat=eof) line_buffer
         if (eof < 0) then
           out_path_value = ""
         else
           ! Parse the line: skip the label (first word), then get the rest
           line_buffer = adjustl(line_buffer)
           idx = index(line_buffer, ' ')
           if (idx > 0) then
             ! Found space after label, get the rest
             out_path_value = adjustl(line_buffer(idx+1:))
           else
             ! No value found
             out_path_value = ""
           end if
         end if
!!!!!output path code
         exit
      enddo
      endif

       close (107)
       
       !! Initialize output path (will use current dir if null/empty)
       call init_output_path(out_path_value)
            
       return

      contains

      subroutine cio_basin_line (line)
      !! Assign in_basin from the tokens actually present on the basin line of file.cio.
      !! Layout: <name> <codes_bas> <parms_bas> [carbon_bsn]
      !! A missing carbon_bsn leaves the component at its default ("carbon.bsn"), which is
      !! correct: carbon_bsn_read returns immediately unless cswat == 2, so a model with
      !! carbon off never touches the file whether or not file.cio names it.
      character(len=*), intent(in) :: line
      character(len=256) :: tok(8)
      integer :: ntok, pos, i0, n

      tok = ""
      ntok = 0
      pos = 1
      n = len_trim(line)
      do
        do while (pos <= n .and. line(pos:pos) == " ")
          pos = pos + 1
        end do
        if (pos > n .or. ntok >= size(tok)) exit
        i0 = pos
        do while (pos <= n .and. line(pos:pos) /= " ")
          pos = pos + 1
        end do
        ntok = ntok + 1
        tok(ntok) = line(i0:pos-1)
      end do

      !! tok(1) is the row label ("basin"); the filenames follow.
      if (ntok >= 2) in_basin%codes_bas  = trim(tok(2))
      if (ntok >= 3) in_basin%parms_bas  = trim(tok(3))
      if (ntok >= 4) in_basin%carbon_bsn = trim(tok(4))
      end subroutine cio_basin_line

      end subroutine readcio_read  
