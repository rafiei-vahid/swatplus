      subroutine path_hru_aqu_read

      use constituent_mass_module
      use input_file_module
      use maximum_data_module

      implicit none

      character (len=80) :: titldum = ""
      character (len=500) :: header = ""
      character (len=500) :: line = ""
      character (len=16) :: name = ""
      integer :: ipathi = 0
      integer :: eof = 0
      integer :: imax = 0
      logical :: i_exist              !none       |check to determine if file exists
      logical :: rows = .false.       !none       |.true. if each set holds one 'pathogen soil plant' row per pathogen
      logical :: at_end = .false.     !none       |.true. once a read has met the end of the file
      logical :: found = .false.      !none       |.true. if a set was read
      integer :: ip = 0
      real :: w = 0.
      real :: b = 0.
      real, dimension(:), allocatable :: soi
      real, dimension(:), allocatable :: plt

      eof = 0

      !read all export coefficient data
      inquire (file=in_init%path_soil, exist=i_exist)
      if (i_exist .or. in_init%path_soil /= "null") then
        do
          open (107,file=in_init%path_soil)
          read (107,*,iostat=eof) titldum
          if (eof < 0) exit
          allocate (soi(cs_db%num_paths), source = 0.)
          allocate (plt(cs_db%num_paths), source = 0.)

          !! the header: some files repeat it before every set, the SWAT+ Editor writes it once; a record that repeats it
          !! before a set is skipped.
          !! the layout, from the first set: the SWAT+ Editor writes one 'pathogen soil plant' row per pathogen. The other
          !! layout is one list-directed read of 'name soil(1:n) plant(1:n)', which may span records. The sets are read as
          !! rows when one of the first set's row-shaped records names a simulated pathogen, unless the record after them
          !! starts with a number (it then continues a list-directed read).
          at_end = .false.
          call next_record
          header = line
          rows = .false.
          if (.not. at_end) then
            call next_record
            if (.not. at_end) then
              do
                read (107,'(a)',iostat=eof) line
                if (eof /= 0) exit
                ip = path_row (line)
                if (ip < 0) exit
                if (ip > 0) rows = .true.
              end do
              if (eof == 0) then
                if (lead_number (line)) rows = .false.
              end if
            end if
          end if

          !! count the sets with the same reads that store them
          call skip_title_header
          imax = 0
          do
            call path_soil_set (found)
            if (.not. found) exit
            imax = imax + 1
          end do

          db_mx%path_ini = imax

          allocate (path_soil_ini(imax))
          allocate (cs_path_solsor(cs_db%num_paths))

          call skip_title_header
          do ipathi = 1, imax
            call path_soil_set (found)
            path_soil_ini(ipathi)%name = name
            allocate (path_soil_ini(ipathi)%soil(cs_db%num_paths), source = soi)
            allocate (path_soil_ini(ipathi)%plt(cs_db%num_paths), source = plt)
          end do
          close (107)
          exit
        end do
      end if

      return

      contains

      !! the next record that is not blank, into line; at_end is .true. if the file ends first
      subroutine next_record
      do
        read (107,'(a)',iostat=eof) line
        if (eof /= 0) then
          at_end = .true.
          line = ""
          return
        end if
        if (len_trim (line) > 0) return
      end do
      end subroutine next_record

      !! back to the first set: past the title and the header record
      subroutine skip_title_header
      rewind (107)
      at_end = .false.
      call next_record
      if (.not. at_end) call next_record
      end subroutine skip_title_header

      !! one set: its name record (after a repeated header, if any), then its values into soi and plt; got is .false. when
      !! the file ends before them
      subroutine path_soil_set (got)
      logical, intent (out) :: got
      integer :: ip
      integer :: ios
      got = .false.
      soi = 0.
      plt = 0.
      if (at_end) return
      call next_record
      if (at_end) return
      if (squeezed (line) == squeezed (header)) then
        call next_record
        if (at_end) return
      end if
      name = ""
      read (line,*,iostat=ios) name
      if (rows) then
        !! rows up to the next record that is not one (a header, the next set's name) or the end of the file;
        !! a row naming a pathogen that is not simulated is skipped
        do
          read (107,'(a)',iostat=eof) line
          if (eof /= 0) then
            at_end = .true.
            exit
          end if
          ip = path_row (line)
          if (ip < 0) then
            backspace (107)
            exit
          end if
          if (ip > 0) then
            soi(ip) = w
            plt(ip) = b
          end if
        end do
      else
        read (107,*,iostat=eof) titldum, soi, plt
        if (eof < 0) then
          at_end = .true.
          return
        end if
      end if
      got = .true.
      end subroutine path_soil_set

      !! str without its blanks, so a repeated header matches whatever its spacing
      function squeezed (str) result (sq)
      character (len=*), intent (in) :: str
      character (len=len(str)) :: sq
      integer :: i
      integer :: k
      sq = ""
      k = 0
      do i = 1, len_trim (str)
        if (str(i:i) /= " " .and. str(i:i) /= char(9)) then
          k = k + 1
          sq(k:k) = str(i:i)
        end if
      end do
      end function squeezed

      !! a 'pathogen soil plant' row: the simulated pathogen's number, with its values in w and b; 0 for a pathogen
      !! that is not simulated; -1 if str does not hold exactly three items, the last two numbers
      integer function path_row (str)
      character (len=*), intent (in) :: str
      character (len=16) :: pname
      character (len=16) :: extra
      integer :: ios
      integer :: i
      path_row = -1
      read (str,*,iostat=ios) pname, w, b
      if (ios /= 0) return
      read (str,*,iostat=ios) pname, w, b, extra
      if (ios == 0) return
      read (str,*,iostat=ios) pname, w, b
      path_row = 0
      do i = 1, cs_db%num_paths
        if (cs_db%paths(i) == pname) path_row = i
      end do
      end function path_row

      !! .true. if the first item of str is a number
      logical function lead_number (str)
      character (len=*), intent (in) :: str
      real :: x
      integer :: ios
      read (str,*,iostat=ios) x
      lead_number = ios == 0
      end function lead_number

      end subroutine path_hru_aqu_read
