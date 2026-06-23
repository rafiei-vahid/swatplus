      subroutine pfas_count_read

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Lightweight EARLY pass over pfas.dat: count the number of simulated PFAS
!!    compounds (npfas) and, if > 0, ACTIVATE the SWAT+ constituent-transport
!!    infrastructure by bumping cs_db%num_tot.  That makes hyd_read_connect
!!    allocate the per-object obcs constituent hydrographs (and set obcs_alloc),
!!    which the in-stream PFAS routing rides -- so PFAS can route on ANY model
!!    without requiring a pesticide / constituents.cs.  The per-type counts
!!    (num_pests / num_paths / num_metals / num_salts / num_cs) stay 0, so every
!!    constituent inner loop is a 0-trip no-op; only the %pfas hydrograph slot
!!    (added by pfas_cha_read) actually carries mass.
!!
!!    cs_db%num_tot is used ONLY as a `> 0` gate everywhere (verified), never as
!!    a loop bound or array size, so bumping it is safe.
!!
!!    Ordering: MUST run AFTER constit_db_read (which sets the base num_tot) and
!!    BEFORE hyd_connect/hyd_read_connect (which allocates obcs).  The full PFAS
!!    database + per-HRU soil pools are read later by pfas_read (which recounts
!!    npfas to the same value -- idempotent).  Missing pfas.dat -> clean no-op.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      use pfas_module, only : npfas
      use constituent_mass_module, only : cs_db

      implicit none

      character (len=80) :: titldum = ""   !! title line
      character (len=80) :: header = ""    !! header line
      integer :: eof = 0                   !! end-of-file flag
      integer :: imax = 0                  !! compound-record count
      integer :: id = 0                    !! PFAS id read from a record
      logical :: i_exist = .false.         !! file-existence flag

      npfas = 0
      inquire (file="pfas.dat", exist=i_exist)
      if (.not. i_exist) return

      open (170, file="pfas.dat")
      read (170,*,iostat=eof) titldum
      if (eof < 0) then
        close (170)
        return
      end if
      read (170,*,iostat=eof) header
      if (eof < 0) then
        close (170)
        return
      end if

      !! count compound records (stop at id<=0 or EOF) -- mirrors pfas_read pass 1
      do
        read (170,*,iostat=eof) id
        if (eof < 0) exit
        if (id <= 0) exit
        imax = imax + 1
      end do
      close (170)

      npfas = imax
      if (npfas > 0) cs_db%num_tot = cs_db%num_tot + npfas

      return
      end subroutine pfas_count_read
