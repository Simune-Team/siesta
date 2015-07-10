!> @brief Consolidates the reading of all allowable types of ps files
!! (vps, psf, psml)
!> @author Alberto Garcia

      module m_ncps_reader

      use m_ncps_froyen_ps_t,    only: pseudopotential_t => froyen_ps_t

      integer, parameter, private :: dp = selected_real_kind(14,100)

      public :: pseudo_read, pseudo_read_from_file

      CONTAINS

        subroutine pseudo_read(label,p,
     $                         psml_handle,has_psml_ps,
     $                         new_grid,a,b,rmax,directory)

        use m_ncps_froyen_reader,  only: pseudo_read_formatted
        use m_ncps_froyen_reader,  only: pseudo_read_unformatted
        use m_ncps_froyen_reader,  only: pseudo_reparametrize
        use m_psml,                only: psml_t => ps_t

        character(len=*), intent(in)   :: label
        type(pseudopotential_t)        :: p
        type(psml_t), intent(inout), target :: psml_handle
        logical, intent(out)           :: has_psml_ps
        logical, intent(in), optional  :: new_grid
        real(dp), intent(in), optional :: a
        real(dp), intent(in), optional :: b
        real(dp), intent(in), optional :: rmax
        character(len=*), intent(in), optional   :: directory

!       PS information can be in a .vps file (unformatted)
!       or in a .psf file (formatted)
!       or in a .psml file 

        character(len=200) fname, prefix
        logical found, reparametrize

        has_psml_ps = .false.

        reparametrize = .false.
        if (present(new_grid)) then
           reparametrize = new_grid
        endif
        if (reparametrize) then
           if (.not. present(a)) call die("New a not present")
           if (.not. present(b)) call die("New b not present")
        endif

        prefix = ""
        if (present(directory)) then
           prefix = trim(directory) // "/"
        endif

        fname  = trim(prefix) // trim(label) // '.vps'
        inquire(file=fname, exist=found)
        if (found) then
           write(6,'(/,a,a,/)') 
     .          'Reading pseudopotential from: ', trim(fname)
           call pseudo_read_unformatted(fname,p)
           if (reparametrize) then
              call pseudo_reparametrize(p,a,b,label,rmax)
           endif
        else
           fname = trim(prefix) // trim(label) // '.psf'
           inquire(file=fname, exist=found)
           if (found) then
              write(6,'(/,a,a,/)') 
     .                'Reading pseudopotential from: ', trim(fname)
              call pseudo_read_formatted(fname,p)
              if (reparametrize) then
                 call pseudo_reparametrize(p,a,b,label,rmax)
              endif
           else
              fname = trim(prefix) // trim(label) // '.psml'
              inquire(file=fname, exist=found)
              if (found) then
                 write(6,'(/,a,a,/)') 
     .                'Reading pseudopotential from: ', trim(fname)
                 call pseudo_read_psml(fname,p,psml_handle,
     $                                 reparametrize,a,b,rmax)
                 has_psml_ps = .true.
              else
                 write(6,'(2a,a)') 'pseudo_read: ERROR: ',
     .                'Pseudopotential file not found: ',
     $                trim(prefix) // trim(label) // '.{psf,vps,psml}'

                 call die("")
              endif
           endif
        endif
        ! Dump locally
        call pseudo_dump(trim(label) // ".psdump",p)
        end subroutine pseudo_read

        subroutine pseudo_read_from_file(filename,p,
     $                                   new_grid,a,b,rmax)

        use m_ncps_froyen_reader,  only: pseudo_read_formatted
        use m_ncps_froyen_reader,  only: pseudo_read_unformatted
        use m_ncps_froyen_reader,  only: pseudo_reparametrize

        character(len=*), intent(in)   :: filename
        type(pseudopotential_t)        :: p

        logical, intent(in), optional  :: new_grid
        real(dp), intent(in), optional :: a
        real(dp), intent(in), optional :: b
        real(dp), intent(in), optional :: rmax

        character(len=30)   :: label, ext
        integer :: status

        logical reparametrize

        reparametrize = .false.
        if (present(new_grid)) then
           reparametrize = new_grid
        endif
        if (reparametrize) then
           if (.not. present(a)) call die("New a not present")
           if (.not. present(b)) call die("New b not present")
        endif

        call get_label_ext(filename,label,ext,status)
        if (status /= 0) call die("Cannot get label and extension")
        if (trim(ext) == ".vps") then
           call pseudo_read_unformatted(filename,p)
           if (reparametrize) then
              call pseudo_reparametrize(p,a,b,label,rmax)
           endif
        else if (trim(ext) == ".psf") then
           call pseudo_read_formatted(filename,p)
           if (reparametrize) then
              call pseudo_reparametrize(p,a,b,label,rmax)
           endif
        else if (trim(ext) == ".psml") then
           call pseudo_read_psml(filename,p,
     $          reparametrize=reparametrize,a=a,b=b,rmax=rmax)
        else
           write(6,'(2a,a)') 'pseudo_read_from_file: ERROR: ',
     .                'Extension not supported: ', trim(ext)
           call die("")
        endif
        ! Dump locally
        call pseudo_dump(trim(label) // ".psdump",p)
        end subroutine pseudo_read_from_file
!
        subroutine pseudo_read_psml(fname,p,
     $                              psml_handle,
     $                              reparametrize,a,b,rmax)

        use m_psml, only: ps_t, ps_destroy, psml_reader
        use m_ncps_translators, only: ncps_xml2froyen_new

        character(len=*), intent(in)              :: fname
        type(pseudopotential_t), intent(out)      :: p
        type(ps_t), intent(inout), optional, target :: psml_handle
        logical, intent(in), optional  :: reparametrize
        real(dp), intent(in), optional :: a
        real(dp), intent(in), optional :: b
        real(dp), intent(in), optional :: rmax

        ! Use the target attribute as per the standard
        ! warning about dangling association...
        type(ps_t), target   :: ps

        if (present(psml_handle)) then
           ! We pass the actual handle to the caller
           call psml_reader(fname,psml_handle)
           call ncps_xml2froyen_new(psml_handle,p,
     $                              reparametrize,a,b,rmax)
        else
           ! We just convert to Froyen form and destroy ps
           call psml_reader(fname,ps)
           call ncps_xml2froyen_new(ps,p,reparametrize,a,b,rmax)
           call ps_destroy(ps)
        endif

        end subroutine pseudo_read_psml
!----
        subroutine pseudo_dump(fname,p)
!
!       Column-oriented output
!
        character(len=*), intent(in) :: fname
        type(pseudopotential_t), intent(in)     :: p

        integer io_ps, i, j

        call get_free_lun(io_ps)
        open(io_ps,file=fname,form='formatted',status='unknown',
     $       action="write",position="rewind")
        write(6,'(3a)') 'Dumping pseudopotential information ',
     $       'in formatted form in ', trim(fname)

 9040    format(i4,7es20.9)
         do j = 1, p%nrval
            write(io_ps,9040) j, p%r(j), (p%vdown(i,j),i=1,p%npotd),
     $                        p%chval(j), p%chcore(j)
         enddo
         close(io_ps)
         end subroutine pseudo_dump

      subroutine get_free_lun(lun)
      integer, intent(out) :: lun

      interface
         subroutine die(str)
         character(len=*), intent(in), optional :: str
         end subroutine die
      end interface

      logical :: used
      integer :: iostat

      do lun= 10,90
         inquire(unit=lun, opened=used, iostat=iostat)
         if (iostat .ne. 0) used = .true.
         if (.not. used) return ! normal return with 'lun' value
      enddo
      call die("No luns available")

      end subroutine get_free_lun

      subroutine get_label_ext(str,label,ext,stat)
      character(len=*), intent(in)   :: str
      character(len=*), intent(out)  :: label
      character(len=*), intent(out)  :: ext
      integer, intent(out)           :: stat

      integer n, i, lo, hi, bar, dot

      n = len_trim(str)
      stat = -1
      dot = -1
      bar = 0
      do i = n, 1, -1
!     print *, "i, c:", i, "|",str(i:i),"|"
         if ( (str(i:i) == ".") .and. (dot == -1) ) then
            dot = i
!     print *, "dot set to: ", dot
         endif
         if ( (str(i:i) == "/") .and. (bar == 0) ) then
            bar = i
!     print *, "bar set to: ", bar
         endif
      enddo

      if ( (dot > 1) .and. (dot>bar)) then
         stat = 0
         label=str(bar+1:dot-1)
         ext=str(dot:n)
      endif

      end subroutine get_label_ext

      end module m_ncps_reader
