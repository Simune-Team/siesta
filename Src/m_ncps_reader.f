!> @brief Consolidates the reading of all allowable types of ps files
!! (vps, psf, psml)
!> @author Alberto Garcia

      module m_ncps_reader

      use m_ncps_froyen_ps_t,    only: pseudopotential_t => froyen_ps_t

      integer, parameter, private :: dp = selected_real_kind(14,100)

      public :: pseudo_read

      CONTAINS

        subroutine pseudo_read(label,p,new_grid,a,b,rmax)
        use m_ncps_froyen_reader,  only: pseudo_read_formatted
        use m_ncps_froyen_reader,  only: pseudo_read_unformatted
        use m_ncps_froyen_reader,  only: pseudo_reparametrize

        character(len=*), intent(in)   :: label
        type(pseudopotential_t)        :: p
        logical, intent(in), optional  :: new_grid
        real(dp), intent(in), optional :: a
        real(dp), intent(in), optional :: b
        real(dp), intent(in), optional :: rmax

!       PS information can be in a .vps file (unformatted)
!       or in a .psf file (formatted)
!       or in a .psml file 

        character(len=40) fname
        logical found, reparametrize

        reparametrize = .false.
        if (present(new_grid)) then
           reparametrize = new_grid
        endif
        if (reparametrize) then
           if (.not. present(a)) call die("New a not present")
           if (.not. present(b)) call die("New b not present")
        endif

        fname  = trim(label) // '.vps'
        inquire(file=fname, exist=found)
        if (found) then
           call pseudo_read_unformatted(fname,p)
           if (reparametrize) then
              call pseudo_reparametrize(p,a,b,label,rmax)
           endif
        else
           fname = trim(label) // '.psf'
           inquire(file=fname, exist=found)
           if (found) then
              call pseudo_read_formatted(fname,p)
              if (reparametrize) then
                 call pseudo_reparametrize(p,a,b,label,rmax)
              endif
           else
              fname = trim(label) // '.psml'
              inquire(file=fname, exist=found)
              if (found) then
                 call pseudo_read_psml(fname,p,reparametrize,a,b,rmax)
              else
                 write(6,'(/,2a,a20,/)') 'read_pseudo: ERROR: ',
     .                'Pseudopotential file not found: ', fname
                 call die
              endif
           endif
        endif
        call pseudo_dump(trim(label) // ".psdump",p)
        end subroutine pseudo_read
!
        subroutine pseudo_read_psml(fname,p,reparametrize,a,b,rmax)

        use m_psml, only: ps_t, ps_destroy, psml_reader
        use m_ncps_translators, only: ncps_xml2froyen_new

        character(len=*), intent(in)              :: fname
        type(pseudopotential_t), intent(out)      :: p
        logical, intent(in), optional  :: reparametrize
        real(dp), intent(in), optional :: a
        real(dp), intent(in), optional :: b
        real(dp), intent(in), optional :: rmax

        ! Use the target attribute as per the standard
        ! warning about dangling association...
        type(ps_t), target   :: ps

        call psml_reader(fname,ps)
        call ncps_xml2froyen_new(ps,p,reparametrize,a,b,rmax)
        call ps_destroy(ps)

        end subroutine pseudo_read_psml
!----
        subroutine pseudo_dump(fname,p)
!
!       Column-oriented output
!
        character(len=*), intent(in) :: fname
        type(pseudopotential_t), intent(in)     :: p

        integer io_ps, i, j

        call io_assign(io_ps)
        open(io_ps,file=fname,form='formatted',status='unknown',
     $       action="write",position="rewind")
        write(6,'(3a)') 'Dumping pseudopotential information ',
     $       'in formatted form in ', trim(fname)

 9040    format(i4,7es20.9)
         do j = 1, p%nrval
            write(io_ps,9040) j, p%r(j), (p%vdown(i,j),i=1,p%npotd),
     $                        p%chval(j), p%chcore(j)
         enddo
         call io_close(io_ps)
         end subroutine pseudo_dump

      end module m_ncps_reader
