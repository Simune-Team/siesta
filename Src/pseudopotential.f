      module pseudopotential

      use sys
      use precision
      
      private

      public :: pseudopotential_t, pseudo_read

      integer, parameter        :: nrmax = 1500


      type pseudopotential_t
        character(len=2)        :: name
        integer                 :: nr
        integer                 :: nrval
        real(dp)                :: zval
        logical                 :: relativistic
        character(len=10)       :: correlation
        character(len=2)        :: icorr
        character(len=3)        :: irel
        character(len=4)        :: nicore
        real(dp)                :: a
        real(dp)                :: b
        character(len=10)       :: method(6)
        character(len=70)       :: text
        integer                 :: npotu
        integer                 :: npotd
        real(dp), pointer       :: r(:)
        real(dp), pointer       :: chcore(:)
        real(dp), pointer       :: chval(:)
        real(dp), pointer       :: vdown(:,:)
        real(dp), pointer       :: vup(:,:)
        integer, pointer        :: ldown(:)
        integer, pointer        :: lup(:)
      end type pseudopotential_t

        CONTAINS

        subroutine pseudo_read(label,p)
        character(len=*), intent(in)   :: label
        type(pseudopotential_t)                    :: p

        character(len=40) fname
        logical found
        integer io_ps, i, j

        call vps_init(p)
        fname  = trim(label) // '.vps'
        inquire(file=fname, exist=found)
        if (.not.found) then
           write(6,'(/,2a,a20)') 'pseudo_read: WARNING: ',
     .          'Pseudopotential file not found: ', fname
           fname = trim(label) // '.psatom.data'
           write(6,'(2a)') 'pseudo_read: WARNING: Looking for ', fname
           inquire(file=fname, exist=found)
           if (.not.found) then
              write(6,'(/,2a,a20,/)') 'read_vps: ERROR: ',
     .             'Pseudopotential file not found: ', fname
              call die
           endif
        endif

        call io_assign(io_ps)
        open(io_ps,file=fname,form='unformatted',status='unknown')

        read(io_ps) p%name, p%icorr, p%irel, p%nicore,
     .       (p%method(i),i=1,6), p%text,
     .       p%npotd, p%npotu, p%nr, p%b, p%a, p%zval

        p%nrval = p%nr + 1
        allocate(p%r(1:p%nrval))
        read(io_ps) (p%r(j),j=2,p%nrval)
        p%r(1) = 0.d0

        if (p%npotd.gt.0) then
           allocate(p%vdown(1:p%npotd,1:p%nrval))
           allocate(p%ldown(1:p%npotd))
        endif
        do i=1,p%npotd
           read(io_ps) p%ldown(i), (p%vdown(i,j), j=2,p%nrval)
           p%vdown(i,1) = p%vdown(i,2)
        enddo

        if (p%npotu.gt.0) then
           allocate(p%vup(1:p%npotu,1:p%nrval))
           allocate(p%lup(1:p%npotu))
        endif
        do i=1,p%npotu
           read(io_ps) p%lup(i), (p%vup(i,j), j=2,p%nrval)
           p%vup(i,1) = p%vup(i,2)
        enddo

        allocate(p%chcore(1:p%nrval))
        allocate(p%chval(1:p%nrval))

        read(io_ps) (p%chcore(j),j=2,p%nrval)
        read(io_ps) (p%chval(j),j=2,p%nrval)
        r2=p%r(2)/(p%r(3)-p%r(2))
        p%chcore(1) = p%chcore(2) - r2*(p%chcore(3)-p%chcore(2))
        p%chval(1) = p%chval(2) - r2*(p%chval(3)-p%chval(2))

        call io_close(io_ps)
        end subroutine pseudo_read

        subroutine vps_init(p)
        type(pseudopotential_t)  :: p
        nullify(p%lup,p%ldown,p%r,p%chcore,p%chval,p%vdown,p%vup)
        end subroutine vps_init

        subroutine vps_print(p)
        type(pseudopotential_t)  :: p
        write(6,*) p%name
        write(6,*) p%nr
        write(6,*) p%lup
        write(6,*) p%ldown
        end subroutine vps_print

        end module pseudopotential



