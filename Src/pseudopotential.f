      module pseudopotential

      use sys
      use precision
      use ionew
      
      private

      public :: pseudopotential_t, pseudo_read, pseudo_header_print

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

!       PS information can be in a .vps file (unformatted)
!       or in a .psf file (formatted)

        character(len=40) fname
        logical found

        call vps_init(p)

        fname  = trim(label) // '.vps'
        inquire(file=fname, exist=found)
        if (found) then
           call pseudo_read_unformatted(fname,p)
        else
           fname = trim(label) // '.psf'
           inquire(file=fname, exist=found)
           if (found) then
              call pseudo_read_formatted(fname,p)
           else
              write(6,'(/,2a,a20,/)') 'read_pseudo: ERROR: ',
     .             'Pseudopotential file not found: ', fname
              call die
           endif
        endif
        end subroutine pseudo_read
!
        subroutine pseudo_read_unformatted(fname,p)
        character(len=*), intent(in) :: fname
        type(pseudopotential_t)                    :: p

        integer io_ps, i, j

        call io_assign(io_ps)
        open(io_ps,file=fname,form='unformatted',status='unknown')
        write(6,'(3a)') 'Reading pseudopotential information ',
     $       'in unformatted form from ', trim(fname)

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
        end subroutine pseudo_read_unformatted
!----
        subroutine pseudo_read_formatted(fname,p)
        character(len=*), intent(in) :: fname
        type(pseudopotential_t)                    :: p

        integer io_ps, i, j
        character(len=70) dummy

        call io_assign(io_ps)
        open(io_ps,file=fname,form='formatted',status='unknown')
        write(6,'(3a)') 'Reading pseudopotential information ',
     $       'in formatted form from ', trim(fname)

 8000   format(1x,i2)
 8005   format(1x,a2,1x,a2,1x,a3,1x,a4)
 8010   format(1x,6a10,/,1x,a70)
 8015   format(1x,2i3,i5,3g20.12)
 8030   format(4(g20.12))
 8040   format(1x,a)

        read(io_ps,8005) p%name, p%icorr, p%irel, p%nicore
        read(io_ps,8010) (p%method(i),i=1,6), p%text
        read(io_ps,8015) p%npotd, p%npotu, p%nr, p%b, p%a, p%zval

        p%nrval = p%nr + 1
        allocate(p%r(1:p%nrval))
        read(io_ps,8040) dummy
        read(io_ps,8030) (p%r(j),j=2,p%nrval)
        p%r(1) = 0.d0

        if (p%npotd.gt.0) then
           allocate(p%vdown(1:p%npotd,1:p%nrval))
           allocate(p%ldown(1:p%npotd))
        endif
        do i=1,p%npotd
           read(io_ps,8040) dummy 
           read(io_ps,8000) p%ldown(i)
           read(io_ps,8030) (p%vdown(i,j), j=2,p%nrval)
           p%vdown(i,1) = p%vdown(i,2)
        enddo

        if (p%npotu.gt.0) then
           allocate(p%vup(1:p%npotu,1:p%nrval))
           allocate(p%lup(1:p%npotu))
        endif
        do i=1,p%npotu
           read(io_ps,8040) dummy 
           read(io_ps,8000) p%lup(i)
           read(io_ps,8030) (p%vup(i,j), j=2,p%nrval)
           p%vup(i,1) = p%vup(i,2)
        enddo

        allocate(p%chcore(1:p%nrval))
        allocate(p%chval(1:p%nrval))

        read(io_ps,8040) dummy
        read(io_ps,8030) (p%chcore(j),j=2,p%nrval)
        read(io_ps,8040) dummy
        read(io_ps,8030) (p%chval(j),j=2,p%nrval)
        r2=p%r(2)/(p%r(3)-p%r(2))
        p%chcore(1) = p%chcore(2) - r2*(p%chcore(3)-p%chcore(2))
        p%chval(1) = p%chval(2) - r2*(p%chval(3)-p%chval(2))

        call io_close(io_ps)
        end subroutine pseudo_read_formatted
!------

        subroutine vps_init(p)
        type(pseudopotential_t)  :: p
        nullify(p%lup,p%ldown,p%r,p%chcore,p%chval,p%vdown,p%vup)
        end subroutine vps_init

!-------
        subroutine pseudo_header_print(lun,p)
        integer, intent(in) :: lun
        type(pseudopotential_t)  :: p

 8005   format(1x,a2,1x,a2,1x,a3,1x,a4)
 8010   format(1x,6a10,/,1x,a70)
 8015   format(1x,2i3,i5,3g20.12)
        
        write(lun,'(a)') '<pseudopotential_header>'
        write(lun,8005) p%name, p%icorr, p%irel, p%nicore
        write(lun,8010) (p%method(i),i=1,6), p%text
        write(lun,'(a)') '</pseudopotential_header>'

        end subroutine pseudo_header_print
!--------
        end module pseudopotential



