      module m_ncps_froyen_reader

      use m_ncps_froyen_ps_t,    only: pseudopotential_t => froyen_ps_t

      public :: pseudo_read_unformatted
      public :: pseudo_read_formatted
      public :: pseudo_reparametrize

      integer, parameter, private :: dp = selected_real_kind(10,100)

      CONTAINS

        subroutine pseudo_read_unformatted(fname,p)
        character(len=*), intent(in) :: fname
        type(pseudopotential_t)                    :: p

        integer io_ps, i, j
        real(dp) :: r2

        call get_free_lun(io_ps)
        open(io_ps,file=fname,form='unformatted',status='unknown')
        write(6,'(3a)') 'Reading pseudopotential information ',
     $       'in unformatted form from ', trim(fname)

        read(io_ps) p%name, p%icorr, p%irel, p%nicore,
     .       (p%method(i),i=1,6), p%text,
     .       p%npotd, p%npotu, p%nr, p%b, p%a, p%zval

!
!       Old style vps files should have the right info in text.
!
        call read_ps_conf(p%irel,p%npotd-1,p%text,p%gen_zval)

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

        close(io_ps)
        end subroutine pseudo_read_unformatted
!----
        subroutine pseudo_read_formatted(fname,p)
        character(len=*), intent(in) :: fname
        type(pseudopotential_t)                    :: p

        integer io_ps, i, j, ios
        character(len=70) dummy
        real(dp) :: r2, gen_zval_inline

        call get_free_lun(io_ps)
        open(io_ps,file=fname,form='formatted',status='unknown')
        write(6,'(3a)') 'Reading pseudopotential information ',
     $       'in formatted form from ', trim(fname)

 8000   format(1x,i2)
 8005   format(1x,a2,1x,a2,1x,a3,1x,a4)
 8010   format(1x,6a10,/,1x,a70)
 8015   format(1x,2i3,i5,4g20.12)
 8030   format(4(g20.12))
 8040   format(1x,a)

        read(io_ps,8005) p%name, p%icorr, p%irel, p%nicore
        read(io_ps,8010) (p%method(i),i=1,6), p%text
        read(io_ps,8015,iostat=ios)
     $       p%npotd, p%npotu, p%nr, p%b, p%a, p%zval,
     $       gen_zval_inline
        if (ios < 0) gen_zval_inline = 0.0_dp
        call read_ps_conf(p%irel,p%npotd-1,p%text,p%gen_zval)
!
!       (Some .psf files contain an extra field corresponding
!       to the ps valence charge at generation time. If that
!       field is not present, the information has to be decoded
!       from the "text" variable.
!
!       "Zero" pseudos have gen_zval = 0, so they need a special case.

        if (p%gen_zval == 0.0_dp) then
           if (gen_zval_inline == 0.0_dp) then
              if (p%method(1) /= "ZEROPSEUDO")
     $             call die("Cannot get gen_zval")
           else
              p%gen_zval = gen_zval_inline
           endif
        endif

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

        close(io_ps)
        end subroutine pseudo_read_formatted
!------

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

      subroutine read_ps_conf(irel,lmax,text,chgvps)
!
!     Attempt to decode the valence configuration used for
!     the generation of the pseudopotential
!     (At least, the valence charge)

      character(len=3), intent(in)  :: irel
      integer, intent(in)           :: lmax
      character(len=70), intent(in) :: text
      real(dp), intent(out)         :: chgvps

      integer  :: l, itext
      real(dp) :: ztot, zup, zdown, rc_read
      character(len=2) :: orb

      chgvps=0.0_dp

            if(irel.eq.'isp') then
               write(6,'(/,2a)')
     .          'Pseudopotential generated from an ',
     .          'atomic spin-polarized calculation'

               write(6,'(/,a)') 'Valence configuration '//
     .                 'for pseudopotential generation:'

               do l=0,min(lmax,3)
                  itext=l*17
                  read(text(itext+1:),err=5000,fmt=8080)
     $                 orb, zdown, zup, rc_read
 8080             format(a2,f4.2,1x,f4.2,1x,f4.2)
                  chgvps = chgvps + zdown + zup
                  write(6,8085) orb, zdown, zup, rc_read
 8085             format(a2,'(',f4.2,',',f4.2,') rc: ',f4.2)
               enddo

            else
               if(irel.eq.'rel') then
                  write(6,'(/,2a)')
     .          'Pseudopotential generated from a ',
     .                 'relativistic atomic calculation'
                  write(6,'(2a)')
     .          'There are spin-orbit pseudopotentials ',
     .                 'available'
                  write(6,'(2a)')
     .          'Spin-orbit interaction is not included in ',
     .                 'this calculation'
               endif

               write(6,'(/,a)') 'Valence configuration '//
     .                 'for pseudopotential generation:'

               do l=0,min(lmax,3)
                  itext=l*17
                  read(text(itext+1:),err=5000,fmt=8090)
     $                 orb, ztot, rc_read
 8090             format(a2,f5.2,4x,f5.2)
                  chgvps = chgvps + ztot
                  write(6,8095) orb, ztot, rc_read
 8095             format(a2,'(',f5.2,') rc: ',f4.2)
               enddo

           endif
           return

 5000    continue       ! Error return: set chgvps to zero

         end subroutine read_ps_conf

         subroutine pseudo_reparametrize(p,a,b,label,new_rmax)
         use flib_spline
!
!        Interpolate values into new grid, given by a and b
!
!        Typical new values:  a = 5x10-4, b=10

         type(pseudopotential_t)          :: p
         real(dp), intent(in)             :: a, b
         character(len=*), intent(in)     :: label
         real(dp), intent(in), optional   :: new_rmax

         real(dp)  :: rmax, rpb, ea, ea2, rr
         integer   :: ir, new_nrval, i, j
         real(dp), dimension(:), pointer   :: func, tmp, new_r
         real(dp), dimension(:,:), pointer :: tmp2

         real(dp), dimension(:), allocatable :: y2 

         if (present(new_rmax)) then
            rmax = new_rmax
            if (rmax < 1.0_dp) rmax = p%r(p%nrval)
         else
            rmax = p%r(p%nrval)
         endif
         print *, "Reparametrization. rmax: ", rmax
         rpb=b
         ea=exp(a)
         ea2=1.0d0
         ir = 0
         do 
            rr = b*(ea2-1.0d0)
            if (rr > rmax) exit
            ir = ir + 1
            rpb=rpb*ea
            ea2=ea2*ea
         enddo
         new_nrval = ir
         allocate(new_r(new_nrval))   ! Will go in derived type
         print *, "Reparametrization. New nrval: ", new_nrval

         rpb=b
         ea=exp(a)
         ea2=1.0d0
         do ir = 1, new_nrval
            new_r(ir) = b*(ea2-1.0d0)
            rpb=rpb*ea
            ea2=ea2*ea
         enddo
         
        allocate(y2(1:p%nrval))
!-----------------------------------------------------------------------
!       Basic idiom to reparametrize
!       Use natural spline (zero second derivative)
!
        func => p%chcore
        call generate_spline(p%r,func,p%nrval,y2,0.0_dp,0.0_dp)
        allocate(tmp(new_nrval))
        do j = 1, new_nrval
           call evaluate_spline(p%r,func,y2,p%nrval,new_r(j),tmp(j))
        enddo
        nullify(func)
        deallocate(p%chcore)    ! Old data
        p%chcore => tmp         ! Point to new memory area
        nullify(tmp)            ! To re-use tmp
!--------------------------------------------------------------------
        func => p%chval
        call generate_spline(p%r,func,p%nrval,y2,0.0_dp,0.0_dp)
        allocate(tmp(new_nrval))
        do j = 1, new_nrval
           call evaluate_spline(p%r,func,y2,p%nrval,new_r(j),tmp(j))
        enddo
        nullify(func)
        deallocate(p%chval)    ! Old data
        p%chval => tmp         ! Point to new memory area
        nullify(tmp)           ! To re-use tmp
        
!
!       Careful with 2D arrays...
!
        allocate(tmp2(p%npotd,new_nrval))
        do i=1,p%npotd
           func => p%vdown(i,:)
           call generate_spline(p%r,func,p%nrval,y2,0.0_dp,0.0_dp)
           do j = 1, new_nrval
            call evaluate_spline(p%r,func,y2,p%nrval,new_r(j),tmp2(i,j))
           enddo
           nullify(func)
        enddo
        deallocate(p%vdown)      ! Old data
        p%vdown => tmp2          ! Point to new memory area
        nullify(tmp2)            ! To re-use tmp

        if (p%npotu > 0) allocate(tmp2(p%npotu,new_nrval))
        do i=1,p%npotu         ! Only executed if npotu > 0 ...
           func => p%vup(i,:)
           call generate_spline(p%r,func,p%nrval,y2,0.0_dp,0.0_dp)
           do j = 1, new_nrval
            call evaluate_spline(p%r,func,y2,p%nrval,new_r(j),tmp2(i,j))
           enddo
           nullify(func)
        enddo
        if (p%npotu > 0) then
           deallocate(p%vup  )  ! Old data
           p%vup => tmp2        ! Point to new memory area
           nullify(tmp2)        ! To re-use tmp
        endif

!
!       Now re-set the values
!
        deallocate(p%r)
        p%r => new_r
        p%nrval = new_nrval
        p%nr    = p%nrval - 1
        p%a     = a
        p%b     = b

        deallocate(y2)

      end subroutine pseudo_reparametrize

      end module m_ncps_froyen_reader
