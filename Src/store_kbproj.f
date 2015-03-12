      subroutine  store_proj(l,ikb,ekb,nrc,erefkb,dkbcos,nrval,rofi,proj)

      use aux_proj, only: a, b, is

      integer, parameter :: dp = selected_real_kind(10,100)
!
! Custom routine to process the information about each projector
! This version for classic Siesta
!
      integer, intent(in)  :: l
      integer, intent(in)  :: ikb
      integer, intent(in)  :: nrc
      real(dp), intent(in) :: ekb
      real(dp), intent(in) :: erefkb
      real(dp), intent(in) :: dkbcos
      integer, intent(in)  :: nrval
      real(dp), intent(in) :: rofi(:)
      real(dp), intent(in) :: proj(:)

      real(dp) :: rc
      rc = rofi(nrc)

!  write(filename,"(a,i1,a,i1)")  "KBproj-rl+1.", l, ".", ikb
!  call file_out(nrc,rofi,proj,trim(filename)) 

! Common block with the information about the  KB projectors

      call comKB(is,a,b,rofi,proj,l,ikb,rc,ekb,nrc)
!
      end subroutine store_proj

      subroutine comKB(is,a,b,rofi,proj,l,ikb,rc,ekb,nrc)
C
C  Creates the common block with all the information about the 
C  Kleinman-Bylander projectors.
C  Written by D. Sanchez-Portal, Aug. 1998.
C  Modified by DSP to allow more than one projector per l, July 1999.
C
        implicit none

        integer l, nrc,is, ikb 

        real(dp) rc, ekb, proj(nrmax), a, b,
     .    rofi(nrmax)  
        character(len=40) filename
C
C Internal variables
C
        integer indx, itb, nr, nmax, nmin, nn, il
        real(dp) delt, r, vphi, dy, yp1, ypn
        integer ir, nrckb
C
C Number of points used by polint for the interpolation
C
        integer npoint
        parameter(npoint=4)
        integer lun
C
        rctb(ikb,l,is)=rc 
C
C Interpolation to generate tables with KB projectors
C
!       This assumes that the routine is going to be called
!       with increasing values of l, and sequentially for
!       the different projectors for each l.
!
        indx=0
        do il=0,l-1
          indx=indx+nkblsave(il,is)
        enddo 
        indx=indx+ikb
        if (ikb.gt.nkblsave(l,is)) then
          write(6,'(/,2a,i3,a,i3)')
     .      'comKB: ERROR: Maximum number of KB projectors',
     .      ' for l=',l,' must be', nkblsave(l,is)
          call die 
        endif 

        delt=rc/(dble(ntbmax-1)+1.0d-20) 
        if (delt.gt.deltmax) then
          write(6,'(a)')
     .    'comKB: WARNING It might be a good idea to increase'
          write(6,'(a)')
     .    'comKB: WARNING parameter ntbmax (in file atmparams.f) '
          write(6,'(a,i6)')
     .    'comKB: WARNING to at least ntbmax = ', 
     .    nint(Rc/deltmax)+2
        endif

        table(1,-indx,is)=delt
        table(2,-indx,is)=ekb

        if (debug_kb_generation) then
          nrckb = nint(log(rc/b+1.0d0)/a)+1
          write(filename,"(a,i1)")
     .          trim(global_label) // "-KBProj.", indx
          call file_out(nrckb,rofi,proj,filename)
        endif
!
        do itb=1,ntbmax-1
          r=delt*(itb-1)
          nr=nint(log(r/b+1.0d0)/a)+1
          nmin=max(1,nr-npoint)
          nmax=min(nrc,nr+npoint)
          nn=nmax-nmin+1
          call polint(rofi(nmin),proj(nmin),nn,r,vphi,dy)
          table(itb+2,-indx,is)=vphi
        enddo 
        table(ntbmax+2,-indx,is)=0.0d0
C
C Table with the second derivative
C
        yp1=0.0d0
        ypn=huge(1.d0)

        call spline(delt, table(3:2+ntbmax,-indx,is), ntbmax,
     &              yp1, ypn, tab2(1:ntbmax,-indx,is) )

        end subroutine comkb


