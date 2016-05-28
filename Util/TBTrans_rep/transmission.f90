! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
! ##################################################################
! ##                                                              ##
! ##        This file is part of the TranSIESTA package.          ##
! ##                                                              ##
! ##           Calculating Full Greens functions of               ## 
! ##              contact coupled to electrodes                   ##
! ##                                                              ##          
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##                                                              ##
! ## Changed to F90 by Jose-Luis Mozos, jlm@kanigo.icmab.es       ##
! ## Optimized by V.V. Maslyuk (2007), maslyuk@gmail.com          ##
! ##################################################################
!
subroutine transmission(UseBulk,nou,Hk,Sk,noD,noL,SFEL,noR,SFER, &
     ZEnergy,GF,GFRGF,TotTrans,tt,ierr)

  use precision,       only : dp
  use sys,             only : die
  use m_ts_aux_rout,   only : csolveg

  implicit none

!     INPUT:
  logical,    intent(in) :: UseBulk          ! true: self-energy only is input else z*S-H-Sigma for bulk is in sfe
  integer,    intent(in) :: nou              ! no. states
  integer,    intent(in) :: noD              ! no. states in contact region
  complex(dp),intent(in) :: Hk(nou,nou)      ! Global Hamiltonian
  complex(dp),intent(in) :: Sk(nou,nou)      ! Global Overlap Matrix
  complex(dp),intent(in) :: ZEnergy          ! contour energy
  integer,    intent(in) :: noL,noR          ! number of orbitals in the left/right region
  complex(dp),intent(in) :: SFEL(noL,noL)    ! Self-energy of Left electrode
  complex(dp),intent(in) :: SFER(noR,noR)    ! Self-energy of Right electrode
! ======================================================================
!     OUTPUT:
  real(dp)                :: TotTrans        ! sum of trace of tt
  complex(dp),intent(out) :: tt(noD,noD)     ! 
  complex(dp),intent(out) :: GF(noD,noD)     ! 1/(Sc*Energy - H)
  complex(dp),intent(out) :: GFRGF(noD,noD)  ! GF^dagger.GammaR.GF
  integer                 :: ierr            ! error in inversion
! ======================================================================
!     Helpers, tempos ...
  complex(dp),dimension(noD,noD)         :: SigmaL, SigmaR !
  complex(dp),dimension(:,:),allocatable :: Hc             !
  complex(dp),dimension(:,:),allocatable :: U,iU           !
  integer,    dimension(:),  allocatable :: ipvt           ! pivoting vector for matrix inv.
  real(dp)                               :: Energy         ! real part of energy
  complex(dp)                            :: csum           !
  integer i,j,ii,jj
! ======================================================================
!     BEGIN
! ======================================================================
  call timer('transmission',1)
  if (noD.ne.nou-(noL+noR)) call die('transmission : noD/nou are wrong') 
  Energy = dreal(ZEnergy)

  allocate(ipvt(nou)) 
  call memory('A','I',nou,'transmis')
  allocate( U(nou,nou))
  call memory('A','Z',nou*nou,'transmis')
  allocate(iU(nou,nou))
  call memory('A','Z',nou*nou,'transmis')

  ipvt = 0
  ierr = 0
!
!     Left Self-energy
!
  U =dcmplx(0.0_dp,0.0_dp)
  iU=dcmplx(0.0_dp,0.0_dp)
  do j = 1, noL
     do i = 1, noL
        if (UseBulk) then
           U(i,j) = SFEL(i,j)                                    
        else
           U(i,j) = Energy*Sk(i,j)-Hk(i,j)-SFEL(i,j)
        end if
     end do
  end do

  do j = 1, noD
     do i = 1, noL
        iU(i,j) = Energy*Sk(i,j+noL)-Hk(i,j+noL)
     end do
  end do

  call csolveg(noL,noD,U(1:noL,1:noL), &
                     iU(1:noL,1:noD),ipvt,ierr)
  if ( ierr /= 0 ) return

  do j = 1, noL
     do i = 1, noD
        U(i,j) = Energy*Sk(i+noL,j)-Hk(i+noL,j)
     end do
  end do

  call zgemm('N','N',noD,noD,noL,(1.0_dp,0.0_dp), &
            U(1:noD,1:noL),noD, &
           iU(1:noL,1:noD),noL,(0.0_dp,0.0_dp), &
       SigmaL(1:noD,1:noD),noD)

!
!     Right Self-energy
!
  U =dcmplx(0.0_dp,0.0_dp)
  iU=dcmplx(0.0_dp,0.0_dp)

  do j = 1, noR
     jj = j + noD+noL
     do i = 1, noR
        ii = i + noD+noL
        if (UseBulk) then
           U(i,j) = SFER(i,j)
        else
           U(i,j) = Energy*Sk(ii,jj)-Hk(ii,jj)-SFER(i,j)
        endif
     end do
  end do

  do j = 1, noD
     jj=j+noL
     do i = 1, noR
        ii = i + noD+noL
        iU(i,j) = Energy*Sk(ii,jj)-Hk(ii,jj)
     end do
  end do

  call csolveg(noR,noD,U(1:noR,1:noR), &
                     iU(1:noR,1:noD),ipvt,ierr)
  if ( ierr /= 0 ) return

  do j = 1, noR
     do i = 1, noD
        U(i,j) = Energy*Sk(i+noL,j+noL+noD) - Hk(i+noL,j+noL+noD)
     end do
  end do

  call zgemm('N','N',noD,noD,noR,(1.0_dp,0.0_dp), &
            U(1:noD,1:noR),noD, &
           iU(1:noR,1:noD),noR,(0.0_dp,0.0_dp), &
       SigmaR(1:noD,1:noD),noD)

!      
!     Calculating inverse matrix
!       
  allocate(Hc(noD,noD))
  call memory('A','Z',noD*noD,'transmis')
  
  do j = 1, noD
     jj=j+noL
     do i = 1, noD
        ii=i+noL
        Hc(i,j) = Energy*Sk(ii,jj)-Hk(ii,jj) - SigmaR(i,j) - SigmaL(i,j)
        GF(i,j) = dcmplx(0.0_dp,0.0_dp)
     end do
     GF(j,j) = dcmplx(1.0_dp,0.0_dp)
  end do

  call csolveg(noD,noD,Hc(1:noD,1:noD), &
                      GF(1:noD,1:noD),ipvt,ierr)
  if ( ierr /= 0 ) return

  call memory('D','Z',noD*noD,'transmis')
  deallocate(Hc)

! 
!     GammaL -> [U], GammaR -> [iU]
!
  do j = 1, noD
     do i = 1, noD
         U(j,i) = SigmaL(j,i) - dconjg(SigmaL(i,j))
        iU(j,i) = SigmaR(j,i) - dconjg(SigmaR(i,j))
     end do
  end do

! --------------------------------------------------------------------
!      TRANSMISSION
!      tt  =  GammaL.(GF*.GammaR.GF):
! --------------------------------------------------------------------

  allocate(Hc(nou,nou))
  call memory('A','Z',nou*nou,'transmis')

!     -GammaR.GF == -[iU].GF -> [Hc]
  call zgemm('N','N',noD,noD,noD,dcmplx(-1.0_dp,0.0_dp), &
       iU,nou, &
       GF,noD,dcmplx(0.0_dp,0.0_dp), &
       Hc,nou)

!     GF*.GammaR.GF == GF*.[Hc] -> [iU]
  call zgemm('C','N',noD,noD,noD,dcmplx(1._dp,0._dp), &
       GF,noD, &
       Hc,nou,dcmplx(0.0_dp,0.0_dp), &
       iU,nou)

  GFRGF = -dcmplx(0.5_dp,0.0_dp) * iU(1:noD,1:noD)

  call memory('D','Z',nou*nou,'transmis')
  deallocate(Hc)

!     GammaL.GF*.GammaR.GF == [U].[iU] -> tt
  call zgemm('N','N',noD,noD,noD,dcmplx(1.0_dp,0.0_dp), &
        U,nou, &
       iU,nou, &
       dcmplx(0.0_dp,0.0_dp), &
       tt,noD)

  call memory('D','Z',size(U),'transmis')
  deallocate(U)
  call memory('D','Z',size(iU),'transmis')
  deallocate(iU)
  call memory('D','I',size(ipvt),'transmis')
  deallocate(ipvt)

!
!     Trace of tt:
!
  csum=dcmplx(0.0_dp,0.0_dp)
  do i = 1 , noD
     csum = csum + tt(i,i)
  end do
  TotTrans = dreal(csum)

  call timer('transmission',2)

end subroutine transmission

