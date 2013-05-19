!
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
!
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

module m_ts_tri_scat

  use precision, only : dp

  implicit none

  private

  public :: calc_GF
  public :: calc_GF_Bias
  public :: calc_GF_Part
  public :: GF_Gamma_GF_Left
  public :: GF_Gamma_GF_Right

  ! Used for BLAS calls (local variables)
  complex(dp), parameter :: z0  = dcmplx( 0._dp, 0._dp)
  complex(dp), parameter :: z1  = dcmplx( 1._dp, 0._dp)
  complex(dp), parameter :: zm1 = dcmplx(-1._dp, 0._dp)
  complex(dp), parameter :: zi  = dcmplx( 0._dp, 1._dp)

contains

! ##################################################################
! ## Calculating Full Greens functions of                         ## 
! ##                                                              ##          
! ##                            By                                ##
! ##            Nick Papior Andersen, nickpapior@gmail.com        ##
! ##                                                              ##
! ##                                                              ##
! ## Completely restructured to be able to handle the             ##
! ## tri-diagonal structure of the Green's function.              ##
! ##                                                              ##
! ##################################################################
  subroutine calc_GF(BiasContour, &
       no_u_TS, GFinv_tri,GF_tri,ierr)
    
    use intrinsic_missing, only: EYE
    use class_zTriMat

    implicit none 

! *********************
! * INPUT variables   *
! *********************
    logical, intent(in) :: BiasContour ! if true we also need b11 and b33
    ! Sizes of the different regions...
    integer, intent(in) :: no_u_TS
    ! Work should already contain Z*S - H
    ! This may seem strange, however, it will clean up this routine extensively
    ! as we dont need to make two different routines for real and complex
    ! Hamiltonian values.
    type(zTriMat), intent(in out) :: GFinv_tri ! the inverted GF
    type(zTriMat), intent(in out) :: GF_tri
    integer,     intent(out) :: ierr              !inversion err

! Local variables
    complex(dp), pointer :: GF(:), GFinv(:)
    complex(dp), pointer :: iGf11(:), iGf12(:)
    complex(dp), pointer :: iGf21(:), iGf22(:), iGf23(:)
    complex(dp), pointer ::           iGf32(:), iGf33(:)
    complex(dp), pointer :: Gf11(:), Gf12(:)
    complex(dp), pointer :: Gf21(:), Gf22(:), Gf23(:)
    complex(dp), pointer ::          Gf32(:), Gf33(:)
    integer :: ipvt(no_u_TS)

    integer :: nL,nC,nR

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE getGF' )
#endif

    call timer('GFT',1) 
    
    ! point the pointers
    GFinv => val(GFinv_tri)
    GF    => val(GF_tri)
    if ( size(GFinv) /= size(GF) ) &
         call die('Size of tri-diagonal matrices are not &
         &consistent.')

    ! Point to the correct tri-diagonal parts
    iGf11 => val(GFinv_tri,1,1)
    iGf12 => val(GFinv_tri,1,2)
    iGf21 => val(GFinv_tri,2,1)
    iGf22 => val(GFinv_tri,2,2)
    iGf23 => val(GFinv_tri,2,3)
    iGf32 => val(GFinv_tri,3,2)
    iGf33 => val(GFinv_tri,3,3)
    Gf11  => val(GF_tri,1,1)
    Gf12  => val(GF_tri,1,2)
    Gf21  => val(GF_tri,2,1)
    Gf22  => val(GF_tri,2,2)
    Gf23  => val(GF_tri,2,3)
    Gf32  => val(GF_tri,3,2)
    Gf33  => val(GF_tri,3,3)

    ierr = 0

! Now we can do MAGIC!!!
    nL = nrows_g(GF_tri,1)
    nC = nrows_g(GF_tri,2)
    nR = nrows_g(GF_tri,3)

! b'11 = a11^-1
    call EYE(nL,Gf11)
    call zgesv(nL,nL,iGf11,nL,ipvt,Gf11,nL,ierr)

! b'33 = a33^-1
    call EYE(nR,Gf33)
    call zgesv(nR,nR,iGf33,nR,ipvt,Gf33,nR,ierr)

! x12 = b'11 * a12 | a11^-1 * a12
    call zgemm('N','N',nL,nC,nL,z1, Gf11,nL, iGf12,nL,z0, Gf21,nL)

! a'22 = a22 - a21 * x12 | a22 - a21 * a11^-1 * a12
    call zgemm('N','N',nC,nC,nL,zm1, iGf21,nC, Gf21,nL,z1, iGf22,nC)

! x32 = b'33 * a32 | a33^-1 * a32
    call zgemm('N','N',nR,nC,nR,z1, Gf33,nR, iGf32,nR,z0, Gf23,nR)

! a''22 = a'22 - a23 * x32 | a22 - a21 * a11^-1 * a12 - a23 * a33^-1 * a32
    call zgemm('N','N',nC,nC,nR,zm1, iGf23,nC, Gf23,nR,z1, iGf22,nC)
    
!*b22 = a''22^-1
    call EYE(nC,Gf22)
    call zgesv(nC,nC,iGf22,nC,ipvt,Gf22,nC,ierr)

!*b12 = - x12 * b22 | b'11 * a12 * a''22^-1
    call zgemm('N','N',nL,nC,nC,zm1, Gf21,nL, Gf22,nC,z0, Gf12,nL)

! x21 = a21 * b'11 | a21 * a11^-1
    call zgemm('N','N',nC,nL,nL,z1, iGf21,nC, Gf11,nL,z0, iGf12,nC)

!*b21 = - b22 * x21 | a''22^-1 * a21 * a11^-1
    call zgemm('N','N',nC,nL,nC,zm1, Gf22,nC, iGf12,nC,z0, Gf21,nC)

    if ( BiasContour ) then
       ! We only need the *full* Green's function for the Bias points...
!*b11 = b'11 - b12 * x21 | a11^-1 - a11^-1 * a12 * a''22^-1 * a21 * a11^-1
       call zgemm('N','N',nL,nL,nC,zm1, Gf12,nL, iGf12,nC,z1, Gf11,nL)
    end if

!*b32 = - x32 * b22 | a33^-1 * a32 * a''22^-1
    call zgemm('N','N',nR,nC,nC,zm1, Gf23,nR, Gf22,nC,z0, Gf32,nR)

! x23 = a23 * b'33 | a23 * a33^-1
    call zgemm('N','N',nC,nR,nR,z1 ,iGf23,nC, Gf33,nR,z0, iGf32,nC)

!*b23 = - b22 * x23 | a''22^-1 * a23 * a33^-1
    call zgemm('N','N',nC,nR,nC,zm1, Gf22,nC, iGf32,nC,z0, Gf23,nC)

    if ( BiasContour ) then
       ! We only need the *full* Green's function for the Bias points...
!*b33 = b'33 - b32 * x23 | a33^-1 - a''22^-1 * a23 * a33^-1 * a''22^-1 * a23 * a33^-1
       call zgemm('N','N',nR,nR,nC,zm1, Gf32,nR, iGf32,nC,z1, Gf33,nR)
    end if

    ! We are done with the inversion of the tri-diagonal case
       
    call timer('GFT',2)  

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS getGF' )
#endif

  end subroutine calc_GF


! ##################################################################
! ## Calculating the Gf(:,[13]) part of the Greens functions of   ## 
! ##                                                              ##          
! ##                            By                                ##
! ##           Nick Papior Andersen, nickpapior@gmail.com         ##
! ##                                                              ##
! ##################################################################
  subroutine calc_GF_Bias(no_u_TS,Gfinv_tri,GF_tri)
    
    use intrinsic_missing, only: EYE
    use class_zTriMat

    implicit none 

! *********************
! * INPUT variables   *
! *********************
    ! Sizes of the different regions...
    integer, intent(in) :: no_u_TS
    ! Work should already contain Z*S - H (and the self-energies)
    ! This may seem strange, however, it will clean up this routine extensively
    ! as we dont need to make two different routines for real and complex
    ! Hamiltonian values.
    type(zTriMat), intent(in out) :: GFinv_tri ! the inverted GF
    type(zTriMat), intent(in out) :: GF_tri

! Local variables
    complex(dp), pointer :: GF(:), GFinv(:)
    complex(dp), pointer :: iGf11(:), iGf12(:)
    complex(dp), pointer :: iGf21(:), iGf22(:), iGf23(:)
    complex(dp), pointer ::           iGf32(:), iGf33(:)
    complex(dp), pointer :: Gf11(:), Gf12(:)
    complex(dp), pointer :: Gf21(:), Gf22(:), Gf23(:)
    complex(dp), pointer ::          Gf32(:), Gf33(:)
    integer :: ipvt(no_u_TS)

    integer :: ierr ! inversion err
    integer :: nL,nC,nR


#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE getGF' )
#endif

    call timer('GFT',1) 
    
    ! point the pointers
    GFinv => val(GFinv_tri)
    GF    => val(GF_tri)
    if ( size(GFinv) /= size(GF) ) &
         call die('Size of tri-diagonal matrices are not &
         &consistent.')

    ! Point to the correct tri-diagonal parts
    iGf11 => val(GFinv_tri,1,1)
    iGf12 => val(GFinv_tri,1,2)
    iGf21 => val(GFinv_tri,2,1)
    iGf22 => val(GFinv_tri,2,2)
    iGf23 => val(GFinv_tri,2,3)
    iGf32 => val(GFinv_tri,3,2)
    iGf33 => val(GFinv_tri,3,3)
    Gf11  => val(GF_tri,1,1)
    Gf12  => val(GF_tri,1,2)
    Gf21  => val(GF_tri,2,1)
    Gf22  => val(GF_tri,2,2)
    Gf23  => val(GF_tri,2,3)
    Gf32  => val(GF_tri,3,2)
    Gf33  => val(GF_tri,3,3)

    ierr = 0

! Now we can do MAGIC!!!
    nL = nrows_g(GF_tri,1)
    nC = nrows_g(GF_tri,2)
    nR = nrows_g(GF_tri,3)

! copy over A3 (we have required that nR <= nC)
    Gf22(1:nR**2) = iGf33(:)

! copy over B2
    Gf32(:) = iGf32(:)

! solve A3x=B2  (X2/C3)
    call zgesv(nR,nC,Gf22,nR,ipvt,Gf32,nR,ierr)
    if ( ierr /= 0 ) call die('error: A3x=B2')

! copy over A2
    Gf22(:) = iGf22(:)
    
! calculate A2-C3 * X2/C3 (from above)
    call zgemm('N','N',nC,nC,nR,zm1, iGf23,nC, Gf32,nR,z1, Gf22,nC)

! copy over B1
    Gf12(:) = iGf21(:)
    
! calculate (A2-X2)^-1B1  (X1/C2)
    call zgesv(nC,nL,Gf22,nC,ipvt,Gf12,nC,ierr)
    if ( ierr /= 0 ) call die('error: (A2-X2)x=B1')

! copy over A1 (we have required that nL <= nC)
    Gf22(1:nL**2) = iGf11(:)

! calculate A1 - X1 = A1-C2*(A2-X2)^-1B1 
    call zgemm('N','N',nL,nL,nC,zm1, iGf12,nL, Gf12,nC,z1, Gf22,nL)

! calculate inv(Gf)11
    call EYE(nL,Gf11)
    call zgesv(nL,nL,Gf22,nL,ipvt,Gf11,nL,ierr)
    if ( ierr /= 0 ) call die('error: (A1-X1)x=Gf11^-1')

! calculate inv(Gf)21
    call zgemm('N','N',nC,nL,nL,zm1, Gf12,nC, Gf11,nL,z0, Gf21,nC)

! calculate inv(Gf)31 
! (untraditionally this we save in Gf12 as the tri-diagonal is not made for full Gf)
    call zgemm('N','N',nR,nL,nC,zm1, Gf32,nR, Gf21,nC,z0, Gf12,nR)

! * Now we have calculated the Gf in the left column *

! We now move to the calculation of Gf in the right column
! Now we can overwrite Gf_inv without problems (we dont need it anymore)

! solve A1x=C2  (Y2/B1)
    call zgesv(nL,nC,iGf11,nL,ipvt,iGf12,nL,ierr)
    if ( ierr /= 0 ) call die('error: A1x=C2')

! calculate A2-B1 * Y2/B1 (from above)
    call zgemm('N','N',nC,nC,nL,zm1, iGf21,nC, iGf12,nL,z1, iGf22,nC)

! calculate [A2-Y2]^-1C3  (Y3/B2)
    call zgesv(nC,nR,iGf22,nC,ipvt,iGf23,nC,ierr)
    if ( ierr /= 0 ) call die('error: (A2-Y2)x=C3')

! calculate A3-Y3 = A3-B2*Y3/B2 = A3-B2 [A2-Y2]^-1C3
    call zgemm('N','N',nR,nR,nC,zm1, iGf32,nR, iGf23,nC, z1, iGf33,nR)

! calculate inv(Gf)33
    call EYE(nR,Gf33)
    call zgesv(nR,nR,iGf33,nR,ipvt,Gf33,nR,ierr)
    if ( ierr /= 0 ) call die('error: (A3-Y3)x=Gf33^-1')

! calculate inv(Gf)23
    call zgemm('N','N',nC,nR,nR,zm1, iGf23,nC, Gf33,nR,z0, Gf23,nC)

! calculate inv(Gf)13
! (untraditionally this we save in Gf32 as the tri-diagonal is not made for full Gf)
    call zgemm('N','N',nL,nR,nC,zm1, iGf12,nL, Gf23,nC,z0, Gf32,nL)

    call timer('GFT',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS getGF' )
#endif

  end subroutine calc_GF_Bias


! ##################################################################
! ## Calculating of part of the Greens functions                  ## 
! ##                                                              ##          
! ##                Nick Papior Andersen, nickpapior@gmail.com    ##
! ##                                                              ##
! ##################################################################
  subroutine calc_GF_Part(no_u_TS, no_L, no_R, GFinv_tri,GF22,ierr)
    
    use intrinsic_missing, only: EYE
    use class_zTriMat
    use alloc

    implicit none 

! *********************
! * INPUT variables   *
! *********************
    ! Sizes of the different regions...
    integer, intent(in) :: no_u_TS, no_L, no_R
    ! Work should already contain Z*S - H
    ! This may seem strange, however, it will clean up this routine extensively
    ! as we dont need to make two different routines for real and complex
    ! Hamiltonian values.
    type(zTriMat), intent(in out) :: GFinv_tri ! the inverted GF (in tri-diagonal form)
    complex(dp), target, intent(inout) :: GF22((no_u_TS-no_R-no_L)**2)
    integer,     intent(out) :: ierr              ! inversion err


! Local variables
    complex(dp), pointer :: GFinv(:)
    complex(dp), pointer :: iGf11(:), iGf12(:)
    complex(dp), pointer :: iGf21(:), iGf22(:), iGf23(:)
    complex(dp), pointer ::           iGf32(:), iGf33(:)
    integer :: ipvt(no_u_TS)
    integer :: nL,nC,nR

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE getGF' )
#endif

    call timer('GFT_P',1) 

    ! point the pointers
    GFinv => val(GFinv_tri)

    ! Point to the parts
    iGf11 => val(GFinv_tri,1,1)
    iGf12 => val(GFinv_tri,1,2)
    iGf21 => val(GFinv_tri,2,1)
    iGf22 => val(GFinv_tri,2,2)
    iGf23 => val(GFinv_tri,2,3)
    iGf32 => val(GFinv_tri,3,2)
    iGf33 => val(GFinv_tri,3,3)

    ierr = 0

! Now we can do MAGIC!!!
    nL = nrows_g(GFinv_tri,1)
    nC = nrows_g(GFinv_tri,2)
    nR = nrows_g(GFinv_tri,3)

! x12 = a11^-1 * a12
    call zgesv(nL,nC,iGf11,nL,ipvt,iGf12,nL,ierr)

! a'22 = a22 - a21 * x12 | a22 - a21 * a11^-1 * a12
    call zgemm('N','N',nC,nC,nL,zm1, iGf21,nC, iGf12,nL,z1, iGf22,nC)

! x32 = a33^-1 * a32
    call zgesv(nR,nC,iGf33,nR,ipvt,iGf32,nR,ierr)

! a''22 = a'22 - a23 * x32 | a22 - a21 * a11^-1 * a12 - a23 * a33^-1 * a32
    call zgemm('N','N',nC,nC,nR,zm1, iGf23,nC, iGf32,nR,z1, iGf22,nC)
    
!*b22 = a''22^-1
    call EYE(nC,Gf22)
    call zgesv(nC,nC,iGf22,nC,ipvt,GF22,nC,ierr)

    call timer('GFT_P',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS getGF' )
#endif

  end subroutine calc_GF_Part

  subroutine GF_Gamma_GF_Left(no_L,Gf_tri,GammaT,GGG_tri)

! ======================================================================
!  This routine returns GGG=GF.Gamma.GF^\dagger, where GF is a tri-diagonal
!  matrix and the states
!  corresponds to the (no_L) Left
!  Gamma is a (no_L)x(no_L) matrix.
! ======================================================================
    use class_zTriMat

    implicit none

! *********************
! * INPUT variables   *
! *********************
    integer, intent(in) :: no_L ! the size of the Gamma
    ! The Green's function (note that it SHOULD be transposed on entry)
    type(zTriMat), intent(inout) :: Gf_tri
    ! i (Sigma - Sigma^dagger)/2
    complex(dp),    intent(in) :: GammaT(no_L,no_L)

! *********************
! * OUTPUT variables  *
! *********************
    type(zTriMat), intent(inout) :: GGG_tri    !GF.GAMMA.GF

    ! local variables
    complex(dp), pointer :: Gf(:), GGG(:), oW(:)
    integer :: nL,nC,nR

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE GFGammaGF' )
#endif

    call timer("GFGGF",1)

    ! should be the same as no_L
    nL = nrows_g(Gf_tri,1)
    nC = nrows_g(Gf_tri,2)
    nR = nrows_g(Gf_tri,3)

    ! This is the full Gf22 array
    oW => val(Gf_tri,2,2)

    Gf  => val(Gf_tri,1,1)
    ! \Gamma Gf^\dagger 11
    call zgemm('T','C',nL,nL,nL,z1, GammaT,nL, Gf,nL,z0, oW,nL)

    ! This is actually not necessary... (we dont have anything to 
    ! update here)
!    GGG => val(GGG_tri,1,1)
!    ! Gf11 \Gamma Gf^\dagger 11 == GGG 11
!    call zgemm('N','N',nL,nL,nL,z1, Gf,nL, oW,nL,z0, GGG,nL)

    Gf  => val(Gf_tri,2,1) ! size nC x nL
    GGG => val(GGG_tri,2,1) ! size nC x nL
    ! Gf21 \Gamma Gf^\dagger 11 == GGG 21
    call zgemm('N','N',nC,nL,nL,z1, Gf,nC, oW,nL,z0, GGG,nC)

    ! > We now have GGG 21 <

    Gf => val(Gf_tri,2,1) ! size nC x nL (note we take the conjugate transpose)
    ! \Gamma Gf^\dagger 12
    call zgemm('T','C',nL,nC,nL,z1, GammaT,nL, Gf,nC,z0, oW,nL)

    Gf  => val(Gf_tri,1,1) ! size nL x nL
    GGG => val(GGG_tri,1,2) ! size nL x nC
    ! Gf 11 \Gamma Gf^\dagger 12 == GGG 12
    call zgemm('N','N',nL,nC,nL,z1, Gf,nL, oW,nL,z0, GGG,nL)

    Gf  => val(Gf_tri,2,1) ! size nC x nL
    GGG => val(GGG_tri,2,2) ! size nC x nC
    ! Gf 21 \Gamma Gf^\dagger 12 == GGG 22
    call zgemm('N','N',nC,nC,nL,z1, Gf,nC, oW,nL,z0, GGG,nC)

    ! NOTICE that Gf contains Gf31 in Gf12 as Gf:2 is not used
    ! for anything
    Gf  => val(Gf_tri,1,2) ! size nR x nL
    GGG => val(GGG_tri,3,2) ! size nR x nC
    ! Gf 31 \Gamma Gf^\dagger 12 == GGG 32
    call zgemm('N','N',nR,nC,nL,z1, Gf,nR, oW,nL,z0, GGG,nR)

    ! > We now have GGG :2 <
    
    ! NOTICE that Gf contains Gf31 in Gf12 as Gf:2 is not used
    ! for anything
    Gf => val(Gf_tri,1,2) ! size nR x nL (note we take the conjugate transpose)
    ! \Gamma Gf^\dagger 13
    call zgemm('T','C',nL,nR,nL,z1, GammaT,nL, Gf,nR,z0, oW,nL)

    Gf  => val(Gf_tri,2,1) ! size nC x nL
    GGG => val(GGG_tri,2,3) ! size nC x nR
    ! Gf 21 \Gamma Gf^\dagger 13 == GGG 23
    call zgemm('N','N',nC,nR,nL,z1, Gf,nC, oW,nL,z0, GGG,nC)

    ! This is actually not necessary... (we dont have anything to 
    ! update here)
    ! NOTICE that Gf contains Gf31 in Gf12 as Gf:2 is not used
    ! for anything
!    Gf  => val(Gf_tri,1,2) ! size nR x nL
!    GGG => val(GGG_tri,3,3) ! size nR x nR
    ! Gf 31 \Gamma Gf^\dagger 13 == GGG 33
!    call zgemm('N','N',nR,nR,nL,z1, Gf,nR, oW,nL,z0, GGG,nR)

    ! > We now have GGG 23 <
    
    call timer('GFGGF',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS GFGammaGF' )
#endif

  end subroutine GF_Gamma_GF_Left


  subroutine GF_Gamma_GF_Right(no_R,Gf_tri,GammaT,GGG_tri)

! ======================================================================
!  This routine returns GGG=GF.Gamma.GF, where GF is a tri-diagonal
!  matrix and the states
!  corresponds to the (no_R) right
!  Gamma is a (no_R)x(no_R) matrix.
! ======================================================================
    use class_zTriMat

    implicit none

! *********************
! * INPUT variables   *
! *********************
    integer, intent(in) :: no_R ! the size of the Gamma
    ! The Green's function (note that it SHOULD be transposed on entry)
    type(zTriMat), intent(inout) :: Gf_tri
    ! i (Sigma - Sigma^dagger)/2
    complex(dp),    intent(in) :: GammaT(no_R,no_R)

! *********************
! * OUTPUT variables  *
! *********************
    type(zTriMat), intent(inout) :: GGG_tri    !GF.GAMMA.GF

    ! local variables
    complex(dp), pointer :: Gf(:), GGG(:), oW(:)
    integer :: nL,nC,nR

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE GFGammaGF' )
#endif

    call timer("GFGGF",1)

    nL = nrows_g(Gf_tri,1)
    nC = nrows_g(Gf_tri,2)
    nR = nrows_g(Gf_tri,3)

    ! We use Gf22 for the calculation array (notice that Gf :2 is not used
    ! for anything, hence this is safe)
    oW => val(Gf_tri,2,2)

    Gf => val(Gf_tri,3,3)
    ! \Gamma Gf^\dagger 33
    call zgemm('T','C',nR,nR,nR,z1, GammaT,nR, Gf,nR,z0, oW,nR)

    ! This is actually not necessary... (we dont have anything to 
    ! update here)
!    GGG => val(GGG_tri,3,3)
!    ! Gf33 \Gamma Gf^\dagger 33 == GGG 33
!    call zgemm('N','N',nR,nR,nR,z1, Gf,nR, oW,nR,z0, GGG,nR)

    Gf  => val(Gf_tri,2,3) ! size nC x nR
    GGG => val(GGG_tri,2,3) ! size nC x nR
    ! Gf23 \Gamma Gf^\dagger 33 == GGG 23
    call zgemm('N','N',nC,nR,nR,z1, Gf,nC, oW,nR,z0, GGG,nC)

    ! > We now have GGG 23 <

    Gf => val(Gf_tri,2,3) ! size nC x nR
    ! \Gamma Gf^\dagger 32
    call zgemm('T','C',nR,nC,nR,z1, GammaT,nR, Gf,nC,z0, oW,nR)

    ! NOTICE that Gf contains Gf13 in Gf32 as Gf:2 is not used
    ! for anything
    Gf  => val(Gf_tri,3,2) ! size nL x nR
    GGG => val(GGG_tri,1,2) ! size nL x nC
    ! Gf 13 \Gamma Gf^\dagger 32 == GGG 12
    call zgemm('N','N',nL,nC,nR,z1, Gf,nL, oW,nR,z0, GGG,nL)

    Gf  => val(Gf_tri,2,3) ! size nC x nR
    GGG => val(GGG_tri,2,2) ! size nC x nC
    ! Gf 23 \Gamma Gf^\dagger 32 == GGG 22
    call zgemm('N','N',nC,nC,nR,z1, Gf,nC, oW,nR,z0, GGG,nC)

    Gf  => val(Gf_tri,3,3) ! size nR x nR
    GGG => val(GGG_tri,3,2) ! size nR x nC
    ! Gf 33 \Gamma Gf^\dagger 32 == GGG 32
    call zgemm('N','N',nR,nC,nR,z1, Gf,nR, oW,nR,z0, GGG,nR)

    ! > We now have GGG :2 <

    ! NOTICE that Gf contains Gf13 in Gf32 as Gf:2 is not used
    ! for anything
    Gf => val(Gf_tri,3,2) ! size nL x nR (note we take the conjugate transpose)
    ! \Gamma Gf^\dagger 31
    call zgemm('T','C',nR,nL,nR,z1, GammaT,nR, Gf,nL,z0, oW,nR)

    ! This is actually not necessary... (we dont have anything to 
    ! update here)
    ! NOTICE that Gf contains Gf13 in Gf32 as Gf:2 is not used
    ! for anything
!    Gf  => val(Gf_tri,3,2) ! size nL x nR
!    GGG => val(GGG_tri,1,1) ! size nL x nL
    ! Gf 13 \Gamma Gf^\dagger 31 == GGG 11
!    call zgemm('N','N',nL,nL,nR,z1, Gf,nL, oW,nR,z0, GGG,nL)

    Gf  => val(Gf_tri,2,3) ! size nC x nR
    GGG => val(GGG_tri,2,1) ! size nC x nL
    ! Gf 23 \Gamma Gf^\dagger 31 == GGG 21
    call zgemm('N','N',nC,nL,nR,z1, Gf,nC, oW,nR,z0, GGG,nC)

    ! > We now have GGG 21 <

    call timer('GFGGF',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS GFGammaGF' )
#endif

  end subroutine GF_Gamma_GF_Right

end module m_ts_tri_scat
