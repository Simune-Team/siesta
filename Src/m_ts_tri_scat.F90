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

  use m_ts_mem_scat, only : UC_expansion_Sigma
  use m_ts_mem_scat, only : UC_expansion_Sigma_Gamma
  use m_ts_mem_scat, only : UC_expansion_Sigma_Bulk

  use m_ts_mem_scat, only : weightDM, weightDMC
  use m_ts_mem_scat, only : read_next_GS, get_scat_region

  implicit none

  public :: calc_GF, calc_GF_Part
  public :: GF_Gamma_GF_Left
  public :: GF_Gamma_GF_Right

  ! Module used, which should be accessible from this module
  public :: weightDM, weightDMC
  public :: read_next_GS, get_scat_region
  public :: UC_expansion_Sigma_Bulk
  public :: UC_expansion_Sigma
  public :: UC_expansion_Sigma_Gamma

  private

  ! Used for BLAS calls
  complex(dp), parameter :: z0  = dcmplx( 0._dp, 0._dp)
  complex(dp), parameter :: z1  = dcmplx( 1._dp, 0._dp)
  complex(dp), parameter :: zm1 = dcmplx(-1._dp, 0._dp)
  complex(dp), parameter :: zi  = dcmplx( 0._dp, 1._dp)
  complex(dp), parameter :: zmi = dcmplx( 0._dp,-1._dp)

contains

! ##################################################################
! ## Calculating Full Greens functions of                         ## 
! ##                                                              ##          
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##                                                              ##
! ##                                                              ##
! ## Completely restructured to be able to handle sparse matrices ##
! ##                                                              ##
! ##                                                              ##
! ##  Modified by Nick Papior Andersen                            ##
! ##################################################################
  subroutine calc_GF(UseBulk, BiasContour, &
       no_u_TS,no_L,no_R, & ! Size of the problem
       SigmaL,SigmaR, & ! Electrode self-energies
       GFinv_tri,GF_tri,ierr)
    
    use intrinsic_missing, only: EYE
    use precision, only: dp
    use class_zTriMat3

    implicit none 

! *********************
! * INPUT variables   *
! *********************
    logical, intent(in) :: UseBulk ! if true self-energy only is input else
!                                    z*S-H-Sigma for bulk is in sfe
    logical, intent(in) :: BiasContour ! if true we also need b11 and b33
    ! Sizes of the different regions...
    integer, intent(in) :: no_u_TS, no_L, no_R
    complex(dp) :: SigmaL(no_L,no_L)   ! Left electrode GF
    complex(dp) :: SigmaR(no_R,no_R)   ! Right electrode GF
    ! Work should already contain Z*S - H
    ! This may seem strange, however, it will clean up this routine extensively
    ! as we dont need to make two different routines for real and complex
    ! Hamiltonian values.
    type(zTriMat3), intent(in out) :: GFinv_tri ! the inverted GF
    type(zTriMat3), intent(in out) :: GF_tri
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
    integer :: i,j,ii


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
    iGf11 => val11(GFinv_tri)
    iGf12 => val12(GFinv_tri)
    iGf21 => val21(GFinv_tri)
    iGf22 => val22(GFinv_tri)
    iGf23 => val23(GFinv_tri)
    iGf32 => val32(GFinv_tri)
    iGf33 => val33(GFinv_tri)
    Gf11  => val11(GF_tri)
    Gf12  => val12(GF_tri)
    Gf21  => val21(GF_tri)
    Gf22  => val22(GF_tri)
    Gf23  => val23(GF_tri)
    Gf32  => val32(GF_tri)
    Gf33  => val33(GF_tri)

    ierr = 0

    ! TODO check that we actually can not use the UseBulk with
    ! Tri-diag... I do not immediately see why not...

    ! Adjust the left electrode part
    ii = 1
    if ( UseBulk ) then 
       ! If we use the bulk electrodes, we will overwrite
       ! the "regular" Hamiltonian elements by the 
       ! Self-energy of the electrode (thereby assuming 
       ! that the electrode states are the perfect bulk system)
       do j = 1, no_L
          do i = 1, no_L
             iGf11(ii) = SigmaL(i,j)
             ii = ii + 1
          end do              !i
       end do                 !j
    else
       ! We don't consider the electrode to be the bulk states,
       ! hence, we will adjust the unit-cell Hamiltonian by the 
       ! self-energy terms.
       do j = 1, no_L
          do i = 1, no_L
             iGf11(ii) = iGF11(ii) - SigmaL(i,j)
             ii = ii + 1
          end do              !i
       end do                 !j
    end if                    !USEBULK

! Adjust the right electrode part
    ii = 1
    if ( UseBulk ) then
       do j = 1 , no_R
          do i = 1 , no_R
             iGf33(ii) = SigmaR(i,j)
             ii = ii + 1
          end do              !i
       end do                 !j
    else
       do j = 1 , no_R
          do i = 1 , no_R
             iGf33(ii) = iGf33(ii) - SigmaR(i,j)
             ii = ii + 1
          end do              !i
       end do                 !j
    end if

! Now we can do MAGIC!!!
    nL = nrows_g_left  (GF_tri)
    nC = nrows_g_center(GF_tri)
    nR = nrows_g_right (GF_tri)

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
    call zgesv(nC,nC,iGf22,nC,ipvt,GF22,nC,ierr)

!*b12 = - x12 * b22 | b'11 * a12 * a''22^-1
    call zgemm('N','N',nL,nC,nC,zm1, Gf21,nL, Gf22,nC,z0, Gf12,nL)

! x21 = a21 * b'11 | a21 * a11^-1
    call zgemm('N','N',nC,nL,nL,z1, iGf21,nC, Gf11,nL,z0, iGf12,nC)

!*b21 = - b22 * x22 | a''22^-1 * a21 * a11^-1
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

! ====================================================================
  END subroutine calc_GF


! ##################################################################
! ## Calculating Full Greens functions of                         ## 
! ##                                                              ##          
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##                                                              ##
! ##                                                              ##
! ## Completely restructured to be able to handle sparse matrices ##
! ##                                                              ##
! ##                                                              ##
! ##  Modified by Nick Papior Andersen                            ##
! ##################################################################
  subroutine calc_GF_Part( &
       no_u_TS,no_L,no_R, & ! Size of the problem
       SigmaL,SigmaR, &             ! Electrode self-energies
       GFinv_tri,GF22,ierr)
    
    use intrinsic_missing, only: EYE
    use precision, only: dp
    use class_zTriMat3
    use alloc

    implicit none 

! *********************
! * INPUT variables   *
! *********************
    ! Sizes of the different regions...
    integer, intent(in) :: no_u_TS, no_L, no_R
    complex(dp), intent(inout) :: SigmaL(no_L,no_L)   ! Left electrode GF
    complex(dp), intent(inout) :: SigmaR(no_R,no_R)   ! Right electrode GF
    ! Work should already contain Z*S - H
    ! This may seem strange, however, it will clean up this routine extensively
    ! as we dont need to make two different routines for real and complex
    ! Hamiltonian values.
    type(zTriMat3), intent(in out) :: GFinv_tri ! the inverted GF (in tri-diagonal form)
    complex(dp), target, intent(inout) :: GF22((no_u_TS-no_R-no_L)**2)
    integer,     intent(out) :: ierr              ! inversion err


! Local variables
    ! TODO change this to a work array...
    complex(dp), pointer :: Gf_OD(:) => null() ! This is for containing off-diagonal products
    complex(dp), pointer :: GFinv(:)
    complex(dp), pointer :: iGf11(:), iGf12(:)
    complex(dp), pointer :: iGf21(:), iGf22(:), iGf23(:)
    complex(dp), pointer ::           iGf32(:), iGf33(:)
    integer :: ipvt(no_u_TS)
    integer :: nL,nC,nR
    integer :: i,j,ii
    logical :: Use_AllocMEM

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE getGF' )
#endif

    call timer('GFT_P',1) 

    ! point the pointers
    GFinv => val(GFinv_tri)

    ! Point to the parts
    iGf11 => val11(GFinv_tri)
    iGf12 => val12(GFinv_tri)
    iGf21 => val21(GFinv_tri)
    iGf22 => val22(GFinv_tri)
    iGf23 => val23(GFinv_tri)
    iGf32 => val32(GFinv_tri)
    iGf33 => val33(GFinv_tri)

    ierr = 0

    ! If we use the bulk electrodes, we will overwrite
    ! the "regular" Hamiltonian elements by the 
    ! Self-energy of the electrode (thereby assuming 
    ! that the electrode states are the perfect bulk system)
    ii = 1
    do j = 1, no_L
       iGf11(ii:ii+no_L-1) = SigmaL(:,j)
       ii = ii + no_L
    end do                 !j

! Adjust the right electrode part
    ii = 1 
    do j = 1 , no_R
       iGf33(ii:ii+no_R-1) = SigmaR(:,j)
       ii = ii + no_R
    end do                 !j

! Now we can do MAGIC!!!
    nL = nrows_g_left  (GFinv_tri)
    nC = nrows_g_center(GFinv_tri)
    nR = nrows_g_right (GFinv_tri)

    ii = max(nL,nR)
    Use_AllocMEM = ii > nC
    ! If the central region is larger than the maximum 
    ! electrode region, we can re-use the memory for 
    ! retaining a temporary calculation array...
    ! This should luckily happen almost always!
    if ( Use_AllocMEM ) then
       call re_alloc(GF_OD,1,ii*nC, &
            name='trimem',routine='transiesta')
    else
       ! This means that the central region is larger
       ! than the electrodes...
       GF_OD => GF22
    end if

! b'11 = a11^-1
    call EYE(nL,SigmaL)
    call zgesv(nL,nL,iGf11,nL,ipvt,SigmaL,nL,ierr)

! b'33 = a33^-1
    call EYE(nR,SigmaR)
    call zgesv(nR,nR,iGf33,nR,ipvt,SigmaR,nR,ierr)

! x12 = b'11 * a12 | a11^-1 * a12
    call zgemm('N','N',nL,nC,nL,z1, SigmaL,nL, iGf12,nL,z0, GF_OD,nL)

! a'22 = a22 - a21 * x12 | a22 - a21 * a11^-1 * a12
    call zgemm('N','N',nC,nC,nL,zm1, iGf21,nC, GF_OD,nL,z1, iGf22,nC)

! x32 = b'33 * a32 | a33^-1 * a32
    call zgemm('N','N',nR,nC,nR,z1, SigmaR,nR, iGf32,nR,z0, GF_OD,nR)

! a''22 = a'22 - a23 * x32 | a22 - a21 * a11^-1 * a12 - a23 * a33^-1 * a32
    call zgemm('N','N',nC,nC,nR,zm1, iGf23,nC, GF_OD,nR,z1, iGf22,nC)
    
!*b22 = a''22^-1
    call EYE(nC,Gf22)
    call zgesv(nC,nC,iGf22,nC,ipvt,GF22,nC,ierr)

    if ( Use_AllocMEM ) then
       call de_alloc(GF_OD, &
            name='trimem',routine='transiesta')
    end if

    call timer('GFT_P',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS getGF' )
#endif

! ====================================================================
  END subroutine calc_GF_Part


  subroutine GF_Gamma_GF_Left(no_L,Gf_tri,Gamma,GGG_tri)

! ======================================================================
!  This routine returns GGG=(-i)*GF.Gamma.GF, where GF is a tri-diagonal
!  matrix and the states
!  corresponds to the (no_L) Left
!  Gamma is a (no_L)x(no_L) matrix.
! ======================================================================
    use precision, only : dp
    use class_zTriMat3

    implicit none

! *********************
! * INPUT variables   *
! *********************
    integer, intent(in) :: no_L ! the size of the Gamma
    ! The Green's function (note that it SHOULD be transposed on entry)
    type(zTriMat3), intent(inout) :: Gf_tri
    ! i (Sigma - Sigma^dagger)/2
    complex(dp),    intent(in) :: Gamma(no_L,no_L)

! *********************
! * OUTPUT variables  *
! *********************
    type(zTriMat3), intent(inout) :: GGG_tri    !GF.GAMMA.GF

    ! local variables
    complex(dp), pointer :: Gf(:), GGG(:), oW(:)
    integer :: nL,nC,nR
    integer :: l_idx, u_idx

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE GFGammaGF' )
#endif

    call timer("GFGGF",1)

    ! should be the same as no_L
    nL = nrows_g_left  (Gf_tri)
    nC = nrows_g_center(Gf_tri)
    nR = nrows_g_right (Gf_tri)

    ! I think this will very rarely happen...
    ! So for now we don't consider this a memory restriction.
    if ( nR * nC * 2 + nR ** 2 < nL * nC ) call die('Your system &
         &has inappropriate sizes')

    ! First we need to point to an empty memory segment of
    ! the tri-diagonal result array...
    GGG => val(GGG_tri)
    l_idx = index(GGG_tri,nL+nC+1,nL+1)
    u_idx = index(GGG_tri,nL+nC+nR,nL+nC+nR)

    ! This is the full GGG 32,23,33 array (note we check the size
    ! requirements so as not to overwrite anything)
    oW => GGG(l_idx:u_idx)

    Gf  => val11(Gf_tri)
    ! \Gamma Gf^\dagger 11
    call zgemm('N','C',nL,nL,nL,z1, Gamma,nL, Gf,nL,z0, oW,nL)

    GGG => val11(GGG_tri)
    ! Gf11 \Gamma Gf^\dagger 11 == GGG 11
    call zgemm('N','N',nL,nL,nL,zmi, Gf,nL, oW,nL,z0, GGG,nL)

    Gf  => val21(Gf_tri) ! size nC x nL
    GGG => val21(GGG_tri) ! size nC x nL
    ! Gf21 \Gamma Gf^\dagger 11 == GGG 21
    call zgemm('N','N',nC,nL,nL,zmi, Gf,nC, oW,nL,z0, GGG,nC)

    ! > We now have GGG :1 <

    Gf => val21(Gf_tri) ! size nC x nL (note we take the conjugate transpose)
    ! \Gamma Gf^\dagger 12
    call zgemm('N','C',nL,nC,nL,z1, Gamma,nL, Gf,nC,z0, oW,nL)

    Gf  => val11(Gf_tri) ! size nL x nL
    GGG => val12(GGG_tri) ! size nL x nC
    ! Gf 11 \Gamma Gf^\dagger 12 == GGG 12
    call zgemm('N','N',nL,nC,nL,zmi, Gf,nL, oW,nL,z0, GGG,nL)

    Gf  => val21(Gf_tri) ! size nC x nL
    GGG => val22(GGG_tri) ! size nC x nC
    ! Gf 21 \Gamma Gf^\dagger 12 == GGG 22
    call zgemm('N','N',nC,nC,nL,zmi, Gf,nC, oW,nL,z0, GGG,nC)

    ! Now we have Gf Gamma Gf^\dagger
    
    ! Only the GGG 1:2,1:2 part!

    ! TODO Check that this is actually not needed...
    ! it resets the value in the right-electrode region which can not
    ! be calculated
    !oW(:) = z0
    
    call timer('GFGGF',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS GFGammaGF' )
#endif

! ====================================================================
  END subroutine GF_Gamma_GF_Left

  subroutine GF_Gamma_GF_Right(no_R,Gf_tri,Gamma,GGG_tri)

! ======================================================================
!  This routine returns GGG=(-i)*GF.Gamma.GF, where GF is a tri-diagonal
!  matrix and the states
!  corresponds to the (no_R) right
!  Gamma is a (no_R)x(no_R) matrix.
! ======================================================================
    use precision, only : dp
    use class_zTriMat3

    implicit none

! *********************
! * INPUT variables   *
! *********************
    integer, intent(in) :: no_R ! the size of the Gamma
    ! The Green's function (note that it SHOULD be transposed on entry)
    type(zTriMat3), intent(inout) :: Gf_tri
    ! i (Sigma - Sigma^dagger)/2
    complex(dp),    intent(in) :: Gamma(no_R,no_R)

! *********************
! * OUTPUT variables  *
! *********************
    type(zTriMat3), intent(inout) :: GGG_tri    !GF.GAMMA.GF

    ! local variables
    complex(dp), pointer :: Gf(:), GGG(:), oW(:)
    integer :: nL,nC,nR
    integer :: l_idx, u_idx

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE GFGammaGF' )
#endif

    call timer("GFGGF",1)

    nL = nrows_g_left  (Gf_tri)
    nC = nrows_g_center(Gf_tri)
    nR = nrows_g_right (Gf_tri)

    ! I think this will very rarely happen...
    ! So for now we don't consider this a memory restriction.
    if ( nL * nC * 2 + nL ** 2 < nR * nC ) call die('Your system &
         &has inappropriate sizes for using the tri-sparse method.')

    ! First we need to point to an empty memory segment of
    ! the tri-diagonal result array...
    GGG   =>  val(GGG_tri)
    l_idx = index(GGG_tri,1,1)
    u_idx = index(GGG_tri,nL,nL+nC)

    ! This is the full GGG 11,21,12 array (note we check the size
    ! requirements so as not to overwrite anything)
    oW => GGG(l_idx:u_idx)

    Gf => val33(Gf_tri)
    ! \Gamma Gf^\dagger 33
    call zgemm('N','C',nR,nR,nR,z1, Gamma,nR, Gf,nR,z0, oW,nR)

    GGG => val33(GGG_tri)
    ! Gf33 \Gamma Gf^\dagger 33 == GGG 33
    call zgemm('N','N',nR,nR,nR,zmi, Gf,nR, oW,nR,z0, GGG,nR)

    Gf  => val23(Gf_tri) ! size nC x nR
    GGG => val23(GGG_tri) ! size nC x nR
    ! Gf23 \Gamma Gf^\dagger 33 == GGG 23
    call zgemm('N','N',nC,nR,nR,zmi, Gf,nC, oW,nR,z0, GGG,nC)

    ! > We now have GGG :3 <

    Gf => val23(Gf_tri) ! size nC x nR
    ! \Gamma Gf^\dagger 32
    call zgemm('N','C',nR,nC,nR,z1, Gamma,nR, Gf,nC,z0, oW,nR)

    !Gf  => val23(Gf_tri) ! size nC x nR
    GGG => val22(GGG_tri) ! size nC x nC
    ! Gf 23 \Gamma Gf^\dagger 32 == GGG 22
    call zgemm('N','N',nC,nC,nR,zmi, Gf,nC, oW,nR,z0, GGG,nC)

    Gf  => val33(Gf_tri) ! size nR x nR
    GGG => val32(GGG_tri) ! size nR x nC
    ! Gf 33 \Gamma Gf^\dagger 32 == GGG 32
    call zgemm('N','N',nR,nC,nR,zmi, Gf,nR, oW,nR,z0, GGG,nR)

    ! Now we have Gf Gamma Gf^\dagger

    ! TODO Check that this is actually not needed...
    ! it resets the value in the right-electrode region which can not
    ! be calculated
    !oW(:) = z0
    
    call timer('GFGGF',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS GFGammaGF' )
#endif

! ====================================================================
  END subroutine GF_Gamma_GF_Right

end module m_ts_tri_scat
