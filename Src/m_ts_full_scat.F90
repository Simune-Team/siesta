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

module m_ts_full_scat

  use precision, only : dp

  implicit none

  private

  public :: calc_GF
  public :: calc_GF_Bias
  public :: calc_GF_Part
  public :: GF_Gamma_GF

contains

  
! Full converted GF.G.GF^\dagger routine for speed.
! This routine is extremely fast compared to any previous implementation.
! It relies on the fact that Gf only contains the left and the right electrode
! columns.
! Furthermore we retain all information by not imposing any symmetry in
! the product (TODO, check that we dont necessarily have this)
  subroutine GF_Gamma_GF(Offset,no_u_TS,no_LR,no_E, &
       GF,GammaT,GGG,nwork,work)

!  This routine returns GGG=GF.Gamma.GF^\dagger, where GF is a (no_u)x(no_L+no_R)
!  matrix and the states
!  corresponds to the (no_E) Left/Right electrode states (decided with Offset)
!  Gamma is a (no_E)x(no_E) matrix.

    use precision, only : dp

    implicit none

! *********************
! * INPUT variables   *
! *********************
    integer, intent(in) :: Offset  ! The offset for where Gamma lives
    integer, intent(in) :: no_u_TS ! no. states in contact region
    integer, intent(in) :: no_LR   ! no. states for both electrodes
    integer, intent(in) :: no_E    ! the size of the Gamma
    ! The Green's function
    complex(dp), intent(inout) :: GF(no_u_TS,no_LR)
    ! i (Sigma - Sigma^dagger)/2
    complex(dp), intent(inout) :: GammaT(no_E*no_E)
    ! A work array for doing the calculation... (nwork has to be larger than no_u_TS)
    integer,     intent(in)    :: nwork
    complex(dp), intent(inout) :: work(nwork)

! *********************
! * OUTPUT variables  *
! *********************
    complex(dp), intent(out) :: GGG(no_u_TS*no_u_TS)    !GF.GAMMA.GF^\dagger

! *********************
! * LOCAL variables   *
! *********************
    complex(dp), parameter :: z0  = dcmplx(0._dp, 0._dp)
    complex(dp), parameter :: z1  = dcmplx(1._dp, 0._dp)

    integer :: i, lE, NB, ind, iB

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE GFGammaGF' )
#endif

    call timer("GFGGF",1)

    lE = Offset + no_E - 1
    ! Number of times we can divide the large matrix
    NB = no_u_TS / no_E

    ! Loop over bottom row matrix 
    do iB = 0 , NB - 1
       
       ! Collect the top row of complex conjugated Gf
       ind = no_u_TS*no_E*iB+1
       do i = Offset , lE
          GGG(ind:ind-1+no_E) = dconjg(Gf(iB*no_E+1:(iB+1)*no_E,i))
          ind = ind + no_E
       end do
       
       ! Do Gamma.Gf^\dagger
       call zgemm('T','T',no_E,no_E,no_E,z1, &
            GammaT, no_E, &
            GGG(no_u_TS*no_E*iB+1),no_E, &
            z0, work,no_E)
       
       ! Calculate the Gf.Gamma.Gf^\dagger product for the entire column
       call zgemm('N','N',no_u_TS,no_E,no_E,z1, &
            Gf(1,Offset),no_u_TS, &
            work        ,   no_E, &
            z0, GGG(no_u_TS*no_E*iB+1),no_u_TS)
    
    end do

    ! in case the block size does not match the matrix order
    if ( NB * no_E /= no_u_TS ) then

       ! The size of the remaining block
       iB = no_u_TS - NB * no_E

       ! Copy over the block
       ind = no_u_TS*no_E*NB+1
       do i = Offset , lE
          ! So this is the complex conjugated of the iB'th block
          GGG(ind:ind-1+iB) = dconjg(Gf(NB*no_E+1:NB*no_E+iB,i))
          ind = ind + iB
       end do

       ! Do Gamma.Gf^\dagger
       call zgemm('T','T',no_E,iB,no_E,z1, &
            GammaT, no_E, &
            GGG(no_u_TS*no_E*NB+1),iB, &
            z0, work,no_E)
       
       ! Calculate the Gf.Gamma.Gf^\dagger product for the entire column
       call zgemm('N','N',no_u_TS,iB,no_E,z1, &
            Gf(1,Offset),no_u_TS, &
            work        ,   no_E, &
            z0, GGG(no_u_TS*no_E*NB+1),no_u_TS)

    end if

#ifdef TRANSIESTA_31
    ! Lets try and impose symmetry...
    call my_symmetrize(no_u_TS,GGG)
#endif

    call timer("GFGGF",2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS GFGammaGF' )
#endif

  end subroutine GF_Gamma_GF

#ifdef TRANSIESTA_31
subroutine my_symmetrize(N,M)
  use parallel, only : IONode
    integer    , intent(in) :: N
    complex(dp), intent(inout) :: M(N,N)
    integer :: i,j
    do j = 1 , N
       do i = 1 , j
!          if(ionode)print *,M(j,i),M(i,j)
          M(j,i) = dimag(M(i,j))
          M(i,j) = M(j,i)

       end do
    end do
  end subroutine my_symmetrize
#endif

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
  subroutine calc_GF(no_u_TS,GFinv,GF,ierr)
    
    use intrinsic_missing, only: EYE
    use precision, only: dp

    implicit none 

! *********************
! * INPUT variables   *
! *********************
    ! Sizes of the different regions...
    integer, intent(in) :: no_u_TS
    ! Work should already contain Z*S - H
    ! This may seem strange, however, it will clean up this routine extensively
    ! as we dont need to make two different routines for real and complex
    ! Hamiltonian values.
    complex(dp), intent(in out) :: GFinv(no_u_TS,no_u_TS) ! the inverted GF
    complex(dp), intent(out) :: GF(no_u_TS**2)
    integer,     intent(out) :: ierr              !inversion err

! Local variables
    integer :: ipvt(no_u_TS)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE getGF' )
#endif

    call timer('GFT',1) 

    ierr = 0

    call EYE(no_u_TS,GF)
    
! Invert directly
    call zgesv(no_u_TS,no_u_TS,GFinv,no_u_TS,ipvt,GF,no_u_TS,ierr)            
       
    call timer('GFT',2)  

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS getGF' )
#endif

  end subroutine calc_GF


! ##################################################################
! ## Calculating Greens functions in the reigon of the electrodes ## 
! ##                                                              ##          
! ##  Fully created by Nick Papior Andersen, nickpapior@gmail.com ##
! ##################################################################
  subroutine calc_GF_Bias(no_u_TS,no_L,no_R,GFinv,GF,ierr)
    
    use precision, only: dp

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
    complex(dp), intent(in out) :: GFinv(no_u_TS,no_u_TS) ! the inverted GF
    ! We only need Gf in the left and right blocks...
    complex(dp), intent(out) :: GF(no_u_TS*(no_L+no_R))
    integer,     intent(out) :: ierr              !inversion err

! Local variables
    integer :: ipvt(no_u_TS)
    integer :: i, o

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE getGF' )
#endif

    call timer('GFTB',1) 

    ierr = 0

    ! Create the RHS for inversion...
    GF(:) = dcmplx(0._dp,0._dp)

    ! Left identity
    do i = 0 , no_L - 1
       GF(i*no_u_TS+i+1) = dcmplx(1._dp,0._dp)
    end do
    ! Right identity
    o = no_L * no_u_TS + 1 + no_u_TS - no_R
    do i = 0 , no_R - 1
       GF(o+i*no_u_TS+i) = dcmplx(1._dp,0._dp)
    end do
    
! Invert directly
    call zgesv(no_u_TS,no_L+no_R,GFinv,no_u_TS,ipvt,GF,no_u_TS,ierr)            
       
    call timer('GFTB',2)  

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS getGF' )
#endif

  end subroutine calc_GF_Bias


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
  subroutine calc_GF_Part(no_u_TS,no_L,no_R, & ! Size of the problem
       GFinv,GF,ierr)
    
    use intrinsic_missing, only: EYE
    use precision, only: dp

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
    complex(dp), intent(in out) :: GFinv(no_u_TS,no_u_TS) ! the inverted GF
    complex(dp), intent(out) :: GF(no_u_TS*(no_u_TS-no_R-no_L))
    integer,     intent(out) :: ierr              ! inversion err

! Local variables
    integer :: ipvt(no_u_TS)
    integer :: i,j,ii

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE getGF' )
#endif

    call timer('GFT_P',1) 

    ierr = 0

! We already know that:
!   UseBulk == UpdateDMCR == .true.
    do j = 1, no_u_TS - no_L - no_R
       ii = (j-1) * no_u_TS
       do i = 1 , no_u_TS
          ii = ii + 1
          if      ( i        <= no_L ) then
             GF(ii) = dcmplx(0._dp,0._dp)
          else if ( i - no_L == j ) then
             GF(ii) = dcmplx(1._dp,0._dp)
          else
             GF(ii) = dcmplx(0._dp,0._dp)
          end if
       end do
    end do

! Invert directly
    call zgesv(no_u_TS,no_u_TS-no_L-no_R,GFinv,no_u_TS,ipvt,GF,no_u_TS,ierr)

    call timer('GFT_P',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS getGF' )
#endif

  end subroutine calc_GF_Part

end module m_ts_full_scat
