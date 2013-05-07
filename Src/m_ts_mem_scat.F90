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

module m_ts_mem_scat

  use precision, only : dp

  implicit none

  public :: calc_GF, calc_GF_Part
  public :: GF_Gamma_GF_Block
  public :: GF_Gamma_GF_Block_Short
  public :: UC_expansion_Sigma_Bulk
  public :: UC_expansion_Sigma
  public :: UC_expansion_Sigma_GammaT
  public :: weightDM, weightDMC
  public :: read_next_GS
  public :: get_scat_region

!  private

contains

  
! Full converted GF.G.GF^\dagger routine for speed.
! This routine is extremely fast compared to any previous implementation.
! The reason is re-use of critical segments.
! Furthermore we retain all information by not imposing any symmetry in
! the product (TODO, check that we dont necessarily have this)
  subroutine GF_Gamma_GF_Block(Offset,no_u_TS,no_E,GF,GammaT,GGG,nwork,work)

!  This routine returns GGG=GF.Gamma.GF^\dagger, where GF is a (no_u)x(no_u)
!  matrix and the states
!  corresponds to the (no_E) Left/Right electrode states (decided with Offset)
!  Gamma is a (no_E)x(no_E) matrix.

    use precision, only : dp

    implicit none

! *********************
! * INPUT variables   *
! *********************
    integer, intent(in) :: Offset ! The offset for where Gamma lives
    integer, intent(in) :: no_u_TS !no. states in contact region
    integer, intent(in) :: no_E ! the size of the Gamma
    ! The Green's function
    complex(dp), intent(inout) :: GF(no_u_TS,no_u_TS)
    ! i (Sigma - Sigma^dagger)/2
    complex(dp), intent(inout) :: GammaT(no_E*no_E)
    ! A work array for doing the calculation... (nwork has to be larger than no_u_TS)
    integer,     intent(in)    :: nwork
    complex(dp), intent(inout) :: work(nwork)

! *********************
! * OUTPUT variables  *
! *********************
    complex(dp), intent(out) :: GGG(no_u_TS*no_u_TS)    !GF.GAMMA.GF

! *********************
! * LOCAL variables   *
! *********************
    complex(dp), parameter :: z0  = dcmplx(0._dp,0._dp)
    complex(dp), parameter :: z1  = dcmplx(1._dp,0._dp)

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
       
       ! Collect the top row of Gf^\dagger
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
          ! So this is the transposed of iB'th block
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

  end subroutine GF_Gamma_GF_Block

! Full converted GF.G.GF^\dagger routine for speed.
! This routine is extremely fast compared to any previous implementation.
! It relies on the fact that Gf only contains the left and the right electrode
! columns.
! Furthermore we retain all information by not imposing any symmetry in
! the product (TODO, check that we dont necessarily have this)
  subroutine GF_Gamma_GF_Block_Short(Offset,no_u_TS,no_LR,no_E, &
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

  end subroutine GF_Gamma_GF_Block_Short

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

! ====================================================================
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
    o =  no_L * no_u_TS + 1
    do i = 0 , no_R - 1
       GF(o+i*no_u_TS+no_u_TS-no_R+i) = dcmplx(1._dp,0._dp)
    end do
    
! Invert directly
    call zgesv(no_u_TS,no_L+no_R,GFinv,no_u_TS,ipvt,GF,no_u_TS,ierr)            
       
    call timer('GFTB',2)  

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS getGF' )
#endif

! ====================================================================
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
    integer :: i,j,ii,o

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

! ====================================================================
  end subroutine calc_GF_Part


  subroutine UC_expansion_Sigma_Bulk(no_u,no_s,NRepA1,NRepA2, &
       na_u,lasto,nq,qb,wq,H,S,GS,Sigma,nwork,work)
    use intrinsic_missing, only : EYE
    use units, only : Pi
! ********************
! * INPUT variables  *
! ********************
    integer,  intent(in) :: no_u, no_s, NRepA1, NRepA2
    integer,  intent(in) :: na_u,lasto(0:na_u)
    integer,  intent(in) :: nq
    real(dp), intent(in) :: qb(3,nq), wq(nq)
    complex(dp), dimension(no_u,no_u,nq), intent(in) :: H, S
    complex(dp), dimension(no_u,no_u,nq), intent(inout) :: GS
! ********************
! * OUTPUT variables *
! ********************
    complex(dp), intent(out) :: Sigma(no_s,no_s)

    integer,  intent(in) :: nwork
    complex(dp), intent(inout) :: work(no_s,no_s)
! ********************
! * LOCAL variables  *
! ********************
    integer :: iq, ierr
    integer :: iow,iau,ia2,ia1,iuo
    integer :: jow,jau,ja2,ja1,juo
    integer :: ipvt(no_s)
    complex(dp), parameter :: zmPi2 = dcmplx(0._dp,-2._dp * Pi)
    complex(dp) :: ph

    call timer('ts_expand',1)
    call timer('ts_expandB',1)

    ! THis should never happen (work is TS-region!)
    if ( nwork < no_s**2 ) call die('Size of work-array is &
         &too small')

    call EYE(no_s,Sigma)

    if ( NRepA1 * NRepA2 == 1 ) then
       if ( no_u /= no_s ) call die('no_E/=no_s')

       ! We will read in a new GS for the following energy
       ! point, hence we do not need GS
       ! GS will be garbage from here on...
       ! We have the matrix to invert in the first no_s**2 values.
       call zgesv(no_s,no_s,GS(1,1,1),no_s,ipvt,Sigma,no_s,ierr)

    else

       ! This is the crucial calcuation.
       ! If we use bulk values in the electrodes
       ! we need not add the expanded H and S values to get the 
       ! electrode \Sigma. Hence, we need only expand
       ! surface Green's function
       iq = 1
       iow = 0
       do iau = 1 , na_u
        do ia2 = 0 , NRepA2-1
         do ia1 = 0 , NRepA1-1
          do iuo = 1 + lasto(iau-1) , lasto(iau)
           iow = iow + 1
           jow = 0
           do jau = 1 , na_u
            do ja2 = 0 , NRepA2-1
             do ja1 = 0 , NRepA1-1
              ph = wq(iq) * cdexp(zmPi2 * &
                   ( (ia1-ja1)*qb(1,iq) + (ia2-ja2)*qb(2,iq) ))
              do juo = 1 + lasto(jau-1) , lasto(jau) 
                 jow = jow + 1
                 
                 work(jow,iow) = ph * GS(juo,iuo,iq)
              end do !juo
             end do !ja1
            end do !ja2
           end do !jau
          end do !iuo
         end do !ia1
        end do !ia2
       end do !iau
       do iq = 2 , nq
        iow = 0
        do iau = 1 , na_u
         do ia2 = 0 , NRepA2-1
          do ia1 = 0 , NRepA1-1
           do iuo = 1 + lasto(iau-1) , lasto(iau)
            iow = iow + 1
            jow = 0
            do jau = 1 , na_u
             do ja2 = 0 , NRepA2-1
              do ja1 = 0 , NRepA1-1
               ph = wq(iq) * cdexp(zmPi2 * &
                    ( (ia1-ja1)*qb(1,iq) + (ia2-ja2)*qb(2,iq) ))

               do juo = 1 + lasto(jau-1) , lasto(jau) 
                  jow = jow + 1
                  
                  work(jow,iow) = work(jow,iow) + ph * GS(juo,iuo,iq)
               end do !juo
              end do !ja1
             end do !ja2
            end do !jau
           end do !iuo
          end do !ia1
         end do !ia2
        end do !iau
       end do !q-points
       
       ! We have the matrix to invert in the first no_s**2 values.
       call zgesv(no_s,no_s,work(1,1),no_s,ipvt,Sigma,no_s,ierr)

    end if
    if ( ierr /= 0 ) THEN
       write(*,*) 'Inversion of surface Greens function failed'
    end if

    call timer('ts_expandB',2)
    call timer('ts_expand',2)

  end subroutine UC_expansion_Sigma_Bulk


  subroutine UC_expansion_Sigma(ZEnergy,no_u,no_s,NRepA1,NRepA2, &
       na_u,lasto,nq,qb,wq,H,S,GS,Sigma,nwork,work)
    use intrinsic_missing, only : EYE
! ********************
! * INPUT variables  *
! ********************
    complex(dp), intent(in) :: ZEnergy
    integer,  intent(in) :: no_u, no_s, NRepA1, NRepA2
    integer,  intent(in) :: na_u,lasto(0:na_u)
    integer,  intent(in) :: nq
    real(dp), intent(in) :: qb(3,nq), wq(nq)
    complex(dp), dimension(no_u,no_u,nq), intent(in) :: H, S
    complex(dp), dimension(no_u,no_u,nq), intent(inout) :: GS
! ********************
! * OUTPUT variables *
! ********************
    complex(dp), intent(out) :: Sigma(no_s,no_s)

    integer,     intent(in)    :: nwork
    complex(dp), intent(inout) :: work(no_s,no_s,2)
! ********************
! * LOCAL variables  *
! ********************
    integer :: ierr
    integer :: io, jo
    integer :: ipvt(no_s)

    call timer('ts_expand',1)

    ! THis should never happen (work is TS-region!)
    if ( nwork < no_s**2*2 ) call die('Size of work-array is &
         &too small')
       
    call update_UC_expansion(ZEnergy,no_u,no_s,NRepA1,NRepA2, &
       na_u,lasto,nq,qb,wq,H,S,GS,nwork,work(1,1,1))

    call EYE(no_s,Sigma)

    if ( NRepA1 * NRepA2 == 1 ) then
       ! GS will be garbage from here on...
       ! We have the matrix to invert in the first no_s**2 values.
       call zgesv(no_s,no_s,GS(1,1,1),no_s,ipvt,Sigma,no_s,ierr)
    else
       ! We have the matrix to invert in the first no_s**2 values.
       call zgesv(no_s,no_s,work(1,1,1),no_s,ipvt,Sigma,no_s,ierr)
    end if
    if ( ierr /= 0 ) then
       write(*,*) 'Inversion of surface Greens function failed'
    end if

    ! Do:
    ! \Sigma = Z*S - H - \Sigma_bulk
    do jo = 1 , no_s
       do io = 1 , no_s
          Sigma(io,jo) = work(io,jo,2) - Sigma(io,jo)
       end do
    end do

    call timer('ts_expand',2)

  end subroutine UC_expansion_Sigma

  subroutine UC_expansion_Sigma_GammaT(UseBulk,ZEnergy,no_u,no_s,NRepA1,NRepA2, &
       na_u,lasto,nq,qb,wq,H,S,GS,Sigma,GammaT,nwork,work)
    use intrinsic_missing, only: EYE
! ********************
! * INPUT variables  *
! ********************
    complex(dp), intent(in) :: ZEnergy
    logical,  intent(in) :: UseBulk
    integer,  intent(in) :: no_u, no_s, NRepA1, NRepA2
    integer,  intent(in) :: na_u,lasto(0:na_u)
    integer,  intent(in) :: nq
    real(dp), intent(in) :: qb(3,nq), wq(nq)
    complex(dp), dimension(no_u,no_u,nq), intent(in) :: H, S
    complex(dp), dimension(no_u,no_u,nq), intent(inout) :: GS
! ********************
! * OUTPUT variables *
! ********************
    complex(dp), intent(out) :: Sigma(no_s,no_s)
!    real(dp), intent(out)    :: Gamma(no_s,no_s)
    complex(dp), intent(out)    :: GammaT(no_s,no_s)

    integer,  intent(in) :: nwork
    complex(dp), intent(inout) :: work(no_s,no_s,2)
! ********************
! * LOCAL variables  *
! ********************
    complex(dp), parameter :: zihalf = dcmplx(0._dp,0.5_dp)
    integer :: ierr
    integer :: io,jo
    integer :: ipvt(no_s)

    call timer('ts_expand',1)
    call timer('ts_expandG',1)
    
    ! THis should never happen (work is TS-region!)
    if ( nwork < no_s**2*2 ) call die('Size of work-array is &
         &too small')

    call update_UC_expansion(ZEnergy,no_u,no_s,NRepA1,NRepA2, &
       na_u,lasto,nq,qb,wq,H,S,GS,nwork,work(1,1,1))

    call EYE(no_s,Sigma)

    if ( NRepA1 * NRepA2 == 1 ) then
       ! We have the matrix to invert in the first no_s**2 values.
       ! NOTICE THAT GS will be rubbish from here on!
       call zgesv(no_s,no_s,GS(1,1,1),no_s,ipvt,Sigma,no_s,ierr)
    else
       ! We have the matrix to invert in the first no_s**2 values.
       call zgesv(no_s,no_s,work(1,1,1),no_s,ipvt,Sigma,no_s,ierr)
    end if
    if ( ierr /= 0 ) then
       write(*,*) 'Inversion of surface Greens function failed'
    end if

    if ( UseBulk ) then

       ! Do:
       ! \Sigma = Z*S - H - \Sigma_bulk
       do jo = 1 , no_s
          do io = 1 , no_s
             work(io,jo,1) = work(io,jo,2) - Sigma(io,jo)
          end do
       end do

       ! Do
       ! \Gamma = i ( \Sigma - \Sigma^\dagger)/2
       do jo = 1 , no_s
          do io = 1 , jo - 1
             GammaT(jo,io) = zihalf * ( &
                  work(io,jo,1)-dconjg(work(jo,io,1)) )
             GammaT(io,jo) = zihalf * ( &
                  work(jo,io,1)-dconjg(work(io,jo,1)) )
          end do
          GammaT(jo,jo) = zihalf * ( &
               work(jo,jo,1)-dconjg(work(jo,jo,1)) )
       end do

    else
       
       ! Do:
       ! \Sigma = Z*S - H - \Sigma_bulk
       do jo = 1 , no_s
          do io = 1 , no_s
             Sigma(io,jo) = work(io,jo,2) - Sigma(io,jo)
          end do
       end do

       ! Do
       ! \Gamma = i ( \Sigma - \Sigma^\dagger)/2
       do jo = 1 , no_s
          do io = 1 , jo - 1
             GammaT(jo,io) = zihalf * ( &
                  Sigma(io,jo)-dconjg(Sigma(jo,io)) )
             GammaT(io,jo) = zihalf * ( &
                  Sigma(jo,io)-dconjg(Sigma(io,jo)) )
          end do
          GammaT(jo,jo) = zihalf * ( &
               Sigma(jo,jo)-dconjg(Sigma(jo,jo)) )
       end do

    end if

    call timer('ts_expandG',2)
    call timer('ts_expand',2)
       
  end subroutine UC_expansion_Sigma_GammaT

  subroutine update_UC_expansion(ZEnergy,no_u,no_s,NRepA1,NRepA2, &
       na_u,lasto,nq,qb,wq,H,S,GS,nwork,work)
    use units, only : Pi
! ********************
! * INPUT variables  *
! ********************
    complex(dp), intent(in) :: ZEnergy
    integer,  intent(in) :: no_u, no_s, NRepA1, NRepA2
    integer,  intent(in) :: na_u,lasto(0:na_u)
    integer,  intent(in) :: nq
    real(dp), intent(in) :: qb(3,nq), wq(nq)
    complex(dp), dimension(no_u,no_u,nq), intent(in) :: H, S, GS
! ********************
! * OUTPUT variables *
! ********************
    integer,     intent(in)    :: nwork
    complex(dp), intent(inout) :: work(no_s,no_s,2)
! ********************
! * LOCAL variables  *
! ********************
    integer :: iq
    integer :: iow,iau,ia2,ia1,iuo
    integer :: jow,jau,ja2,ja1,juo
    complex(dp), parameter :: zmPi2 = dcmplx(0._dp,-2.0_dp * Pi)
    complex(dp) :: ph

    if ( NRepA1 * NRepA2 == 1 ) then
       if ( no_u /= no_s ) call die('no_E/=no_s')

       ! We do not need this...
       !work(:,:,1) = GS(:,:,1)
       work(:,:,2) = ZEnergy * S(:,:,1) - H(:,:,1)

    else

       ! This is the crucial calcuation.
       ! If we use bulk values in the electrodes
       ! we need not add the expanded H and S values to get the 
       ! electrode \Sigma. Hence, we need only expand
       ! surface Green's function
       iq = 1
       iow = 0
       do iau = 1 , na_u
        do ia2 = 0 , NRepA2-1
         do ia1 = 0 , NRepA1-1
          do iuo = 1 + lasto(iau-1) , lasto(iau)
           iow = iow + 1
           jow = 0
           do jau = 1 , na_u
            do ja2 = 0 , NRepA2-1
             do ja1 = 0 , NRepA1-1
              ! TODO qb(:,1) == 0.0_dp in any case currently used
              ! So we could do without.
              ph = wq(iq) * cdexp(zmPi2 * &
                   ( (ia1-ja1)*qb(1,iq) + (ia2-ja2)*qb(2,iq) ))

              do juo = 1 + lasto(jau-1) , lasto(jau)
                jow = jow + 1
                
                work(jow,iow,1) = ph * GS(juo,iuo,iq)
                
                work(jow,iow,2) = ph * (ZEnergy*S(juo,iuo,iq)-H(juo,iuo,iq))
                
              end do !juo
             end do !ja1
            end do !ja2
           end do !jau
          end do !iuo
         end do !ia1
        end do !ia2
       end do !iau
       do iq = 2 , nq
        iow = 0
        do iau = 1 , na_u
         do ia2 = 0 , NRepA2-1
          do ia1 = 0 , NRepA1-1
           do iuo = 1 + lasto(iau-1) , lasto(iau)
            iow = iow + 1
            jow = 0
            do jau = 1 , na_u
             do ja2 = 0 , NRepA2-1
              do ja1 = 0 , NRepA1-1
               ph = wq(iq) * cdexp(zmPi2 * &
                    ( (ia1-ja1)*qb(1,iq) + (ia2-ja2)*qb(2,iq) ))

               do juo = 1 + lasto(jau-1) , lasto(jau)
                  jow = jow + 1
                  
                  work(jow,iow,1) = work(jow,iow,1) + ph * GS(juo,iuo,iq)
   
                  work(jow,iow,2) = work(jow,iow,2) + ph * &
                       (ZEnergy*S(juo,iuo,iq)-H(juo,iuo,iq))
   
               end do !juo
              end do !ja1
             end do !ja2
            end do !jau
           end do !iuo
          end do !ia1
         end do !ia2
        end do !iau
       end do !q-points

    end if

  end subroutine update_UC_expansion



! ##################################################################
! ## Subroutine which read-in GS for the 1x1 surface cell         ##
! ##    The expansion of the arrays are performed else-where.     ##
! ##                                                              ##
! ##    Nick Papior Andersen, nickpapior@gmail.com                ##
! ## Fully recoded to conform with the memory reduced TranSIESTA  ##
! ##################################################################

  subroutine read_next_GS(uGF,NEReqs,ikpt,no_GS,nq,HAA,SAA,GAA,ZEnergy, &
       nwork,work)

    use precision, only : dp
    use parallel,  only : IONode, Node, Nodes

#ifdef MPI
    use mpi_siesta, only : MPI_Bcast, MPI_Isend, MPI_Irecv,  &
         MPI_Sum, MPI_Integer, MPI_Status_Size, MPI_Comm_World
    use mpi_siesta, only : MPI_double_complex

#endif
! *********************
! * INPUT variables   *
! *********************
    ! file-unit, and k-point index
    integer, intent(in) :: uGF, ikpt
    integer, intent(in) :: NEReqs
    ! Size of the arrays we will read etc.
    integer, intent(in) :: no_GS, nq
    ! Hamiltonian, overlap, and GS
    complex(dp), intent(inout) :: HAA(no_GS,no_GS,nq)
    complex(dp), intent(inout) :: SAA(no_GS,no_GS,nq)
    complex(dp), intent(inout) :: GAA(no_GS,no_GS,nq)
    complex(dp), intent(in) :: ZEnergy
    ! The work array passed, this means we do not have
    ! to allocate anything down here.
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(no_GS,no_GS,nq)

! *********************
! * LOCAL variables   *
! *********************
    real(dp), parameter :: EPS = 1.d-7
    integer :: read_Size

#ifdef MPI
    integer :: MPIerror, Request, Status(MPI_Status_Size)
#endif

    complex(dp) :: ZE_cur, wGFi
    integer :: iEni, iNode, ikGS

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE read_next_GS' )
#endif

#ifdef MPI
    ! This will make the timings be more realistic
    call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif

    call timer('rn_GS',1)

    read_Size = no_GS * no_GS * nq

    ! Check if the number of energy points requested are 
    ! inconsistent
    if ( NEReqs <= 0 ) then
       if(IONode) &
            write(*,'(a,i0,a)') 'ERROR GetSFE: Requested E-points=', &
            NEReqs,'< 0'  
       call die('ERROR in reading GF file')
    else if ( Nodes < NEReqs ) then
       if(IONode) &
            write(*,'(2(a,i0))')  &
            'ERROR GetSFE: Requested E-points= ', &
            NEReqs,' > Nodes = ', Nodes
       call die('ERROR in reading GF file')
    end if

    if ( nwork < read_Size ) then
       write(*,*) 'Size of work array while reading GS was not large &
            &enough. Something went wrong.'
       call die('ERROR in reading GF file')
    end if
    
! Loop over nodes. If root send stuff out to relevant nodes, if
! not root receive information from root.

    ! We only loop over the requested energy points...
    do iNode = 0, NEReqs - 1

       if ( IONode ) then
          ! read in header of GF-point
          read(uGF) iEni,ZE_cur,wGFI,ikGS
       endif

#ifdef MPI
       ! distribute the energy index (in case of iEni == 1)
       call MPI_Bcast(iEni,1,MPI_integer,0, &
            MPI_Comm_World,MPIerror)
       ! distribute the current energy point
       if ( IONode .and. Node == iNode ) then
          ! do nothing
       else if ( Node == iNode ) then
          ! recieve from the host node
          call MPI_iRecv(ZE_cur,1,MPI_Double_Complex,    0,iNode, &
               MPI_Comm_World,Request,MPIerror)
          call MPI_Wait(Request,Status,MPIerror)
       else if ( IONode ) then
          call MPI_iSend(ZE_cur,1,MPI_Double_Complex,iNode,iNode, &
               MPI_Comm_World,Request,MPIerror)
       endif
#endif


       ! The test of the energy-point is performed on
       ! the calculating node...
       if ( Node == iNode ) then
          if ( cdabs(ZEnergy-ZE_cur) > EPS ) then
             call die('Energy point in GF file does &
                  not match the internal energy-point in transiesta. &
                  &Please correct your GF files.')
          end if
       end if
       
       ! If the k-point does not match what we expected...
       if ( IONode .and. ikpt /= ikGS ) then
          call die('Read k-point in GF file does not match &
               &the requested k-point. Please correct your &
               &GF files.')
       end if

       if ( iEni == 1 ) then

          ! We have to read in electrode Hamiltonian and overlap...
          if ( IONode ) then
             read(uGF) HAA
             read(uGF) SAA
          endif

#ifdef MPI
          call MPI_Bcast(HAA(1,1,1),read_Size,MPI_Double_Complex, &
               0,MPI_Comm_World,MPIerror)
          call MPI_Bcast(SAA(1,1,1),read_Size,MPI_Double_Complex, &
               0,MPI_Comm_World,MPIerror)
#endif
       end if 


       
#ifdef MPI
       if ( IONode .and. iNode == Node ) then
#endif
          ! Read in point directly to correct array...
          read(uGF) GAA
#ifdef MPI
       else if ( IONode ) then
          ! We need to ensure that the energy point has been sent
          ! I.e. this is for the ZE_cur message
          call MPI_Wait(Request,Status,MPIerror)

          read(uGF) work(1:no_GS,1:no_GS,1:nq)
          
          call MPI_ISend(work(1,1,1),read_Size,MPI_Double_Complex, &
               iNode,1,MPI_Comm_World,Request,MPIerror) 
          call MPI_Wait(Request,Status,MPIerror)
       else if ( Node ==  iNode ) then
          call MPI_IRecv(GAA(1,1,1),read_Size,MPI_Double_Complex, &
               0    ,1,MPI_Comm_World,Request,MPIerror)
          call MPI_Wait(Request,Status,MPIerror)
       end if
#endif

    end do

    call timer('rn_GS',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS read_next_GS' )
#endif

  end subroutine read_next_GS

! ##################################################################
! ##   Mixing the Density matrixes according to the smallest      ##
! ##    realspace integral                                        ##
! ##                                                              ##
! ##  Version 011200  Kurt Stokbro, ks@mic.dtu.dk                 ##
! ##  Heavily edited by Nick Papior Andersen to be used with the  ##
! ##  full sparsity pattern in transiesta                         ##
! ##################################################################
  subroutine weightDM(no_C_L, no_C_R, &
       SpArrDML , SpArrDMR , SpArrDMneqL, SpArrDMneqR, &
       SpArrEDML, SpArrEDMR)
!  This routine find weight for the DM integral. On output
!  DML := w (DML+DMneqR) + (1-w) (DMR+DMneqL)
!  EDML:= w EDML +(1-w) EDMR
!  In left electrode w=1 and in right electrode w=0

    use precision, only: dp
    use parallel,  only: IONode
    use class_Sparsity
    use class_dSpData1D

    implicit none

! *********************
! * INPUT variables   *
! *********************
    integer, intent(in) :: no_C_L, no_C_R
! *********************
! * OUTPUT variables  *
! *********************
    ! Contour part of DM integration
    type(dSpData1D), intent(inout) :: SpArrDML, SpArrDMR
    ! Real-axis part of DM integration
    type(dSpData1D), intent(inout) :: SpArrDMneqL, SpArrDMneqR
    ! L-R estimates of EDM
    type(dSpData1D), intent(inout) :: SpArrEDML, SpArrEDMR

! *********************
! * LOCAL variables   *
! *********************
    real(dp) :: wL,wR,wSUM

    ! arrays for looping in the sparsity pattern
    type(Sparsity), pointer :: sp
    real(dp), pointer :: DML(:), DMR(:)
    real(dp), pointer :: DMneqL(:), DMneqR(:)
    real(dp), pointer :: EDML(:), EDMR(:)
    integer, pointer :: l_ncol(:)
    integer, pointer :: l_ptr(:)
    integer, pointer :: l_col(:)
    integer :: nr
    integer :: io, jo, ind, j
    ! For error estimation
    integer :: eM_i,eM_j,neM_i,neM_j
    real(dp) :: eM, neM, tmp

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE weightDM' )
#endif

    ! TODO Enforce that sparsity is the same
    ! (however, we know that they are the same.
    sp => spar(SpArrDML)
    l_ncol => n_col   (sp)
    l_ptr  => list_ptr(sp)
    l_col  => list_col(sp)
    nr = nrows(sp)
    ! Obtain the values in the arrays...
    DML => val(SpArrDML)
    DMR => val(SpArrDMR)
    DMneqL => val(SpArrDMneqL)
    DMneqR => val(SpArrDMneqR)
    EDML => val(SpArrEDML)
    EDMR => val(SpArrEDMR)

    ! initialize the errors
    eM  = 0._dp
    neM = 0._dp

    do io = 1 , nr
       ! We are in a buffer region...
       if ( l_ncol(io) == 0 ) cycle
       do j = 1 , l_ncol(io)

          ind = l_ptr(io)+j
          ! Retrive the connecting orbital
          jo = l_col(ind)

          wL = DMneqL(ind)*DMneqL(ind)
          wR = DMneqR(ind)*DMneqR(ind)
          wSUM = wL + wR

          ! The weights
          if ( wSUM > 0._dp ) then
             wL = wL / wSUM
             wR = wR / wSUM
          else
             wL = 0.5_dp
             wR = 0.5_dp
             wSUM = 1._dp
          end if
          
          ! Do error estimation (capture before update)
          tmp = (DML(ind) + DMneqR(ind) - DMR(ind) - DMneqL(ind))**2

          DML(ind)  = wL * (DML(ind) + DMneqR(ind)) &
                    + wR * (DMR(ind) + DMneqL(ind))
          EDML(ind) = wL * EDML(ind) + wR * EDMR(ind)

          ! this is absolute error
          if ( tmp > eM ) then
             eM   = tmp
             eM_i = io
             eM_j = jo
          end if
          ! this is normalized absolute error
          tmp = tmp * wL * wR
          if ( tmp > neM ) then
             neM   = tmp
             neM_i = io
             neM_j = jo
          end if

       end do
    end do

    call print_error_estimate(IONode,eM,eM_i,eM_j,neM,neM_i,neM_j)    

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS weightDM' )
#endif

  end subroutine weightDM


! ##################################################################
! ##   Mixing the Density matrixes according to the smallest      ##
! ##    realspace integral *** COMPLEX VERSION ***                ##
! ##                                                              ##
! ##  Version 011200  Kurt Stokbro, ks@mic.dtu.dk                 ##
! ##  Heavily edited by Nick Papior Andersen to be used with the  ##
! ##  full sparsity pattern in transiesta                         ##
! ##################################################################
  subroutine weightDMC(no_C_L, no_C_R, &
       SpArrDML , SpArrDMR , SpArrDMneqL, SpArrDMneqR, &
       SpArrEDML, SpArrEDMR)
!  This routine find weight for the DM integral. On output
!  DML := w (DML+DMneqR) + (1-w) (DMR+DMneqL)
!  EDML:= w EDML +(1-w) EDMR
!  In left electrode w=1 and in right electrode w=0

    use precision, only: dp
    use parallel,  only: IONode
    use class_Sparsity
    use class_zSpData1D

    implicit none

! *********************
! * INPUT variables   *
! *********************
    integer, intent(in) :: no_C_L, no_C_R
! *********************
! * OUTPUT variables  *
! *********************
    ! Contour part of DM integration
    type(zSpData1D), intent(inout) :: SpArrDML, SpArrDMR
    ! Real-axis part of DM integration
    type(zSpData1D), intent(inout) :: SpArrDMneqL, SpArrDMneqR
    ! L-R estimates of EDM
    type(zSpData1D), intent(inout) :: SpArrEDML, SpArrEDMR

! *********************
! * LOCAL variables   *
! *********************
    real(dp) :: wL,wR,wSUM

    ! arrays for looping in the sparsity pattern
    type(Sparsity), pointer :: sp
    complex(dp), pointer :: DML(:), DMR(:)
    complex(dp), pointer :: DMneqL(:), DMneqR(:)
    complex(dp), pointer :: EDML(:), EDMR(:)
    integer, pointer :: l_ncol(:)
    integer, pointer :: l_ptr(:)
    integer, pointer :: l_col(:)
    integer :: nr
    integer :: io, jo, ind, j
    ! For error estimation
    integer :: eM_i,eM_j,neM_i,neM_j
    complex(dp) :: ztmp
    real(dp) :: eM, neM, rtmp

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE weightDMC' )
#endif

    ! TODO Enforce that sparsity is the same
    ! (however, we know that they are the same.
    sp => spar(SpArrDML)
    l_ncol => n_col   (sp)
    l_ptr  => list_ptr(sp)
    l_col  => list_col(sp)
    nr = nrows(sp)
    ! Obtain the values in the arrays...
    DML => val(SpArrDML)
    DMR => val(SpArrDMR)
    DMneqL => val(SpArrDMneqL)
    DMneqR => val(SpArrDMneqR)
    EDML => val(SpArrEDML)
    EDMR => val(SpArrEDMR)

    ! initialize the errors
    eM  = 0._dp
    neM = 0._dp

    do io = 1 , nr
       ! We are in a buffer region...
       if ( l_ncol(io) == 0 ) cycle
       do j = 1 , l_ncol(io)

          ind = l_ptr(io)+j
          ! Retrive the connecting orbital
          jo = l_col(ind)

          ! It is weighted in the density (not the imaginary part of 
          ! \rho_L). Note that here \rho_L\equiv -i\rho_L!
          wL = aimag(DMneqL(ind)) ** 2
          wR = aimag(DMneqR(ind)) ** 2
          wSUM = wL + wR
 
          ! The weights (in any case we always have the full Gf.G.Gf, hence
          ! it is safe to use this method always.
          ! No need to force either correction term in the left/right regions
          if ( wSUM > 0._dp ) then
             wL = wL / wSUM
             wR = wR / wSUM
          else
             wL = 0.5_dp
             wR = 0.5_dp
             wSUM = 1._dp
          end if

          ! We need to capture the error before the update...
          ztmp = DML(ind) + DMneqR(ind) - DMR(ind) - DMneqL(ind)

          DML(ind)  = wL * (DML(ind) + DMneqR(ind)) &
                    + wR * (DMR(ind) + DMneqL(ind))
          EDML(ind) = wL * EDML(ind) + wR * EDMR(ind)

          ! Do error estimation (we are only interested in 
          ! error in the density...)
          
          rtmp = aimag(ztmp) * aimag(ztmp)
          ! this is absolute error
          if ( rtmp > eM ) then
             eM = rtmp
             eM_i = io
             eM_j = jo
          end if
          ! this is normalized absolute error
          rtmp = rtmp * wL * wR
          if ( rtmp > neM ) then
             neM = rtmp
             neM_i = io
             neM_j = jo
          end if

       end do
    end do

    call print_error_estimate(IONode,eM,eM_i,eM_j,neM,neM_i,neM_j)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS weightDMC' )
#endif

  end subroutine weightDMC

  subroutine print_error_estimate(IONode,eM,eM_i,eM_j,neM,neM_i,neM_j)
    logical, intent(in) :: IONode
    real(dp), intent(in) :: eM, neM
    integer, intent(in) :: eM_i,eM_j, neM_i,neM_j

    if ( IONode ) then
       write(*,'(a,'' |('',i5,'','',i5,'')| = '',g9.4,&
            &'' , |('',i5,'','',i5,'')|~ = '',g9.4)') &
            'ts: integration EE.:',&
            eM_i,eM_j,sqrt(eM), &
            neM_i,neM_j,sqrt(neM)
    end if

  end subroutine print_error_estimate
    

!     Function for calculating which region in the scattering matrix the
!     hamiltonian elements adhere.
  pure function get_scat_region(io,noL,jo, noR,no_u)
    integer, intent(in) :: io, noL, jo, noR, no_u
    integer :: get_scat_region

    if ( io < 1 .and. jo < 1 ) then
       get_scat_region = 1 ! Left buffer
    else if ( (io < 1 .and. jo <= noL) .or. &
         (io <= noL .and. jo < 1) ) then
       get_scat_region = 2 ! Left buffer / left electrode
    else if ( (1 <= io .and. 1 <= jo) .and. &
         (io <= noL .and. jo <= noL) ) then
       get_scat_region = 3 ! Left electrode
    else if ( ((1 <= io .and. noL < jo) .and. &
         (io <= noL .and. jo <= no_u-noR)) .or. &
         ((noL < io .and. 1 <= jo) .and. &
         (io <= no_u-noR .and. jo <= noL)) ) then
       get_scat_region = 4 ! Left electrode / contact region
    else if ( (noL < io .and. noL < jo) .and. &
         (io <= no_u-noR .and. jo <= no_u-noR) ) then
       get_scat_region = 5 ! Contact region
    else if ( ((noL < io .and. no_u-noR < jo) .and. &
         (io <= no_u-noR .and. jo <= no_u)) .or. &
         ((no_u-noR < io .and. noL < jo) .and. &
         (io <= no_u .and. jo <= no_u-noR)) ) then
       get_scat_region = 6 ! Contact region / right electrode
    else if ( (no_u-noR < io .and. no_u-noR < jo) .and. &
         (io <= no_u .and. jo <= no_u) ) then
       get_scat_region = 7 ! Right electrode
    else if ( ((no_u-noR < io .and. no_u < jo) .and. &
         io <= no_u) .or. &
         ((no_u < io .and. no_u-noR < jo) .and. &
         jo <= no_u) ) then
       get_scat_region = 8 ! Right electrode / right buffer
    else if ( no_u < io .and. no_u < jo ) then
       get_scat_region = 9 ! Right buffer
    else
       get_scat_region = 0 ! Everything else (LE/RE ... etc.)
    end if

  end function get_scat_region

end module m_ts_mem_scat
