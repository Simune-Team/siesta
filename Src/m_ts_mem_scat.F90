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
  public :: GF_Gamma_GF
  public :: UC_expansion_Sigma_Bulk
  public :: UC_expansion_Sigma
  public :: UC_expansion_Sigma_Gamma
  public :: weightDM, weightDMC
  public :: read_next_GS
  public :: get_scat_region

  private

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
  subroutine calc_GF(UseBulk, &
       no_u_TS,no_L,no_R, & ! Size of the problem
       SigmaL,SigmaR, & ! Electrode self-energies
       GFinv,GF,ierr)
    
    use intrinsic_missing, only: EYE
    use precision, only: dp

    implicit none 

! *********************
! * INPUT variables   *
! *********************
    logical, intent(in) :: UseBulk           ! if true self-energy only is input else
!                                 z*S-H-Sigma for bulk is in sfe
    ! Sizes of the different regions...
    integer, intent(in) :: no_u_TS, no_L, no_R
    complex(dp) :: SigmaL(no_L,no_L)   ! Left electrode GF
    complex(dp) :: SigmaR(no_R,no_R)   ! Right electrode GF
    ! Work should already contain Z*S - H
    ! This may seem strange, however, it will clean up this routine extensively
    ! as we dont need to make two different routines for real and complex
    ! Hamiltonian values.
    complex(dp), intent(in out) :: GFinv(no_u_TS**2) ! the inverted GF
    complex(dp), intent(out) :: GF(no_u_TS**2)
    integer,     intent(out) :: ierr              !inversion err

! Local variables
    integer :: ipvt(no_u_TS)
    integer :: i,j,ii,o

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE getGF' )
#endif

    call timer('GFT',1) 

    ierr = 0

! We need to check the input to the routine.
! We can calculate either the full GF, if .not. (UseBulk .and

    call EYE(no_u_TS,GF)
    
! Adjust the left electrode part
    if ( UseBulk ) then 
       ! If we use the bulk electrodes, we will overwrite
       ! the "regular" Hamiltonian elements by the 
       ! Self-energy of the electrode (thereby assuming 
       ! that the electrode states are the perfect bulk system)
       do j = 1, no_L
          ii = (j-1) * no_u_TS
          do i = 1, no_L
             ii = ii + 1
             GFinv(ii) = SigmaL(i,j)
          end do              !i
       end do                 !j
    else
       ! We don't consider the electrode to be the bulk states,
       ! hence, we will adjust the unit-cell Hamiltonian by the 
       ! self-energy terms.
       do j = 1, no_L
          ii = (j-1) * no_u_TS
          do i = 1, no_L
             ii = ii + 1
             GFinv(ii) = GFinv(ii) - SigmaL(i,j)
          end do              !i
       end do                 !j
    end if                    !USEBULK

! Adjust the right electrode part
    o = no_u_TS - no_R
    if ( UseBulk ) then
       do j = 1 , no_R
          ii = (o + j - 1) * no_u_TS + o
          do i = 1 , no_R
             ii = ii + 1
             GFinv(ii) = SigmaR(i,j)
          end do              !i
       end do                 !j
    else
       do j = 1 , no_R
          ii = (o + j - 1) * no_u_TS + o
          do i = 1 , no_R
             ii = ii + 1
             GFinv(ii) = GFinv(ii) - SigmaR(i,j)
          end do              !i
       end do                 !j
    end if

! Invert directly
    call zgesv(no_u_TS,no_u_TS,GFinv,no_u_TS,ipvt,GF,no_u_TS,ierr)            
       
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
       GFinv,GF,ierr)
    
    use intrinsic_missing, only: EYE
    use precision, only: dp

    implicit none 

! *********************
! * INPUT variables   *
! *********************
    ! Sizes of the different regions...
    integer, intent(in) :: no_u_TS, no_L, no_R
    complex(dp) :: SigmaL(no_L,no_L)   ! Left electrode GF
    complex(dp) :: SigmaR(no_R,no_R)   ! Right electrode GF
    ! Work should already contain Z*S - H
    ! This may seem strange, however, it will clean up this routine extensively
    ! as we dont need to make two different routines for real and complex
    ! Hamiltonian values.
    complex(dp), intent(in out) :: GFinv(no_u_TS**2) ! the inverted GF
    complex(dp), intent(out) :: GF(no_u_TS*(no_u_TS-no_R-no_L))
    integer,     intent(out) :: ierr              ! inversion err

! Local variables
    integer :: ipvt(no_u_TS)
    integer :: i,j,ii,jj,o

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

    ! If we use the bulk electrodes, we will overwrite
    ! the "regular" Hamiltonian elements by the 
    ! Self-energy of the electrode (thereby assuming 
    ! that the electrode states are the perfect bulk system)
    do j = 1, no_L
       ii = (j-1) * no_u_TS
       do i = 1, no_L
          ii = ii + 1
          GFinv(ii) = SigmaL(i,j)
       end do              !i
    end do                 !j

! Adjust the right electrode part
    o = no_u_TS - no_R
    do j = 1 , no_R
       ii = (o + j - 1) * no_u_TS + o
       do i = 1 , no_R
          ii = ii + 1
          GFinv(ii) = SigmaR(i,j)
       end do              !i
    end do                 !j


! Invert directly
    call zgesv(no_u_TS,no_u_TS-no_L-no_R,GFinv,no_u_TS,ipvt,GF,no_u_TS,ierr)

    call timer('GFT_P',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS getGF' )
#endif

! ====================================================================
  END subroutine calc_GF_Part


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
    complex(dp), dimension(no_u,no_u,nq), intent(in) :: H, S, GS
! ********************
! * OUTPUT variables *
! ********************
    complex(dp), intent(out) :: Sigma(no_s,no_s)

    integer,  intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)
! ********************
! * LOCAL variables  *
! ********************
    integer :: iq, ierr
    integer :: iow,iau,ia2,ia1,iuo
    integer ::     jau,ja2,ja1,juo
    integer :: ipvt(no_s)
    real(dp), parameter :: Pi2 = 2._dp * Pi
    complex(dp) :: ph

    call timer('ts_expand',1)
    call timer('ts_expandB',1)

    ! THis should never happen (work is TS-region!)
    if ( nwork < no_s**2*2 ) call die('Size of work-array is &
         &too small')

    do iow = 1 , no_s*no_s
       work(iow) = dcmplx(0._dp,0._dp)
    end do

    ! This is the crucial calcuation.
    ! If we use bulk values in the electrodes
    ! we need not add the expanded H and S values to get the 
    ! electrode \Sigma. Hence, we need only expand
    ! surface Green's function
    do iq = 1 , nq
     iow = 0
     do iau = 1 , na_u
      do ia2 = 0 , NRepA2-1
       do ia1 = 0 , NRepA1-1
        do iuo = 1 + lasto(iau-1) , lasto(iau)
         do jau = 1 , na_u
          do ja2 = 0 , NRepA2-1
           do ja1 = 0 , NRepA1-1
            do juo = 1 + lasto(jau-1) , lasto(jau) 
               iow = iow + 1
               ph = wq(iq) * cdexp(dcmplx(0.0_dp,-Pi2) * &
                    ( (ia1-ja1)*qb(1,iq) + (ia2-ja2)*qb(2,iq) ))
               
               work(iow) = work(iow) + ph * GS(juo,iuo,iq)
            end do !juo
           end do !ja1
          end do !ja2
         end do !jau
        end do !iuo
       end do !ia1
      end do !ia2
     end do !iau
    end do !q-points

    call EYE(no_s,Sigma)

    ! We have the matrix to invert in the first no_s**2 values.
    call zgesv(no_s,no_s,work(1),no_s,ipvt,Sigma,no_s,ierr)

    if ( ierr /= 0 ) THEN
       write(*,*) 'Inversion of surface Greens function failed'
    end if

    call timer('ts_expandB',2)
    call timer('ts_expand',2)

  end subroutine UC_expansion_Sigma_Bulk


  subroutine UC_expansion_Sigma(ZEnergy,no_u,no_s,NRepA1,NRepA2, &
       na_u,lasto,nq,qb,wq,H,S,GS,Sigma,nwork,work)
    use units, only : Pi
    use intrinsic_missing, only : EYE
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
    complex(dp), intent(out) :: Sigma(no_s,no_s)

    integer,     intent(in)    :: nwork
    complex(dp), intent(inout) :: work(nwork)
! ********************
! * LOCAL variables  *
! ********************
    integer :: iq, ierr
    integer :: iow1,iow2,iau,ia2,ia1,iuo
    integer ::           jau,ja2,ja1,juo
    integer :: ipvt(no_s)
    real(dp), parameter :: Pi2 = 2.0_dp * Pi
    complex(dp) :: ph

    call timer('ts_expand',1)

    ! THis should never happen (work is TS-region!)
    if ( nwork < no_s**2*2 ) call die('Size of work-array is &
         &too small')

    do iow1 = 1 , no_s*no_s*2
       work(iow1) = dcmplx(0._dp,0._dp)
    end do

    ! This is the crucial calcuation.
    ! If we use bulk values in the electrodes
    ! we need not add the expanded H and S values to get the 
    ! electrode \Sigma. Hence, we need only expand
    ! surface Green's function
    do iq = 1 , nq
     iow1 = 0
     iow2 = no_s**2
     do iau = 1 , na_u
      do ia2 = 0 , NRepA2-1
       do ia1 = 0 , NRepA1-1
        do iuo = 1 + lasto(iau-1) , lasto(iau) 
         do jau = 1 , na_u
          do ja2 = 0 , NRepA2-1
           do ja1 = 0 , NRepA1-1
            do juo = 1 + lasto(jau-1) , lasto(jau) 
               iow1 = iow1 + 1
               iow2 = iow2 + 1
               ph = wq(iq)*cdexp(dcmplx(0.0_dp,-Pi2) * &
                    ( (ia1-ja1)*qb(1,iq) + (ia2-ja2)*qb(2,iq) ))
               
               work(iow1) = work(iow1) + ph * GS(juo,iuo,iq)
               work(iow2) = work(iow2) + ph * &
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

    call EYE(no_s,Sigma)
    ! We have the matrix to invert in the first no_s**2 values.
    call zgesv(no_s,no_s,work(1),no_s,ipvt,Sigma,no_s,ierr)
    if ( ierr /= 0 ) then
       write(*,*) 'Inversion of surface Greens function failed'
    end if

    ! Do:
    ! \Sigma = Z*S - H - \Sigma_bulk
    iow2 = no_s**2
    do juo = 1 , no_s
       do iuo = 1 , no_s
          iow2 = iow2+1
          Sigma(iuo,juo) = work(iow2) - Sigma(iuo,juo)
       end do
    end do

    call timer('ts_expand',2)

  end subroutine UC_expansion_Sigma

  subroutine UC_expansion_Sigma_Gamma(UseBulk,ZEnergy,no_u,no_s,NRepA1,NRepA2, &
       na_u,lasto,nq,qb,wq,H,S,GS,Sigma,Gamma,nwork,work)
    use intrinsic_missing, only: EYE
    use units, only : Pi
! ********************
! * INPUT variables  *
! ********************
    complex(dp), intent(in) :: ZEnergy
    logical,  intent(in) :: UseBulk
    integer,  intent(in) :: no_u, no_s, NRepA1, NRepA2
    integer,  intent(in) :: na_u,lasto(0:na_u)
    integer,  intent(in) :: nq
    real(dp), intent(in) :: qb(3,nq), wq(nq)
    complex(dp), dimension(no_u,no_u,nq), intent(in) :: H, S, GS
! ********************
! * OUTPUT variables *
! ********************
    complex(dp), intent(out) :: Sigma(no_s,no_s)
    real(dp), intent(out)    :: Gamma(no_s,no_s)

    integer,  intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)
! ********************
! * LOCAL variables  *
! ********************
    integer :: iq, ierr
    integer :: iow1, iow2
    integer :: io,iau,ia2,ia1,iuo
    integer :: jo,jau,ja2,ja1,juo
    integer :: ipvt(no_s)
    real(dp), parameter :: Pi2 = 2._dp * Pi
    complex(dp) :: ph

    call timer('ts_expand',1)
    call timer('ts_expandG',1)
    
    ! THis should never happen (work is TS-region!)
    if ( nwork < no_s**2*2 ) call die('Size of work-array is &
         &too small')

    do iow1 = 1 , no_s*no_s*2
       work(iow1)      = dcmplx(0._dp,0._dp)
    end do
    do jo = 1 , no_s
       do io = 1 , no_s
          Gamma(io,jo) = 0._dp
       end do
    end do

    ! This is the crucial calcuation.
    ! If we use bulk values in the electrodes
    ! we need not add the expanded H and S values to get the 
    ! electrode \Sigma. Hence, we need only expand
    ! surface Green's function
    do iq = 1 , nq
     iow1 = 0
     iow2 = no_s**2
     do iau = 1 , na_u
      do ia2 = 0 , NRepA2-1
       do ia1 = 0 , NRepA1-1
        do iuo = 1 + lasto(iau-1) , lasto(iau)   
         do jau = 1 , na_u
          do ja2 = 0 , NRepA2-1
           do ja1 = 0 , NRepA1-1
            do juo = 1 + lasto(jau-1) , lasto(jau)  
               iow1 = iow1 + 1
               iow2 = iow2 + 1
               ph = wq(iq) * cdexp(dcmplx(0.0_dp,-Pi2) * &
                    ( (ia1-ja1)*qb(1,iq) + (ia2-ja2)*qb(2,iq) ))
               
               work(iow1) = work(iow1) + ph * GS(juo,iuo,iq)

               work(iow2) = work(iow2) + ph * &
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

    call EYE(no_s,Sigma)
    ! We have the matrix to invert in the first no_s**2 values.
    call zgesv(no_s,no_s,work(1),no_s,ipvt,Sigma,no_s,ierr)
    if ( ierr /= 0 ) then
       write(*,*) 'Inversion of surface Greens function failed'
    end if

    if ( UseBulk ) then

       ! Do:
       ! \Sigma = Z*S - H - \Sigma_bulk
       iow1 = 0 
       iow2 = no_s**2
       do jo = 1 , no_s
          do io = 1 , no_s
             iow1 = iow1 + 1
             work(iow1) = work(iow2+iow1) - Sigma(io,jo)
          end do
       end do

       ! Do
       ! \Gamma = -\Im \Sigma
       do jo = 1 , no_s
          iow1 = (jo-1) * no_s
          do io = 1 , jo
             iow2 = (io-1) * no_s
             Gamma(io,jo) = -0.5_dp * dimag( &
                  work(io+iow1)-dconjg(work(jo+iow2)) )
             Gamma(jo,io) = -0.5_dp * dimag( &
                  work(jo+iow2)-dconjg(work(io+iow1)) )
          end do
       end do

    else
       
       ! Do:
       ! \Sigma = Z*S - H - \Sigma_bulk
       iow2 = no_s**2
       do jo = 1 , no_s
          do io = 1 , no_s
             iow2 = iow2 + 1
             Sigma(io,jo) = work(iow2) - Sigma(io,jo)
          end do
       end do

       ! Do
       ! \Gamma = \Im \Sigma
       do jo = 1 , no_s
          do io = 1 , jo
             Gamma(io,jo) = -0.5_dp * dimag( &
                  Sigma(io,jo)-dconjg(Sigma(jo,io)) )
             Gamma(jo,io) = -0.5_dp * dimag( &
                  Sigma(jo,io)-dconjg(Sigma(io,jo)) )
          end do
       end do

    end if

    call timer('ts_expandG',2)
    call timer('ts_expand',2)
       
  end subroutine UC_expansion_Sigma_Gamma


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
    complex(dp), dimension(no_GS,no_GS,nq), intent(inout) :: HAA, SAA, GAA
    complex(dp), intent(in) :: ZEnergy
    ! The work array passed, this means we do not have
    ! to allocate anything down here.
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork) ! or is this a 2D-array?

! *********************
! * LOCAL variables   *
! *********************
    real(dp), parameter :: EPS = 1.d-7
    logical, save :: first_call = .true.
    integer, save :: icall = 1
    integer :: read_Size

#ifdef MPI
    integer :: MPIerror, Request, Status(MPI_Status_Size)
#endif

    complex(dp) :: ZE_cur, wGFi
    integer :: iEni, iNode, ikGS

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE read_next_GS' )
#endif

    call timer('rn_GS',1)

    if ( first_call ) then
       icall = 1
       first_call = .false.
    else
       icall = icall + 1
    end if

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

          ! We will ensure that we have passed on the information here
          ! This is for ZE_cur
          call MPI_Wait(Request,Status,MPIerror)
          
          read(uGF) work(1:read_Size)
          
          call MPI_ISend(work(1),read_Size,MPI_Double_Complex, &
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


  subroutine GF_Gamma_GF(Offset,no_u_TS,no_E,GF,Gamma,GGG,nwork,work)

! ======================================================================
!  This routine returns GGG=(-i)*GF.Gamma.GF, where GF is a (no_u)x(no_u)
!  matrix and the first[LEFTFLAG true]/last[LEFTFLAG false] states
!  corresponds to the (no_E) Left/Right electrode states.
!  Gamma is a (no_E)x(no_E) matrix.
! ======================================================================
    use precision, only : dp

    implicit none

! *********************
! * INPUT variables   *
! *********************
    integer, intent(in) :: Offset ! this will help leverage the loop construct and
!                                 ! for the right electrode this is the no_u_TS - no_E + 1
    integer, intent(in) :: no_u_TS !no. states in contact region
    integer, intent(in) :: no_E ! the size of the Gamma
    complex(dp), intent(in) :: GF(no_u_TS,no_u_TS)
    ! i (Sigma - Sigma^dagger)/2
    real(dp),    intent(in) :: Gamma(Offset:Offset+no_E-1,Offset:Offset+no_E-1)
    ! A work array for doing the calculation... (nwork has to be larger than no_u_TS)
    integer,     intent(in)    :: nwork
    complex(dp), intent(inout) :: work(nwork)

! *********************
! * OUTPUT variables  *
! *********************
    complex(dp), intent(out) :: GGG(no_u_TS,no_u_TS)    !GF.GAMMA.GF

    integer :: i, j, ie, je, lE

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE GFGammaGF' )
#endif

    call timer("GFGGF",1)

    ! calculate the end of the bound for the Gamma
    lE = Offset + no_E - 1
    if ( lE /= no_u_TS .and. lE /= no_E ) call die('Wrong indices')

    do i = 1 , no_u_TS

       do ie = Offset , lE

          work(ie) = dcmplx(0._dp,0._dp)
          do je = Offset , lE
             work(ie) = work(ie) + Gamma(je,ie) * GF(i,je)
          end do

       end do
       
       do j = 1 , i

          GGG(j,i) = dcmplx(0._dp,0._dp)

          do ie = Offset , lE
             GGG(j,i) = GGG(j,i) + work(ie) * dconjg(GF(j,ie))
          end do

          ! We only take the real part once (faster)
          GGG(j,i) = dreal(GGG(j,i)) * dcmplx(0._dp,-1._dp)
          GGG(i,j) = GGG(j,i)

       end do
       
    end do

    call timer('GFGGF',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS GFGammaGF' )
#endif

! ====================================================================
  END subroutine GF_Gamma_GF

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
    use class_dSpArr1D

    implicit none

! *********************
! * INPUT variables   *
! *********************
    integer, intent(in) :: no_C_L, no_C_R
! *********************
! * OUTPUT variables  *
! *********************
    ! Contour part of DM integration
    type(dSpArr1D), intent(inout) :: SpArrDML, SpArrDMR
    ! Real-axis part of DM integration
    type(dSpArr1D), intent(inout) :: SpArrDMneqL, SpArrDMneqR
    ! L-R estimates of EDM
    type(dSpArr1D), intent(inout) :: SpArrEDML, SpArrEDMR

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

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE weightDM' )
#endif

    call timer('weightDM',1)

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

    do io = 1 , nr
       ! We are in a buffer region...
       if ( l_ncol(io) == 0 ) cycle
       do j = 1 , l_ncol(io)

          ind = l_ptr(io)+j
          ! Retrive the connecting orbital
          jo = l_col(ind)

! For now we do not allow tri-diagonalization in the memory reduced method
! Otherwise this has to be changed.
          if ( io < no_C_L .or. jo < no_C_L ) then
             wL = 0._dp
             wR = 1._dp
          else if ( no_C_R < io .or. no_C_R < jo ) then
             wL = 1._dp
             wR = 0._dp
          else
             ! The weights
             wL = DMneqL(ind)*DMneqL(ind)
             wSUM = wL + DMneqR(ind)*DMneqR(ind)
             if ( wSUM > 0._dp ) then
                wL = wL / wSUM
                wR = 1._dp - wL
             else
                wL = 0.5_dp
                wR = 0.5_dp
             end if
          end if
          DML(ind) = wL * (DML(ind) + DMneqR(ind)) &
               + wR * (DMR(ind) + DMneqL(ind))
          EDML(ind) = wL * EDML(ind) + wR * EDMR(ind)
       end do
    end do

    call timer('weightDM',2)

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
    use class_zSpArr1D

    implicit none

! *********************
! * INPUT variables   *
! *********************
    integer, intent(in) :: no_C_L, no_C_R
! *********************
! * OUTPUT variables  *
! *********************
    ! Contour part of DM integration
    type(zSpArr1D), intent(inout) :: SpArrDML, SpArrDMR
    ! Real-axis part of DM integration
    type(zSpArr1D), intent(inout) :: SpArrDMneqL, SpArrDMneqR
    ! L-R estimates of EDM
    type(zSpArr1D), intent(inout) :: SpArrEDML, SpArrEDMR

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

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE weightDMC' )
#endif

    call timer('weightDMC',1)

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

    do io = 1 , nr
       ! We are in a buffer region...
       if ( l_ncol(io) == 0 ) cycle
       do j = 1 , l_ncol(io)

          ind = l_ptr(io)+j
          ! Retrive the connecting orbital
          jo = l_col(ind)
! For now we do not allow tri-diagonalization in the memory reduced method
! Otherwise this has to be changed.
          if ( io < no_C_L .or. jo < no_C_L ) then
             wL = 0._dp
             wR = 1._dp
          else if ( no_C_R < io .or. no_C_R < jo ) then
             wL = 1._dp
             wR = 0._dp
          else
             ! The weights
             wL = dconjg(DMneqL(ind))*DMneqL(ind)
             wSUM = wL + dconjg(DMneqR(ind))*DMneqR(ind)
             if ( wSUM > 0._dp ) then
                wL = wL/wSUM
                wR = 1._dp - wL
             else
                wL = 0.5_dp
                wR = 0.5_dp
             end if
          end if
          DML(ind) = wL * (DML(ind) + DMneqR(ind)) &
               + wR * (DMR(ind) + DMneqL(ind))
          EDML(ind) = wL * EDML(ind) + wR * EDMR(ind)
       end do
    end do

    call timer('weightDMC',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS weightDMC' )
#endif

  end subroutine weightDMC

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
