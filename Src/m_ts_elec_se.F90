! This module will control what is used of the electrodes in the transiesta SCF
!
! Hence we here collect the routines for reading and expanding the self-energies
! in the GF-files.

module m_ts_elec_se

  use precision, only : dp

  implicit none

  private

  public :: UC_expansion
  public :: UC_expansion_Sigma_Bulk
  public :: UC_expansion_Sigma
  public :: UC_expansion_Sigma_GammaT
  public :: read_next_GS
  public :: read_next_GS_LR

contains

  subroutine UC_expansion(non_Eq,UseBulk,ZEnergy, &
       no_u,no_s,NRepA1,NRepA2, &
       na_u,lasto,nq,qb,wq, &
       H,S,GS,Sigma,Gamma,nwork,work)
! ********************
! * INPUT variables  *
! ********************
    logical,  intent(in) :: non_Eq, UseBulk
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
    complex(dp), intent(inout) :: Sigma(no_s,no_s)
    complex(dp), intent(inout) :: Gamma(no_s,no_s)

    integer,  intent(in) :: nwork
    complex(dp), intent(inout) :: work(no_s,no_s,2)

    if ( non_Eq ) then
       call UC_expansion_Sigma_GammaT(UseBulk,ZEnergy, &
            no_u,no_s,NRepA1,NRepA2, &
            na_u,lasto,nq,qb,wq,H,S,GS,Sigma,Gamma,nwork,work)
    else
       if ( UseBulk ) then
          call UC_expansion_Sigma_Bulk(no_u,no_s,NRepA1,NRepA2, &
               na_u,lasto,nq,qb,wq,H,S,GS,Sigma,nwork,work)
       else
          call UC_expansion_Sigma(ZEnergy,no_u,no_s,NRepA1,NRepA2, &
               na_u,lasto,nq,qb,wq,H,S,GS,Sigma,nwork,work)
       end if
    end if

  end subroutine UC_expansion

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
    complex(dp), intent(inout) :: Sigma(no_s,no_s)

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
    complex(dp), intent(inout) :: Sigma(no_s,no_s)

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
    complex(dp), intent(inout) :: Sigma(no_s,no_s)
    complex(dp), intent(inout) :: GammaT(no_s,no_s)

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

  subroutine read_next_GS_LR(uGF,NEReqs,ikpt,no_GS,nq,HAA,SAA,GAA,ZEnergy, &
       nwork,work)

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

  end subroutine read_next_GS_LR


  ! Subroutine for reading in both the left and right next energy point
  
  subroutine read_next_GS(iPE,cNEn,Z,ikpt, &
       uGFL, no_L_HS, nqL, HAAL, SAAL, GAAL, &
       uGFR, no_R_HS, nqR, HAAR, SAAR, GAAR, &
       nzwork, zwork)

    use parallel, only : Node, Nodes
    use m_ts_cctype
      
    integer, intent(in) :: iPE, cNEn, ikpt
    complex(dp), intent(in) :: Z
    integer, intent(in) :: uGFL, no_L_HS, nqL
    complex(dp), intent(inout), dimension(no_L_HS,no_L_HS,nqL) :: HAAL, SAAL, GAAL
    integer, intent(in) :: uGFR, no_R_HS, nqR
    complex(dp), intent(inout), dimension(no_R_HS,no_R_HS,nqR) :: HAAR, SAAR, GAAR
    integer, intent(in) :: nzwork
    complex(dp), intent(inout) :: zwork(nzwork)

    integer :: iE, NEReqs

    ! obtain a valid energy point (truncate at NEn)
    iE = min(iPE,cNEn)
    
    ! save the current weight of the point
    ! This is where we include the factor-of-two for spin and
    ! and the (1/Pi) from DM = Im[G]/Pi
    ! Furthermore we include the weight of the k-point

    ! the number of points we wish to read in this segment
    NEReqs = min(Nodes, cNEn-(iPe-1-Node))

    ! TODO Move reading of the energy points
    ! directly into the subroutines which need them
    ! In this way we can save both GAA, Sigma AND Gamma arrays!!!!
    ! However, this will probably come at the expense 
    ! of doing the same "repetition" expansion twice, we can live with
    ! that!

    ! Read in the left electrode
    call read_next_GS_LR(uGFL, NEReqs, &
         ikpt,no_L_HS,nqL, HAAL, SAAL, &
         GAAL, Z, nzwork, zwork)
    
    ! Read in the right electrode
    call read_next_GS_LR(uGFR, NEReqs, &
         ikpt,no_R_HS,nqR, HAAR, SAAR, &
         GAAR, Z, nzwork, zwork)

  end subroutine read_next_GS

end module m_ts_elec_se
