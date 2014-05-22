! This module will control what is used of the electrodes in the transiesta SCF
!
! Hence we here collect the routines for reading and expanding the self-energies
! in the GF-files.

module m_ts_elec_se

  use precision, only : dp

  use m_ts_electype
  use m_ts_cctype

  implicit none

  private

  public :: UC_expansion
  !public :: UC_expansion_Sigma_Bulk
  !public :: UC_expansion_Sigma
  !public :: UC_expansion_Sigma_GammaT
  public :: read_next_GS
!  public :: read_next_GS_LR

contains

  subroutine UC_expansion(cE,El,nwork,work, &
       non_Eq)
! ********************
! * INPUT variables  *
! ********************
    type(ts_c_idx), intent(in) :: cE
    type(Elec), intent(in out) :: El

! ********************
! * WORK variables   *
! ********************
    integer,  intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)

    logical,  intent(in), optional :: non_Eq

    logical :: lnon_Eq

    if ( cE%fake ) return

    call timer('ts_expand',1)

    lnon_Eq = .false.
    if ( present(non_Eq) ) lnon_Eq = non_Eq

    if ( lnon_Eq ) then
       call UC_expansion_Sigma_GammaT(cE%e, &
            El%no_used,TotUsedOrbs(El),El, &
            El%na_used,El%lasto_used,Rep(El), &
            El%HA,El%SA,El%GA,El%Sigma,El%Gamma,nwork,work)
    else
       if ( El%Bulk ) then
          call UC_expansion_Sigma_Bulk(El%no_used,TotUsedOrbs(El),El, &
               El%na_used,El%lasto_used,Rep(El), &
               El%HA,El%SA,El%GA,El%Sigma,nwork,work)
       else
          call UC_expansion_Sigma(cE%e,El%no_used,TotUsedOrbs(El),El, &
               El%na_used,El%lasto_used,Rep(El), &
               El%HA,El%SA,El%GA,El%Sigma,nwork,work)
       end if
    end if

    call timer('ts_expand',2)

  end subroutine UC_expansion

  subroutine UC_expansion_Sigma_Bulk(no_u,no_s,El, &
       na_u,lasto,nq,H,S,GS,Sigma,nwork,work)
    use intrinsic_missing, only : EYE
    use units, only : Pi
! ********************
! * INPUT variables  *
! ********************
    integer,  intent(in) :: no_u, no_s
    type(Elec), intent(in) :: El
    integer,  intent(in) :: na_u,lasto(0:na_u)
    integer,  intent(in) :: nq
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
    integer :: iow,iau,ia3,ia2,ia1,iuo
    integer :: jow,jau,ja3,ja2,ja1,juo
    integer :: ipvt(no_s)
    complex(dp), parameter :: zmPi2 = dcmplx(0._dp,-2._dp * Pi)
    complex(dp), parameter :: zPi2  = dcmplx(0._dp, 2._dp * Pi)
    complex(dp) :: ph
    real(dp) :: qmPi(3,nq), wq

    ! THis should never happen (work is TS-region!)
    if ( nwork < no_s**2 ) call die('Size of work-array is &
         &too small')

    if ( nq == 1 ) then
       if ( no_u /= no_s ) call die('no_E/=no_s')

       ! When no repetition we save it "as is"
       Sigma(:,:) = GS(:,:,1)

    else

       do iq = 1 , nq 
          qmPi(1:3,iq) = - 2._dp * Pi * q_exp(El,iq)
       end do
       wq = 1._dp / nq

       ! This is the crucial calcuation.
       ! If we use bulk values in the electrodes
       ! we need not add the expanded H and S values to get the 
       ! electrode \Sigma. Hence, we need only expand
       ! surface Green's function
       ! We do the equivalent of :
       ! ph = wq(iq) * cdexp(zmPi2 * &
       !      ( (ia1-ja1)*qb(1,iq) + (ia2-ja2)*qb(2,iq) ))
       iq = 1
       iow = 0
       do iau = 1 , na_u
        do ia3 = 1 , El%RepA3
        do ia2 = 1 , El%RepA2
        do ia1 = 1 , El%RepA1
          do iuo = 1 + lasto(iau-1) , lasto(iau)
           iow = iow + 1
           jow = 0
           do jau = 1 , na_u
            do ja3 = 1 , El%RepA3
            do ja2 = 1 , El%RepA2
            do ja1 = 1 , El%RepA1
              ph = wq * cdexp(dcmplx(0._dp, &
                   (ia1-ja1)*qmPi(1,iq) + (ia2-ja2)*qmPi(2,iq) + (ia3-ja3)*qmPi(3,iq) ) )
              do juo = 1 + lasto(jau-1) , lasto(jau) 
                 jow = jow + 1
                 
                 work(jow,iow) = ph * GS(juo,iuo,iq)
              end do !juo
            end do !ja1
            end do !ja2
            end do !ja3
           end do !jau
          end do !iuo
        end do !ia1
        end do !ia2
        end do !ia3
       end do !iau
       do iq = 2 , nq
        iow = 0
        do iau = 1 , na_u
         do ia3 = 1 , El%RepA3
         do ia2 = 1 , El%RepA2
         do ia1 = 1 , El%RepA1
           do iuo = 1 + lasto(iau-1) , lasto(iau)
            iow = iow + 1
            jow = 0
            do jau = 1 , na_u
             do ja3 = 1 , El%RepA3
             do ja2 = 1 , El%RepA2
             do ja1 = 1 , El%RepA1
               ph = wq * cdexp(dcmplx(0._dp, &
                    (ia1-ja1)*qmPi(1,iq) + (ia2-ja2)*qmPi(2,iq) + (ia3-ja3)*qmPi(3,iq) ) )
               do juo = 1 + lasto(jau-1) , lasto(jau) 
                  jow = jow + 1
                  
                  work(jow,iow) = work(jow,iow) + ph * GS(juo,iuo,iq)
               end do !juo
             end do !ja1
             end do !ja2
             end do !ja3
            end do !jau
           end do !iuo
         end do !ia1
         end do !ia2
         end do !ia3
        end do !iau
       end do !q-points
       
       call EYE(no_s,Sigma)

       ! We have the matrix to invert in the first no_s**2 values.
       call zgesv(no_s,no_s,work(1,1),no_s,ipvt,Sigma,no_s,ierr)
       if ( ierr /= 0 ) &
            write(*,*) 'Inversion of surface Greens function failed'
       
    end if

  end subroutine UC_expansion_Sigma_Bulk


  subroutine UC_expansion_Sigma(ZEnergy,no_u,no_s,El, &
       na_u,lasto,nq,H,S,GS,Sigma,nwork,work)
    use intrinsic_missing, only : EYE
! ********************
! * INPUT variables  *
! ********************
    complex(dp), intent(in) :: ZEnergy
    integer,  intent(in) :: no_u, no_s
    type(Elec), intent(in) :: El
    integer,  intent(in) :: na_u,lasto(0:na_u)
    integer,  intent(in) :: nq
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

    ! THis should never happen (work is TS-region!)
    if ( nwork < no_s**2*2 ) call die('Size of work-array is &
         &too small')
       
    call update_UC_expansion(ZEnergy,no_u,no_s,El, &
         na_u,lasto,nq,H,S,GS,nwork,work(1,1,1))

    if ( nq == 1 ) then

       ! When no repetition we save it "as is"
       Sigma(:,:) = GS(:,:,1)

    else
       call EYE(no_s,Sigma)

       ! We have the matrix to invert in the first no_s**2 values.
       call zgesv(no_s,no_s,work(1,1,1),no_s,ipvt,Sigma,no_s,ierr)

       if ( ierr /= 0 ) &
            write(*,*) 'Inversion of surface Greens function failed'

    end if

    ! Do:
    ! \Sigma = Z*S - H - \Sigma_bulk
    do jo = 1 , no_s
       do io = 1 , no_s
          Sigma(io,jo) = work(io,jo,2) - Sigma(io,jo)
       end do
    end do

  end subroutine UC_expansion_Sigma

  subroutine UC_expansion_Sigma_GammaT(ZEnergy,no_u,no_s,El, &
       na_u,lasto,nq,H,S,GS,Sigma,GammaT,nwork,work)
    use intrinsic_missing, only: EYE
! ********************
! * INPUT variables  *
! ********************
    complex(dp), intent(in) :: ZEnergy
    integer,  intent(in) :: no_u, no_s
    type(Elec), intent(in) :: El
    integer,  intent(in) :: na_u,lasto(0:na_u)
    integer,  intent(in) :: nq
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

    ! THis should never happen (work is TS-region!)
    if ( nwork < no_s**2*2 ) call die('Size of work-array is &
         &too small')

    call update_UC_expansion(ZEnergy,no_u,no_s,El, &
         na_u,lasto,nq,H,S,GS,nwork,work(1,1,1))

    if ( nq == 1 ) then

       ! When no repetition we save it "as is"
       Sigma(:,:) = GS(:,:,1)
       
    else
       call EYE(no_s,Sigma)

       ! We have the matrix to invert in the first no_s**2 values.
       call zgesv(no_s,no_s,work(1,1,1),no_s,ipvt,Sigma,no_s,ierr)

       if ( ierr /= 0 ) then
          write(*,*) 'Inversion of surface Greens function failed'
       end if

    end if

    if ( El%Bulk ) then

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

  end subroutine UC_expansion_Sigma_GammaT

  subroutine update_UC_expansion(ZEnergy,no_u,no_s,El, &
       na_u,lasto,nq,H,S,GS,nwork,work)
    use units, only : Pi
! ********************
! * INPUT variables  *
! ********************
    complex(dp), intent(in) :: ZEnergy
    integer,  intent(in) :: no_u, no_s
    type(Elec), intent(in) :: El
    integer,  intent(in) :: na_u,lasto(0:na_u), nq
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
    integer :: iow,iau,ia3,ia2,ia1,iuo
    integer :: jow,jau,ja3,ja2,ja1,juo
    complex(dp), parameter :: zmPi2 = dcmplx(0._dp,-2.0_dp * Pi)
    complex(dp), parameter :: zPi2  = dcmplx(0._dp, 2.0_dp * Pi)
    complex(dp) :: ph
    real(dp) :: qmPi(3,nq), wq

    if ( nq == 1 ) then
       if ( no_u /= no_s ) call die('no_E/=no_s')

       ! We do not need this...
       !work(:,:,1) = GS(:,:,1)
       work(:,:,2) = ZEnergy * S(:,:,1) - H(:,:,1)

    else

       do iq = 1 , nq 
          qmPi(1:3,iq) = - 2._dp * Pi * q_exp(El,iq)
       end do
       wq = 1._dp / real(nq,dp)

       ! This is the crucial calcuation.
       ! If we use bulk values in the electrodes
       ! we need not add the expanded H and S values to get the 
       ! electrode \Sigma. Hence, we need only expand
       ! surface Green's function
       iq = 1
       iow = 0
       do iau = 1 , na_u
        do ia3 = 1 , El%RepA3
        do ia2 = 1 , El%RepA2
        do ia1 = 1 , El%RepA1
          do iuo = 1 + lasto(iau-1) , lasto(iau)
           iow = iow + 1
           jow = 0
           do jau = 1 , na_u
            do ja3 = 1 , El%RepA3
            do ja2 = 1 , El%RepA2
            do ja1 = 1 , El%RepA1
              ph = wq * cdexp(dcmplx(0._dp, &
                   (ia1-ja1)*qmPi(1,iq) + (ia2-ja2)*qmPi(2,iq) + (ia3-ja3)*qmPi(3,iq) ) )
              do juo = 1 + lasto(jau-1) , lasto(jau)
                jow = jow + 1
                
                work(jow,iow,1) = ph * GS(juo,iuo,iq)
                
                work(jow,iow,2) = ph * (ZEnergy*S(juo,iuo,iq)-H(juo,iuo,iq))
                
              end do !juo
            end do !ja1
            end do !ja2
            end do !ja3
           end do !jau
          end do !iuo
        end do !ia1
        end do !ia2
        end do !ia3
       end do !iau
       do iq = 2 , nq
        iow = 0
        do iau = 1 , na_u
         do ia3 = 1 , El%RepA3
         do ia2 = 1 , El%RepA2
         do ia1 = 1 , El%RepA1
           do iuo = 1 + lasto(iau-1) , lasto(iau)
            iow = iow + 1
            jow = 0
            do jau = 1 , na_u
             do ja3 = 1 , El%RepA3
             do ja2 = 1 , El%RepA2
             do ja1 = 1 , El%RepA1
               ph = wq * cdexp(dcmplx(0._dp, &
                    (ia1-ja1)*qmPi(1,iq) + (ia2-ja2)*qmPi(2,iq) + (ia3-ja3)*qmPi(3,iq) ) )
               do juo = 1 + lasto(jau-1) , lasto(jau)
                  jow = jow + 1
                  
                  work(jow,iow,1) = work(jow,iow,1) + ph * GS(juo,iuo,iq)
   
                  work(jow,iow,2) = work(jow,iow,2) + ph * &
                       (ZEnergy*S(juo,iuo,iq)-H(juo,iuo,iq))
   
               end do !juo
             end do !ja1
             end do !ja2
             end do !ja3
            end do !jau
           end do !iuo
         end do !ia1
         end do !ia2
         end do !ia3
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

  subroutine read_next_GS_Elec(uGF,NEReqs,ikpt,El,cE, &
       nwork,work, &
       forward)

    use parallel, only : IONode, Node, Nodes
    use units,    only : eV

#ifdef MPI
    use mpi_siesta, only : MPI_Bcast, MPI_Isend, MPI_Irecv
    use mpi_siesta, only : MPI_Sum, MPI_Integer, MPI_Double_Complex
    use mpi_siesta, only : MPI_Status_Size, MPI_Comm_World
#endif

! *********************
! * INPUT variables   *
! *********************
    ! file-unit, and k-point index
    integer, intent(in) :: uGF, ikpt
    integer, intent(in) :: NEReqs
    ! The electrode also contains the arrays
    type(Elec), intent(in out) :: El
    type(ts_c_idx), intent(in)     :: cE
    ! The work array passed, this means we do not have
    ! to allocate anything down here.
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)
    ! If 'forward' is .true. we will read consecutively 
    ! and distribute from Node = 0, to Node = Nodes - 1
    ! (default)
    ! 
    ! if 'forward' is .false. we will read consecutively 
    ! and distribute from Node = Nodes - 1 to Node = 0
    logical, intent(in), optional :: forward

! *********************
! * LOCAL variables   *
! *********************
    real(dp), parameter :: EPS = 1.d-7
    integer :: read_Size

#ifdef MPI
    integer :: MPIerror, Request, Status(MPI_Status_Size)
#endif

    complex(dp) :: ZE_cur
    integer :: iNode, ikGS, iEni
    integer :: iNodeS, iNodeE, iNodeStep
    logical :: lforward

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE read_next_GS' )
#endif

    lforward = .true.
    if ( present(forward) ) lforward = forward
    if ( lforward ) then
       iNodeS = 0
       iNodeE = NEReqs - 1
       iNodeStep = 1
    else
       iNodeS = Nodes - 1
       iNodeE = Nodes - NEReqs
       iNodeStep = -1
    end if

    read_Size = El%no_used ** 2 * Rep(El) ! no_GS * no_GS * nq

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
    do iNode = iNodeS, iNodeE, iNodeStep

       if ( IONode ) then
          ! read in header of GF-point
          read(uGF) ikGS, iEni, ZE_cur
       end if

#ifdef MPI
       call MPI_Bcast(iEni,1,MPI_Integer, &
            0,MPI_Comm_World,MPIerror)

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
          call MPI_Wait(Request,Status,MPIerror)
       end if
#endif

       ! The test of the energy-point is performed on
       ! the calculating node...
       if ( Node == iNode ) then
          if ( cdabs(cE%e-ZE_cur) > EPS ) then
             write(*,*) 'GF-file: '//trim(GFFile(El))
             write(*,'(2(a,2(tr1,g12.5)))') 'Energies, TS / Gf:', &
                  cE%e / eV, ' /', ZE_cur / eV
             call die('Energy point in GF file does &
                  not match the internal energy-point in transiesta. &
                  &Please correct your GF files.')
          end if
       end if
       
       ! If the k-point does not match what we expected...
       if ( IONode .and. ikpt /= ikGS ) then
          write(*,*) 'GF-file: '//trim(GFFile(El))
          write(*,'(2(a,i0))') 'k-point, TS / Gf: ', &
               ikpt, ' / ', ikGS
          call die('Read k-point in GF file does not match &
               &the requested k-point. Please correct your &
               &GF files.')
       end if

       if ( iEni == 1 ) then

          ! read in the electrode Hamiltonian and overlap...
          if ( IONode ) then
             read(uGF) El%HA
             read(uGF) El%SA
          end if

#ifdef MPI
          call MPI_Bcast(El%HA(1,1,1),read_Size,MPI_Double_Complex, &
               0,MPI_Comm_World,MPIerror)
          call MPI_Bcast(El%SA(1,1,1),read_Size,MPI_Double_Complex, &
               0,MPI_Comm_World,MPIerror)
#endif
       end if 


#ifdef MPI
       if ( IONode .and. iNode == Node ) then
#endif
          ! read in surface Green's function
          read(uGF) El%GA
#ifdef MPI
       else if ( IONode ) then

          read(uGF) work(1:read_Size)
          
          call MPI_ISend(work(1),read_Size,MPI_Double_Complex, &
               iNode,1,MPI_Comm_World,Request,MPIerror) 
          call MPI_Wait(Request,Status,MPIerror)
       else if ( Node ==  iNode ) then
          call MPI_IRecv(El%GA(1,1),read_Size,MPI_Double_Complex, &
               0    ,1,MPI_Comm_World,Request,MPIerror)
          call MPI_Wait(Request,Status,MPIerror)
       end if
#endif

    end do

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS read_next_GS' )
#endif

  end subroutine read_next_GS_Elec


  ! Subroutine for reading in both the left and right next energy point
  
  subroutine read_next_GS(ispin,ikpt, kpt, cE, &
       NElecs, uGF, Elecs, &
       nzwork, zwork, RemUCellDistance, reread, &
       forward )

    use parallel, only : Node, Nodes, IONode

#ifdef MPI
    use mpi_siesta, only : MPI_AllReduce, MPI_Sum, MPI_Integer
    use mpi_siesta, only : MPI_Comm_World
#endif

    use m_ts_cctype
    use m_ts_electrode, only: calc_next_GS_Elec
      
    integer, intent(in) :: ispin, ikpt
    real(dp), intent(in) :: kpt(3)
    type(ts_c_idx), intent(in) :: cE
    integer, intent(in) :: NElecs, uGF(NElecs)
    type(Elec), intent(inout) :: Elecs(NElecs)
    integer, intent(in) :: nzwork
    complex(dp), intent(inout), target :: zwork(nzwork)
    logical, intent(in) :: RemUCellDistance
    logical, intent(in), optional :: reread, forward

    integer :: iE, NEReqs, i, j
#ifdef MPI
    integer :: MPIerror
#endif

    ! save the current weight of the point
    ! This is where we include the factor-of-two for spin and
    ! and the (1/Pi) from DM = Im[G]/Pi
    ! Furthermore we include the weight of the k-point

    ! the number of points we wish to read in this segment
    NEReqs = 1
#ifdef MPI
    i = 0
    if ( cE%exist .and. .not. cE%fake) i = 1
    call MPI_AllReduce(i,NEReqs,1,MPI_Integer, MPI_Sum, &
         MPI_Comm_World, MPIerror)
#endif

    if ( present(reread) ) then
       if ( IONode .and. reread ) then
          do j = 1 , NElecs
             if ( .not. Elecs(j)%out_of_core ) cycle
             do i = 1 , NEReqs * 2
                backspace(unit=uGF(j))
             end do
          end do
! Currently the equilibrium energy points are just after
! the k-point, hence we will never need to backspace behind the
! HAA and SAA reads
! however, when we add kpoints for the bias contour, then it might be
! necessary!
!           if ( new_kpt ) then
!             do i = 1 , 2
!                backspace(unit=uGFL)
!                backspace(unit=uGFR)
!             end do
!          end if
       end if
    end if

    ! TODO Move reading of the energy points
    ! directly into the subroutines which need them
    ! In this way we can save both GAA, Sigma AND Gamma arrays!!!!
    ! However, this will probably come at the expense 
    ! of doing the same "repetition" expansion twice, we can live with
    ! that!
    do i = 1 , NElecs
       if ( Elecs(i)%out_of_core ) then
          call read_next_GS_Elec(uGF(i), NEReqs, &
               ikpt, Elecs(i), cE, &
               nzwork, zwork, forward = forward)
       else
          ! the electrode already have the correct 
          ! variables which decides whether it is
          ! RemUCellDistance or not!
          call calc_next_GS_Elec(Elecs(i),ispin,kpt,cE%e, &
               nzwork, zwork, RemUCellDistance)
       end if
    end do

  end subroutine read_next_GS

end module m_ts_elec_se
