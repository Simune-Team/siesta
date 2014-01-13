module m_ts_electrode
!
! Routines that are used for Electrodes GFs calculations
! Heavily updated by Nick Papior Andersen, 2012
!
!=============================================================================
! CONTAINS:
!          1) surface_Green
!          2) create_Green
!          3) init_electrode_HS
!          4) set_electrode_HS_Transfer

! TODO
! Remove all references to Gamma for the electrode.
! The TranSIESTA routine will never create the surface Green's
! function without having a transfer matrix.
! Thus we ENFORCE Gamma == .false. and the program should die if
! the electrode calculation was a Gamma calculation!

  use precision, only : dp

  implicit none

  public :: create_Green
  public :: init_Electrode_HS

  private

  ! Accuracy of the surface-Green's function
  real(dp), parameter :: accur = 1.e-15_dp
  ! BLAS parameters
  complex(dp), parameter :: z_1  = dcmplx(1._dp,0._dp)
  complex(dp), parameter :: z_m1 = dcmplx(-1._dp,0._dp)
  complex(dp), parameter :: z_0  = dcmplx(0._dp,0._dp)

contains


  ! Calculates the surface Green's function for the electrodes
  ! Handles both the left and right one
  ! this is the Sancho, Sancho and Rubio algorithm
  subroutine SSR_sGreen_DOS(no,ZE,H00,S00,H01,S01,GS, &
       zDOS, &
       nwork, zwork, &
       iterations, final_invert)
       
! ***************** INPUT **********************************************
! integer     no      : Number of orbitals in the electrode
! complex(dp) ZE      : The energy of the Green's function evaluation
! complex(dp) H00     : Hamiltonian within the first unit cell (discarding T-direction)
! complex(dp) S00     : Overlap matrix within the first unit cell (discarding T-direction)
! complex(dp) H01     : Transfer matrix from H00 to the neighbouring cell (in T-direction)
! complex(dp) S01     : Transfer matrix from S00 to the neighbouring cell (in T-direction)
! ***************** OUTPUT *********************************************
! complex(dp) GS      : Surface Green's function of the electrode
! **********************************************************************
    use m_mat_invert
    use parallel,  only: IONode
    use precision, only: dp

! ***********************
! * INPUT variables     *
! ***********************
    integer,     intent(in) :: no
    complex(dp), intent(in) :: ZE 
    complex(dp), intent(in) :: H00(no*no),S00(no*no)
    complex(dp), intent(in) :: H01(no*no),S01(no*no)

    integer,     intent(in) :: nwork

    logical, intent(in), optional :: final_invert

! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), intent(out), target :: GS(no*no)
    complex(dp), pointer :: zwork(:)
    complex(dp), intent(out) :: zDOS

    integer, intent(out), optional :: iterations

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: nom1, no2, nosq
    integer :: ierr             !error in inversion
    integer :: i,j,ic,ic2
    logical :: as_first

    real(dp) :: ro

    ! on the stack...
    integer :: ipvt(no)
    complex(dp), dimension(:), pointer :: rh,rh1,w,alpha,beta,gb
    complex(dp), dimension(:), pointer :: gsL,gsR

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE surface_Green' )
#endif

    ! Initialize counter
    if ( present(iterations) ) iterations = 0

    call timer('ts_GS',1)

    nom1 = no - 1
    no2  = 2 * no
    nosq = no * no

    if ( nwork < 9 * nosq ) call die('surface_Green: &
         &Not enough work space')
    i = 0
    rh  => zwork(i+1:i+2*nosq) 
    i = i + 2*nosq
    rh1 => zwork(i+1:i+2*nosq) 
    i = i + 2*nosq
    alpha => zwork(i+1:i+nosq) 
    i = i + nosq
    beta => zwork(i+1:i+nosq) 
    i = i + nosq
    w => zwork(i+1:i+nosq)
    i = i + nosq
    gb => zwork(i+1:i+nosq) 
    i = i + nosq

    gsL => zwork(i+1:i+nosq) 
    gsR => GS


! gb    =   Z*S00-H00
! alpha = -(Z*S01-H01)
    do i = 1 , nosq
       gb(i)    = ZE * S00(i) - H00(i)
       alpha(i) = H01(i) - ZE * S01(i)
    end do

! gs  = Z*S00-H00
    do i = 1 , nosq
       gsL(i) = gb(i)
       gsR(i) = gb(i)
    end do

! beta = -(Z*S10-H10)
    do j = 1 , no
       ic = no * (j-1)
       do i = 1 , no
          ic2 = j + no*(i-1)
          beta(ic+i) = dconjg(H01(ic2)) - ZE * dconjg(S01(ic2))
       end do
    end do

    ! Initialize loop
    ro = accur + 1._dp
    as_first = .false.
    do while ( ro > accur ) 

       ! Increment iterations
       if ( present(iterations) ) &
            iterations = iterations + 1

! rh = -(Z*S01-H01) ,j<no
! rh = -(Z*S10-H10) ,j>no
       do i = 1, nosq
          rh(i)       = alpha(i)
          rh(nosq+i)  = beta(i)
       end do

! w = Z*S00-H00
       w(:) = gb(:)

! rh =  rh1^(-1)*rh
! rh =  t0
       call zgesv(no, no2, w, no, ipvt, rh, no, ierr)

       if ( ierr /= 0 ) then
          write(*,*) 'ERROR: calc_green 1 MATRIX INVERSION FAILED'
          write(*,*) 'ERROR: LAPACK INFO = ',ierr
       end if

       ! switch pointers instead of copying elements
       call switch_alpha_beta_rh1(as_first)

! alpha = -(Z*S01-H01)*t0
       call zgemm('N','N',no,no,no,z_1,rh1(1),no,rh(1),no,z_0,alpha,no)
! beta  = -(Z*S10-H10)*t0 ??
       call zgemm('N','N',no,no,no,z_1,rh1(nosq+1),no,rh(nosq+1),no,z_0,beta,no)

! ba    = (Z*S10-H10)*t0b
       call zgemm('N','N',no,no,no,z_m1,rh1(nosq+1),no,rh(1),no,z_0,w,no)
       do i = 1 , nosq
          gb(i)  = gb(i) + w(i)
          gsL(i) = gsL(i) + w(i)
       end do

! ab    = (Z*S01-H01)*t0
       call zgemm('N','N',no,no,no,z_m1,rh1(1),no,rh(nosq+1),no,z_0,w,no)
       ro = -1._dp
       do i = 1 , nosq
          gb(i)  = gb(i) + w(i)
          gsR(i) = gsR(i) + w(i)

          ! also update the criteria
          ro = max(ro,abs(w(i)))
       end do
       
    end do

    if ( present(final_invert) ) then
       ! If we do not need to invert it, save it for later.
       if ( .not. final_invert ) then
          rh1(1:nosq) = GS(:)
       end if
    end if

    ! Invert to get the Surface Green's function
    call mat_invert(gsL,w,no,MI_IN_PLACE_LAPACK, ierr=ierr)

    if ( ierr /= 0 ) then
       write(*,*) 'ERROR: calc_green GSL MATRIX INVERSION FAILED'
       write(*,*) 'ERROR: LAPACK INFO = ',ierr
    end if

    ! Invert to get the Surface Green's function
    call mat_invert(gsR,w,no,MI_IN_PLACE_LAPACK, ierr=ierr)

    if ( ierr /= 0 ) then
       write(*,*) 'ERROR: calc_green GSR MATRIX INVERSION FAILED'
       write(*,*) 'ERROR: LAPACK INFO = ',ierr
    end if

    ! Invert to obtain the bulk Green's function
    call mat_invert(GB,w,no,MI_IN_PLACE_LAPACK, ierr=ierr)

    if ( ierr /= 0 ) then
       write(*,*) 'ERROR: calc_green GB MATRIX INVERSION FAILED'
       write(*,*) 'ERROR: LAPACK INFO = ',ierr
    end if

    ! We now calculate the density of states...
    do i = 1 , nosq
       alpha(i) = H01(i) -        ZE  * S01(i)
       ! notice, we utilize the relation (H10-z*S10) = (H01-conjg(z)*S01)^H
       beta(i)  = H01(i) - dconjg(ZE) * S01(i)
    end do

    i = 1
    j = nosq + 1 
    ! zDOS = Tr{ G_b * S00 + 
    !            G_l * (H01 - E * S01 ) * G_b * S10 +
    !            G_r * (H10 - E * S10 ) * G_b * S01   }
    call zgemm('N','N',no,no,no,z_1,gsL  ,no,alpha,no,z_0,w    ,no)
    call zgemm('N','N',no,no,no,z_1,w    ,no,gb   ,no,z_0,rh(i),no)
    call zgemm('C','N',no,no,no,z_1,beta ,no,gb   ,no,z_0,w    ,no)
    call zgemm('N','N',no,no,no,z_1,gsR  ,no,w    ,no,z_0,rh(j),no)
    call zgemm('N','N',no,no,no,z_1,gb   ,no,s00  ,no,z_0,w    ,no)
    call zgemm('N','C',no,no,no,z_1,rh(i),no,s01  ,no,z_1,w    ,no)
    call zgemm('N','N',no,no,no,z_1,rh(j),no,s01  ,no,z_1,w    ,no)

    zDOS = 0.0_dp
    do j = 0 , nom1
       zDOS = zDOS + w(1+j*(no+1))
    end do
    ! Normalize DOS
    zDOS = zDOS / real(no,dp)

    if ( present(final_invert) ) then
       ! If we do not need to invert it, return the value
       if ( .not. final_invert ) then
          GS(:) = rh1(1:nosq)
       end if
    end if

    call timer('ts_GS',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS surface_Green' )
#endif

  contains

    ! We supply a routine to switch the pointer position of alpha,beta / rh1
    subroutine switch_alpha_beta_rh1(as_first)
      logical, intent(inout) :: as_first
      integer :: i 
      ! start
      i = 2 * nosq

      if ( as_first ) then
         rh1 => zwork(i+1:i+2*nosq) 
         i = i + 2*nosq
         alpha => zwork(i+1:i+nosq) 
         i = i + nosq
         beta => zwork(i+1:i+nosq) 
      else
         alpha => zwork(i+1:i+nosq) 
         i = i + nosq
         beta => zwork(i+1:i+nosq) 
         i = i + nosq
         rh1 => zwork(i+1:i+2*nosq) 
      end if
      as_first = .not. as_first

    end subroutine switch_alpha_beta_rh1

  end subroutine SSR_sGreen_DOS

  ! Calculates the surface Green's function for the electrodes
  ! Handles both the left and right one
  ! this is the Sancho, Sancho and Rubio algorithm
  subroutine SSR_sGreen_NoDOS(no,ZE,H00,S00,H01,S01,GS, &
       nwork, zwork, &
       iterations, final_invert)
       
! ***************** INPUT **********************************************
! integer     no      : Number of orbitals in the electrode
! complex(dp) ZE      : The energy of the Green's function evaluation
! complex(dp) H00     : Hamiltonian within the first unit cell (discarding T-direction)
! complex(dp) S00     : Overlap matrix within the first unit cell (discarding T-direction)
! complex(dp) H01     : Transfer matrix from H00 to the neighbouring cell (in T-direction)
! complex(dp) S01     : Transfer matrix from S00 to the neighbouring cell (in T-direction)
! ***************** OUTPUT *********************************************
! complex(dp) GS      : Surface Green's function of the electrode
! **********************************************************************
    use m_mat_invert
    use parallel,  only: IONode
    use precision, only: dp

! ***********************
! * INPUT variables     *
! ***********************
    integer,     intent(in) :: no
    complex(dp), intent(in) :: ZE 
    complex(dp), intent(in) :: H00(no*no),S00(no*no)
    complex(dp), intent(in) :: H01(no*no),S01(no*no)

    integer,     intent(in) :: nwork

    logical, intent(in), optional :: final_invert

! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), target :: GS(no*no)
    complex(dp), target :: zwork(nwork)

    integer, intent(out), optional :: iterations

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: nom1, no2, nosq
    integer :: ierr             !error in inversion
    integer :: i,j,ic,ic2
    logical :: as_first

    real(dp) :: ro

    ! on the stack...
    integer :: ipvt(no)
    complex(dp), dimension(:), pointer :: rh,rh1,w,alpha,beta,gb

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE surface_Green' )
#endif

    ! Initialize counter
    if ( present(iterations) ) iterations = 0

    call timer('ts_GS',1)

    nom1 = no - 1
    no2  = 2 * no
    nosq = no * no

    if ( nwork < 8 * nosq ) call die('surface_Green: &
         &Not enough work space')
    i = 0
    rh  => zwork(i+1:i+2*nosq) 
    i = i + 2*nosq
    rh1 => zwork(i+1:i+2*nosq) 
    i = i + 2*nosq
    alpha => zwork(i+1:i+nosq) 
    i = i + nosq
    beta => zwork(i+1:i+nosq) 
    i = i + nosq
    w => zwork(i+1:i+nosq)
    i = i + nosq
    gb => zwork(i+1:i+nosq) 

! gb    =   Z*S00-H00
! alpha = -(Z*S01-H01)
! gs  = Z*S00-H00
    do i = 1 , nosq
       gb(i)    = ZE * S00(i) - H00(i)
       GS(i)    = gb(i)
       alpha(i) = H01(i) - ZE * S01(i)
    end do

! beta = -(Z*S10-H10)
    do j = 1 , no
       ic = no * (j-1) + 1
       do i = 0 , nom1
          ic2 = j + no*i
          beta(ic+i) = dconjg(H01(ic2)) - ZE * dconjg(S01(ic2))
       end do
    end do

    ! Initialize loop
    ro = accur + 1._dp
    as_first = .false.
    do while ( ro > accur ) 

       ! Increment iterations
       if ( present(iterations) ) &
            iterations = iterations + 1

! rh = -(Z*S01-H01) ,j<no
! rh = -(Z*S10-H10) ,j>no
       do i = 1, nosq
          rh(i)       = alpha(i)
          rh(nosq+i)  = beta(i)
       end do

! w = Z*S00-H00
       w(:) = gb(:)

! rh =  rh1^(-1)*rh
! rh =  t0
       call zgesv(no, no2, w, no, ipvt, rh, no, ierr)

       if ( ierr /= 0 ) then
          write(*,*) 'ERROR: calc_green 1 MATRIX INVERSION FAILED'
          write(*,*) 'ERROR: LAPACK INFO = ',ierr
       end if

       ! switch pointers instead of copying elements
       call switch_alpha_beta_rh1(as_first)

! alpha = -(Z*S01-H01)*t0
       call zgemm('N','N',no,no,no,z_1,rh1(1),no,rh(1),no,z_0,alpha,no)
! beta  = -(Z*S10-H10)*t0 ??
       call zgemm('N','N',no,no,no,z_1,rh1(nosq+1),no,rh(nosq+1),no,z_0,beta,no)

! gb = gb + [ba    = (Z*S10-H10)*t0b]
       call zgemm('N','N',no,no,no,z_m1,rh1(nosq+1),no,rh(1),no,z_1,gb,no)

! ab    = (Z*S01-H01)*t0
       call zgemm('N','N',no,no,no,z_m1,rh1(1),no,rh(nosq+1),no,z_0,w,no)

       ro = -1._dp
       do i = 1 , nosq
          gb(i) = gb(i) + w(i)
          gs(i) = gs(i) + w(i)

          ! also do the accuracy calculation
          ro = max(ro,abs(w(i)))
       end do

    end do

    if ( present(final_invert) ) then
       if ( final_invert ) then
          ! Invert to get the Surface Green's function
          call mat_invert(GS,w,no,MI_IN_PLACE_LAPACK, ierr=ierr)
       end if
    else
       call mat_invert(GS,w,no,MI_IN_PLACE_LAPACK, ierr=ierr)
    end if

    if ( ierr /= 0 ) then
       write(*,*) 'ERROR: calc_green GS MATRIX INVERSION FAILED'
       write(*,*) 'ERROR: LAPACK INFO = ',ierr
    end if

    call timer('ts_GS',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS surface_Green' )
#endif

  contains

    ! We supply a routine to switch the pointer position of alpha,beta / rh1
    subroutine switch_alpha_beta_rh1(as_first)
      logical, intent(inout) :: as_first
      integer :: i 
      ! start
      i = 2 * nosq

      if ( as_first ) then
         rh1 => zwork(i+1:i+2*nosq) 
         i = i + 2*nosq
         alpha => zwork(i+1:i+nosq) 
         i = i + nosq
         beta => zwork(i+1:i+nosq) 
      else
         alpha => zwork(i+1:i+nosq) 
         i = i + nosq
         beta => zwork(i+1:i+nosq) 
         i = i + nosq
         rh1 => zwork(i+1:i+2*nosq) 
      end if
      as_first = .not. as_first

    end subroutine switch_alpha_beta_rh1

  end subroutine SSR_sGreen_NoDOS

!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------



! ##################################################################
! ## Driver subroutine for calculating the (ideal)                ##
! ## Handles both Left and Right surface Greens function.         ## 
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##                 Updated by : Nick Papior Andersen            ##
! ## It has now been parallelized to speed up electrode           ##
! ## surface Green's function generation.                         ##
! ## It generates the surface Green's function by handling        ##
! ## repetition as well.                                          ##
! ##################################################################

  subroutine create_Green(El, &
       ucell,nkpnt,kpoint,kweight, &
       RemUCellDistance, &
       NEn,ce, &
       CalcDOS,ZBulkDOS)

    use precision,  only : dp
    use parallel  , only : Node, Nodes, IONode
    use units,      only : eV
    use sys ,       only : die
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World
    use mpi_siesta, only : MPI_Bcast,MPI_ISend,MPI_IRecv
    use mpi_siesta, only : MPI_Sum, MPI_Max, MPI_integer
    use mpi_siesta, only : MPI_Wait,MPI_Status_Size
    use mpi_siesta, only : MPI_double_complex
    use mpi_siesta, only : MPI_double_precision
#endif
    use m_ts_electype
    use m_mat_invert

    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

    use m_iterator

! ***********************
! * INPUT variables     *
! ***********************
    type(Elec), intent(inout)     :: El  ! The electrode 
    integer, intent(in)           :: nkpnt ! Number of k-points
    real(dp),intent(in)           :: kpoint(3,nkpnt) ! k-points
    real(dp),intent(in)           :: kweight(nkpnt) ! weights of kpoints
    logical, intent(in)           :: RemUCellDistance ! Whether to remove the unit cell distance in the Hamiltonian.
    real(dp), dimension(3,3)      :: ucell ! The unit cell of the CONTACT
    integer, intent(in)           :: NEn ! Number of energy points
    complex(dp), intent(in)       :: ce(NEn) ! the energy points
    logical, intent(in)           :: CalcDOS ! whether or not to calculate the bulk-density of states.

! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), intent(out)      :: ZBulkDOS(NEn,Spin(El)) 

! ***********************
! * LOCAL variables     *
! ***********************
    ! Array for holding converted k-points
    real(dp), allocatable :: kE(:,:)
    real(dp) :: kpt(3), qpt(3), ktmp(3), wq
    
    ! Dimensions
    integer :: nq, nspin
    integer :: nuo_E, nS, nuou_E, nuS

    ! Electrode transfer and hamiltonian matrix
    complex(dp), pointer :: H00(:) => null()
    complex(dp), pointer :: S00(:) => null()
    complex(dp), pointer :: H01(:) => null()
    complex(dp), pointer :: S01(:) => null()
    complex(dp), pointer :: zwork(:) => null()
    complex(dp), pointer :: zHS(:) => null()

    ! Green's function variables
    complex(dp), pointer :: GS(:)
    complex(dp), pointer :: Hq(:), Sq(:), Gq(:)
    complex(dp) :: ZEnergy, zDOS

    ! In order to print information about the recursize algorithm
    integer, allocatable :: iters(:,:)
    real(dp) :: i_mean, i_std

    integer :: ierror,uGF
    ! Big loop counters
    type(itt2) :: it2
    integer, pointer :: ispin, ikpt
    integer :: iEn, iqpt
    ! Counters
    integer :: i, j, k, io, jo, off

    logical :: is_left, Gq_allocated, final_invert

#ifdef MPI
    integer :: MPIerror, curNode
    integer :: req, status(MPI_Status_Size)
    integer, allocatable :: reqs(:)
#endif
    
#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE create_Green' )
#endif

    call timer('TS_SE',1)

    ! Check input for what to do
    if( El%inf_dir == INF_NEGATIVE ) then
       is_left = .true.
    else if( El%inf_dir == INF_POSITIVE ) then
       is_left = .false.
    else
       call die("init electrode has received wrong job ID [L,R].")
    endif

    ! capture information from the electrode
    nspin  = Spin(El)
    nuo_E  = Orbs(El)
    nS     = nuo_E ** 2
    nuou_E = UsedOrbs(El)
    nuS    = nuou_E ** 2
    ! create expansion q-points (weight of q-points)
    nq     = Rep(El)
    wq     = 1._dp / real(nq,dp)
    final_invert = nq /= 1

    if (IONode) then
       write(*,'(/,2a,/,2a)') &
            "Creating Green's function file for: ",trim(name(El)), &
            "Green's function file title: ",trim(GFTitle(El))

       write(*,*) "Electrodes with transport k-points &
            & (Bohr**-1) and weights:"
       do i = 1 , nkpnt
          ! From CONTACT to electrode k-point
          ! First convert to units of reciprocal vectors
          ! Then convert to 1/Bohr in the electrode unit cell coordinates
          call kpoint_convert(ucell,kpoint(:,i),ktmp,1)
          if ( RepA1(El) > 1 ) ktmp(1) = ktmp(1)/real(RepA1(El),dp)
          if ( RepA2(El) > 1 ) ktmp(2) = ktmp(2)/real(RepA2(El),dp)
          if ( RepA3(El) > 1 ) ktmp(3) = ktmp(3)/real(RepA3(El),dp)
          call kpoint_convert(UnitCell(El),ktmp,kpt,-1)
          write(*,'(i4,2x,4(E14.5))') i, kpt,kweight(i)
       end do

       ! Show the number of used atoms and orbitals
       write(*,'(a,i6,'' / '',i6)') ' Atoms available    / used atoms   : ', &
            Atoms(El),UsedAtoms(El)
       write(*,'(a,i6,'' / '',i6)') ' Orbitals available / used orbitals: ', &
            Orbs(El),UsedOrbs(El)


       ! We show them in units of Bohr**-1
       write(*,'(a)') ' q-points for expanding electrode (Bohr**-1):'
       do i = 1 , nq
          call kpoint_convert(unitcell(El),q_exp(El,i),qpt,-1)
          write(*,'(i4,2x,4(E14.5))') i,qpt,wq
       end do
       write(*,'(a,f14.5,1x,a)') &
            " Fermi level shift in electrode (chemical potential) : ",El%mu%mu/eV,' eV'
    end if

    ! Initialize Green's function and Hamiltonian arrays
    nullify(GS)
    if ( nS /= nuS ) then
       allocate(GS(nS))
       call memory('A','Z',nS,'create_green')
    !else
    !  the regions are of same size, so we can just point
    !  to the correct memory segment
    end if

    ! Allocate work array
    allocate(zwork(max(nS*9,nuS*nq*2)))
    call memory('A','Z',size(zwork),'create_Green')

    ! Point the hamiltonian and the overlap to the work array
    ! The work-array is only used for calculation the surface
    ! Green's function and
    Hq => zwork(1:nuS*nq)
    Sq => zwork(nuS*nq+1:nuS*nq*2)
    if ( size(zwork) >= nS * 9 + nuS*nq ) then
       Gq => zwork(nS*9+1:nS*9+nuS*nq)
       Gq_allocated = .false.
    else
       nullify(Gq)
       allocate(Gq(nuS*nq))
       call memory('A','Z',size(Gq),'create_green')
       Gq_allocated = .true.
    end if

    ! all the Hamiltonian and overlaps
    allocate(zHS(nS * nq * 4))
    call memory('A','Z',nS * nq * 4,'create_Green')

    ! Prepare for the inversion
    call init_mat_inversion(nuo_E)

    ! Reset bulk DOS
    if ( CalcDOS ) then
       ZBulkDOS(:,:) = dcmplx(0._dp,0._dp)
    end if

!******************************************************************
!           Start Green's function calculation
!******************************************************************
    
    if (IONode) then
       call io_assign(uGF)
       open(FILE=GFfile(El),UNIT=uGF,FORM='UNFORMATTED')

       ! Initial header for file
       write(uGF) GFTitle(El)
       ! Electrode information
       write(uGF) Spin(El), UnitCell(El)
       write(uGF) UsedAtoms(El),UsedOrbs(El)
       write(uGF) El%xa_used, El%lasto_used
       write(uGF) RepA1(El),RepA2(El),RepA3(El)
       write(uGF) El%mu%mu

       ! Write out explicit information about this content
       write(uGF) RemUCellDistance
       write(uGF) nkpnt
       ! Notice that we write the k-points for the ELECTRODE
       ! Do a conversion here
       allocate(kE(3,nkpnt))
       call memory('A','D',nkpnt*3,'create_green')
       do i = 1 , nkpnt
          ! Init kpoint, in reciprocal vector units ( from CONTACT ucell)
          call kpoint_convert(ucell,kpoint(:,i),ktmp,1)
          ktmp(1) = ktmp(1)/real(RepA1(El),dp)
          ktmp(2) = ktmp(2)/real(RepA2(El),dp)
          ktmp(3) = ktmp(3)/real(RepA3(El),dp)
          ! Convert back to reciprocal units (to electrode ucell_E)
          call kpoint_convert(UnitCell(El),ktmp,kE(:,i),-1)
       end do
       write(uGF) kE,kweight
       call memory('D','D',nkpnt*3,'create_green')
       deallocate(kE)

       ! write out the contour information
       write(uGF) NEn
       write(uGF) ce ! energy points
    end if
    
#ifdef MPI
    if ( IONode ) then
       allocate(reqs(Nodes-1))
       call memory('A','I',Nodes-1,'create_green')
       ! Create request handles for communication
       ! This is a rather new feature which enhances communication times.
       ! However, this is perhaps overkill as we never have VERY many 
       ! contour points. Say NEn > 1000
       ! Look in the loop for MPI_Start(...) for where this is used
       do i = 1 , Nodes - 1
          call MPI_Recv_Init(Gq(1),nuS*nq,MPI_double_complex, &
               i,i,MPI_Comm_World,reqs(i),MPIerror)
       end do
    else
       ! Create request handles for communication
       call MPI_Send_Init(Gq(1),nuS*nq,MPI_double_complex, &
            0,Node,MPI_Comm_World,req,MPIerror)
    end if
#endif

    ! prepare the iteration counter
    allocate(iters(NEn,2))
    if ( IONode ) then
       ! TODO when adding new surface-Green's functions schemes, please update here
       write(*,'(1x,a)') 'Lopez Sancho, Lopez Sancho & Rubio recursive &
            &surface self-energy calculation...'
    end if

    ! start up the iterators
    call itt_init  (it2,end1=nspin,end2=nkpnt)
    call itt_attach(it2,cur1=ispin,cur2=ikpt)

    ! do spin and k-point loop in one go...
    do while ( .not. itt_step(it2) )
       
       if ( itt_stepped(it2,1) ) then
          ! Number of iterations
          iters(:,:) = 0
       end if
       
       ! Init kpoint, in reciprocal vector units ( from CONTACT ucell)
       call kpoint_convert(ucell,kpoint(:,ikpt),ktmp,1)
       ktmp(1) = ktmp(1)/real(RepA1(El),dp)
       ktmp(2) = ktmp(2)/real(RepA2(El),dp)
       ktmp(3) = ktmp(3)/real(RepA3(El),dp)
       ! Convert back to reciprocal units (to electrode)
       call kpoint_convert(UnitCell(El),ktmp,kpt,-1)
       
       ! loop over the repeated cell...
       HSq_loop: do iqpt = 1 , nq
             
          ! point to the correct segment of memory
          H00 => zHS((     iqpt-1)*nS+1:      iqpt *nS)
          S00 => zHS((  nq+iqpt-1)*nS+1:(  nq+iqpt)*nS)
          H01 => zHS((2*nq+iqpt-1)*nS+1:(2*nq+iqpt)*nS)
          S01 => zHS((3*nq+iqpt-1)*nS+1:(3*nq+iqpt)*nS)

          ! init qpoint in reciprocal lattice vectors
          call kpoint_convert(UnitCell(El),q_exp(El,iqpt),qpt,-1)

          ! Setup the transfer matrix and the intra cell at the k-point and q-point
          if ( RemUCellDistance ) then
             call die('Not working yet')
          end if

          ! Calculate transfer matrices @Ef (including the chemical potential)
          call set_electrode_HS_Transfer(ispin, El, kpt, qpt, &
               nS, H00,S00,H01,S01, RemUCellDistance=RemUCellDistance)

          i = (iqpt-1)*nuS
          if ( nuo_E /= nuou_E ) then
             if( is_left ) then
                ! Left, we use the last orbitals
                off = nuo_E - nuou_E + 1
                do jo = off - 1 , nuo_E - 1
                   do io = off , nuo_E
                      i = i + 1
                      Sq(i) = S00(io+nuo_E*jo)
                      Hq(i) = H00(io+nuo_E*jo)
                   end do
                end do
             else
                ! Right, the first orbitals
                do jo = 0 , nuou_E - 1
                   do io = 1 , nuou_E
                      i = i + 1
                      Sq(i) = S00(io+nuo_E*jo)
                      Hq(i) = H00(io+nuo_E*jo)
                   end do   ! io
                end do      ! jo
             end if
          end if
          
       end do HSq_loop
       
       if ( IONode ) then
          write(uGF) ikpt, 1, ce(1) ! k-point and energy point
          if ( nuo_E /= nuou_E ) then
             write(uGF) Hq
             write(uGF) Sq
          else
             H00 => zHS(      1:nq*nS  )
             S00 => zHS(nq*nS+1:nq*nS*2)
             write(uGF) H00
             write(uGF) S00
          end if
       end if
       
       Econtour_loop: do iEn = 1, NEn
          
#ifdef MPI
          ! Every node takes one energy point
          ! This asserts that IONode = Node == 0 will have iEn == 1
          ! Important !
          curNode = MOD(iEn-1,Nodes)
          E_Nodes: if ( curNode == Node ) then
#endif
             ! as we already have shifted H,S to Ef + mu, and ZEnergy is
             ! wrt. mu, we don't need to subtract mu again
             ZEnergy = ce(iEn)
             
! loop over the repeated cell...
             q_loop: do iqpt = 1 , nq

                H00 => zHS((     iqpt-1)*nS+1:      iqpt *nS)
                S00 => zHS((  nq+iqpt-1)*nS+1:(  nq+iqpt)*nS)
                H01 => zHS((2*nq+iqpt-1)*nS+1:(2*nq+iqpt)*nS)
                S01 => zHS((3*nq+iqpt-1)*nS+1:(3*nq+iqpt)*nS)
                if ( nS == nuS ) then
                   ! instead of doing a copy afterward, we can
                   ! put it the correct place immediately
                   GS => Gq((    iqpt-1)*nS+1:      iqpt *nS)
                end if

                ! Calculate the surface Green's function
                ! Zenergy is wrt. to the system Fermi-level
                if ( CalcDOS ) then
                   call SSR_sGreen_DOS(nuo_E,ZEnergy,H00,S00,H01,S01,GS, &
                        zDOS,9*nS,zwork, &
                        iterations=iters(iEn,1), final_invert = final_invert)
                   
                   ! We also average the k-points.
                   ZBulkDOS(iEn,ispin) = ZBulkDOS(iEn,ispin) + &
                        wq * zDOS * kweight(ikpt)

                else
                   call SSR_sGreen_NoDos(nuo_E,ZEnergy,H00,S00,H01,S01,GS, &
                        8*nS,zwork, &
                        iterations=iters(iEn,1), final_invert = final_invert)
                   
                end if

                  
                ! Copy over surface Green's function
                i = (iqpt-1)*nuS
                if ( nS /= nuS ) then
                   if ( is_left ) then
                      ! Left, we use the last orbitals
                      do jo = nuo_E - nuou_E , nuo_E - 1
                         do io = nuo_E - nuou_E + 1 , nuo_E
                            i = i + 1
                            Gq(i) = GS(io+nuo_E*jo)
                         end do           ! io
                      end do              ! jo
                   else
                      ! Right, the first orbitals
                      do jo = 0,nuou_E-1
                         do io = 1,nuou_E
                            i = i + 1
                            Gq(i) = GS(io+nuo_E*jo)
                         end do           ! io
                      end do              ! jo
                   end if
                end if
                
             end do q_loop
             
             if (IONode) then
                ! Write out calculated information at E point

                if ( iEn /= 1 ) write(uGF) ikpt, iEn, ce(iEn)
                write(uGF) Gq

             end if

#ifdef MPI
             ! If not IONode we should send message
             ! This message parsing is directly connected to 
             ! a predefined size of the message, see right before
             ! spin loop.
             ! It communicates the Gq array to the Gq array
             if ( .not. IONode ) then
                call MPI_Start(req,MPIerror)
                call MPI_Wait(req,status,MPIerror)
             end if
             
          end if E_Nodes

          ! If IONode, we should receive in each energy point
          ! There is no need to create a buffer array for the Gq
          ! We will not use it until we are in the loop again
          if ( IONode .and. curNode /= Node ) then
             call MPI_Start(reqs(curNode),MPIerror)
             if ( iEn /= 1 ) write(uGF) ikpt, iEn, ce(iEn)
             call MPI_Wait(reqs(curNode),status,MPIerror)
             write(uGF) Gq
          end if

#endif

       end do Econtour_loop
          

       if ( itt_last(it2,2) ) then
#ifdef MPI
          call MPI_Reduce(iters(1,1), iters(1,2), NEn, MPI_Integer, MPI_Sum, &
               0, MPI_Comm_World, MPIerror)
#else
          iters(:,2) = iters(:,1)
#endif
          if ( IONode ) then
             i_mean = sum(iters(:,2)) / real(NEn,dp)
             i_std = 0._dp
             do i = 1 , NEn
                i_std = i_std + ( iters(i,2) - i_mean ) ** 2
             end do
             i_std = sqrt(i_std/real(NEn,dp))
             ! TODO if new surface-Green's function scheme is implemented, fix here
             write(*,'(1x,a,f10.4,'' / '',f10.4)') 'Lopez Sancho, Lopez Sancho & Rubio: &
                  &Mean/std iterations: ', i_mean             , i_std
             write(*,'(1x,a,i10,'' / '',i10)')     'Lopez Sancho, Lopez Sancho & Rubio: &
                  &Min/Max iterations : ', minval(iters(:,2)) , maxval(iters(:,2))
          end if
       end if

    end do
!*******************************************************************
!         Green's function calculation is done
!*******************************************************************

    deallocate(iters)

#ifdef MPI
    ! Free requests made for the communications
    if ( IONode ) then
       do i = 1 , Nodes - 1 
          call MPI_Request_Free(reqs(i),MPIerror)
       end do
       call memory('D','I',Nodes-1,'create_green')
       deallocate(reqs)
    else
       call MPI_Request_Free(req,MPIerror)
    end if
#endif

    ! Close file
    if ( IONode ) then
       call io_close(uGF)
       write(*,'(a)') "Done creating '"//trim(GFFile(El))//"'."  
    end if
    
    ! Clean up computational arrays
    if ( nS /= nuS ) then
       call memory('D','Z',size(GS),'create_green')
       deallocate(GS)
    end if

    if ( Gq_allocated ) then
       call memory('D','Z',size(Gq),'create_green')
       deallocate(Gq)
    end if

    ! Work-arrays
    call memory('D','Z',size(zwork),'create_green')
    deallocate(zwork)
    call memory('D','Z',size(zHS),'create_green')
    deallocate(zHS)

    call itt_destroy(it2)

    call clear_mat_inversion()

#ifdef MPI
    if ( CalcDOS ) then
       ! Sum the bulkdensity of states
       ! Here we can safely use the array as temporary (Gq)
       allocate(Gq(NEn*nspin))
       call memory('A','Z',NEn*nspin,'create_green')
       Gq = 0.0_dp
       call MPI_AllReduce(ZBulkDOS(1,1),Gq(1),NEn*nspin, MPI_double_complex, &
            MPI_Sum,MPI_Comm_World,MPIerror)
       i = 0
       do j = 1 , nspin
          do io = 1 , NEn
             i = i + 1
             ZBulkDOS(io,j) = Gq(i)
          end do
       end do
       call memory('D','Z',NEn*nspin,'create_green')
       deallocate(Gq)
    end if
#endif

    call timer('TS_SE',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS create_Green' )
#endif
    
  end subroutine create_Green

  subroutine init_Electrode_HS(El,RemUCellDistance)
    use m_ts_electype
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

    type(Elec), intent(inout) :: El
    logical, intent(in) :: RemUCellDistance
    
    ! Read-in and create the corresponding transfer-matrices
    call delete(El) ! ensure clean electrode
    call read_Elec(El,Bcast=.true.)

    if ( .not. initialized(El%xij) ) then
       call die('An electrode file needs to be a non-Gamma calculation. &
            &Ensure at least two k-points in the T-direction.')
    end if

    ! print out the precision of the electrode (whether it extends
    ! beyond first principal layer)
    call check_Connectivity(El)
    
    call create_sp2sp01(El,calc_xijo=Rep(El)/=1 .or. RemUCellDistance)
    ! Clean-up, we will not need these!
    ! we should not be very memory hungry now, but just in case...
    call delete(El%H)
    call delete(El%S)
   
    ! We do not accept onlyS files
    if ( .not. initialized(El%H00) ) then
       call die('An electrode file must contain the Hamiltonian')
    end if

    call delete(El%xij)
    call delete(El%sp)
    if ( RemUCellDistance ) then
       call delete(El%xij00)
       call delete(El%xij01)
    end if

  end subroutine init_Electrode_HS


!**********
! Create the Hamiltonian for the electrode as well
! as creating the transfer matrix.
!**********
  subroutine set_electrode_HS_Transfer(ispin,El,k,q,nS,Hk,Sk,Hk_T,Sk_T,RemUCellDistance)
    use sys, only : die
    use precision, only : dp
    use m_ts_electype
    use geom_helper, only : ucorb
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

! ***********************
! * INPUT variables     *
! ***********************
    integer, intent(in)    :: ispin, nS
    type(Elec), intent(inout) :: El
    real(dp), intent(in)   :: k(3)   ! k-point in [1/Bohr]
    real(dp), intent(in)   :: q(3)   ! expansion k-point in [1/Bohr]
! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), dimension(nS) :: Hk,Sk,Hk_T,Sk_T

    logical, intent(in), optional :: RemUCellDistance

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: no_u
    real(dp) :: kqxij, Ef
    complex(dp) :: cphase
    integer :: i,j,iuo,juo,ind
    integer, pointer :: ncol00(:), l_ptr00(:),l_col00(:)
    integer, pointer :: ncol01(:), l_ptr01(:),l_col01(:)
    real(dp), pointer :: xij00(:,:), xij01(:,:), xijo00(:,:), xijo01(:,:) 
    real(dp), pointer :: H00(:,:) , S00(:), H01(:,:), S01(:)
    integer :: t_dir
    logical :: has_o

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE elec_HS_Transfer' )
#endif

    t_dir = El%t_dir
    ! we need to subtract as the below code shifts to Ef
    Ef    = El%Ef - El%mu%mu
    no_u  = Orbs(El) ! this has to be the total size of the electrode

    if ( no_u ** 2 /= nS ) call die('Wrong size of the electrode array')

    ! retrieve values
    call attach(El%sp00,n_col=ncol00,list_ptr=l_ptr00,list_col=l_col00)
    call attach(El%sp01,n_col=ncol01,list_ptr=l_ptr01,list_col=l_col01)
    ! point to the data-segments...
    H00 => val(El%H00)
    H01 => val(El%H01)
    S00 => val(El%S00)
    S01 => val(El%S01)
    if ( RemUCellDistance ) then
       xij00 => val(El%xijo00)
       xij01 => val(El%xijo01)
    else
       xij00 => val(El%xij00)
       xij01 => val(El%xij01)
    end if
    xijo00 => val(El%xijo00)
    xijo01 => val(El%xijo01)
    has_o = initialized(El%xijo00)


    ! Initialize arrays
    do i = 1, nS
       Hk(i)   = dcmplx(0._dp,0._dp)
       Sk(i)   = dcmplx(0._dp,0._dp)
       Hk_T(i) = dcmplx(0._dp,0._dp)
       Sk_T(i) = dcmplx(0._dp,0._dp)
    enddo

    do iuo = 1 , no_u

       ! Create the 00
       do j = 1 , ncol00(iuo)
          ind = l_ptr00(iuo) + j
          juo = Ucorb(l_col00(ind),no_u)
          kqxij = 0._dp
          if ( has_o ) then
             do i = 1 , 3 
                if ( i == t_dir ) cycle
                kqxij = kqxij + k(i) * xij00(i,ind) + q(i) * xijo00(i,ind)
             end do
          else
             do i = 1 , 3 
                if ( i == t_dir ) cycle
                kqxij = kqxij + k(i) * xij00(i,ind)
             end do
          end if

          cphase = cdexp(dcmplx(0._dp,kqxij) )
          
          i = iuo+(juo-1)*no_u
          Hk(i) = Hk(i) + H00(ind,ispin) * cphase
          Sk(i) = Sk(i) + S00(ind)       * cphase
       enddo

       ! Create the 01
       do j = 1 , ncol01(iuo)
          ind = l_ptr01(iuo) + j
          juo = ucorb(l_col01(ind),no_u)
          kqxij = 0._dp
          if ( has_o ) then
             do i = 1 , 3 
                if ( i == t_dir ) cycle
                kqxij = kqxij + k(i) * xij01(i,ind) + q(i) * xijo01(i,ind)
             end do
          else
             do i = 1 , 3 
                if ( i == t_dir ) cycle
                kqxij = kqxij + k(i) * xij00(i,ind)
             end do
          end if

          cphase = cdexp(dcmplx(0._dp,kqxij) )
          
          i = iuo+(juo-1)*no_u
          Hk_T(i) = Hk_T(i) + H01(ind,ispin) * cphase
          Sk_T(i) = Sk_T(i) + S01(ind)       * cphase
       end do
    end do

    ! Symmetrize and make EF the energy-zero
    do iuo = 1,no_u
       do juo = 1,iuo-1
          i = iuo+(juo-1)*no_u
          j = juo+(iuo-1)*no_u

          Sk(j) = 0.5_dp*( Sk(j) + dconjg(Sk(i)) )
          Sk(i) = dconjg(Sk(j))

          Hk(j) = 0.5_dp*( Hk(j) + dconjg(Hk(i)) ) - Ef * Sk(j)
          Hk(i) = dconjg(Hk(j))

          ! Transfer matrix is not symmetric
          Hk_T(i) = Hk_T(i) - Ef * Sk_T(i)
          Hk_T(j) = Hk_T(j) - Ef * Sk_T(j)

       end do
       
       i = iuo+(iuo-1)*no_u
       Sk(i) = Sk(i) - dcmplx(0._dp,dimag(Sk(i)))
       
       Hk(i) = Hk(i) - dcmplx(0._dp,dimag(Hk(i))) - Ef * Sk(i)
       
       ! Transfer matrix
       Hk_T(i) = Hk_T(i) - Ef * Sk_T(i)
    end do

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS elec_HS_Transfer' )
#endif

  end subroutine set_electrode_HS_Transfer

end module m_ts_electrode
