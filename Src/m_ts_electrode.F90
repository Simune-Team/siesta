! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
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
! The TranSIESTA routine can at the moment not create the surface Green's
! function without having a transfer matrix.
! Thus we ENFORCE Gamma == .false. and the program should die if
! the electrode calculation was a Gamma calculation!

  implicit none

  public :: create_Green

  private !


  ! This is the integer telling which direction is the propagation direction
  ! DO NOT CHANGE THIS! 
  ! At the moment all code is enforced to use PropDir == 3
  ! TODO move this to an option in the m_ts_options module.
  ! It should be user specific (however, it requires a lot of
  ! recoding)
  integer, parameter :: PropDir = 3

contains

  ! Calculates the surface Green's function for the electrodes
  ! Handles both the left and right one
  subroutine surface_Green(tjob,nv,Zenergy,h00,s00,h01,s01, &
       gs,zdos)
! ***************** INPUT **********************************************
! character   tjob    : Specifies the left or the right electrode
! integer     nv      : Number of orbitals in the electrode
! complex(dp) Zenergy : The energy of the Green's function evaluation
! complex(dp) H00     : Hamiltonian within the first unit cell
! complex(dp) S00     : Overlap matrix within the first unit cell
! complex(dp) H01     : Transfer matrix from H00 to the neighbouring cell
! complex(dp) S01     : Transfer matrix from S00 to the neighbouring cell
! ***************** OUTPUT *********************************************
! complex(dp) gs      : Surface Green's function of the electrode
! complex(dp) zdos    : Density of energy point
! **********************************************************************
    use precision, only : dp
    use units, only : Pi
    use fdf, only : leqi
    use m_ts_aux_rout, only : csolveg

! ***********************
! * INPUT variables     *
! ***********************
    character(len=1) :: tjob
    integer :: nv
    complex(dp) :: ZEnergy 
    complex(dp) :: h00(0:nv*nv-1),s00(0:nv*nv-1)
    complex(dp) :: h01(0:nv*nv-1),s01(0:nv*nv-1)
! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp) :: gs(0:nv*nv-1)
    complex(dp) :: zdos

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: nv2,iter
    integer :: ierr             !error in inversion
    integer :: i,j,ic,ic2

    complex(dp) :: a,b
    real(dp)    :: ro
    real(dp), parameter :: accur=1.d-15

    integer, dimension(:), allocatable :: ipvt
    complex(dp), dimension(:), allocatable :: &
         rh,rh1,rh3,alpha,beta,ab,ba,gb,gs2

#ifdef TRANSIESTA_DEBUG
      call write_debug( 'PRE surface_Green' )
#endif


    allocate(ipvt(nv))
    allocate(rh(0:2*nv*nv))
    allocate(rh1(0:2*nv*nv))
    allocate(rh3(0:4*nv*nv))
    allocate(alpha(0:nv*nv-1))
    allocate(beta(0:nv*nv-1))
    allocate(ba(0:nv*nv-1))
    allocate(ab(0:nv*nv-1))
    allocate(gb(0:nv*nv-1))
    allocate(gs2(0:nv*nv-1))
    call memory('A','I',nv,'calc_green')
    call memory('A','Z',14*nv*nv+3,'calc_green')

    a=(1.d0,0.d0)
    b=(0.d0,0.d0)
    nv2 =2*nv

! FDN
! gb = Z*S00-H00
! alpha = -(Z*S01-H01)        
    do i=0,nv*nv-1
       gb(i) = zenergy*s00(i)-h00(i)
       alpha(i) = h01(i)-zenergy*s01(i)
    end do

! FDN
! gs = Z*S00-H00
! gs2 = Z*S00-H00
    do i=0,nv*nv-1
       gs(i) = gb(i)
       gs2(i) = gb(i)
    end do

! FDN
! beta = -(Z*S10-H10)
    do j=0,nv-1
       do i=0,nv-1
          ic = i + nv*j
          ic2 = j + nv*i
          beta(ic) = dconjg(h01(ic2))-zenergy*dconjg(s01(ic2))
       end do
    end do


    iter=0
1000 continue
    iter=iter+1


! FDN
! nv2=2*nv
! rh = -(Z*S01-H01) ,j<nv
! rh = -(Z*S10-H10) ,j>nv
    do j=0,nv2-1
       do i=0,nv-1
          ic =i + j*nv
          ic2=i + (j - nv)*nv
          if(j.lt.nv)then
             rh(ic) = alpha(ic)
          else
             rh(ic) = beta(ic2)
          endif
       end do
    end do

! FDN
! rh3 = Z*S00-H00
    do i=0,nv*nv-1
       rh3(i) = gb(i)
    end do

! FDN
! rh =  rh3^(-1)*rh
! rh = t0
    call csolveg(nv,nv2,rh3,rh,ipvt,ierr) 

    if(IERR.ne.0) then
       write(*,*) 'ERROR: calc_green 1 MATRIX INVERSION FAILED'
       write(*,*) 'ERROR: LAPACK INFO = ',IERR
    end if


! FDN
! nv2=2*nv
! rh1 = -(Z*S01-H01) ,j<nv
! rh1 = -(Z*S10-H10) ,j>nv
! a = 1
! b = 0
    do j=0,nv-1
       do i=0,nv2-1
          ic =i + j*nv
          ic2 =i-nv + j*nv
          if(i.lt.nv)then
             rh1(i + nv2*j) = alpha(ic)
          else
             rh1(i + nv2*j) = beta(ic2)
          end if
       end do
    end do
! FDN
! rh3 = 1*rh1*rh + 0*rh3                   
! rh3 = -(Z*S01-H01)*t0
    call zgemm('N','N',nv2,nv2,nv,a,rh1,nv2,rh,nv,b,rh3,nv2)

! FDN
! aplha = -(Z*S01-H01)*t0
! ba = -(Z*S10-H10)*t0b
    do j=0,nv-1
       do i=0,nv2-1
          ic =i + j*nv
          ic2 =i-nv + j*nv
          if(i.lt.nv)then
             alpha(ic) = rh3(i + nv2*j) 
          else
             ba(ic2) = -rh3(i + nv2*j) 
          end if
       end do
    end do
    do j=nv,nv2-1
       do i=0,nv2-1
          ic=i + (j - nv)*nv
          ic2=i - nv + (j - nv)*nv
          if(i.lt.nv)then
             ab(ic)= -rh3(i + nv2*j) 
          else
             beta(ic2)= rh3(i + nv2*j) 
          end if
       end do
    end do

    do i=0,nv*nv-1
       gb(i) =  gb(i) + ba(i) + ab(i)
       gs(i) =  gs(i) + ab(i) 
       gs2(i) =  gs2(i) + ba(i) 
    end do

    ro =-1.0d0
    do j =0,nv*nv-1
       ro =max(ro,dreal(ab(j))**2+dimag(ab(j))**2)
    end do
    ro =dsqrt(ro)

    if(ro.gt.accur) go to 1000


    do i=0,nv*nv-1
       rh3(i) = gs(i)
       rh(i) = 0.0d0
    end do

    do j=0,nv-1
       rh(j*(nv+1)) = 1.d0
    end do

    call csolveg(nv,nv,rh3,rh,ipvt,ierr)

    if(IERR.ne.0) then
       write(*,*) 'ERROR: calc_green 2 MATRIX INVERSION FAILED'
       write(*,*) 'ERROR: LAPACK INFO = ',IERR
    end if

    do i=0,nv*nv-1
       gs(i) = rh(i)
    end do



    do i=0,nv*nv-1
       rh3(i) = gs2(i)
       rh(i) = 0.0d0
    end do

    do j=0,nv-1
       rh(j*(nv+1)) = 1.d0
    end do

    call csolveg(nv,nv,rh3,rh,ipvt,ierr)


    if(IERR.ne.0) then
       write(*,*) 'ERROR: calc_green 3 MATRIX INVERSION FAILED'
       write(*,*) 'ERROR: LAPACK INFO = ',IERR
    end if

    do i=0,nv*nv-1
       gs2(i) = rh(i)
    end do


    do i=0,nv*nv-1
       rh3(i) = gb(i)
       rh(i) = 0.0d0
    end do

    do j=0,nv-1
       rh(j*(nv+1)) = 1.d0
    end do

    call csolveg(nv,nv,rh3,rh,ipvt,ierr)

    if(IERR.ne.0) then
       write(*,*) 'ERROR: calc_green 4 MATRIX INVERSION FAILED'
       write(*,*) 'ERROR: LAPACK INFO = ',IERR
    end if

    do i=0,nv*nv-1
       gb(i) = rh(i)
    end do


!      ----      DOS     -----

    do i=0,nv*nv-1
       alpha(i) = h01(i)-zenergy*s01(i)
    end do
    call zgemm('N','N',nv,nv,nv,a,gs2,nv,alpha,nv,b,ab,nv)
    do i=0,nv*nv-1
       alpha(i) = ab(i)
    end do
    call zgemm('N','N',nv,nv,nv,a,alpha,nv,gb,nv,b,ab,nv)


    do j=0,nv-1
       do i=0,nv-1
          ic = i + nv*j
          ic2 = j + nv*i
          beta(ic) = dconjg(h01(ic2))-zenergy*dconjg(s01(ic2))
       end do
    end do
    call zgemm('N','N',nv,nv,nv,a,gs,nv,beta,nv,b,ba,nv)
    do i=0,nv*nv-1
       beta(i) = ba(i)
    end do
    call zgemm('N','N',nv,nv,nv,a,beta,nv,gb,nv,b,ba,nv)


    do i=0,nv*nv-1
       rh3(i) = 0.0d0
    end do

    call zgemm('N','N',nv,nv,nv,a,gb,nv,s00,nv,b,rh3,nv)
    call zgemm('N','C',nv,nv,nv,a,ab,nv,s01,nv,a,rh3,nv)
    call zgemm('N','N',nv,nv,nv,a,ba,nv,s01,nv,a,rh3,nv)


    zdos = 0.0d0

    do j=0,nv-1
       zdos = zdos + (rh3(j*(nv+1)))
    end do



    if( leqi(tjob,'L') ) then
       do i=0,nv*nv-1
          gs(i) =  gs2(i) 
       end do
    endif

    call memory('D','Z',14*nv*nv+3,'calc_green')
    call memory('D','I',nv,'calc_green')
    deallocate(ipvt)
    deallocate(rh,rh1,rh3)
    deallocate(alpha,beta)
    deallocate(ba,ab)
    deallocate(gb,gs2)

#ifdef TRANSIESTA_DEBUG
      call write_debug( 'POS surface_Green' )
#endif

  end subroutine surface_Green

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

  subroutine create_Green(tElec, HSFile, GFFile, GFTitle, &
       ElecValenceBandBot, &
       nkpnt,kpoint,kweight, &
       NBufAt,NUsedAtoms,NA1,NA2, RemUCellDistance, &
       ucell,xa,nua,NEn,contour,chem_shift,ZBulkDOS,nspin)

    use precision,  only : dp
    use fdf,        only : leqi
    use parallel  , only : Node, Nodes, IONode
    use units,      only : eV
    use sys ,       only : die
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World
    use mpi_siesta, only : MPI_Bcast,MPI_ISend,MPI_IRecv
    use mpi_siesta, only : MPI_Sum
    use mpi_siesta, only : MPI_Wait,MPI_Status_Size
    use mpi_siesta, only : DAT_dcomplex => MPI_double_complex, &
                           DAT_double => MPI_double_precision
#endif
    use m_hs_matrix,only : set_HS_matrix, matrix_symmetrize
    use m_ts_cctype
! ***********************
! * INPUT variables     *
! ***********************
    character(len=1), intent(in) :: tElec   ! 'L' for Left electrode, 'R' for right
    character(len=*), intent(in) :: HSFile  ! The electrode TSHS file
    character(len=*), intent(in) :: GFFile  ! The electrode GF file to be saved to
    character(len=*), intent(in) :: GFTitle ! The title to be written in the GF file
    logical, intent(in)          :: ElecValenceBandBot ! Whether or not to calculate electrodes valence bandbottom
    integer, intent(in)            :: nkpnt ! Number of k-points
    real(dp),dimension(3,nkpnt),intent(in) :: kpoint ! k-points
    real(dp),dimension(nkpnt),intent(in) :: kweight ! weights of kpoints
    integer, intent(in)            :: NBufAt,NA1,NA2 ! Buffer/Rep a1/Rep a2
    logical, intent(in)            :: RemUCellDistance ! Whether to remove the unit cell distance in the Hamiltonian.
    integer, intent(in)            :: NUsedAtoms ! Needs update here
    integer, intent(in)            :: nua ! Full system count of atoms in unit cell
    real(dp), dimension(3,3)       :: ucell ! The unit cell of the CONTACT
    real(dp), intent(in)           :: xa(3,nua) ! Coordinates in the system for the TranSIESTA routine
    integer, intent(in)            :: nspin ! spin in system
    integer, intent(in)            :: NEn ! Number of energy points
    type(ts_ccontour), intent(in)  :: contour(NEn) ! contours path for GF
    real(dp), intent(in)           :: chem_shift ! the Fermi-energy we REQUIRE the electrode
! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), intent(out)       :: ZBulkDOS(NEn,nspin) 

! ***********************
! * LOCAL variables     *
! ***********************
    logical :: Gamma
    character(len=5)   :: GFjob ! Contains either 'Left' or 'Right'
    logical :: exist ! Checking for file existance
    
! >>>>>>>>>> Electrode TSHS variables
! We suffix with _E to distinguish from CONTACT
    integer                           :: nua_E,nuo_E,maxnh_E ! Unit cell atoms / orbitals / Hamiltonian size
    integer                           :: nuou_E ! # used orbitals from electrode (if NUsedAtoms < nua_E)
    integer                           :: notot_E ! Total number of orbitals
    real(dp), dimension(:,:), pointer :: H_E,xij_E,xijo_E ! Hamiltonian, differences with unitcell, differences without unitcell
    real(dp), dimension(:,:), pointer :: xa_E ! atomic coordinats
    real(dp), dimension(:),   pointer :: S_E ! Overlap
    integer,  dimension(:),   pointer :: zconnect_E ! Has 0 values for indices where there are no z-connection.
    integer,  dimension(:),   pointer :: numh_E,listhptr_E,listh_E,indxuo_E,lasto_E
    real(dp)                          :: Ef_E ! Efermi
    real(dp), dimension(3,3)          :: ucell_E ! The unit cell of the electrode

    ! Array for holding eigen values
    real(dp), dimension(:), allocatable :: eig
! <<<<<<<<<< Electrode TSHS variables

    integer :: nq ! number of q-points, set 'ts_mkqgrid'
    real(dp), dimension(:,:), pointer :: qb => null() 
                                            ! q points for repetition, in units
                                            ! of reciprocal lattice vectors (hence the b)
    real(dp), dimension(:), pointer   :: wq => null() 
                                            ! weights for q points for repetition
    real(dp) :: kpt(3), qpt(3), ktmp(3)
    
    ! Electrode transfer and hamiltonian matrix
    complex(dp), dimension(:), pointer :: H00 => null()
    complex(dp), dimension(:), pointer :: S00 => null()
    complex(dp), dimension(:), pointer :: H01 => null()
    complex(dp), dimension(:), pointer :: S01 => null()

    ! Green's function variables
    complex(dp), dimension(:), allocatable :: GS
    complex(dp), dimension(:,:), allocatable :: Hq,Sq,Gq
    complex(dp) :: ZEnergy, ZSEnergy, zdos

    integer :: ierror,uGF
    ! Big loop counters
    integer :: iEn,ispin, iqpt,ikpt
    ! Counters
    integer :: i,j,ia,ja,io,jo

#ifdef MPI
    integer :: MPIerror,curNode
    integer :: req, status(MPI_Status_Size)
    integer, allocatable :: reqs(:)
#endif
    
#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE create_Green' )
#endif

    ! Should we read Gamma in TSHS file?
    ! Must be .false., otherwise we cannot access Transfer matrix
    ! This Gamma is to be used for the remaining part of the tests
    ! We cannot use the ts_gamma in case of Gamma point in kxy direction
    Gamma = .false.

    ! Check input for what to do
    if( leqi(tElec,'L') ) then
       GFjob = 'left '
    else if ( leqi(tElec,'R') ) then
       GFjob = 'right'
    else
       call die("init electrode has received wrong job ID [L,R]: "//tElec)
    endif
    
    call timer('genGreen',1)

    if (IONode) then
       write(*,'(/,2a)') &
            "Creating Green's function file for: ",GFjob
    end if
    
    ! Read in all variables from the TSHS electrode file.
    ! Broadcasting within the routine is performed in MPI run
    call init_electrode_HS(tElec,NUsedAtoms,Gamma,xa,nua,nspin, &
         NBufAt, NA1, NA2, HSFile, &
         nua_E,nuo_E,maxnh_E,notot_E,xa_E,H_E,S_E,xij_E, &
         xijo_E,zconnect_E,numh_E,listhptr_E,listh_E,indxuo_E,  &
         lasto_E, &
         Ef_E,ucell_E)

    ! Calculate the k-points used in the electrode
    if (IONode) &
         write(*,*) "Electrodes with transport k-points &
         & (Bohr**-1) and weights:"
    do i = 1 , nkpnt
       ! From CONTACT to electrode k-point
       ! First convert to units of reciprocal vectors
       ! Then convert to 1/Bohr in the electrode unit cell coordinates
       call kpoint_convert(ucell,kpoint(:,i),ktmp,1)
       if ( NA1 > 1 ) ktmp(1) = ktmp(1)/real(NA1,dp)
       if ( NA2 > 1 ) ktmp(2) = ktmp(2)/real(NA2,dp)
       call kpoint_convert(ucell_E,ktmp,kpt,-1)
       if (IONode) write(*,'(i4,2x,3(E14.5))') i,&
            kpt(1),kpt(2),kweight(i)
    end do

! >>>>>>>>>>>>>>>>>>> Valence Band bottom Calculation <<<<<<<<<<<<<<<<< 
    if ( ElecValenceBandBot ) then
       ! Calculate the Valence Band Bottom for the electrode
       allocate(H00(nuo_E*nuo_E),S00(nuo_E*nuo_E))
       call memory('A','Z',nuo_E*nuo_E*2,'create_green')
       kpt = 0.0_dp
       call set_HS_matrix(Gamma,ucell_E,nua_E,nuo_E,notot_E,maxnh_E, &
            xij_E,numh_E,listhptr_E,listh_E,indxuo_E,H_E(:,1),S_E, &
            kpt,H00,S00)
       call matrix_symmetrize(nuo_E,H00,S00,Ef_E)
       
       ! We "reuse" H01 here. No need to create a temporary array.
       ! H01 is not used until later
       ! H01 => eigenvectors
       ! eig => eigenvalues
       allocate(H01(nuo_E*nuo_E))
       call memory('A','Z',nuo_E*nuo_E,'create_green')
       allocate(eig(nuo_E))
       call memory('A','D',nuo_E,'create_green')
       call cdiag(H00,S00,nuo_E,nuo_E,nuo_E,eig,H01,nuo_E,10,ierror)
       if ( IONode ) then
          if ( ierror == 0 ) then
             write(*,'(a,f10.4,a)')' Valence Band Bottom: ',eig(1)/eV,' eV'
          else
             write(*,'(a,i6)')' Error in calculating the Band Bottom: ',ierror
          end if
       end if
       call memory('D','Z',nuo_E*nuo_E*3,'create_green')
       deallocate(H00,S00,H01)
       nullify(H00,S00,H01)
       call memory('D','D',nuo_E,'create_green')
       deallocate(eig)
    end if
! >>>>>>>>>>>>>>>>>>>>> End of Valence band bottom calculation <<<<<<<


    ! Count number of used orbitals in the electrode
    ! Here we count the orbitals by using the option variable:
    ! TS.NumUsedAtoms[Left|Right]
    nuou_E = 0
    if( leqi(tElec,'L') ) then
       ! Left, we use the last atoms in the list
       do ia = nua_E - NUsedAtoms + 1, nua_E
          nuou_E = nuou_E + lasto_E(ia) - lasto_E(ia-1)
       end do
    else
       ! Right, the first atoms in the list
       do ia = 1 , NUsedAtoms
          nuou_E = nuou_E + lasto_E(ia) - lasto_E(ia-1)
       end do
    end if

    ! Show the number of used atoms and orbitals
    if ( IONode ) then
       write(*,'(a,i6,'' / '',i6)') ' Atoms available    / used atoms   : ', &
            nua_E,NUsedAtoms
       write(*,'(a,i6,'' / '',i6)') ' Orbitals available / used orbitals: ', &
            lasto_E(nua_E),nuou_E
    end if
    ! Clean up what can be cleaned up
    call memory('D','I',nua_E+1,'create_green')
    deallocate(lasto_E)

! FDN cell,kscell,kdispl added as dummys
! q,wq:
!  this is WHERE the initial q and wq points are generated.
!  they are read in by 'read_green' later on.
! They are in units of the reciprocal lattice vectors (hence suffix b)
    call mkqgrid(NA1,NA2,nq,qb,wq)

    if ( IONode ) then
! We show them in units of Bohr**-1
       write(*,'(a)') &
            ' q-points for expanding electrode (Bohr**-1):'
       do i=1,nq
          call kpoint_convert(ucell_E,qb(:,i),qpt,-1)
          write(*,'(i4,2x,3(E14.5))') i,qpt(1),qpt(2),wq(i)
       end do
       write(*,'(a,f14.5,1x,a)') &
            " Fermi level shift in electrode : ",chem_shift/eV,' eV'
    end if

    ! Initialize Green's function and Hamiltonian arrays
    allocate(GS(nuo_E*nuo_E))
    call memory('A','Z',nuo_E*nuo_E,'create_green')
    allocate(Hq(nuou_E*nuou_E,nq))
    allocate(Sq(nuou_E*nuou_E,nq))
    allocate(Gq(nuou_E*nuou_E,nq))
    call memory('A','Z',nuou_E*nuou_E*nq*3,'create_green')

    ! Allocate all H00,S00,H01 and S01 arrays
    allocate(H00(nuo_E*nuo_E),S00(nuo_E*nuo_E))
    allocate(H01(nuo_E*nuo_E),S01(nuo_E*nuo_E))
    call memory('A','Z',4*nuo_E*nuo_E,'create_Green')

    ! Reset bulk DOS
    do ispin = 1 , nspin
       do iEn = 1 , NEn
          ZBulkDOS(iEn,ispin) = dcmplx(0.d0,0.d0)
       end do
    end do

!******************************************************************
!           Start Green's function calculation
!******************************************************************
    
    if (IONode) then
       call io_assign(uGF)
       open(FILE=GFfile,UNIT=uGF,FORM='UNFORMATTED')

       ! Initial header for file
       write(uGF) GFTitle
       write(uGF) chem_shift,NEn
       write(uGF) RemUCellDistance
       write(uGF) NUsedAtoms,NA1,NA2,nkpnt,nq
       ! Write spin, ELECTRODE unit-cell
       write(uGF) nspin, ucell_E
       ! Write out the atomic coordinates of the used electrode
       if( leqi(tElec,'L') ) then
          ! Left, we use the last atoms in the list
          write(uGF) xa_E(:,nua_E-NUsedAtoms+1:nua_E)
       else
          ! Right, the first atoms in the list
          write(uGF) xa_E(:,1:NUsedAtoms)
       end if
       ! Notice that we write the k-points for the ELECTRODE
       ! Do a conversion here
       allocate(eig(nkpnt*3))
       call memory('A','D',nkpnt*3,'create_green')
       i = 0
       do ikpt = 1 , nkpnt
          ! Init kpoint, in reciprocal vector units ( from CONTACT ucell)
          call kpoint_convert(ucell,kpoint(:,ikpt),ktmp,1)
          ktmp(1) = ktmp(1)/real(NA1,dp)
          ktmp(2) = ktmp(2)/real(NA2,dp)
          ! Convert back to reciprocal units (to electrode ucell_E)
          call kpoint_convert(ucell_E,ktmp,kpt,-1)
          do j = 1 , 3
             i = i + 1
             eig(i) = kpt(j)
          end do
       end do
       write(uGF) contour(:)%c,contour(:)%w,eig,kweight,qb,wq
       call memory('D','D',nkpnt*3,'create_green')
       deallocate(eig)
       ! Write the number of USED orbitals for the calculation
       write(uGF) nuou_E
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
          call MPI_Recv_Init(Gq(1,1),nuou_E*nuou_E*nq,DAT_dcomplex, &
               i,i,MPI_Comm_World,reqs(i),MPIerror)
       end do
    else
       ! Create request handles for communication
       call MPI_Send_Init(Gq(1,1),nuou_E*nuou_E*nq,DAT_dcomplex, &
            0,Node,MPI_Comm_World,req,MPIerror)
    end if
#endif


! Spin loop .............................................. Spin loop
    spin_loop: do ispin = 1 , nspin
       
! TS k-point loop ........................................ k-point loop
       kpoint_loop: do ikpt = 1 , nkpnt
            
          ! Init kpoint, in reciprocal vector units ( from CONTACT ucell)
          call kpoint_convert(ucell,kpoint(:,ikpt),ktmp,1)
          ktmp(1) = ktmp(1)/real(NA1,dp)
          ktmp(2) = ktmp(2)/real(NA2,dp)
          ! Convert back to reciprocal units (to electrode ucell_E)
          call kpoint_convert(ucell_E,ktmp,kpt,-1)
          
          ! No need for reseting computational arrays
          ! They are completely overwritten

! Energy contour loop .................................... Energy loop
          Econtour_loop: do iEn = 1, NEn
             
#ifdef MPI
             ! Every node takes one energy point
             ! This asserts that IONode = Node == 0 will have iEn == 1
             ! Important !
             curNode = MOD(iEn-1,Nodes)
             E_Nodes: if ( curNode == Node ) then
#endif
             ZEnergy  = contour(iEn)%c
             ZSEnergy = ZEnergy-dcmplx(chem_shift,0.0_dp)
             
! loop over the repeated cell...
             q_loop: do iqpt = 1 , nq

                ! init qpoint in reciprocal lattice vectors
                call kpoint_convert(ucell_E,qb(:,iqpt),qpt,-1)

                ! Setup the transfer matrix and the intra cell at the k-point and q-point
                if ( RemUCellDistance ) then
                   call set_electrode_HS_Transfer(Gamma,nuo_E,maxnh_E, &
                        notot_E,nspin,H_E,S_E,xijo_E,xijo_E,zconnect_E,numh_E, &
                        listhptr_E,listh_E,indxuo_E,Ef_E, &
                        ispin, kpt, qpt, &
                        H00,S00,H01,S01)
                else
                   call set_electrode_HS_Transfer(Gamma,nuo_E,maxnh_E, &
                        notot_E,nspin,H_E,S_E,xij_E,xijo_E,zconnect_E,numh_E, &
                        listhptr_E,listh_E,indxuo_E,Ef_E, &
                        ispin, kpt, qpt, &
                        H00,S00,H01,S01)
                end if
                   
                
                ! This requires IONode = Node == 0 !
                if ( iEn .eq. 1 ) then
                   ! Copy over the Hamiltonian and overlap
                   if( leqi(tElec,'L') ) then
                      ! Left, we use the last orbitals
                      i = 0
                      do jo = 1 , nuou_E
                         do io = 1 , nuou_E
                            i = i + 1
                            Sq(i,iqpt) = S00(io+(nuo_E-nuou_E)+&
                                 nuo_E*(jo+(nuo_E-nuou_E)-1))
                            Hq(i,iqpt) = H00(io+(nuo_E-nuou_E)+&
                                 nuo_E*(jo+(nuo_E-nuou_E)-1)) +&
                                 chem_shift * Sq(i,iqpt)
                         end do  ! io
                      end do     ! jo
                   else
                      ! Right, the first orbitals
                      i=0
                      do jo = 1 , nuou_E
                         do io = 1 , nuou_E
                            i = i + 1
                            Sq(i,iqpt) = S00(io+nuo_E*(jo-1))
                            Hq(i,iqpt) = H00(io+nuo_E*(jo-1)) +&
                                 chem_shift * Sq(i,iqpt)
                         end do   ! io
                      end do      ! jo
                   end if
                   
                end if
                  
                ! Calculate the surface Green's function
                ! ZSenergy is Zenergy together with the chemical shift
                call surface_Green(tElec,nuo_E,ZSEnergy,H00,S00,H01,S01, &
                     GS,zdos)

                ! We also average the k-points.
                ZBulkDOS(iEn,ispin) = ZBulkDOS(iEn,ispin) + &
                     wq(iqpt)*zdos * kweight(ikpt)
                  
                ! Copy over surface Green's function
                if( leqi(tElec,'L') ) then
                   ! Left, we use the last orbitals
                   i = 0
                   do jo = 1 , nuou_E
                      do io = 1 , nuou_E
                         i = i + 1
                         Gq(i,iqpt) = GS(io+(nuo_E-nuou_E)+ &
                              nuo_E*(jo+(nuo_E-nuou_E)-1) )
                      end do           ! io
                   end do              ! jo
                else
                   ! Right, the first orbitals
                   i = 0
                   do jo = 1,nuou_E
                      do io = 1,nuou_E
                         i = i + 1
                         Gq(i,iqpt) = GS(io+nuo_E*(jo-1))
                      end do           ! io
                   end do              ! jo
                end if
                
             end do q_loop
             
             if (IONode) then
                ! Write out calculated information at E point

                write(uGF) iEn,contour(iEn)%c,contour(iEn)%w,ikpt
                
                if(iEn .eq. 1) then ! This will only occur
                   write(uGF) Hq
                   write(uGF) Sq
                end if
                
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
             write(uGF) iEn,contour(iEn)%c,contour(iEn)%w,ikpt
             call MPI_Start(reqs(curNode),MPIerror)
             call MPI_Wait(reqs(curNode),status,MPIerror)
             write(uGF) Gq
          end if
#endif             

          end do Econtour_loop
          
       end do kpoint_loop
       
    end do spin_loop
!*******************************************************************
!         Green's function calculation is done
!*******************************************************************


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
       write(*,'(a)') "Done creating '"//trim(GFFile)//"'."  
    end if
    
    ! Clean up computational arrays
    call memory('D','Z',nuou_E*nuou_E*nq*3+nuo_E*nuo_E,'create_green')
    deallocate(GS,Hq,Sq,Gq)

    ! Hamiltonians
    call memory('D','Z',4*nuo_E*nuo_E,'create_green')
    deallocate(H00,S00,H01,S01)

    if ( .not. Gamma ) then
       call memory('D','I',notot_E,'create_green')
       deallocate(indxuo_E)
       call memory('D','D',3*maxnh_E,'create_green')
       deallocate(xij_E)
       call memory('D','D',3*maxnh_E,'create_green')
       deallocate(xijo_E)
       call memory('D','I',maxnh_E,'create_green')
       deallocate(zconnect_E)
    end if

    call memory('D','I',nuo_E,'create_green')
    deallocate(numh_E)
    call memory('D','I',nuo_E,'create_green')
    deallocate(listhptr_E)
    call memory('D','I',maxnh_E,'create_green')    
    deallocate(listh_E)
    call memory('D','D',maxnh_E*nspin,'create_green')
    deallocate(H_E)
    call memory('D','D',maxnh_E,'create_green')
    deallocate(S_E)

#ifdef MPI
    ! Sum the bulkdensity of states
    ! Here we can safely use the array as temporary (Gq)
    allocate(Gq(NEn,nspin))
    call memory('A','Z',NEn*nspin,'create_green')
    Gq = 0.0_dp
    call MPI_AllReduce(ZBulkDOS(1,1),Gq(1,1),NEn*nspin, DAT_dComplex, MPI_Sum, &
         MPI_Comm_World,MPIerror)
    ZBulkDOS = Gq
    call memory('D','Z',NEn*nspin,'create_green')
    deallocate(Gq)
#endif

    call timer('genGreen',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS create_Green' )
#endif
    
  end subroutine create_Green


!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
  subroutine init_electrode_HS(tElec,NUsedAtoms,Gamma,xa_sys,nua_sys,nspin_sys, &
       NBufAt,NA1,NA2, &
       HSfile,nua,nuo,maxnh,notot,xa,H,S,xij,xijo,zconnect, &
       numh,listhptr,listh,indxuo,lasto, &
       Ef,ucell)

    use precision, only: dp
    use fdf, only : leqi
    use sys, only : die
    use units, only : Ang
    use parallel, only : IONode, Node,Nodes
    use m_ts_io, only  : ts_iohs
    use files, only: label_length
#ifdef MPI
    use mpi_siesta, only: MPI_Comm_World, MPI_LOR
    use mpi_siesta, only: MPI_Bcast,MPI_Barrier
    use mpi_siesta, only: DAT_double => MPI_double_precision
    use mpi_siesta, only: MPI_Logical,MPI_Integer
    use mpi_siesta, only: MPI_Reduce
#endif
#ifdef TBTRANS
    use m_tbt_kpoints, only : ts_kscell => kscell
    use m_tbt_kpoints, only : ts_kdispl => kdispl
    use m_tbt_kpoints, only : ts_gamma_scf => Gamma
#else
    use m_ts_kpoints, only : ts_gamma_scf, ts_kscell, ts_kdispl
#endif


! ***********************
! * INPUT variables     *
! ***********************
    character(len=1)     :: tElec   ! 'L' for Left electrode, 'R' for right
    integer, intent(in)  :: NUsedAtoms ! The number of atoms used...
    integer, intent(in)  :: nua_sys ! Full system count of atoms in unit cell
    real(dp), intent(in) :: xa_sys(3,nua_sys) ! Coordinates in the system for the TranSIESTA routine
    integer, intent(in)  :: nspin_sys ! spin in system
    character(len=label_length+5),intent(in) :: HSfile !H,S parameter file 
    integer, intent(in)  :: NBufAt,NA1,NA2 ! Buffer atoms, and repetitions 
! ***********************
! * OUTPUT variables    *
! ***********************
    logical                           :: Gamma
    integer                           :: nua,nuo,maxnh ! Unit cell atoms / orbitals / Hamiltonian size
    real(dp), dimension(:,:), pointer :: xa ! The atomic coordinates
    real(dp), dimension(:,:), pointer :: H,xij,xijo ! Hamiltonian, differences with unitcell, differences without unitcell
    real(dp), dimension(:),   pointer :: S ! Overlap
    integer,  dimension(:),   pointer :: zconnect ! Has 0 values for indices where there are no z-connection.
    integer,  dimension(:),   pointer :: numh,listhptr,listh,indxuo,lasto
    real(dp) :: Ef ! Efermi
    real(dp), intent(inout)           :: ucell(3,3) ! The unit cell

! ***********************
! * LOCAL variables     *
! ***********************
! >>>> Related to the Electrode TSHS
    character(len=5) :: GFjob
    integer :: notot  ! Total orbitals in all supercells
    integer :: nspin  ! The spin polarization
    integer,  dimension(:),   pointer :: isa ! atomic species
    logical :: ts_gamma ! Read gamma from file
    real(dp) :: kdispl(3)
    integer  :: kscell(3,3)
    real(dp) :: qtot,Temp ! total charges, electronic temperature
    integer :: dummy1,dummy2 ! dummy variables


    logical :: eXa ! Errors when testing for position of electrodes
    real(dp),dimension(3) :: xa_o,xa_sys_o ! Origo coordinates of the electrodes
    real(dp) :: zc
    real(dp), dimension(:,:), allocatable :: xo
    real(dp), dimension(3,3) :: recell ! Reciprocal cell without 2Pi, used for zconnect
    integer :: sysElec ! the first atom of the electrode in the SYSTEM setup
    integer :: elecElec ! The first atom in the electrode in the ELECTRODE setup
    integer :: i,j,ia,iuo,juo,ind,iaa !Loop counters
    character(8) :: iotask

#ifdef MPI
    logical :: eXa_buff
    integer :: MPIerror
#endif
!=======================================================================
    real(dp), parameter :: EPS = 1.0d-4
!=======================================================================

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE init_elec_HS' )
#endif

    nullify(H,S)
    nullify(xij,xijo,xa,isa)
    nullify(numh,listhptr,listh,indxuo)
    nullify(lasto)
    nullify(zconnect)

    if ( IONode ) then
       iotask = 'READ'
       call ts_iohs(iotask,Gamma, .false., nuo, notot, nspin, &
            indxuo, maxnh, numh, listhptr, listh, H, S, qtot, Temp, &
            xij, label_length+5, HSfile, nua, lasto, isa, Ef, &
            ucell, kscell, kdispl, ts_gamma, xa, dummy1,dummy2, &
            check_kcell=.false.) 
       ! Checking the electrode k-grid against the system k-grid...
       ! This will also check for the repetition...
       ! This check is not advisable in the TBtrans utility.
       ! This is because the k-grid is often more fine when calculating
       ! the transmission function.
       ! Perhaps this should only check that the electrode k-sampling
       ! is AT LEAST as good as the system k-grid for a TranSiesta run?
#ifndef TBTRANS
       eXa = .false.
       do j = 1 , 3
          do i = 1 , 3
             if ( j == 1 ) then
                eXa = eXa .or. ( kscell(i,j) /= ts_kscell(i,j)*NA1 )
             else if ( j == 2 ) then
                eXa = eXa .or. ( kscell(i,j) /= ts_kscell(i,j)*NA2 )
             else
                eXa = eXa .or. ( kscell(i,j) /= ts_kscell(i,j) )
             end if
          end do
          eXa = eXa .or. ( kdispl(j) /= ts_kdispl(j) )
       end do
       if ( eXa ) then
          write(*,'(a)') 'Incompatible k-grids...'
          write(*,'(a)') 'Electrode file k-grid:'
          do j = 1 , 3
             write(*,'(3(i4,tr1),f8.4)') (kscell(i,j),i=1,3),kdispl(j)
          end do
          write(*,'(a)') 'System k-grid:'
          do j = 1 , 3
             write(*,'(3(i4,tr1),f8.4)') (ts_kscell(i,j),i=1,3),ts_kdispl(j)
          end do
          write(*,'(a)') 'Electrode file k-grid should be:'
          kscell(:,1) = ts_kscell(:,1) * NA1
          kscell(:,2) = ts_kscell(:,2) * NA2
          kscell(:,3) = ts_kscell(:,3)
          do j = 1 , 3
             write(*,'(3(i4,tr1),f8.4)') (kscell(i,j),i=1,3),ts_kdispl(j)
          end do
          call die('Incompatible electrode k-grids') 
       end if
#endif

       ! Deallocate isa, not needed anymore
       call memory('D','I',nua,'elec_HS')
       deallocate(isa)

       if( nspin_sys /= nspin ) then
          write(*,*) "Spin of electrode '"//trim(HSfile) &
               //"' is different."
          write(*,*) "Spin electrode / Spin system",nspin,nspin_sys
          call die('Wrong spin! Check system and electrode!')
       end if
    end if

! Do communication of variables
#ifdef MPI
    call MPI_Bcast(nua,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(nuo,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(notot,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(nspin,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(maxnh,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(Ef,1,DAT_double,0, MPI_Comm_World,MPIerror)
    call MPI_Bcast(ucell(1,1),3*3,DAT_double,0, &
         MPI_Comm_World,MPIerror)
    if ( .not. IONode ) then
       if ( .not. Gamma ) then 
! they behave as dummy arrays in case of Gamma == .true.
! However, the electrode must be Gamma == .false. as we need the transfer matrix
! TODO in options start-up check that the electrode is in fact a Gamma == .false. calculation
! This can be done together with retrieving the nua_E variable.
          allocate(indxuo(notot))
          call memory('A','I',notot,'elec_HS')
          allocate(xij(3,maxnh))
          call memory('A','D',maxnh*3,'elec_HS')
       end if
       allocate(xa(3,nua))
       call memory('A','D',3*nua,'elec_HS')
       allocate(lasto(0:nua))
       call memory('A','I',1+nua,'elec_HS')
       allocate(numh(nuo))
       call memory('A','I',nuo,'elec_HS')
       allocate(listhptr(nuo))
       call memory('A','I',nuo,'elec_HS')
       allocate(listh(maxnh))
       call memory('A','I',maxnh,'elec_HS')
       allocate(H(maxnh,nspin))
       call memory('A','D',maxnh*nspin,'elec_HS')
       allocate(S(maxnh))
       call memory('A','D',maxnh,'elec_HS')
    end if
    call MPI_Bcast(xa(1,1),3*nua,DAT_Double,0, &
         MPI_Comm_World,MPIerror)
    if ( .not. Gamma ) then
       call MPI_Bcast(indxuo,notot,MPI_Integer,0, MPI_Comm_World,MPIerror)
       call MPI_Bcast(xij(1,1),3*maxnh,DAT_Double,0, &
            MPI_Comm_World,MPIerror)
    end if
    call MPI_Bcast(lasto(0),1+nua,MPI_Integer,0, MPI_Comm_World,MPIerror)
    call MPI_Bcast(numh,nuo,MPI_Integer,0, MPI_Comm_World,MPIerror)
    call MPI_Bcast(listhptr,nuo,MPI_Integer,0, MPI_Comm_World,MPIerror)
    call MPI_Bcast(listh,maxnh,MPI_Integer,0, MPI_Comm_World,MPIerror)
    call MPI_Bcast(H(1,1),maxnh*nspin,DAT_Double,0, &
         MPI_Comm_World,MPIerror)
    call MPI_Bcast(S,maxnh,DAT_Double,0, MPI_Comm_World,MPIerror)
#endif


    ! Do a recheck if the electrode file has been overwritted or??
    if ( NUsedAtoms > nua ) then
       write(*,*) "# of requested atoms is larger than available."
       write(*,*) "Requested: ",NUsedAtoms
       write(*,*) "Available: ",nua
       call die("Error on requested atoms.")
    end if


! Create reciprocal cell, without 2Pi
    call reclat(ucell,recell,0)

    if( leqi(tElec,'L') ) then
       GFjob = 'Left'
       sysElec = NbufAt + 1
       elecElec = nua - NUsedAtoms + 1
    else if ( leqi(tElec,'R') ) then
       GFjob = 'Right'
       sysElec = nua_sys - NbufAt - NUsedAtoms * NA1 * NA2 + 1
       elecElec = 1
    else
       call die("init electrode has received wrong job ID [L,R]: "//tElec)
    endif

! Print out structural information of the system versus the electrode
    struct_info: if ( IONode ) then
       write(*,*) trim(GFjob)//' unit cell (Ang):'
       do j=1,3
          write(*,'(3F8.4)') (ucell(i,j)/Ang,i=1,3)
       end do

       write(*,'(a,t35,a)') &
            " Structure of the "//trim(GFjob)//" electrode","| System electrode:"
       write(*,'(t3,3a10,''  |'',3a10)') &
            "X (Ang)","Y (Ang)","Z (Ang)", &
            "X (Ang)","Y (Ang)","Z (Ang)"

       ! Save origo of System electrode
       xa_sys_o(:) = xa_sys(:,sysElec)
       ! Save origo of electrode 
       xa_o(:) = xa(:,elecElec)

       ! Initialize error parameter
       eXa = .false.
       iaa = sysElec
       do ia = elecElec , elecElec + NUsedAtoms - 1
          do j=0,NA2-1
             do i=0,NA1-1
                write(*,'(t3,3f10.5,''  |'',3f10.5)') &
                     (xa(1,ia)-xa_o(1)+ucell(1,1)*i+ucell(1,2)*j)/Ang, &
                     (xa(2,ia)-xa_o(2)+ucell(2,1)*i+ucell(2,2)*j)/Ang, &
                     (xa(3,ia)-xa_o(3))/Ang, &
                     (xa_sys(1,iaa)-xa_sys_o(1))/Ang, &
                     (xa_sys(2,iaa)-xa_sys_o(2))/Ang, &
                     (xa_sys(3,iaa)-xa_sys_o(3))/Ang
                eXa=eXa.or.abs(xa(1,ia)-xa_o(1) + &
                     ucell(1,1)*i+ucell(1,2)*j - &
                     xa_sys(1,iaa)+xa_sys_o(1)) > EPS
                eXa=eXa.or.abs(xa(2,ia)-xa_o(2) + &
                     ucell(2,1)*i+ucell(2,2)*j - &
                     xa_sys(2,iaa)+xa_sys_o(2)) > EPS
                eXa=eXa.or.abs(xa(3,ia)-xa_o(3) - &
                     xa_sys(3,iaa)+xa_sys_o(3)) > EPS
                iaa = iaa + 1
             end do
          end do
       end do
       if ( eXa ) then
          write(*,'(a)') "Coordinates from the electrode repeated out to an FDF file"
          write(*,'(t3,3a20)') &
               "X (Ang)","Y (Ang)","Z (Ang)"
          iaa = sysElec
          do ia = elecElec , elecElec + NUsedAtoms - 1
             do j=0,NA2-1
                do i=0,NA1-1
                   write(*,'(t2,3(tr1,f20.12))') &
                        (xa(1,ia)+ucell(1,1)*i+ucell(1,2)*j)/Ang, &
                        (xa(2,ia)+ucell(2,1)*i+ucell(2,2)*j)/Ang, &
                        (xa(3,ia))/Ang
                end do
             end do
          end do
          call die("The electrodes are not situated in the same coordinates. Please correct.")
       end if

    end if struct_info

! Create them for passing to other routines 
! in such case, they are dummy arrays (this occurs only if Gamma .eq. .true.
    allocate(xijo(3,maxnh))
    call memory('A','D',3*maxnh,'elec_HS')
    allocate(zconnect(maxnh))
    call memory('A','I',maxnh,'elec_HS')
    
    ! Initialize the error parameter
    eXa = .false.

    ! Create the z-connect array
    zconnect_gamma: if ( .not. Gamma ) then

       ! temporary array in this part
       allocate(xo(3,nuo))
       call memory('A','D',3*nuo,'elec_HS')

! We now create zconnect to contain 0 for interconnects without
! z-direction
! This needs to be the full electrode, no matter NUsedAtoms!

! Create xo array (orbital coordinates)
       do ia = 1 , nua
          do i = lasto(ia-1)+1 , lasto(ia)
             xo(1,i) = xa(1,ia)
             xo(2,i) = xa(2,ia)
             xo(3,i) = xa(3,ia)
          end do           !i
       end do              !ia in uc

       do iuo = 1 , nuo
          do j = 1 , numh(iuo)
             ind = listhptr(iuo) + j
             juo = indxuo(listh(ind))
             xijo(:,ind) = xij(:,ind)-(xo(:,juo)-xo(:,iuo))
             zc = 0.0_dp
             do i = 1 , 3
! recell is already without 2*Pi
                zc = zc + xijo(i,ind) * recell(i,PropDir)
             end do
             zconnect(ind) = nint(zc)

             if ( abs(zconnect(ind)) > 1 ) then
                eXa = .true.
             end if
          end do
       end do
       call memory('D','D',3*nuo,'elec_HS')
       deallocate(xo)
    end if zconnect_gamma
    
#ifdef MPI
    eXa_buff = eXa
    call MPI_Reduce(eXa_buff,eXa,1, MPI_LOGICAL,MPI_LOR, &
         0,MPI_Comm_World,MPIerror)
#endif

    if ( IONode .and. eXa ) then
       write(0,*) "WARNING: Connections across 2 unit cells or more &
            &in the transport direction."
       write(0,*) "WARNING: This is inadvisable."
       write(0,*) "WARNING: Please increase the electrode size &
            &in the transport direction."
       write(0,*) "WARNING: Will proceed without further notice."
       write(*,*) "WARNING: Connections across 2 unit cells or more &
            &in the transport direction."
       write(*,*) "WARNING: This is inadvisable."
       write(*,*) "WARNING: Please increase the electrode size &
            &in the transport direction."
       write(*,*) "WARNING: Will proceed without further notice."
    end if

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS init_elec_HS' )
#endif

  end subroutine init_electrode_HS


!**********
! Create the Hamiltonian for the electrode as well
! as creating the transfer matrix.
!**********
  subroutine set_electrode_HS_Transfer(Gamma,nuo,maxnh,notot,nspin, &
       H,S,xij,xijo,zconnect,numh, &
       listhptr,listh,indxuo,Ef, &
       ispin,k,q,Hk,Sk,Hk_T,Sk_T)
    use sys, only : die
    use precision, only : dp
! ***********************
! * INPUT variables     *
! ***********************
    logical                           :: Gamma ! Is it a Gamma Calculation?
    integer                           :: nuo ! Unit cell orbitals
    integer                           :: maxnh,notot,nspin ! Hamiltonian size / total orbitals / spins
    real(dp)                          :: H(maxnh,nspin) ! Hamiltonian
    real(dp)                          :: S(maxnh) ! Overlap
    real(dp)                          :: xij(3,maxnh) ! differences with unitcell, differences with unitcell
    real(dp)                          :: xijo(3,maxnh) ! differences with unitcell, differences without unitcell
    integer                           :: zconnect(maxnh) ! 0 for no connection in z-direction
    integer                           :: numh(nuo),listhptr(nuo)
    integer                           :: listh(maxnh),indxuo(notot)
    real(dp)                          :: Ef ! Efermi
    integer                           :: ispin
    real(dp), dimension(3)            :: k   ! k-point in [1/Bohr]
    real(dp), dimension(3)            :: q   ! expansion k-point in [1/Bohr]
! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), dimension(nuo*nuo)   :: Hk,Sk,Hk_T,Sk_T

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: kqxij
    complex(dp) :: cphase
    integer :: i,j,iuo,juo,ind

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE elec_HS_Transfer' )
#endif

    if ( Gamma ) then
       write(*,*) 'Transfer matrix not possible with Gamma-calculation.'
       call die("Transfer matrix not possible with Gamma-calculation")
    end if
    
! Initialize arrays
    do i = 1,nuo*nuo
       Hk(i) = dcmplx(0.d0,0.d0)
       Sk(i) = dcmplx(0.d0,0.d0)
       Hk_T(i) = dcmplx(0.d0,0.d0)
       Sk_T(i) = dcmplx(0.d0,0.d0)
    enddo

    do iuo = 1 , nuo
       do j = 1 , numh(iuo)
          ind = listhptr(iuo) + j
          juo = indxuo(listh(ind))
          kqxij = &
               k(1)       * xij(1,ind) + &
               k(2)       * xij(2,ind) + &
               k(3)       * xij(3,ind) - &
               k(PropDir) * xij(PropDir,ind) + &
               q(1)       * xijo(1,ind) + &
               q(2)       * xijo(2,ind) + &
               q(3)       * xijo(3,ind) - &
               q(PropDir) * xijo(PropDir,ind)

          cphase = cdexp(dcmplx(0d0,1d0)*kqxij )
          
          i = iuo+(juo-1)*nuo
          if (zconnect(ind).eq.0) then
             Hk(i) = Hk(i)+H(ind,ispin)*cphase
             Sk(i) = Sk(i)+S(ind)*cphase
          else if (zconnect(ind).eq.1) then
             Hk_T(i) = Hk_T(i)+H(ind,ispin)*cphase
             Sk_T(i) = Sk_T(i)+S(ind)*cphase
          endif
          
       enddo
    enddo

!
!     Symmetrize and make EF the energy-zero*!!!
!
    do iuo = 1,nuo
       do juo = 1,iuo-1
          i = iuo+(juo-1)*nuo
          j = juo+(iuo-1)*nuo

          Sk(j) = 0.5d0*( Sk(j) + dconjg(Sk(i)) )
          Sk(i) = dconjg(Sk(j))

          Hk(j) = 0.5d0*( Hk(j) + dconjg(Hk(i)) ) &
               - Ef*Sk(j)
          Hk(i) = dconjg(Hk(j))

          ! Transfer matrix
          Hk_T(i) = Hk_T(i) - Ef*Sk_T(i)
          Hk_T(j) = Hk_T(j) - Ef*Sk_T(j)

       enddo
       
       i = iuo+(iuo-1)*nuo
       Sk(i) = Sk(i) - dcmplx(0d0,1d0)*dimag(Sk(i))
       
       Hk(i) = Hk(i) - dcmplx(0d0,1d0)*dimag(Hk(i)) &
            - Ef*Sk(i)
       
       ! Transfer matrix
       Hk_T(i) = Hk_T(i) - Ef*Sk_T(i)
    enddo

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS elec_HS_Transfer' )
#endif

!-----------------------------------------------------------------
  end subroutine set_electrode_HS_Transfer


! ******************************************************
! * ts_mkqgrid creates a q-point grid which is used    *
! * for expanding a smaller unit cell into a larger    *
! * coinciding one.                                    *
! *                                                    *
! * Example:                                           *
! *   Electrode of 1 atom, can be repeated 3x3 times   *
! *   to form a 3x3 atom layer in the TranSIESTA       *
! *   calculation.                                     *
! *                                                    *
! * (Re-)Introduced by Nick Papior Andersen            *
! ******************************************************
  subroutine mkqgrid(NA1,NA2,nq,q,wq)
    use precision, only : dp

! ***********************
! * INPUT variables     *
! ***********************
    integer , intent(in)     :: NA1,NA2  ! no. repetitions of simple unitcell in A1,A2 directions   

! ***********************
! * OUTPUT variables    *
! ***********************
    integer , intent(out)          :: nq      ! no. q-points (<= NA1*NA2 for gamma)
    real(dp), pointer              :: q(:,:)  ! q-points
    real(dp), pointer              :: wq(:)   ! weight of q-points (k_||)

! ***********************
! * LOCAL variables     *
! ***********************
    integer                  :: i,j,iq

    nq = NA1*NA2                !initial value

! To comply with new standard 3-dimension regime
    allocate(q(3,nq))
    call memory('A','D',3*nq,'mkqgrid')
    allocate(wq(nq))
    call memory('A','D',nq,'mkqgrid')
    ! Initialize to 0.0
    q(:,:)  = 0.0_dp
    wq(:)   = 0.0_dp
    iq = 0
    do i=1,NA1
       do j=1,NA2
          iq = iq+1
          q(1,iq)= 1.0_dp*(i-1)/real(NA1,dp)
          q(2,iq)= 1.0_dp*(j-1)/real(NA2,dp)
          q(3,iq)= 0.0_dp
          wq(iq) = 1.0_dp/real(NA1*NA2,dp)
       end do
    end do
    
  end subroutine mkqgrid


end module m_ts_electrode
