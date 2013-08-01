subroutine amn( ispin )

  use precision,          only: dp                  ! Real double precision type
  use parallel,           only: Nodes               ! Total number of Nodes
  use parallel,           only: Node                ! Local Node
  use parallel,           only: IONode              ! Input/output node
  use atomlist,           only: rmaxo               ! Max. cutoff atomic orbital
  use siesta_geom,        only: scell               ! Lattice vector of the
                                                    !   supercell in real space
  use siesta_geom,        only: na_s                ! Number of atoms in the
                                                    !   supercell
  use siesta_geom,        only: xa                  ! Atomic positions
  use m_siesta2wannier90, only: latvec              ! Lattice vectors in real 
                                                    !   space
  use m_siesta2wannier90, only: numproj             ! Total number of projectors
  use m_siesta2wannier90, only: projections         ! Trial projection functions
  use m_siesta2wannier90, only: numkpoints          ! Total number of k-points
                                                    !   for which the overlap 
                                                    !   of the projection
                                                    !   function with the
                                                    !   eigenstate at k will be
                                                    !   computed
  use m_siesta2wannier90, only: kpointsfrac         ! List of k points relative
                                                    !   to the reciprocal 
                                                    !   lattice vectors.
                                                    !   First  index: component
                                                    !   Second index: k-point 
                                                    !      index in the list
  use m_siesta2wannier90, only: Amnmat              ! Matrix of the overlaps of 
                                                    !   trial projector funtions
                                                    !   with Eigenstates of the
                                                    !   Hamiltonian
  use trialorbitalclass,  only: trialorbital        ! Derived type to define the
                                                    !    localized trial
                                                    !    orbitals
  use m_matel_registry,   only: register_in_tf_pool ! Subroutine that assigns a
                                                    !    global index to the 
                                                    !    trial projection
                                                    !    functions
  use m_new_matel,        only: new_matel           ! New MATEL implementation 
                                                    !   with the global indices 
                                                    !   of the radial functions 
                                                    !   as inputs
  use atmfuncs,           only: orb_gindex          ! Subroutine that gives
                                                    !   the global index of an
                                                    !   atomic orbital

!
! Variables for the diagonalization
!
  use m_spin,             only: nspin        ! Number of spin components
  use atomlist,           only: no_s         ! Number of orbitals in supercell
                                             ! NOTE: When running in parallel,
                                             !   this is core independent
  use atomlist,           only: no_l         ! Number of orbitals in local node
                                             ! NOTE: When running in parallel,
                                             !   this is core dependent
                                             !   Sum_{cores} no_l = no_u
  use atomlist,           only: no_u         ! Number of orbitals in unit cell
                                             ! NOTE: When running in parallel,
                                             !   this is core independent
  use atomlist,           only: indxuo       ! Index of equivalent orbital in 
                                             !   the unit cell
  use sparse_matrices,    only: maxnh        ! Maximum number of orbitals
                                             !   interacting
                                             ! NOTE: While running in parallel,
                                             !   maxnh changes from one core to 
                                             !   the other
  use sparse_matrices,    only: numh         ! Number of nonzero element of each
                                             !   row of the hamiltonian matrix 
  use sparse_matrices,    only: listh        ! Nonzero hamiltonian-matrix elemen
  use sparse_matrices,    only: listhptr     ! Pointer to start of each row 
                                             !   of the hamiltonian matrix
  use sparse_matrices,    only: H            ! Hamiltonian matrix in sparse form
  use sparse_matrices,    only: S            ! Overlap matrix in sparse form
  use sparse_matrices,    only: xijo         ! Vectors between orbital centers
  use densematrix,        only: Haux         ! Hamiltonian matrix in dense form
  use densematrix,        only: Saux         ! Overlap matrix in dense form
  use densematrix,        only: psi          ! Coefficients of the wave function
                                             ! NOTE: while running in parallel,
                                             ! each core knows about the 
                                             ! wave functions of no_l bands
  use alloc,              only: re_alloc     ! Reallocation routines
  use m_noccbands,        only: noccupied    ! Number of occupied bands for a 
                                             !   given spin direction

  use parallelsubs,         only: LocalToGlobalOrb ! Converts an orbital index
                                                 !   in the local frame 
                                                 !   to the global frame
  use parallelsubs,         only: GetNodeOrbs    ! Calculates the number of
                                                 !   orbitals stored on the 
                                                 !   local Node.


#ifdef MPI
      use mpi_siesta
#endif


  implicit none

  type orbitallinkedlist
    real(dp),dimension(3)           :: center
    integer                         :: specie
    integer                         :: specieindex
    integer                         :: globalindex
    type(orbitallinkedlist),pointer :: nextitem
  end type

! Passed arguments
  integer,  intent(in) :: ispin                   ! Spin component

! Local variables
  integer  :: numproj_l       ! Number of projections to be computed locally
                              !   in this node

  integer  :: minpernode      ! Minimum number of projections per node
  integer  :: remainder       ! Remainder of unassigned projections
  integer  :: iproj           ! Counter for loop on projections
  integer  :: ik              ! Counter for loop on kpoints
  integer  :: io              ! Counter for loop on orbital 
  integer  :: iuo             ! Counter for loop on orbital 
  integer  :: imu             ! Counter for loop on orbital 
  integer  :: iband           ! Counter for loop on bands
  integer  :: gindex          ! Global index of the trial projector function
                              !   in the list of functions that will be
                              !   evaluated in Matel
  integer  :: globalindexorbital ! Global index of the neighbour atomic orbital
                              !   in the list of functions that will be
                              !   evaluated in Matel
  integer  :: globalindexproj ! Global index of the neighbour atomic orbital
                              !   in the list of functions that will be
                              !   evaluated in Matel
  integer  :: indexproj       ! Index of the projector 
                              !   This index runs from 1 to the total number
                              !   of projections
  integer  :: nhs             ! Variable to dimension the Hamiltonian / Overlap
  integer  :: npsi            ! Variable to dimension the coefficient vector
  integer  :: nbands          ! Number of occupied bands
  real(dp) :: trialcenter(3)  ! Position where the trial function is centered
                              !   (in Bohr)
  real(dp) :: trialrcut       ! Cutoff radius of the trial function
  real(dp) :: r12(3)          ! Relative position of the trial function with
                              !   respect to the neighbour orbital
  real(dp) :: overlap         ! Overlap between the trial function and a 
                              !   given atomic orbital
  real(dp) :: gradient(3)     ! Grandient of the overlap between 
                              !   the trial function and a given atomic orbital
  real(dp) :: kvector(3)      ! k-point vector for which the Overlap matrix 
                              !   between the projection function and the 
                              !   eigenvector of the Hamiltonian will be
                              !   computed
  real(dp) :: phase           ! Product of the k-vector with the position
                              !   where the neighbour orbital is centered

  real(dp), dimension(:), pointer :: epsilon ! Eigenvalues of the Hamiltonian
  complex(dp), dimension(:,:), pointer :: psiloc ! Coefficients of the wave
                                             !   function
  complex(dp) :: exponential                 ! Exponential of exp( i * phase )
  complex(dp) :: cstar                       ! Conjugate of the coefficient
                                             !   of the wave function      
  complex(dp), parameter :: iu = cmplx(0.0_dp,1.0_dp,kind=dp) ! Imaginary unit

  type(trialorbital) :: tf                   ! Projetion function

  type(orbitallinkedlist),pointer   :: item  ! Variable that contains all the
                                             !   information of the atomic
                                             !   orbitals neighbours of a given
                                             !   trial projection function
  integer, dimension(:), allocatable :: projector_gindex 
                                             ! Array that gives the global index
                                             !   of a projector in the list of
                                             !   functions that will be
                                             !   evaluated by MATEL

#ifdef MPI
  integer     :: MPIerror
  integer     :: inode                       ! Counter for the loop on nodes
  integer     :: nbands_max_loc              ! Maximum number of bands stored
                                             !   on a node
  integer     :: nbands_loc                  ! Number of bands stored
                                             !   on the local node
  integer     :: iband_global                ! Global index for a band
  complex(dp), dimension(:,:), pointer :: psitmp ! Temporal array used to
                                             !   broadcast the coefficients 
                                             !   of the wave function
  complex(dp), dimension(:,:,:), pointer :: auxloc ! Temportal array for the
                                             !   the global reduction of Amnmat
#endif 



  external :: timer

! Start time counter
  call timer( 'amn', 1 )

! Allocate memory related with the overlap matrix between the trial projection
! function and the Hamiltonian eigenstate.
! These matrix will depend on three indices (see paragraph after Eq. (16) 
! of the review by
! N. Marzari et al., Review of Modern Physics 84, 1419 (2012),
! or Eq. (1.8) of the Wannier Users Guide, Version 1.2
! $A_{m n}^{(\vec{k})} = 
!    \langle \psi_{m \vec{k}} \vbar g_{n} \rangle =
!    \sum_{\mu} \sum_{lcell} c_{\mu m}^{\ast} (\vec{k}) \times
!    exp^{i \vec{k} \cdot \left( \tau_{\mu} + \vec{R}_{lcell} \right)} \times
!    \langle \phi_{\mu} (\vec{r}-\tau_{\mu} - \vec{r}_{lcell} \vbar g_{n}\rangle
! where $m$ runs between 1 and the number of occupied bands,
! $n$ runs between 1 and the number of projection functions, 
! $\vec{k}$ runs over all the k-points in the first BZ where these matrices 
! will be computed 
  nbands = noccupied( ispin )
  nullify( Amnmat )
  call re_alloc( Amnmat,            &
 &               1, nbands,         &
 &               1, numproj,        &
 &               1, numkpoints,     &
 &               'Amnmat',          &
 &               'amn'           )


! Allocate memory related with the dense matrices that will be used in
! the diagonalization routines.
! (Hamiltonian, Overlap and Coefficient vectors)
! These matrices are defined in the module dense matrix
  nhs  = 2 * no_u * no_l
  npsi = 2 * no_u * no_l

  call re_alloc( Haux, 1, nhs,  'Haux', 'Amn' )
  call re_alloc( Saux, 1, nhs,  'Saux', 'Amn' )
  call re_alloc( psi,  1, npsi, 'psi',  'Amn' )

! Initialise psi
  do io = 1, npsi
    psi(io)     = 0.0_dp
  enddo

! Allocate memory related with the eigenvalues of the Hamiltonian (epsilon)
! and with a local variable where the coefficients of the eigenvector at the
! k-point will be stored
  nullify( epsilon, psiloc )
  call re_alloc( epsilon,      1, no_u,            'epsilon', 'Amn' )
  call re_alloc( psiloc,       1, no_u, 1, no_l,   'psiloc',  'Amn' )

! Allocate the temporal variable that will be used to broadcast
! the coefficients of the wave function
! Compute the number of bands stored on Node 0
! We assume that this is the maximum number of orbitals that will be stored
! on any node.
  call GetNodeOrbs( no_u, 0, Nodes, nbands_max_loc)
  nullify( psitmp )
  call re_alloc( psitmp, 1, no_u, 1, nbands_max_loc, &
 &               name='psitmp', routine='amn' )

! Assign a global index to the trial functions
! The same global index is assigned to a given trial function in all the nodes
  allocate ( projector_gindex(numproj) )
  do iproj = 1, numproj
    tf = projections( iproj )
    call register_in_tf_pool( tf, gindex )
    projector_gindex( iproj ) = gindex
  enddo

! Calculate the minimum number of projections per node
  minpernode = numproj / Nodes

! Find the remainder of unassigned projections
  remainder = numproj - minpernode * Nodes 

! Workout the local number of projectors
  numproj_l = minpernode 
  if (Node .lt. remainder) numproj_l = numproj_l + 1

kpoints:                                                        &
  do ik = 1, numkpoints
!   Compute the wave vector in bohr^-1 for every vector in the list
!   (done in the subroutine getkvector).
!   Remember that kpointsfrac are read from the .nnkp file in reduced units,
!   so we have to multiply then by the reciprocal lattice vector.
    call getkvector( kpointsfrac(:,ik), kvector )
    call diagpol( ispin, nspin, no_l, no_s, no_u,                    &
 &                maxnh, numh, listhptr, listh, H, S, xijo, indxuo,  &
 &                kvector, epsilon, psi, 2, Haux, Saux )
    iuo = 0
    do iband = 1, no_l
      do io = 1, no_u
        psiloc(io,iband) = cmplx(psi(iuo+1), psi(iuo+2),kind=dp)
        iuo = iuo + 2
      enddo
    enddo

! Loop on the projections that will be computed in the local node
    do iproj = 1, numproj_l 
!     Identify the global index 
      if ( remainder .eq. 0 ) then
        indexproj = minpernode * Node + iproj
      else 
        if ( Node .lt. remainder ) then 
          indexproj = ( minpernode + 1 ) * Node + iproj 
        else
          indexproj = ( minpernode + 1 ) * remainder + &
 &                    ( Node - remainder )+ 1
        endif
      endif

!     Find the global index of the projector in the list of radial functions
!     that will be evaluated by MATEL
      globalindexproj = projector_gindex(indexproj) 

!     Find where the trial function is centered
      trialcenter = projections(indexproj)%center
    
!     Find the cutoff radius of the trial function
      trialrcut   = projections(indexproj)%rcut

!!     For debugging
!      write(6,'(a,3i5,4f12.5)')' Node, iproj, indexproj = ',  &
!                  Node, iproj, indexproj, trialcenter, trialrcut 
!!     End debugging

!     Find the atomic orbitals that ovelap with our radial orbital
!     centered at trialcenter and with range trialrcut
      item => get_overlapping_orbitals( scell, rmaxo, trialrcut, na_s, &
 &                                      xa, trialcenter )

OrbitalQueue:                                                        &
      do while (associated(item))
        r12 = trialcenter - item%center
        globalindexorbital = orb_gindex( item%specie, item%specieindex )
        call new_matel('S',                & ! Compute the overlap
 &                     globalindexorbital, & ! Between orbital with globalinde
 &                     globalindexproj,    & ! And projector with globalindex
 &                     r12,                & 
 &                     overlap,            & 
 &                     gradient )

        phase = -1.0_dp * dot_product( kvector, item%center )
        exponential = exp( iu * phase )

!!       For debugging
!        write(6,'(/,a,2i5)')               &
! &        'amn: Node, globalindex   = ', Node, item%globalindex
!        write(6,'(a,2i5)')                 &
! &        'amn: Node, globalindexorb= ', Node, globalindexorbital
!        write(6,'(a,2i5)')                 &
! &        'amn: Node, globalindexpro= ', Node, globalindexproj
!        write(6,'(a,i5,3f12.5)')           &
! &        'amn: Node, r12           = ', Node, r12
!        write(6,'(a,i5,f12.5)')            &
! &        'amn: Node, overlap       = ', Node, overlap
!        write(6,'(a,i5,3f12.5)')           &
! &        'amn: Node, gradient      = ', Node, gradient
!!       End debugging

!       Loop on all the nodes:
!       This is required because a given Node, for instance Node = 0,
!       knows the coefficients of the wave function 
!       for all the atomic orbitals (mu = 1, nuotot), 
!       but only for no_l bands

        do inode = 0, Nodes-1
        
!         Compute the number of occupied bands stored on Node inode 
          call GetNodeOrbs( no_u, inode, Nodes, nbands_loc )

!         Copy the coefficients for inode to a temporal variable
          if( Node .eq. inode ) then
            do iband = 1, nbands_loc
              do imu = 1, no_u
                psitmp(imu,iband) = psiloc(imu,iband)
              enddo
            enddo
          endif

!         Broadcast the auxiliary matrix from node inode to all the other nodes
          call MPI_Bcast( psitmp(1,1), no_u*nbands_loc, &
 &                  MPI_double_complex, inode, MPI_Comm_World, MPIerror )

!         Loop over occupied bands
          do iband = 1, nbands_loc  

!         Identify the global index of the occupied band
            call LocalToGlobalOrb( iband, inode, Nodes, iband_global )

            if( iband_global .le. nbands ) then
              cstar = conjg( psitmp(item%globalindex,iband) )
              Amnmat(iband_global,indexproj,ik) =   & 
 &              Amnmat(iband_global,indexproj,ik) + &
 &              exponential * cstar * overlap
            endif
          enddo  ! End loop on the bands

        enddo    ! End loop on the nodes
  

        item => item%nextitem
      enddo OrbitalQueue

    enddo   ! Loop on projections on the local node

  enddo kpoints

#ifdef MPI
! Allocate workspace array for global reduction
  nullify( auxloc )
  call re_alloc( auxloc, 1, nbands, 1, numproj, 1, numkpoints,     &
 &               name='auxloc', routine='Amn' )
! Global reduction of auxloc matrix
  auxloc(:,:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
  call MPI_AllReduce( Amnmat(1,1,1), auxloc(1,1,1), nbands*numproj*numkpoints, &
 &                    MPI_double_complex,MPI_sum,MPI_Comm_World,MPIerror )
  Amnmat(:,:,:) = auxloc(:,:,:)
#endif

  deallocate ( projector_gindex )

! Write the Amn overlap matrices in a file, in the format required
! by Wannier90
  if(IOnode) call writeamn( ispin )

! End time counter
  call timer( 'amn', 2 )

  return

  contains

  function get_overlapping_orbitals( latvec, atomrcut, trialrcut, numatoms, &
 &                                   atomcoords, trialcenter )              &
 &  result(firstitem)

    use precision,          only: dp            ! Real double precision type
    use sys,                only: die                 
    use neighbour,          only: maxnna        ! Maximum number of neighbours
    use neighbour,          only: jan           ! Atom-index of neighbours
    use neighbour,          only: r2ij          ! Squared distances to neighbors
    use neighbour,          only: xij           ! Vector from a given atom
                                                !   to its neighbours
    use neighbour,          only: mneighb       ! Subroutine to compute the
                                                !   number of neighbours
    use siesta_geom,        only: isa           ! Species index of each atom
    use atomlist,           only: lasto         ! Position of last orbital 
                                                !   of each atom
    use atomlist,           only: iphorb        ! Orbital index of each  orbital
                                                !   in its atom
    use atomlist,           only: indxuo        ! Index of equivalent orbital  
                                                !   in "u" cell
    use atmfuncs,           only: rcut          ! Function that determines the
                                                !   cutoff radius of a given
                                                !   orbital of a given specie

!   This subroutine yields a list of basis orbitals that overlap 
!   with a given trial orb.
!   This subroutine bridges siesta's mneighb() with our purposes.

    implicit none
  
!   Passing variables
    real(dp),dimension(3,3),intent(in) :: latvec      ! Lattice vectors of the 
                                                      !   supercell in 
                                                      !   real space.
    real(dp)               ,intent(in) :: atomrcut    ! Maximum cutoff radius in
                                                      !   the atomic orbital 
                                                      !   basis
    integer                ,intent(in) :: numatoms    ! Number of atoms in 
                                                      !   the supercell
    real(dp),dimension(3,numatoms),intent(in) :: atomcoords  ! Atomic positions 
    real(dp),dimension(3)  ,intent(in) :: trialcenter ! Center of the 
                                                      !  trial function 
    real(dp)               ,intent(in) :: trialrcut   ! Cutoff of the 
                                                      !  trial function 

!   Passing variables
    real(dp), dimension(:,:), allocatable,save :: coords 
!   A new array containing the coordinates of the na_s atoms in the unit cell
!   plus the position of the trial function will be required to call mneighb.
!   This new array, coords, will have (3,na_s+1) dimensions

    integer          :: jneig   ! Counter for the loop on neighbours
    integer          :: nneig   ! Number of neighbors of a given trial projector
    integer          :: atom    ! Atomic index of the neighbour
    integer          :: specie  ! Atomic species of the neighbour
    integer          :: iorb    ! Counter for loop on neighbour orbitals
    integer          :: joa     ! Index for the atomic orbital within 
                                !    a given atom
    real(dp)         :: radius  ! Radius that determine the scope of the search
    real(dp)         :: rij     ! Squared distance to neighbours

    type(orbitallinkedlist), pointer     :: firstitem
    type(orbitallinkedlist), pointer     :: newitem

! Initialize mneighb
    if( .not. allocated(coords) ) then !set up static "save" variables
      allocate( coords(3,numatoms+1) )
      coords(1:3,1:numatoms) = atomcoords(1:3,1:numatoms)
    endif

    coords(:,numatoms+1) = trialcenter(:)
    radius = trialrcut + atomrcut
    call mneighb( latvec, radius, numatoms+1, coords, 0 , 0, nneig )

!
! Look for atoms that overlap with the trial orbital
!
    call mneighb( latvec, radius, numatoms+1, coords, numatoms+1, 0, nneig )
    if (nneig.gt.maxnna)                                             &
 &    call die("amn: insufficient array shapes; see mneighb(..)")

!
! Prepare output list 
!

    firstitem => null()
    do jneig = 1, nneig
      atom = jan(jneig)
      if ( atom .gt. numatoms ) cycle !it is a trial projection function
      specie = isa(atom)              !it is an atom
      rij = dsqrt( r2ij(jneig) )
      do iorb = lasto( atom-1 )+1, lasto( atom )
        joa = iphorb( iorb )
        if ( (rcut(specie,joa) + trialrcut) .gt. rij ) then
          allocate(newitem)
          newitem%center(1:3) = xij(1:3,jneig) + trialcenter(1:3)
          newitem%specie      = specie
          newitem%specieindex = joa
          newitem%globalindex = indxuo( iorb )
          newitem%nextitem    => firstitem
          firstitem => newitem
        endif
      enddo
    enddo    ! End loop on number of neighbours
  end function get_overlapping_orbitals
end subroutine amn

