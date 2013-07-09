
subroutine Mmn( ispin ) 

  use precision,          only: dp           ! Real double precision type
  use m_siesta2wannier90, only: numkpoints   ! Total number of k-points
  use m_siesta2wannier90, only: kpointsfrac  ! List of k points relative to the 
                                             !   reciprocal lattice vectors.
                                             !   First  index: component
                                             !   Second index: k-point index 
                                             !      in the list
  use m_siesta2wannier90, only: nncount      ! Number of nearest k-pnt neighbors
  use m_siesta2wannier90, only: nnlist       ! nnlist(ikp,inn) is the index of
                                             !   the inn-neighbour of ikp-point
                                             !   in the Monkhorst-Pack grid 
                                             !   folded to first Brillouin zone
  use m_siesta2wannier90, only: nnfolding    ! nnfolding(i,ikp,inn) is the 
                                             !   i-component of the reciprocal 
                                             !   lattice vector 
                                             !   (in reduced units) that brings
                                             !   the inn-neighbour specified 
                                             !   in nnlist
                                             !   (which is in the first BZ) to 
                                             !   the actual \vec{k} + \vec{b} 
                                             !   that we need.
  use m_siesta2wannier90, only: numbands     ! Number of bands for wannierizatio
  use m_siesta2wannier90, only: bvectorsfrac ! Vectors that connect each 
                                             !   mesh k-point to its 
                                             !   nearest neighbours.
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
  use atomlist,           only: iaorb        ! Atomic index of each orbital
  use siesta_geom,        only: xa           ! Atomic positions
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
  use m_digest_nnkp,      only: getdelkmatgenhandle

! For debugging
  use parallel,           only: Node, Nodes, IOnode
  use sys,                only: die
! End debugging

  implicit none

  integer,  intent(in) :: ispin                   ! Spin component
! Internal variables

  integer  :: ik             ! Counter for the k-point loop
  integer  :: io             ! Counter for the orbital loop
  integer  :: jo             ! Counter for the orbital loop
  integer  :: iuo            ! Counter for the orbital loop
  integer  :: iband          ! Counter for the bands loop
  integer  :: inn            ! Counter for the k-point neighbor loop
  integer  :: nhs            ! Variable to dimension the Hamiltonian and Overlap
  integer  :: npsi           ! Variable to dimension the coefficient vector
  integer  :: indexneig      ! Index of the neighbour k-point in the list
                             !   of k-points
  integer  :: fold           ! Should we fold the coefficients of the 
                             !   wavefunction at the neighbour k-point
                             !   out of the first Brillouin zone?
  integer  :: ia             ! Index of the atom to which an atomic orbital
                             !   belongs
  integer  :: handle         ! Given a k-point vector and a neighbor, 
                             !   separated by a vector \vec{b},
                             !   it gives the position of the delkmatgen 
                             !   array where exp^(i \vec{b} \cdot \vec{r})
                             !   will be stored
  real(dp) :: kvector(3)     ! Wave vector where the Overlap matrix between
                             ! the periodic part of the wave functions will be
                             ! computed
  real(dp) :: kvectorneig(3) ! Wave vector of the neighbor k-point
  real(dp) :: gfold(3)       ! Reciprocal lattice vector that brings
                             !   the inn-neighbour specified in nnlist
                             !   (which is in the first BZ) to the
                             !   actual \vec{k} + \vec{b} that we need.
  real(dp) :: bvectoraux(3)  ! Auxiliary vector
  real(dp) :: bvector(3)     ! Vector that connects a given k-point with its 
                             !   neighbour (in Bohr^-1)
  real(dp) :: gxij           ! Dot product of the reciprocal lattice vector
                             !   gfold with an atomic position
  real(dp) :: foldfrac(3)    ! Auxiliar vector to compute the folding vector
  complex(dp) :: eigx        ! Exponential exp^( i * gxij )
  complex(dp) :: psiaux      ! Auxiliar complex values

  real(dp), dimension(:), pointer :: epsilon ! Eigenvalues of the Hamiltonian
  real(dp), dimension(:), pointer :: psiatk  ! Coefficients of the eigenvector
                                             !   at the k-point

! Start time counter
  call timer( 'Mmn', 1 )

! Allocate memory related with the dense matrices
! (Hamiltonian, Overlap and Coefficient vectors)
! These matrices are defined in the module dense matrix
  nhs  = 2 * no_u * no_l
  npsi = 2 * no_u * no_l

  call re_alloc( Haux, 1, nhs,  'Haux', 'Mmn' )
  call re_alloc( Saux, 1, nhs,  'Saux', 'Mmn' )
  call re_alloc( psi,  1, npsi, 'psi',  'Mmn' )

! Allocate memory related with the eigenvalues of the Hamiltonian (epsilon)
! and with a local variable where the coefficients of the eigenvector for the
! k-point will be stored
  nullify( epsilon, psiatk )
  call re_alloc( epsilon,      1, no_u, 'epsilon',      'Mmn' )
  call re_alloc( psiatk,       1, npsi, 'psiatk',       'Mmn' )

! Initialise psi
  do io = 1, npsi
    psi(io)     = 0.0_dp
    psiatk(io)  = 0.0_dp
  enddo

  do ik = 1, numkpoints
!   Compute the wave vector in bohr^-1 for every vector in the list
!   (done in the subroutine getkvector).
!   Remember that kpointsfrac are read in reduced units, so we have
!   to multiply then by the reciprocal lattice vector.
    call getkvector( kpointsfrac(:,ik), kvector )

!   For debugging
    if( IOnode ) then
      write(6,'(/,a,3f12.5)')         &
 &      'mmn: kvector = ', kvector
!      do io = 1, maxnh
!        write(6,'(i5,2f12.5)')io, H(io,ispin), S(io)
!      enddo 
    endif
!   End debugging

    call diagpol( ispin, nspin, no_l, no_s, no_u,                             &
 &                maxnh, numh, listhptr, listh, H, S, xijo, indxuo, kVector,  &
 &                epsilon, psi, 2, Haux, Saux )

!!   For debugging
!!   NOTE OF CAUTION: beware when comparing the coefficients of the
!!   wave function obtained in differente machines, specially for 
!!   degenerate states. 
!!   There is a phase that is arbitrary and might change from one machine
!!   to the other. Also, any linear combination of eigenvectors with 
!!   the same eigenvalue is also a solution of the Hamiltonian,
!!   and the coefficients of the linear combination might be different.
!!    if( Node .eq. 1 ) then
!    if( IOnode ) then
!      iuo = 0
!      do iband = 1, no_l
!        write(6,'(a,i5,f12.5)')         &
! &        'mmn: iband, epsilon = ', iband, epsilon(iband)
!        do io = 1, no_u
!          write(6,'(a,i5,2f12.5)')         &
! &          'mmn: io, coeff = ', io, psi(iuo+1), psi(iuo+2)
!          iuo = iuo + 2
!        enddo 
!      enddo 
!    endif
!!   End debugging

!   Store the wavefunction at the k-point. 
    call savepsi( psiatk, psi, no_l, no_u, numbands(ispin) )

    do inn = 1, nncount

! Initialise psi
      do io = 1, npsi
        psi(io)     = 0.0_dp
      enddo 

!     Get the coordinates of the neighbor k-point.
      indexneig = nnlist(ik,inn)
      call getkvector( kpointsfrac(:,indexneig), kvectorneig )
!   For debugging
      if( IOnode ) then
        write(6,'(a,3f12.5)')         &
 &        'mmn: kvectorneig = ', kvectorneig
      endif
!   End debugging

!     Find the coefficients of the wave function for the neighbour k-point
      call diagpol( ispin, nspin, no_l, no_s, no_u,                   &
 &                  maxnh, numh, listhptr, listh, H, S, xijo, indxuo, &
 &                  kvectorneig, epsilon, psi, 2, Haux, Saux )

      
!     The neighbour k-point, as specified in the nnlist, 
!     is always in the first Brillouin zone.
!     To find the actual k-point, we might have to add a vector of the
!     reciprocal lattice, as specified in the nnfolding matrix
      fold =  nnfolding(1,ik,inn)**2 +  &
 &            nnfolding(2,ik,inn)**2 +  &
 &            nnfolding(3,ik,inn)**2

!     The coeffs of the wave function, as obtained in diagpol, 
!     are for a wave vector in the first Brillouin zone. 
!     If the actual neighbour is out,
!     we apply a simple transformation that holds in the periodic gauge
!     c_{i\mu}(k+b) = c_{i\mu}(k+b-G)exp(-iG\cdot r_\mu)
!     where -G brings the k+b to the first BZ
      if ( fold .gt. 0 ) then
        foldfrac(:) = nnfolding(:,ik,inn) * 1.0_dp
        call getkvector( foldfrac, gfold )  
!        do iband = 1, numincbands(ispin)
!!       For debugging
!        if ( IOnode )
!          write(6,'(a,i5,3i5,3f12.5)')         &
! &          'mmn: inn, gfold = ', inn, nnfolding(:,ik,inn), gfold
!        endif
!!       End debugging
        iuo = 0
        do iband = 1, no_l
          do jo = 1, no_u
!           Localize the position where the atom is centered
            ia = iaorb(jo)
!           Compute exp( i G tau_{mu} )
            gxij = dot_product( gfold,xa(:,ia) )
            eigx = cmplx( dcos(gxij), dsin(gxij), kind=dp )
!           Multiply the coefficient times the gauge
            psiaux = cmplx(psi(iuo+1), psi(iuo+2), kind=dp) * conjg( eigx )
            psi(iuo+1) = real(psiaux)  * 1.0_dp 
            psi(iuo+2) = aimag(psiaux) * 1.0_dp
            iuo = iuo + 2
          enddo
        enddo
      endif

!!     For debugging
!!      if( Node .eq. 1 ) then
!      if( IOnode ) then
!        write(6,'(a,2i5)')         &
! &        'mmn: inn, fold = ', inn, fold
!        iuo = 0
!        do iband = 1, no_l
!          write(6,'(a,i5,f12.5)')         &
! &          'mmn: iband, epsilon = ', iband, epsilon(iband)
!          do io = 1, no_u
!            write(6,'(a,i5,2f12.5)')         &
! &            'mmn: io, coeff = ', io, psi(iuo+1), psi(iuo+2)
!            iuo = iuo + 2
!          enddo 
!        enddo 
!      endif
!!     End debugging

!     Now we have to determine the position in the delkmatgen array,
!     where exp^(i \vec{b} \cdot \vec{r}) is stored,
!     where \vec{b} is the vector connecting kvector and kvectorneig.
      bvectoraux(:) = kpointsfrac(:,nnlist(ik,inn)) +               &
 &                     nnFolding(:,ik,inn) - kPointsFrac(:,ik)
      handle = getdelkmatgenhandle( bvectoraux, nncount, bvectorsfrac )
      call getkvector( bvectoraux, bvector )
!! For debugging
!      if ( IOnode ) then
!        write(6,'(a,2i5)')         &
! &        'mmn: inn, handle  = ', inn, handle
!        write(6,'(a,i5,3f12.5)')   &
! &        'mmn: inn, bvector = ', inn, bvector
!      endif 
!! End debugging

    enddo ! End loop on neighbour k-points

  enddo  ! End loop on the number of k-points

! End time counter
  call timer( 'Mmn', 2 )
 
end subroutine Mmn
