
subroutine Mmn( ispin ) 
!
!     In this subroutine we compute the overlaps between Bloch orbitals
!     at neighboring k points:
!
!     $M_{m n}^{(\vec{k}, \vec{b})} = 
!        \langle u_{m \vec{k}} \vbar u_{n \vec{k} + \vec{b} \rangle$
!
!     Eq. (27) of the paper by N. Marzari et al. 
!     Review of Modern Physics 84, 1419 (2012)
!
!     In the previous formula, only the periodic part of the wave functions
!     at neighbour k-points enter in the Equation.
!     We have to adapt this equation to the input provided by Siesta
!     (the coefficients of the whole wave function, not only of the periodic
!     part). 
!     This is done following Eq. (5) of the paper by 
!     D. Sanchez-Portal et al., Fundamental Physics for Ferroelectrics 
!     (AIP Conf. Proc. Vol 535) ed R. Cohen (Melville, AIP) pp 111-120 (2000).
!
!     OUTPUT: 
!     File called seedname.mmn, where the overlap matrices are written in 
!     the format required by Wannier90
!
!     File called seedname.eig, where the eigenvalues of the Hamiltonian
!     at the required k-points are written according to the format required
!     by Wannier90.
!
!     Implemented by J. Junquera and R. Korytar, July 2013
!

  use precision,          only: dp           ! Real double precision type
  use m_siesta2wannier90, only: numkpoints   ! Total number of k-points
                                             !   for which the overlap of the
                                             !   periodic part of the wavefunct
                                             !   with a neighbour k-point will
                                             !   be computed
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
                                             !   before excluding bands         
  use m_siesta2wannier90, only: numincbands  ! Number of bands for wannierizatio
                                             !   after excluding bands         
  use m_siesta2wannier90, only: numproj      ! Total number of projectors
  use m_siesta2wannier90, only: bvectorsfrac ! Vectors that connect each 
                                             !   mesh k-point to its 
                                             !   nearest neighbours.
  use m_siesta2wannier90, only: Mmnkb        ! Matrix of the overlaps of 
                                             !   periodic parts of Bloch waves.
                                             !   <u_{n,k}|u_{m,k+b}>
  use m_siesta2wannier90, only: Amnmat       ! Matrix of the overlaps of 
                                             !   trial projector funtions
                                             !   with Eigenstates of the
                                             !   Hamiltonian
  use m_siesta2wannier90, only: eo           ! Eigenvalues of the Hamiltonian 
                                             !   at the numkpoints introduced in
                                             !   kpointsfrac 
                                             !   First  index: band index
                                             !   Second index: k-point index
  use m_siesta2wannier90, only: isexcluded   ! Masks excluded bands
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
  use siesta_options, only: towriteamn       ! Write the Amn matrix for the
                                             !   interface with Wannier
  use siesta_options, only: towriteunk       ! Write the periodic part of
                                             !   the wave function at the point
                                             !   of a real space mesh
  use siesta_options, only: towriteeig       ! Write the eigenvalues for the
                                             !   interface with Wannier
! Notes about the use of psi:
! the array psi is allocated with a size equal to 2 * no_u * no_l
! 2 is for the complex nature of the coefficients (we have to store the
! real and the imaginary parts).
! no_u comes from the fact that for every band we have as many coefficients
! as basis orbitals in the unit cell.
! no_l is due to the fact that each node knows the coefficients of no_l bands.
! When running in serial, no_l = no_u.
! When running in parallel, no_l < no_u, and sum_{nodes} no_l = no_u.
! In mmn.F90, psi is a one dimensional array.
! 
! In diagpol, as in cdiag, this array has three dimensions,
! (2, nuotot=no_u, nuo=no_l).
! The first index refers to the real or imaginary part.
! The second index refers to the atomic orbital.
! The third index is the band index.
! 
! Since only a few bands might be used for wannierization, 
! and those bands might not be the numincbands lowest bands,
! a first reordering of the coefficients is done in the subroutine reordpsi.
! After calling it, the bands whose number ranges from 1 to numincbands(ispin)
! correspond to those bands that will be included for wannierization.
! The rest of the coefficients are set to zero.
! 
! Then, the array with the coefficients is reallocated in savepsi, 
! when we want to store it
! to save the coefficients for a given k-point.
! There, a new variable is defined, psiatk (psiprev in ksv.f), 
! with a dimension of 
! (2, nuo=no_l, nuotot=no_u).
! Again, the first index refers to the real or imaginary part.
! The second index refers to the atomic orbital.
! The third index is the band index.
!
! But now, while running in parallel, regarding psiatk each node knows the 
! coefficients of no_l orbitals of all the OCCUPIED BANDS.
! For the unoccupied bands the coefficients are set to zero.

  use alloc,              only: re_alloc     ! Reallocation routines
  use alloc,              only: de_alloc     ! Deallocation routines
  use parallel,           only: IOnode       ! Input/output node
  use m_digest_nnkp,      only: getdelkmatgenhandle
  use m_noccbands,        only: noccupied    ! Number of occupied bands for a 
                                             !   given spin direction
  use units,              only: eV           ! Conversion factor from Ry to eV
  use m_planewavematrixvar, only: delkmat    ! Matrix elements of a plane wave
                                             !   Only for one \vec{b} vector
  use m_planewavematrixvar, only: delkmatgen ! Matrix elements of a plane wave
                                             !   (array that contains the 
                                             !   matrices for all the \vec{b}
                                             !   vectors

! For debugging
  use parallel,           only: Node, Nodes
! End debugging

  implicit none

  integer,  intent(in) :: ispin              ! Spin component
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
  integer  :: nincbands      ! Number of occupied bands
  real(dp) :: kvector(3)     ! k-point vector for which the Overlap matrix 
                             !   between the periodic part of the 
                             !   wave functions will be computed
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
  real(dp), dimension(:), pointer :: psitmp  ! Temporal array used to 
                                             !   store the coefficients of the
                                             !   wave functions

! Start time counter
  call timer( 'Mmn', 1 )

! Allocate memory related with the dense matrices that will be used in
! the diagonalization routines.
! (Hamiltonian, Overlap and Coefficient vectors)
! These matrices are defined in the module dense matrix
  nhs  = 2 * no_u * no_l
  npsi = 2 * no_u * no_l

  call re_alloc( Haux, 1, nhs,  'Haux', 'Mmn' )
  call re_alloc( Saux, 1, nhs,  'Saux', 'Mmn' )
  call re_alloc( psi,  1, npsi, 'psi',  'Mmn' )

! Allocate memory related with the eigenvalues of the Hamiltonian (epsilon)
! and with a local variable where the coefficients of the eigenvector at the
! k-point will be stored
  nullify( epsilon, psiatk, psitmp )
  call re_alloc( epsilon, 1, no_u, 'epsilon', 'Mmn' )
  call re_alloc( psiatk,  1, npsi, 'psiatk',  'Mmn' )
  call re_alloc( psitmp,  1, npsi, 'psitmp',  'Mmn' )

! Allocate memory related with the overlap matrix between periodic parts
! of Bloch functions at neighbour k-points.
! These matrix will depend on four indices (see Eq. (27) of the review by
! N. Marzari et al., Review of Modern Physics 84, 1419 (2012):
! $M_{m n}^{(\vec{k}, \vec{b})} = 
!    \langle u_{m \vec{k} \vbar u_{n \vec{k} + \vec{b}} \rangle
! where $m$ and $n$ run between 1 and the number of occupied bands
! $\vec{k}$ runs over all the k-points in the first BZ where these matrices 
! will be computed and 
! $\vec{b}$ runs over all the neighbours of the k-point.
! These last two variables are read from the .nnkp file.
  nincbands = numincbands( ispin )
  nullify( Mmnkb )
  call re_alloc( Mmnkb,          &
 &               1, nincbands,   &
 &               1, nincbands,   &
 &               1, numkpoints,  &
 &               1, nncount,     &
 &               'Mmnkb',        &
 &               'Mmn' )

! Allocate memory related with the overlap matrix between the trial projection
! function and the Hamiltonian eigenstate.
! These matrices will depend on three indices (see paragraph after Eq. (16) 
! of the review by
! N. Marzari et al., Review of Modern Physics 84, 1419 (2012),
! or Eq. (1.8) of the Wannier Users Guide, Version 1.2
! $A_{m n}^{(\vec{k})} = 
!    \langle \psi_{m \vec{k}} \vbar g_{n} \rangle =
!    \sum_{\mu} \sum_{lcell} c_{\mu m}^{\ast} (\vec{k}) \times
!    exp^{i \vec{k} \cdot \left( \tau_{\mu} + \vec{R}_{lcell} \right)} \times
!    \langle \phi_{\mu} (\vec{r}-\tau_{\mu} - \vec{r}_{lcell} \vbar g_{n}\rangle
! where $m$ runs between 1 and the number of bands considered for wannierization
! $n$ runs between 1 and the number of projection functions, 
! $\vec{k}$ runs over all the k-points in the first BZ where these matrices 
! will be computed 

  nullify( Amnmat )
  call re_alloc( Amnmat,         &
 &               1, nincbands,   &
 &               1, numproj,     &
 &               1, numkpoints,  &
 &               'Amnmat',       &
 &               'Mmn'           )

! Allocate memory related with the eigenvalues of the Hamiltonian
  nullify( eo )
  call re_alloc( eo,             &
 &               1, nincbands,   &
 &               1, numkpoints,  &
 &               'eo',           &
 &               'Mmn' )

! Initialise psi
  do io = 1, npsi
    psi(io)     = 0.0_dp
    psiatk(io)  = 0.0_dp
    psitmp(io)  = 0.0_dp
  enddo

  do ik = 1, numkpoints
!   Compute the wave vector in bohr^-1 for every vector in the list
!   (done in the subroutine getkvector).
!   Remember that kpointsfrac are read from the .nnkp file in reduced units, 
!   so we have to multiply then by the reciprocal lattice vector.
    call getkvector( kpointsfrac(:,ik), kvector )

!!   For debugging
!    if( IOnode ) then
!      write(6,'(/,a,3f12.5)')         &
! &      'mmn: kvector = ', kvector
!!      do io = 1, maxnh
!!        write(6,'(i5,2f12.5)')io, H(io,ispin), S(io)
!!      enddo 
!    endif
!!   End debugging

!   Diagonalize the Hamiltonian for the k-point.
!   Here, we obtain $\psi_{n} (\vec{k})$, where n runs between 1 and no_l
    call diagpol( ispin, nspin, no_l, no_s, no_u,                             &
 &                maxnh, numh, listhptr, listh, H, S, xijo, indxuo, kvector,  &
 &                epsilon, psitmp, 2, Haux, Saux )

!!   For debugging
!!   NOTE OF CAUTION: beware when comparing the coefficients of the
!!   wave function obtained in differente machines, specially for 
!!   degenerate states. 
!!   There is a phase that is arbitrary and might change from one machine
!!   to the other. Also, any linear combination of eigenvectors with 
!!   the same eigenvalue is also a solution of the Hamiltonian,
!!   and the coefficients of the linear combination might be different.
!    if( Node .eq. 1 ) then
!!    if( IOnode ) then
!      do iband = 1, no_u
!        write(6,'(i5,f12.5)') iband, epsilon(iband)/eV
!      enddo
!      iuo = 0
!      do iband = 1, no_l
!        write(6,'(a,i5,f12.5)')         &
! &        'mmn: iband, epsilon = ', iband, epsilon(iband)/eV
!        do io = 1, no_u
!          write(6,'(a,i5,2f12.5)')         &
! &          'mmn: io, coeff = ', io, psitmp(iuo+1), psitmp(iuo+2)
!          iuo = iuo + 2
!        enddo 
!      enddo 
!    endif
!!   End debugging

!   Store the eigenvalues, while skipping the excluded bands
    eo(1:numincbands(ispin),ik) = pack(epsilon/eV,.not.isexcluded)

!   Keep only the coefficients of the included bands for wannierization,
!   while skipping the excluded eigenvectors.
!   The coefficients of the included eigenvectors will be stored in psi.
!   Bands whose band index ranges from 1 to numbands(ispin) correspond
!   to the bands included for wannierization.
    call reordpsi( psi, psitmp, no_l, no_u, numbands(ispin) )

!!   For debugging
!    if( IOnode ) then
!      iuo = 0
!      do iband = 1, no_l
!        write(6,'(a,i5,f12.5)')         &
! &        'mmn: iband, eo = ', iband, eo(iband,ik)
!        do io = 1, no_u
!          write(6,'(a,i5,2f12.5)')         &
! &          'mmn: io, psi   = ', io, psi(iuo+1), psi(iuo+2)
!          iuo = iuo + 2
!        enddo 
!      enddo 
!    endif
!!   End debugging

!   Compute the overlaps between Bloch states onto trial localized orbitals
    if( towriteamn ) call amn( ik, kvector, nincbands, npsi, psi )

!   Compute the overlaps between Bloch states onto trial localized orbitals
    if( towriteunk ) call writeunk( ispin, ik, kvector, npsi, psi )

!   Store the wavefunction at the k-point. 
!   Take into account the previous comments about how psi is stored.
!   Only the coefficients for the bands considered for wannierization are kept.
!   The coefficients for the non-considered bands are set to 0.
    call savepsi( psiatk, psi, no_l, no_u, nincbands )

!!   For debugging
!    if( IOnode ) then
!      do iband = 1, numincbands(ispin)
!        write(6,*) ik, iband, eo(iband,ik)
!      enddo
!      do iband = 1, no_u
!        write(6,*) ik, iband, epsilon(iband)/eV
!      enddo
!    endif
!    if( Node .eq. 1 ) then
!!    if( IOnode ) then
!      iuo = 0
!      do iband = 1, no_u
!        do io = 1, no_l
!          write(6,'(a,2i5,2f12.5)')         &
! &          'savepsi: io, coeff = ', io, iband, psiatk(iuo+1), psiatk(iuo+2)
!          iuo = iuo + 2
!        enddo 
!      enddo 
!    endif
!!   End debugging

!   Loop on the neighbour k-points for a given k.
    do inn = 1, nncount
!     Initialise psi
      do io = 1, npsi
        psi(io)     = 0.0_dp
        psitmp(io)  = 0.0_dp
      enddo 

!     Get the coordinates of the neighbor k-point.
      indexneig = nnlist(ik,inn)
      call getkvector( kpointsfrac(:,indexneig), kvectorneig )
!!     For debugging
!      if( IOnode ) then
!        write(6,'(a,3f12.5)')         &
! &        'mmn: kvectorneig = ', kvectorneig
!      endif
!!     End debugging

!     Find the coefficients of the wave function for the neighbour k-point
!     Here, we obtain $\psi_{m} (\vec{k} + \vec{b})$, 
!     where m runs between 1 and no_l
      call diagpol( ispin, nspin, no_l, no_s, no_u,                   &
 &                  maxnh, numh, listhptr, listh, H, S, xijo, indxuo, &
 &                  kvectorneig, epsilon, psitmp, 2, Haux, Saux )

!     Keep only the coefficients of the included bands for wannierization,
!     while skipping the excluded eigenvectors.
!     The coefficients of the included eigenvectors will be stored in psi.
!     Bands whose band index ranges from 1 to numincbands(ispin) correspond
!     to the bands included for wannierization, while the rest of the
!     coefficients are set to zero.
      call reordpsi( psi, psitmp, no_l, no_u, numbands(ispin) )

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
!     c_{i \mu} (k+b) = c_{i\mu}(k+b-G) exp(-iG\cdot r_\mu)
!     where -G brings the k+b to the first BZ
      if ( fold .gt. 0 ) then
        foldfrac(:) = nnfolding(:,ik,inn) * 1.0_dp
        call getkvector( foldfrac, gfold )  
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

      delkmat(1:maxnh) = delkmatgen(handle,1:maxnh)

!! For debugging
!      if ( IOnode ) then
!        write(6,'(a,2i5)')         &
! &        'mmn: inn, handle  = ', inn, handle
!        write(6,'(a,i5,3f12.5)')   &
! &        'mmn: inn, bvector = ', inn, bvector
!!        do io = 1, maxnh
!!          write(6,'(a,i5,2f12.5)')         &
!! &          'mmn: io, delkmat  = ', io, delkmat(io)
!!        enddo
!      endif 
!! End debugging

!     Initialize Mmnkb
      Mmnkb(:,:,ik,inn) = cmplx(0.0_dp, 0.0_dp, kind=dp)

!     Compute Mmnkb, following Eq. (5) of the paper by 
!     D. Sanchez-Portal et al., Fundamental Physics for Ferroelectrics 
!     (AIP Conf. Proc. Vol 535) ed R. Cohen (Melville, AIP) pp 111-120 (2000).
      call overkkneig( kvector, bvector, no_l, no_u, psiatk, psi,   &
 &                     maxnh, delkmat, nincbands, Mmnkb(:,:,ik,inn) )

    enddo ! End loop on neighbour k-points (inn)

  enddo  ! End loop on the number of k-points (ik)


! Write the Mmn overlap matrices in a file, in the format required
! by Wannier90
  if( IOnode ) call writemmn( ispin )

! Write the Amn overlap matrices in a file, in the format required
! by Wannier90
  if( IOnode .and. towriteamn ) call writeamn( ispin )

! Write the eigenvalues in a file, in the format required
! by Wannier90
  if( IOnode .and. towriteeig ) call writeeig( ispin )

! Deallocate some of the arrays
  call de_alloc( Haux,    'Haux',    'mmn' )
  call de_alloc( Saux,    'Saux',    'mmn' )
  call de_alloc( psi,     'psi',     'mmn' )
  call de_alloc( epsilon, 'epsilon', 'mmn' )
  call de_alloc( psiatk,  'psiatk',  'mmn' )
  call de_alloc( psitmp,  'psitmp',  'mmn' )

! End time counter
  call timer( 'Mmn', 2 )
 
end subroutine Mmn
