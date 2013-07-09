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
subroutine amn( ispin )

  use precision,          only: dp                  ! Real double precision type
  use parallel,           only: Nodes               ! Total number of Nodes
  use parallel,           only: Node                ! Local Node
  use atomlist,           only: rmaxo               ! Max. cutoff atomic orbital
  use siesta_geom,        only: na_u                ! Number of atoms in the
                                                    !   unit cell
  use siesta_geom,        only: xa                  ! Atomic positions
  use m_siesta2wannier90, only: latvec              ! Lattice vectors in real 
                                                    !   space
  use m_siesta2wannier90, only: numproj             ! Total number of projectors
  use m_siesta2wannier90, only: projections         ! Trial projection functions
  use neighbour,          only: xij

  implicit none

! Passed arguments
  integer,  intent(in) :: ispin                   ! Spin component

! Local variables
  integer  :: numproj_l       ! Number of projections to be computed locally
                              !   in this node

  integer  :: minpernode      ! Minimum number of projections per node
  integer  :: remainder       ! Remainder of unassigned projections
  integer  :: iproj           ! Counter for loops on projetions
  integer  :: globalindexproj ! Global index of the projector 
                              !   This index runs from 1 to the total number
                              !   of projections

  real(dp) :: trialcentre(3)  ! Position where the trial function is centered
                              !   (in Bohr)
  real(dp) :: trialrcut       ! Cutoff radius of the trial function

  external :: timer

! Start time counter
  call timer( 'amn', 1 )

! Calculate the minimum number of projections per node
  minpernode = numproj / Nodes

! Find the remainder of unassigned projections
  remainder = numproj - minpernode * Nodes 

! Workout the local number of orbitals
  numproj_l = minpernode 
  if (Node .lt. remainder) numproj_l = numproj_l + 1

! Loop on the projections that will be computed in the local node
  do iproj = 1, numproj_l 
!   Identify the global index 
    if ( remainder .eq. 0 ) then
      globalindexproj = minpernode * Node + iproj
    else 
      if ( Node .lt. remainder ) then 
         globalindexproj = ( minpernode + 1 ) * Node + iproj 
      else
         globalindexproj = ( minpernode + 1 ) * remainder + &
 &                         ( Node - remainder )+ 1
      endif
    endif

!   Find where the trial function is centered
    trialcentre = projections(globalindexproj)%center
    
!   Find the cutoff radius of the trial function
    trialrcut   = projections(globalindexproj)%rcut

!!   For debugging
!    write(6,'(a,3i5,4f12.5)')' Node, iproj, globalindexproj = ',  &
!                Node, iproj, globalindexproj, trialcentre, trialrcut 
!!   End debugging

!   FIND THE CORRESPONDING INDEX IN THE LIST OF RADIAL FUNCTIONS TO CALL MATEL

!   Find the atomic orbitals that ovelap with our radial orbital
!   centered at trialcentre and with range trialrcut
    call get_overlapping_orbitals( latvec, rmaxo, trialrcut, na_u, &
 &                                 xa, trialcentre )

!   For debugging
    write(6,'(/,a,2i5)')         &
 &      'amn: Node, nna     = ', Node, xij
!   End debugging



  enddo   ! Loop on projections on the local node

! End time counter
  call timer( 'amn', 2 )

end subroutine amn

subroutine get_overlapping_orbitals( latvec, atomrcut, trialrcut, numatoms, &
 &                                   atomcoords, trialcentre )

  use precision,          only: dp                  ! Real double precision type
  use sys,                only: die                 
  use neighbour

! This subroutine yields a list of basis orbitals that overlap 
! with a given trial orb.
! This subroutine bridges siesta's mneighb() with our purposes.

  implicit none
  
! Passing variables
  real(dp),dimension(3,3),intent(in) :: latvec      ! Lattice vectors in 
                                                    !   real space.
                                                    !   Ordered as read from 
                                                    !   the nnkp file: 
                                                    !   First  index: vector
                                                    !   Second index: component
  real(dp)               ,intent(in) :: atomrcut    ! Maximum cutoff radius in
                                                    !   the atomic orbital basis
  integer                ,intent(in) :: numatoms    ! Number of atoms in 
                                                    !   the unit cell
  real(dp),dimension(:,:),intent(in) :: atomcoords  ! Atomic positions (3,na_u)
  real(dp),dimension(3)  ,intent(in) :: trialcentre ! Center of the trial funct 
  real(dp)               ,intent(in) :: trialrcut   ! Cutoff of the trial funct

! Passing variables
  real(dp), dimension(:,:), allocatable,save :: coords 
! A new array containing the coordinates of the na_u atoms in the unit cell
! plus the position of the trial function will be required to call mneighb.
! This new array, coords, will have (3,na_u+1) dimensions

  real(dp), dimension(3,3), save             :: latvect   
! neighb requires the introduction of the lattice vector in real space
! in the transpose format introduced in nnkp  

! Iterators: neighbors,orbitals...
  integer                            :: nneig,orb,specie,atom,jneig
  integer                            :: joa
! Scope of the search
  real(dp)                           :: radius

!
! Initialize mneighb
!
  if( .not. allocated(coords) ) then !set up static "save" variables
    allocate( coords(3,numatoms+1) )
    coords(1:3,1:numatoms) = atomcoords(1:3,1:numatoms)
    latvect = transpose(latvec)
  endif

  coords(:,numatoms+1) = trialcentre(:)
  radius = trialrcut + atomrcut
  call mneighb( latvect, radius, numAtoms+1, coords, 0 , 0, nneig )

!
! Look for atoms that overlap with the trial orbital
!
  call mneighb(latVecT,radius,numAtoms+1,coords,numAtoms+1,0,nneig)
  if (nneig.gt.maxnna)                                             &
 &  call die("amn: insufficient array shapes; see mneighb(..)")

end subroutine get_overlapping_orbitals

