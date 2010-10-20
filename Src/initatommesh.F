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
!***********************************************************************
!
! subroutine InitAtomMesh
! ----------------------------------------------------------------------
! Initialises the atom information relating to the mesh for a given data
! distribution. It creates and initialise ipa. This arrays is allocated
! inside a meshDisType structure and can be accesed by pointers of the
! mesh module.
!
! Only necesary for those data distributions who needs to use the
! extended mesh.
!
! ----------------------------------------------------------------------
! Input :
! ----------------------------------------------------------------------
! integer  distr        : Used data distribution
! integer  na           : number of atoms
! integer  xa(3,na)     : Atomic coordinates
! ----------------------------------------------------------------------
! Output :
! ----------------------------------------------------------------------
! All output quantities are in the mesh and moreMeshSubs modules.
! The ipa array is stored in moreMeshSubs but accessed using the mesh
! module pointer.
!
!***********************************************************************
      subroutine InitAtomMesh( distr, na, xa )
      use precision, only: dp
      use mesh, only: ipa      ! atoms' mesh points within my Extended mesh
      use mesh, only: dxa      ! atom's position relative to mesh point
      use mesh, only: cmesh    ! Mesh cell vectors
      use mesh, only: meshLim  ! My processor's box of mesh points
      use mesh, only: ne       ! Points in rmax along each lat. vector
      use mesh, only: nmsc     ! Mesh points in each supercell vector
      use mesh, only: rcmesh   ! Reciprocal mesh cell vectors
      use alloc,only: re_alloc ! (Re)allocation routine
      use moreMeshSubs, only : allocIpaDistr
      implicit none
!     Input parameters
      integer,       intent(in) :: na        ! Number of atoms in supercell
      real(dp),      intent(in) :: xa(3,na)  ! Atomic coordinates
      integer,       intent(in) :: distr     ! Used data distribution
!     Local variables
      character(len=*),parameter:: myName = 'InitAtomMesh '
      integer,             save :: lastna = 0  ! Number of atoms on last call
      integer                   :: ia, ixa(3), jsc(3), jxa(3), nem(3),
     &                             j1, j2, j3, myBox(2,3), myExtBox(2,3)
      real(dp)                  :: cxa(3)

#ifdef DEBUG
      call write_debug( '      PRE InitAtomMesh' )
#endif
!     (Re)allocate atom-mesh arrays
      call allocIpaDistr( distr, na )
      if (na.ne.lastna) then
        call re_alloc( dxa, 1,3, 1,na, myName//'dxa' )
      endif
      lastna = na

      myBox(:,:) = meshLim(:,:) - 1

!     Add 'wings extensions' to myBox of mesh points, containing all points 
!     within the orbital spheres that may intersect myBox. Thus, wings must
!     be one diameter wide.
      myExtBox(1,:) = myBox(1,:) - 2*ne(:)
      myExtBox(2,:) = myBox(2,:) + 2*ne(:)

!     Find number of extended-box points.
      nem(:) = myExtBox(2,:) - myExtBox(1,:) + 1

!     Find atomic positions relative to mesh
      do ia = 1,na

!       Find atomic coordinates in mesh basis vectors
        cxa(:) = matmul( xa(:,ia), rcmesh(:,:) )
        ! fix: Make sure we consider equivalent positions
        !      inside the cell, even if some atoms are outside.
        cxa(:) = modulo( cxa(:), real(nmsc(:),dp) )

!       Find indexes of mesh cell in which atom is
        ixa(:) = floor(cxa(:))

!       Find atom position within mesh cell
        cxa(:) = cxa(:) - ixa(:)                  ! Mesh coordinates
        dxa(:,ia) = matmul( cmesh(:,:), cxa(:) )  ! Cartesian coords

!       Find atom's mesh index within myExtBox, but only if the
!       (periodic) atomic sphere (of width ne) can intersect myBox.
!       Otherwise, make ipa=0. nem is the width of myExtBox.
!       nmsc is the width of the supercell.
!       Notice that, even if several periodic images of the same 
!       atom may fit within myExtBox, only one must be considered,
!       since the points within its sphere, once folded into myBox,
!       will be equivalent to those of any other periodic image.
        ipa(ia) = 0
        sc: do j3 = -1,1   ! Loop on neighbor supercells
            do j2 = -1,1
            do j1 = -1,1
              jsc(:) = (/j1,j2,j3/)
              jxa(:) = ixa(:) + jsc(:)*nmsc(:)
              if ( all(jxa(:)>=myBox(1,:)-ne(:)) .and. 
     &             all(jxa(:)<=myBox(2,:)+ne(:)) ) then
                jxa(:) = jxa(:) - myExtBox(1,:)
                ipa(ia) = 1 + jxa(1) + nem(1)*jxa(2)
     &                      + nem(1)*nem(2)*jxa(3)
                
                exit sc
              end if
            end do    ! j1
          end do      ! j2
        end do sc     ! j3
      enddo ! ia

#ifdef DEBUG
      call write_debug( '      POS InitAtomMesh' )
#endif
      end subroutine InitAtomMesh
