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
! subroutine setupExtMesh
! ----------------------------------------------------------------------
! Setup the extended mesh variables for a given data distribution.
! It creates and initialise indexp, xdop and idop. These arrays are
! allocated inside a meshDisType structure and can be accesed by pointers
! of the mesh module.
!
! Only necesary for those data distributions who needs to use the
! extended mesh.
! ----------------------------------------------------------------------
! Input :
! ----------------------------------------------------------------------
! integer  distr        : distribution number
! real(dp) rmax         : Maximum orbital radius
!
! ----------------------------------------------------------------------
! Output :
! ----------------------------------------------------------------------
! All output quantities are in the mesh and moreMeshSubs modules.
! Data is stored in moreMeshSubs but accessed using the mesh module
! pointers.
!
!***********************************************************************
      subroutine setupExtMesh( distr, rmax )
      use precision, only : dp
      use mesh,      only : meshLim  ! My processor's box of mesh points
      use mesh,      only : ne       ! Points in rmax along each lat. vector
      use mesh,      only : nmsc     ! Mesh points in each supercell vector
      use mesh,      only : mop      ! Accumulated num. of orbital points
      use mesh,      only : nsp      ! Number of sub-points of each mesh point
      use mesh,      only : cmesh    ! Mesh cell vectors
      use mesh,      only : xdsp     ! Vector to mesh sub-points
      use mesh,      only : indexp   ! ranslation from extended to
                                     ! normal mesh index
      use mesh,      only : idop     ! Mesh-index span of points 
                                     ! within an atomic sphere
      use mesh,      only : xdop     ! Coordinates of (super)points 
                                     ! within an atomic sphere
      use moreMeshSubs, only : allocExtMeshDistr
      implicit none
C     Input parameters
      integer              :: distr
      real(dp)             :: rmax
C     Local parameters
      integer              :: i1, i2, i3, j1, j2, j3, k1, k2, k3,
     &                        i, j, k, nep, boxWidth(3), extWidth(3),
     &                        isp, myBox(2,3), nem(3), myExtBox(2,3)
      real(dp)             :: r, dxp(3), dx(3)
      logical              :: within
C     Functions
      real(dp)             :: dismin

      myBox(:,:) = meshLim(:,:) - 1

! Add 'wings extensions' to myBox of mesh points, containing all points 
! within the orbital spheres that may intersect myBox. Thus, wings must
! be one diameter wide.
      myExtBox(1,:) = myBox(1,:) - 2*ne(:)
      myExtBox(2,:) = myBox(2,:) + 2*ne(:)

! Find number of extended-box points.
      nem(:) = myExtBox(2,:) - myExtBox(1,:) + 1
      nep = nem(1) * nem(2) * nem(3)

!     Allocate memory for indexp(nep), idop(mop) and xdop(3,mop) of
!     the current data distribution
      call allocExtMeshDistr( distr, nep, mop )

! Find relationship between extended and normal box points
      boxWidth(:) =    myBox(2,:)    - myBox(1,:) + 1
      extWidth(:) = myExtBox(2,:) - myExtBox(1,:) + 1
      do i3 = myExtBox(1,3),myExtBox(2,3)
        do i2 = myExtBox(1,2),myExtBox(2,2)
          do i1 = myExtBox(1,1),myExtBox(2,1)

            ! Find periodic indexes within first supercell
            j1 = modulo( i1, nmsc(1) )
            j2 = modulo( i2, nmsc(2) )
            j3 = modulo( i3, nmsc(3) )

            ! Find indexes relative to box origins
            j1 = j1 -    myBox(1,1)
            j2 = j2 -    myBox(1,2)
            j3 = j3 -    myBox(1,3)
            k1 = i1 - myExtBox(1,1)
            k2 = i2 - myExtBox(1,2)
            k3 = i3 - myExtBox(1,3)

            ! Find combined mesh indexes.
            j = 1 + j1 + boxWidth(1)*j2 + boxWidth(1)*boxWidth(2)*j3
            k = 1 + k1 + extWidth(1)*k2 + extWidth(1)*extWidth(2)*k3

            ! Find myExtBox -> myBox index translation
            if (j1>=0 .and. j1<boxWidth(1) .and.
     .          j2>=0 .and. j2<boxWidth(2) .and.
     .          j3>=0 .and. j3<boxWidth(3)) then
              indexp(k) = j   ! Point is within myBox
            else
              indexp(k) = 0   ! Point is outside myBox
            endif ! (j>=myBox(1) and j<=myBox(2))

          enddo ! i1
        enddo ! i2
      enddo ! i3

! Find points within rmax (orbital points)
      mop = 0
      do i3 = -ne(3), ne(3)
        do i2 = -ne(2), ne(2)
          do i1 = -ne(1), ne(1)
            dxp(:) = cmesh(:,1) * i1 +
     &               cmesh(:,2) * i2 +
     &               cmesh(:,3) * i3 
            within = .false.
            do isp = 1,nsp
              dx(1:3) = dxp(1:3) + xdsp(1:3,isp) 
              r = dismin( cmesh, dx )
              if ( r .lt. rmax ) within = .true.
            enddo
            if ( within ) then
              mop = mop + 1
              ! Store index-distance and vector-distance to point.
              idop(mop)     = i1 + nem(1)*i2 + nem(1)*nem(2)*i3
              xdop(1:3,mop) = dxp(1:3)
            endif ! (within)
          enddo ! i1
        enddo ! i2
      enddo ! i3

      end subroutine setupExtMesh
