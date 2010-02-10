      subroutine setupExtMesh( distr, rmax )
      use precision,    only : dp
      use mesh,         only : meshLim, myExtBox, ne, nem, nmsc, mop,
     &                         indexp, idop, xdop, cmesh, xdsp, nsp
      use moreMeshSubs, only : allocExtMeshDistr
      implicit none
C     Input parameters
      integer              :: distr
      real(dp)             :: rmax
C     Local parameters
      integer              :: i1, i2, i3, j1, j2, j3, k1, k2, k3,
     &                        i, j, k, nep, boxWidth(3), extWidth(3),
     &                        isp, myBox(2,3)
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
