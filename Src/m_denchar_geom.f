! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module m_denchar_geom

      use precision

      public :: atompla, length, matvect, colinear, plane, crossv

      private

      CONTAINS


      SUBROUTINE ATOMPLA( NA, ORIGIN, XA, MROT, IDIMEN, ISCALE,
     .                   NATINPLA, INDICES, XAPLANE )
C **********************************************************************
C Find the coordinates of some selected atoms on the 
C plane-reference-frame. Atoms in plane will have the third coordinate 
C equal to zero.
C For 3D-grid, all the atoms are selected
C Coded by J. Junquera May '99
C Modified by P. Ordejon to include 3D capability; June 2003
C **********************************************************************

      USE FDF
      USE PARSE
      USE SYS

      IMPLICIT NONE

      INTEGER
     .  NA, NATINPLA, INDICES(NA), ISCALE

      DOUBLE PRECISION
     .  ORIGIN(3), XA(3,NA), MROT(3,3), XAPLANE(3,NA)

C ******* INPUT ********************************************************
C INTEGER NA             : Number of atoms
C REAL*8  ORIGIN(3)      : Origin of the plane reference frame
C REAL*8  XA(3,NA)       : Atomic coordinates in lattice reference frame
C REAL*8  MROT(3,3)      : Rotation matrix from the lattice reference frame 
C                          to the in-plane reference frame
C INTEGER IDIMEN         : Determines if run is plane or 3D-grid
C INTEGER ISCALE         : Units of the atomic coordinates
C                          ISCALE = 1 => bohrs, ISCALE = 2 => Ang
C ******* OUTPUT *******************************************************
C INTEGER NATINPLA       : Number of atoms whose coordinates will be rotated
C INTEGER INDICES(NA)    : Atomic indices of the atoms whose coordinates
C                          will be rotated
C REAL*8  XAPLANE(3,NA)  : Atomic coordinates in plane reference frame
C **********************************************************************

      CHARACTER  LINE*150

      INTEGER  NATINPL_DEFECT, IDIMEN, IUNIT, IAT, IX
 
      DOUBLE PRECISION  VAUX1(3), VAUX2(3)

      LOGICAL  ATINPLA

      TYPE(PARSED_LINE), POINTER :: P
      TYPE(BLOCK), POINTER       :: BP


C **********************************************************************
C INTEGER NATINPLA       : Number of atoms whose coordinates will be 
C                          rotated from the lattice reference frame to 
C                          the in-plane reference frame
C REAL*8  VAUX1(3)       : Auxiliar vector
C REAL*8  VAUX2(3)       : Auxiliar vector
C LOGICAL ATINPLA        : Is block Denchar.AtomsInPlane present in fdf?
C **********************************************************************

      IF (IDIMEN .EQ. 2) THEN

C Read fdf data block 'Denchar.AtomsInPlane' ---------------------------
        NATINPL_DEFECT = 0
        NATINPLA = 0
  
        NULLIFY(BP)
        IF ( .NOT. FDF_BLOCK('Denchar.AtomsInPlane',BP) )  GOTO 2000

        LOOP: DO
          IF (.NOT. FDF_BLINE(BP,LINE)) EXIT LOOP
          P=>DIGEST(LINE)
          IF (.NOT. MATCH(P,"I") ) 
     .         CALL DIE("Wrong format in Denchar.AtomsInPlane")
          NATINPLA = NATINPLA + 1
          INDICES(NATINPLA) = INTEGERS(P,1) 
          CALL DESTROY(P)
        ENDDO LOOP
        CALL DESTROY(BP) 
 2000   CONTINUE

      ELSE IF (IDIMEN .EQ. 3) THEN
        NATINPLA = NA
        DO IAT = 1, NA
          INDICES(IAT) = IAT
        ENDDO
      ELSE
        CALL DIE("Wrong IDIMEN in ATOMPLA")
      ENDIF
     

C Rotate the coordinates -----------------------------------------------
        DO IAT = 1, NATINPLA

          DO IX = 1,3
            VAUX1(IX) = XA(IX,INDICES(IAT)) - ORIGIN(IX)
          ENDDO

          CALL MATVECT(MROT, VAUX1, VAUX2)

          DO IX = 1,3
            XAPLANE(IX,INDICES(IAT)) = VAUX2(IX)
          ENDDO
  
C        IF (ISCALE .EQ. 2) THEN
C          DO IX = 1,3
C           XAPLANE(IX,INDICES(IAT))=XAPLANE(IX,INDICES(IAT))*0.529177D0
C          ENDDO
C        ENDIF

        ENDDO

        END subroutine atompla
C!------------------------------------------------------------------

      SUBROUTINE COLINEAR( COORPO, COLIN )

C **********************************************************************
C Checks if three points lye in the same straight line (if they are 
C colinear).
C Coded by J. Junquera November'98
C **********************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) ::  COORPO(3,3)
 
      LOGICAL, INTENT(OUT) :: COLIN

C **********************************************************************
C REAL*8 COORPO(3,3)   : Coordinates of three points. COORPO(POINT,IX)
C LOGICAL COLIN        : True => Three points are colinear
C                        False=> Three points NOT colinear
C **********************************************************************

C Internal variables ---------------------------------------------------

      INTEGER 
     .  IX
 
      DOUBLE PRECISION
     .  VEC1(3), VEC2(3), NORMAL(3), EPS

      DATA EPS /1.D-12/

      DO IX = 1,3
        VEC1(IX) = COORPO(2,IX) - COORPO(1,IX)
        VEC2(IX) = COORPO(3,IX) - COORPO(1,IX)
      ENDDO

      CALL CROSSV( VEC1, VEC2, NORMAL )
     
      IF( (ABS(NORMAL(1)) .LT. EPS) .AND. 
     .    (ABS(NORMAL(2)) .LT. EPS) .AND.
     .    (ABS(NORMAL(3)) .LT. EPS) ) COLIN = .TRUE.

      END subroutine colinear

C!-------------------------------------------------------------------------
      
      function length(a)
      real*8 length
*
*     This function calculates the length of a vector defined in
*    cartesian coordinates.
*
      real*8 a(3)
*
*     a(3) : cartesian coordinates of the vector.  
*
      length = sqrt(a(1)**2 + a(2)**2 + a(3)**2)

      return
      end function length

C!---------------------------------------------------------------------

      subroutine matvect(matrix,vector,product)
*
*     This subroutine calculates the vector result of the product of a 
*    matrix (3,3) and a vector.
*
      implicit none
      integer i,j
      real*8 matrix(3,3)
      real*8 vector(3),product(3)

      do i = 1,3
	 product(i) = 0.d0
	 do j = 1,3
	    product(i) = product(i) + matrix(i,j) * vector(j)
         enddo
      enddo

      end subroutine matvect

C!--------------------------------------------------------------------

      SUBROUTINE PLANE( NA, NPLAMAX, IDIMEN, OPTION, 
     .                  XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, 
     .                  NPX, NPY, NPZ, COORPO, NORMALP, 
     .                  DIRVER1, DIRVER2,
     .                  XA, NAPLA, INDICES, ISCALE,
     .                  LATPOINT, PLAPOINT, XAPLANE )

C **********************************************************************
C This subroutine calculates coordinates of the points of a plane 
C or in the 3D-grid in two reference frames: in the lattice reference 
C frame (lrf) and in another one centered on the plane in which the 
C third coordiante is always cero (prf)
C **********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::
     .   NPLAMAX, OPTION, NPX, NPY, NPZ, ISCALE, NA, NAPLA, INDICES(NA),
     .   IDIMEN

      DOUBLE PRECISION, INTENT(IN) ::
     .   XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .   NORMALP(3), COORPO(3,3), XA(3,NA)

      DOUBLE PRECISION, INTENT(IN) ::  DIRVER1(3), DIRVER2(3)

      DOUBLE PRECISION, INTENT(OUT) ::
     .   PLAPOINT(NPLAMAX,3), LATPOINT(NPLAMAX,3),
     .   XAPLANE(3,NA)

C **** INPUTS **********************************************************
C INTEGER NA         : Number of atoms in supercell
C INTEGER NPLAMAX    : Maximum number of points in the plane
C INTEGER IDIMEN     : Specify if the run is to plot quantities
C                      in a plane or in a 3D grid (2 or 3, respect)
C INTEGER OPTION     : Choose one of these options to define a plane:
C                      1 = components of the normal vector
C                      2 = equation of two straigth lines
C                      3 = coordinates of three points
C                      4 = indices of three atoms
C INTEGER NPX,NPY,NPZ: Number of points generated along x and y
C                      directions (and z for 3D-grids) 
C REAL*8 XMIN, XMAX  : Limits of the plane in the PRF for x-direction
C REAL*8 YMIN, YMAX  : Limits of the plane in the PRF for y-direction
C REAL*8 ZMIN, ZMAX  : Limits of the plane in the PRF for z-direction
C REAL*8  NORMAL(3)  : Components of the normal vector used to define
C                      the plane
C REAL*8  DIRVER1(3) : Components of the first vector contained
C                      in the plane
C                      (Only used if ioption = 2)
C REAL*8  DIRVER2(3) : Components of the second vector contained
C                       in the plane
C                       (Only used if ioption = 2)
C REAL*8  COORPO(3,3): Coordinates of the three points used to define
C                      the plane (Only used if ioption = 3)
C INTEGER NAPLA      : Number of atoms whose coordinates will be rotated
C INTEGER INDICES(NA): Indices of the atoms whose coordinates will
C                      be rotated from the lattice reference frame
C                      to the in-plane reference frame
C INTEGER ISCALE     : Unit of the points of the plane
C REAL*8  XA(3,NA)   : Atomic positions in cartesian coordinates
C                      (in bohr)
C **** OUTPUT **********************************************************
C REAL*8  LATPOINT(NPLAMAX,3) : Coordinates of the points 
C                               (lattice ref. frame)
C REAL*8  PLAPOINT(NPLAMAX,3) : Coordinates of the points 
C                               (plane refer. frame)
C REAL*8  XAPLANE(3,NA): Atomic coordinates in plane reference frame
C **********************************************************************

C Internal variables ---------------------------------------------------

      integer i,j
      real*8 modnor,moddi1
      real*8 normal(3)
      real*8 xplane(3),yplane(3)
      real*8 origin(3)
      real*8 mrot(3,3),inmrot(3,3)

      real*8 aux1(3), aux2(3)

C **********************************************************************
C REAL*8 MODNOR, MODDI1       : Length of the vector normal and dirver1
C REAL*8 XPLANE(3), YPLANE(3) 
C        ZPLANE(3)            : Plane reference frame
C REAL*8 ORIGIN(3)            : Coordinates of a point of the plane 
C                               that will be considered the origin of the 
C                               reference frame fixed in the plane. 
C REAL*8 MROT(3,3)            : matrix that relates the in plane
C                               reference frame with the lattice reference
C                               frame.
C REAL*8 INMROT(3,3)          : Inverse of the previous matrix
C **********************************************************************

      if(option .eq. 1) then

         do i = 1,3
            xplane(i) = coorpo(2,i) - coorpo(1,i)
         enddo

         modnor = length(normalp)
         moddi1 = length(xplane)

         do i = 1,3 
            normal(i) = (1/modnor)*normalp(i)
            xplane(i) = (1/moddi1)*xplane(i)
         enddo

c         write(*,*)(normal(i),i=1,3)
c         write(*,*)(xplane(i),i=1,3)
*=======================================================================
      else if(option .eq. 2) then     

         call crossv(dirver1,dirver2,normal)
         
         modnor = length(normal)
         moddi1 = length(dirver1)

         do i = 1,3
            normal(i) = (1/modnor) * normal(i)
            xplane(i) = (1/moddi1) * dirver1(i)
         enddo

*=======================================================================
      else if( (option .eq. 3) .or. (option .eq. 4) ) then

         do i = 1,3
            aux1(i) = coorpo(2,i) - coorpo(1,i)
            aux2(i) = coorpo(3,i) - coorpo(1,i)
         enddo
        
         call crossv(aux1,aux2,normal)
        
         modnor = length(normal)
         moddi1 = length(aux1)
        
         do i = 1,3
            normal(i) = (1/modnor)*normal(i)
            xplane(i) = (1/moddi1)*aux1(i)
         enddo
          
      endif

      call crossv(normal,xplane,yplane)

      do i = 1,3
         origin(i) = coorpo(1,i)
      enddo


*
*       We define now the matrix that describe the rotation.
*
      do j = 1,3
         mrot(j,1) = xplane(j) 
         mrot(j,2) = yplane(j)
         mrot(j,3) = normal(j)
      enddo
        
      do i = 1,3
         do j = 1,3
            inmrot(i,j) = mrot(j,i)
         enddo
      enddo

      call popla(nplamax,npx,npy,npz,xmin,xmax,ymin,ymax,zmin,zmax,
     .          idimen,plapoint)
        
      call rotation(nplamax,npx,npy,npz,mrot,inmrot,origin,
     .              plapoint,latpoint)

      call atompla( na, origin, xa, inmrot, idimen, iscale, 
     .              napla, indices, xaplane )

      end subroutine plane


*=======================================================================
*=======================================================================
      subroutine popla(nplamax,npx,npy,npz,xmin,xmax,ymin,ymax,
     .           zmin,zmax,idimen,plapoint)
*
*     This subroutine generates the coordinates of the points of the plane in
*    a reference frame fixed in the plane and where the z-axis is the normal
*    vector (in this way the third coordinate of the points is always zero).
*
*-----------------------------------------------------------------------
*     VARIABLES
*
        implicit none
    
        integer nplamax
        integer npx,npy,npz,i,k,l,m,idimen
        real*8 xmin,xmax,ymin,ymax,zmin,zmax,deltax,deltay,deltaz
        real*8 plapoint(nplamax,3)
*-----------------------------------------------------------------------
*     INICIALIZATION
*
        deltax=(xmax-xmin)/(npx - 1)
        deltay=(ymax-ymin)/(npy - 1)
        if (npz .ne. 1) then
          deltaz=(zmax-zmin)/(npz - 1)
        else
          deltaz=0.0d0
        endif

*-----------------------------------------------------------------------
*     COMPUTATIONAL BLOCK
*

C Do the loops differently if 2D or 3D. 
C This is necessary to comply with the ordering
C expected by Cube, and to maintain the former 
C ordering in the 2D case from former versions

        i = 1
        if (idimen .eq. 2) then
          do k = 1,npx
            do l = 1,npy
               plapoint(i,1) = xmin + (k-1)*deltax
               plapoint(i,2) = ymin + (l-1)*deltay
               plapoint(i,3) = 0.0d0
               i = i+1
            enddo
          enddo
        else if (idimen .eq. 3) then
          do m = 1,npz
             do l = 1,npy
               do k = 1,npx
                  plapoint(i,1) = xmin + (k-1)*deltax
                  plapoint(i,2) = ymin + (l-1)*deltay
                  plapoint(i,3) = zmin + (m-1)*deltaz
                  i = i+1
                enddo
            enddo
          enddo
        endif

        end subroutine popla

*=======================================================================
*=======================================================================       

      subroutine rotation(nplamax,npx,npy,npz,mrot,inmrot,origin,
     .                    plapoint,latpoint)
*
*     This subroutine makes the transformation of the coordinates of the points
*    of the plane in the plane reference frame to the lattice reference frame.
*    
*-----------------------------------------------------------------------
*     VARIABLES
*
      implicit none

      integer nplamax
      integer npx,npy,npz,i,j
      real*8 mrot(3,3),inmrot(3,3)
      real*8 origin(3),vaux1(3),vaux2(3)
      real*8 plapoint(nplamax,3),latpoint(nplamax,3)

*    
*     vaux1,vaux2 = auxiliar vectors needed in the computation of latpoint.
*-----------------------------------------------------------------------
*     COMPUTATIONAL BLOCK
*

      do i = 1, npx*npy*npz
         do j = 1,3
           vaux1(j) = plapoint(i,j)
         enddo

         call matvect(mrot,vaux1,vaux2)

         do j = 1,3
            latpoint(i,j) = vaux2(j) + origin(j)
         enddo
      enddo

      end subroutine rotation
*=======================================================================
*=======================================================================
      subroutine dotv(a,b,c)
*
*     This subroutine calculates the dot product of two vectors defined
*    in cartesian coordinates.
*
      real*8 a(3),b(3)
      real*8 c
*
*     a(3), b(3) : cartesian coordintes of the two vectors.
*     c : dot product of the two vectors.
*

      c = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

      end subroutine dotv

*=======================================================================
*=======================================================================
      SUBROUTINE CROSSV( A, B, C )

C **********************************************************************
C This subroutine calculates the cross product of two vectors defined 
C in cartesian coordinates.
C **********************************************************************

      DOUBLE PRECISION, INTENT(IN) ::
     .   A(3), B(3)

      DOUBLE PRECISION, INTENT(OUT) ::
     .   C(3)

C **** INPUT ***********************************************************
C REAL*8 A(3)  : Cartesian components of the first vector
C REAL*8 B(3)  : Cartesian components of the second vector
C **** OUTPUT **********************************************************
C REAL*8 C(3)  : Cartesian components of cross product vector
C **********************************************************************

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)

      END SUBROUTINE CROSSV

      end module m_denchar_geom

