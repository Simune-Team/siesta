      SUBROUTINE PLANE( NA, NPLAMAX, OPTION, XMIN, XMAX, YMIN, YMAX,
     .                  NPX, NPY, COORPO, NORMALP, DIRVER1, DIRVER2,
     .                  XA, NAPLA, INDICES, ISCALE,
     .                  LATPOINT, PLAPOINT, XAPLANE )

C **********************************************************************
C This subroutine calculates coordinates of the points of a plane in 
C two reference frames: in the lattice reference frame (lrf) and in another 
C one centered on the plane in which the third coordiante 
C is always cero (prf)
C **********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::
     .   NPLAMAX, OPTION, NPX, NPY, ISCALE, NA, NAPLA, INDICES(NA)

      DOUBLE PRECISION, INTENT(IN) ::
     .   XMIN, XMAX, YMIN, YMAX,
     .   NORMALP(3), COORPO(3,3), XA(3,NA)

      DOUBLE PRECISION, INTENT(INOUT) ::
     .   DIRVER1(3), DIRVER2(3)

      DOUBLE PRECISION, INTENT(OUT) ::
     .   PLAPOINT(NPLAMAX,3), LATPOINT(NPLAMAX,3),
     .   XAPLANE(3,NA)

C **** INPUTS **********************************************************
C INTEGER NA         : Number of atoms in supercell
C INTEGER NPLAMAX    : Maximum number of points in the plane
C INTEGER OPTION     : Choose one of these options to define a plane:
C                      1 = components of the normal vector
C                      2 = equation of two straigth lines
C                      3 = coordinates of three points
C                      4 = indices of three atoms
C INTEGER NPX,NPY    : Number of points generated along x and y
C                      direction in a system of reference in which
C                      the third components od the points of the plane is
C                      zero (Plane Reference Frame; PRF)
C REAL*8 XMIN, XMAX  : Limits of the plane in the PRF for x-direction
C REAL*8 YMIN, YMAX  : Limits of the plane in the PRF for y-direction
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
C INTEGER ISCALE     : Unit if the points of the plane
C REAL*8  XA(3,NA)   : Atomic positions in cartesian coordinates
C                      (in bohr)
C **** OUTPUT **********************************************************
C REAL*8  LATPOINT(NPLAMAX,3) : Coordinates of the points of the plane
C                               (lattice ref. frame)
C REAL*8  PLAPOINT(NPLAMAX,3) : Coordinates of the points of the plane
C                               (plane refer. frame)
C REAL*8  XAPLANE(3,NA): Atomic coordinates in plane reference frame
C **********************************************************************

C Internal variables ---------------------------------------------------

      integer i,j
      real*8 length,modnor,moddi1
      real*8 normal(3)
      real*8 xplane(3),yplane(3)
      real*8 origin(3)
      real*8 mrot(3,3),inmrot(3,3)

      external length, atompla


C **********************************************************************
C REAL*8 MODNOR, MODDI1       : Length of the vector normal and dirver1
C REAL*8 XPLANE(3), YPLANE(3) : Plane reference frame
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
            dirver1(i) = coorpo(2,i) - coorpo(1,i)
            dirver2(i) = coorpo(3,i) - coorpo(1,i)
         enddo
        
         call crossv(dirver1,dirver2,normal)
        
         modnor = length(normal)
         moddi1 = length(dirver1)
        
         do i = 1,3
            normal(i) = (1/modnor)*normal(i)
            xplane(i) = (1/moddi1)*dirver1(i)
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

      call popla(nplamax,npx,npy,xmin,xmax,ymin,ymax,plapoint)
        
      call rotation(nplamax,npx,npy,mrot,inmrot,origin,
     .              plapoint,latpoint)

      call atompla( na, origin, xa, inmrot, napla, indices, iscale,
     .              xaplane )

      return
      end


*=======================================================================
*=======================================================================
      subroutine popla(nplamax,npx,npy,xmin,xmax,ymin,ymax,plapoint)
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
        integer npx,npy,i,k,l
        real*8 xmin,xmax,ymin,ymax,deltax,deltay
        real*8 plapoint(nplamax,3)
*-----------------------------------------------------------------------
*     INICIALIZATION
*
        i = 1
        deltax=(xmax-xmin)/(npx - 1)
        deltay=(ymax-ymin)/(npy - 1)

*-----------------------------------------------------------------------
*     COMPUTATIONAL BLOCH
*
        do k = 1,npx
           do l = 1,npy
             plapoint(i,1) = xmin + (k-1)*deltax
             plapoint(i,2) = ymin + (l-1)*deltay
             plapoint(i,3) = 0.d0
             i = i+1
          enddo
        enddo

	return
	end
*=======================================================================
*=======================================================================       

      subroutine rotation(nplamax,npx,npy,mrot,inmrot,origin,
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
      integer npx,npy,i,j
      real*8 mrot(3,3),inmrot(3,3)
      real*8 origin(3),vaux1(3),vaux2(3)
      real*8 plapoint(nplamax,3),latpoint(nplamax,3)

      external matvect
*    
*     vaux1,vaux2 = auxiliar vectors needed in the computation of latpoint.
*-----------------------------------------------------------------------
*     COMPUTATIONAL BLOCH
*

      do i = 1, npx*npy
         do j = 1,3
           vaux1(j) = plapoint(i,j)
         enddo

         call matvect(mrot,vaux1,vaux2)

         do j = 1,3
            latpoint(i,j) = vaux2(j) + origin(j)
         enddo
      enddo

      return
      end
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

      return
      end
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

      RETURN
  
      END SUBROUTINE CROSSV

