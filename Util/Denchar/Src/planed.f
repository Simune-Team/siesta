      subroutine plane(nplamax, option,xmin,xmax,ymin,ymax,
     .                 npx,npy,coorpo,normalp,dirver1,dirver2,
     .                 latpoint,plapoint,na,napla,indices,iscale,
     .                 xa,xaplane)
*
*     This subroutine calculates coordinates of the points of one plane in 
*    two reference frames: in the lattice reference frame (lrf) and in another 
*    one centered on the plane in which the thid coordiante is always cero (prf)
*  INPUTS:  
*    option: Choose one of these options to define a plane:
*                   1 = components of the normal vector
*                   2 = equation of two straigth lines
*                   3 = coordinates of three points
*                   4 = indexes of three atoms
*    xmin,xmax,ymin,ymax = limits of the portion of the plane simulated (prf).
*    npx,npy = number of points generated in the directions x and y in prf 
*    coorpo(3,3) = coordinates of the three points utilized to define the plane
*                  (if option = 3) (lrf)
*    normal(3) = components of the normal vector.
*    dirver(2,3) = components of the directions vectors. 
*    na = number of atoms
*    indices(na) = indices ofthe atoms whose coordinates will be rotated
*                  from the lattice reference frame to the in-plane 
*                  reference frame
*    iscale = units of the points of the plane
*    xa(3,na) = atomic coordinates in lattice reference frame
*    xaplane(3,na) = atomic coordinates in plane reference frame
*
*    Coded by J.Junquera (March/97)
*-----------------------------------------------------------------------
*     VARIABLES
*
      implicit none

      integer i,j,option,npx,npy,nplamax, na, napla, indices(na), iscale
      real*8 xmin,xmax,ymin,ymax,length,modnor,moddi1
      real*8 normal(3),normalp(3), xa(3,na), xaplane(3,na)
      real*8 plapoint(nplamax,3)
      real*8 latpoint(nplamax,3)
      real*8 coorpo(3,3)
      real*8 dirver1(3),dirver2(3)
      real*8 xplane(3),yplane(3)
      real*8 origin(3)
      real*8 mrot(3,3),inmrot(3,3)
      external length, atompla

*
*    modnor,moddi1 = length of the vecotr normal and dirver1
*    xplane,yplane = plane reference frame. 
*    plapoint(3) = coordinates of some points of the plane(plane refer. frame)
*    latpoint(3) = coordinates of some points of the plane(lattice ref. frame)
*    origin(3) = coordinates of a point of the plane that will be considered 
*                the origin of the reference frame fixed in the plane. 
*    mrot(3,3) = matriz del cambio de base de las coordenadas en el sistema
*                de referencia del plano al sistema de referencia de la red.
*    inmrot(3,3) = inversa de la matriz anterior.
*

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
      subroutine crossv(a,b,c)
*
*     This subroutine calculates the cross product of two vectors defined in
*    cartesian coordinates.
*
      real*8 a(3),b(3),c(3)
*
*     a(3),b(3),c(3) :cartesian coordiantes of the three vectors.
*

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)

      return
      end

