      subroutine rhooda( NO, endC,  listC,  C, NSPmax, NSP, NP, Datm,
     .                   RHOatm )

C ********************************************************************
C Finds the Harris density at the mesh points from the atomic
C occupations.
C Written by P.Ordejon and J.M.Soler. May'95.
C *********************** INPUT **************************************
C integer NO              : Number of basis orbitals
C integer endC(NO)        : Accumulated umber of nonzero elements in
C                           each row of C
C integer listC(*)        : List of nonzero elements in each row of C
C real*4  C(NSPmax,*)     : Values of nonzero elements in C
C integer NSPmax          : First dimension of array C
C integer NSP             : Number of sub-points of each mesh point
C integer NP              : Number of mesh points
C real*8  Datm(NO)        : Occupations of basis orbitals in free atom
C *********************** OUTPUT **************************************
C real*4  RHOatm(NSP,NP)  : Harris (sum of atoms) density at mesh points
C *********************************************************************

      implicit none

      integer          NO, NP, NSP, NSPmax, endC(0:NO), listC(*)
      real             C(NSPmax,*), RHOatm(NSP,NP)
      double precision Datm(NO)

      integer          i, in, ip, isp
      double precision Ci

C     Initialize RHOatm
      do 20 ip = 1, NP
        do 10 isp = 1, NSP
          RHOatm(isp,ip) = 0.d0
   10   continue
   20 continue

C     Loop on orbitals
      do 50 i = 1, NO

C       Loop on mesh points of orbital i
        do 40 in = 1+endC(i-1), endC(i)
          ip = listC(in)

C         Loop over sub-points
          do 30 isp = 1, NSP
            Ci = C(isp,in)
            RHOatm(isp,ip) = RHOatm(isp,ip) + Datm(i) * Ci * Ci
   30     continue

   40   enddo

   50 enddo

      end
