      subroutine rhoofd( NO, endC,  listC,  C, NSPmax, NSP, 
     .                   NP, endCt, listCt, CtoCt,
     .                   NDmax, numDs, listDs, Dscf,
     .                   RHOscf, MaxDi, Di )

C ********************************************************************
C Finds the SCF density at the mesh points from the density matrix.
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
C integer endCt(NP)       : Accumulated number of nonzero elements in
C                           each column of C
C integer listCt(*)       : List of nonzero elements in each column of C
C integer CtoCt(*)        : Index of matrix C where the nonzero
C                           elements in each column are stored
C integer NDmax           : First dimension of listD and Dscf, and
C                           maximum number of nonzero elements in
C                           any row of Dscf
C integer numDs(NO)       : Number of nonzero elemts in each row of Dscf
C integer listDs(NVmax,NO): List of nonzero elements in each row of Dscf
C real*8  Dscf(NVmax,NO)  : Value of nonzero elemens in each row of Dscf
C integer MaxDi           : Dimension of auxiliary array Di
C *********************** OUTPUT **************************************
C real*4  RHOscf(NSP,NP)  : SCF density at mesh points
C ********************* AUXILIARY *************************************
C real*8  Di(MaxDi)       : Space that can be used freely between calls
C                           MaxDi must be at least NO
C *********************************************************************

      implicit none

      integer
     .   NO, NP, NDmax, NSP, NSPmax, MaxDi,
     .   endC(0:NO),  listC(*),
     .   endCt(0:NP), listCt(*), CtoCt(*),
     .   numDs(NO), listDs(NDmax,NO)
      real
     .   C(NSPmax,*), RHOscf(NSP,NP)
      double precision
     .   Dscf(NDmax,NO), Di(NO)

      integer i, in, ip, isp, j, jn, kn
      double precision Ci, Cj, Dij
      
      call timer('rhoofd',1)

C     Check size of auxiliary array
      if (NO .gt. MaxDi) stop 'RHOOFD: dimension MaxDi too small'

C     Initialize RHOscf
      do 80 ip = 1, NP
        do 70 isp = 1, NSP
          RHOscf(isp,ip) = 0.d0
   70   continue
   80 continue

C     Full initialization of array Di done only once
      do 90 j = 1, NO
        Di(j) = 0.d0
   90 continue

C     Loop on orbitals
      do 160 i = 1, NO

C       Copy row i of Dij from sparse format in Dscf to
C         full-row format in Di
        do 100 in = 1, numDs(i)
          j = listDs(in,i)
          Di(j) = Dscf(in,i)
  100   continue

C       Loop on mesh points of orbital i
        do 120 in = 1+endC(i-1), endC(i)
          ip = listC(in)

C         Loop over non-zero orbitals in mesh point
          do 110 kn = 1+endCt(ip-1), endCt(ip)
            j = listCt(kn)
            jn = CtoCt(kn)
            Dij = Di(j)
*           write(6,'(a,i3,5i6,i3,i6)')
*    .        'rhoofd: i,endC(i-1),in,ip,endCt(ip-1),kn,j,jn=',
*    .        i, endC(i-1), in, ip, endCt(ip-1), kn, j, jn

C           Loop over sub-points
            do 105 isp = 1, NSP
              Ci = C(isp,in)
              Cj = C(isp,jn)
              RHOscf(isp,ip) = RHOscf(isp,ip) + Dij * Ci * Cj
  105       continue

  110     enddo

  120   enddo

C       Restore Di(j) for next orbital i
        do 150 in = 1, numDs(i)
          j = listDs(in,i)
          Di(j) = 0.d0
  150   continue

  160 enddo

      call timer('rhoofd',2)
      end
