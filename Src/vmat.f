      subroutine vmat( NO, endC,  listC, C, NSPmax, NSP,
     .                 NP, endCt, listCt, CtoCt,
     .                 VolCel, V, NVmax, numVs, listVs, Vs,
     .                 MaxVi, Vi )

C ********************************************************************
C Finds the matrix elements of the potential.
C Written by P.Ordejon.
C Name and interfase modified by J.M.Soler. May'95.
C *********************** INPUT **************************************
C integer NO              : Number of basis orbitals
C integer endC(NO)        : Accumulated umber of nonzero elements in
C                           each row of C
C integer listC(*)        : List of nonzero elements in each row of C
C real*4  C(NSPmax,*)     : Values of nonzero elements in C
C integer NSPmax          : First dimension of array C
C integer NSP             : Number of sub-points of each mesh point
C integer NP              : Number of columns in C (num. of mesh points)
C integer endCt(NP)       : Accumulated number of nonzero elements in
C                           each column of C
C integer listCt(*)       : List of nonzero elements in each column of C
C integer CtoCt(*)        : Index of matrix C where the nonzero
C                           elements in each column are stored
C real*8  VolCel          : Unit cell volume
C real*4  V(NP)           : Value of the potential at the mesh points
C integer NVmax           : First dimension of listV and Vs, and maxim.
C                           number of nonzero elements in any row of Vs
C integer numVs(NO)       : Number of nonzero elements in each row of Vs
C integer listVs(NVmax,NO): List of nonzero elements in each row of Vs
C integer MaxVi           : Dimension of auxiliary array Vi
C ******************** INPUT AND OUTPUT *******************************
C real*8  Vs(NVmax,NO)    : Value of nonzero elements in each row of Vs
C                           to which the potential matrix elements are
C                           summed up
C ********************* AUXILIARY *************************************
C real*8  Vi(MaxVi)       : Space that can be used freely between calls
C                           MaxVi must be at least NO
C *********************************************************************

      implicit none

      integer
     .   NO, NP, NVmax, NSP, NSPmax, MaxVi,
     .   endC(0:NO),  listC(*),
     .   endCt(0:NP), listCt(*), CtoCt(*),
     .   numVs(NO), listVs(NVmax,NO)

      real
     .  C(NSPmax,*), V(NSPmax,NP)

      double precision
     .   VolCel, Vs(NVmax,NO), Vi(NO)

      integer i, in, ip, isp, j, jn, kn
      double precision dVol
      
      call timer('vmat',1)

C     Check size of auxiliary array
      if (NO .gt. MaxVi) stop 'VMAT: dimension MaxVi too small'

C     Find volume per mesh point
      dVol = VolCel / (NP*NSP)

C     Full initialization of array Vi done only once (care!)
      do 90 j = 1, NO
        Vi(j) = 0.d0
   90 continue

C     Loop on orbitals
      do 160 i = 1, NO

C       Loop on mesh points of orbital i
        do 120 in = 1+endC(i-1), endC(i)
          ip = listC(in)

C         Loop on other non-zero orbitals in mesh point
          do 110 kn = 1+endCt(ip-1), endCt(ip)
            j = listCt(kn)
            jn = CtoCt(kn)

C           Loop over sub-points
            do 100 isp = 1, NSP
C             Accumulate <Ci|V|Cj> in prod(j) in full-row format
              Vi(j) = Vi(j) + V(isp,ip) * C(isp,in) * C(isp,jn)
  100       enddo

  110     enddo

  120   enddo

C       Collect <Ci|V|Cj> and store in sparse matrix format
C       Notice that in this version Vs is NOT initialized
        do 130 in = 1, numVs(i)
          j = listVs(in,i)
*         Vs(in,i) = Vi(j) * dVol
          Vs(in,i) = Vs(in,i) + Vi(j) * dVol
C         Restore Vi(j) for next orbital i
C         It is critical that listVs and listC are fully consistent!
*         Vi(j) = 0.d0
  130   enddo

C       A safer (but longer) restoration of Vi(j)
        do 150 in = 1+endC(i-1), endC(i)
          ip = listC(in)
          do 140 kn = 1+endCt(ip-1), endCt(ip)
            j = listCt(kn)
            Vi(j) = 0.d0
  140     enddo
  150   enddo

  160 enddo

      call timer('vmat',2)
      end

