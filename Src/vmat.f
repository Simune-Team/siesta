C $Id: vmat.f,v 1.10.2.3 1999/06/09 09:24:07 emilio Exp $

      subroutine vmat( NO, indxuo, endC, listC, C, NSPmax, NSP,
     .                 NP, endCt, listCt, CtoCt,
     .                 VolCel, V, NVmax, numVs, listVs, Vs,
     .                 MaxAux, ilocal )

C ********************************************************************
C Finds the matrix elements of the potential.
C First version written by P.Ordejon.
C Name and interface modified by J.M.Soler. May'95.
C Re-ordered so that mesh is the outer loop and the orbitals are
C handled as lower-half triangular. J.D.Gale and J.M.Soler, Feb'99
C *********************** INPUT **************************************
C integer NO              : Number of basis orbitals
C integer indxuo(NO)      : Index of equivalent atom in unit cell
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
C integer MaxAux          : Dimension of auxiliary array ilocal
C ******************** INPUT AND OUTPUT *******************************
C real*8  Vs(NVmax,NO)    : Value of nonzero elements in each row of Vs
C                           to which the potential matrix elements are
C                           summed up
C ********************* AUXILIARY *************************************
C integer ilocal(MaxAux)  : Space that can be used freely between calls
C                           MaxAux must be at least NO
C *********************************************************************

      implicit none

      integer
     .   NO, NP, NVmax, NSP, NSPmax, MaxAux,
     .   endC(0:NO),  listC(*),
     .   endCt(0:NP), listCt(*), CtoCt(*),
     .   numVs(NO), listVs(NVmax,NO), ilocal(MaxAux), indxuo(NO)

      real
     .  C(NSPmax,*), V(NSPmax,NP)

      double precision
     .   VolCel, Vs(NVmax,NO)

      integer maxloc
      parameter ( maxloc = 300 )

      integer i, ii, il, imp, in, iorb(0:maxloc), ip, isp, iu,
     .        j, jl, jmp, jn, last, nlocal
      double precision Ci, Cj, dVol, Vij, Vlocal(0:maxloc,0:maxloc)
      
      call timer('vmat',1)

C  Check size of auxiliary array
      if (NO .gt. MaxAux) then
        stop 'VMAT: dimension MaxAux too small'
      endif

C  Find volume per mesh point
      dVol = VolCel / (NP*NSP)

C  Initialize iorb and Vlocal
      do il = 0,maxloc
        iorb(il) = 0
        do jl = 0,maxloc
          Vlocal(jl,il) = 0.d0
        enddo
      enddo
      last = 0

C  Full initialization of array ilocal done only once
      do j = 1, NO
        ilocal(j) = 0
      enddo

C  Loop over grid points
      do ip = 1,np

C  Check that Vlocal can contain Vij for this point.
        if (endCt(ip)-endCt(ip-1) .gt. maxloc)
     .    stop 'VMAT: parameter maxloc too small'

C  Find new required size of Vlocal
        nlocal = last
        do imp = 1+endCt(ip-1), endCt(ip)
          i = listCt(imp)
          if (ilocal(i) .eq. 0) nlocal = nlocal + 1
        enddo

C  If overflooded, add Vlocal to Vs and reinitialize it
        if (nlocal .gt. maxloc) then
          do il = 1,last
            i = iorb(il)
            iu = indxuo(i)
            do ii = 1, numVs(i)
              j = listVs(ii,i)
              jl = ilocal(j)
              Vs(ii,iu) = Vs(ii,iu) + dVol * Vlocal(jl,il)
            enddo
          enddo
          do il = 1,last
            i = iorb(il)
            ilocal(i) = 0
            iorb(il) = 0
            do jl = 1,last
              Vlocal(jl,il) = 0.d0
            enddo
          enddo
          last = 0
        endif

C  Look for required orbitals not yet in Vlocal
        if (nlocal .gt. last) then
          do imp = 1+endCt(ip-1), endCt(ip)
            i = listCt(imp)
            if (ilocal(i) .eq. 0) then
              last = last + 1
              ilocal(i) = last
              iorb(last) = i
            endif
          enddo
        endif

C  Loop on first orbital of mesh point
        do imp = 1+endCt(ip-1), endCt(ip)
          i = listCt(imp)
          il = ilocal(i)
          iu = indxuo(i)
          in = CtoCt(imp)

C  Loop on second orbital of mesh point
C  Notice that the loop runs only for jmp.le.imp
          do jmp = 1+endCt(ip-1), imp
            j = listCt(jmp)
            jl = ilocal(j)
            jn = CtoCt(jmp)

C  Loop over sub-points
            Vij = 0
            do isp = 1, nsp
              Ci = C(isp,in)
              Cj = C(isp,jn)
              Vij = Vij + V(isp,ip) * Ci * Cj
            enddo
            Vlocal(jl,il) = Vlocal(jl,il) + Vij
            if (imp .ne. jmp) then
              Vlocal(il,jl) = Vlocal(il,jl) + Vij
            endif

          enddo
        enddo
      enddo

C  Add final Vlocal to Vs
      do il = 1,last
        i = iorb(il)
        iu = indxuo(i)
        do ii = 1, numVs(i)
          j = listVs(ii,i)
          jl = ilocal(j)
          Vs(ii,iu) = Vs(ii,iu) + dVol * Vlocal(jl,il)
        enddo
      enddo

      call timer('vmat',2)
      end

