C $Id: rhoofd.f,v 1.10.2.3 1999/06/09 09:24:08 emilio Exp $

      subroutine rhoofd( NO, indxuo, endC,  listC,  C, NSPmax, NSP, 
     .                   NP, endCt, listCt, CtoCt,
     .                   NDmax, numDs, listDs, Dscf,
     .                   RHOscf, MaxAux, ilocal )

C ********************************************************************
C Finds the SCF density at the mesh points from the density matrix.
C Written by P.Ordejon and J.M.Soler. May'95.
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
C integer MaxAux          : Dimension of auxiliary array ilocal
C *********************** OUTPUT **************************************
C real*4  RHOscf(NSP,NP)  : SCF density at mesh points
C ********************* AUXILIARY *************************************
C integer ilocal(MaxAux)  : Space that can be used freely between calls
C                           MaxAux must be at least NO
C *********************************************************************

      implicit none

      integer
     .   NO, NP, NDmax, NSP, NSPmax, MaxAux, 
     .   endC(0:NO),  listC(*), 
     .   endCt(0:NP), listCt(*), CtoCt(*),
     .   numDs(NO), listDs(NDmax,NO), ilocal(MaxAux), indxuo(NO)
      real
     .   C(NSPmax,*), RHOscf(NSP,NP)
      double precision
     .   Dscf(NDmax,*)


      integer maxloc
      parameter ( maxloc = 300 )

      integer i, ii, il, imp, in, iorb(0:maxloc), ip, isp, iu,
     .        j, jl, jn, jmp, last

      double precision Ci, Cj, Dij, Dlocal(0:maxloc,0:maxloc)
      
      call timer('rhoofd',1)

C  Check size of auxiliary array
      if (NO .gt. MaxAux) then
        stop 'RHOOFD: dimension MaxAux too small'
      endif

C  Full initialization of arrays ilocal and iorb done only once
      do j = 1, NO
        ilocal(j) = 0
      enddo
      do il = 0,maxloc
        iorb(il) = 0
      enddo
      last = 0

C  Loop over grid points
      do ip = 1,np

C  Initialise RHOatm
        do isp = 1, nsp
          RHOscf(isp,ip) = 0.d0
        enddo

C  Check that Dlocal can contain Dij for this point.
        if (endCt(ip)-endCt(ip-1) .gt. maxloc)
     .    stop 'RHOOFD: parameter maxloc too small'

C  iorb(il)>0 means that row il of Dlocal must not be overwritten
C  iorb(il)=0 means that row il of Dlocal is empty
C  iorb(il)<0 means that row il of Dlocal contains a valid row of 
C             Dscf, but which is not required at this point
        do imp = 1+endCt(ip-1), endCt(ip)
          i = listCt(imp)
          il = ilocal(i)
          if (il.gt.0) iorb(il) = i
        enddo

C  Look for required rows of Dscf not yet stored in Dlocal
        do imp = 1+endCt(ip-1), endCt(ip)
          i = listCt(imp)
          if (ilocal(i) .eq. 0) then

C           Look for an available row in Dlocal
            do il = 1,maxloc
C             last runs circularly over rows of Dlocal
              last = last + 1
              if (last .gt. maxloc) last = 1
              if (iorb(last) .le. 0) goto 10
            enddo
              stop 'rhoofd: no slot available in Dlocal'
   10       continue

C  Copy row i of Dscf into row last of Dlocal
            j = abs(iorb(last))
            if (j.ne.0) ilocal(j) = 0
            ilocal(i) = last
            iorb(last) = i
            il = last
            iu = indxuo(i)
            do ii = 1, numDs(i)
              j = listDs(ii,i)
              jl = ilocal(j)
              Dlocal(il,jl) = Dscf(ii,iu)
              Dlocal(jl,il) = Dscf(ii,iu)
            enddo
          endif
        enddo

C  Loop on first orbital of mesh point
        do imp = 1+endCt(ip-1), endCt(ip)
          i = listCt(imp)
          il = ilocal(i)
          iu = indxuo(i)
          in = CtoCt(imp)

C  Loop on second orbital of mesh point
          do jmp = 1+endCt(ip-1), imp
            j = listCt(jmp)
            jl = ilocal(j)
            jn = CtoCt(jmp)
            if (imp .eq. jmp) then
              Dij = Dlocal(il,jl)
            else
              Dij = 2*Dlocal(il,jl)
            endif

C  Loop over sub-points
            do isp = 1, nsp
              Ci = C(isp,in)
              Cj = C(isp,jn)
              RHOscf(isp,ip) = RHOscf(isp,ip) + Dij * Ci * Cj
            enddo

          enddo
        enddo

C  Restore iorb for next point
        do imp = 1+endCt(ip-1), endCt(ip)
          i = listCt(imp)
          il = ilocal(i)
          iorb(il) = -i
        enddo

      enddo

      call timer('rhoofd',2)
      end
