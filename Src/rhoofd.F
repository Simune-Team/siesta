      subroutine rhoofd( NO, indxuo, C, maxC, NSP, 
     .                   NP, endCt, listCt, CtoCt,
     .                   NDmax, numDs, listDs, Dscf,
     .                   RHOscf, nuotot, NeedDscfL,
     .                   nrowsDscfL, iaorb, iphorb,
     .                   isa, listP2, dxa, xdop, xdsp,
     .                   DirectPhi )
C ********************************************************************
C Finds the SCF density at the mesh points from the density matrix.
C Written by P.Ordejon and J.M.Soler. May'95.
C Re-ordered so that mesh is the outer loop and the orbitals are
C handled as lower-half triangular. J.D.Gale and J.M.Soler, Feb'99
C Version of rhoofd that optionally uses a direct algorithm to save 
C memory. Modified by J.D.Gale, November'99
C *********************** INPUT **************************************
C integer NO              : Number of basis orbitals
C integer indxuo(NO)      : Index of equivalent atom in unit cell
C real*4  C(NSPmax,*)     : Values of nonzero elements in C
C integer maxC            : Second dimension of C
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
C integer NeedDscfL(NO)   : Pointer to local elements of Dscf needed
C integer nrowsDscfL      : Number of local rows of Dscf
C integer iaorb(*)        : Pointer to atom to which orbital belongs
C integer iphorb(*)       : Orbital index within each atom
C integer isa(*)          : Species index of all atoms
C integer listP2(*)       : Mesh point pointer for orbital
C real*8  dxa(3,*)        : Atom position within mesh-cell
C real*8  xdop(3,*)       : Vector to mesh points within rmax
C real*8  xdsp(3,*)       : Vector to mesh sub-points
C logical DirectPhi       : if .true. then C is not to be used and phi
C                         : values must be calculated on the fly
C *********************** OUTPUT **************************************
C real*4  RHOscf(NSP,NP)  : SCF density at mesh points
C *********************************************************************

C
C  Modules
C
      use precision
      use atmfuncs, only: rcut, phiatm
      use listsc_module, only: listsc
#ifdef MPI
      use mpi
      use parallel
#endif

      implicit none

      integer
     .   NO, NP, NDmax, NSP, nuotot, maxC,
     .   endCt(0:NP), listCt(*), CtoCt(*),
     .   numDs(NO), listDs(NDmax,NO), indxuo(NO),
     .   NeedDscfL(*), nrowsDscfL, iaorb(*),
     .   iphorb(*), isa(*), listP2(*)
      real
     .   C(NSP,maxC), RHOscf(NSP,NP)
      double precision
     .   Dscf(NDmax,*), dxa(3,*), xdop(3,*), xdsp(3,*)
      logical 
     .   DirectPhi
      external
     .   memory

      integer maxloc

      integer i, ii, il, imp, in, ip, isp, iu, iop, is, iphi, ia,
     .        j, jl, jn, jmp, last, ic, jc, ix

      integer, dimension(:), allocatable, save :: 
     .  ilocal, iorb

#ifdef MPI
      integer 
     .   MPIerror, BNode, Node, Nodes
      double precision, dimension(:), allocatable, save :: 
     .  buffer
      double precision, dimension(:,:), allocatable, save :: 
     .  DscfL
#endif

      double precision Ci, Cj, Dij, r2o, dxsp(3), r2sp, phip
      double precision, dimension(:,:), allocatable, save :: 
     .  Dlocal, gradCi
      real, dimension(:,:), allocatable, save :: 
     .  Clocal
      
      call timer('rhoofd',1)

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#endif

C  Find value of maxloc
      maxloc = 100
      do ip = 1,np
        maxloc = max(maxloc,endCt(ip)-endCt(ip-1))
      enddo

C  Allocate local memory
      allocate(ilocal(NO))
      call memory('A','I',NO,'rhoofd')
      allocate(iorb(0:maxloc))
      call memory('A','I',maxloc+1,'rhoofd')
      allocate(Dlocal(0:maxloc,0:maxloc))
      call memory('A','D',(maxloc+1)*(maxloc+1),'rhoofd')
      if (DirectPhi) then
        allocate(gradCi(3,NSP))
        call memory('A','D',3*NSP,'rhoofd')
        allocate(Clocal(NSP,maxloc))
        call memory('A','S',NSP*maxloc,'rhoofd')
      endif

C  Full initialization of arrays ilocal and iorb done only once
      do j = 1, no
        ilocal(j) = 0
      enddo
      do il = 0,maxloc
        iorb(il) = 0
      enddo
      last = 0

#ifdef MPI
C  Allocate local Dscf storage array
      allocate(DscfL(NDmax,nrowsDscfL))
      call memory('A','D',NDmax*nrowsDscfL,'rhoofd')
      allocate(buffer(NDmax))
      call memory('A','D',NDmax,'rhoofd')

C  Broadcast Dscf
      do i = 1,nuotot
        call WhichNodeOrb(i,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(i,Node,Nodes,iu)
          do ii = 1,numDs(i)
            buffer(ii) = Dscf(ii,iu)
          enddo
        endif
        call MPI_Bcast(buffer,numDs(i),DAT_double,BNode,
     .    MPI_Comm_World,MPIerror)
        if (NeedDscfL(i).gt.0) then
          do ii = 1,numDs(i)
            DscfL(ii,NeedDscfL(i)) = buffer(ii)
          enddo
        endif
      enddo
      call memory('D','D',size(buffer),'rhoofd')
      deallocate(buffer)
#endif

C  Just as a check zero Dlocal
      do i = 0,maxloc
        do j = 0,maxloc
          Dlocal(j,i) = 0.0
        enddo
      enddo

C  Loop over grid points
      do ip = 1,np

C  Initialise RHOatm
        do isp = 1, nsp
          RHOscf(isp,ip) = 0.d0
        enddo

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
            if (i .eq. iu) then
              do ii = 1, numDs(iu)
                j = listDs(ii,iu)
                jl = ilocal(j)
#ifdef MPI
                Dij = DscfL(ii,NeedDscfL(iu))
#else
                Dij = Dscf(ii,iu)
#endif
                Dlocal(il,jl) = Dij
                Dlocal(jl,il) = Dij
              enddo
            else
              do ii = 1, numDs(iu)
                j = listsc( i, iu, listDs(ii,iu) )
                jl = ilocal(j)
#ifdef MPI
                Dij = DscfL(ii,NeedDscfL(iu))
#else
                Dij = Dscf(ii,iu)
#endif
                Dlocal(il,jl) = Dij
                Dlocal(jl,il) = Dij
              enddo
            endif
          endif
        enddo

C  Loop on first orbital of mesh point
        ic = 0
        do imp = 1+endCt(ip-1), endCt(ip)
          ic = ic + 1
          i = listCt(imp)
          il = ilocal(i)
          iu = indxuo(i)
          in = CtoCt(imp)

          if (DirectPhi) then
C  Generate phi value
            iphi = iphorb(i)
            ia = iaorb(i)
            is = isa(ia)
            r2o = rcut(is,iphi)**2
            iop = listP2(in)
            do isp = 1,nsp
              do ix = 1,3
                dxsp(ix) = xdop(ix,iop) + xdsp(ix,isp) - dxa(ix,ia)
              enddo
              r2sp = dxsp(1)**2 + dxsp(2)**2 + dxsp(3)**2
              if (r2sp.lt.r2o) then
                call phiatm(is,iphi,dxsp,phip,gradCi(1,isp))
                Clocal(isp,ic) = phip
              else
                Clocal(isp,ic) = 0.0
              endif
            enddo
          endif

C  Loop on second orbital of mesh point
          jc = 0
          do jmp = 1+endCt(ip-1), imp
            jc = jc + 1
            j = listCt(jmp)
            jl = ilocal(j)
            jn = CtoCt(jmp)
            if (imp .eq. jmp) then
              Dij = Dlocal(il,jl)
            else
              Dij = 2*Dlocal(il,jl)
            endif

C  Loop over sub-points
            if (DirectPhi) then
              do isp = 1, nsp
                Ci = Clocal(isp,ic)
                Cj = Clocal(isp,jc)
                RHOscf(isp,ip) = RHOscf(isp,ip) + Dij * Ci * Cj
              enddo
            else
              do isp = 1, nsp
                Ci = C(isp,in)
                Cj = C(isp,jn)
                RHOscf(isp,ip) = RHOscf(isp,ip) + Dij * Ci * Cj
              enddo
            endif

          enddo
        enddo

C  Restore iorb for next point
        do imp = 1+endCt(ip-1), endCt(ip)
          i = listCt(imp)
          il = ilocal(i)
          iorb(il) = -i
        enddo

      enddo

#ifdef MPI
C  Deallocate local Dscf storage array
      call memory('D','D',size(DscfL),'rhoofd')
      deallocate(DscfL)
#endif

C  Free local memory
      call memory('D','I',size(ilocal),'rhoofd')
      deallocate(ilocal)
      call memory('D','I',size(iorb),'rhoofd')
      deallocate(iorb)
      call memory('D','D',size(Dlocal),'rhoofd')
      deallocate(Dlocal)
      if (DirectPhi) then
        call memory('D','D',size(gradCi),'rhoofd')
        deallocate(gradCi)
        call memory('D','S',size(Clocal),'rhoofd')
        deallocate(Clocal)
      endif

      call timer('rhoofd',2)
      end
