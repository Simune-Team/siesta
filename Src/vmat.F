      subroutine vmat( NO, indxuo, C, maxC, NSP,
     .                 NP, NPF, endCt, listCt, CtoCt,
     .                 VolCel, V, NVmax, numVs, listVs, Vs,
     .                 nuo, nuotot, NeedVsL, nrowsVsL,
     .                 iaorb, iphorb, isa, listP2,
     .                 dxa, xdop, xdsp, DirectPhi )

C ********************************************************************
C Finds the matrix elements of the potential.
C First version written by P.Ordejon.
C Name and interface modified by J.M.Soler. May'95.
C Re-ordered so that mesh is the outer loop and the orbitals are
C handled as lower-half triangular. J.D.Gale and J.M.Soler, Feb'99
C Version of vmat that use a direct algorithm to save memory.
C Modified by J.D.Gale, November'99
C *********************** INPUT **************************************
C integer NO              : Number of basis orbitals
C integer indxuo(NO)      : Index of equivalent atom in unit cell
C real*4  C(NSPmax,*)     : Values of nonzero elements in C
C integer maxC            : Second dimension of C
C integer NSP             : Number of sub-points of each mesh point
C integer NP              : Number of columns in C (local)
C integer NPF             : Number of columns in C (global)
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
C integer NeedVsL(NO)     : Pointer to local elements of Vs needed
C integer nrowsVsL        : Number of rows of Vs needed
C integer iaorb(*)        : Pointer to atom to which orbital belongs
C integer iphorb(*)       : Orbital index within each atom
C integer isa(*)          : Species index of all atoms
C integer listP2(*)       : Mesh point pointer for orbital
C real*8  dxa(3,*)        : Atom position within mesh-cell
C real*8  xdop(3,*)       : Vector to mesh points within rmax
C real*8  xdsp(3,*)       : Vector to mesh sub-points
C logical DirectPhi       : if .true. then C is not to be used and phi
C                         : values must be calculated on the fly
C ******************** INPUT AND OUTPUT *******************************
C real*8  Vs(NVmax,NO)    : Value of nonzero elements in each row of Vs
C                           to which the potential matrix elements are
C                           summed up
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
     .   NO, NP, NPF, NVmax, NSP, NSPmax, nuo, nuotot,
     .   maxC, endCt(0:NP), listCt(*), CtoCt(*),
     .   numVs(NO), listVs(NVmax,NO), indxuo(NO),
     .   NeedVsL(*), nrowsVsL, iaorb(*),
     .   iphorb(*), isa(*), listP2(*)

      real
     .  C(NSP,maxC), V(NSP,NP)

      double precision
     .   VolCel, Vs(NVmax,NO), dxa(3,*), xdop(3,*), xdsp(3,*)

      logical
     .   DirectPhi

      integer maxloc

      integer i, ii, il, imp, in, ip, isp, iu, ic, jc, iop,
     .        j, jl, jmp, jn, last, nlocal, ia, iphi, is, ix

      integer, dimension(:), allocatable, save :: ilocal, iorb

#ifdef MPI
      integer MPIerror, nVsL, nVsLmax, Node, Nodes, iul
      double precision, dimension(:,:), allocatable, save :: 
     .  Vi, VsL
      integer, dimension(:), allocatable, save :: 
     .  nVsLPtr
#endif
      double precision Ci, Cj, dVol, Vij, r2o, r2sp, dxsp(3), phip
      double precision, dimension(:,:), allocatable, save :: 
     .  Vlocal, gradCi
      real, dimension(:,:), allocatable, save ::
     .  Clocal
      
      call timer('vmat',1)

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
      call memory('A','I',NO,'vmat')
      allocate(iorb(0:maxloc))
      call memory('A','I',maxloc+1,'vmat')
      allocate(Vlocal(0:maxloc,0:maxloc))
      call memory('A','D',(maxloc+1)*(maxloc+1),'vmat')
      if (DirectPhi) then
        allocate(gradCi(3,NSP))
        call memory('A','D',3*NSP,'vmat')
        allocate(Clocal(NSP,maxloc))
        call memory('A','S',NSP*maxloc,'vmat')
      endif

C  Find volume per mesh point
      dVol = VolCel / (NPF*NSP)

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

#ifdef MPI
C  Allocate local Vs storage array
      allocate(VsL(NVmax,nrowsVsL))
      call memory('A','D',NVmax*nrowsVsL,'vmat')

C  Initialise VsL
      do i = 1,nrowsVsL
        do ii = 1,NVmax
          VsL(ii,i) = 0.0d0
        enddo
      enddo
#endif

C  Loop over grid points
      do ip = 1,np

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
#ifdef MPI
            iul = NeedVsL(iu)
#endif
            if (i .eq. iu) then
              do ii = 1, numVs(iu)
                j = listVs(ii,iu)
                jl = ilocal(j)
#ifdef MPI
                VsL(ii,iul) = VsL(ii,iul) + dVol * Vlocal(jl,il) 
#else
                Vs(ii,iu) = Vs(ii,iu) + dVol * Vlocal(jl,il) 
#endif
              enddo
            else
              do ii = 1, numVs(iu)
                j = listsc( i, iu, listVs(ii,iu) )
                jl = ilocal(j)
#ifdef MPI
                VsL(ii,iul) = VsL(ii,iul) + dVol * Vlocal(jl,il) 
#else
                Vs(ii,iu) = Vs(ii,iu) + dVol * Vlocal(jl,il) 
#endif
              enddo
            endif
          enddo
          do il = 1,last
            i = iorb(il)
            ilocal(i) = 0
            iorb(il) = 0
            do jl = 1,last
              Vlocal(jl,il) = 0.0d0
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
C  Notice that the loop runs only for jmp.le.imp
          jc = 0
          do jmp = 1+endCt(ip-1), imp
            jc = jc + 1
            j = listCt(jmp)
            jl = ilocal(j)
            jn = CtoCt(jmp)

C  Loop over sub-points
            Vij = 0.0d0
            if (DirectPhi) then
              do isp = 1, nsp
                Ci = Clocal(isp,ic)
                Cj = Clocal(isp,jc)
                Vij = Vij + V(isp,ip) * Ci * Cj
              enddo
            else
              do isp = 1, nsp
                Ci = C(isp,in)
                Cj = C(isp,jn)
                Vij = Vij + V(isp,ip) * Ci * Cj
              enddo
            endif
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
#ifdef MPI
        iul = NeedVsL(iu)
#endif
        if (i .eq. iu) then
          do ii = 1, numVs(iu)
            j = listVs(ii,iu)
            jl = ilocal(j)
#ifdef MPI
            VsL(ii,iul) = VsL(ii,iul) + dVol * Vlocal(jl,il) 
#else
            Vs(ii,iu) = Vs(ii,iu) + dVol * Vlocal(jl,il) 
#endif
          enddo
        else
          do ii = 1, numVs(iu)
            j = listsc( i, iu, listVs(ii,iu) )
            jl = ilocal(j)
#ifdef MPI
            VsL(ii,iul) = VsL(ii,iul) + dVol * Vlocal(jl,il) 
#else
            Vs(ii,iu) = Vs(ii,iu) + dVol * Vlocal(jl,il) 
#endif
          enddo
        endif
      enddo

#ifdef MPI
      call MPI_AllReduce(nrowsVsL,nVsLmax,1,MPI_integer,
     .  MPI_max,MPI_Comm_World,MPIerror)
      allocate(Vi(NVmax,nVsLmax))
      call memory('A','D',NVmax*nVsLmax,'vmat')
      allocate(nVsLPtr(no))
      call memory('A','I',no,'vmat')

C  Reduce VsL back to global storage
      do in = 0,Nodes-1
        if (in.eq.Node) nVsL = nrowsVsL
        call MPI_Bcast(nVsL,1,MPI_integer,in,MPI_Comm_World,MPIerror)
        if (in.eq.Node) then
          do i = 1,no
            nVsLPtr(i) = NeedVsL(i)
          enddo
        endif
        call MPI_Bcast(nVsLPtr,no,MPI_integer,in,MPI_Comm_World,
     .    MPIerror)
        if (in.eq.Node) then
          do i = 1,nrowsVsL
            do ii = 1,NVmax
              Vi(ii,i) = VsL(ii,i)
            enddo
          enddo
        endif
        call MPI_Bcast(Vi(1,1),NVmax*nVsL,DAT_double,in,
     .    MPI_Comm_World,MPIerror)
        do i = 1,nuo
          call LocalToGlobalOrb(i,Node,Nodes,iu)
          iul = nVsLPtr(iu)
          if (iul.gt.0) then
            do ii = 1,numVs(iu)
              Vs(ii,i) = Vs(ii,i) + Vi(ii,iul)
            enddo
          endif
        enddo
      enddo
      call memory('D','I',size(nVsLPtr),'vmat')
      deallocate(nVsLPtr)
      call memory('D','D',size(Vi),'vmat')
      deallocate(Vi)

C  Deallocate local Vs storage array
      call memory('D','D',size(VsL),'vmat')
      deallocate(VsL)
#endif

C  Free local memory
      call memory('D','D',size(Vlocal),'vmat')
      deallocate(Vlocal)
      call memory('D','I',size(iorb),'vmat')
      deallocate(iorb)
      call memory('D','I',size(ilocal),'vmat')
      deallocate(ilocal)
      if (DirectPhi) then
        call memory('D','D',size(gradCi),'vmat')
        deallocate(gradCi)
        call memory('D','S',size(Clocal),'vmat')
        deallocate(Clocal)
      endif

      call timer('vmat',2)
      return
      end

