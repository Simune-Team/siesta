      subroutine vmat( no, indxuo, np, dvol, V, nvmax, 
     .                 numVs, listVsptr, listVs, Vs, 
     .                 nuo, nuotot, iaorb, iphorb, isa )

C ********************************************************************
C Finds the matrix elements of the potential.
C First version written by P.Ordejon.
C Name and interface modified by J.M.Soler. May'95.
C Re-ordered so that mesh is the outer loop and the orbitals are
C handled as lower-half triangular. J.D.Gale and J.M.Soler, Feb'99
C Version of vmat that use a direct algorithm to save memory.
C Modified by J.D.Gale, November'99
C *********************** INPUT **************************************
C integer no              : Number of basis orbitals
C integer indxuo(no)      : Index of equivalent atom in unit cell
C integer np              : Number of columns in C (local)
C real*8  dvol            : Volume per mesh point
C real*4  V(np)           : Value of the potential at the mesh points
C integer nvmax           : First dimension of listV and Vs, and maxim.
C                           number of nonzero elements in any row of Vs
C integer numVs(nuo)      : Number of non-zero elements in a row of Vs
C integer listVsptr(nuo)  : Pointer to the start of rows in listVs
C integer listVs(nvmax)   : List of non-zero elements of Vs
C integer iaorb(*)        : Pointer to atom to which orbital belongs
C integer iphorb(*)       : Orbital index within each atom
C integer isa(*)          : Species index of all atoms
C ******************** INPUT AND OUTPUT *******************************
C real*8  Vs(nvmax)       : Value of nonzero elements in each row 
C                           of Vs to which the potential matrix 
C                           elements are summed up
C *********************************************************************

C  Modules
      use precision
      use atmfuncs, only: rcut, phiatm, all_phi, nsmax=>nspecies
      use listsc_module, only: listsc
      use mesh, only: dxa, nsp, xdop, xdsp
      use meshdscf
      use meshphi

      implicit none

C Argument types and dimensions
      integer
     .   no, np, nvmax, nuo, nuotot, indxuo(no), iaorb(*),
     .   iphorb(*), isa(*), numVs(nuo), listVsptr(nuo), listVs(nvmax)
      real
     .   V(nsp,np)
      real*8
     .   dvol, Vs(nvmax)

C Internal variables and arrays

      integer, parameter ::
     .  minloc = 100,  ! Min buffer size for local copy of Dscf
     .  maxoa  = 100   ! Max # of orbitals per atom
      integer
     .  i, ia, ic, ii, il, imp, in, ind, iop, ip, iphi, io,
     .  is, isp, iu, iul, ix, j, jc, jl, jmp, jn, 
     .  last, lasta, lastop, maxloc, maxloc2, nc, nlocal, 
     .  nphiloc, nvmaxl
      integer, dimension(:), allocatable, save ::
     .  ilc, ilocal, iorb
      logical
     .  Parallel
      real*8
     .  Vij, dxsp(3), phia(maxoa,nsp), r2cut(nsmax), r2sp
      real*8, dimension(:),   allocatable, save :: 
     .  VClocal
      real*8, dimension(:,:), allocatable, save :: 
     .  Clocal, Vlocal

C  Start time counter
      call timer('vmat',1)

C  Set algorithm logical
      Parallel = (nuo .ne. nuotot)

C  Find value of maxloc

      maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
      maxloc = maxloc2 + minloc
      maxloc = min( maxloc, no )

C  Allocate local memory
      allocate(ilocal(no))
      call memory('A','I',no,'vmat')
      allocate(ilc(maxloc2))
      call memory('A','I',maxloc2,'vmat')
      allocate(iorb(0:maxloc))
      call memory('A','I',maxloc+1,'vmat')
      allocate(Vlocal(0:maxloc,0:maxloc))
      call memory('A','D',(maxloc+1)*(maxloc+1),'vmat')
      allocate(Clocal(nsp,maxloc2))
      call memory('A','D',nsp*maxloc2,'vmat')
      allocate(VClocal(nsp))
      call memory('A','D',nsp,'vmat')

      if (Parallel) then
        nvmaxl = listdlptr(nrowsDscfL) + numdl(nrowsDscfL)
        allocate(DscfL(nvmaxl,1))
        call memory('A','D',nvmaxl,'meshdscf')
        DscfL(1:nvmaxl,1) = 0.0d0
      endif

C  Full initializations done only once
      ilocal(1:no) = 0
      iorb(0:maxloc) = 0
      Vlocal(:,:) = 0.0
      last = 0

C  Find atomic cutoff radii
      r2cut(:) = 0.0d0
      do i = 1,nuotot
        ia = iaorb(i)
        is = isa(ia)
        io = iphorb(i)
        r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
      enddo

C  Loop over grid points
      do ip = 1,np

C  Find number of nonzero orbitals at this point
        nc = endpht(ip) - endpht(ip-1)

C  Find new required size of Vlocal
        nlocal = last
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          if (ilocal(i) .eq. 0) nlocal = nlocal + 1
        enddo

C  If overflooded, add Vlocal to Vs and reinitialize it
        if (nlocal .gt. maxloc) then
          do il = 1,last
            i = iorb(il)
            iu = indxuo(i)
            if (Parallel) then
              iul = NeedDscfL(iu)
              if (i .eq. iu) then
                do ii = 1, numdl(iul)
                  ind = listdlptr(iul)+ii
                  j = listdl(ind)
                  jl = ilocal(j)
                  DscfL(ind,1) = DscfL(ind,1) + dVol * Vlocal(jl,il) 
                enddo
              else
                do ii = 1, numdl(iul)
                  ind = listdlptr(iul)+ii
                  j = listsc( i, iu, listdl(ind) )
                  jl = ilocal(j)
                  DscfL(ind,1) = DscfL(ind,1) + dVol * Vlocal(jl,il) 
                enddo
              endif
            else
              if (i .eq. iu) then
                do ii = 1, numVs(iu)
                  ind = listVsptr(iu)+ii
                  j = listVs(ind)
                  jl = ilocal(j)
                  Vs(ind) = Vs(ind) + dVol * Vlocal(jl,il) 
                enddo
              else
                do ii = 1, numVs(iu)
                  ind = listVsptr(iu)+ii
                  j = listsc( i, iu, listVs(ind) )
                  jl = ilocal(j)
                  Vs(ind) = Vs(ind) + dVol * Vlocal(jl,il) 
                enddo
              endif
            endif
          enddo
          ilocal(iorb(1:last)) = 0
          iorb(1:last) = 0
          Vlocal(1:last,1:last) = 0.0d0
          last = 0
        endif

C  Look for required orbitals not yet in Vlocal
        if (nlocal .gt. last) then
          do ic = 1,nc
            imp = endpht(ip-1) + ic
            i = lstpht(imp)
            if (ilocal(i) .eq. 0) then
              last = last + 1
              ilocal(i) = last
              iorb(last) = i
            endif
          enddo
        endif

C  Loop on first orbital of mesh point
        lasta=0
        lastop=0
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          il = ilocal(i)
          iu = indxuo(i)
          ia = iaorb(i)
          is = isa(ia)
          iop = listp2(imp)
          ilc(ic) = il

C  Generate or retrieve phi values
          if (DirectPhi) then
            if (ia.ne.lasta .or. iop.ne.lastop) then
              lasta = ia
              lastop = iop
              do isp = 1,nsp
                do ix = 1,3
                  dxsp(ix) = xdsp(ix,isp) + xdop(ix,iop) - dxa(ix,ia)
                enddo
                r2sp = sum(dxsp**2)
                if (r2sp.lt.r2cut(is)) then
                  call all_phi( is, +1, dxsp, nphiloc, phia(:,isp) )
                else
                  phia(:,isp) = 0.d0
                endif
              enddo
            endif
            iphi = iphorb(i)
            do isp = 1,nsp
              Clocal(isp,ic) = phia(iphi,isp)
            enddo
          else
            do isp = 1,nsp
              Clocal(isp,ic) = phi(isp,imp)
            enddo
          endif

C  Pre-multiply V and Clocal(,ic)
          do isp = 1,nsp
            VClocal(isp) = V(isp,ip) * Clocal(isp,ic)
          enddo

C  Loop on second orbital of mesh point (only for jc.le.ic)
          do jc = 1,ic
            jl = ilc(jc)

C           Loop over sub-points
            Vij = 0.0d0
            do isp = 1,nsp
              Vij = Vij + VClocal(isp) * Clocal(isp,jc)
            enddo

            Vlocal(jl,il) = Vlocal(jl,il) + Vij
            if (ic .ne. jc) then
              Vlocal(il,jl) = Vlocal(il,jl) + Vij
            endif

          enddo
        enddo
      enddo

C  Add final Vlocal to Vs
      do il = 1,last
        i = iorb(il)
        iu = indxuo(i)
        if (Parallel) then
          iul = NeedDscfL(iu)
          if (i .eq. iu) then
            do ii = 1, numdl(iul)
              ind = listdlptr(iul)+ii
              j = listdl(ind)
              jl = ilocal(j)
              DscfL(ind,1) = DscfL(ind,1) + dVol * Vlocal(jl,il) 
            enddo
          else
            do ii = 1, numdl(iul)
              ind = listdlptr(iul)+ii
              j = listsc( i, iu, listdl(ind) )
              jl = ilocal(j)
              DscfL(ind,1) = DscfL(ind,1) + dVol * Vlocal(jl,il) 
            enddo
          endif
        else
          if (i .eq. iu) then
            do ii = 1, numVs(iu)
              ind = listVsptr(iu)+ii
              j = listVs(ind)
              jl = ilocal(j)
              Vs(ind) = Vs(ind) + dVol * Vlocal(jl,il) 
            enddo
          else
            do ii = 1, numVs(iu)
              ind = listVsptr(iu)+ii
              j = listsc( i, iu, listVs(ind) )
              jl = ilocal(j)
              Vs(ind) = Vs(ind) + dVol * Vlocal(jl,il) 
            enddo
          endif
        endif
      enddo

C  Free local memory
      call memory('D','D',size(Vlocal),'vmat')
      deallocate(Vlocal)
      call memory('D','I',size(iorb),'vmat')
      deallocate(iorb)
      call memory('D','I',size(ilocal),'vmat')
      deallocate(ilocal)
      call memory('D','I',size(ilc),'vmat')
      deallocate(ilc)
      call memory('D','D',size(Clocal),'vmat')
      deallocate(Clocal)
      call memory('D','D',size(VClocal),'vmat')
      deallocate(VClocal)

      if (Parallel) then
C Redistribute Hamiltonian from mesh to orbital based distribution
        call matrixMtoO( nvmaxl, nvmax, numVs, listVsptr, nuo, 
     .    nuotot, DscfL, Vs )
C Free memory 
        call memory('D','D',size(DscfL),'meshdscf')
        deallocate(DscfL)
      endif

      call timer('vmat',2)
      return
      end
