      subroutine vmatsp( no, indxuo, np, dvol, nspin, V, nvmax, 
     .                 numVs, listVsptr, listVs, Vs, 
     .                 nuo, nuotot, iaorb, iphorb, isa, q )

C ********************************************************************
C Finds the matrix elements of the potential.
C First version written by P.Ordejon.
C Name and interface modified by J.M.Soler. May'95.
C Re-ordered so that mesh is the outer loop and the orbitals are
C handled as lower-half triangular. J.D.Gale and J.M.Soler, Feb'99
C Version of vmat that use a direct algorithm to save memory.
C Modified by J.D.Gale, November'99
C Spiral version written by V. M. Garcia-Suarez. June 2002.
C *********************** INPUT **************************************
C integer no              : Number of basis orbitals
C integer indxuo(no)      : Index of equivalent atom in unit cell
C integer np              : Number of columns in C (local)
C real*8  dvol            : Volume per mesh point
C integer nspin           : Number of spin components
C real*4  V(np,nspin)     : Value of the potential at the mesh points
C integer nvmax           : First dimension of listV and Vs, and maxim.
C                           number of nonzero elements in any row of Vs
C integer numVs(nuo)      : Number of non-zero elements in a row of Vs
C integer listVsptr(nuo)  : Pointer to the start of rows in listVs
C integer listVs(nvmax)   : List of non-zero elements of Vs
C integer iaorb(*)        : Pointer to atom to which orbital belongs
C integer iphorb(*)       : Orbital index within each atom
C integer isa(*)          : Species index of all atoms
C real*8 q(3)             : Wave vector for spiral configuration.
C ******************** INPUT AND OUTPUT *******************************
C real*8  Vs(nvmax,nspin) : Value of nonzero elements in each row 
C                           of Vs to which the potential matrix 
C                           elements are summed up
C *********************************************************************

C  Modules
      use precision
      use atmfuncs, only: rcut, phiatm, all_phi, nsmax=>nspecies
      use listsc_module, only: listsc
      use mesh, only: dxa, nsp, xdop, xdsp, cmesh, nmeshg, nsm
      use meshdscf
      use meshphi

      implicit none

C Argument types and dimensions
      integer
     .   no, np, nvmax, nuo, nuotot, indxuo(no), iaorb(*), nspin,
     .   iphorb(*), isa(*), numVs(nuo), listVsptr(nuo), listVs(nvmax)
      real
     .   V(nsp,np,nspin)
      real*8
     .   dvol, Vs(nvmax,nspin), q(3)
      external
     .   memory, timer, ipack

C Internal variables and arrays

      integer, parameter ::
     .  minloc = 100,  ! Min buffer size for local copy of Dscf
     .  maxoa  = 100   ! Max # of orbitals per atom
      integer
     .  i, ia, iua, ic, ii, ijl, il, imp, ind, iop, ip, iphi, io,
     .  is, isp, ispin, iu, iul, ix, j, jc, jl, iv, nsd,
     .  last, lasta, lastop, maxloc, maxloc2, nc, nlocal, 
     .  nphiloc, nvmaxl, iii(3)
      integer, dimension(:), allocatable, save ::
     .  ilc, ilocal, iorb
      logical
     .  Parallel
      real*8
     .  Vij, Vij1, Vij2, Vij3, Vij4, dxsp(3), phia(maxoa,nsp),
     .  r2cut(nsmax), r2sp, mod, xr(3), Rdi(3), qRdi, cqRdi, sqRdi
      real*8, dimension(:),   allocatable, save :: 
     .  VClocal, VClocal1, VClocal2, VClocal3, VClocal4    
      real*8, dimension(:,:), allocatable, save :: 
     .  Clocal, Vlocal
      real*8, dimension(:,:,:), allocatable, save ::
     .  VlocalSp

C  Start time counter
      call timer('vmatsp',1)

C  Set algorithm logical
      Parallel = (nuo .ne. nuotot)

C  Find value of maxloc

      maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
      maxloc = maxloc2 + minloc
      maxloc = min( maxloc, no )

C  If spiral, the diagonal elements of Vlocal do not change
      nsd = 2

C  Allocate local memory
      allocate(ilocal(no))
      call memory('A','I',no,'vmatsp')
      allocate(ilc(maxloc2))
      call memory('A','I',maxloc2,'vmatsp')
      allocate(iorb(0:maxloc))
      call memory('A','I',maxloc+1,'vmatsp')
      ijl = (maxloc+1)*(maxloc+2)/2
      allocate(Vlocal(ijl,nsd))
      call memory('A','D',ijl*nsd,'vmatsp')
      allocate(VlocalSp(0:maxloc,0:maxloc,nsd))
      call memory('A','D',(maxloc+1)*(maxloc+1)*nsd,'vmatsp')
      allocate(Clocal(nsp,maxloc2))
      call memory('A','D',nsp*maxloc2,'vmatsp')
      allocate(VClocal(nsp))
      call memory('A','D',nsp,'vmatsp')
      allocate(VClocal1(nsp))
      call memory('A','D',nsp,'vmatsp')
      allocate(VClocal2(nsp))
      call memory('A','D',nsp,'vmatsp')
      allocate(VClocal3(nsp))
      call memory('A','D',nsp,'vmatsp')
      allocate(VClocal4(nsp))
      call memory('A','D',nsp,'vmatsp')

      if (Parallel) then
        nvmaxl = listdlptr(nrowsDscfL) + numdl(nrowsDscfL)
        allocate(DscfL(nvmaxl,nsd+2))
        call memory('A','D',nvmaxl*(nsd+2),'meshdscf')
        DscfL(1:nvmaxl,1:nsd+2) = 0.0d0
      endif

C  Full initializations done only once
      ilocal(1:no) = 0
      iorb(0:maxloc) = 0
      Vlocal(:,:) = 0.0
      VlocalSp(:,:,:) = 0.0
      last = 0

C  Find atomic cutoff radiae
      r2cut(:) = 0.0d0
      do i = 1,nuotot
        ia = iaorb(i)
        is = isa(ia)
        io = iphorb(i)
        r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
      enddo

C  Loop over grid points
      do ip = 1,np

C  Find point coordinates
        call ipack(-1,3,nmeshg/nsm,iii,ip)
        do ix = 1,3
          xr(ix) = iii(1) * cmesh(ix,1) + iii(2)*cmesh(ix,2) +
     .             iii(3) * cmesh(ix,3)
        enddo

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
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nsd
                    DscfL(ind,ispin) = DscfL(ind,ispin) + dVol * 
     .                Vlocal(ijl,ispin)
                    DscfL(ind,ispin+2) = DscfL(ind,ispin+2) + dVol * 
     .                VlocalSp(jl,il,ispin)
                  enddo
                enddo
              else
                do ii = 1, numdl(iul)
                  ind = listdlptr(iul)+ii
                  j = listsc( i, iu, listdl(ind) )
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nsd
                    DscfL(ind,ispin) = DscfL(ind,ispin) + dVol * 
     .                Vlocal(ijl,ispin)
                    DscfL(ind,ispin+2) = DscfL(ind,ispin+2) + dVol *
     .                VlocalSp(jl,il,ispin)
                  enddo
                enddo
              endif
            else
              if (i .eq. iu) then
                do ii = 1, numVs(iu)
                  ind = listVsptr(iu)+ii
                  j = listVs(ind)
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nsd
                    Vs(ind,ispin) = Vs(ind,ispin) + dVol * 
     .                Vlocal(ijl,ispin)
                    Vs(ind,ispin+2) = Vs(ind,ispin+2) + dVol *
     .                VlocalSp(jl,il,ispin)
                  enddo
                enddo
              else
                do ii = 1, numVs(iu)
                  ind = listVsptr(iu)+ii
                  j = listsc( i, iu, listVs(ind) )
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nsd
                    Vs(ind,ispin) = Vs(ind,ispin) + dVol * 
     .                Vlocal(ijl,ispin)
                    Vs(ind,ispin+2) = Vs(ind,ispin+2) + dVol *
     .                VlocalSp(jl,il,ispin)
                  enddo
                enddo
              endif
            endif
          enddo
          ilocal(iorb(1:last)) = 0
          iorb(1:last) = 0
          ijl = (last+1)*(last+2)/2
          Vlocal(1:ijl,1:nsd) = 0.0d0
          VlocalSp(1:last,1:last,1:nsd) = 0.0d0
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

        lasta=0
        lastop=0

C  Generate or retrieve phi values for all orbitals up to nc
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          il = ilocal(i)
          iu = indxuo(i)
          ia = iaorb(i)
          is = isa(ia)
          iop = listp2(imp)
          ilc(ic) = il

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
        enddo

C  Loop on first orbital of mesh point
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          il = ilocal(i)
          iu = indxuo(i)
          ia = iaorb(i)
          is = isa(ia)
          iop = listp2(imp)

C  Calculate spiral phase
          Rdi(1:3) =  xr(1:3) - xdop(1:3,iop) + dxa(1:3,ia)
          qRdi = q(1) * Rdi(1) + q(2) * Rdi(2) + q(3) * Rdi(3)
          cqRdi = cos(qRdi)
          sqRdi = sin(qRdi)

C  Pre-multiply V and Clocal(,ic)
          do isp = 1,nsp
            VClocal1(isp) = V(isp,ip,1) * Clocal(isp,ic)
            VClocal2(isp) = V(isp,ip,2) * Clocal(isp,ic)
            VClocal3(isp) = (V(isp,ip,3)*cqRdi + V(isp,ip,4)*sqRdi)
     .                      *Clocal(isp,ic)
            VClocal4(isp) = (V(isp,ip,4)*cqRdi - V(isp,ip,3)*sqRdi)
     .                      *Clocal(isp,ic)
          enddo

C  Loop on second orbital of mesh point (NOT only for jc.le.ic)
          do jc = 1,nc
            jl = ilc(jc)
  
C           Loop over sub-points
            Vij1 = 0.0d0
            Vij2 = 0.0d0
            Vij3 = 0.0d0
            Vij4 = 0.0d0
            do isp = 1,nsp
              Vij1 = Vij1 + VClocal1(isp) * Clocal(isp,jc)
              Vij2 = Vij2 + VClocal2(isp) * Clocal(isp,jc)
              Vij3 = Vij3 + VClocal3(isp) * Clocal(isp,jc)
              Vij4 = Vij4 + VClocal4(isp) * Clocal(isp,jc)
            enddo

            if (il.gt.jl) then
              ijl = il*(il+1)/2 + jl + 1
            else
              ijl = jl*(jl+1)/2 + il + 1
            endif

            if (jc.le.ic) then
              if (ic.ne.jc.and.il.eq.jl) then
                Vlocal(ijl,1) = Vlocal(ijl,1) + 2.0*Vij1
                Vlocal(ijl,2) = Vlocal(ijl,2) + 2.0*Vij2
              else
                Vlocal(ijl,1) = Vlocal(ijl,1) + Vij1
                Vlocal(ijl,2) = Vlocal(ijl,2) + Vij2
              endif
            endif

            VlocalSp(jl,il,1) = VlocalSp(jl,il,1) + Vij3
            VlocalSp(jl,il,2) = VlocalSp(jl,il,2) + Vij4

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
              if (il.gt.jl) then
                ijl = il*(il+1)/2 + jl + 1
              else
                ijl = jl*(jl+1)/2 + il + 1
              endif
              do ispin = 1,nsd
                DscfL(ind,ispin) = DscfL(ind,ispin) + dVol * 
     .            Vlocal(ijl,ispin) 
                DscfL(ind,ispin+2) = DscfL(ind,ispin+2) + dVol *
     .            VlocalSp(jl,il,ispin)
              enddo
            enddo
          else
            do ii = 1, numdl(iul)
              ind = listdlptr(iul)+ii
              j = listsc( i, iu, listdl(ind) )
              jl = ilocal(j)
              if (il.gt.jl) then
                ijl = il*(il+1)/2 + jl + 1
              else
                ijl = jl*(jl+1)/2 + il + 1
              endif
              do ispin = 1,nsd
                DscfL(ind,ispin) = DscfL(ind,ispin) + dVol * 
     .            Vlocal(ijl,ispin) 
                DscfL(ind,ispin+2) = DscfL(ind,ispin+2) + dVol *
     .            VlocalSp(jl,il,ispin)
              enddo
            enddo
          endif
        else
          if (i .eq. iu) then
            do ii = 1, numVs(iu)
              ind = listVsptr(iu)+ii
              j = listVs(ind)
              jl = ilocal(j)
              if (il.gt.jl) then
                ijl = il*(il+1)/2 + jl + 1
              else
                ijl = jl*(jl+1)/2 + il + 1
              endif
              do ispin = 1,nsd
                Vs(ind,ispin) = Vs(ind,ispin) + dVol * 
     .            Vlocal(ijl,ispin)
                Vs(ind,ispin+2) = Vs(ind,ispin+2) + dVol *
     .            VlocalSp(jl,il,ispin)
              enddo
            enddo
          else
            do ii = 1, numVs(iu)
              ind = listVsptr(iu)+ii
              j = listsc( i, iu, listVs(ind) )
              jl = ilocal(j)
              if (il.gt.jl) then
                ijl = il*(il+1)/2 + jl + 1
              else
                ijl = jl*(jl+1)/2 + il + 1
              endif
              do ispin = 1,nsd
                Vs(ind,ispin) = Vs(ind,ispin) + dVol * 
     .            Vlocal(ijl,ispin)
                Vs(ind,ispin+2) = Vs(ind,ispin+2) + dVol *
     .            VlocalSp(jl,il,ispin)
              enddo
            enddo
          endif
        endif
      enddo

C  Free local memory
      call memory('D','D',size(Vlocal),'vmatsp')
      deallocate(Vlocal)
      call memory('D','D',size(VlocalSp),'vmatsp')
      deallocate(VlocalSp)
      call memory('D','I',size(iorb),'vmatsp')
      deallocate(iorb)
      call memory('D','I',size(ilocal),'vmatsp')
      deallocate(ilocal)
      call memory('D','I',size(ilc),'vmatsp')
      deallocate(ilc)
      call memory('D','D',size(Clocal),'vmatsp')
      deallocate(Clocal)
      call memory('D','D',size(VClocal),'vmatsp')
      deallocate(VClocal)
      call memory('D','D',size(VClocal1),'vmatsp')
      deallocate(VClocal1)
      call memory('D','D',size(VClocal2),'vmatsp')
      deallocate(VClocal2)
      call memory('D','D',size(VClocal3),'vmatsp')
      deallocate(VClocal3)
      call memory('D','D',size(VClocal4),'vmatsp')
      deallocate(VClocal4)

      if (Parallel) then
C Redistribute Hamiltonian from mesh to orbital based distribution
        call matrixMtoO( nvmaxl, nvmax, numVs, listVsptr, nuo, 
     .      nuotot, nspin, DscfL, Vs )
C Free memory 
        call memory('D','D',size(DscfL),'meshdscf')
        deallocate(DscfL)
      endif

      call timer('vmatsp',2)
      return
      end
