      subroutine rhoofd( no, indxuo, np, maxnd, numd, listdptr, 
     .                   listd, Dscf, rhoscf, nuo, nuotot, 
     .                   iaorb, iphorb, isa )
C ********************************************************************
C Finds the SCF density at the mesh points from the density matrix.
C Written by P.Ordejon and J.M.Soler. May'95.
C Re-ordered so that mesh is the outer loop and the orbitals are
C handled as lower-half triangular. J.D.Gale and J.M.Soler, Feb'99
C Version of rhoofd that optionally uses a direct algorithm to save 
C memory. Modified by J.D.Gale, November'99
C *********************** InpUT **************************************
C integer no              : Number of basis orbitals
C integer indxuo(no)      : Index of equivalent atom in unit cell
C integer np              : Number of mesh points
C integer maxnd           : First dimension of listD and Dscf, and
C                           maximum number of nonzero elements in
C                           any row of Dscf
C integer numd(nuo)       : Number of nonzero elemts in each row of Dscf
C integer listdptr(nuo)   : Pointer to start of rows in listd
C integer listd(maxnd)    : List of nonzero elements in each row of Dscf
C real*8  Dscf(maxnd)     : Rows of Dscf that are non-zero 
C integer nuo             : Number of orbitals in unit cell locally
C integer nuotot          : Number of orbitals in unit cell in total
C integer iaorb(*)        : Pointer to atom to which orbital belongs
C integer iphorb(*)       : Orbital index within each atom
C integer isa(*)          : Species index of all atoms
C *********************** OUTPUT **************************************
C real    rhoscf(nsp,np)  : SCF density at mesh points
C *********************************************************************

C  Modules
      use precision
      use atmfuncs, only: rcut, phiatm, all_phi, nsmax=>nspecies
      use listsc_module, only: listsc
      use mesh,     only: nsp, dxa, xdop, xdsp
      use meshdscf
      use meshphi

      implicit none

C Argument types and dimensions
      integer
     .   no, np, maxnd, nuo, nuotot, indxuo(no), iaorb(*),
     .   iphorb(*), isa(*), numd(nuo), listdptr(nuo), listd(maxnd)

      real
     .   rhoscf(nsp,np)

      real*8
     .   Ci, Cj, Dscf(maxnd)

      external
     .   memory, timer

C Internal variables and arrays
      integer, parameter ::
     .  minloc = 100,  ! Min buffer size for local copy of Dscf
     .  maxoa  = 100   ! Max # of orbitals per atom

      integer
     .  i, ia, ic, ii, il, imp, in, ind,
     .  iop, ip, iphi, is, isp, iu, iul, ix,
     .  j, jc, jl, jn, jmp,
     .  last, lasta, lastop, maxloc, maxloc2, maxndl, nc, nphiloc

      integer, dimension(:), allocatable, save :: 
     .  ilc, ilocal, iorb

      logical
     .  Parallel

      real*8
     .  Dij, dxsp(3), phia(maxoa,nsp), r2cut(nsmax), r2sp

      real*8, dimension(:,:), allocatable ::
     .  Clocal, Dlocal

C  Start time counter
      call timer('rhoofd',1)

C  Set algorithm logical
      Parallel = (nuo .ne. nuotot)

C  Find size of buffers to store partial copies of Dscf and C
      maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
      maxloc = maxloc2 + minloc
      maxloc = min( maxloc, no )

C  Allocate local memory
      allocate(ilocal(no))
      call memory('A','I',no,'rhoofd')
      allocate(ilc(maxloc2))
      call memory('A','I',maxloc2,'rhoofd')
      allocate(iorb(0:maxloc))
      call memory('A','I',maxloc+1,'rhoofd')
      allocate(Dlocal(0:maxloc,0:maxloc))
      call memory('A','D',(maxloc+1)*(maxloc+1),'rhoofd')
      allocate(Clocal(nsp,maxloc2))
      call memory('A','D',nsp*maxloc2,'rhoofd')

      if (Parallel) then
        maxndl = listdlptr(nrowsDscfL)+numdl(nrowsDscfL)
        allocate(DscfL(maxndl,1))
        call memory('A','D',maxndl,'meshdscf')
C Redistribute Dscf to DscfL form
        call matrixOtoM( maxnd, numd, listdptr, maxndl, nuo, nuotot,
     .    Dscf, DscfL )
      endif

C  Find atomic cutoff radiae
      do i = 1,no
        ia = iaorb(i)
        is = isa(ia)
        r2cut(is) = rcut(is,0)**2
      enddo

C  Initializations
      rhoscf(:,:) = 0.0
      Dlocal(:,:) = 0.0d0
      ilocal(:) = 0
      iorb(:) = 0
      last = 0

C  Loop over grid points
      do ip = 1,np

C  Find number of nonzero orbitals at this point
        nc = endpht(ip) - endpht(ip-1)

C  iorb(il)>0 means that row il of Dlocal must not be overwritten
C  iorb(il)=0 means that row il of Dlocal is empty
C  iorb(il)<0 means that row il of Dlocal contains a valid row of 
C             Dscf, but which is not required at this point
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          il = ilocal(i)
          if (il.gt.0) iorb(il) = i
        enddo

C  Look for required rows of Dscf not yet stored in Dlocal
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
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
            if (Parallel) then
              iul = NeedDscfL(iu)
              if (i .eq. iu) then
                do ii = 1, numdl(iul)
                  ind = listdlptr(iul)+ii
                  j = listdl(ind)
                  jl = ilocal(j)
                  Dij = DscfL(ind,1)
                  Dlocal(il,jl) = Dij
                  Dlocal(jl,il) = Dij
                enddo
              else
                do ii = 1, numdl(iul)
                  ind = listdlptr(iul)+ii
                  j = listsc( i, iu, listdl(ind) )
                  jl = ilocal(j)
                  Dij = DscfL(ind,1)
                  Dlocal(il,jl) = Dij
                  Dlocal(jl,il) = Dij
                enddo
              endif
            else
              if (i .eq. iu) then
                do ii = 1, numd(iu)
                  ind = listdptr(iu)+ii
                  j = listd(ind)
                  jl = ilocal(j)
                  Dij = Dscf(ind)
                  Dlocal(il,jl) = Dij
                  Dlocal(jl,il) = Dij
                enddo
              else
                do ii = 1, numd(iu)
                  ind = listdptr(iu)+ii
                  j = listsc( i, iu, listd(ind) )
                  jl = ilocal(j)
                  Dij = Dscf(ind)
                  Dlocal(il,jl) = Dij
                  Dlocal(jl,il) = Dij
                enddo
              endif
            endif
          endif
        enddo

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
                dxsp(:) = xdsp(:,isp) + xdop(:,iop) - dxa(:,ia)
                r2sp = sum(dxsp**2)
                if (r2sp.lt.r2cut(is)) then
                  call all_phi( is, +1, dxsp, nphiloc, phia(:,isp) )
                else
                  phia(:,isp) = 0.0d0
                endif
              enddo
            endif
            iphi = iphorb(i)
            Clocal(:,ic) = phia(iphi,:)
          else
            Clocal(:,ic) = phi(:,imp)
          endif

C  Loop on second orbital of mesh point
          do jc = 1,ic
            jl = ilc(jc)
            if (ic .eq. jc) then
              Dij = Dlocal(il,jl)
            else
              Dij = 2*Dlocal(il,jl)
            endif

C  Loop over sub-points
            do isp = 1, nsp
              rhoscf(isp,ip) = rhoscf(isp,ip) +
     .                         Dij * Clocal(isp,ic) * Clocal(isp,jc)
            enddo

          enddo

        enddo

C  Restore iorb for next point
        do imp = 1+endpht(ip-1), endpht(ip)
          i = lstpht(imp)
          il = ilocal(i)
          iorb(il) = -i
        enddo

      enddo

C  Free local memory
      call memory('D','I',size(ilocal),'rhoofd')
      deallocate(ilocal)
      call memory('D','I',size(ilc),'rhoofd')
      deallocate(ilc)
      call memory('D','I',size(iorb),'rhoofd')
      deallocate(iorb)
      call memory('D','D',size(Dlocal),'rhoofd')
      deallocate(Dlocal)
      call memory('D','D',size(Clocal),'rhoofd')
      deallocate(Clocal)

      if (Parallel) then
        call memory('D','D',size(DscfL),'meshdscf')
        deallocate(DscfL)
      endif

      call timer('rhoofd',2)
      end
