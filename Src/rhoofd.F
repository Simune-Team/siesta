! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine rhoofd( no, np, maxnd, numd, listdptr, listd, nspin, 
     &                   Dscf, rhoscf, nuo, nuotot, iaorb, iphorb, isa )
C ********************************************************************
C Finds the SCF density at the mesh points from the density matrix.
C Written by P.Ordejon and J.M.Soler. May'95.
C Re-ordered so that mesh is the outer loop and the orbitals are
C handled as lower-half triangular. J.D.Gale and J.M.Soler, Feb'99
C Version of rhoofd that optionally uses a direct algorithm to save 
C memory. Modified by J.D.Gale, November'99
C *********************** InpUT **************************************
C integer no              : Number of basis orbitals
C integer np              : Number of mesh points
C integer maxnd           : First dimension of listD and Dscf, and
C                           maximum number of nonzero elements in
C                           any row of Dscf
C integer numd(nuo)       : Number of nonzero elemts in each row of Dscf
C integer listdptr(nuo)   : Pointer to start of rows in listd
C integer listd(maxnd)    : List of nonzero elements in each row of Dscf
C integer nspin           : Number of spin components
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

      use precision,     only: dp, grid_p
      use atmfuncs,      only: rcut, all_phi
      use atm_types,     only: nsmax=>nspecies
      use atomlist,      only: indxuo
      use listsc_module, only: LISTSC
      use mesh,          only: nsp, dxa, xdop, xdsp, meshLim
      use meshdscf,      only: matrixOtoM
      use meshdscf,      only: nrowsDscfL, listdl, listdlptr, NeedDscfL,
     &                         numdl, DscfL
      use meshphi,       only: DirectPhi, endpht, lstpht, listp2, phi
      use parallel,      only: Nodes, node
      use sys,           only: die
      use alloc,         only: re_alloc, de_alloc
      use parallelsubs,  only: GlobalToLocalOrb
#ifdef MPI
      use mpi_siesta
#endif
      implicit none

C     Argument types and dimensions
      integer,       intent(in) :: no, np, maxnd, nuo, numd(nuo), nspin,
     &                             nuotot, iaorb(*), iphorb(*), 
     &                             isa(*), listdptr(nuo), listd(maxnd)
      real(dp),      intent(in) :: Dscf(maxnd,nspin)
      real(grid_p), intent(out) :: rhoscf(nsp,np,nspin)
      external                  :: memory, timer
C     Internal variables and arrays
      integer,        parameter :: minloc = 1000,  ! Min buffer size for local copy of Dscf
     &                             maxoa  = 100   ! Max # of orbitals per atom
      logical                   :: ParallelLocal
      integer                   :: i, ia, ic, ii, ijl, il, imp, ind,
     &                             ispin, io, iop, ip, iphi, is,
     &                             isp, iu, iul, j, jc, jl, last, lasta,
     &                             lastop, maxloc, maxloc2, triang, nc,
     &                             maxndl, nphiloc, lenx, leny, lenxy,
     &                             lenz
      real(dp)                  :: Dij, r2sp, dxsp(3), Cij(nsp)
      integer,          pointer :: ilc(:), ilocal(:), iorb(:)
      real(dp),         pointer :: r2cut(:), Clocal(:,:),
     &                             Dlocal(:,:), phia(:,:)
#ifdef MPI
      integer                   :: MPIerror
#endif

#ifdef DEBUG
      call write_debug( '    PRE rhoofd' )
#endif
#ifdef _TRACE_
      call MPI_Barrier( MPI_Comm_World, MPIerror )
      call MPItrace_event( 1000, 1 )
#endif
C     Start time counter
      call timer('rhoofd',1)

C     Set algorithm logical
      ParallelLocal = (Nodes.gt.1)

      if (ParallelLocal) then
        if (nrowsDscfL.gt.0) then
          maxndl = listdlptr(nrowsDscfL) + numdl(nrowsDscfL)
        else
          maxndl = 1
        endif
        nullify(DscfL)
        call re_alloc( DscfL, 1, maxndl, 1, nspin, 'DscfL', 'rhoofd' )
C       Redistribute Dscf to DscfL form
        call matrixOtoM( maxnd, numd, listdptr, maxndl, nuo,
     &                   nspin, Dscf, DscfL )
      endif

C     Find atomic cutoff radii
      nullify(r2cut)
      call re_alloc( r2cut, 1, nsmax, 'r2cut', 'rhoofd' )
      r2cut = 0.0
      do i = 1,nuotot
        ia = iaorb(i)
        is = isa(ia)
        io = iphorb(i)
        r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
      enddo

C     Find size of buffers to store partial copies of Dscf and C
      maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
      maxloc = maxloc2 + minloc
      maxloc = min( maxloc, no )
      triang = (maxloc+1)*(maxloc+2)/2

      lenx  = meshLim(2,1) - meshLim(1,1) + 1
      leny  = meshLim(2,2) - meshLim(1,2) + 1
      lenz  = meshLim(2,3) - meshLim(1,3) + 1
      lenxy = lenx*leny

!$OMP PARALLEL DEFAULT(NONE) SHARED(indxuo,nsp,dxa,xdop,xdsp,lenxy,lenx,
!$OMP&  listdl,listdlptr,NeedDscfL,numdl,DscfL,directphi,endpht,lstpht,
!$OMP&  listp2,phi,no,np,nspin,iaorb,iphorb,isa,numd,listdptr,listd,
!$OMP&  dscf,maxloc,maxloc2,triang,ParallelLocal,r2cut,rhoscf,
!$OMP&  nodes,node)
!$OMP&  PRIVATE(ip,i,ia,ic,ii,ijl,il,imp,ind,ispin,io,iop,iphi,is,isp,
!$OMP&  iu,iul,j,jc,jl,last,lasta,lastop,nc,nphiloc,Dij,r2sp,ilc,ilocal,
!$OMP&  iorb,Cij,Clocal,Dlocal,phia,dxsp)

!$OMP CRITICAL(MEMORY_LOCK)
C     Allocate local memory
      nullify ( ilocal, ilc, iorb, Dlocal, Clocal, phia )
!      call re_alloc( ilocal, 1, no, 'ilocal', 'rhoofd' )
!      call re_alloc( ilc, 1, maxloc2, 'ilc', 'rhoofd' )
!      call re_alloc( iorb, 1, maxloc, 'iorb', 'rhoofd' )
!      call re_alloc( Cij, 1, nsp, 'Cij', 'rhoofd' )
!      call re_alloc( Dlocal, 1, triang, 1, nspin, 'Dlocal', 'rhoofd' )
!      call re_alloc( Clocal, 1, nsp, 1, maxloc2, 'Clocal', 'rhoofd' )
!      call re_alloc( phia, 1, , maxoa, 1, nsp, 'phia', 'rhoofd' )
!      call re_alloc( dxsp, 1, , 3, 'dxsp', 'rhoofd' )
      allocate( ilocal(no), ilc(maxloc2), iorb(1:maxloc),
     &          Dlocal(triang,nspin), Clocal(nsp,maxloc2),
     &          phia(maxoa,nsp) )
!$OMP END CRITICAL(MEMORY_LOCK)


C     Initializations
      Dlocal = 0.0_dp
      ilocal = 0
      iorb   = 0
      last   = 0

!$OMP DO SCHEDULE(GUIDED,LENXY)
      do ip = 1,np
C       Initializations
        rhoscf(:,ip,:) = 0.0_grid_p

C       Find number of nonzero orbitals at this point
        nc = endpht(ip) - endpht(ip-1)
C       iorb(il)>0 means that row il of Dlocal must not be overwritten
C       iorb(il)=0 means that row il of Dlocal is empty
C       iorb(il)<0 means that row il of Dlocal contains a valid row of 
C             Dscf, but which is not required at this point
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          il = ilocal(i)
          if (il.gt.0) iorb(il) = i
        enddo

C       Look for required rows of Dscf not yet stored in Dlocal
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
            call die('rhoofd: no slot available in Dlocal')
   10       continue
C           Copy row i of Dscf into row last of Dlocal
            j = abs(iorb(last))
            if (j.ne.0) ilocal(j) = 0
            ilocal(i) = last
            iorb(last) = i
            il = last
            iu = indxuo(i)
            if (ParallelLocal) then
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
                  do ispin = 1,nspin
                    Dij = DscfL(ind,ispin)
                    Dlocal(ijl,ispin) = Dij
                  enddo
                enddo
              else
                do ii = 1, numdl(iul)
                  ind = listdlptr(iul)+ii
                  j = LISTSC( i, iu, listdl(ind) )
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nspin
                    Dij = DscfL(ind,ispin)
                    Dlocal(ijl,ispin) = Dij
                  enddo
                enddo
              endif
            else
              call GlobalToLocalOrb( iu, Node, Nodes, iul )
              if (i .eq. iu) then
                do ii = 1, numd(iul)
                  ind = listdptr(iul)+ii
                  j = listd(ind)
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nspin
                    Dij = Dscf(ind,ispin)
                    Dlocal(ijl,ispin) = Dij
                  enddo
                enddo
              else
                do ii = 1, numd(iul)
                  ind = listdptr(iul)+ii
                  j = LISTSC( i, iu, listd(ind) )
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nspin
                    Dij = Dscf(ind,ispin)
                    Dlocal(ijl,ispin) = Dij
                  enddo
                enddo
              endif
            endif
          endif
        enddo

C       Loop on first orbital of mesh point
        lasta = 0
        lastop = 0
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          il = ilocal(i)
          iu = indxuo(i)
          ia = iaorb(i)
          is = isa(ia)
          iop = listp2(imp)
          ilc(ic) = il

C         Generate or retrieve phi values
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
                  phia(:,isp) = 0.0_dp
                endif
              enddo
            endif
            iphi = iphorb(i)
            Clocal(:,ic) = phia(iphi,:)
          else
            Clocal(:,ic) = phi(:,imp)
          endif

C         Loop on second orbital of mesh point
          do jc = 1,ic
            jl = ilc(jc)
            do isp = 1, nsp
              Cij(isp) = Clocal(isp,ic) * Clocal(isp,jc)
            enddo
            if (il.gt.jl) then
              ijl = il*(il+1)/2 + jl + 1
            else
              ijl = jl*(jl+1)/2 + il + 1
            endif
            do ispin = 1,nspin
              if (ic .eq. jc) then
                Dij = Dlocal(ijl,ispin)
              else
                Dij = 2*Dlocal(ijl,ispin)
              endif

C             Loop over sub-points
              do isp = 1,nsp
                rhoscf(isp,ip,ispin) = rhoscf(isp,ip,ispin) + 
     &            Dij*Cij(isp)
              enddo

            enddo

          enddo

        enddo

C       Restore iorb for next point
        do imp = 1+endpht(ip-1), endpht(ip)
          i = lstpht(imp)
          il = ilocal(i)
          iorb(il) = -i
        enddo
      enddo
!$OMP END DO

!$OMP CRITICAL(MEMORY_LOCK)
C     Free local memory
!      call de_alloc( phia, 'phia', 'rhoofd' )
!      call de_alloc( Clocal, 'Clocal', 'rhoofd' )
!      call de_alloc( Dlocal, 'Dlocal', 'rhoofd' )
!      call de_alloc( iorb, 'iorb', 'rhoofd' )
!      call de_alloc( Cij, 'Cij', 'rhoofd' )
!      call de_alloc( ilc, 'ilc', 'rhoofd' )
!      call de_alloc( ilocal, 'ilocal', 'rhoofd' )
!      call de_alloc( dxsp, 'dxsp', 'rhoofd' )
      deallocate( ilocal, ilc, iorb, Dlocal, Clocal, phia )
!$OMP END CRITICAL(MEMORY_LOCK)

!$OMP END PARALLEL
      call de_alloc( r2cut, 'r2cut', 'rhoofd' )
      if (ParallelLocal) then
        call de_alloc( DscfL, 'DscfL', 'rhoofd' )
      endif

#ifdef _TRACE_
      call MPI_Barrier( MPI_Comm_World, MPIerror )
      call MPItrace_event( 1000, 0 )
#endif
      call timer('rhoofd',2)

#ifdef DEBUG
      call write_debug( '    POS rhoofd' )
#endif
      end subroutine rhoofd
