! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_vmat
      private
      public :: vmat
      CONTAINS
      subroutine vmat( no, np, dvol, nspin, V, nvmax, 
     &                 numVs, listVsptr, listVs, Vs, 
     &                 nuo, nuotot, iaorb, iphorb, isa )

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
C ******************** INPUT AND OUTPUT *******************************
C real*8  Vs(nvmax,nspin) : Value of nonzero elements in each row 
C                           of Vs to which the potential matrix 
C                           elements are summed up
C *********************************************************************

C  Modules
      use precision,     only: dp, grid_p
      use atmfuncs,      only: rcut, all_phi
      use atm_types,     only: nsmax=>nspecies
      use atomlist,      only: indxuo
      use listsc_module, only: LISTSC
      use mesh,          only: dxa, nsp, xdop, xdsp, meshLim
      use meshdscf,      only: matrixMtoO
      use meshdscf,      only: needdscfl, listdl, numdl, nrowsdscfl,
     &                         listdlptr
      use meshphi,       only: directphi, endpht, lstpht, listp2, phi
      use parallel,      only: Nodes, node
      use alloc,         only: re_alloc, de_alloc
      use parallelsubs,  only: GlobalToLocalOrb
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

C Argument types and dimensions
      integer                  :: no, np, nvmax, nuo, nuotot, iaorb(*),
     &                            nspin, iphorb(*), isa(*), numVs(nuo),
     &                            listVsptr(nuo), listVs(nvmax)
      real(grid_p), intent(in) :: V(nsp,np,nspin)
      real(dp)                 :: dvol
      real(dp),         target :: Vs(nvmax,nspin)
C Internal variables and arrays
      integer,       parameter :: minloc = 1000,  ! Min buffer size
     &                            maxoa  = 100   ! Max # of orb/atom
      integer                  :: i, ia, ic, ii, ijl, il, imp, ind, iop,
     &                            ip, iphi, io, is, isp, ispin, iu, iul,
     &                            ix, j, jc, jl, last, lasta, lastop,
     &                            maxloc, maxloc2, nc, nlocal, nphiloc,
     &                            nvmaxl, triang, lenx, leny, lenz,lenxy
      logical                  :: ParallelLocal
      real(dp)                 :: Vij, r2sp, dxsp(3), VClocal(nsp)

      integer,         pointer :: ilc(:),    omp_ilc(:,:),
     &                            ilocal(:), omp_ilocal(:,:),
     &                            iorb(:),   omp_iorb(:,:)
      real(dp),        pointer :: DscfL(:,:), omp_DscfL(:,:,:),
     &                            Vss(:,:), omp_Vss(:,:,:),
     &                            Clocal(:,:), omp_Clocal(:,:,:),
     &                            Vlocal(:,:), omp_Vlocal(:,:,:),
     &                            phia(:,:), omp_phia(:,:,:),
     &                            r2cut(:)
      integer                  :: NTH, TID
#ifdef _OMP_
      integer                  :: omp_get_num_threads,
     &                            omp_get_thread_num
#endif
#ifdef MPI
      integer                  :: MPIerror
#endif

#ifdef DEBUG
        call write_debug( '    PRE vmat' )
#endif

#ifdef _TRACE_
      call MPI_Barrier( MPI_Comm_World, MPIerror )
      call MPItrace_event( 1000, 4 )
#endif

C       Start time counter
      call timer('vmat',1)

C     Find atomic cutoff radii
      nullify(r2cut)
      call re_alloc( r2cut, 1, nsmax, 'r2cut', 'vmat' )
      r2cut = 0.0_dp
      do i = 1,nuotot
        ia = iaorb(i)
        is = isa(ia)
        io = iphorb(i)
        r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
      enddo

C     Set algorithm logical
      ParallelLocal = (Nodes.gt.1)
      lenx  = meshLim(2,1) - meshLim(1,1) + 1
      leny  = meshLim(2,2) - meshLim(1,2) + 1
      lenz  = meshLim(2,3) - meshLim(1,3) + 1
      lenxy = lenx*leny

C     Find value of maxloc
      maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
      maxloc  = maxloc2 + minloc
      maxloc  = min( maxloc, no )
      triang  = (maxloc+1)*(maxloc+2)/2
      if (ParallelLocal) then
        if (nrowsDscfL.gt.0) then
          nvmaxl = listdlptr(nrowsDscfL) + numdl(nrowsDscfL)
        else
          nvmaxl = 1
        endif
      endif

C     Allocate local memory
!$OMP   PARALLEL DEFAULT(NONE) SHARED(no,np,nvmax,iaorb,nspin,iphorb,
!$OMP&    isa,numVs,listVsptr,listVs,V,dVol,Vs,maxloc,maxloc2,
!$OMP&    nvmaxl,triang,lenx,leny,lenxy,parallellocal,omp_ilc,
!$OMP&    omp_ilocal,omp_iorb,omp_DscfL,omp_Vss,omp_Clocal,omp_Vlocal,
!$OMP&    omp_phia,r2cut,indxuo,dxa,nsp,xdop,xdsp,NeedDscfL,listdl,
!$OMP&    numdl,listdlptr,directphi,endpht,lstpht,listp2,phi,nodes,node)
!$OMP&    PRIVATE(i,ia,ic,ii,ijl,il,imp,ind,iop,ip,iphi,is,isp,ispin,
!$OMP&    iu,iul,ix,j,jc,jl,last,lasta,lastop,nc,nlocal,nphiloc,Vij,
!$OMP&    r2sp,dxsp,VClocal,ilc,ilocal,iorb,DscfL,Vss,Clocal,Vlocal,
!$OMP&    phia,nth,tid)
#ifdef _OMP_
      NTH = omp_get_num_threads( )
      TID = omp_get_thread_num( ) + 1
#else
      NTH = 1
      TID = 1
#endif

!$OMP SINGLE
      if (ParallelLocal) then
        nullify ( omp_DscfL )
        call re_alloc( omp_DscfL, 1, nvmaxl, 1, nspin, 1, NTH, 
     &                 'DscfL',  'vmat' )
      else
        if (NTH.gt.1) then
          nullify ( omp_Vss )
          call re_alloc( omp_Vss, 1, nvmax, 1, nspin, 1, NTH, 
     &                   'Vss',  'vmat' )
        endif
      endif

      nullify ( omp_ilocal, omp_ilc, omp_iorb, omp_Vlocal,
     &          omp_Clocal, omp_phia )
      call re_alloc( omp_ilocal, 1, no, 1, NTH, 'ilocal', 'vmat' )
      call re_alloc( omp_ilc, 1, maxloc2, 1, NTH, 'ilc', 'vmat' )
      call re_alloc( omp_iorb, 1, maxloc, 1, NTH, 'iorb', 'vmat' )
      call re_alloc( omp_Vlocal, 1, triang, 1, nspin, 1, NTH,
     &               'Vlocal', 'vmat' )
      call re_alloc( omp_Clocal, 1, nsp, 1, maxloc2, 1, NTH,
     &               'Clocal', 'vmat' )
      call re_alloc( omp_phia, 1, maxoa, 1, nsp, 1, NTH,
     &               'phia', 'vmat' )
!$OMP END SINGLE

      if (ParallelLocal) then
        DscfL => omp_DscfL(1:nvmaxl,1:nspin,TID)
        DscfL(1:nvmaxl,1:nspin) = 0.0_dp
      else
        if (NTH.gt.1) then
          Vss => omp_Vss(1:nvmax,1:nspin,TID)
        else
          Vss => Vs
        endif
      endif
      ilocal   => omp_ilocal(1:no,TID)
      ilc      => omp_ilc(1:maxloc2,TID)
      iorb     => omp_iorb(1:maxloc,TID)
      Vlocal   => omp_Vlocal(1:triang,1:nspin,TID)
      Clocal   => omp_Clocal(1:nsp,1:maxloc2,TID)
      phia     => omp_phia(1:maxoa,1:nsp,TID)

C     Full initializations done only once
      ilocal(1:no)             = 0
      iorb(1:maxloc)           = 0
      Vlocal(1:triang,1:nspin) = 0.0_dp
      last                     = 0

C     Loop over grid points
!$OMP DO SCHEDULE(GUIDED,LENXY)
      do ip = 1,np
C       Find number of nonzero orbitals at this point
        nc = endpht(ip) - endpht(ip-1)
C       Find new required size of Vlocal
        nlocal = last
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          if (ilocal(i) .eq. 0) nlocal = nlocal + 1
        enddo

C       If overflooded, add Vlocal to Vs and reinitialize it
        if (nlocal .gt. maxloc) then
          do il = 1,last
            i = iorb(il)
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
                    DscfL(ind,ispin) = DscfL(ind,ispin) + dVol * 
     &                Vlocal(ijl,ispin) 
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
                    DscfL(ind,ispin) = DscfL(ind,ispin) + dVol * 
     &                Vlocal(ijl,ispin) 
                  enddo
                enddo
              endif
            else
              call GlobalToLocalOrb( iu, Node, Nodes, iul )
              if (i .eq. iu) then
                do ii = 1, numVs(iul)
                  ind = listVsptr(iul)+ii
                  j = listVs(ind)
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nspin
                    Vss(ind,ispin) = Vss(ind,ispin) + dVol * 
     &                Vlocal(ijl,ispin) 
                  enddo
                enddo
              else
                do ii = 1, numVs(iul)
                  ind = listVsptr(iul)+ii
                  j = LISTSC( i, iu, listVs(ind) )
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nspin
                    Vss(ind,ispin) = Vss(ind,ispin) + dVol * 
     &                Vlocal(ijl,ispin) 
                  enddo
                enddo
              endif
            endif
          enddo
!         Reset local arrays
          do ii= 1, last
            ilocal(iorb(ii)) = 0
          enddo
          iorb(1:last) = 0
          ijl = (last+1)*(last+2)/2
          Vlocal(1:ijl,1:nspin) = 0.0_dp
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
                do ix = 1,3
                  dxsp(ix) = xdsp(ix,isp) + xdop(ix,iop) - dxa(ix,ia)
                enddo
                r2sp = sum(dxsp**2)
                if (r2sp.lt.r2cut(is)) then
                  call all_phi( is, +1, dxsp, nphiloc, phia(:,isp) )
                else
                  phia(:,isp) = 0.0_dp
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

C         Pre-multiply V and Clocal(,ic)
          do ispin = 1,nspin
            do isp = 1,nsp
              VClocal(isp) = V(isp,ip,ispin) * Clocal(isp,ic)
            enddo

C           Loop on second orbital of mesh point (only for jc.le.ic)
            do jc = 1,ic
              jl = ilc(jc)

C             Loop over sub-points
              Vij = 0.0_dp
              do isp = 1,nsp
                Vij = Vij + VClocal(isp) * Clocal(isp,jc)
              enddo

              if (il.gt.jl) then
                ijl = il*(il+1)/2 + jl + 1
              else
                ijl = jl*(jl+1)/2 + il + 1
              endif
              if (ic.ne.jc.and.il.eq.jl) then
                Vlocal(ijl,ispin) = Vlocal(ijl,ispin) + 2.0_dp*Vij
              else
                Vlocal(ijl,ispin) = Vlocal(ijl,ispin) + Vij
              endif

            enddo
          enddo
        enddo
      enddo
!$OMP END DO

C  Add final Vlocal to Vs
      do il = 1,last
        i = iorb(il)
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
                DscfL(ind,ispin) = DscfL(ind,ispin) + dVol * 
     &            Vlocal(ijl,ispin) 
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
                DscfL(ind,ispin) = DscfL(ind,ispin) + dVol * 
     &            Vlocal(ijl,ispin) 
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
              do ispin = 1,nspin
                Vss(ind,ispin) = Vss(ind,ispin) + dVol * 
     &            Vlocal(ijl,ispin) 
              enddo
            enddo
          else
            do ii = 1, numVs(iu)
              ind = listVsptr(iu)+ii
              j = LISTSC( i, iu, listVs(ind) )
              jl = ilocal(j)
              if (il.gt.jl) then
                ijl = il*(il+1)/2 + jl + 1
              else
                ijl = jl*(jl+1)/2 + il + 1
              endif
              do ispin = 1,nspin
                Vss(ind,ispin) = Vss(ind,ispin) + dVol * 
     &            Vlocal(ijl,ispin) 
              enddo
            enddo
          endif
        endif
      enddo

!$OMP BARRIER

      if (ParallelLocal) then
!$OMP   DO SCHEDULE(STATIC)
        do ind = 1, nvmaxl
          do ispin = 1,nspin
            do ii = 2, NTH
              omp_DscfL(ind,ispin,1) = omp_DscfL(ind,ispin,1) +
     &                                 omp_DscfL(ind,ispin,ii)
            enddo
          enddo
        enddo
!$OMP   END DO
      else if (NTH.gt.1) then
!$OMP   DO SCHEDULE(STATIC)
        do ind = 1, nvmax
          do ispin = 1,nspin
            do ii = 1, NTH
              Vs(ind,ispin) = Vs(ind,ispin) + omp_Vss(ind,ispin,ii)
            enddo
          enddo
        enddo
!$OMP   END DO
      endif

!$OMP END PARALLEL

C     Free memory
      call de_alloc( omp_phia, 'phia', 'vmat' )
      call de_alloc( omp_Clocal, 'Clocal', 'vmat' )
      call de_alloc( omp_Vlocal, 'Vlocal', 'vmat' )
      call de_alloc( omp_iorb, 'iorb', 'vmat' )
      call de_alloc( omp_ilc, 'ilc', 'vmat' )
      call de_alloc( omp_ilocal, 'ilocal', 'vmat' )

      if (ParallelLocal) then
C       Redistribute Hamiltonian from mesh to orbital based distribution
        DscfL => omp_DscfL(1:nvmaxl,1:nspin,1)
        call matrixMtoO( nvmaxl, nvmax, numVs, listVsptr, nuo, 
     &                   nspin, DscfL, Vs )
        call de_alloc( omp_DscfL, 'DscfL', 'vmat' )
      else if (NTH.gt.1) then
        call de_alloc( omp_Vss, 'Vss', 'vmat' )
      endif

      call de_alloc( r2cut, 'r2cut', 'vmat' )

#ifdef _TRACE_
      call MPI_Barrier( MPI_Comm_World, MPIerror )
      call MPItrace_event( 1000, 0 )
#endif

      call timer('vmat',2)

#ifdef DEBUG
      call write_debug( '    POS vmat' )
#endif
      return
      end subroutine vmat
      end module m_vmat
