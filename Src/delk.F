! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_delk
      private
      public :: delk
      CONTAINS
      subroutine delk( iexpikr, no, np, dvol, nvmax, 
     &                 numVs, listVsptr, listVs, 
     &                 nuo, nuotot, iaorb, iphorb, isa )

C ********************************************************************
C Finds the matrix elements of a plane wave
C < \phi_{\mu} | exp^( iexpikr * i * \vec{k} \cdot \vec{r} ) | \phi_{\nu} >
C First version written by J. Junquera in Feb. 2008
C Adapted from an existing version of vmat after the parallelization 
C designed by BSC in July 2011.
C *********************** INPUT **************************************
C integer iexpikr         : Prefactor of the dot product between the
C                           the k-vector and the position-vector in exponent.
C                           it might be +1 or -1
C integer no              : Number of basis orbitals
C integer np              : Number of columns in C (local)
C real*8  dvol            : Volume per mesh point
C integer nvmax           : First dimension of listV and Vs, and maxim.
C                           number of nonzero elements in any row of delkmat
C integer numVs(nuo)      : Number of non-zero elements in a row of delkmat
C integer listVsptr(nuo)  : Pointer to the start of rows in listVs
C integer listVs(nvmax)   : List of non-zero elements of delkmat
C integer iaorb(*)        : Pointer to atom to which orbital belongs
C integer iphorb(*)       : Orbital index within each atom
C integer isa(*)          : Species index of all atoms
C *****************************  OUTPUT *******************************
C complex*16  delkmat(nvmax) : Value of nonzero elements in each row
C                              of the matrix elements of exp(i*\vec{k}\vec{r})
C                              this variable is defined in the module
C                              m_dimensionsefield (file m_dimefield.f90)
C *********************************************************************

C  Modules
      use precision,     only: dp, grid_p
      use atmfuncs,      only: rcut, all_phi
      use atm_types,     only: nsmax=>nspecies
      use atomlist,      only: indxuo, indxua
      use siesta_geom,   only: na_s, xa
      use listsc_module, only: listsc
      use mesh,          only: dxa, nsp, xdop, xdsp, ne, nem,
     &                         cmesh, ipa, idop, nmsc, iatfold,
     &                         meshLim
      use meshdscf,      only: matrixMtoOC
      use meshdscf,      only: needdscfl, listdl, numdl, nrowsdscfl,
     &                         listdlptr
      use meshphi,       only: directphi, endpht, lstpht, listp2, phi
      use parallel,      only: Nodes, node
      use alloc,         only: re_alloc, de_alloc
      use parallelsubs,  only: GlobalToLocalOrb
      use m_planewavematrixvar, only: delkmat, wavevector
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

C Argument types and dimensions
      integer                  :: iexpikr, no, np, nvmax, nuo, nuotot, 
     &                            iaorb(*), iphorb(*), isa(*), 
     &                            numVs(nuo),     
     &                            listVsptr(nuo), listVs(nvmax)
      real(dp)                 :: dvol
C Internal variables and arrays
      integer,       parameter :: minloc = 100,  ! Min buffer size
     &                            maxoa  = 100   ! Max # of orb/atom
      integer                  :: i, ia, ic, ii, ijl, il, imp, ind, iop,
     &                            ip, iphi, io, is, isp, iu, iul,
     &                            ix, j, jc, jl, last, lasta, lastop,
     &                            maxloc, maxloc2, nc, nlocal, nphiloc,
     &                            nvmaxl, triang, lenx, leny, lenz,lenxy
      logical                  :: ParallelLocal
      real(dp)                 :: Vij, r2sp, dxsp(3), VClocal(nsp)

      integer,         pointer :: ilc(:),          omp_ilc(:,:),
     &                            ilocal(:),       omp_ilocal(:,:),
     &                            iorb(:),         omp_iorb(:,:)
      real(dp),        pointer :: DscfL(:,:),      omp_DscfL(:,:,:),
     &                            Clocal(:,:),     omp_Clocal(:,:,:),
     &                            Vlocal(:,:),     omp_Vlocal(:,:,:),
     &                            phia(:,:),       omp_phia(:,:,:),
     &                            r2cut(:)
      integer                  :: NTH, TID

C Variables to compute the matrix element of the plane wave
C (not included in the original vmat subroutine)
      integer                  :: irel, iua, irealim, inmp(3)
      real(dp)                 :: kxij, aux(2), dist(3), kpoint(3),
     &                            dxp(3), displaat(3)
      real(dp)                 :: dxpgrid(3,nsp), comkxij(2,nsp)
      complex(dp),     pointer :: delkmats(:),   omp_delkmats(:,:)


#ifdef _OMP_
      integer                  :: omp_get_num_threads,
     &                            omp_get_thread_num
#endif
#ifdef MPI
      integer                  :: MPIerror
#endif

#ifdef DEBUG
        call write_debug( '    PRE delk' )
#endif

#ifdef _TRACE_
      call MPI_Barrier( MPI_Comm_World, MPIerror )
      call MPItrace_event( 1000, 4 )
#endif

C     Start time counter
      call timer('delk',1)

C  Initialize the matrix elements of exp(i*\vec{k} \vec{r})
      delkmat(:)      = cmplx(0.0_dp,0.0_dp,kind=dp)
      kpoint(:)       = dble(iexpikr) * wavevector(:)

!! For debugging
!      kpoint(:) = 0.0_dp
!! End debugging

C     Find atomic cutoff radii
      nullify(r2cut)
      call re_alloc( r2cut, 1, nsmax, 'r2cut', 'delk' )
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
!$OMP   PARALLEL DEFAULT(NONE) SHARED(no,np,nvmax,iaorb,iphorb,
!$   &    isa,numVs,listVsptr,listVs,dVol,maxloc,maxloc2,
!$   &    nvmaxl,triang,lenx,leny,lenxy,parallellocal,omp_ilc,
!$   &    omp_ilocal,omp_iorb,omp_DscfL,omp_delkmats,omp_Clocal,
!$   &    omp_Vlocal,
!$   &    omp_phia,r2cut,indxuo,dxa,nsp,xdop,xdsp,NeedDscfL,listdl,
!$   &    numdl,listdlptr,directphi,endpht,lstpht,listp2,phi,nodes,node)
!$   &    PRIVATE(i,ia,ic,ii,ijl,il,imp,ind,iop,ip,iphi,is,isp,irealim,
!$   &    iu,iul,ix,j,jc,jl,last,lasta,lastop,nc,nlocal,nphiloc,Vij,
!$   &    r2sp,dxsp,VClocal,ilc,ilocal,iorb,DscfL,delkmats,Clocal,
!$   &    Vlocal,      
!$   &    phia,nth,tid)
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
        call re_alloc( omp_DscfL, 1, nvmaxl, 1, 2, 1, NTH, 
     &                 'DscfL',  'delk' )
      else
        if (NTH.gt.1) then
          nullify ( omp_delkmats )
          call re_alloc( omp_delkmats, 1, nvmax, 1, NTH, 
     &                   'delkmats',  'delk' )
        endif
      endif

      nullify ( omp_ilocal,  omp_ilc, omp_iorb, omp_Vlocal,
     &          omp_Clocal,  omp_phia )
      call re_alloc( omp_ilocal, 1, no, 1, NTH, 'ilocal', 'delk' )
      call re_alloc( omp_ilc, 1, maxloc2, 1, NTH, 'ilc', 'delk' )
      call re_alloc( omp_iorb, 1, maxloc, 1, NTH, 'iorb', 'delk' )
      call re_alloc( omp_Vlocal, 1, triang, 1, 2, 1, NTH,
     &               'Vlocal', 'delk' )
      call re_alloc( omp_Clocal, 1, nsp, 1, maxloc2, 1, NTH,
     &               'Clocal', 'delk' )
      call re_alloc( omp_phia, 1, maxoa, 1, nsp, 1, NTH,
     &               'phia', 'delk' )
!$OMP END SINGLE

      if (ParallelLocal) then
        DscfL => omp_DscfL(1:nvmaxl,1:2,TID)
        DscfL(1:nvmaxl,1:2) = 0.0_dp
      else
        if (NTH.gt.1) then
          delkmats => omp_delkmats(1:nvmax,TID)
        else
          delkmats => delkmat
        endif
      endif
      ilocal     => omp_ilocal(1:no,TID)
      ilc        => omp_ilc(1:maxloc2,TID)
      iorb       => omp_iorb(1:maxloc,TID)
      Vlocal     => omp_Vlocal(1:triang,1:2,TID)
      Clocal     => omp_Clocal(1:nsp,1:maxloc2,TID)
      phia       => omp_phia(1:maxoa,1:nsp,TID)

C     Full initializations done only once
      ilocal(1:no)             = 0
      iorb(1:maxloc)           = 0
      Vlocal(1:triang,1:2)     = 0.0_dp
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

C       If overflooded, add Vlocal to delkmat and reinitialize it
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
! The variables we want to compute in this subroutine are complex numbers
! Here, when irealim =1 we refer to the real part, and
! when irealim = 2 we refer to the imaginary part
                  do irealim = 1, 2
                    DscfL(ind,irealim) = DscfL(ind,irealim) + dVol * 
     &                Vlocal(ijl,irealim) 
                  enddo
                enddo
              else
                ia  = iaorb(i)
                iua = indxua(ia)
                do ix = 1, 3
                  displaat(ix) =
     .             (iatfold(1,ia)*nmsc(1))*cmesh(ix,1)+
     .             (iatfold(2,ia)*nmsc(2))*cmesh(ix,2)+
     .             (iatfold(3,ia)*nmsc(3))*cmesh(ix,3)
                enddo
                dist(:) = xa(:,iua) - xa(:,ia) - displaat(:)
                kxij    = kpoint(1) *  dist(1)  +
     .                    kpoint(2) *  dist(2)  +
     .                    kpoint(3) *  dist(3)
                do ii = 1, numdl(iul)
                  ind = listdlptr(iul)+ii
                  j   = listsc( i, iu, listdl(ind) )
                  jl  = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  aux(1) = Vlocal(ijl,1) * dcos(kxij) -
     .                     Vlocal(ijl,2) * dsin(kxij)
                  aux(2) = Vlocal(ijl,2) * dcos(kxij) +
     .                     Vlocal(ijl,1) * dsin(kxij)
                  do irealim = 1, 2
                    DscfL(ind,irealim) = DscfL(ind,irealim) + dVol * 
     &                aux(irealim)
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

                  delkmats(ind) = delkmats(ind) +
     &              dVol * cmplx(Vlocal(ijl,1), Vlocal(ijl,2), kind=dp)

                enddo
              else
                ia  = iaorb(i)
                iua = indxua(ia)
                do ix = 1, 3
                  displaat(ix) =
     .             (iatfold(1,ia)*nmsc(1))*cmesh(ix,1)+
     .             (iatfold(2,ia)*nmsc(2))*cmesh(ix,2)+
     .             (iatfold(3,ia)*nmsc(3))*cmesh(ix,3)
                enddo
                dist(:) = xa(:,iua) - xa(:,ia) - displaat(:)
                kxij    = kpoint(1) *  dist(1)  +
     .                    kpoint(2) *  dist(2)  +
     .                    kpoint(3) *  dist(3)
                do ii = 1, numVs(iul)
                  ind = listVsptr(iul)+ii
                  j = listsc( i, iu, listVs(ind) )
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  aux(1) = Vlocal(ijl,1) * dcos(kxij) -
     .                     Vlocal(ijl,2) * dsin(kxij)
                  aux(2) = Vlocal(ijl,2) * dcos(kxij) +
     .                     Vlocal(ijl,1) * dsin(kxij)

                  delkmats(ind) = delkmats(ind) +
     .              dVol * cmplx(aux(1),aux(2),kind=dp)

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
          Vlocal(1:ijl,1:2) = 0.0_dp
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
          imp     = endpht(ip-1) + ic
          i       = lstpht(imp)
          il      = ilocal(i)
          iu      = indxuo(i)
          ia      = iaorb(i)
          iua     = indxua(ia)
          is      = isa(ia)
          iop     = listp2(imp)
          ilc(ic) = il

C  Localize the position of the mesh point
          irel = idop(iop) + ipa(ia)
          call ipack( -1, 3, nem, inmp, irel )
          inmp(:) = inmp(:) + (meshLim(1,:) - 1)
          inmp(:) = inmp(:) - 2 * ne(:)

          do ix = 1, 3
           dxp(ix) = cmesh(ix,1) * inmp(1) +
     .               cmesh(ix,2) * inmp(2) +
     .               cmesh(ix,3) * inmp(3)
          enddo

          do isp = 1, nsp
            do ix = 1, 3
              dxpgrid(ix,isp) = dxp(ix) + xdsp(ix,isp)
            enddo
          enddo

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
          do irealim = 1, 2
            do isp = 1,nsp
              kxij           = kpoint(1) * dxpgrid(1,isp) +
     .                         kpoint(2) * dxpgrid(2,isp) +
     .                         kpoint(3) * dxpgrid(3,isp)
              comkxij(1,isp) = dcos(kxij)
              comkxij(2,isp) = dsin(kxij)
              VClocal(isp)   = comkxij(irealim,isp) * Clocal(isp,ic)
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
                Vlocal(ijl,irealim) = Vlocal(ijl,irealim) + 2.0_dp*Vij
              else
                Vlocal(ijl,irealim) = Vlocal(ijl,irealim) + Vij
              endif

            enddo
          enddo
        enddo
      enddo
!$OMP END DO

C  Add final Vlocal to delkmat
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
              do irealim = 1, 2
                DscfL(ind,irealim) = DscfL(ind,irealim) + dVol * 
     &            Vlocal(ijl,irealim) 
              enddo
            enddo
          else
            ia  = iaorb(i)
            iua = indxua(ia)
            do ix = 1, 3
              displaat(ix) =
     .          (iatfold(1,ia)*nmsc(1))*cmesh(ix,1)+
     .          (iatfold(2,ia)*nmsc(2))*cmesh(ix,2)+
     .          (iatfold(3,ia)*nmsc(3))*cmesh(ix,3)
            enddo
            dist(:) = xa(:,iua) - xa(:,ia) - displaat(:)
            kxij    = kpoint(1) * dist(1) +
     .                kpoint(2) * dist(2) +
     .                kpoint(3) * dist(3)
            do ii = 1, numdl(iul)
              ind = listdlptr(iul)+ii
              j = listsc( i, iu, listdl(ind) )
              jl = ilocal(j)
              if (il.gt.jl) then
                ijl = il*(il+1)/2 + jl + 1
              else
                ijl = jl*(jl+1)/2 + il + 1
              endif
              aux(1) = Vlocal(ijl,1) * dcos(kxij) -
     .                 Vlocal(ijl,2) * dsin(kxij)
              aux(2) = Vlocal(ijl,2) * dcos(kxij) +
     .                 Vlocal(ijl,1) * dsin(kxij)
              do irealim = 1, 2
                DscfL(ind,irealim) = DscfL(ind,irealim) + dVol * 
     &            aux(irealim) 
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
              delkmats(ind) = delkmats(ind) +
     .           dVol * cmplx(Vlocal(ijl,1), Vlocal(ijl,2),kind=dp)
            enddo
          else
            ia  = iaorb(i)
            iua = indxua(ia)
            do ix = 1, 3
              displaat(ix) =
     .          (iatfold(1,ia)*nmsc(1))*cmesh(ix,1)+
     .          (iatfold(2,ia)*nmsc(2))*cmesh(ix,2)+
     .          (iatfold(3,ia)*nmsc(3))*cmesh(ix,3)
            enddo
            dist(:) = xa(:,iua) - xa(:,ia) - displaat(:)
            kxij    = kpoint(1) * dist(1) +
     .                kpoint(2) * dist(2) +
     .                kpoint(3) * dist(3)
            do ii = 1, numVs(iu)
              ind = listVsptr(iu)+ii
              j = listsc( i, iu, listVs(ind) )
              jl = ilocal(j)
              if (il.gt.jl) then
                ijl = il*(il+1)/2 + jl + 1
              else
                ijl = jl*(jl+1)/2 + il + 1
              endif
              aux(1) = Vlocal(ijl,1) * dcos(kxij) -
     .                 Vlocal(ijl,2) * dsin(kxij)
              aux(2) = Vlocal(ijl,2) * dcos(kxij) +
     .                 Vlocal(ijl,1) * dsin(kxij)
              delkmats(ind) = delkmats(ind) +
     .           dVol * cmplx(aux(1),aux(2),kind=dp)
            enddo
          endif
        endif
      enddo

!$OMP BARRIER

      if (ParallelLocal) then
!$OMP   DO SCHEDULE(STATIC)
        do ind = 1, nvmaxl
          do irealim = 1, 2
            do ii = 2, NTH
              omp_DscfL(ind,irealim,1) = omp_DscfL(ind,irealim,1) +
     &                                   omp_DscfL(ind,irealim,ii)
            enddo
          enddo
        enddo
!$OMP   END DO
      else if (NTH.gt.1) then
!$OMP   DO SCHEDULE(STATIC)
        do ind = 1, nvmax
          do ii = 1, NTH
            delkmat(ind) = delkmat(ind) + 
     &                     omp_delkmats(ind,ii)
          enddo
        enddo
!$OMP   END DO
      endif

!$OMP END PARALLEL

C     Free memory
      call de_alloc( omp_phia, 'phia', 'delk' )
      call de_alloc( omp_Clocal, 'Clocal', 'delk' )
      call de_alloc( omp_Vlocal, 'Vlocal', 'delk' )
      call de_alloc( omp_iorb, 'iorb', 'delk' )
      call de_alloc( omp_ilc, 'ilc', 'delk' )
      call de_alloc( omp_ilocal, 'ilocal', 'delk' )

      if (ParallelLocal) then
C       Redistribute Hamiltonian from mesh to orbital based distribution
        DscfL => omp_DscfL(1:nvmaxl,1:2,1)
        call matrixMtoOC( nvmaxl, nvmax, numVs, listVsptr, nuo, 
     &                    DscfL, delkmat )
        call de_alloc( omp_DscfL, 'DscfL', 'delk' )
      else if (NTH.gt.1) then
        call de_alloc( omp_delkmats, 'delkmats', 'delk' )
      endif

      call de_alloc( r2cut, 'r2cut', 'delk' )

#ifdef _TRACE_
      call MPI_Barrier( MPI_Comm_World, MPIerror )
      call MPItrace_event( 1000, 0 )
#endif

      call timer('delk',2)

#ifdef DEBUG
      call write_debug( '    POS delk' )
#endif
      return
      end subroutine delk
      end module m_delk
