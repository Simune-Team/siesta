! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
module m_delk

  implicit none

  private

  public :: delk

contains

  subroutine delk( iexpikr, no, np, dvol, nvmax, &
       numVs, listVsptr, listVs, &
       nuo, nuotot, iaorb, iphorb, isa )

! ********************************************************************
! Finds the matrix elements of a plane wave
! < \phi_{\mu} | exp^( iexpikr * i * \vec{k} \cdot \vec{r} ) | \phi_{\nu} >
! First version written by J. Junquera in Feb. 2008
! Adapted from an existing version of vmat after the parallelization 
! designed by BSC in July 2011.
! *********************** INPUT **************************************
! integer iexpikr         : Prefactor of the dot product between the
!                           the k-vector and the position-vector in exponent.
!                           it might be +1 or -1
! integer no              : Number of basis orbitals
! integer np              : Number of columns in C (local)
! real*8  dvol            : Volume per mesh point
! integer nvmax           : First dimension of listV and Vs, and maxim.
!                           number of nonzero elements in any row of delkmat
! integer numVs(nuo)      : Number of non-zero elements in a row of delkmat
! integer listVsptr(nuo)  : Pointer to the start of rows in listVs
! integer listVs(nvmax)   : List of non-zero elements of delkmat
! integer iaorb(*)        : Pointer to atom to which orbital belongs
! integer iphorb(*)       : Orbital index within each atom
! integer isa(*)          : Species index of all atoms
! *****************************  OUTPUT *******************************
! complex*16  delkmat(nvmax) : Value of nonzero elements in each row
!                              of the matrix elements of exp(i*\vec{k}\vec{r})
!                              this variable is defined in the module
!                              m_dimensionsefield (file m_dimefield.f90)
! *********************************************************************

!  Modules
    use precision,     only: dp, grid_p
    use atmfuncs,      only: rcut, all_phi
    use atm_types,     only: nsmax=>nspecies
    use atomlist,      only: indxuo, indxua
    use siesta_geom,   only: xa
    use listsc_module, only: listsc
    use mesh,          only: dxa, nsp, xdop, xdsp, ne, nem
    use mesh,          only: cmesh, ipa, idop, nmsc, iatfold
    use mesh,          only: meshLim
    use meshdscf,      only: matrixMtoOC
    use meshdscf,      only: needdscfl, listdl, numdl, nrowsdscfl, listdlptr
    use meshphi,       only: directphi, endpht, lstpht, listp2, phi
    use parallel,      only: Nodes, node
    use alloc,         only: re_alloc, de_alloc
    use parallelsubs,  only: GlobalToLocalOrb
    use m_planewavematrixvar, only: delkmat, wavevector
#ifdef MPI
    use mpi_siesta
#endif
#ifdef _OPENMP
    use omp_lib
#endif

! Argument types and dimensions
    integer :: iexpikr, no, np, nvmax, nuo, nuotot
    integer :: iaorb(*), iphorb(*), isa(*), numVs(nuo)
    integer :: listVsptr(nuo), listVs(nvmax)
    real(dp) :: dvol
! Internal variables and arrays
    integer, parameter :: minloc = 1000 ! Min buffer size
    integer, parameter :: maxoa  = 100  ! Max # of orb/atom
    integer :: i, ia, ic, ii, ijl, il, imp, ind, iop
    integer :: ip, iphi, io, is, isp, iu, iul
    integer :: ix, j, jc, jl, last, lasta, lastop
    integer :: maxloc, maxloc2, nc, nlocal, nphiloc
    integer :: nvmaxl, triang, lenx, leny, lenz,lenxy
    logical :: ParallelLocal
    real(dp) :: Vij, r2sp, dxsp(3), VClocal(nsp)

    integer, pointer :: ilc(:), ilocal(:), iorb(:)
    real(dp), pointer :: DscfL(:,:),  t_DscfL(:,:,:), Clocal(:,:)
    real(dp), pointer :: Vlocal(:,:), phia(:,:), r2cut(:)
    integer :: NTH, TID

! Variables to compute the matrix element of the plane wave
! (not included in the original vmat subroutine)
    integer :: irel, iua, irealim, inmp(3)
    real(dp) :: kxij, aux(2), dist(3), kpoint(3)
    real(dp) :: dxp(3), displaat(3)
    real(dp) :: dxpgrid(3,nsp), comkxij(2,nsp)
    complex(dp), pointer :: delkmats(:), t_delkmats(:,:)

#ifdef _TRACE_
    integer :: MPIerror
#endif

#ifdef DEBUG
    call write_debug( '    PRE delk' )
#endif

#ifdef _TRACE_
    call MPI_Barrier( MPI_Comm_World, MPIerror )
    call MPItrace_event( 1000, 4 )
#endif

!   Start time counter
    call timer('delk',1)

!   Initialize the matrix elements of exp(i*\vec{k} \vec{r})
    delkmat(:) = 0.0_dp
    kpoint(:)  = dble(iexpikr) * wavevector(:)

!! For debugging
!      kpoint(:) = 0.0_dp
!! End debugging

!   Find atomic cutoff radii
    nullify(r2cut)
    call re_alloc( r2cut, 1, nsmax, 'r2cut', 'delk' )
    r2cut = 0.0_dp
    do i = 1,nuotot
       ia = iaorb(i)
       is = isa(ia)
       io = iphorb(i)
       r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
    enddo

!   Set algorithm logical
    ParallelLocal = (Nodes > 1)
    lenx  = meshLim(2,1) - meshLim(1,1) + 1
    leny  = meshLim(2,2) - meshLim(1,2) + 1
    lenz  = meshLim(2,3) - meshLim(1,3) + 1
    lenxy = lenx*leny

!   Find value of maxloc
    maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
    maxloc  = maxloc2 + minloc
    maxloc  = min( maxloc, no )
    triang  = (maxloc+1)*(maxloc+2)/2
    if ( ParallelLocal ) then
       if ( nrowsDscfL > 0 ) then
          nvmaxl = listdlptr(nrowsDscfL) + numdl(nrowsDscfL)
       else
          nvmaxl = 1
       endif
    endif

!   Allocate local memory
!$OMP parallel default(shared), &
!$OMP&shared(NTH,t_DscfL,t_delkmats), &
!$OMP&private(TID,last,delkmats,irealim), &
!$OMP&private(ip,nc,nlocal,ic,imp,i,il,iu,iul,ii,ind,j,ijl,jl), &
!$OMP&private(lasta,lastop,ia,is,iop,isp,ix,dxsp,r2sp,nphiloc,iphi,jc), &
!$OMP&private(Vij,VClocal,DscfL,ilocal,ilc,iorb,Vlocal,Clocal,phia), &
!$OMP&private(irel,inmp,dxp,dxpgrid,dist,kxij,iua,displaat,aux)

!$OMP single
#ifdef _OPENMP
    NTH = omp_get_num_threads( )
#else
    NTH = 1
#endif
!$OMP end single ! implicit barrier, IMPORTANT

#ifdef _OPENMP
    TID = omp_get_thread_num( ) + 1
#else
    TID = 1
#endif

    nullify(Clocal,phia,ilocal,ilc,iorb,Vlocal)
!$OMP critical
! Perhaps the critical section is not needed,
! however it "tells" the OS to allocate per
! thread, possibly waiting for each thread to
! place the memory in the best position.
    allocate( Clocal(nsp,maxloc2) )
    allocate( phia(maxoa,nsp) )
    allocate( ilocal(no) , ilc(maxloc2) , iorb(maxloc) )
    allocate( Vlocal(triang,2) )
!$OMP end critical

!$OMP single
    if ( ParallelLocal ) then
       nullify( t_DscfL )
       call re_alloc( t_DscfL, 1, nvmaxl, 1, 2, 1, NTH, &
            'DscfL',  'delk' )
    else
       if ( NTH > 1 ) then
          nullify( t_delkmats )
          call re_alloc( t_delkmats, 1, nvmax, 1, NTH, &
               'delkmats',  'delk' )
       endif
    endif
!$OMP end single

    if ( ParallelLocal ) then
       DscfL => t_DscfL(1:nvmaxl,1:2,TID)
       DscfL(1:nvmaxl,1:2) = 0.0_dp
    else
       if ( NTH > 1 ) then
          delkmats => t_delkmats(1:nvmax,TID)
       else
          delkmats => delkmat
       endif
    endif

!   Full initializations done only once
    ilocal(1:no)             = 0
    iorb(1:maxloc)           = 0
    Vlocal(1:triang,1:2)     = 0.0_dp
    last                     = 0

!   Loop over grid points
!$OMP do
    do ip = 1,np
!      Find number of nonzero orbitals at this point
       nc = endpht(ip) - endpht(ip-1)
!      Find new required size of Vlocal
       nlocal = last
       do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          if (ilocal(i) .eq. 0) nlocal = nlocal + 1
       enddo

!      If overflooded, add Vlocal to delkmat and reinitialize it
       if (nlocal > maxloc .and. last > 0) then
          do il = 1,last
             i = iorb(il)
             iu = indxuo(i)
             if ( ParallelLocal ) then
                iul = NeedDscfL(iu)
                if (i .eq. iu) then
                   do ii = 1, numdl(iul)
                      ind = listdlptr(iul)+ii
                      j = listdl(ind)
                      jl = ilocal(j)
                      ijl = idx_ijl(il,jl)
! The variables we want to compute in this subroutine are complex numbers
! Here, when irealim =1 we refer to the real part, and
! when irealim = 2 we refer to the imaginary part
                      do irealim = 1, 2
                         DscfL(ind,irealim) = DscfL(ind,irealim) + dVol * &
                              Vlocal(ijl,irealim) 
                      enddo
                   enddo
                else
                   ia  = iaorb(i)
                   iua = indxua(ia)
                   do ix = 1, 3
                      displaat(ix) = &
                           (iatfold(1,ia)*nmsc(1))*cmesh(ix,1)+ &
                           (iatfold(2,ia)*nmsc(2))*cmesh(ix,2)+ &
                           (iatfold(3,ia)*nmsc(3))*cmesh(ix,3)
                   enddo
                   dist(:) = xa(:,iua) - xa(:,ia) - displaat(:)
                   kxij    = kpoint(1) *  dist(1)  + &
                        kpoint(2) *  dist(2)  + &
                        kpoint(3) *  dist(3)
                   do ii = 1, numdl(iul)
                      ind = listdlptr(iul)+ii
                      j   = listsc( i, iu, listdl(ind) )
                      jl  = ilocal(j)
                      ijl = idx_ijl(il,jl)
                      aux(1) = Vlocal(ijl,1) * dcos(kxij) - &
                           Vlocal(ijl,2) * dsin(kxij)
                      aux(2) = Vlocal(ijl,2) * dcos(kxij) + &
                           Vlocal(ijl,1) * dsin(kxij)
                      do irealim = 1, 2
                         DscfL(ind,irealim) = DscfL(ind,irealim) + dVol * &
                              aux(irealim)
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
                      ijl = idx_ijl(il,jl)
                      delkmats(ind) = delkmats(ind) + &
                           dVol * cmplx(Vlocal(ijl,1), Vlocal(ijl,2), kind=dp)

                   enddo
                else
                   ia  = iaorb(i)
                   iua = indxua(ia)
                   do ix = 1, 3
                      displaat(ix) = &
                           (iatfold(1,ia)*nmsc(1))*cmesh(ix,1)+ &
                           (iatfold(2,ia)*nmsc(2))*cmesh(ix,2)+ &
                           (iatfold(3,ia)*nmsc(3))*cmesh(ix,3)
                   enddo
                   dist(:) = xa(:,iua) - xa(:,ia) - displaat(:)
                   kxij    = kpoint(1) *  dist(1)  + &
                        kpoint(2) *  dist(2)  + &
                        kpoint(3) *  dist(3)
                   do ii = 1, numVs(iul)
                      ind = listVsptr(iul)+ii
                      j = listsc( i, iu, listVs(ind) )
                      jl = ilocal(j)
                      ijl = idx_ijl(il,jl)
                      aux(1) = Vlocal(ijl,1) * dcos(kxij) - &
                           Vlocal(ijl,2) * dsin(kxij)
                      aux(2) = Vlocal(ijl,2) * dcos(kxij) + &
                           Vlocal(ijl,1) * dsin(kxij)

                      delkmats(ind) = delkmats(ind) + &
                           dVol * cmplx(aux(1),aux(2),kind=dp)

                   enddo
                endif
             endif
          enddo
!         Reset local arrays
          do ii = 1, last
             ilocal(iorb(ii)) = 0
          enddo
          iorb(1:last) = 0
          ijl = (last+1)*(last+2)/2
          Vlocal(1:ijl,1:2) = 0.0_dp
          last = 0
       endif

!      Look for required orbitals not yet in Vlocal
       if (nlocal > last) then
          do ic = 1,nc
             imp = endpht(ip-1) + ic
             i = lstpht(imp)
             if (ilocal(i) == 0) then
                last = last + 1
                ilocal(i) = last
                iorb(last) = i
             endif
          enddo
       endif

!      Loop on first orbital of mesh point
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

!         Localize the position of the mesh point
          irel = idop(iop) + ipa(ia)
          call ipack( -1, 3, nem, inmp, irel )
          inmp(:) = inmp(:) + (meshLim(1,:) - 1)
          inmp(:) = inmp(:) - 2 * ne(:)

          do ix = 1, 3
             dxp(ix) = cmesh(ix,1) * inmp(1) + &
                  cmesh(ix,2) * inmp(2) + &
                  cmesh(ix,3) * inmp(3)
          enddo

          do isp = 1, nsp
             do ix = 1, 3
                dxpgrid(ix,isp) = dxp(ix) + xdsp(ix,isp)
             enddo
          enddo

!         Generate or retrieve phi values
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

!         Pre-multiply V and Clocal(,ic)
          do irealim = 1, 2
             do isp = 1,nsp
                kxij           = kpoint(1) * dxpgrid(1,isp) + &
                     kpoint(2) * dxpgrid(2,isp) + &
                     kpoint(3) * dxpgrid(3,isp)
                comkxij(1,isp) = dcos(kxij)
                comkxij(2,isp) = dsin(kxij)
                VClocal(isp)   = comkxij(irealim,isp) * Clocal(isp,ic)
             enddo

!            Loop on second orbital of mesh point (only for jc.le.ic)
             do jc = 1,ic
                jl = ilc(jc)

!               Loop over sub-points
                Vij = 0.0_dp
                do isp = 1,nsp
                   Vij = Vij + VClocal(isp) * Clocal(isp,jc)
                enddo

                ijl = idx_ijl(il,jl)
                if (ic.ne.jc.and.il.eq.jl) then
                   Vlocal(ijl,irealim) = Vlocal(ijl,irealim) + 2.0_dp*Vij
                else
                   Vlocal(ijl,irealim) = Vlocal(ijl,irealim) + Vij
                endif

             enddo
          enddo
       enddo
    enddo
!$OMP end do nowait

!   Add final Vlocal to delkmat
    if ( ParallelLocal .and. last > 0 ) then
       do il = 1 , last
          i = iorb(il)
          iu = indxuo(i)
          if (ParallelLocal) then
             iul = NeedDscfL(iu)
             if (i .eq. iu) then
                do ii = 1, numdl(iul)
                   ind = listdlptr(iul)+ii
                   j = listdl(ind)
                   jl = ilocal(j)
                   ijl = idx_ijl(il,jl)
                   do irealim = 1, 2
                      DscfL(ind,irealim) = DscfL(ind,irealim) + &
                           dVol * Vlocal(ijl,irealim) 
                   enddo
                enddo
             else
                ia  = iaorb(i)
                iua = indxua(ia)
                do ix = 1, 3
                   displaat(ix) = &
                        (iatfold(1,ia)*nmsc(1))*cmesh(ix,1)+ &
                        (iatfold(2,ia)*nmsc(2))*cmesh(ix,2)+ &
                        (iatfold(3,ia)*nmsc(3))*cmesh(ix,3)
                enddo
                dist(:) = xa(:,iua) - xa(:,ia) - displaat(:)
                kxij    = kpoint(1) * dist(1) + &
                     kpoint(2) * dist(2) + &
                     kpoint(3) * dist(3)
                do ii = 1, numdl(iul)
                   ind = listdlptr(iul)+ii
                   j = listsc( i, iu, listdl(ind) )
                   jl = ilocal(j)
                   ijl = idx_ijl(il,jl)
                   aux(1) = Vlocal(ijl,1) * dcos(kxij) - &
                        Vlocal(ijl,2) * dsin(kxij)
                   aux(2) = Vlocal(ijl,2) * dcos(kxij) + &
                        Vlocal(ijl,1) * dsin(kxij)
                   do irealim = 1, 2
                      DscfL(ind,irealim) = DscfL(ind,irealim) + &
                           dVol * aux(irealim) 
                   enddo
                enddo
             endif
          else
             if (i .eq. iu) then
                do ii = 1, numVs(iu)
                   ind = listVsptr(iu)+ii
                   j = listVs(ind)
                   jl = ilocal(j)
                   ijl = idx_ijl(il,jl)
                   delkmats(ind) = delkmats(ind) + &
                        dVol * cmplx(Vlocal(ijl,1), Vlocal(ijl,2),kind=dp)
                enddo
             else
                ia  = iaorb(i)
                iua = indxua(ia)
                do ix = 1, 3
                   displaat(ix) = &
                        (iatfold(1,ia)*nmsc(1))*cmesh(ix,1)+ &
                        (iatfold(2,ia)*nmsc(2))*cmesh(ix,2)+ &
                        (iatfold(3,ia)*nmsc(3))*cmesh(ix,3)
                enddo
                dist(:) = xa(:,iua) - xa(:,ia) - displaat(:)
                kxij    = kpoint(1) * dist(1) + &
                     kpoint(2) * dist(2) + &
                     kpoint(3) * dist(3)
                do ii = 1, numVs(iu)
                   ind = listVsptr(iu)+ii
                   j = listsc( i, iu, listVs(ind) )
                   jl = ilocal(j)
                   ijl = idx_ijl(il,jl)
                   aux(1) = Vlocal(ijl,1) * dcos(kxij) - &
                        Vlocal(ijl,2) * dsin(kxij)
                   aux(2) = Vlocal(ijl,2) * dcos(kxij) + &
                        Vlocal(ijl,1) * dsin(kxij)
                   delkmats(ind) = delkmats(ind) + &
                        dVol * cmplx(aux(1),aux(2),kind=dp)
                enddo
             endif
          endif
       enddo

    end if

!$OMP barrier

    if ( ParallelLocal .and. NTH > 1 ) then
!$OMP do collapse(2)
       do ind = 1, nvmaxl
          do irealim = 1, 2
             do ii = 2, NTH
                t_DscfL(ind,irealim,1) = t_DscfL(ind,irealim,1) + &
                     t_DscfL(ind,irealim,ii)
             enddo
          enddo
       enddo
!$OMP end do
    else if ( NTH > 1 ) then
!$OMP do
       do ind = 1, nvmax
          do ii = 1, NTH
             delkmat(ind) = delkmat(ind) + t_delkmats(ind,ii)
          enddo
       enddo
!$OMP end do
    endif

!   Free local memory
    deallocate(Clocal,phia,ilocal,ilc,iorb,Vlocal)

!$OMP master
    if ( ParallelLocal ) then
!      Redistribute Hamiltonian from mesh to orbital based distribution
       DscfL => t_DscfL(1:nvmaxl,1:2,1)
       call matrixMtoOC( nvmaxl, nvmax, numVs, listVsptr, nuo, DscfL, delkmat )
       call de_alloc( t_DscfL, 'DscfL', 'delk' )
    else if ( NTH > 1 ) then
       call de_alloc( t_delkmats, 'delkmats', 'delk' )
    endif
!$OMP end master

!$OMP end parallel

    call de_alloc( r2cut, 'r2cut', 'delk' )

#ifdef _TRACE_
    call MPI_Barrier( MPI_Comm_World, MPIerror )
    call MPItrace_event( 1000, 0 )
#endif

    call timer('delk',2)

#ifdef DEBUG
    call write_debug( '    POS delk' )
#endif

  contains

!   In any case will the compiler most likely inline this
!   small routine. So it should not pose any problem.
    pure function idx_ijl(i,j) result(ij)
      integer, intent(in) :: i,j
      integer :: ij
      if ( i > j ) then
         ij = i * (i + 1)/2 + j + 1
      else
         ij = j * (j + 1)/2 + i + 1
      end if
    end function idx_ijl

  end subroutine delk

end module m_delk
