      subroutine PhiOnMesh( nmpl, norb, iaorb, iphorb, isa, numphi )
C
C Calculates the values of the orbitals at the mesh points
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer nmpl          : Number of mesh points in unit cell locally
C integer norb          : Total number of basis orbitals in supercell
C integer iaorb(norb)   : Atom to which each orbital belongs
C integer iphorb(norb)  : Orbital index (within atom)
C integer isa(na)       : Species index of all atoms in supercell
C ----------------------------------------------------------------------
C Output :
C ----------------------------------------------------------------------
C All output quantities are in the module meshphi
C ----------------------------------------------------------------------

C Modules
      use precision, only : dp, grid_p
      use atmfuncs,  only : rcut, phiatm
      use mesh,      only : dxa, idop, indexp, ipa, mop, nsp, xdop,
     &                      xdsp, meshLim
      use meshphi,   only : directphi, nphi, phi, lstpht, listp2, endpht
      use fdf,       only : fdf_get
      use parallel,  only : IONode
      use alloc,     only : re_alloc

      implicit none
C Passed arguments
      integer,    intent(in) :: nmpl, norb, iaorb(norb),
     &                          iphorb(norb), isa(*)

      ! This array was set in distriPhiOnMesh, and holds
      ! the pre-computed count of orbitals touching a given
      ! point. It is re-used and overwritten here.
      integer, intent(inout) :: numphi(nmpl)

C     Local variables
      integer                :: ia, io, iop, ip, ip0, iphi, is, isp, n,
     &                          nlist, nliste, nmax

      logical                :: within
      logical,          save :: firsttime = .true.
      real(dp)               :: dxsp(3,nsp), grphi(3), phip, r2o,
     &                          r2sp(nsp)

!--------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '      PRE PhiOnMesh' )
#endif

      call timer( 'PHION', 1 )
C     On first call set the logical DirectPhi
      if (firsttime) then
        DirectPhi = fdf_get('DirectPhi', .false. )
        firsttime = .false.
      else
        ! reset the data structures in module meshphi
        call resetMeshPhi( )
      endif

C     Allocate memory related to nmpl
      nullify( endpht, phi, lstpht, listp2 )
      call re_alloc( endpht, 0, nmpl, 'endpht', 'PhiOnMesh' )

C     Initialise pointer array
      endpht(0) = 0
      do ip = 1,nmpl
        endpht(ip) = endpht(ip-1) + numphi(ip)
      enddo

C     Allocate phi if this is not a direct calculation
      nlist = endpht(nmpl)
      if (DirectPhi) then
        nphi = 1
      else
        nphi = nlist
        if (IONode) then
           write(6,"(a,i20)")
     $          "PhiOnMesh: Number of (b)points on node 0 = ", nmpl
           write(6,"(a,i20)")
     $          "PhiOnMesh: nlist on node 0 = ", nlist
        endif
      endif

C     Add an extra margin of error to nlist to minimise reallocations
      nliste = 1.01*nlist

C     Adjust dimensions of phi if necessary
      if (associated(phi)) then
        nmax = size(phi,2)
      else
        nmax = 0
      end if
      if (nphi.gt.nmax) then
        if (DirectPhi) then
          call re_alloc( phi, 1, nsp, 1, nphi, 'phi', 'PhiOnMesh' )
       else
          call re_alloc( phi, 1, nsp, 1, nliste, 'phi', 'PhiOnMesh' )
       endif
      endif

C     Adjust dimensions of list arrays if necessary
      if (associated(lstpht) .and. associated(listp2)) then
        nmax = min( size(lstpht), size(listp2) )
      else
        nmax = 0
      end if
      if (nlist.gt.nmax) then
        call re_alloc( lstpht, 1, nliste, 'lstpht', 'PhiOnMesh' )
        call re_alloc( listp2, 1, nliste, 'listp2', 'PhiOnMesh' )
      endif

C     Find indexes and values of atomic orbitals at mesh points
      numphi = 0
      do io = 1,norb
        ia = iaorb(io)
        if (ipa(ia)==0) cycle ! Atomic sphere does not intersect my Box
        iphi = iphorb(io)
        is = isa(ia)
        r2o = rcut(is,iphi)**2
        do iop = 1, mop
          ip0 = indexp( ipa(ia) + idop(iop) )
          if (ip0 .gt. 0) then
            within = .false.
            do isp = 1,nsp
              dxsp(1:3,isp) = xdop(1:3,iop)+xdsp(1:3,isp)-dxa(1:3,ia)
              r2sp(isp) = sum( dxsp(1:3,isp)**2 )
              if (r2sp(isp) .lt. r2o) within = .true.
            enddo
            if (within) then
              numphi(ip0) = numphi(ip0) + 1
              n = endpht(ip0-1) + numphi(ip0)
              lstpht(n) = io
              listp2(n) = iop
              if (.not.DirectPhi) then
                do isp = 1,nsp
                  if (r2sp(isp) .lt. r2o) then
                    call phiatm( is, iphi, dxsp(1,isp),
     &                           phip, grphi )
                    phi(isp,n) = phip
                  else
                    phi(isp,n) = 0.0_grid_p
                  endif
                enddo
              endif
            endif
          endif
        enddo
      enddo

      call timer( 'PHION', 2 )
#ifdef DEBUG
      call write_debug( '      POS PhiOnMesh' )
#endif
!------------------------------------------------------------------------ END
      end subroutine PhiOnMesh
