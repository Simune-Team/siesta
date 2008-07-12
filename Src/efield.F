! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine efield( cell, na, isa, xa, meshG, nsm, v, field )
c **********************************************************************
c Adds the potential created by an external electric field, whose value
c is readed from the FDF block ExternalElectricField.
c Written by J.M.Soler. Feb. 1998.
c Modified to operate only on the sub-matrix of the potential stored
c locally. J.D.Gale March 1999.
c ********* Input ******************************************************
c real*8  cell(3,3) : Unit cell vectors
c integer na        : Number of atoms
c integer isa(na)   : Atomic species indexes
c real*8  xa(3,na)  : Atomic positions (cartesian coordinates)
c integer meshG(3)  : Number of mesh divisions in each cell direction
c integer nsm       : Number of sub-mesh points along each axis
c ********* Input and output *******************************************
c real    v(*)      : Electron potential, to which that created by the
c                     electric field is added. Notice single precision.
c ********* Output *****************************************************
c real*8  field(3)  : Electric field
c ********* Units ******************************************************
c Distances in Bohr radiae
c Energies in Rydbergs
c Electric field in Ry/Bohr
c ********* Behaviour **************************************************
c The sign of the potential is that for electrons (v=+E*x), i.e. 
c  opposite to that of the conventional electrostatic potential.
c Notice that the potential is not initialized.
c Bulk electric fields are not allowed. If the specified electric field
c  is not orthogonal to all bulk directions, it is orthogonalized, and
c  a warning message is printed.
c The electric field produces a discontinuity of the potential in the
c  periodic cell, which is automatically placed in the middle of the
c  vacuum region.
c The output electric field is obtained even if mesh=0 (so that no
c  potential may be calculated)
c ********* Usage ******************************************************
c Sample FDF electric field specification:
c    %block ExternalElectricField
c        0.000  0.000  3.000  V/Ang
c    %endblock ExternalElectricField
c **********************************************************************

C
C  Modules
C
      use precision, only: dp, grid_p

      use parallel,     only : IOnode, Node, Nodes, ProcessorY
      use parallelsubs, only : HowManyMeshPerNode
      use atmfuncs,     only : rcut
      use fdf
      use parsing
      use siesta_cml
      use sys,          only : die
      use mesh,         only : meshLim
#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

      integer           na, isa(na), meshG(3), nsm
      real(grid_p), intent(inout)  ::  v(*)
      real(dp)          cell(3,3), ddot, field(3), xa(3,na)
      external          cross, ddot, reclat, shaper

C Internal parameters
C tol : tolerance for bulk components of the electric field
      real(dp), save :: tol = 1.0d-12

C Internal variables
      logical           found, isfield, orthog
      logical,  save :: frstme = .true.
      character         eunits*10, shape*8, line*130, names*20
      integer           i0(3), i1, i2, i3, ia, imesh, int(1),
     &                  is, iu, ix,
     &                  j1, j2, j3, last, lc(0:1), 
     &                  nbcell, ni, nn, nr, nv, meshl(3),
     &                  Xoffset, Yoffset, Zoffset, i30, i20, i10
#ifdef MPI
      integer           MPIerror, npl
#endif
      real(dp)          b1xb2(3), bcell(3,3), cfactor, dplane(3),
     &                  e(3), e0(3), eb1, eb2, eb3,
     &                  f(3), rc, rcell(3,3), v0,
     &                  xfrac, xmax(3), xmean, xmin(3)
      save              e, f, isfield, i0, v0

C Find and store the electric field only the first time
      if (frstme) then
        frstme = .false.
        isfield = .false.
        e = 0.0_dp

        if (ionode) then
C Read the electric field block from the fdf input file
          found = fdf_block('ExternalElectricField',iu)
          if (found) then
            read(iu,'(a)') line
            last = index(line,'#') - 1
            if (last .le. 0) last = len(line)
            call parse( line(1:last), nn, lc, names, nv, e,
     &                ni, int, nr, e0 )
            eunits = names(lc(0)+1:lc(1))
            cfactor = fdf_convfac(eunits,'Ry/Bohr/e')
            do ix = 1,3
              if (e(ix) .ne. 0.0_dp) isfield = .true.
              e(ix) = e(ix) * cfactor
              e0(ix) = e(ix)
            enddo
          endif
        endif
#ifdef MPI
        call MPI_Bcast(e,3,MPI_double_precision,0,MPI_Comm_World,
     &    MPIerror)
        call MPI_Bcast(isfield,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

C Check that the field is orthogonal to the bulk directions
        if (isfield) then
          call shaper( cell, na, isa, xa, shape, nbcell, bcell )
          orthog = .true.
          if (nbcell .eq. 1) then
            eb1 = ddot(3,e,1,bcell,1) / ddot(3,bcell,1,bcell,1)
            if (abs(eb1) .gt. tol) then
              orthog = .false.
              do ix = 1,3
                e(ix) = e(ix) - eb1 * bcell(ix,1)
              enddo
            endif
          elseif (nbcell .eq. 2) then
            eb1 = ddot(3,e,1,bcell(1,1),1)/
     &            ddot(3,bcell(1,1),1,bcell(1,1),1)
            eb2 = ddot(3,e,1,bcell(1,2),1)/
     &            ddot(3,bcell(1,2),1,bcell(1,2),1)
            if (abs(eb1).gt.tol .or. abs(eb2).gt.tol) then
              orthog = .false.
              call cross( bcell(1,1), bcell(1,2), b1xb2 )
              eb3 = ddot(3,e,1,b1xb2,1)/ddot(3,b1xb2,1,b1xb2,1)
              do ix = 1,3
                e(ix) = eb3 * b1xb2(ix)
              enddo
            endif
          elseif (nbcell .eq. 3) then
            orthog = .false.
            do ix = 1,3
              e(ix) = 0.0_dp
            enddo
          endif
          if (ionode) then
            if (orthog) then
              write(6,'(/,a,3f12.6,a)')
     &          'efield: Electric field =', e, ' Ry/Bohr/e'
              if (cml_p) call cmlAddProperty(xf=mainXML, property=e, 
     &             dictref='siesta:elfield', 
     &             units='siestaUnits:Ry_Bohr_e')
            else
              write(6,'(a,(/,a,3f12.6))')
     &          'efield: ERROR: Non zero bulk electric field.',
     &          'efield: Input field (Ry/Bohr/e) =', e0,
     &          'efield: Orthogonalized field    =', e
            endif
          endif
        endif
      endif

C Find the origin of a shited cell, with the system centered in it
C This is done at every call, because of possible atomic movements
      if (isfield) then

C Find reciprocal unit cell and distance between lattice planes
        call reclat( cell, rcell, 0 )
        do ix = 1,3
          dplane(ix) = sqrt(ddot(3,rcell(1,ix),1,rcell(1,ix),1))
        enddo

C Find the geometric center of the system
        do ix = 1,3
          xmin(ix) =  1.0e30_dp
          xmax(ix) = -1.0e30_dp
        enddo
        do ia = 1,na
          is = isa(ia)
          rc = rcut(is,0)
          do ix = 1,3
            xfrac = ddot(3,xa(1,ia),1,rcell(1,ix),1)
            xmin(ix) = min( xmin(ix), xfrac-rc/dplane(ix) )
            xmax(ix) = max( xmax(ix), xfrac+rc/dplane(ix) )
          enddo
        enddo

C Find the mesh index of the origin of the shifted cell
        do ix = 1,3
          xmean = (xmin(ix) + xmax(ix)) / 2
          i0(ix) = nint( (xmean-0.5_dp) * meshG(ix) )
        enddo

C Find the electric field in mesh coordinates, so that
C v = e*x = f*index
        do ix = 1,3
          f(ix) = ddot(3,e,1,cell(1,ix),1) / max( meshG(ix), 1 )
        enddo

C Find the potential at the origin of the shifted cell, so that
C the potential is zero at the center of the cell
        v0 = (-0.5_dp)*(f(1)*meshG(1) + f(2)*meshG(2) + f(3)*meshG(3))
      endif        

      if (isfield) then
        meshl(1) = (meshLim(2,1) - meshLim(1,1)+1)*nsm
        meshl(2) = (meshLim(2,2) - meshLim(1,2)+1)*nsm
        meshl(3) = (meshLim(2,3) - meshLim(1,3)+1)*nsm

C Calculate starting point for grid
        Xoffset = (meshLim(1,1)-1)*nsm
        Yoffset = (meshLim(1,2)-1)*nsm
        Zoffset = (meshLim(1,3)-1)*nsm

C Add the electric field potential to the input potential
        imesh = 0
        i30 = Zoffset - 1
        do i3 = 0,meshl(3)-1
          i30 = i30 + 1
          i20 = Yoffset - 1
          do i2 = 0,meshl(2)-1
            i20 = i20 + 1
            i10 = Yoffset - 1
            do i1 = 0,meshl(1)-1
              i10 = i10 + 1
              imesh = imesh + 1
              j1 = mod( i10-i0(1)+10*meshG(1), meshG(1) )
              j2 = mod( i20-i0(2)+10*meshG(2), meshG(2) )
              j3 = mod( i30-i0(3)+10*meshG(3), meshG(3) )
              v(imesh) = v(imesh) + v0 + f(1)*j1 + f(2)*j2 + f(3)*j3
            enddo
          enddo
        enddo
      endif

C Copy the electric field to the output array
      do ix = 1,3
        field(ix) = e(ix)
      enddo
      end
