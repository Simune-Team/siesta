c $Id: efield.f,v 1.4 1999/01/31 11:14:51 emilio Exp $

      subroutine efield( cell, na, isa, xa, mesh, v, field )

c **********************************************************************
c Adds the potential created by an external electric field, whose value
c is readed from the FDF block ExternalElectricField.
c Written by J.M.Soler. Feb. 1998.
c ********* Input ******************************************************
c real*8  cell(3,3) : Unit cell vectors
c integer na        : Number of atoms
c integer isa(na)   : Atomic species indexes
c real*8  xa(3,na)  : Atomic positions (cartesian coordinates)
c integer mesh(3)   : Number of mesh divisions in each cell direction
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

      implicit          none
      integer           na, isa(na), mesh(3)
      real              v(*)
      double precision  cell(3,3), dot, field(3), rcut, xa(3,na)
      external          cross, dot, parse, rcut, reclat, shaper

c Internal parameters
c tol : tolerance for bulk components of the electric field
      double precision tol
      parameter ( tol = 1.d-12 )

c Internal variables
      logical           found, frstme, isfield, orthog
      character         eunits*10, shape*8, line*130, names*20
      integer           i0(3), i1, i2, i3, ia, imesh, int, is, iu, ix,
     .                  j1, j2, j3, last, lc(0:1),
     .                  nbcell, ni, nn, nr, nv
      double precision  b1xb2(3), bcell(3,3), cfactor, dplane(3),
     .                  e(3), e0(3), eb1, eb2, eb3,
     .                  f(3), rc, rcell(3,3), v0,
     .                  xfrac, xmax(3), xmean, xmin(3)
      save              e, f, frstme, isfield, i0, v0

c FDF-related declarations
      double precision fdf_convfac
      external fdf_convfac
      include 'fdf/fdfdefs.h'

      data frstme /.true./
*     data eunits /'          '/

c Find and store the electric field only the first time
      if (frstme) then
        frstme = .false.
        isfield = .false.
        do ix = 1,3
          e(ix) = 0.d0
        enddo

c       Read the electric field block from the fdf input file
        found = fdf_block('ExternalElectricField',iu)
        if (found) then
          read(iu,'(a)') line
          last = index(line,'#') - 1
          if (last .le. 0) last = len(line)
          call parse( line(1:last), nn, lc, names, nv, e,
     .                ni, int, nr, e0 )
          eunits = names(lc(0)+1:lc(1))
          cfactor = fdf_convfac(eunits,'Ry/Bohr/e')
          do ix = 1,3
            if (e(ix) .ne. 0.d0) isfield = .true.
            e(ix) = e(ix) * cfactor
            e0(ix) = e(ix)
          enddo
        endif

c       Check that the field is orthogonal to the bulk directions
        if (isfield) then
          call shaper( cell, na, isa, xa, shape, nbcell, bcell )
          orthog = .true.
          if (nbcell .eq. 1) then
            eb1 = dot(e,bcell,3) / dot(bcell,bcell,3)
            if (abs(eb1) .gt. tol) then
              orthog = .false.
              do ix = 1,3
                e(ix) = e(ix) - eb1 * bcell(ix,1)
              enddo
            endif
          elseif (nbcell .eq. 2) then
            eb1 = dot(e,bcell(1,1),3) / dot(bcell(1,1),bcell(1,1),3)
            eb2 = dot(e,bcell(1,2),3) / dot(bcell(1,2),bcell(1,2),3)
            if (abs(eb1).gt.tol .or. abs(eb2).gt.tol) then
              orthog = .false.
              call cross( bcell(1,1), bcell(1,2), b1xb2 )
              eb3 = dot(e,b1xb2,3) / dot(b1xb2,b1xb2,3)
              do ix = 1,3
                e(ix) = eb3 * b1xb2(ix)
              enddo
            endif
          elseif (nbcell .eq. 3) then
            orthog = .false.
            do ix = 1,3
              e(ix) = 0.d0
            enddo
          endif
          if (orthog) then
            write(6,'(/,a,3f12.6,a))')
     .        'efield: Electric field =', e, ' Ry/Bohr/e'
          else
            write(6,'(a,(/,a,3f12.6))')
     .        'efield: ERROR: Non zero bulk electric field.',
     .        'efield: Input field (Ry/Bohr/e) =', e0,
     .        'efield: Orthogonalized field    =', e
          endif
        endif
      endif

c Find the origin of a shited cell, with the system centered in it
c This is done at every call, because of possible atomic movements
      if (isfield) then

c       Find reciprocal unit cell and distance between lattice planes
        call reclat( cell, rcell, 0 )
        do ix = 1,3
          dplane(ix) = sqrt( dot(rcell(1,ix),rcell(1,ix),3) )
        enddo

c       Find the geometric center of the system
        do ix = 1,3
          xmin(ix) =  1.d30
          xmax(ix) = -1.d30
        enddo
        do ia = 1,na
          is = isa(ia)
          rc = rcut(is,0)
          do ix = 1,3
            xfrac = dot( xa(1,ia), rcell(1,ix), 3 )
            xmin(ix) = min( xmin(ix), xfrac-rc/dplane(ix) )
            xmax(ix) = max( xmax(ix), xfrac+rc/dplane(ix) )
          enddo
        enddo

c       Find the mesh index of the origin of the shifted cell
        do ix = 1,3
          xmean = (xmin(ix) + xmax(ix)) / 2
          i0(ix) = nint( (xmean-0.5d0) * mesh(ix) )
        enddo

c       Find the electric field in mesh coordinates, so that
c       v = e*x = f*index
        do ix = 1,3
          f(ix) = dot( e, cell(1,ix), 3 ) / max( mesh(ix), 1 )
        enddo

c       Find the potential at the origin of the shifted cell, so that
c       the potential is zero at the center of the cell
        v0 = (- 0.5d0) * (f(1)*mesh(1) + f(2)*mesh(2) + f(3)*mesh(3))
      endif        

c Add the electric field potential to the input potential
      if (isfield) then
        imesh = 0
        do i3 = 0,mesh(3)-1
        do i2 = 0,mesh(2)-1
        do i1 = 0,mesh(1)-1
          imesh = imesh + 1
          j1 = mod( i1-i0(1)+10*mesh(1), mesh(1) )
          j2 = mod( i2-i0(2)+10*mesh(2), mesh(2) )
          j3 = mod( i3-i0(3)+10*mesh(3), mesh(3) )
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


