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
      subroutine efield( cell, na, isa, xa, mesh, nsm, dipole, 
     .                   v, field )
c **********************************************************************
c Adds the potential created by an external electric field, whose value
c is readed from the FDF block ExternalElectricField.
c Written by J.M.Soler. Feb. 1998.
c Modified to operate only on the sub-matrix of the potential stored
c locally. J.D.Gale March 1999.
c Modified to calculate the electric field needed to compensate
c the system dipole in slab calculations.   P. Ordejon, July 2009
c ********* Input ******************************************************
c real*8  cell(3,3) : Unit cell vectors
c integer na        : Number of atoms
c integer isa(na)   : Atomic species indexes
c real*8  xa(3,na)  : Atomic positions (cartesian coordinates)
c integer mesh(3)   : Number of mesh divisions in each cell direction
c integer nsm       : Number of sub-mesh points along each axis
c integer dipole(3) : Dipole moment of system
c ********* Input and output *******************************************
c real    v(*)      : Electron potential, to which that created by the
c                     electric field is added. Notice single precision.
c ********* Output *****************************************************
c real*8  field(3)  : Electric field
c ********* Units ******************************************************
c Distances in Bohr radiae
c Energies in Rydbergs
c Electric field in Ry/Bohr
c Dipoles in electrons*Bohr
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
c The electric field required can be added in two possible ways:
c 1) FDF block 'ExternalElectricField'. This routine reads the
c    electric field input by the user, and adds the corresponding
c    potential to the grid.
c    Sample FDF electric field specification:
c       %block ExternalElectricField
c           0.000  0.000  3.000  V/Ang
c       %endblock ExternalElectricField
c 2) FDF Option 'SlabDipoleCorrection'. This options calculates
c    the electric field required to compensate the dipole of the
c    system. The potential added to the grid corresponds to that
c    of a dipole layer at the middle of the vacuum layer. For slabs,
c    this exactly compensates the electric field at the vacuum, created
c    by the dipole moment of the system, thus allowing to deal with
c    asymmetric slabs (and compute properties such as the work funcion
c    of each of the surfaces).
c    Sample FDF input:   
c       SlabDipoleCorrection   .true.
c    The dipole moment of the system is input as an argument of the
c    routine.
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
      use sys,          only: die
#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

      integer           na, isa(na), mesh(3), nsm
      real(grid_p), intent(inout)  ::  v(*)
      real(dp)          cell(3,3), ddot, field(3), xa(3,na), dipole(3),
     .                  volcel
      external          cross, ddot, reclat, shaper, volcel

C Internal parameters
C tol : tolerance for bulk components of the electric field
      real(dp), save :: tol = 1.0d-12

C Internal variables
      logical,  save :: found1   
      logical           found2, found, isfield, orthog
      logical,  save :: frstme = .true., lwrite
      character         eunits*10, shape*8, line*130, names*20
      integer           i0(3), i1, i2, i3, ia, imesh, int(1),
     $                  is, iu, ix,
     .                  j1, j2, j3, last, lc(0:1), 
     .                  nbcell, ni, nn, nr, nv, meshl(3),
     .                  ProcessorZ, BlockSizeY, BlockSizeZ, Yoffset,
     .                  Zoffset, Py, Pz, i30, i20, meshnsm(3),
     .                  NRemY, NRemZ
#ifdef MPI
      integer           MPIerror, npl
#endif
      real(dp)          b1xb2(3), bcell(3,3), cfactor, dplane(3),
     .                  e(3), e0(3), esave(3), eb1, eb2, eb3,
     .                  f(3), pi, rc, rcell(3,3), v0,
     .                  xfrac, xmax(3), xmean, xmin(3)
      save              e, f, isfield, i0, v0, found2, esave


      pi =4.*atan(1.0d0)
C Find and store the electric field only the first time
      lwrite = .false.
      if (frstme) then
        lwrite = .true.
        frstme = .false.
        isfield = .false.
        e = 0.0_dp

        if (ionode) then
C Read the electric field block from the fdf input file
          found1 = fdf_block('ExternalElectricField',iu)
          if (found1) then
            read(iu,'(a)') line
            last = index(line,'#') - 1
            if (last .le. 0) last = len(line)
            call parse( line(1:last), nn, lc, names, nv, e,
     .                ni, int, nr, e0 )
            eunits = names(lc(0)+1:lc(1))
            cfactor = fdf_convfac(eunits,'Ry/Bohr/e')
            do ix = 1,3
              e(ix) = e(ix) * cfactor
              e0(ix) = e(ix)
              esave(ix) = e(ix)
            enddo
          endif
        endif
      else
        do ix = 1,3
          e(ix) = esave(ix)
        enddo
      endif

      if (ionode) then
        found2=fdf_boolean('SlabDipoleCorrection',.false.)
        if (found2) then
          if (lwrite) then
            write(6,'(/,(a))')
     .        'efield: SlabDipoleCorrection = .true. in input file',
     .        'efield: A dipole layer will be introduced in the vacuum',
     .        'efield: region to compensate the system dipole'
          endif
          do ix = 1,3
            e(ix) = e(ix) - 4.0*pi*dipole(ix)/volcel(cell)*2.0
            e0(ix) = e(ix)
            if (e(ix) .ne. 0.0_dp) isfield = .true.
          enddo
        endif
        found = (found1 .or. found2)
      endif

#ifdef MPI
      call MPI_Bcast(e,3,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
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
     .          ddot(3,bcell(1,1),1,bcell(1,1),1)
          eb2 = ddot(3,e,1,bcell(1,2),1)/
     .          ddot(3,bcell(1,2),1,bcell(1,2),1)
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
            if (lwrite) then
              write(6,'(/,a,3f12.6,a)')
     .          'efield: External electric field =', e, ' Ry/Bohr/e'
              if (cml_p) call cmlAddProperty(xf=mainXML, value=e, 
     .           dictref='siesta:elfield', 
     .           units='siestaUnits:Ry_Bohr_e')
            endif
          else
            write(6,'(a,(/,a,3f12.6))')
     .        'efield: ERROR: Non zero bulk electric field.',
     .        'efield: Input field (Ry/Bohr/e) =', e0,
     .        'efield: Orthogonalized field    =', e
          endif
        endif
      endif

C Find the origin of a shifted cell, with the system centered in it
C This is done at every call, because of possible atomic movements
      if (isfield) then

C Find reciprocal unit cell and distance between lattice planes
        call reclat( cell, rcell, 0 )
        do ix = 1,3
          dplane(ix) = 1.0_dp /
     $                   sqrt(ddot(3,rcell(1,ix),1,rcell(1,ix),1))
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
          i0(ix) = nint( (xmean-0.5_dp) * mesh(ix) )
        enddo

C Find the electric field in mesh coordinates, so that
C v = e*x = f*index
        do ix = 1,3
          f(ix) = ddot(3,e,1,cell(1,ix),1) / max( mesh(ix), 1 )
        enddo

C Find the potential at the origin of the shifted cell, so that
C the potential is zero at the center of the cell
        v0 = (- 0.5_dp) * (f(1)*mesh(1) + f(2)*mesh(2) + f(3)*mesh(3))
      endif        

      if (isfield) then

C Find local number of mesh points
        meshnsm(1) = mesh(1)/nsm
        meshnsm(2) = mesh(2)/nsm
        meshnsm(3) = mesh(3)/nsm
#ifdef MPI
        call HowManyMeshPerNode(meshnsm,Node,Nodes,npl,meshl)
        meshl(1) = meshl(1)*nsm
        meshl(2) = meshl(2)*nsm
        meshl(3) = meshl(3)*nsm
#else
        meshl(1) = mesh(1)
        meshl(2) = mesh(2)
        meshl(3) = mesh(3)
#endif

C Check that ProcessorY is a factor of the number of processors
      if (mod(Nodes,ProcessorY).gt.0)
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')
        ProcessorZ = Nodes/ProcessorY

C Calculate blocking sizes
        BlockSizeY = (meshnsm(2)/ProcessorY)*nsm
        NRemY = (mesh(2) - ProcessorY*BlockSizeY)/nsm
        BlockSizeZ = (meshnsm(3)/ProcessorZ)*nsm
        NRemZ = (mesh(3) - ProcessorZ*BlockSizeZ)/nsm

C Calculate coordinates of current node in processor grid
        Py = (Node/ProcessorZ)+1
        Pz = Node - (Py - 1)*ProcessorZ + 1

C Calculate starting point for grid
        Yoffset = (Py-1)*BlockSizeY + nsm*min(Py-1,NRemY)
        Zoffset = (Pz-1)*BlockSizeZ + nsm*min(Pz-1,NRemZ)

C Add the electric field potential to the input potential
        imesh = 0
        i30 = Zoffset - 1
        do i3 = 0,meshl(3)-1
          i30 = i30 + 1
          i20 = Yoffset - 1
          do i2 = 0,meshl(2)-1
            i20 = i20 + 1
            do i1 = 0,meshl(1)-1
              imesh = imesh + 1
              j1 = mod( i1-i0(1)+10*mesh(1), mesh(1) )
              j2 = mod( i20-i0(2)+10*mesh(2), mesh(2) )
              j3 = mod( i30-i0(3)+10*mesh(3), mesh(3) )
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
