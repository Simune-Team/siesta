      subroutine coor(na,cell)

c Reads atomic coordinates and format in which they are given, to be
c transformed into Bohr cartesian for internal handling. 
c It also shifts them all according to AtomicCoordinatesOrigin.
c Written by E. Artacho, December 1997, on the original piece of the
c redata subroutine written by P. Ordejon in December 1996.
c Modified by J.M.Soler. August 1998.

!! Modified by Alberto Garcia, May 16, 2000

c integer na  : Number of atoms in unit cell (or supercell, if defined)

      use precision
      use atomlist, only: xa, isa
      use fdf2
      use alloc

      implicit          none

      integer, intent(out)          :: na
      double precision, intent(out) :: cell(3,3)

c Internal variables and arrays
      character         acf*22, acf_default*22
      logical           isdiag, leqi
      integer           i, ia, ic, iscale, iua, iunit, ix,
     .                  mscell(3,3), ncells, nsc(3), nua
      double precision  alat, dcell(3,3), dscell(3,3), origin(3),
     .                  scell(3,3), ucell(3,3), volcel, xac(3)
      external          digcel, redcel, superx, volcel

      data origin /3*0.d0/

C Read unit cell and supercell
      call redcel( alat, ucell, scell, mscell )
!
!     Set cell argument to "supercell" for use by the program
!
      cell = scell

C Format of atomic coordinates

      acf_default = 'NotScaledCartesianBohr'
      acf = fdf_get('AtomicCoordinatesFormat',acf_default)
      if (leqi(acf,'NotScaledCartesianBohr') .or.
     .    leqi(acf,'Bohr') ) then
        iscale = 0
        write(6,'(a,a)')
     .   'coor: Atomic-coordinates input format  = ',
     .   '    Cartesian coordinates'
        write(6,'(a,a)')
     .   'coor:                                    ',
     .   '    (in Bohr units)'
      else if (leqi(acf,'NotScaledCartesianAng') .or.
     .         leqi(acf,'Ang') ) then
        iscale = 1
        write(6,'(a,a)')
     .   'coor: Atomic-coordinates input format  = ',
     .   '    Cartesian coordinates'
        write(6,'(a,a)')
     .   'coor:                                    ',
     .   '    (in Angstroms)'
      else if (leqi(acf,'ScaledCartesian')) then
        if (alat.eq.0.d0) then
           write(6,"(/,2a)") 'coor: ERROR: Explicit lattice ',
     .       'constant is needed for ScaledCartesian format'
           call die
        endif
        iscale = 2
        write(6,'(a,a)')
     .   'coor: Atomic-coordinates input format  = ',
     .   '    Cartesian coordinates'
        write(6,'(a,a)')
     .   'coor:                                    ',
     .   '    (in units of alat)'
      else if (leqi(acf,'ScaledByLatticeVectors') .or. 
     .         leqi(acf,'Fractional') ) then
        if (alat.eq.0.d0) then
           write(6,"(/,2a)") 'coor: ERROR: Explicit lattice ',
     .       'constant is needed for Fractional format'
           call die
        endif
        iscale = 3
        write(6,'(a,a)')
     .   'coor: Atomic-coordinates input format  = ',
     .   '    Fractional'
      else
        write(6,"(/,'coor: ',72(1h*))")
        write(6,"('coor:                  INPUT ERROR')")
        write(6,'(a)') 'coor: '
        write(6,'(2a)') 'coor: You must use one of the following',
     .                            ' coordinate scaling options:'
        write(6,'(a)') 'coor:     - NotScaledCartesianBohr (or Bohr)'
        write(6,'(a)') 'coor:     - NotScaledCartesianAng (or Ang) '
        write(6,'(a)') 'coor:     - ScaledCartesian                '
        write(6,'(2a)') 'coor:     - ScaledByLatticeVectors ',
     .                                               '(or Fractional)'
        write(6,"('coor: ',72(1h*))")

        call die

      endif


c Read atomic coordinates and species

      na = fdf_get('NumberOfAtoms',0)
      if (na.eq.0) call die("Must specify number of atoms!")
!
!     Check if we need more space to accomodate supercell
!     (still a "real" supercell, specified for convenience!)
!
      if (volcel(ucell) .lt. 1.d-8) then
        ncells = 1
      else
        ncells = nint( volcel(scell) / volcel(ucell) )
      endif

      nua = na
      na  = na * ncells

      nullify(isa,xa)
      call realloc(isa,1,na,name='isa',routine='coor')
      call realloc(xa,1,3,1,na,name='xa',routine='coor')

      if ( fdf_block('AtomicCoordinatesAndAtomicSpecies',iunit) )
     .     then
         do ia = 1, nua
            read(iunit,*) (xa(i,ia), i=1,3), isa(ia)
         enddo
      else
         call die("coo: You must specify the atomic coordinates")
      endif

      if ( fdf_block('AtomicCoordinatesOrigin',iunit) ) then
         read(iunit,*) (origin(i),i=1,3)
         do ia = 1,nua
            do i = 1,3
               xa(i,ia) = xa(i,ia) + origin(i)
            enddo
         enddo
      endif

c       Scale atomic coordinates
c       Coord. option = 0 => Do nothing
c       Coord. option = 1 => Multiply by 1./0.529177 (Ang --> Bohr)
c       Coord. option = 2 => Multiply by lattice constant
c       Coord. option = 3 => Multiply by lattice vectors

      if (iscale .eq. 1) then
         xa = xa / 0.529177d0 
      elseif (iscale .eq. 2) then
         xa = xa * alat
      elseif (iscale .eq. 3) then
         do ia = 1,nua
            do ix = 1,3
               xac(ix) = xa(ix,ia)
            enddo
            do ix = 1,3
               xa(ix,ia) = ucell(ix,1) * xac(1) +
     .              ucell(ix,2) * xac(2) +
     .              ucell(ix,3) * xac(3)
            enddo
         enddo
      endif

c Expand the coordinates to whole supercell

      if (ncells.gt.1) then

C         Find equivalent diagonal combination of cell/supercell
          call digcel( ucell, mscell, dcell, dscell, nsc, isdiag )

C         Expand coordinates
          call superx( dcell, nsc, nua, na, xa, dscell )

C         Expand index isa
          ia = 0
          do ic = 1,ncells
            do iua = 1,nua
              ia = ia + 1
              isa(ia) = isa(iua)
            enddo
          enddo

        endif

      end
