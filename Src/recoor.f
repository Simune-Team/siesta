c $Id: recoor.f,v 1.6 1999/04/13 17:03:25 emilio Exp $

      subroutine recoor( maxa, na, isa, xa ) 

c *******************************************************************
c Reads atomic coordinates and format in which they are given, to be
c transformed into Bohr cartesian for internal handling. 
c It also shifts them all according to AtomicCoordinatesOrigin.
c Written by E. Artacho, December 1997, on the original piece of the
c redata subroutine written by P. Ordejon in December 1996.
c Modified by J.M.Soler. August 1998.
c ********* INPUT ***************************************************
c integer maxa : Dimension of arrays xa and isa
c ********* OUTPUT **************************************************
c integer na  : Number of atoms in unit cell (or supercell, if defined)
c integer isa(maxa)  : Atomic species of different atoms
c real*8  xa(3,maxa) : Atomic catesian coordinates in Bohr
c *******************************************************************

      implicit          none
      integer           maxa, na
      integer           isa(maxa)
      double precision  xa(3,maxa)

c Internal variables and arrays
      character         acf*22, acf_default*22
      logical           isdiag, leqi
      integer           i, ia, ic, iscale, iua, iunit, ix,
     .                  mscell(3,3), ncells, nsc(3), nua
      double precision  alat, dcell(3,3), dscell(3,3), origin(3),
     .                  scell(3,3), ucell(3,3), volcel, xac(3)
      external          digcel, redcel, superx, volcel

c Enable FDF input/output
      include 'fdf/fdfdefs.h'

      data origin /3*0.d0/

C Read unit cell and supercell
      call redcel( alat, ucell, scell, mscell )

C Format of atomic coordinates
      acf_default = 'NotScaledCartesianBohr'
      acf = fdf_string('AtomicCoordinatesFormat',acf_default)
      if (leqi(acf,'NotScaledCartesianBohr') .or.
     .    leqi(acf,'Bohr') ) then
        iscale = 0
        write(6,'(a,a)')
     .   'recoor: Atomic-coordinates input format  = ',
     .   '    Cartesian coordinates'
        write(6,'(a,a)')
     .   'recoor:                                    ',
     .   '    (in Bohr units)'
      else if (leqi(acf,'NotScaledCartesianAng') .or.
     .         leqi(acf,'Ang') ) then
        iscale = 1
        write(6,'(a,a)')
     .   'recoor: Atomic-coordinates input format  = ',
     .   '    Cartesian coordinates'
        write(6,'(a,a)')
     .   'recoor:                                    ',
     .   '    (in Angstroms)'
      else if (leqi(acf,'ScaledCartesian')) then
        if (alat.eq.0.d0) then
           write(6,"(/,2a)") 'recoor: ERROR: Explicit lattice ',
     .       'constant is needed for ScaledCartesian format'
           stop 'recoor: ERROR: Explicit lattice constant needed'
        endif
        iscale = 2
        write(6,'(a,a)')
     .   'recoor: Atomic-coordinates input format  = ',
     .   '    Cartesian coordinates'
        write(6,'(a,a)')
     .   'recoor:                                    ',
     .   '    (in units of alat)'
      else if (leqi(acf,'ScaledByLatticeVectors') .or. 
     .         leqi(acf,'Fractional') ) then
        if (alat.eq.0.d0) then
           write(6,"(/,2a)") 'recoor: ERROR: Explicit lattice ',
     .       'constant is needed for Fractional format'
           stop 'recoor: ERROR: Explicit lattice constant needed'
        endif
        iscale = 3
        write(6,'(a,a)')
     .   'recoor: Atomic-coordinates input format  = ',
     .   '    Fractional'
      else
        write(6,"(/,'recoor: ',72(1h*))")
        write(6,"('recoor:                  INPUT ERROR')")
        write(6,'(a)') 'recoor: '
        write(6,'(2a)') 'recoor: You must use one of the following',
     .                            ' coordinate scaling options:'
        write(6,'(a)') 'recoor:     - NotScaledCartesianBohr (or Bohr)'
        write(6,'(a)') 'recoor:     - NotScaledCartesianAng (or Ang) '
        write(6,'(a)') 'recoor:     - ScaledCartesian                '
        write(6,'(2a)') 'recoor:     - ScaledByLatticeVectors ',
     .                                               '(or Fractional)'
        write(6,"('recoor: ',72(1h*))")
        stop 'recoor: ERROR: Wrong atomic-coordinate input format'
      endif


c Read atomic coordinates and species

      na = fdf_integer('NumberOfAtoms',0)
      if (na.gt.0 .and. na.le.maxa) then

        if ( fdf_block('AtomicCoordinatesAndAtomicSpecies',iunit) )
     .    then
          do ia = 1,na
            read(iunit,*) (xa(i,ia), i=1,3), isa(ia)
          enddo
        else
          write(6,"(/,'recoor: ',72(1h*))")
          write(6,"('recoor:                  INPUT ERROR')")
          write(6,'(a)')
     .    'recoor:   You must specify the atomic coordinates'
          write(6,"('recoor: ',72(1h*))")
          stop 'recoor: ERROR: Atomic coordinates missing'
        endif

        if ( fdf_block('AtomicCoordinatesOrigin',iunit) ) then
          read(iunit,*) (origin(i),i=1,3)
          do ia = 1,na
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
          do ia = 1,na
            do ix = 1,3
              xa(ix,ia) = 1.d0 / 0.529177d0 * xa(ix,ia)
            enddo
          enddo
        elseif (iscale .eq. 2) then
          do ia = 1,na
            do ix = 1,3
              xa(ix,ia) = alat * xa(ix,ia)
            enddo
          enddo
        elseif (iscale .eq. 3) then
          do ia = 1,na
            do ix = 1,3
              xac(ix) = xa(ix,ia)
            enddo
            do ix = 1,3
              xa(ix,ia) = ucell(ix,1) * xac(1) +
     .                    ucell(ix,2) * xac(2) +
     .                    ucell(ix,3) * xac(3)
            enddo
          enddo
        endif

      endif

c Expand the coordinates to whole supercell
      if (volcel(ucell) .lt. 1.d-8) then
        ncells = 1
      else
        ncells = nint( volcel(scell) / volcel(ucell) )
      endif
      if (ncells .gt. 1) then
        nua = na
        na  = na * ncells
        if (na.gt.0 .and. na.le.maxa) then

C         Find equivalent diagonal combination of cell/supercell
          call digcel( ucell, mscell, dcell, dscell, nsc, isdiag )

C         Expand coordinates
          call superx( dcell, nsc, nua, maxa, xa, dscell )

C         Expand index isa
          ia = 0
          do ic = 1,ncells
            do iua = 1,nua
              ia = ia + 1
              isa(ia) = isa(iua)
            enddo
          enddo

        endif
      endif

c Dump to output the coordinates is now done by siesta after ioxv

c     if (na .le. maxa) then
c       write(6,'(a)') 'recoor: Atomic coordinates (Bohr) and species'
c       do ia = 1,na
c         write(6,"('recoor: ',i4,2x,3f10.5,i3)")
c    .                    ia,(xa(ix,ia),ix=1,3),isa(ia)
c       enddo
c     endif

      return
      end
