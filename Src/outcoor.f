      subroutine outcoor(cell, xa, isa, na, atm_label, cohead) 

c *******************************************************************
c Writes atomic coordinates in format given by fdf:AtomCoorFormatOut
c Default: same format as was used for the reading.
c
c Written by E. Artacho, December 1997, on the original piece of the
c redata subroutine written by P. Ordejon in December 1996.
c ********* INPUT ***************************************************
c integer na                : number of atoms
c double precision cell(3,3): Lattice (supercell) vectors
c double precision alat     : Lattice constant (in Bohr) 
c double precision xa(3,na) : atomic coordinates in Bohr cartesian
c integer isa(na)           : atomic species of different atoms
c character(ns)*20          : atomic labels
c character*(*) cohead      : special phrase to include in heading
c *******************************************************************
c inversion of lattice vectors for fractional coordinates is done
c through subroutine reclat.
c *******************************************************************

      implicit          none
      integer           na
      integer           isa(na)
      character         atm_label(*)*20, cohead*(*)
      double precision  xa(3,na), cell(3,3)

c Internal variables and arrays

      integer namax

      parameter (namax=2000)

      character         acf*22, acf_defect*22, acfout*22, 
     .                  pieceh*20, titl*60, pasteb*60
      logical           leqi
      integer           ia, ix
      double precision  xac(3), xap(3,namax), recell(3,3), alat

      external          pasteb

c enable FDF input/output

      include 'fdf/fdfdefs.h'


C Check dimensions ..........................................................
      if (na .gt. namax) then
        write(6,*) 'coceri: Wrong namax; Must be at least ',na
        stop
      endif
C ..................

C format for output of atomic coordinates

      acf_defect = 'NotScaledCartesianBohr'
      acf = fdf_string('AtomicCoordinatesFormat',acf_defect)
      acfout = fdf_string('AtomCoorFormatOut',acf)

c Scale atomic coordinates
c   Coord. option Bohr       => Do nothing
c   Coord. option Ang        => Multiply by 0.529177 (Bohr --> Ang)
c   Coord. option Scaled     => Divide by lattice constant
c   Coord. option Fractional => Multiply by inverse of lattice vectors

      alat = fdf_physical('LatticeConstant',0.d0,'Bohr')
      if (alat.eq.0.d0 .and. leqi(acfout,'ScaledCartesian')) then
         write(6,"(/,2a)") 'recoor: WARNING: Explicit lattice ',
     .       'constant is needed for ScaledCartesian output.'
         write(6,"(2a)")   'recoor:          NotScaledCartesianAng',
     .       'format being used instead.'
         acfout = 'NotScaledCartesianAng'
      endif

      if (leqi(acfout,'NotScaledCartesianBohr')) then
        pieceh = '(Bohr):'
        do ia = 1,na
          do ix = 1,3
            xap(ix,ia) = xa(ix,ia)
          enddo
        enddo
      else if (leqi(acfout,'NotScaledCartesianAng')) then
        pieceh = '(Ang):'
        do ia = 1,na
          do ix = 1,3
            xap(ix,ia) = 0.529177d0 * xa(ix,ia)
          enddo
        enddo
      else if (leqi(acfout,'ScaledCartesian')) then
        pieceh = '(scaled):'
        do ia = 1,na
          do ix = 1,3
            xap(ix,ia) = xa(ix,ia) / alat
          enddo
        enddo
      else if (leqi(acfout,'ScaledByLatticeVectors') .or. 
     .         leqi(acfout,'Fractional') ) then
        pieceh = '(fractional):'
        call reclat(cell, recell, 0)
        do ia = 1,na
          do ix = 1,3
            xac(ix) = xa(ix,ia)
          enddo
          do ix = 1,3
            xap(ix,ia) = recell(ix,1) * xac(1) +
     .                   recell(ix,2) * xac(2) +
     .                   recell(ix,3) * xac(3)
          enddo
        enddo
      else
        write(6,"(/,'outcoor: ',72(1h*))")
        write(6,"('outcoor:                  INPUT ERROR')")
        write(6,'(a)') 'outcoor: '
        write(6,'(2a)') 'outcoor: You must use one of the following',
     .                            ' coordinate output options:'
        write(6,'(a)') 'outcoor:     - NotScaledCartesianBohr         '
        write(6,'(a)') 'outcoor:     - NotScaledCartesianAng          '
        write(6,'(a)') 'outcoor:     - ScaledCartesian                '
        write(6,'(2a)') 'outcoor:     - ScaledByLatticeVectors ',
     .                                               '(or Fractional)'
        write(6,"('outcoor: ',72(1h*))")
        stop 'outcoor: ERROR: Wrong atomic-coordinate output format'
      endif

c writing a heading for the coordinates

      if (cohead .eq. ' ') then
         titl = pasteb( 'outcoor: Atomic coordinates', pieceh )
      else
         titl = pasteb( 'outcoor:', cohead )
         titl = pasteb( titl, 'atomic coordinates')
         titl = pasteb( titl, pieceh)
      endif

      write(6,'(/,a)') titl

c writing the coordinates

      write(6,'(3f14.8,i4,2x,a6,i5)')
     .  ((xap(ix,ia),ix=1,3),isa(ia),atm_label(isa(ia)),ia,ia=1,na)

      return
      end
