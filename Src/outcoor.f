      subroutine outcoor(cell, xa, isa, na, cohead, writec) 

c *******************************************************************
c Writes atomic coordinates in format given by fdf:AtomCoorFormatOut
c Default: same format as was used for the reading.
c
c Written by E. Artacho, December 1997, on the original piece of the
c redata subroutine written by P. Ordejon in December 1996. 
C Modified by DSP, Aug. 1998.
c ********* INPUT ***************************************************
c integer na                : number of atoms
c double precision cell(3,3): Lattice (supercell) vectors
c double precision alat     : Lattice constant (in Bohr) 
c double precision xa(3,na) : atomic coordinates in Bohr cartesian
c integer isa(na)           : atomic species of different atoms
c character*(*) cohead      : special phrase to include in heading
c logical writec            : writing coor only if true
c *******************************************************************
c inversion of lattice vectors for fractional coordinates is done
c through subroutine reclat.
c *******************************************************************

      use atmfuncs, only: labelfis
      use fdf

      implicit          none
      integer           na
      integer           isa(na)
      logical           writec
      character         cohead*(*)
      double precision  xa(3,na), cell(3,3)

c Internal variables and arrays

      character         acf*22, acf_defect*22, acfout*22, 
     .                  pieceh*20, titl*60, pasteb*60
      logical           leqi, frstme
      integer           ia, ix
      double precision  xac(3), recell(3,3), alat

      double precision, dimension(:,:), allocatable ::
     .                  xap

      external          pasteb, memory
      save              frstme, acfout, alat

      data              frstme /.true./
c--------------------------------------------------------------------

      if (frstme) then

C read format for output of atomic coordinates

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
           write(6,"(/,2a)") 'outcoor: WARNING: Explicit lattice ',
     .         'constant is needed for ScaledCartesian output.'
           write(6,"(2a)")   'outcoor:          NotScaledCartesianAng',
     .         'format being used instead.'
           acfout = 'NotScaledCartesianAng'
        endif

        frstme = .false.
      endif

C Allocate local memory
      allocate(xap(3,na))
      call memory('A','D',3*na,'outcoor')

c Write coordinates at every time or relaxation step?

      if ( (cohead .eq. ' ') .and. ( .not. writec) ) return

c write coordinates according to format 

      if (leqi(acfout,'NotScaledCartesianBohr') .or. 
     .    leqi(acfout,'Bohr') ) then
        pieceh = '(Bohr):'
        do ia = 1,na
          do ix = 1,3
            xap(ix,ia) = xa(ix,ia)
          enddo
        enddo
      else if (leqi(acfout,'NotScaledCartesianAng') .or.
     .         leqi(acfout,'Ang') ) then
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
            xap(ix,ia) = recell(1,ix) * xac(1) +
     .                   recell(2,ix) * xac(2) +
     .                   recell(3,ix) * xac(3)
          enddo
        enddo
      else
        write(6,"(/,'outcoor: ',72(1h*))")
        write(6,"('outcoor:                  INPUT ERROR')")
        write(6,'(a)') 'outcoor: '
        write(6,'(2a)') 'outcoor: You must use one of the following',
     .                            ' coordinate output options:'
        write(6,'(a)') 'outcoor:     - NotScaledCartesianBohr (or Bohr)'
        write(6,'(a)') 'outcoor:     - NotScaledCartesianAng (or Ang) '
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
     .  ((xap(ix,ia),ix=1,3),isa(ia),labelfis(isa(ia)),ia,ia=1,na)

C Deallocate local memory
      call memory('D','D',size(xap),'outcoor')
      deallocate(xap)

      return
      end
