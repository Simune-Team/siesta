      subroutine coceri(iza, xa, cell, na, sname, slabel)

c *******************************************************************
c Writes coordinates in format to be read by CERIUS
c
c It implies atomic symbols, atomic coordinates (fractional format)
c and lattice parameters (modules (in Ang) and angles)
c
c Written by E. Artacho. December 1997.
c ********* INPUT ***************************************************
c integer   iza(na)   : Atomic numbers of different atoms
c double    xa(3,na)  : Atom coordinates (in Bohr)
c double    cell(3,3) : Lattice vectors (in Bohr)
c integer   na        : Number of atoms
c character slabel*20 : Label for file naming
c character sname*150 : Label for the title  
c ******************************************************************

      implicit          none
      character         slabel*20, sname*150, paste*24
      integer           na
      integer           iza(na)
      double precision  xa(3,na), cell(3,3)
      external          io_assign, io_close, paste, symbol

c Internal variables and arrays
 
      character         fname*24, symbol*2
      integer           unit,ix, iv,  i, ia
      double precision  celang(3), cellm(3), recell(3,3),
     .                  xac(3), pi, Ang 

      double precision, dimension(:,:), allocatable ::
     .                  xap

      data pi, Ang      / 3.1415926d0, 0.529177d0 /

C Allocate local memory
      allocate(xap(3,na))
      call memory('A','D',3*na,'coceri')
C ..................

c Find lattice parameters out of lattice vectors: first modules:

      do iv = 1, 3
         cellm(iv) = 0.d0
         do ix = 1, 3
            cellm(iv) = cellm(iv) + cell(ix,iv)*cell(ix,iv)
         enddo
         cellm(iv) = sqrt(cellm(iv))
      enddo

c and angles

      celang(1) = 0.d0
      do ix = 1, 3
         celang(1) = celang(1) + cell(ix,2)*cell(ix,3)
      enddo
      celang(1) = acos(celang(1)/(cellm(2)*cellm(3)))*180.d0/pi
      celang(2) = 0.d0
      do ix = 1, 3
         celang(2) = celang(2) + cell(ix,1)*cell(ix,3)
      enddo
      celang(2) = acos(celang(2)/(cellm(1)*cellm(3)))*180.d0/pi
      celang(3) = 0.d0
      do ix = 1, 3
         celang(3) = celang(3) + cell(ix,1)*cell(ix,2)
      enddo
      celang(3) = acos(celang(3)/(cellm(1)*cellm(2)))*180.d0/pi

c Obtain fractional coordinates (reclat inverts matrix)

      call reclat(cell, recell, 0)
      do ia = 1,na
        do ix = 1,3
          xac(ix) = xa(ix,ia)
        enddo
        do ix = 1,3
          xap(ix,ia) = recell(1,ix) * xac(1) +
     .                 recell(2,ix) * xac(2) +
     .                 recell(3,ix) * xac(3)
        enddo
      enddo


c Find file name

      fname = paste(slabel,'.xtl')

      write(6,'(/,2a)')'coceri: Writing CERIUS coordinates into file ',
     .                  fname

      call io_assign(unit)
      open( unit, file=fname, form = 'formatted', status='unknown')
      rewind(unit)

c Write file

      write(unit,'(a,a70)') 'TITLE ', sname
      write(unit,'(a)')  'DIMENSION 3'
      write(unit,'(a,6f11.5)') 
     .          'CELL', (cellm(iv)*Ang,iv=1,3), (celang(i),i=1,3)
      write(unit,'(a)') 'SYMMETRY  NUMBER 1  LABEL P1'
      write(unit,'(3a)') 
     .       'SYM MAT  1.000000  0.000000  0.000000  0.000000',
     .              '  1.000000  0.000000  0.000000  0.000000',
     .              '  1.000000 0.0000 0.0000 0.0000'
      write(unit,'(/,a)') 'ATOMS'
      write(unit,'(2a)') 'NAME       X          Y          Z     ',
     .                             'CHARGE   TEMP    OCCUP   SCAT'
      write(unit,'(2x,a2,3f11.5,3f8.4,3x,a2)')
     .  ( symbol(iza(ia)), (xap(i,ia),i=1,3),
     .                0.d0, 0.d0, 1.d0, symbol(iza(ia)), ia=1,na )
      write(unit,'(a)') 'EOF'

      call io_close(unit)
      
C Deallocate local memory
      call memory('D','D',size(xap),'coceri')
      deallocate(xap)

      return
      end

