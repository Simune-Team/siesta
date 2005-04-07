      subroutine outcell(cell) 

c *******************************************************************
c Writes lattice vectors and lattice parameters
c
c Written by E. Artacho, December 1997.
c ********* INPUT ***************************************************
c double precision cell(3,3): Lattice (supercell) vectors
c *******************************************************************

      use precision, only : dp
      use siesta_cml
      use units, only : Ang, pi

      implicit none

      real(dp)   :: cell(3,3)

c Internal variables and arrays

      integer    :: iv, ix
      real(dp)   :: cellm(3), celang(3), 
     .              volume, volcel

      external volcel

c Writing cell vectors

      write(6,'(/,a,3(/,a,3f12.6))')
     . 'outcell: Unit cell vectors (Ang):',
     . ('    ', (cell(ix,iv)/Ang,ix=1,3), iv =1,3)

c Cell-vector modules

      do iv = 1,3
        cellm(iv) = 0.0_dp
        do ix = 1,3
          cellm(iv) = cellm(iv) + cell(ix,iv)*cell(ix,iv)
        enddo
        cellm(iv) = sqrt(cellm(iv))
      enddo

c Cell-vector angles

      celang(1) = 0.0_dp
      do ix = 1,3
         celang(1) = celang(1) + cell(ix,1)*cell(ix,2)
      enddo
      celang(1) = acos(celang(1)/(cellm(1)*cellm(2)))*180._dp/pi
      celang(2) = 0.0_dp
      do ix = 1,3
         celang(2) = celang(2) + cell(ix,1)*cell(ix,3)
      enddo
      celang(2) = acos(celang(2)/(cellm(1)*cellm(3)))*180._dp/pi
      celang(3) = 0.0_dp
      do ix = 1,3
         celang(3) = celang(3) + cell(ix,2)*cell(ix,3)
      enddo
      celang(3) = acos(celang(3)/(cellm(2)*cellm(3)))*180._dp/pi
      write(6,'(/,a,3f12.6)')
     . 'outcell: Cell vector modules (Ang)   :',
     .          (cellm(iv)/Ang,iv=1,3)
      write(6,'(a,3f12.4)')
     . 'outcell: Cell angles (23,13,12) (deg):',
     .          (celang(iv),iv=3,1,-1)

c Cell volume

      volume = volcel( cell )
      write(6,'(a,f12.4)')
     .  'outcell: Cell volume (Ang**3)        :', volume/Ang**3

      if (cml_p) then
        call cmlAddCrystal(xf=mainXML, 
     .       title='Lattice Parameters', fmt='(f12.6)', 
     .       alpha=celang(1), beta=celang(2), gamma=celang(3),
     .       a=cellm(1)/Ang, b=cellm(2)/Ang, c=cellm(3)/Ang )
      endif

      end subroutine outcell
