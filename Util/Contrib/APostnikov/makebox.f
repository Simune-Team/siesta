C
      subroutine makebox(obox,rbox)
C
C     asks for origin ans spanning vectors, 
C     constructs the output box and its inversion.
C
C     Input:   none (interactive in/out)
C     Output:  obox - origin of output box,
C              rbox - its three spanning vectors
C
      implicit none
      integer ii
      double precision obox(3),rbox(3,3),b2ang
      parameter (b2ang=0.529177)   !  Bohr to Angstroem
      character unitlab*1,labunit*9
      logical unitb

      write (6,702)
  702 format(" Now define the grid cell for your XCrysDen plot.",/
     .       " Note that it can be arbitrarily chosen",
     .       " with respect to the Siesta simulation cell,",
     .       " and it needs not to be orthogonal.",
     .       " We'll define it by the origin point and three",
     .       " spanning vectors. They can be given in Bohr or Ang.")
  101 write (6,703,advance="no")
  703 format (" Would you use Bohr (B) or Ang (A) ? ")
      read (5,*) unitlab 
      if (unitlab.eq.'B'.or.unitlab.eq.'b') then
        unitb = .true.
        labunit = 'Bohr'
      elseif (unitlab.eq.'A'.or.unitlab.eq.'a') then
        unitb = .false.
        labunit = 'Angstroem'
      else 
        write (6,*) ' Sorry, no third choice.' 
        goto 101
      endif
      write (6,*) ' Enter origin point in ',labunit,':'
      read (5,*) (obox(ii),ii=1,3)
      write (6,*) ' Enter first spanning vector in ',
     .              labunit(1:len_trim(labunit)),':'
      read (5,*) (rbox(ii,1),ii=1,3)
      write (6,*) ' Enter second spanning vector in ',
     .              labunit(1:len_trim(labunit)),':'
      read (5,*) (rbox(ii,2),ii=1,3)
      write (6,*) ' Enter third spanning vector in ',
     .              labunit(1:len_trim(labunit)),':'
      read (5,*) (rbox(ii,3),ii=1,3)
      if (unitb) then
C       transform anything into Ang, as is standard in XCrysden:
        obox=obox*b2ang
        rbox=rbox*b2ang
      endif
      return
      end
