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
      SUBROUTINE BONDS( CELL, NA, ISA, XA, RMAX, filename )

C **********************************************************************
C Finds the bond lengths (actually the distances to neighbors)
C Based on shaper (by J.M.Soler)
C Alberto Garcia, Oct 2006
C ************ INPUT ***************************************************
C REAL*8  CELL(3,3) : Lattice vectors CELL(Ixyz,Ivector)
C INTEGER NA        : Number of atoms
C INTEGER ISA(NA)   : Species index of each atom
C REAL*8  XA(3,NA)  : Cartesian atomic coordinates
C REAL*8  RMAX      : Maximum distance considered
C CHARACTER(LEN=*)  filename : File name for output
C ************ UNITS ***************************************************
C CELL and XA must be in the same units
C **********************************************************************

      use precision, only : dp
      use atmfuncs,  only : labelfis
      use units,     only : Ang
      use m_recipes, only : sort
      use sorting,   only : order, iorder

      implicit          none

      INTEGER, intent(in)  ::     NA, ISA(NA)
      real(dp), intent(in) ::     CELL(3,3), XA(3,NA), RMAX
      character(len=*), intent(in) :: filename

      EXTERNAL             ::     NEIGHB, io_assign, io_close


C Internal variables and arrays
      integer, save ::
     .  MAXNA

      integer 
     .  IA, IN, IS, JA, JS, NNA, maxnain, iu

      integer, dimension(:), allocatable, save :: 
     .  JAN, index

      real(dp)
     .  RI, RIJ, RJ

      real(dp), dimension(:), allocatable, save :: 
     .  R2IJ

      real(dp), dimension(:,:), allocatable, save :: 
     .  XIJ

      logical
     .  overflow

      real(dp) :: tol = 1.0e-8_dp

      data MAXNA / 1000 /

C Allocate local memory
  100 if (.not.allocated(jan)) then
        allocate(jan(maxna))
        call memory('A','I',maxna,'shaper')
      endif
      if (.not.allocated(index)) then
        allocate(index(maxna))
        call memory('A','I',maxna,'shaper')
      endif
      if (.not.allocated(r2ij)) then
        allocate(r2ij(maxna))
        call memory('A','D',maxna,'shaper')
      endif
      if (.not.allocated(xij)) then
        allocate(xij(3,maxna))
        call memory('A','D',3*maxna,'shaper')
      endif

C Initialize neighbour-locater routine
      NNA = MAXNA
      CALL NEIGHB( CELL, RMAX, NA, XA, 0, 0, NNA, JAN, XIJ, R2IJ )
      overflow = (nna.gt.maxna)
      if (overflow) then
        call memory('D','I',size(jan),'shaper')
        deallocate(jan)
        call memory('D','I',size(index),'shaper')
        deallocate(index)
        call memory('D','D',size(r2ij),'shaper')
        deallocate(r2ij)
        call memory('D','D',size(xij),'shaper')
        deallocate(xij)
        maxna = nna + nint(0.1*nna)
        goto 100
      endif

C Main loop

      call io_assign( iu )
      open(iu,file=filename,form='formatted',
     $     status='replace', action="write")      
      maxnain = maxna
      do IA = 1,NA

C Find neighbours of atom IA
        NNA = maxnain
        CALL NEIGHB( CELL, RMAX, NA, XA, IA, 0,
     .               NNA, JAN, XIJ, R2IJ )
        if (.not.overflow) overflow = (nna.gt.maxna)
        if (overflow) then
          maxna = max(maxna,nna)
        endif

C Don't do the actual work if the neighbour arrays are too small
        if (.not.overflow) then
           if (nna < 2 ) then
              write(iu , *) "Atom ", ia, 
     $             ": no bonds for rmax specified: ", rmax/Ang,
     $             " Ang."
              call pxfflush(iu)
              cycle   ! loop over ia
           endif
           call sort( nna, r2ij, index )
           call iorder( jan, 1, nna, index )
           call order(  r2ij, 1, nna, index )

           write(iu,fmt="(a,i3,1x,a,3f8.4)")
     $       "Neighbors of: ",
     $          ia, trim(labelfis(isa(ia))) // " at: ", xa(:,ia)

          do IN = 1, NNA
            JA = JAN(IN)
            JS = ISA(JA)
            RIJ = SQRT(R2IJ(IN))
            if (rij > 0.0001_dp)
     $        write(iu,fmt="(i3,1x,a,f8.4,2x,a,3f8.4)")
     $           ja, labelfis(js), rij/Ang, "Ang. at: ", xa(:,ja)
          enddo
        endif

      enddo

C If maxna dimension was no good then reset arrays and start again
      if (overflow) then
        call memory('D','I',size(jan),'shaper')
        deallocate(jan)
        call memory('D','I',size(index),'shaper')
        deallocate(index)
        call memory('D','D',size(r2ij),'shaper')
        deallocate(r2ij)
        call memory('D','D',size(xij),'shaper')
        deallocate(xij)
        maxna = maxna + nint(0.1*maxna)
        goto 100
      endif

C Deallocate local memory
      call memory('D','I',size(jan),'shaper')
      deallocate(jan)
      call memory('D','I',size(index),'shaper')
      deallocate(index)
      call memory('D','D',size(r2ij),'shaper')
      deallocate(r2ij)
      call memory('D','D',size(xij),'shaper')
      deallocate(xij)

      call io_close(iu)

      end
