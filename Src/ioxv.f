c $Id: ioxv.f,v 1.5 1999/01/31 11:14:58 emilio Exp $

      subroutine ioxv( task, cell, vcell,
     .                 na, isa, iza, xa, va, found )

c *******************************************************************
c Saves positions and velocities.
c J.M.Soler. July 1997.
c ********** INPUT **************************************************
c character task*(*) : 'read' or 'write' (or 'READ' or 'WRITE')
c ********** INPUT OR OUTPUT (depending of task) ********************
c real*8  cell(3,3)  : Unit cell vectors
c real*8  vcell(3,3) : Unit cell vector velocities (Parrinello-Rahman)
c integer na         : Number of atoms
c integer isa(na)    : Atomic species index
c integer iza(na)    : Atomic numbers
c real*8  xa(3,na)   : Atomic positions
c real*8  va(3,na)   : Atomic velocities
c ********** OUTPUT *************************************************
c logical found      : Has input file been found
c                      (only for task='read')
c ********** UNITS **************************************************
c Units are arbitrary, but the use with task='write' and task='read'
c must be consistent
c *******************************************************************

      implicit          none
      character         task*(*), paste*33
      logical           found
      integer           na, isa(na), iza(na)
      double precision  cell(3,3), va(3,na), vcell(3,3), xa(3,na)
      external          io_assign, io_close, paste
      include          'fdf/fdfdefs.h'

c Internal variables and arrays
      character  sname*30, fname*33
      integer    ia, iu, iv, ix
      logical    frstme
      save       frstme, fname
      data frstme /.true./


c Find name of file
      if (frstme) then
        sname = fdf_string( 'SystemLabel', 'siesta' )
        fname = paste( sname, '.XV' )
        frstme = .false.
      endif

c Choose between read or write
      if (task.eq.'read' .or. task.eq.'READ') then

c       Check if input file exists
        inquire( file=fname, exist=found )
        if (found) then

c         Open file
          call io_assign( iu )
          open( iu, file=fname, status='old' )      

c         Read data
          write(6,'(/,a)') 
     .     'ioxv: Reading coordinates and velocities from file'
          do iv = 1,3
            read(iu,*) (cell(ix,iv),ix=1,3),(vcell(ix,iv),ix=1,3)
          enddo
          read(iu,*) na
          do ia = 1,na
            read(iu,*)
     .        isa(ia),iza(ia),(xa(ix,ia),ix=1,3),(va(ix,ia),ix=1,3)
          enddo

c         Close file
          call io_close( iu )

        else
c         If input file not found, go to exit point
          goto 999
        endif

      elseif (task.eq.'write' .or. task.eq.'WRITE') then

c       Open file
        call io_assign( iu )
        open( iu, file=fname, form='formatted', status='unknown' )      

c       Write data on file
        write(iu,'(2(3x,3f18.9))')
     .    ((cell(ix,iv),ix=1,3),(vcell(ix,iv),ix=1,3),iv=1,3)
        write(iu,*) na
        do ia = 1,na
          write(iu,'(i3,i6,3f18.9,3x,3f18.9)')
     .      isa(ia),iza(ia),(xa(ix,ia),ix=1,3),(va(ix,ia),ix=1,3)
        enddo

c       Close file
        call io_close( iu )

      endif

c Exit point
  999 continue
      end

