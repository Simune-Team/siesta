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

C
C  Modules
C
      use precision
      use fdf
#ifdef MPI
      use mpi
#endif

      implicit          none

      character         task*(*), paste*33
      logical           found
      integer           na, isa(na), iza(na)
      double precision  cell(3,3), va(3,na), vcell(3,3), xa(3,na)
      external          io_assign, io_close, paste

c Internal variables and arrays
      character  sname*30, fname*33
      integer    ia, iu, iv, ix, Node
#ifdef MPI
      integer    MPIerror
#endif
      logical    frstme
      save       frstme, fname
      data frstme /.true./

c Set Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

c Only do reading and writing for IOnode
      if (Node.eq.0) then

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
     .       'ioxv: Reading coordinates and velocities from file'
            do iv = 1,3
              read(iu,*) (cell(ix,iv),ix=1,3),(vcell(ix,iv),ix=1,3)
            enddo
            read(iu,*) na
            do ia = 1,na
              read(iu,*)
     .          isa(ia),iza(ia),(xa(ix,ia),ix=1,3),(va(ix,ia),ix=1,3)
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
     .      ((cell(ix,iv),ix=1,3),(vcell(ix,iv),ix=1,3),iv=1,3)
          write(iu,*) na
          do ia = 1,na
            write(iu,'(i3,i6,3f18.9,3x,3f18.9)')
     .        isa(ia),iza(ia),(xa(ix,ia),ix=1,3),(va(ix,ia),ix=1,3)
          enddo

c       Close file
          call io_close( iu )

        endif
      endif

  999 continue

c If data has been read in then broadcast the values to all Nodes
#ifdef MPI
      call MPI_Bcast(found,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      if (found.and.(task.eq.'read' .or. task.eq.'READ')) then
        call MPI_Bcast(na,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#ifdef NODAT
        call MPI_Bcast(cell(1,1),9,MPI_double_precision,0,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(vcell(1,1),9,MPI_double_precision,0,
     .    MPI_Comm_World,MPIerror)
#else
        call MPI_Bcast(cell(1,1),9,DAT_double,0,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(vcell(1,1),9,DAT_double,0,
     .    MPI_Comm_World,MPIerror)
#endif
        call MPI_Bcast(isa,na,MPI_integer,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(iza,na,MPI_integer,0,MPI_Comm_World,MPIerror)
#ifdef NODAT
        call MPI_Bcast(xa(1,1),3*na,MPI_double_precision,0,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(va(1,1),3*na,MPI_double_precision,0,
     .    MPI_Comm_World,MPIerror)
#else
        call MPI_Bcast(xa(1,1),3*na,DAT_double,0,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(va(1,1),3*na,DAT_double,0,
     .    MPI_Comm_World,MPIerror)
#endif
      endif
#endif

      return
      end

