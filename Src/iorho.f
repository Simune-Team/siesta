C $Id: iorho.f,v 1.5 1999/01/31 11:14:58 emilio Exp $

      subroutine iorho( task, fname, cell, mesh, maxp, nspin, rho,
     .                  found )

C *********************************************************************
C Saves the electron density at the mesh points.
C Writen by J.Soler July 1997.
C *************************** INPUT **********************************
C character*(*) task      : 'read'/'READ' or 'write'/'WRITE'
C character*(*) fname     : File name for input or output
C integer maxp            : First dimension of array rho
C integer nspin           : Number of spin polarizations (1 or 2)
C ******************** INPUT or OUTPUT (depending on task) ***********
C real*8  cell(3,3)       : Lattice vectors
C integer mesh(3)         : Number of mesh divisions of each
C                           lattice vector
C real    rho(maxp,nspin) : Electron density
C                           Notice single precision in this version
C ******************** OUTPUT *****************************************
C logical found           : Were data found? (only when task='read')
C *************************** UNITS ***********************************
C Units should be consistent between task='read' and 'write'
C *********************************************************************

      implicit          none
      character*(*)     fname, task
      integer           maxp, mesh(3), nspin
      real              rho(maxp,nspin)
      double precision  cell(3,3)
      external          io_assign, io_close, paste
      include          'fdf/fdfdefs.h'

c Internal variables and arrays
      character  paste*33, ffor*9, fform*(*)
c     character  sname*30
      integer    i, ip, iu, is, j, ju, np, ns
      logical    baddim, found

c Fix whether formatted or unformatted files wil be used
      parameter ( fform = 'unformatted' )

c Choose between read or write
      if (task.eq.'read' .or. task.eq.'READ') then

c       Check if input file exists
        inquire( file=fname, exist=found )
        if (found) then

c         Open file
          call io_assign( iu )
          open( iu, file=fname, form=fform, status='old' )      

c         Read cell vectors and number of points
          ffor = 'formatted'
          if (fform .eq. ffor) then
            read(iu) cell
            read(iu) mesh, ns
          else
            read(iu,*) cell
            read(iu,*) mesh, ns
          endif
          np = mesh(1) * mesh(2) * mesh(3)

c         Check dimensions
          baddim = .false.
          if (ns .ne. nspin) baddim = .true.
          if (np .gt. maxp)  baddim = .true.
          if (baddim) then
            call io_assign( ju )
            open( ju, file='iorho.h', status='unknown' )
            write(ju,'(a)') 'C Dimensions for input to iorho'
            write(ju,'(6x,a,i8,a)') 'parameter ( nspin =', ns, ' )'
            write(ju,'(6x,a,i8,a)') 'parameter ( maxp  =', np, ' )'
            call io_close( ju )
            write(6,'(a)') 'iorho: ERROR: BAD DIMENSIONS'
            stop           'iorho: ERROR: BAD DIMENSIONS'
          endif

c         Read data
          ffor = 'formatted'
          if (fform .eq. ffor) then
            do is = 1,ns
              read(iu,*) (rho(ip,is),ip=1,np)
            enddo
          else
            do is = 1,ns
              read(iu) (rho(ip,is),ip=1,np)
            enddo
          endif

c         Close file
          call io_close( iu )

        endif

      elseif (task.eq.'write' .or. task.eq.'WRITE') then

c       Open file
        call io_assign( iu )
        open( iu, file=fname, form=fform, status='unknown' )      

c       Write data
        np = mesh(1) * mesh(2) * mesh(3)
        ffor = 'formatted'
        if (fform .eq. ffor) then
          do i = 1,3
            write(iu,*) (cell(j,i),j=1,3)
          enddo
          write(iu,*) mesh, nspin
          do is = 1,nspin
            write(iu,'(e15.6)') (rho(ip,is),ip=1,np)
          enddo
        else
          write(iu) cell
          write(iu) mesh, nspin
          do is = 1,nspin
            write(iu) (rho(ip,is),ip=1,np)
          enddo
        endif

c       Close file
        call io_close( iu )
      endif
      end


