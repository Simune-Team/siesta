      subroutine iohs( task, no, nspin, maxo, maxno,
     .                 numh, listh, H, S, qtot, temp, xij )
C *********************************************************************
C Saves the hamiltonian and overlap matrices, and other data required
C to obtain the bands and density of states
C Writen by J.Soler July 1997.
C *************************** INPUT **********************************
C character*(*) task          : 'read'/'READ' or 'write'/'WRITE'
C ******************** INPUT or OUTPUT (depending on task) ***********
C integer no                  : Number of basis orbitals
C integer nspin               : Spin polarization (1 or 2)
C integer maxo                : Maximum number of basis  orbitals
C integer maxno               : Maximum number of orbitals interacting
C                               with any orbital
C integer numh(no)            : Number of nonzero elements of each row
C                               of hamiltonian matrix
C integer listh(maxno,no)     : Nonzero hamiltonian-matrix element column
C                               indexes for each matrix row
C real*8  H(maxno,maxo,nspin) : Hamiltonian in sparse form
C real*8  S(maxno,maxo)       : Overlap in sparse form
C real*8  qtot                : Total number of electrons
C real*8  temp                : Electronic temperature for Fermi smearing
C real*8  xij(3,maxno,maxo)   : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C *************************** UNITS ***********************************
C Units should be consistent between task='read' and 'write'
C *********************************************************************

      implicit          none
      character         task*(*), paste*33
      integer           maxo, maxno, no, nspin
      integer           listh(maxno,no), numh(no)
      double precision  H(maxno,maxo,nspin), S(maxno,maxo),
     .                  qtot, temp, xij(3,maxno,maxo)
      external          io_assign, io_close, paste
      include          'fdf/fdfdefs.h'

c Internal variables and arrays
      character  sname*30, fname*33
      integer    iu, ju, mno, mo, ns
      logical    baddim, found, frstme
      save       frstme, fname
      data frstme /.true./

c Find name of file
      if (frstme) then
        sname = fdf_string( 'SystemLabel', 'siesta' )
        fname = paste( sname, '.HS' )
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

c         Read dimensions
          read(iu) no, ns, mo, mno

c         Check dimensions
          baddim = .false.
          if (ns  .ne. nspin) baddim = .true.
          if (mo  .ne. maxo)  baddim = .true.
          if (mno .ne. maxno) baddim = .true.
          if (baddim) then
            call io_assign( ju )
            open( ju, file='iohs.h', status='unknown' )
            write(ju,'(a)') 'C Dimensions for input to iohs'
            write(ju,'(6x,a,i8,a)') 'parameter ( nspin =', ns,  ' )'
            write(ju,'(6x,a,i8,a)') 'parameter ( maxo  =', maxo,  ' )'
            write(ju,'(6x,a,i8,a)') 'parameter ( maxno =', maxno, ' )'
            call io_close( ju )
            stop 'iohs: BAD DIMENSIONS'
          endif

c         Read data
          read(iu) numh, listh, H, S, qtot, temp, xij

c         Close file
          call io_close( iu )

        else
          write(6,*) 'iohs: ERROR: file not found: ', fname
          stop 'iohs: ERROR: file not found'
        endif

      elseif (task.eq.'write' .or. task.eq.'WRITE') then

c       Open file
        call io_assign( iu )
        open( iu, file=fname, form='unformatted', status='unknown' )      

c       Write data
        write(iu) no, nspin, maxo, maxno
        write(iu) numh, listh, H, S, qtot, temp, xij

c       Close file
        call io_close( iu )
      endif
      end

