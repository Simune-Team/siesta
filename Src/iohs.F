      subroutine iohs( task, gamma, no, nspin, maxo, maxno,
     .                 numh, listh, H, S, qtot, temp, xij )
C *********************************************************************
C Saves the hamiltonian and overlap matrices, and other data required
C to obtain the bands and density of states
C Writen by J.Soler July 1997.
C *************************** INPUT **********************************
C character*(*) task          : 'read'/'READ' or 'write'/'WRITE'
C logical       gamma         : Is only gamma point used?
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
C                               (not read/written if only gamma point)
C *************************** UNITS ***********************************
C Units should be consistent between task='read' and 'write'
C ********* PARALLEL **************************************************
C In order to preserve the order of the density matrix this
C routine reads by reading in strict order the part of the
C density matrix local to the current processor.
C *********************************************************************

C
C  Modules
C
      use precision
      use parallel
      use fdf
#ifdef MPI
      use mpi
#endif

      implicit          none

      character         task*(*), paste*33
      logical           gamma
      integer           maxo, maxno, no, nspin
      integer           listh(maxno,no), numh(no)
      double precision  H(maxno,maxo,nspin), S(maxno,maxo),
     .                  qtot, temp, xij(3,maxno,maxo)
      external          io_assign, io_close, paste

c Internal variables and arrays
      character  sname*30, fname*33
      integer    iu, ju, mno, mo, ns, j, k, is
      integer    Node,Nodes,ih,hl,maxotot,notot,maxoloc,noloc
#ifdef MPI
      integer    MPIerror
#endif
      logical    baddim, found, frstme
      save       frstme, fname
      data frstme /.true./

c Get the Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif

c Find name of file
      if (frstme) then
        if (Node.eq.0) then
          sname = fdf_string( 'SystemLabel', 'siesta' )
        endif
#ifdef MPI
        call MPI_Bcast(sname,30,MPI_character,0,MPI_Comm_World,
     .    MPIerror)
#endif
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
          if (Node.eq.0) then
            read(iu) noloc, ns, mo, mno
#ifdef MPI
            call MPI_Bcast(noloc,1,MPI_integer,0,MPI_Comm_World,
     .        MPIerror)
            call GetNodeOrbs(noloc,Node,Nodes,no)
            call MPI_Bcast(ns,1,MPI_integer,0,MPI_Comm_World,MPIerror)
            call MPI_Bcast(mo,1,MPI_integer,0,MPI_Comm_World,MPIerror)
            call MPI_Bcast(mno,1,MPI_integer,0,MPI_Comm_World,MPIerror)
            call MPI_AllReduce(maxo,maxoloc,1,MPI_integer,MPI_sum,
     .        MPI_Comm_World,MPIerror)
#else
            no = noloc
            maxoloc = maxo
#endif
          endif

c         Check dimensions
          baddim = .false.
          if (ns  .ne. nspin) baddim = .true.
          if (mo  .ne. maxoloc)  baddim = .true.
          if (mno .ne. maxno) baddim = .true.
          if (baddim) then
            if (Node.eq.0) then
              call io_assign( ju )
              open( ju, file='iohs.h', status='unknown' )
              write(ju,'(a)') 'C Dimensions for input to iohs'
              write(ju,'(6x,a,i8,a)') 'parameter ( nspin =', ns,  ' )'
              write(ju,'(6x,a,i8,a)') 'parameter ( maxo  =', maxo,  ' )'
              write(ju,'(6x,a,i8,a)') 'parameter ( maxno =', maxno, ' )'
              call io_close( ju )
              stop 'iohs: BAD DIMENSIONS'
            else
              stop
            endif
          endif

c Read data in a processor independent fashion
          do ih = 1,no
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            if (hl.gt.0) then
              read(iu) numh(ih)
            endif
#ifdef MPI
            call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
          enddo
          do ih = 1,no
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            if (hl.gt.0) then
              do j=1,numh(ih)
                read(iu) listh(j,ih)
              enddo
            endif
#ifdef MPI
            call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
          enddo
          do is = 1,nspin
            do ih = 1,no
              call GlobalToLocalOrb(ih,Node,Nodes,hl)
              if (hl.gt.0) then
                do j=1,numh(ih)
                  read(iu) H(j,hl,is)
                enddo
              endif
#ifdef MPI
              call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
            enddo
          enddo
          do ih = 1,no
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            if (hl.gt.0) then
              do j=1,numh(ih)
                read(iu) S(j,hl)
              enddo
            endif
#ifdef MPI
            call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
          enddo
          
          if (Node.eq.0) then
            read(iu) qtot,temp
#ifdef MPI
            call MPI_Bcast(qtot,1,DAT_double,0,
     .        MPI_Comm_World,MPIerror)
            call MPI_Bcast(temp,1,DAT_double,0,
     .        MPI_Comm_World,MPIerror)
#endif
          endif

          if (.not.gamma) then
            do ih = 1,no
              call GlobalToLocalOrb(ih,Node,Nodes,hl)
              if (hl.gt.0) then
                do j=1,numh(ih)
                  read(iu) (xij(k,j,hl),k=1,3)
                enddo
              endif
#ifdef MPI
              call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
            enddo
          endif

c         Close file
          call io_close( iu )

        else
          if (Node.eq.0) then
            write(6,*) 'iohs: ERROR: file not found: ', fname
            stop 'iohs: ERROR: file not found'
          else
            stop 
          endif
        endif

      elseif (task.eq.'write' .or. task.eq.'WRITE') then

c       Open file
        call io_assign( iu )
        open( iu, file=fname, form='unformatted', status='unknown' )      

c Find total numbers over all Nodes
#ifdef MPI
        call MPI_AllReduce(no,notot,1,MPI_integer,MPI_sum,
     .    MPI_Comm_World,MPIerror)
        call MPI_AllReduce(maxo,maxotot,1,MPI_integer,MPI_sum,
     .    MPI_Comm_World,MPIerror)
#else
        notot = no
        maxotot = maxo
#endif

c Write overall data
        if (Node.eq.0) then
          write(iu) notot, nspin, maxotot, maxno
        endif

c Write data in a processor independent fashion
        do ih = 1,no
          call GlobalToLocalOrb(ih,Node,Nodes,hl)
          if (hl.gt.0) then
            write(iu) numh(ih)
          endif
#ifdef MPI
          call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
        enddo
        do ih = 1,no
          call GlobalToLocalOrb(ih,Node,Nodes,hl)
          if (hl.gt.0) then
            do j=1,numh(ih)
              write(iu) listh(j,ih)
            enddo
          endif
#ifdef MPI
          call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
        enddo
        do is = 1,nspin
          do ih = 1,no
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            if (hl.gt.0) then
              do j=1,numh(ih)
                write(iu) H(j,hl,is)
              enddo
            endif
#ifdef MPI
            call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
          enddo
        enddo
        do ih = 1,no
          call GlobalToLocalOrb(ih,Node,Nodes,hl)
          if (hl.gt.0) then
            do j=1,numh(ih)
              write(iu) S(j,hl)
            enddo
          endif
#ifdef MPI
          call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
        enddo

        write(iu) qtot,temp

        if (.not.gamma) then
          do ih = 1,no
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            if (hl.gt.0) then
              do j=1,numh(ih)
                write(iu) (xij(k,j,hl),k=1,3)
              enddo
            endif
#ifdef MPI
            call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
          enddo
        endif

c       Close file
        call io_close( iu )
      endif
      end

