      subroutine iodm( task, maxnd, maxo, nbasis, nspin,
     .                 numd, listd, dm, found )
C *******************************************************************
C Reads/writes density matrix from/to file
C Written by P.Ordejon and J.M.Soler. May 1997.
C ********* INPUT ***************************************************
C character task*3 : 'read' or 'write'
C integer   maxnd  : First dimension of listd and dm
C integer   maxo   : Second dimension of dm (must be maxo.ge.nbasis)
C integer   nbasis : Number of atomic orbitals
C integer   nspin  : Number of spins (1 or 2)
C ********* INPUT OR OUTPUT (depending on task) *********************
C integer numd(nbasis)        : Control vector of DM matrix
C                               (number of nonzero elements of each row)
C integer listd(maxnd,nbasis) : Control vector of DM matrix
C                               (list of nonzero elements of each row)
C real*8  dm(maxnd,maxo,nspin): Density matrix
C ********* OUTPUT *************************************************
C logical found : Has DM been found in disk? (Only when task='read')
C ********* PARALLEL ***********************************************
C In order to preserve the order of the density matrix this
C routine reads by reading in strict order the part of the
C density matrix local to the current processor.
C ******************************************************************
      
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
      logical           found
      integer           maxnd, maxo, nbasis, nspin
      integer           listd(maxnd,*), numd(*)
      double precision  dm(maxnd,maxo,nspin)
      external          chkdim, io_assign, io_close, paste, timer,
     .                  memory

c     Internal variables and arrays
      character fname*33, sname*30
      logical   exist1, exist2, exist3, frstme
      integer   im, is, unit1, unit2, m, nb, ndmax, ns
      integer   Node,Nodes,nbasistot,ml
#ifdef MPI
      integer   MPIerror, Request, Status(MPI_Status_Size)
      integer   BNode
      double precision, dimension(:), allocatable :: buffer
#endif
      save      frstme, fname
      data      frstme /.true./

*     call timer( 'iodm', 1 )

c Get the Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif
c     Find file name
      if (frstme) then
        if (Node.eq.0) then
          sname = fdf_string('SystemLabel','siesta')
        endif
#ifdef MPI
        call MPI_Bcast(sname,30,MPI_character,0,MPI_Comm_World,
     .    MPIerror)
#endif
        fname = paste(sname,'.DM')
        frstme = .false.
      endif

c Find total number of basis functions over all Nodes
#ifdef MPI
      if (Nodes.gt.1) then
        call MPI_AllReduce(nbasis,nbasistot,1,MPI_integer,MPI_sum,
     .    MPI_Comm_World,MPIerror)
      else
        nbasistot = nbasis
      endif
#else
      nbasistot = nbasis
#endif

#ifdef MPI
c Create buffer array for transfering density matrix between nodes
      allocate(buffer(maxnd))
      call memory('A','D',maxnd,'iodm')
#endif

      if (task.eq.'read' .or. task.eq.'READ') then
        if (Node.eq.0) then
          inquire (file='SAVE.DM',     exist=exist1)
          inquire (file='SAVE.ctrlDM', exist=exist2)
          inquire (file=fname,         exist=exist3)
        endif

#ifdef MPI
c       Broadcast logicals so that all processors take the same route
        call MPI_Bcast(exist1,1,MPI_logical,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(exist2,1,MPI_logical,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(exist3,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

c       First look for old-format files
        if (exist1 .and. exist2) then
          if (Node.eq.0) then
            write(6,'(/,a)') 'iodm: Reading Density Matrix from files'
            call io_assign(unit1)
            call io_assign(unit2)
            open( unit1, file='SAVE.DM', status='unknown')
            open( unit2, file='SAVE.ctrlDM', status='unknown')
            rewind(unit1)
            rewind(unit2)
          endif

          ndmax = 0
          if (Node.eq.0) then
            do m = 1,nbasistot
              read(unit2,*) numd(m)
              ndmax = max( ndmax, numd(m) )
            enddo
          endif
#ifdef MPI
          call MPI_Bcast(numd,nbasistot,MPI_integer,0,MPI_Comm_World,
     .      MPIerror)
          call MPI_Bcast(ndmax,1,MPI_integer,0,MPI_Comm_World,
     .      MPIerror)
#endif
          call chkdim( 'iodm', 'maxnd', maxnd, ndmax, 1 )
          do m = 1,nbasistot
            call GlobalToLocalOrb(m,Node,Nodes,ml)
            if (ml.gt.0) then
              do im = 1,numd(m)
                read(unit2,*) listd(im,m)
              enddo
            endif
#ifdef MPI
            call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
          enddo
          do is = 1,nspin
            do m = 1,nbasistot
#ifdef MPI
              call WhichNodeOrb(m,Nodes,BNode)
              if (BNode.eq.0.and.Node.eq.BNode) then
                call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
                ml = m
#endif
                do im = 1,numd(m)
                  read(unit1,*) dm(im,ml,is)
                enddo
#ifdef MPI
              elseif (Node.eq.0) then
                do im = 1,numd(m)
                  read(unit1,*) buffer(im)
                enddo
                call MPI_ISend(buffer,numd(m),DAT_double,
     .            BNode,1,MPI_Comm_World,Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)
              elseif (Node.eq.BNode) then
                call GlobalToLocalOrb(m,Node,Nodes,ml)
                call MPI_IRecv(dm(1,ml,is),numd(m),
     .            DAT_double,0,1,
     .            MPI_Comm_World,Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)
              endif
              if (BNode.ne.0) then
                call MPI_Barrier(MPI_Comm_World,MPIerror)
              endif
#endif
            enddo
          enddo

          if (Node.eq.0) then
            call io_close(unit1)
            call io_close(unit2)
          endif
          found = .true.

c       Look now for new-format files
        elseif (exist3) then
          if (Node.eq.0) then
            write(6,'(/,a)') 'iodm: Reading Density Matrix from files'
            call io_assign(unit1)
            open( unit1, file=fname,
     .          form='unformatted', status='unknown' )
            rewind(unit1)
            read(unit1) nb, ns
          endif

c Communicate the values to all Nodes and adjust to allow for
c distributed memory before checking the dimensions
#ifdef MPI
          call MPI_Bcast(nb,1,MPI_integer,0,MPI_Comm_World,MPIerror)
          call MPI_Bcast(ns,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

c Check dimensions
          call chkdim( 'iodm', 'nbasis', nbasistot, nb, 0 )
          call chkdim( 'iodm', 'nspin',  nspin,  ns, 0 )

          ndmax = 0
          if (Node.eq.0) then
            do m = 1,nbasistot
              read(unit1) numd(m)
              ndmax = max( ndmax, numd(m) )
            enddo
          endif
#ifdef MPI
          call MPI_Bcast(numd,nbasistot,MPI_integer,0,MPI_Comm_World,
     .      MPIerror)
          call MPI_Bcast(ndmax,1,MPI_integer,0,MPI_Comm_World,
     .      MPIerror)
#endif
          call chkdim( 'iodm', 'maxnd', maxnd, ndmax, 1 )
          do m = 1,nbasistot
            if (Node.eq.0) then
              do im = 1,numd(m)
                read(unit1) listd(im,m)
              enddo
            endif
#ifdef MPI
            call MPI_Bcast(listd(1,m),numd(m),MPI_integer,0,
     .        MPI_Comm_World,MPIerror)
#endif
          enddo
          do is = 1,nspin
            do m = 1,nbasistot
#ifdef MPI
              call WhichNodeOrb(m,Nodes,BNode)
              if (BNode.eq.0.and.Node.eq.BNode) then
                call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
                ml = m
#endif
                do im = 1,numd(m)
                  read(unit1) dm(im,ml,is)
                enddo
#ifdef MPI
              elseif (Node.eq.0) then
                do im = 1,numd(m)
                  read(unit1) buffer(im)
                enddo
                call MPI_ISend(buffer,numd(m),DAT_double,
     .            BNode,1,MPI_Comm_World,Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)
              elseif (Node.eq.BNode) then
                call GlobalToLocalOrb(m,Node,Nodes,ml)
                call MPI_IRecv(dm(1,ml,is),numd(m),
     .            DAT_double,0,1,MPI_Comm_World,
     .            Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)
              endif
              if (BNode.ne.0) then
                call MPI_Barrier(MPI_Comm_World,MPIerror)
              endif
#endif
            enddo
          enddo

          if (Node.eq.0) then
            call io_close(unit1)
          endif
          found = .true.

        else
          found = .false.
        endif

      elseif (task.eq.'write' .or. task.eq.'WRITE') then

        if (Node.eq.0) then
          call io_assign(unit1)
          open( unit1, file=fname,
     .      form='unformatted', status='unknown' )
          rewind(unit1)
          write(unit1) nbasistot, nspin
        endif

        if (Node.eq.0) then
          do m=1,nbasistot
            write(unit1) numd(m)
          enddo
        endif

        if (Node.eq.0) then
          do m=1,nbasistot
            do im=1,numd(m)
              write(unit1) listd(im,m)
            enddo
          enddo
        endif

        do is=1,nspin
          do m=1,nbasistot
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
              ml = m
#endif
              do im=1,numd(m)
                write(unit1) dm(im,ml,is)
              enddo
#ifdef MPI
            elseif (Node.eq.0) then
              call MPI_IRecv(buffer,numd(m),DAT_double,
     .          BNode,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
              call MPI_ISend(dm(1,ml,is),numd(m),DAT_double,
     .          0,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
              if (Node.eq.0) then
                do im=1,numd(m)
                  write(unit1) buffer(im)
                enddo
              endif
            endif
#endif

          enddo
        enddo

        if (Node.eq.0) then
          call io_close(unit1)
        endif

      else
        if (Node.eq.0) then
          stop 'iodm: incorrect task'
        else
          stop
        endif
      endif

#ifdef MPI
c Free buffer array 
      call memory('D','D',size(buffer),'iodm')
      deallocate(buffer)
#endif

*     call timer( 'iodm', 2 )
      end




