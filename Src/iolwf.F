      subroutine iolwf( task, maxc, ncmax, maxoloc, nbasis, nspin,
     .                  numc, listc, c, found )

C *******************************************************************
C Reads/writes localized wave functions from/to file
C Written by P.Ordejon and J.M.Soler. May 1997.
C ********* INPUT ***************************************************
C character task*3 : 'read' or 'write'
C integer   maxc   : First dimension of listd and c
C integer   maxoloc: Second dimension of c/listc 
C                  : (must be maxoloc.ge.nbasisloc)
C integer   nbasis : Number of atomic orbitals
C integer   nspin  : Number of spins (1 or 2)
C ********* INPUT OR OUTPUT (depending on task) *********************
C integer numc(nbasis)       : Control vector of c matrix
C                              (number of nonzero elements of each row)
C integer listc(maxc,maxoloc): Control vector of c matrix
C                              (list of nonzero elements of each row)
C real*8  c(maxc,maxoloc,nspin) : Density matrix
C ********* OUTPUT *************************************************
C integer ncmax : actual value of maxc needed
C logical found : Have LWF's been found in disk? (Only when task='read')
C ******************************************************************

C
C  Modules
C
      use precision
      use parallel
      use fdf
#ifdef MPI
      use mpi_siesta
#endif

      implicit  none

      character task*(*), paste*33
      logical   found
      integer   maxc, maxoloc, nbasis, nspin, ncmax
      integer   listc(maxc,maxoloc), numc(nbasis)
      real*8    c(maxc,maxoloc,nspin)
      external  chkdim, io_assign, io_close, paste, timer,
     .          memory

C Internal variables and arrays
      character fname*33, sname*30
      logical   exist1, exist2, exist3, frstme
      integer   im, is, unit1, unit2, m, nb, ns, Node, Nodes,
     .          ml, BNode
#ifdef MPI
      integer   MPIerror, Request, Status(MPI_Status_Size)
      real*8,  dimension(:), allocatable, save :: bdens
      integer, dimension(:), allocatable, save :: idens
#endif

      save      frstme, fname
      data      frstme /.true./

*     call timer( 'iolwf', 1 )

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif

#ifdef MPI
C Allocate local buffer array
      allocate(bdens(maxc))
      call memory('A','D',maxc,'iolwf')
      allocate(idens(maxc))
      call memory('A','I',maxc,'iolwf')
#endif
      
C Find file name
      if (frstme) then
        if (Node.eq.0) then
          sname = fdf_string('SystemLabel','siesta')
          fname = paste(sname,'.LWF')
        endif
        frstme = .false.
      endif

      if (task.eq.'read' .or. task.eq.'READ') then
        if (Node.eq.0) then
          inquire (file='SAVE.LWF',     exist=exist1)
          inquire (file='SAVE.ctrlLWF', exist=exist2)
          inquire (file=fname,          exist=exist3)
        endif

#ifdef MPI
C Broadcast logicals
        call MPI_Bcast(exist1,1,MPI_logical,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(exist2,1,MPI_logical,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(exist3,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

C First look for old-format files
        if (exist1 .and. exist2) then
          if (Node.eq.0) then
            write(6,'(/,a)') 'iolwf: Reading LWFs from file'

            call io_assign(unit1)
            call io_assign(unit2)
            open( unit1, file='SAVE.LWF', status='unknown')
            open( unit2, file='SAVE.ctrlLWF', status='unknown')
            rewind(unit1)
            rewind(unit2)
          endif
          ncmax = 0
          if (Node.eq.0) then
            do m = 1,nbasis
              read(unit2,*) numc(m)
              ncmax = max( ncmax, numc(m) )
            enddo
          endif
#ifdef MPI
          call MPI_Bcast(numc,nbasis,MPI_integer,0,MPI_Comm_World,
     .      MPIerror)
          call MPI_Bcast(ncmax,1,MPI_integer,0,MPI_Comm_World,
     .      MPIerror)
#endif
C Check dimensions are OK
          if (ncmax.le.maxc) then
            do m = 1,nbasis
              call WhichNodeOrb(m,Nodes,BNode)

#ifdef MPI
              if (Node.eq.0.and.BNode.eq.Node) then
C Local read
#endif
                call GlobalToLocalOrb(m,Node,Nodes,ml)
                do im = 1,numc(m)
                  read(unit2,*) listc(im,ml)
                enddo

#ifdef MPI
              elseif (Node.eq.0) then
C Read then send to node where density matrix row is stored
                do im = 1,numc(m)
                  read(unit2,*) idens(im)
                enddo
                call MPI_ISend(idens,numc(m),MPI_integer,
     .            BNode,1,MPI_Comm_World,Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)

              elseif (Node.eq.BNode) then
C Receive density matrix row
                call GlobalToLocalOrb(m,Node,Nodes,ml)
                call MPI_IRecv(listc(1,ml),numc(m),MPI_integer,
     .            0,1,MPI_Comm_World,Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)

              endif

              if (BNode.ne.0) then
                call MPI_Barrier(MPI_Comm_World,MPIerror)
              endif
#endif
            enddo

            do is = 1,nspin
              do m = 1,nbasis
                call WhichNodeOrb(m,Nodes,BNode)

#ifdef MPI
                if (Node.eq.0.and.BNode.eq.Node) then
C Local read
#endif
                  call GlobalToLocalOrb(m,Node,Nodes,ml)
                  do im = 1,numc(m)
                    read(unit1,*) c(im,ml,is)
                  enddo

#ifdef MPI
                elseif (Node.eq.0) then
C Read then send to node where density matrix row is stored
                  do im = 1,numc(m)
                    read(unit1,*) bdens(im)
                  enddo
#ifdef NODAT
                  call MPI_ISend(bdens,numc(m),MPI_double_precision,
     .              BNode,1,MPI_Comm_World,Request,MPIerror)
#else
                  call MPI_ISend(bdens,numc(m),DAT_double,
     .              BNode,1,MPI_Comm_World,Request,MPIerror)
#endif
                  call MPI_Wait(Request,Status,MPIerror)

                elseif (Node.eq.BNode) then
C Receive density matrix row
                  call GlobalToLocalOrb(m,Node,Nodes,ml)
#ifdef NODAT
                  call MPI_IRecv(c(1,ml,is),numc(m),
     .              MPI_double_precision,0,1,MPI_Comm_World,
     .              Request,MPIerror)
#else
                  call MPI_IRecv(c(1,ml,is),numc(m),DAT_double,
     .              0,1,MPI_Comm_World,Request,MPIerror)
#endif
                  call MPI_Wait(Request,Status,MPIerror)

                endif

                if (BNode.ne.0) then
                  call MPI_Barrier(MPI_Comm_World,MPIerror)
                endif
#endif

              enddo
            enddo
          endif
          if (Node.eq.0) then
            call io_close(unit1)
            call io_close(unit2)
          endif
          found = .true.

c       Look now for new-format files
        elseif (exist3) then
          if (Node.eq.0) then
            write(6,'(/,a)') 'iolwf: Reading LWFs from files'
            call io_assign(unit1)
            open( unit1, file=fname,
     .          form='unformatted', status='unknown' )
            rewind(unit1)
            read(unit1) nb, ns
          endif
#ifdef MPI
          call MPI_Bcast(nb,1,MPI_integer,0,MPI_Comm_World,MPIerror)
          call MPI_Bcast(ns,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
          call chkdim( 'iolwf', 'nbasis', nbasis, nb, 0 )
          call chkdim( 'iolwf', 'nspin',  nspin,  ns, 0 )
          if (Node.eq.0) then
            read(unit1) (numc(m), m=1,nbasis)
          endif
#ifdef MPI
          call MPI_Bcast(numc,nbasis,MPI_integer,0,MPI_Comm_World,
     .      MPIerror)
#endif
          ncmax = 0
          do m = 1,nbasis
            ncmax = max( ncmax, numc(m) )
          enddo
          if (ncmax.le.maxc) then
            do m = 1,nbasis
              call WhichNodeOrb(m,Nodes,BNode)

#ifdef MPI
              if (Node.eq.0.and.BNode.eq.Node) then
C Local read
#endif
                call GlobalToLocalOrb(m,Node,Nodes,ml)
                read(unit1) (listc(im,ml), im=1,numc(m))

#ifdef MPI
              elseif (Node.eq.0) then
C Read density matrix row and send to node where it is to be stored
                read(unit1) (idens(im), im=1,numc(m))
                call MPI_ISend(idens,numc(m),MPI_integer,
     .            BNode,1,MPI_Comm_World,Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)

              elseif (Node.eq.BNode) then
C Receive density matrix row
                call GlobalToLocalOrb(m,Node,Nodes,ml)
                call MPI_IRecv(listc(1,ml),numc(m),MPI_integer,
     .            0,1,MPI_Comm_World,Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)

              endif

              if (BNode.ne.0) then
                call MPI_Barrier(MPI_Comm_World,MPIerror)
              endif
#endif
            enddo

            do is = 1,nspin
              do m = 1,nbasis
                call WhichNodeOrb(m,Nodes,BNode)

#ifdef MPI
                if (Node.eq.0.and.BNode.eq.Node) then
C Local read
#endif
                  call GlobalToLocalOrb(m,Node,Nodes,ml)
                  read(unit1) (c(im,ml,is), im=1,numc(m))

#ifdef MPI
                elseif (Node.eq.0) then
C Read density matrix row and send to node where it is to be stored
                  read(unit1) (bdens(im), im=1,numc(m))
#ifdef NODAT
                  call MPI_ISend(bdens,numc(m),MPI_double_precision,
     .              BNode,1,MPI_Comm_World,Request,MPIerror)
#else
                  call MPI_ISend(bdens,numc(m),DAT_double,
     .              BNode,1,MPI_Comm_World,Request,MPIerror)
#endif
                  call MPI_Wait(Request,Status,MPIerror)

                elseif (Node.eq.BNode) then
C Receive density matrix row
                  call GlobalToLocalOrb(m,Node,Nodes,ml)
#ifdef NODAT
                  call MPI_IRecv(c(1,ml,is),numc(m),
     .              MPI_double_precision,0,1,MPI_Comm_World,
     .              Request,MPIerror)
#else
                  call MPI_IRecv(c(1,ml,is),numc(m),DAT_double,
     .              0,1,MPI_Comm_World,Request,MPIerror)
#endif
                  call MPI_Wait(Request,Status,MPIerror)

                endif

                if (BNode.ne.0) then
                  call MPI_Barrier(MPI_Comm_World,MPIerror)
                endif
#endif
              enddo
            enddo

          endif

          if (Node.eq.0) then
            call io_close(unit1)
          endif
          found = .true.

        else

          found = .false.
          ncmax = 0

        endif

      elseif (task.eq.'write' .or. task.eq.'WRITE') then

        if (Node.eq.0) then
          call io_assign(unit1)
          open( unit1, file=fname,
     .        form='unformatted', status='unknown' )
          rewind(unit1)
          write(unit1) nbasis, nspin
          write(unit1) (numc(m), m=1,nbasis)
        endif

        do m = 1,nbasis
          call WhichNodeOrb(m,Nodes,BNode)

#ifdef MPI
          if (Node.eq.0.and.BNode.eq.Node) then
C Local write
#endif
            call GlobalToLocalOrb(m,Node,Nodes,ml)
            write(unit1)(listc(im,ml), im=1,numc(m))

#ifdef MPI
          elseif (Node.eq.0) then
C Receive density matrix row
            call MPI_IRecv(idens,numc(m),MPI_integer,
     .        BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)

          elseif (Node.eq.BNode) then
C Send density matrix row
            call GlobalToLocalOrb(m,Node,Nodes,ml)
            call MPI_ISend(listc(1,ml),numc(m),MPI_integer,
     .        0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)

          endif

          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
            if (Node.eq.0) then
              write(unit1) (idens(im),im=1,numc(m))
            endif
          endif
#endif
        enddo

        do is = 1,nspin
          do m = 1,nbasis
            call WhichNodeOrb(m,Nodes,BNode)

#ifdef MPI
            if (Node.eq.0.and.BNode.eq.Node) then
C Local write
#endif
              call GlobalToLocalOrb(m,Node,Nodes,ml)
              write(unit1)(c(im,ml,is), im=1,numc(m))

#ifdef MPI
            elseif (Node.eq.0) then
C Receive density matrix row
#ifdef NODAT
              call MPI_IRecv(bdens,numc(m),MPI_double_precision,
     .          BNode,1,MPI_Comm_World,Request,MPIerror)
#else
              call MPI_IRecv(bdens,numc(m),DAT_double,
     .          BNode,1,MPI_Comm_World,Request,MPIerror)
#endif
              call MPI_Wait(Request,Status,MPIerror)

            elseif (Node.eq.BNode) then
C Send density matrix row
              call GlobalToLocalOrb(m,Node,Nodes,ml)
#ifdef NODAT
              call MPI_ISend(c(1,ml,is),numc(m),MPI_double_precision,
     .          0,1,MPI_Comm_World,Request,MPIerror)
#else
              call MPI_ISend(c(1,ml,is),numc(m),DAT_double,
     .          0,1,MPI_Comm_World,Request,MPIerror)
#endif
              call MPI_Wait(Request,Status,MPIerror)

            endif

            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
              if (Node.eq.0) then
                write(unit1) (bdens(im),im=1,numc(m))
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
          stop 'iolwf: incorrect task'
        else
          stop 
        endif
      endif

#ifdef MPI
C Free local buffer memory
      call memory('D','D',size(bdens),'iolwf')
      deallocate(bdens)
      call memory('D','I',size(idens),'iolwf')
      deallocate(idens)
#endif

*     call timer( 'iolwf', 2 )
      end

