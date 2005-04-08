      subroutine iodm( task, maxnd, nbasis, nspin, numd, 
     .                 listdptr, listd, dm, found )
C *******************************************************************
C Reads/writes density matrix from/to file
C Written by P.Ordejon and J.M.Soler. May 1997.
C ********* INPUT ***************************************************
C character task*(*) : 'read' or 'write'
C integer   maxnd    : First dimension of listd and dm
C integer   nbasis   : Number of atomic orbitals
C integer   nspin    : Number of spins (1 or 2)
C ********* INPUT OR OUTPUT (depending on task) *********************
C integer numd(nbasis)     : Control vector of DM matrix
C                            (number of nonzero elements of each row)
C integer listdptr(nbasis) : Control vector of DM matrix
C                            (pointer to the start of each row)
C integer listd(maxnd)     : Control vector of DM matrix
C                            (list of nonzero elements of each row)
C real*8  dm(maxnd,nspin)  : Density matrix
C ********* OUTPUT *************************************************
C logical found : Has DM been found in disk? (Only when task='read')
C ******************************************************************
      
C
C  Modules
C
      use precision,    only : dp
      use parallel,     only : Node, Nodes
      use parallelsubs, only : LocalToGlobalOrb
      use sys,          only : die
      use fdf,          only : fdf_boolean, fdf_string
      use ionew
#ifdef MPI
      use mpi_siesta
      use parallel_subs, only : WhichNodeOrb, GlobalToLocalOrb
#endif

      implicit  none

      character(len=*), intent(in) :: task
      integer, intent(in) :: maxnd
      integer, intent(in) :: nbasis
      integer, intent(in) :: nspin

      integer, intent(inout) :: numd(nbasis)
      integer, intent(inout) :: listdptr(nbasis)
      integer, intent(inout) :: listd(maxnd)
      real(dp), intent(inout) :: dm(maxnd, nspin)

      logical, intent(out) :: found

C Internal variables
      logical   exist1, exist2, exist3
      integer   im, is, unit1, unit2, m, nb, ndmax, ns
      integer   nbasistot, ml, ndmaxg
      integer, dimension(:), allocatable, save :: numdg
#ifdef MPI
      integer   MPIerror, Request, Status(MPI_Status_Size)
      integer   BNode
      real(dp), dimension(:), allocatable, save :: buffer
      integer,  dimension(:), allocatable, save :: ibuffer
#endif

      

C Saved internal variables:
      logical,           save :: frstme = .true.
      character(len=30), save :: sname
      character(len=33), save :: fnameu
      character(len=34), save :: fnamef
      character(len=34), save :: fnamei, fnameo
      logical,           save :: fmti, fmto
      character(len=11), save :: formin, formout

! Character formats for formatted I/O
! We assume that we have no integers greater than (1e10-1)
! (which is bigger than the largest 32-bit integer)
!
! and that floating point accuracy is sufficiently represented
! by 16 decimal digits of precision (which suffices for 64-bit
! IEEE floats)

      character(len=*), parameter :: intfmt = '(I11)'
      character(len=*), parameter :: floatfmt = '(ES22.14)'

C External procedures
      external          chkdim, memory

C Find file name
      if (frstme) then
        if (node == 0) then
          fmto = fdf_boolean('DM.FormattedFiles', .false.)
          fmti = fdf_boolean('DM.FormattedInput', fmto)
          fmto = fdf_boolean('DM.FormattedOutput', fmto)
          sname = fdf_string('SystemLabel', 'siesta')
          fnameu = trim(sname) // '.DM'
          fnamef = trim(sname) // '.DMF'
          if (fmto) then
            formout = 'formatted'
            fnameo = fnamef
          else
            formout = 'unformatted'
            fnameo = fnameu
          endif
          if (fmti) then
            formin = 'formatted'
            fnamei = fnamef
          else
            formin = 'unformatted'
            fnamei = fnameu
          endif
        endif
      endif

C Find total number of basis functions over all Nodes
#ifdef MPI
      call MPI_AllReduce(nbasis,nbasistot,1,MPI_integer,MPI_sum,
     .     MPI_Comm_World,MPIerror)
#else
      nbasistot = nbasis
#endif

C Allocate local buffer array for globalised numd
      allocate(numdg(nbasistot))
      call memory('A','I',nbasistot,'iodm')

      if (task.eq.'read' .or. task.eq.'READ') then
        if (Node.eq.0) then
          inquire (file='SAVE.DM',     exist=exist1)
          inquire (file='SAVE.ctrlDM', exist=exist2)
          inquire (file=fnamei,        exist=exist3)
        endif

#ifdef MPI
C Broadcast logicals so that all processors take the same route
        call MPI_Bcast(exist1,1,MPI_logical,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(exist2,1,MPI_logical,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(exist3,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

        if (exist1 .and. exist2) then
C Old-format files
          if (Node.eq.0) then
            write(6,'(/,a)') 'iodm: Reading Density Matrix from files'
            call io_assign(unit1)
            call io_assign(unit2)
            open( unit1, file='SAVE.DM', status='unknown')
            open( unit2, file='SAVE.ctrlDM', status='unknown')
            rewind(unit1)
            rewind(unit2)
          endif

          if (Node.eq.0) then
            read(unit2,*) (numdg(m),m=1,nbasistot)
          endif
#ifdef MPI
          call MPI_Bcast(numdg,nbasistot,MPI_integer,0,MPI_Comm_World,
     .      MPIerror)
#endif

C Convert global numd pointer to local form and generate listdptr
          ndmax = 0
          do m = 1,nbasis
            call LocalToGlobalOrb(m,Node,Nodes,ml)
            numd(m) = numdg(ml)
            ndmax = ndmax + numd(m)
            if (m .eq. 1) then
              listdptr(1) = 0
            else
              listdptr(m) = listdptr(m-1) + numd(m-1)
            endif
          enddo
          ndmaxg = 0
          do m = 1,nbasistot
            ndmaxg = max(ndmaxg,numdg(m))
          enddo

C Check size of first dimension of dm
          call chkdim( 'iodm', 'maxnd', maxnd, ndmax, 1 )

#ifdef MPI
C Create buffer arrays for transfering density matrix between nodes and lists
          allocate(buffer(ndmaxg))
          call memory('A','D',ndmaxg,'iodm')
          allocate(ibuffer(ndmaxg))
          call memory('A','I',ndmaxg,'iodm')
#endif

          do m = 1,nbasistot
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
              ml = m
#endif
              read(unit2,*) (listd(listdptr(ml)+im),im=1,numd(ml))
#ifdef MPI
            elseif (Node.eq.0) then
              do im = 1,numdg(m)
                read(unit2,*) ibuffer(im)
              enddo
              call MPI_ISend(ibuffer,numdg(m),MPI_integer,
     .          BNode,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
              call MPI_IRecv(listd(listdptr(ml)+1),numd(ml),
     .          MPI_integer,0,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
            endif
#endif
          enddo

#ifdef MPI
          call memory('D','I',size(ibuffer),'iodm')
          deallocate(ibuffer)
#endif

          do is = 1,nspin
            do m = 1,nbasistot
#ifdef MPI
              call WhichNodeOrb(m,Nodes,BNode)
              if (BNode.eq.0.and.Node.eq.BNode) then
                call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
                ml = m
#endif
                do im = 1,numd(ml)
                  read(unit1,*) dm(listdptr(ml)+im,is)
                enddo
#ifdef MPI
              elseif (Node.eq.0) then
                do im = 1,numdg(m)
                  read(unit1,*) buffer(im)
                enddo
#ifdef NODAT
                call MPI_ISend(buffer,numdg(m),MPI_double_precision,
     .            BNode,1,MPI_Comm_World,Request,MPIerror)
#else
                call MPI_ISend(buffer,numdg(m),DAT_double,
     .            BNode,1,MPI_Comm_World,Request,MPIerror)
#endif
                call MPI_Wait(Request,Status,MPIerror)
              elseif (Node.eq.BNode) then
                call GlobalToLocalOrb(m,Node,Nodes,ml)
#ifdef NODAT
                call MPI_IRecv(dm(listdptr(ml)+1,is),numd(ml),
     .            MPI_double_precision,0,1,MPI_Comm_World,Request,
     .            MPIerror)
#else
                call MPI_IRecv(dm(listdptr(ml)+1,is),numd(ml),
     .            DAT_double,0,1,MPI_Comm_World,Request,MPIerror)
#endif
                call MPI_Wait(Request,Status,MPIerror)
              endif
              if (BNode.ne.0) then
                call MPI_Barrier(MPI_Comm_World,MPIerror)
              endif
#endif
            enddo
          enddo

#ifdef MPI
C Free buffer array
          call memory('D','D',size(buffer),'iodm')
          deallocate(buffer)
#endif
          if (Node.eq.0) then
            call io_close(unit1)
            call io_close(unit2)
          endif

          found = .true.

        elseif (exist3) then
C New-format files
          if (Node.eq.0) then
            write(6,'(/,a)') 'iodm: Reading Density Matrix from files'
            call io_assign(unit1)
            open( unit1, file=fnamei, form=formin, status='old' )
            rewind(unit1)
            if (fmti) then
              read(unit1, intfmt) nb, ns
            else
              read(unit1) nb, ns
            endif
          endif

C Communicate the values to all Nodes and adjust to allow for
C distributed memory before checking the dimensions
#ifdef MPI
          call MPI_Bcast(nb,1,MPI_integer,0,MPI_Comm_World,MPIerror)
          call MPI_Bcast(ns,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

C Check dimensions
          call chkdim( 'iodm', 'nbasis', nbasistot, nb, 0 )
          call chkdim( 'iodm', 'nspin',  nspin,  ns, 0 )

          if (Node.eq.0) then
            if (fmti) then
              read(unit1, intfmt) (numdg(m),m=1,nbasistot)
            else
              read(unit1) (numdg(m),m=1,nbasistot)
            endif
          endif
#ifdef MPI
          call MPI_Bcast(numdg,nbasistot,MPI_integer,0,MPI_Comm_World,
     .      MPIerror)
#endif

C Convert global numd pointer to local form and generate listdptr
          ndmax = 0
          do m = 1,nbasis
            call LocalToGlobalOrb(m,Node,Nodes,ml)
            numd(m) = numdg(ml)
            ndmax = ndmax + numd(m)
            if (m .eq. 1) then
              listdptr(1) = 0
            else
              listdptr(m) = listdptr(m-1) + numd(m-1)
            endif
          enddo
          ndmaxg = 0
          do m = 1,nbasistot
            ndmaxg = max(ndmaxg,numdg(m))
          enddo

C Check size of first dimension of dm
          call chkdim( 'iodm', 'maxnd', maxnd, ndmax, 1 )

#ifdef MPI
C Create buffer arrays for transfering density matrix between nodes and lists
          allocate(buffer(ndmaxg))
          call memory('A','D',ndmaxg,'iodm')
          allocate(ibuffer(ndmaxg))
          call memory('A','I',ndmaxg,'iodm')
#endif

          do m = 1,nbasistot
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
              ml = m
#endif
              if (fmti) then
                read(unit1, intfmt) 
     .               (listd(listdptr(ml)+im),im=1,numd(ml))
              else
                read(unit1) (listd(listdptr(ml)+im),im=1,numd(ml))
              endif
#ifdef MPI
            elseif (Node.eq.0) then
              if (fmti) then
                read(unit1, intfmt) (ibuffer(im),im=1,numdg(m))
              else
                read(unit1) (ibuffer(im),im=1,numdg(m))
              endif
              call MPI_ISend(ibuffer,numdg(m),MPI_integer,
     .          BNode,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
              call MPI_IRecv(listd(listdptr(ml)+1),numd(ml),
     .          MPI_integer,0,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
            endif
#endif
          enddo

#ifdef MPI
          call memory('D','I',size(ibuffer),'iodm')
          deallocate(ibuffer)
#endif

          do is = 1,nspin
            do m = 1,nbasistot
#ifdef MPI
              call WhichNodeOrb(m,Nodes,BNode)
              if (BNode.eq.0.and.Node.eq.BNode) then
                call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
                ml = m
#endif
                if (fmti) then
                  read(unit1, floatfmt)
     .                 (dm(listdptr(ml)+im,is),im=1,numd(ml))
                else
                  read(unit1) (dm(listdptr(ml)+im,is),im=1,numd(ml))
                endif
#ifdef MPI
              elseif (Node.eq.0) then
                if (fmti) then
                  read(unit1, floatfmt) (buffer(im),im=1,numdg(m))
                else
                  read(unit1) (buffer(im),im=1,numdg(m))
                endif
#ifdef NODAT
                call MPI_ISend(buffer,numdg(m),MPI_double_precision,
     .            BNode,1,MPI_Comm_World,Request,MPIerror)
#else
                call MPI_ISend(buffer,numdg(m),DAT_double,
     .            BNode,1,MPI_Comm_World,Request,MPIerror)
#endif
                call MPI_Wait(Request,Status,MPIerror)
              elseif (Node.eq.BNode) then
                call GlobalToLocalOrb(m,Node,Nodes,ml)
#ifdef NODAT
                call MPI_IRecv(dm(listdptr(ml)+1,is),numd(ml),
     .            MPI_double_precision,0,1,MPI_Comm_World,Request,
     .            MPIerror)
#else
                call MPI_IRecv(dm(listdptr(ml)+1,is),numd(ml),
     .            DAT_double,0,1,MPI_Comm_World,Request,MPIerror)
#endif
                call MPI_Wait(Request,Status,MPIerror)
              endif
              if (BNode.ne.0) then
                call MPI_Barrier(MPI_Comm_World,MPIerror)
              endif
#endif
            enddo
          enddo

#ifdef MPI
C Free buffer array
          call memory('D','D',size(buffer),'iodm')
          deallocate(buffer)
#endif
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
          open( unit1, file=fnameo, form=formout, status='unknown' )
          rewind(unit1)
          if (fmto) then
            write(unit1, intfmt) nbasistot, nspin
          else
            write(unit1) nbasistot, nspin
          endif
        endif

C Create globalised numd
        do m = 1,nbasistot
#ifdef MPI
          call WhichNodeOrb(m,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
            ml = m
#endif
            numdg(m) = numd(ml)
#ifdef MPI
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
            call MPI_ISend(numd(ml),1,MPI_integer,
     .        0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.0) then
            call MPI_IRecv(numdg(m),1,MPI_integer,
     .        BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
          endif
#endif
        enddo

C Write out numd array
        if (Node.eq.0) then
          ndmaxg = 0
          do m = 1,nbasistot
            ndmaxg = max(ndmaxg,numdg(m))
          enddo
          if (fmto) then
            write(unit1, intfmt) (numdg(m),m=1,nbasistot)
          else
            write(unit1) (numdg(m),m=1,nbasistot)
          endif
#ifdef MPI
          allocate(buffer(ndmaxg))
          call memory('A','D',ndmaxg,'iodm')
          allocate(ibuffer(ndmaxg))
          call memory('A','I',ndmaxg,'iodm')
#endif
        endif

C Write out listd array
        do m = 1,nbasistot
#ifdef MPI
          call WhichNodeOrb(m,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
            ml = m
#endif
            if (fmto) then
              write(unit1, intfmt)
     .             (listd(listdptr(ml)+im),im=1,numd(ml))
            else
              write(unit1) (listd(listdptr(ml)+im),im=1,numd(ml))
            endif
#ifdef MPI
          elseif (Node.eq.0) then
            call MPI_IRecv(ibuffer,numdg(m),MPI_integer,BNode,1,
     .        MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
            call MPI_ISend(listd(listdptr(ml)+1),numd(ml),MPI_integer,
     .        0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
            if (Node.eq.0) then
              if (fmto) then
                write(unit1, intfmt) (ibuffer(im),im=1,numdg(m))
              else
                write(unit1) (ibuffer(im),im=1,numdg(m))
              endif
            endif
          endif
#endif
        enddo

#ifdef MPI
        if (Node.eq.0) then
          call memory('D','I',size(ibuffer),'iodm')
          deallocate(ibuffer)
        endif
#endif

C Write density matrix
        do is=1,nspin
          do m=1,nbasistot
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
              ml = m
#endif
              if (fmto) then
                write(unit1, floatfmt) 
     .               (dm(listdptr(ml)+im,is),im=1,numd(ml))
              else
                write(unit1) (dm(listdptr(ml)+im,is),im=1,numd(ml))
              endif
#ifdef MPI
            elseif (Node.eq.0) then
#ifdef NODAT
              call MPI_IRecv(buffer,numdg(m),MPI_double_precision,
     .          BNode,1,MPI_Comm_World,Request,MPIerror)
#else
              call MPI_IRecv(buffer,numdg(m),DAT_double,
     .          BNode,1,MPI_Comm_World,Request,MPIerror)
#endif
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
#ifdef NODAT
              call MPI_ISend(dm(listdptr(ml)+1,is),numd(ml),
     .          MPI_double_precision,0,1,MPI_Comm_World,Request,
     .          MPIerror)
#else
              call MPI_ISend(dm(listdptr(ml)+1,is),numd(ml),DAT_double,
     .          0,1,MPI_Comm_World,Request,MPIerror)
#endif
              call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
              if (Node.eq.0) then
                if (fmto) then
                  write(unit1, floatfmt) (buffer(im),im=1,numdg(m))
                else
                  write(unit1) (buffer(im),im=1,numdg(m))
                endif
              endif
            endif
#endif
          enddo
        enddo

        if (Node.eq.0) then
#ifdef MPI
          call memory('D','D',size(buffer),'iodm')
          deallocate(buffer)
#endif
          call io_close(unit1)
        endif

      else
        if (Node.eq.0) then
          call die('iodm: incorrect task')
        endif
      endif

C Deallocate local buffer array for globalised numd
      call memory('D','I',size(numdg),'iodm')
      deallocate(numdg)

      end subroutine iodm
