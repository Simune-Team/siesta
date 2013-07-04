   MODULE m_matrixwrite
   implicit none
   PUBLIC

   CONTAINS
      subroutine write_global_matrix_singlenodewrite( &
          no_u, no_s, maxnh, numh, listhptr, listh,   &
          matrix_data, fname)

!
!  Modules
!
      use precision, only: dp
      use parallel,     only : Node, Nodes
      use parallelsubs, only : WhichNodeOrb, LocalToGlobalOrb, &
                               GlobalToLocalOrb, GetNodeOrbs
      use files,        only : slabel, label_length
      use sys,          only : die
#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

      integer, intent(in)                :: maxnh, no_u, no_s
      integer, intent(in), dimension(:)  :: listh, numh, listhptr
      real(dp), intent(in), dimension(:) :: matrix_data
      character(*), intent(in)           :: fname

      external          io_assign, io_close

! Internal variables and arrays
      integer    im, iu
      integer    ih,hl,maxnhtot,maxhg
      integer, dimension(:), allocatable :: numhg, hg_ptr

#ifdef MPI
      integer    MPIerror, Request, Status(MPI_Status_Size), BNode
      integer,  dimension(:),   allocatable :: ibuffer
      real(dp), dimension(:),   allocatable :: buffer
#endif


! Find total numbers over all Nodes
#ifdef MPI
      call MPI_AllReduce(maxnh,maxnhtot,1,MPI_integer,MPI_sum,MPI_Comm_World,MPIerror)
#else
      maxnhtot = maxnh
#endif

        if (Node.eq.0) then
! Open file
          call io_assign( iu )
!          open( iu, file=fname, form='unformatted', status='unknown' )      
          open( unit=iu, file=fname, status='unknown' )

! Write overall data
          write(iu,*) no_s, no_s, maxnhtot

! Allocate local array for global numh
          allocate(numhg(no_u))
          call memory('A','I',no_u,'iohs')

        endif

! Create globalised numh
        do ih = 1,no_u
#ifdef MPI
          call WhichNodeOrb(ih,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
            hl = ih
#endif
            numhg(ih) = numh(hl)
#ifdef MPI
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            call MPI_ISend(numh(hl),1,MPI_integer,0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.0) then
            call MPI_IRecv(numhg(ih),1,MPI_integer,BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
          endif
#endif
        enddo

        if (Node.eq.0) then
! Write row pointers
          allocate(hg_ptr(no_u+1))
          maxhg = 0
          hg_ptr(1) = 1
          do ih = 1,no_u
            maxhg = max(maxhg,numhg(ih))
            hg_ptr(ih+1) = hg_ptr(ih) + numhg(ih)
          enddo
          write(iu,*) (hg_ptr(ih),ih=1,no_u+1)
          deallocate(hg_ptr)

#ifdef MPI
          allocate(buffer(maxhg))
          call memory('A','D',maxhg,'iohs')
          allocate(ibuffer(maxhg))
          call memory('A','I',maxhg,'iohs')
#endif
        endif

! Write listh
        do ih = 1,no_u
#ifdef MPI
          call WhichNodeOrb(ih,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
            hl = ih
#endif
            write(iu,*) (listh(listhptr(hl)+im),im = 1,numh(hl))
#ifdef MPI
          elseif (Node.eq.0) then
            call MPI_IRecv(ibuffer,numhg(ih),MPI_integer,BNode,1, &
              MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            call MPI_ISend(listh(listhptr(hl)+1),numh(hl),MPI_integer, &
              0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
            if (Node.eq.0) then
               write(iu,*) (ibuffer(im),im = 1,numhg(ih))
            endif
          endif
#endif
        enddo

#ifdef MPI
        if (Node.eq.0) then
          call memory('D','I',size(ibuffer),'iohs')
          deallocate(ibuffer)
        endif
#endif

! Write Hamiltonian
!        do is=1,nspin
          do ih=1,no_u
#ifdef MPI
            call WhichNodeOrb(ih,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
              hl = ih
#endif
              write(iu,*) (matrix_data(listhptr(hl)+im), im=1,numh(hl))
#ifdef MPI
            elseif (Node.eq.0) then
              call MPI_IRecv(buffer,numhg(ih),MPI_double_precision, &
                BNode,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(ih,Node,Nodes,hl)
              call MPI_ISend(matrix_data(listhptr(hl)+1),numh(hl), &
                MPI_double_precision,0,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
              if (Node.eq.0) then
                 write(iu,*) (buffer(im),im=1,numhg(ih))
              endif
            endif
#endif
          enddo
!        enddo


#ifdef MPI
          if (Node .eq. 0) then
! Free buffer array
             call memory('D','D',size(buffer),'iohs')
             deallocate(buffer)
          endif
#endif


        if (Node.eq.0) then


! Deallocate local array for global numh
          call memory('D','I',size(numhg),'iohs')
          deallocate(numhg)
! Close file
          call io_close( iu )
        endif

      end subroutine write_global_matrix_singlenodewrite

end module m_matrixwrite
 
