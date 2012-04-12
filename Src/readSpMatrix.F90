
 subroutine readSpMatrix (filename, SpM, found, ref_dist )

   use class_SpMatrix
   use class_Sparsity
   use class_Array2D
   use class_OrbitalDistribution

#ifdef MPI
   use mpi_siesta
#endif

   character(len=*), intent(in) :: filename

   ! Note: inout is essential to avoid memory leaks
   type(SpMatrix), intent(inout)  :: SpM
   logical, intent(out)         :: found
   type(OrbitalDistribution), intent(in) :: ref_dist

      logical   exist3
      integer   im, is, lun, m, nb
      integer   ml
      integer :: Node, Nodes, Node_io, Comm

      integer, parameter :: dp = selected_real_kind(10,100)

      integer, allocatable, dimension(:) :: numdg, numd, listdptr, listd
      real(dp), allocatable              :: dm(:,:)

      type(Sparsity)    :: sp_read
      type(Array2D)     :: a2d_read

#ifdef MPI
      integer   MPIerror, Request, Status(MPI_Status_Size)
      integer   BNode, ndmaxg

      real(dp), dimension(:), allocatable :: buffer
      integer,  dimension(:), allocatable :: ibuffer
#endif

      Node = ref_dist%data%node
      Nodes = ref_dist%data%nodes
      Node_io = ref_dist%data%node_io
      Comm = ref_dist%data%comm

!     Find file name
      if (Node.eq. Node_io) then
         inquire (file=filename, exist=exist3)
      endif
#ifdef MPI
      call MPI_Bcast(exist3,1,MPI_logical,Node_io,Comm,MPIerror)
#endif

      found = .false.
      if ( .not. exist3) RETURN

      if (Node.eq. Node_io) then
         write(6,'(/,a)') 'Reading SpMatrixfrom file '// trim(filename)
         lun = 88
         open( lun, file=filename, form="unformatted", status='old' )
         rewind(lun)
         read(lun) no_u, nspin
      endif

      
!     Communicate the values to all Nodes and adjust to allow for
!     distributed memory before checking the dimensions
#ifdef MPI
      call MPI_Bcast(no_u,1,MPI_integer,Node_io,Comm,MPIerror)
      call MPI_Bcast(nspin,1,MPI_integer,Node_io,Comm,MPIerror)
#endif

!     Allocate local buffer array for globalised numd
      allocate(numdg(1:no_u))
      if (Node.eq.Node_io) then
         read(lun) (numdg(m),m=1,no_u)
      endif
#ifdef MPI
      call MPI_Bcast(numdg,no_u,MPI_integer,Node_io,Comm,MPIerror)
#endif

      no_l = num_local_elements(ref_dist,no_u,Node)
!     Convert global numd pointer to local form and generate listdptr
      allocate(numd(1:no_l),listdptr(1:no_l))
      maxnd = 0
      do m = 1,no_l
#ifdef MPI
         mg = index_local_to_global(ref_dist,m,Node)
#else
         mg = m
#endif
         numd(m) = numdg(mg)
         maxnd = maxnd + numdg(mg)
         if (m .eq. 1) then
            listdptr(1) = 0
         else
            listdptr(m) = listdptr(m-1) + numd(m-1)
         endif
      enddo

      allocate(listd(1:maxnd))
      allocate(dm(1:maxnd,1:nspin))

#ifdef MPI
!  Create buffer arrays for transfering density matrix between nodes and lists
      ndmaxg = maxval(numdg(1:no_u))
      allocate(buffer(1:ndmaxg))
      allocate(ibuffer(1:ndmaxg))
#endif

      do m = 1,no_u
#ifdef MPI
         Bnode = node_handling_element(ref_dist,m)
         if (BNode==Node_io .and. Node==BNode) then
            ml = index_global_to_local(ref_dist,m,Node)
#else
            ml = m
#endif
            read(lun) (listd(listdptr(ml)+im),im=1,numd(ml))
#ifdef MPI
         elseif (Node == Node_io) then
            read(lun) (ibuffer(im),im=1,numdg(m))
            call MPI_ISend(ibuffer,numdg(m),MPI_integer,  &
                           BNode,1,Comm,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)

         elseif (Node == BNode) then
            ml = index_global_to_local(ref_dist,m,Node)
            call MPI_IRecv(listd(listdptr(ml)+1),numd(ml),   &
                MPI_integer,Node_io,1,Comm,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
         endif
         if (BNode.ne.Node_io) then
            call MPI_Barrier(Comm,MPIerror)
         endif
#endif
      enddo

#ifdef MPI
      deallocate(ibuffer)
#endif

      do is = 1,nspin
         do m = 1,no_u
#ifdef MPI
            Bnode = node_handling_element(ref_dist,m)
            if (BNode==Node_io .and. Node==BNode) then
               ml = index_global_to_local(ref_dist,m,Node)
#else
               ml = m
#endif
               read(lun) (dm(listdptr(ml)+im,is),im=1,numd(ml))
#ifdef MPI
            elseif (Node == Node_io) then
               read(lun) (buffer(im),im=1,numdg(m))
               call MPI_ISend(buffer,numdg(m),MPI_double_precision,  &
                   BNode,1,Comm,Request,MPIerror)
               call MPI_Wait(Request,Status,MPIerror)

            elseif (Node == BNode) then
               ml = index_global_to_local(ref_dist,m,Node)
               call MPI_IRecv(dm(listdptr(ml)+1,is),numd(ml),   &
                    MPI_double_precision,Node_io,1,Comm,Request, &
                    MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.Node_io) then
               call MPI_Barrier(Comm,MPIerror)
            endif
#endif
         enddo
      enddo

#ifdef MPI
      deallocate(buffer)
#endif
      deallocate(numdg)

      if (Node == Node_io) then
         close(lun)
      endif

      found = .true.

      call newSparsity(sp_read,no_l,no_u,  &
                       maxnd,numd,listdptr,listd,  &
                       "(read from " // trim(filename) // ")")
      call newArray2D(a2d_read,dm,name="(Array from read dm)")
      call newSpMatrix(sp_read,a2d_read,ref_dist,SpM, &
                       "(read from " // trim(filename) // ")")

      call delete(sp_read)
      call delete(a2d_read)
      deallocate(numd,listdptr,listd,dm)

    end subroutine readSpMatrix
