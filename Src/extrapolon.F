      subroutine extrapolon(istep,iord,nspin,nrow,nmrow,nmrowl,nmax,
     .                      num,list,numold,listold,mm2,mnew)
C ******************************************************************************
C Subroutine to extrapolate a given matrix M (like the coefficients of the
C wave functions, or the density matrix) for the next MD step.
C The matrix M is given in sparse form.
C Order N version.
C
C Writen by P.Ordejon, November'96.
C ******************************* INPUT ***************************************
C integer istep                : Time step of the simulation
C integer iord                 : Extrapolation order (0 or 1)
C                                0 = 0th order;  1 = 1st order
C integer nspin                : Number of spin polarizations (1 or 2)
C integer nrow                 : Number of rows of matrix M
C integer nmrow                : Maximum number of rows of matrix M (dimension)
C integer nmrowl               : Maximum local number of rows of matrix M 
C integer nmax                 : First dimension of M matrix, and maximum
C                                number of nonzero elements of each column of M
C integer num(nmax)            : Control vector 1 of M matrix at t
C integer list(nmax,nmrow)     : Control vector 2 of M matrix at t
C ************************** INPUT AND OUTPUT *********************************
C integer numold(nmax)         : Input: Control vector 1 of M matrix at t-dt
C                                       (if istep .ne. 1)
C                                Output: Control vector 1 of M matrix at t
C integer listold(nmax,nmrowl) : Input: Control vector 2 of M matrix at t-dt
C                                       (if istep .ne. 1)
C                                Output: Control vector 2 of M matrix at t
C real*8 mm2(nmax,nmrowl,nspin) : Input: matrix M at t-2dt
C                                Output: matrix M at t-dt
C real*8 mnew(nmax,nmrowl,nspin): New matrix M (extrapolated)
C                                Input: matrix at t-dt
C                                Output: matrix at t
C                                If istep = 1, mnew returned uncahanged
C **************************** BEHAVIOUR **************************************
C The routine allows for the sparse structure of the matrix M to change
C between MD time steps. On input, the matrices of former steps (mnew and mm2) 
C have the structure of last step (t-dt): numold and listold; whereas the new
C (extrapolated) matrix has the structure of the current time step (which
C must be determined before calling this routine!!): num and list.
C On output, the routine updates the structure of mnew and mm2, to that
C at the current (t) time steps respectively. Same with numold and listold
C 
C For the first MD time step (istep = 1), there is no extrapolation. 
C In that case, mnew is returned unchanged.
C Also, in that case numold and listold are only an output, and are set equal
C to num and list
C *****************************************************************************
C
C  Modules
C
      use precision
      use parallel
#ifdef MPI
      use mpi
#endif

      implicit none

#ifdef MPI
      integer
     .  MPIerror
      logical
     .  lbuffer
#endif

      integer 
     .  iord,istep,nmax,nmrow,nmrowl,nrow,nspin

      integer 
     .  list(nmax,nmrow),listold(nmax,nmrowl),num(nmrow),numold(nmrow)

      double precision
     .  mm2(nmax,nmrowl,nspin),mnew(nmax,nmrowl,nspin)
 
C  Internal variables .......................................................

      integer
     .  i,in,ispin,j,il,Node,Nodes
      logical
     .  changed
      double precision
     .  msave
      double precision, dimension(:,:), allocatable, save ::
     .  aux

      external
     .  memory
C ...........................................................................

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif

      if (iord .ne. 0 .and. iord . ne. 1) then
        if (Node.eq.0) then
          write(6,*) 'extrapol: Wrong iord: only 0 and 1 order ',
     .      'available'
        endif
        stop
      endif

C Just initialize numold and listold if istep = 1 ...........................
      if (istep .eq. 1) then
        do i = 1,nmrow
          numold(i) = num(i)
          call GlobalToLocalOrb(i,Node,Nodes,il)
          if (il.gt.0) then
            do in = 1,num(i)
              listold(in,il) = list(in,i)
              do ispin = 1,nspin
                mm2(in,il,ispin) = 0.d0
              enddo
            enddo
          endif
        enddo
        return
C .....................

      else

C Check if sparse structure has changed .....................................
        changed = .false.
        do i = 1,nmrow
          if (numold(i) .ne. num(i)) changed = .true.
        enddo
C Check list arrays if num arrays are unchanged
        if (.not.changed) then
          do il = 1,nmrowl
            call LocalToGlobalOrb(il,Node,Nodes,i)
            do in = 1,num(i)
              if (listold(in,il) .ne. list(in,i)) changed = .true.
            enddo
          enddo
        endif

#ifdef MPI
C Globalise changed flag
        call MPI_AllReduce(changed,lbuffer,1,MPI_logical,MPI_lor,
     .    MPI_Comm_World,MPIerror)
        changed = lbuffer
#endif

        if (changed) then
C .....................

C Allocate local scratch array
          allocate(aux(2,nrow))
          call memory('A','D',2*nrow,'extrapol')

C If sparse structure has changed, re-order mnew and mm2 
C and change numold and listold to current ones .............................

          do i = 1,nrow
            do j = 1,2
              aux(j,i) = 0.0d0
            enddo
          enddo
  
          do i = 1,nmrowl
            call LocalToGlobalOrb(i,Node,Nodes,il)
            do ispin = 1,nspin
              do in = 1,numold(il)
                j = listold(in,i)
                aux(1,j) = mnew(in,i,ispin)
                aux(2,j) = mm2(in,i,ispin)
              enddo
              do in = 1,num(il)
                j = list(in,il)
                mnew(in,i,ispin) = aux(1,j)
                mm2(in,i,ispin) = aux(2,j)
              enddo
              do in = 1,numold(il)
                j = listold(in,i)
                aux(1,j) = 0.0d0
                aux(2,j) = 0.0d0
              enddo
            enddo
            numold(il) = num(il)
            do in = 1,num(il)
              listold(in,i) = list(in,il)
            enddo
          enddo

#ifdef MPI
C  Globalise num/numold
          do i=1,nmrow
            numold(i) = 0
          enddo
          do i=1,nmrowl
            call LocalToGlobalOrb(i,Node,Nodes,il)
            if (il.gt.0) then
              numold(il) = num(il)
            endif
          enddo
          call MPI_AllReduce(numold,num,nmrow,MPI_integer,
     .      MPI_sum,MPI_Comm_World,MPIerror)
          do i=1,nmrow
            numold(i) = num(i)
          enddo
#endif

C Deallocate local scratch array
          call memory('D','D',size(aux),'extrapol')
          deallocate(aux)

        endif
C ..................

C Extrapolate matrix M ......................................................

        do ispin = 1,nspin
          do i = 1,nmrowl
            call LocalToGlobalOrb(i,Node,Nodes,il)
            do in = 1,num(il)
              msave = mnew(in,i,ispin)
              if (iord .eq. 1 .and. mm2(in,i,ispin) .ne. 0.0d0) then
                mnew(in,i,ispin) = 2.0d0 * mnew(in,i,ispin) -
     .                             mm2(in,i,ispin)
              endif
              mm2(in,i,ispin) = msave
            enddo
          enddo
        enddo
C ....................

      endif

      return
      end
