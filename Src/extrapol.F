      subroutine extrapol(istep,iord,nspin,nrow,nmrow,nmrowl,nmax,num,
     .                    listptr,list,numold,listptrold,listold,
     .                    mm2,mnew)
C ******************************************************************************
C Subroutine to extrapolate a given matrix M (like the coefficients of the
C wave functions, or the density matrix) for the next MD step.
C The matrix M is given in sparse form.
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
C                                number of nonzero elements of M
C integer num(nmrowl)          : Control vector 1 of M matrix at t
C integer listptr(nmrowl)      : Control vector 2 of M matrix at t
C integer list(nmax)           : Control vector 3 of M matrix at t
C ************************** INPUT AND OUTPUT *********************************
C integer numold(nmrowl)       : Input: Control vector 1 of M matrix at t-dt
C                                       (if istep .ne. 1)
C                                Output: Control vector 1 of M matrix at t
C integer listptrold(nmrowl)   : Input: Control vector 1 of M matrix at t-dt
C                                       (if istep .ne. 1)
C                                Output: Control vector 1 of M matrix at t
C integer listold(nmax)        : Input: Control vector 3 of M matrix at t-dt
C                                       (if istep .ne. 1)
C                                Output: Control vector 3 of M matrix at t
C real*8 mm2(nmax,nspin)       : Input: matrix M at t-2dt
C                                Output: matrix M at t-dt
C real*8 mnew(nmax,nspin)      : New matrix M (extrapolated)
C                                Input: matrix at t-dt
C                                Output: matrix at t
C                                If istep = 1, mnew returned unchanged
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
     .  list(nmax),listold(nmax),num(nmrowl),numold(nmrowl),
     .  listptr(nmrowl),listptrold(nmrowl)

      double precision
     .  mm2(nmax,nspin),mnew(nmax,nspin)
 
C  Internal variables .......................................................

      integer
     .  i,in,ind,ispin,j,Node,nh
      logical
     .  changed
      double precision
     .  msave
      double precision, dimension(:), allocatable, save ::
     .  aux, mtmp

      external
     .  memory
C ...........................................................................

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

C Find the maximum number of non-zero elements in matrices
      nh = listptr(nmrowl) + num(nmrowl)

      if (iord .ne. 0 .and. iord . ne. 1) then
        if (Node.eq.0) then
          write(6,*) 'extrapol: Wrong iord: only 0 and 1 order ',
     .      'available'
        endif
        stop
      endif

C Just initialize numold and listold if istep = 1 ...........................
      if (istep .eq. 1) then
        numold(1:nmrowl) = num(1:nmrowl)
        listptrold(1:nmrowl) = listptr(1:nmrowl)
        listold(1:nh) = list(1:nh)
        mm2(1:nmax,1:nspin) = 0.0d0
        return
C .....................

      else

C Check if sparse structure has changed .....................................
        changed = .false.
        do i = 1,nmrowl
          if (numold(i) .ne. num(i)) changed = .true.
        enddo
C Check list arrays if num arrays are unchanged
        if (.not.changed) then
          do i = 1,nh
            if (listold(i) .ne. list(i)) changed = .true.
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

C If sparse structure has changed, re-order mnew and mm2 
C and change numold and listold to current ones .............................
C It is now necessary to use a buffer array during changes to mnew and mm2
C because boundaries between rows may change in either direction

C Allocate local scratch arrays
          allocate(aux(nrow))
          call memory('A','D',nrow,'extrapol')
          allocate(mtmp(nh))
          call memory('A','D',nh,'extrapol')

          aux(1:nrow) = 0.0d0
  
          do ispin = 1,nspin
            do i = 1,nmrowl
              do in = 1,numold(i)
                ind = listptrold(i)+in
                j = listold(ind)
                aux(j) = mnew(ind,ispin)
              enddo
              do in = 1,num(i)
                ind = listptr(i)+in
                j = list(ind)
                mtmp(ind) = aux(j)
              enddo
              do in = 1,numold(i)
                j = listold(listptrold(i)+in)
                aux(j) = 0.0d0
              enddo
            enddo
            mnew(1:nh,ispin) = mtmp(1:nh)
          enddo
          do ispin = 1,nspin
            do i = 1,nmrowl
              do in = 1,numold(i)
                ind = listptrold(i)+in
                j = listold(ind)
                aux(j) = mm2(ind,ispin)
              enddo
              do in = 1,num(i)
                ind = listptr(i)+in
                j = list(ind)
                mtmp(ind) = aux(j)
              enddo
              do in = 1,numold(i)
                j = listold(listptrold(i)+in)
                aux(j) = 0.0d0
              enddo
            enddo
            mm2(1:nh,ispin) = mtmp(1:nh)
          enddo

          numold(1:nmrowl) = num(1:nmrowl)
          listptrold(1:nmrowl) = listptr(1:nmrowl)
          listold(1:nh) = list(1:nh)

C Deallocate local scratch arrays
          call memory('D','D',size(mtmp),'extrapol')
          deallocate(mtmp)
          call memory('D','D',size(aux),'extrapol')
          deallocate(aux)

        endif
C ..................

C Extrapolate matrix M ......................................................

        do ispin = 1,nspin
          do i = 1,nh
            msave = mnew(i,ispin)
            if (iord .eq. 1 .and. mm2(i,ispin) .ne. 0.0d0) then
              mnew(i,ispin) = 2.0d0 * mnew(i,ispin) - mm2(i,ispin)
            endif
            mm2(i,ispin) = msave
          enddo
        enddo
C ....................

      endif

      return
      end
