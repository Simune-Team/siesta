! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine extrapol(istep,iord,nspin,nrow,nmrowl,nmax,num,
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
      use precision,    only : dp
      use parallel,     only : node
      use sys,          only : die
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      integer, intent(in) :: istep
      integer, intent(in) :: iord
      integer, intent(in) :: nspin
      integer, intent(in) :: nrow
      integer, intent(in) :: nmrowl
      integer, intent(in) :: nmax
      integer, intent(in) :: num(nmrowl)
      integer, intent(in) :: listptr(nmrowl)
      integer, intent(in) :: list(nmax)

      integer, intent(inout) :: numold(nmrowl)
      integer, intent(inout) :: listptrold(nmrowl)
      integer, intent(inout) :: listold(nmax)
      real(dp), intent(inout) :: mm2(nmax, nspin)
      real(dp), intent(inout) :: mnew(nmax, nspin)

C  Internal variables .......................................................

      integer
     .  i,in,ind,ispin,j,nh,nrmax
      logical
     .  changed
      real(dp)
     .  msave
      real(dp), dimension(:), allocatable :: aux, mtmp

      external
     .  memory

#ifdef MPI
      integer
     .  MPIerror
      logical
     .  lbuffer
#endif

C ...........................................................................

C Find the maximum number of non-zero elements in matrices
      if (nmrowl>0) then
        nh = listptr(nmrowl) + num(nmrowl)
      else
        nh = 0
        ! in this case, none of the loops will be executed; we'll
        ! fall straight through everything, but still execute the
        ! MPI_AllReduce.
      endif

      if (iord /= 0 .and. iord /= 1) then
        if (node == 0) then
          call die('extrapol: Wrong iord: only 0 and 1 order available')
        endif
      endif

C Just initialize numold and listold if istep = 1 ...........................
      if (istep .eq. 1) then
        numold(1:nmrowl) = num(1:nmrowl)
        listptrold(1:nmrowl) = listptr(1:nmrowl)
        listold(1:nh) = list(1:nh)
        mm2(1:nmax,1:nspin) = 0.0_dp
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
          nrmax = nrow
C Determine size of auxiliary matrix
          do i = 1,nmrowl
            do in = 1,numold(i)
              ind = listptrold(i)+in
              j = listold(ind)
              nrmax = max(nrmax,j)
            enddo
            do in = 1,num(i)
              ind = listptr(i)+in
              j = list(ind)
              nrmax = max(nrmax,j)
            enddo
          enddo

          allocate(aux(nrmax))
          call memory('A','D',nrmax,'extrapol')
          allocate(mtmp(nh))
          call memory('A','D',nh,'extrapol')

          aux(1:nrow) = 0.0_dp
  
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
                aux(j) = 0.0_dp
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

        endif !changed
C ..................

C Extrapolate matrix M ......................................................

        do ispin = 1,nspin
          do i = 1,nh
            msave = mnew(i,ispin)
            if (iord == 1) then
              mnew(i,ispin) = 2.0_dp * mnew(i,ispin) - mm2(i,ispin)
            endif
            mm2(i,ispin) = msave
          enddo
        enddo
C ....................

      endif !istep==1

      end
