      subroutine extrapolon(istep,iord,nspin,nbasis,nbasisCloc,maxnc,
     .                      numc,listc,aux,numcold,listcold,cold,c,Node)
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
C integer nbasis               : Number of rows of matrix M
C integer nbasisCloc           : Maximum number of rows of matrix M (dimension)
C integer maxnc                : First dimension of M matrix, and maximum
C                                number of nonzero elements of each column of M
C integer numc(maxnc)          : Control vector 1 of M matrix at t
C integer listc(maxnc,nbasisCloc) : Control vector 2 of M matrix at t
C real*8 aux(2,nbasis)         : Auxiliary storage array
C ************************** INPUT AND OUTPUT *********************************
C integer numcold(maxnc)         : Input: Control vector 1 of M matrix at t-dt
C                                       (if istep .ne. 1)
C                                Output: Control vector 1 of M matrix at t
C integer listcold(maxnc,nbasisCloc) : Input: Control vector 2 of M matrix at t-dt
C                                       (if istep .ne. 1)
C                                Output: Control vector 2 of M matrix at t
C real*8 cold(maxnc,nbasisCloc,nspin): Input: matrix M at t-2dt
C                                Output: matrix M at t-dt
C real*8 c(maxnc,nbasisCloc,nspin): New matrix M (extrapolated)
C                                Input: matrix at t-dt
C                                Output: matrix at t
C                                If istep = 1, c returned uncahanged
C **************************** BEHAVIOUR **************************************
C The routine allows for the sparse structure of the matrix M to change
C between MD time steps. On input, the matrices of former steps (c and cold) 
C have the structure of last step (t-dt): numold and listold; whereas the new
C (extrapolated) matrix has the structure of the current time step (which
C must be determined before calling this routine!!): num and list.
C On output, the routine updates the structure of c and cold, to that
C at the current (t) time steps respectively. Same with numold and listold
C 
C For the first MD time step (istep = 1), there is no extrapolation. 
C In that case, c is returned unchanged.
C Also, in that case numold and listold are only an output, and are set equal
C to num and list
C *****************************************************************************

      use precision
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      integer 
     .  iord,istep,maxnc,nbasis,nbasisCloc,nspin,Node

      integer 
     .  listc(maxnc,nbasisCloc),listcold(maxnc,nbasisCloc),
     .  numc(nbasisCloc),numcold(nbasisCloc)

      real(dp) ::
     .  aux(2,nbasis),cold(maxnc,nbasisCloc,nspin),
     .  c(maxnc,nbasisCloc,nspin)
 
C  Internal variables .......................................................

      integer
     .  i,in,ispin,j

      real(dp) ::
     .  msave
C ...........................................................................

      if (iord .ne. 0 .and. iord . ne. 1) then
        if (Node.eq.0) then
      write(6,*) 'extrapol: Wrong iord: only 0 and 1 order available'
        endif
        stop
      endif

C Just initialize numcold and listcold if istep = 1 ...........................
      if (istep .eq. 1) then
        do i = 1,nbasisCloc
          numcold(i) = numc(i)
          do in = 1,numc(i)
            listcold(in,i) = listc(in,i)
            do ispin = 1,nspin
              cold(in,i,ispin) = 0.0d0
            enddo
          enddo
        enddo
        return

      else

C Check if sparse structure has changed .....................................
        do i = 1,nbasisCloc
          if (numcold(i).ne.numc(i)) goto 10
          do in = 1,numc(i)
            if (listcold(in,i).ne.listc(in,i)) goto 10
          enddo
        enddo
        goto 20

C If sparse structure has changed, re-order c and cold 
C and change numcold and listcold to current ones .............................
10      continue

        do i = 1,nbasis
          do j = 1,2
            aux(j,i) = 0.0 
          enddo
        enddo
  
        do i = 1,nbasisCloc
          do ispin = 1,nspin
            do in = 1,numcold(i)
              j = listcold(in,i)
              aux(1,j) = c(in,i,ispin)
              aux(2,j) = cold(in,i,ispin)
            enddo
            do in = 1,numc(i)
              j = listc(in,i)
              c(in,i,ispin) = aux(1,j)
              cold(in,i,ispin) = aux(2,j)
            enddo
            do in = 1,numcold(i)
              j = listcold(in,i)
              aux(1,j) = 0.0d0
              aux(2,j) = 0.0d0
            enddo
          enddo
          numcold(i) = numc(i)
          do in = 1,numc(i)
            listcold(in,i) = listc(in,i)
          enddo
        enddo

C Extrapolate matrix M ......................................................
20      continue

        do ispin = 1,nspin
          do i = 1,nbasisCloc
            do in = 1,numc(i)
              msave = c(in,i,ispin)
              if (iord.eq.1 .and. cold(in,i,ispin).ne.0.0d0) then
                c(in,i,ispin) = 2.0d0*c(in,i,ispin) - cold(in,i,ispin)
              endif
              cold(in,i,ispin) = msave
            enddo
          enddo
        enddo

      endif

      return
      end
