C $Id: extrapol.f,v 1.3 1999/01/31 11:14:53 emilio Exp $

      subroutine extrapol(istep,iord,nspin,nrow,nmrow,nmax,num,list,aux,
     .                numold,listold,mm2,mnew)
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
C integer nmax                 : First dimension of M matrix, and maximum
C                                number of nonzero elements of each column of M
C integer num(nmax)            : Control vector 1 of M matrix at t
C integer list(nmax,nmrow)     : Control vector 2 of M matrix at t
C real*8 aux(2,nmrow)          : Auxiliary storage array
C ************************** INPUT AND OUTPUT *********************************
C integer numold(nmax)         : Input: Control vector 1 of M matrix at t-dt
C                                       (if istep .ne. 1)
C                                Output: Control vector 1 of M matrix at t
C integer listold(nmax,nmrow)  : Input: Control vector 2 of M matrix at t-dt
C                                       (if istep .ne. 1)
C                                Output: Control vector 2 of M matrix at t
C real*8 mm2(nmax,nmrow,nspin) : Input: matrix M at t-2dt
C                                Output: matrix M at t-dt
C real*8 mnew(nmax,nmrow,nspin): New matrix M (extrapolated)
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
      implicit none

      integer 
     .  iord,istep,nmax,nmrow,nrow,nspin

      integer 
     .  list(nmax,nmrow),listold(nmax,nmrow),num(nmrow),numold(nmrow)

      double precision
     .  aux(2,nmrow),mm2(nmax,nmrow,nspin),mnew(nmax,nmrow,nspin)
 
C  Internal variables .......................................................

      integer
     .  i,in,ispin,j

      double precision
     .  msave
C ...........................................................................

      if (iord .ne. 0 .and. iord . ne. 1) then
        write(6,*) 'extrapol: Wrong iord: only 0 and 1 order available'
        stop
      endif

C Just initialize numold and listold if istep = 1 ...........................
      if (istep .eq. 1) then
        do i = 1,nrow
          numold(i) = num(i)
          do in = 1,num(i)
            listold(in,i) = list(in,i)
            do ispin = 1,nspin
              mm2(in,i,ispin) = 0.d0
            enddo
          enddo
        enddo
        return
C .....................

      else

C Check if sparse structure has changed .....................................
        do i = 1,nrow
          if (numold(i) .ne. num(i)) goto 10
          do in = 1,num(i)
            if (listold(in,i) .ne. list(in,i)) goto 10
          enddo
        enddo
        goto 20
C .....................

C If sparse structure has changed, re-order mnew and mm2 
C and change numold and listold to current ones .............................
10      continue

        do i = 1,nrow
          do j = 1,2
            aux(j,i) = 0.0 
          enddo
        enddo
  
        do i = 1,nrow
          do ispin = 1,nspin
            do in = 1,numold(i)
              j = listold(in,i)
              aux(1,j) = mnew(in,i,ispin)
              aux(2,j) = mm2(in,i,ispin)
            enddo
            do in = 1,num(i)
              j = list(in,i)
              mnew(in,i,ispin) = aux(1,j)
              mm2(in,i,ispin) = aux(2,j)
            enddo
            do in = 1,numold(i)
              j = listold(in,i)
              aux(1,j) = 0.0
              aux(2,j) = 0.0
            enddo
          enddo
          numold(i) = num(i)
          do in = 1,num(i)
            listold(in,i) = list(in,i)
          enddo
        enddo
C ..................

C Extrapolate matrix M ......................................................
20      continue

        do ispin = 1,nspin
          do i = 1,nrow
            do in = 1,num(i)
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
