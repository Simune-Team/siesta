      include '../poison.f'

      SUBROUTINE CPUTIM (TIME)

C  RETURNS CPU TIME IN SECONDS SINCE PROGRAM START
C  WRITEN BY J.SOLER (JSOLER AT EMDUAM11)

      DOUBLE PRECISION TIME
      REAL TIMES(2)

C  NEXT LINES FOR IBM SYSTEMS-370 (ASSEMBLE ROUTINE TASKTM REQUIRED)
*     CALL TASKTM (ITIME)
*     TIME = 1.D-4 * ITIME

C  NEXT LINES FOR IBM-3090
*     CALL CPUTIME (TIME,RCODE)
*     TIME = 1.D-6 * TIME

C  NEXT LINE FOR CRAY
*     TIME = SECOND()

C  NEXT LINE FOR SUN OR DEC WORKSTATIONS
      TIME = ETIME(TIMES)

C  NEXT LINE FOR IBM RS/6000 WORKSTATIONS
*     TIME = MCLOCK()*0.01D0

C  NEXT LINE FOR ANYTHING ELSE
*     TIME = 0.D0

      END

      subroutine rdiag (h,s,n,nm,w,z,fv)
C ***************************************************************************
C Subroutine  to solve all eigenvalues and eigenvectors of the
C real general eigenvalue problem  H z = w S z,  with H and S
C real symmetric matrices, by calling the  LAPACK routine DSYGV
C Writen by G.Fabricius and J.Soler, March 1998
C ************************** INPUT ******************************************
C real*8 h(nm,nm)                 : Symmetric H matrix
C real*8 s(nm,nm)                 : Symmetric S matrix
C integer n                       : Order of  the generalized  system
C integer nm                      : Dimension of H and S matrices
C ************************** OUTPUT *****************************************
C real*8 w(nm)                    : Eigenvalues
C real*8 z(nm,nm)                 : Eigenvectros
C ************************* AUXILIARY ***************************************
C real*8  fv(nm,2)                : Auxiliary storage array
C ***************************************************************************
      implicit          none
      integer           i1,i2,INFO,LWORK,n,nm
      character         JOBZ,UPLO
      double precision  h(nm,nm), s(nm,nm), w(nm),  z(nm,nm), fv(nm,3)
C ......................

C  start time count
      call timer('rdiag',1)

       JOBZ='V'
       UPLO='U'
       LWORK=3*n
       call dsygv(1,JOBZ,UPLO,n,h,nm,s,nm,w,fv,LWORK,INFO)
       do i1=1,n
         do i2=1,n
           z(i1,i2)=h(i1,i2)
         enddo
       enddo

      if (info .ne. 0) then
        write(6,*) 'rdiag: ERROR in routine dsygv:'
        if (info.lt.0) then
          write(6,*) -info,'-th argument had an illegal value'
        elseif (info.le.n)   then
          write(6,*) 'DSYEV failed to converge. ',info,
     .               ' off-diagonal elements of an intermediate '
          write(6,*) 'tridiagonal form did not converge to zero'
        else
          write(6,*)'The leading minor of order ',info-n,
     .              ' of B is not positive definite.'
          write(6,*)'The factorization of B could not be completed',
     .              ' and no eigenvalues or eigenvectors were computed.'
        endif
        stop 'rdiag: ERROR in routine dsygv'
      endif

c  stop time count
      call timer('rdiag',2)
      end

