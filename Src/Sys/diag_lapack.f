c
c  Diagonalization routines using LAPACK
c
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




      subroutine cdiag (h,nh,s,ns,n,w,z,nz,fv)
C ***************************************************************************
C Subroutine  to solve all eigenvalues and eigenvectors of the
C real general eigenvalue problem  H z = w S z,  with H and S
C complex hermitian matrices, by calling the  LAPACK routine ZHEGV
C Writen by G.Fabricius and J.Soler, August 1998
C ************************** INPUT ******************************************
C complex*16 h(nh,n)   : Hermitian H matrix (destroyed on output)
C integer    nh        : First dimension of H matrix
C complex*16 s(ns,n)   : Hermitian S matrix (destroyed on output)
C integer    ns        : First dimension of S matrix
C integer    n         : Order of the generalized system
C integer    nz        : First dimension of Z matrix
C ************************** OUTPUT *****************************************
C real*8     w(n)      : Eigenvalues
C complex*16 z(nz,n)   : Eigenvectors
C ************************* AUXILIARY ***************************************
C real*8     fv(n,2)   : Auxiliary storage array
C ***************************************************************************
      implicit          none
      integer           i1,i2,INFO,LWORK,n,nh,ns,nz
      character         JOBZ,UPLO
      double precision  h(2,nh,n), s(2,ns,n), w(n), z(2,nz,n), fv(n,3)
C ......................

C  start time count
      call timer('cdiag',1)

       JOBZ='V'
       UPLO='U'
       LWORK=nz*n
C      Matrix z used for auxiliary space within zhegv
       call zhegv(1,JOBZ,UPLO,n,h,nh,s,ns,w,z,LWORK,fv,INFO)
       do i2=1,n
         do i1=1,n
           z(1,i1,i2)=h(1,i1,i2)
           z(2,i1,i2)=h(2,i1,i2)
         enddo
       enddo

      if (info .ne. 0) then
        write(6,*) 'cdiag: ERROR in routine zhegv:'
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
        stop 'cdiag: ERROR in routine zhegv'
      endif

c  stop time count
      call timer('cdiag',2)
      end




