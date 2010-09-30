      module jacobi_davidson
      use matrixop
      use precond
      use precision
      use alloc
      use gmres
      use cg
      implicit none

      PUBLIC :: jd_diag
      PRIVATE

      real(dp),  pointer :: jd_v0(:) 
      integer            :: Kmax
      integer            :: Mmin
      integer            :: Mmax
      integer            :: JDreortho
      real               :: JDepsilon
      real               :: ITepsilon

      contains

      subroutine jd_diag( nuotot, nuo, nspin, neig, fixspin, temp, qtot,
     &                    qs, occtol, numh, listh, H, S, ef, efs,
     &                    Entropy, Dnew, Enew )
      use printMat
      use domain_decom, only : dd_cperm, dd_ncolum

      use parallel,     only : node, nodes
      use fdf
      implicit none
C     Input variables
      integer,   intent(in) :: nuotot, nuo, nspin
      logical,   intent(in) :: fixspin
      real(dp),  intent(in) :: temp, qtot, qs(nspin), occtol
      integer,   intent(in) :: numh(*), listh(*)
      real(dp),  intent(in) :: H(*), S(*)
C     Output variables
      integer,  intent(out) :: neig
      real(dp), intent(out) :: ef, efs(nspin), Entropy
      real(dp), intent(out) :: Dnew(*), Enew(*)
C     Local variables
      integer              :: io, nnz, i, j
      integer,     pointer :: ind(:), col(:)
      character*128        :: filename
      type(matrix_t)       :: A, B
      real(dp),    pointer :: w(:), Z(:,:), qo(:)
      real(dp)             :: tt
      logical              :: getD, getPSI

      call write_debug( '    PRE Jacobi' )
!      neig = min(neig,200)
!      neig = (qtot*6)/5
!      neig = min(neig,100)
!      neig = min(neig,nuotot)

#ifdef KKKKKKK
      Kmax = neig
!      Kmax = min(50,neig)
      Mmax      = 50
      Mmin      = 20
      JDepsilon = 1.0e-5
      ITepsilon = 1.0e-5
      ITmaxiter = min(50,nuotot)
      JDreortho = 100
      gmres_dim = min(max(257,nuotot/2),nuotot)

#else
      Kmax      = fdf_integer( 'jd.kmax',  neig )
      Mmax      = fdf_integer( 'jd.mmax',  50 )
      Mmin      = fdf_integer( 'jd.mmin',  20 )
      JDepsilon = fdf_double(  'jd.eps',   1.0e-5_dp )
      ITepsilon = fdf_double(  'jd.iteps', 1.0e-5_dp )
      ITmaxiter = fdf_integer( 'jd.maxit', min(50,nuotot) )
      JDreortho = fdf_integer( 'jd.reort', 100 )
#endif
      Kmax = min(Kmax,nuotot)
      neig = Kmax

#ifdef DEBUG
      write(23,*) '=========  Inside jacobi  ========='
      write(23,*) '     n. orbitales global:', nuotot
      write(23,*) '     n. orbitales local :', nuo
      write(23,*) '     qtot               :', qtot
      write(23,*) '     n. de autovalores  :', Kmax 
      write(23,*) '     JD max krylov dim  :', Mmax
      write(23,*) '     JD min krylov dim  :', Mmin
      write(23,*) '     JD epsilon         :', JDepsilon
      write(23,*) '     Iter Epsilon       :', ITepsilon
      write(23,*) '     Iter Max Iter      :', ITmaxiter
      flush(23)
#endif

      nullify( ind )
      call re_alloc( ind, 1, nuo+1, 'ind', 'jd_diag' )
      ind(1) = 1
      do io= 1, nuo
        ind(io+1) = ind(io) + numh(io)
      enddo
      nnz = ind(nuo+1)-1

      nullify( col )
      call re_alloc( col, 1, nnz, 'col', 'jd_diag' )
      do io= 1, nnz
        col(io) = dd_cperm(listh(io))
      enddo

#ifdef DEBUG_PS
      if (Nodes.lt.10) then
        write(filename,'(A7,I1,A3)') 'jacobi.', Node, '.ps'
      else if (Nodes.lt.100) then
        write(filename,'(A7,I2,A3)') 'jacobi.', Node, '.ps'
      else if (Nodes.lt.1000) then
        write(filename,'(A7,I3,A3)') 'jacobi.', Node, '.ps'
      else if (Nodes.lt.10000) then
        write(filename,'(A7,I4,A3)') 'jacobi.', Node, '.ps'
      else
        write(filename,'(A7,I5,A3)') 'jacobi.', Node, '.ps'
      endif

      call printMatrix( nuo, nuotot, col, ind, filename, 'MyPart' )
#endif

      call create_matrix( nuo, dd_ncolum, nnz, ind, col, H,  A )
      call create_matrix( nuo, dd_ncolum, nnz, ind, col, S,  B )
      write(23,*) '     n. columnas local  :', A%nc
      flush(23)

      nullify( w, Z, qo )
      call re_alloc( w, 1, neig, 'w', 'jd_diag' )
      call re_alloc( Z, 1, nuo, 1, neig, 'Z', 'jd_diag' )
      call re_alloc( qo, 1, neig, 'qo', 'jd_diag' )

C     Set jacobi options
      nullify( jd_v0 )
      call re_alloc( jd_v0, 1, nuo, 'v0', 'jd_diag' )
      tt    = nuotot
      tt    = 1.0/sqrt(tt)
      jd_v0 = tt

      getD   = .true.
      getPSI = .false.

      call jd_diagg( nuo, nspin, neig, nnz, getD, getPSI, fixspin, temp,
     &               qtot, qs, occtol, A, B, ef, efs, Entropy,
     &               w, Z, qo, Dnew, Enew )

      call de_alloc( jd_v0, 'v0', 'jd_diag' )

      call de_alloc( qo, 'qo', 'jd_diag' )
      call de_alloc( Z, 'Z', 'jd_diag' )
      call de_alloc( w, 'w', 'jd_diag' )

      call de_alloc( ind, 'ind', 'jd_diag' )
      call de_alloc( col, 'col', 'jd_diag' )

      call write_debug( '    POS Jacobi' )
      end subroutine jd_diag

      subroutine jd_diagg( n, nspin, neig, nnz, getD, getPSI, fixspin,
     &                     temp, qtot, qs, occtol, A, B, ef, efs,
     &                     Entropy, w, Z, qo, Dnew, Enew )
      use m_fermid,      only : fermid, fermispin
      implicit none
C     Input variables
      integer,        intent(in) :: n, nspin, neig, nnz
      logical,        intent(in) :: getD, getPSI, fixspin
      real(dp)                   :: temp, qtot, qs(nspin), occtol
      type(matrix_t), intent(in) :: A, B
C     Output variables
      real(dp),      intent(out) :: ef, efs(nspin), Entropy
      real(dp),      intent(out), target :: Z(n,neig)
      real(dp),      intent(out) :: w(neig), qo(neig),
     &                              Dnew(nnz), Enew(nnz)
C     Local variables
      integer,         parameter :: nk = 1
      real(dp),        parameter :: wk(1) = (/ 1.0_dp /)
      logical                    :: notfound
      integer                    :: io, ie, j, jo, uneig, nc
      real(dp)                   :: qe, ee, qei, eei
      real(dp),          pointer :: ZZ(:,:)

      call medirTiempo( 1, 7 )
      call jd_solver( n, neig, A, B, W, Z )
      call medirTiempo( 2, 7 )
      write(23,*) '=== AUTOVALORES de salida ==='
      flush(23)
      do j= 1, min(neig-5,5)
        write(23,*) 'EV(', j, '):', W(j)   ! , Z(:,j)
        flush(23)
      enddo
      do j= max(neig-5,1), neig
        write(23,*) 'EV(', j, '):', W(j)   ! , Z(:,j)
        flush(23)
      enddo





      call medirTiempo( 2, 6 )
      write(23,*) ' '
      write(23,*) '========    TIMES    ========'
      flush(23)
      call showTiempo( 23, 1,  'SetupHam         :' )
      call showTiempo( 23, 2,  '  Rhoofd         :' )
      call showTiempo( 23, 3,  '  Poison         :' )
      call showTiempo( 23, 4,  '  Cellxc         :' )
      call showTiempo( 23, 5,  '  Vmat           :' )
      call showTiempo( 23, 7,  'jd_solver        :' )
!      call printmemory( 23, 1 )
!      flush(23)

      call MPI_Finalize( io )
      STOP






C Check if we are done
      if (.not.getD) return
      write(23,*) '=== PRE fermid ==='
      flush(23)

C Find new Fermi energy and occupation weights
      if (fixspin) then
        call fermispin( nspin, nspin, nk, wk, neig, neig, 
     &                  w, temp, qs, qo, efs, Entropy )
      else
        call fermid( nspin, nspin, nk, wk, neig, neig, 
     &               w, temp, qtot, qo, ef, Entropy )
      endif
      write(23,*) '=== POS fermid ==='
      flush(23)

C     Count the number of usefull eigenvalues/vectors
      uneig    = 1
      notfound = .true.
      do while(notfound)
        if (qo(uneig).gt.occtol) then
          if (uneig+1.gt.neig) then
            notfound = .false.
          else
            uneig = uneig + 1
          endif
        else
          notfound = .false.
          uneig    = uneig - 1
        endif
      enddo

#ifdef MPI
      nc = A%nc
      nullify(ZZ)
      call re_alloc( ZZ, 1, nc, 1, uneig, 'ZZ', 'jd_diagg' )
      call reduce_matrix( n, nc, uneig, Z, ZZ )
#else
      ZZ => Z
#endif

      Dnew = 0.0
      Enew = 0.0
      do ie = 1, uneig
        qe = qo(ie)
        ee = qe*w(ie)
        do io = 1, n
          qei = qe*ZZ(io,ie)
          eei = ee*ZZ(io,ie)
          do j = A%ind(io), A%ind(io+1)-1
            jo = A%col(j)
            Dnew(j) = Dnew(j) + qei*ZZ(jo,ie)
            Enew(j) = Enew(j) + eei*ZZ(jo,ie)
          enddo
        enddo
      enddo
      write(23,*) '=== FIN jd_diagg ==='
      flush(23)

#ifdef MPI
      call de_alloc( ZZ, 'ZZ', 'jd_diagg' )
#endif
      end subroutine jd_diagg

      subroutine jd_solver( n, neig, A, B, lambda, Q )
      implicit none
C     Input variables
      type(matrix_t), intent(in) :: A, B
      integer,        intent(in) :: n, neig
C     Output variables
      real(dp),      intent(out) :: lambda(neig), Q(n,neig)
C     Local variables
      integer              :: it, i, j, k, k2, lastk, m, m2,
     &                        offset, gm_it, gm_it_tot, lwork,
     &                        reortho, nreort
      real(dp)             :: dotp, rteta, nrm2, taut, dvdot
      logical              :: nomin
      real(dp),    pointer :: t(:), V(:,:), Va(:,:),
     &                        Vb(:,:), VAV(:,:), VAUX(:,:),
     &                        teta(:), S(:), MM(:,:),
     &                        vi(:), via(:), vib(:),
     &                        vmb(:), s1(:), u(:), p(:),
     &                        ua(:), r(:), si(:), siplus1(:), Z(:,:),
     &                        work(:)
      real(dp),  parameter :: ONE=1.0, ZERO=0.0, EPS_IT_MIN=1.0E-1
#ifdef MPI
      integer              :: MPIerror
#endif

      t   => jd_v0

      k         = 0
      m         = 0
      gm_it_tot = 0

      m2    = ((Mmax+1)*Mmax)/2
      k2    = ((Kmax+1)*Kmax)/2
      lwork = max(Kmax*2+5*n+A%nc,max(3*mmax,2*(mmax+kmax)))
      nullify( V, Va, Vb, VAUX, VAV, teta, MM,
     &         S, u, p, ua, r, Z, work )
      call re_alloc( V, 1, n, 1, Mmax, 'V', 'jd_solver' )
      call re_alloc( Va, 1, n, 1, Mmax, 'Va', 'jd_solver' )
      call re_alloc( Vb, 1, n, 1, Mmax, 'Vb', 'jd_solver' )
      call re_alloc( VAUX, 1, n, 1, Mmax, 'VAUX', 'jd_solver' )
      call re_alloc( VAV, 1, m2, 1, 2, 'VAV', 'jd_solver' )
      call re_alloc( MM, 1, k2, 1, 2, 'MM', 'jd_solver' )
      call re_alloc( teta, 1, Mmax, 'teta', 'jd_solver' )
      call re_alloc( S, 1, Mmax*Mmax, 'S', 'jd_solver' )
      call re_alloc( u, 1, n, 'u', 'jd_solver' )
      call re_alloc( p, 1, n, 'p', 'jd_solver' )
      call re_alloc( ua, 1, n, 'ua', 'jd_solver' )
      call re_alloc( r, 1, n, 'r', 'jd_solver' )
      call re_alloc( Z, 1, n, 1, Kmax, 'Z', 'jd_solver' )
      call re_alloc( work, 1, lwork, 'work', 'jd_solver' )

      ! (1)
      V       = 0.0
      Va      = 0.0
      Vb      = 0.0

      it = 1
      reortho = 0
      nreort  = 0
      ITeps   = EPS_IT_MIN
      nomin   = .true.
      taut    = 0.0

      ! (2)
      do while ( k .lt. Kmax )
        write(23,*) ' '
        write(23,*) '  NEW JD_ITER:', it, 'ITSOL_IT:', gm_it_tot,
     &             'JD_dim:', m+1, '  EV found:', k, ' nreort:', nreort
        flush(23)

        ! (3) - (5)
        call GramSchmidt( N, M, Vb, V, t, work )

        ! (6) - (10)
        m = m + 1
        call dcopy( N, t, 1, V(1,m), 1 )        ! Vm    = t
        call projectionMatrix( n, m, A, B, V, VA, VB, VAV, work )

        ! (11)
        call computeEigenpairs( m, m2, VAV, teta, S, work )
!        write(23,*) '      CANDIDATE EIGENVALUES:', teta(1:min(4,m))

!        call selectEigenpair( nomin, taut, m, teta, S, work )
!        write(23,*) '      SELECTED EIGENVALUE:', teta(1)

        ! (12)
        call DGEMV( 'N', n, m, ONE, V, n, S, 1, ZERO, u, 1 )
        call DGEMV( 'N', n, m, ONE, Vb, n, S, 1, ZERO, p, 1 )
        call DGEMV( 'N', n, m, ONE, Va, n, S, 1, ZERO, ua, 1 )
        do i= 1, n
          r(i) = ua(i) - teta(1)*p(i)
        enddo

        lastk     = k
        ! (13) 
        do while ( residual_norm( n, r, u) .le. JDepsilon )
          write(23,*) '      GOOD EIGENVALUE:', teta(1)

          ! (14)
          call dcopy( n, u, 1, Q(1,k+1), 1 )

          call dcopy( n, p, 1, Z(1,k+1), 1 )
          lambda(k+1)   = teta(1)
          k             = k + 1

          ! (15)
          if ( k .eq. Kmax ) goto 100

          ! (16)
          m = m - 1

          ! (17)
          if ( m.gt.0 ) then
            ! (18)
            call DGEMM( 'N', 'N', n, m, m+1, ONE, V, n, S(m+2), m+1,
     &                  ZERO, VAUX, n )
            call dcopy( n*m, VAUX, 1, V, 1 )
            call DGEMM( 'N', 'N', n, m, m+1, ONE, Va, n, S(m+2), m+1,
     &                  ZERO, VAUX, n )
            call dcopy( n*m, VAUX, 1, Va, 1 )
            call DGEMM( 'N', 'N', n, m, m+1, ONE, Vb, n, S(m+2), m+1,
     &                  ZERO, VAUX, n )
            call dcopy( n*m, VAUX, 1, Vb, 1 )

            reortho = reortho + m + 1

            ! (19) (20)
            offset = 0
            do i = 1, m
              do j= 1, i-1
                offset = offset + 1
                VAV(offset,1) = 0.0
              enddo
              offset = offset + 1
              VAV(offset,1) = teta(i+1)
            enddo

            offset = 0
            do i = 1, m
              do j= 1, m
                offset = offset + 1
                if (i.eq.j) then
                  S(offset) = 1.0
                else
                  S(offset) = 0.0
                endif
              enddo
            enddo

            do i = 1, m
              teta(i)     = teta(i+1)
            enddo
          endif
          ! (21)
          call dcopy( n, V,  1, u, 1 )
          call dcopy( n, Vb, 1, p, 1 )
          do i= 1, n
            r(i) = Va(i,1) - teta(1)*p(i)
          enddo
!          write(23,*) '      CANDIDATE EIGENVALUE:', teta(1)

          ITeps = EPS_IT_MIN
        ! (22)  
        enddo

        ! (23)
        if ( m .ge. Mmax ) then
!          write(23,*) '    RESTART JACOBI'

          ! (24)
C          VAV = 0.0
          ! (25) - (27)
          call DGEMM( 'N', 'N', n, Mmin-1, m, ONE, V, n, S(m+1),
     &                m, ZERO, VAUX, n )
          call dcopy( n, u, 1, V, 1 )
          call dcopy( n*(Mmin-1), VAUX, 1, V(1,2), 1 )

          call DGEMM( 'N', 'N', n, Mmin-1, m, ONE, Va, n, S(m+1),
     &                m, ZERO, VAUX, n )
          call dcopy( n, ua, 1, Va, 1 )
          call dcopy( n*(Mmin-1), VAUX, 1, Va(1,2), 1 )

          call DGEMM( 'N', 'N', n, Mmin-1, m, ONE, Vb, n, S(m+1),
     &                m, ZERO, VAUX, n )
          call dcopy( n, p, 1, Vb, 1 )
          call dcopy( n*(Mmin-1), VAUX, 1, Vb(1,2), 1 )

          reortho = reortho + m

          ! (28)
          m      = Mmin
          offset = 0
          do i = 1, m
            do j= 1, i-1
              offset = offset + 1
              VAV(offset,1) = 0.0
            enddo
            offset = offset + 1
            VAV(offset,1) = teta(i)
          enddo

        ! (29)
        endif

        if (reortho.ge.JDreortho) then
!          write(23,*) '    REORTOGONALIZE BASE'
          nreort = nreort + 1
          do i= 1, m
            call GramSchmidt2( N, k, i-1, Z, Q, Vb, V, V(1,i), work )
            call projectionMatrix( n, i, A, B, V, VA, VB, VAV, work )
          enddo

          reortho = 0
        endif

        ! (30)
        rteta = teta(1)
        call dcopy( n, p, 1, Z(1,k+1), 1 )
        call dcopy( n, u, 1, Q(1,k+1), 1 )

        ! (31)
        t  = 0.0
        call jd_correction( n, k+1, lastk, Kmax, Z, MM, A, B, rteta,
     &                      r, gm_it, t, work, lwork )

        ITeps = ITeps*0.5
        if (ITeps.lt.ITepsilon) ITeps = ITepsilon

        gm_it_tot = gm_it_tot + gm_it
        it = it + 1

!         call showTiempo( 23,  8,  'GramSchmidt of V   :' )
!         call showTiempo( 23,  9,  'M=V*A*V            :' )
!         call showTiempo( 23, 10,  'EV of M            :' )
!         call showTiempo( 23, 11,  'Check EV           :' )
!         call showTiempo( 23, 12,  'Restart V          :' )
!         call showTiempo( 23, 13,  'Reortogonalize     :' )
!         call showTiempo( 23, 14,  'JD_Correction      :' )
!         call showTiempo( 23, 15,  '  DDMTM            :' )
!         call showTiempo( 23, 16,  '  choleskyFact     :' )
!         call showTiempo( 23, 20,  '  SYMLQ            :' )
!         call showTiempo( 23, 21,  '    DEFLATED MxV   :' )
!         call showTiempo( 23, 22,  '     y = (A-t*B)*v :' )
!         call showTiempo( 23, 17,  '     DDMTV         :' )
!         call showTiempo( 23, 18,  '     choleskySolve :' )
!         call showTiempo( 23, 19,  '     DGEMV         :' )


      ! (32)  
      end do


100   continue

      call showTiempo( 23, 7,  'Partial Time JacobiDa:' )
      flush(23)

      call de_alloc( work, 'work', 'jd_solver' )
      call de_alloc( Z, 'Z', 'jd_solver' )
      call de_alloc( r, 'r', 'jd_solver' )
      call de_alloc( ua, 'ua', 'jd_solver' )
      call de_alloc( p, 'p', 'jd_solver' )
      call de_alloc( u, 'u', 'jd_solver' )
      call de_alloc( S, 'S', 'jd_solver' )
      call de_alloc( teta, 'teta', 'jd_solver' )
      call de_alloc( MM, 'MM', 'jd_solver' )
      call de_alloc( VAV, 'VAV', 'jd_solver' )
      call de_alloc( Vb, 'Vb', 'jd_solver' )
      call de_alloc( Va, 'Va', 'jd_solver' )
      call de_alloc( V, 'V', 'jd_solver' )
      end subroutine jd_solver

      subroutine computeEigenpairs( m, m2, AP, W, VR, work )
      use sys,       only : die
      implicit none
C     Input parameters
      integer,    intent(in) :: m, m2
C     Output parameters
      real(dp)               :: AP(m2,2), W(m), VR(m*m), work(3*M)
C     Local Variables
      integer                :: info

      ! LAPACK
      !
      ! DSPEV computes all the eigenvalues and, optionally, eigenvectors
      ! of a complex Hermitian matrix in packed storage
      !
      ! The right eigenvector v(j) of A satisfies
      !                  A * v(j) = lambda(j) * v(j)
      ! where lambda(j) is its eigenvalue.
      ! The left eigenvector u(j) of A satisfies
      !               u(j)**H * A = lambda(j) * u(j)**H
      ! where u(j)**H denotes the conjugate transpose of u(j).
      !
      ! The computed eigenvectors are normalized to have Euclidean norm
      ! equal to 1 and largest component real.
      call dcopy( m2, AP, 1, AP(1,2), 1 )
      call DSPEV( 'V', 'U', m, AP(1,2), W, VR, m, work, info )

      if (info.ne.0) call die( 'ERROR in computeEigenpairs->DSPEV' )

      end subroutine computeEigenpairs


      subroutine selectEigenpair( nomin, taut, m, teta, S, work )
      use sys,       only : die
      implicit none
C     Input parameters
      logical,   intent(in) :: nomin
      integer,   intent(in) :: m
      real(dp),  intent(in) :: taut
C     Output parameters
      real(dp), intent(out) :: teta(m), S(m,m), work(m)
C     Local variables
      integer               :: ind, i, j
      real(dp)              :: rr

      if (.not.nomin) then
        ind = 1
        do while (teta(1).lt.taut)
          rr = teta(1)
          do i= 1, m-1
            teta(i) = teta(i+1)
          enddo
          teta(m) = rr

          do j= 1, m
            work(j) = S(j,1)
          enddo
          do i= 1, m-1
            do j= 1, m
              S(j,i) = S(j,i+1)
            enddo
          enddo
          do j= 1, m
            S(j,m) = work(j)
          enddo

          ind = ind + 1
          if (ind.eq.m) call die( 'ERROR in selectEigenpair' )
        end do
      endif

      end subroutine selectEigenpair


      subroutine jd_correction( n, k, lastk, Kmax, Z, MM, A, B, teta,
     &                          r, gm_it, t, work, lwork )
      implicit none
C     Input parameters
      integer               :: n, k, lastk, Kmax, lwork
      real(dp)              :: teta, r(n)
      real(dp),      target :: Z(n*k), work(lwork)
      type(matrix_t)        :: A, B
C     Output parameters
      real(dp)              :: t(n)
      real(dp),      target :: MM(((Kmax+1)*Kmax)/2,2)
      integer               :: gm_it
C     Local variables
      integer               :: nc, o1, o2
      real(dp),     pointer :: M(:), Maux(:), LU(:), gama(:),
     &                         gmaux(:)
      real*8,     parameter :: ONE=1.0, ZERO=0.0, M_ONE=-1.0

      nc    = A%nc

C     Set local pointers
      M    => MM(:,1)
      Maux => MM(:,2)
      LU   => MM(:,2)

      o1 = 1
      o2 = o1+Kmax
      gama => work(o1:o2-1)

      o1 = o2
      o2 = o1+Kmax
      gmaux => work(o1:o2-1)

      ! (1) - (2)
      call ddmtmSym( n, k, lastk, Z, M, Maux )

      ! (3)
C     computes the cholesky factorization of a symetric N-by-N matrix.
      call choleskyFact( k, lastk, M, LU )

      ! (4) - (6)
C     gama = Z^t * r
      call ddmtv( n, k, Z, r, gama, gmaux )

      ! (7) (8)
C     gama = M^(-1)*gama using LU from choleskyFact
      call choleskySolve( k, LU, gama )

      ! (9) r = Z*gama - r
      call DGEMV( 'N', n, k, ONE, Z, n, gama, 1, M_ONE, r, 1 )

      ! (10)
      call com_symmlq( n, k, kmax, A, B, teta, r, Z,
     &                 LU, t, gm_it, work, 5*n+2*kmax+nc )
      end subroutine jd_correction

      subroutine projectionMatrix( n, m, A, B, V, VA, VB, VAV, work )
      implicit none
C     Input variables
      integer,           intent(in) :: n, m
      type(matrix_t),    intent(in) :: A, B
C     Output variables
      real(dp), target, intent(out) :: V(n,m), VA(n,m), VB(n,m),
     &                                 VAV((m*(m+1))/2), work(*)
C     Local variables
      integer                       :: nc, offset
      real(dp),             pointer :: t(:), ta(:), tb(:)
      real(dp)                      :: dvdot, nrm

      nc    = A%nc

      t  => V(1:n,m)
      ta => VA(1:n,m)
      tb => VB(1:n,m)

      ! (6)
      call dsmv( n, nc, B, t, work, tb )
      nrm = dvdot( n, t, tb )
      nrm = 1.0/sqrt(nrm)

      ! (7)
      call dscal( n, nrm, t, 1 )         ! Vm  = t/nrm
      call dsmv( n, nc, A, t, work, ta )   ! VAm = A*Vm 
      call dscal( n, nrm, tb, 1 )        ! VBm = tb/nrm

      ! (8) (9) (10)
      offset = (m-1)*m/2
      call ddmtv( n, m, V, ta, VAV(offset+1), work )
      end subroutine projectionMatrix

      function residual_norm( n, x, y )
      use precision,    only : dp
#ifdef MPI
      use m_mpi_utils,  only : globalize_sum
#endif
      implicit none
C     Input variables
      integer,  intent(in) :: n
      real(dp), intent(in) :: x(n), y(n)
C     Output variables
      real(dp)             :: residual_norm
C     Local variables
      integer              :: i
      real(dp)             :: local(2), global(2)
      real(dp),   external :: ddot

      local(1) = ddot( n, x, 1, x, 1 )
      local(2) = ddot( n, y, 1, y, 1 )
#ifdef MPI
      call globalize_sum( local, global )
      residual_norm = sqrt(global(1)/global(2))
#else
      residual_norm = sqrt(local(1)/local(2))
#endif
!      write(23,*) '      RESIDUAL NORM        :',
!     &            residual_norm, JDepsilon
      end function residual_norm

      end module jacobi_davidson      
