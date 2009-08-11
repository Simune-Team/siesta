      module jacobi_davidson
      use matrixop
      use precond,      only : precond_t, setIdPrecond,
     &                         precond_solver_mat,
     &                         precond_solver_vec
      use precision
      use alloc
      use gmres
      use cg
      implicit none

      PUBLIC :: jd_diag
      PRIVATE

      real(dp),  pointer :: jd_v0(:) 
      integer            :: jd_reorthogonalize
      integer            :: jd_kmax
      integer            :: jd_mmin
      integer            :: jd_mmax
      real               :: jd_eps

      type(precond_t)    :: PC


      contains

      subroutine jd_diag( nuotot, nuo, nspin, neig, fixspin, temp, qtot,
     &                    qs, occtol, numh, listh, H, S, ef, efs,
     &                    Entropy, Dnew, Enew )
      use printMat
      use domain_decom, only : dd_cperm, dd_ncolum

      use parallel,     only : node, nodes

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

      neig = (qtot*7)/10

      write(23,*) '=========  Inside jacobi  ========='
      write(23,*) '     n. orbitales global:', nuotot
      write(23,*) '     n. orbitales local :', nuo
      write(23,*) '     n. de autovalores  :', neig
      write(23,*) '     qtot               :', qtot
      flush(23)

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

      jd_kmax = neig
!      jd_kmax = min(50,neig)
      jd_mmax            = 50
      jd_mmin            = 20
      jd_eps             = 1.0e-4
      jd_reorthogonalize = 2
      write(23,*) '     JD max krylov dim  :', jd_mmax
      write(23,*) '     JD min krylov dim  :', jd_mmin
      flush(23)

C     Set MINRES, GMRES or CG options and Preconditioner
#ifdef KKKKKKKKKK
      gmres_maxiter = nuotot
      gmres_dim     = min(max(257,nuotot/2),nuotot)
      cg_maxiter    = nuotot
      call setIdPrecond( PC, nuo )
#else
      minres_eps     = 1.0e-4
      minres_maxiter = nuotot

      minres_maxiter = 50

#endif
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
      write(23,*) 'POS jd_diag'
      flush(23)

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


      call medirTiempo( 1, 6 )
      call jd_solver( n, neig, A, B, W, Z )
      call medirTiempo( 2, 6 )
      write(23,*) '=== AUTOVALORES de salida ==='
      do j= 1, neig
        write(23,*) 'EV(', j, '):', W(j), Z(:,j)
      enddo
      write(23,*) ' '
      write(23,*) '========    TIMES    ========'
      flush(23)

      call showTiempo( 23, 1,  'SetupHam         :' )
      call showTiempo( 23, 2,  '  Rhoofd         :' )
      call showTiempo( 23, 3,  '  Poison         :' )
      call showTiempo( 23, 4,  '  Cellxc         :' )
      call showTiempo( 23, 5,  '  Vmat           :' )
      call showTiempo( 23, 6,  'JacobiDa         :' )
      call showTiempo( 23, 12, '  Eigenval       :' )
      call showTiempo( 23, 13, '  GramSchB       :' )
      call showTiempo( 23, 14, '  A·x            :' )
      call showTiempo( 23, 7,  '  Correcti       :' )
      call showTiempo( 23, 8,  '    GMRES        :' )
!      call showTiempo( 23, 11, '      GramSchm   :' )
      call showTiempo( 23, 9,  '      gmres_mv   :' )
      call showTiempo( 23, 10, '        (A-o·B)x :' )
      call showTiempo( 23, 11, '        rd_array :' )
      flush(23)
!      call MPI_finalize( io )
!      STOP

C Check if we are done
      if (.not.getD) return

C Find new Fermi energy and occupation weights
      if (fixspin) then
        call fermispin( nspin, nspin, nk, wk, neig, neig, 
     .                  w, temp, qs, qo, efs, Entropy )
      else
        call fermid( nspin, nspin, nk, wk, neig, neig, 
     .               w, temp, qtot, qo, ef, Entropy )
      endif

      write(23,*) 'Fermi energy', ef, Entropy
      write(23,*) 'Occupations', qo
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

      write(23,*) 'uneig:', uneig
      flush(23)

      nc = A%nc
#ifdef MPI
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

      call de_alloc( ZZ, 'ZZ', 'jd_diagg' )
      end subroutine jd_diagg

      subroutine jd_solver( n, neig, A, B, lambda, Q )
      implicit none
C     Input variables
      type(matrix_t), intent(in) :: A, B
      integer,        intent(in) :: n, neig
C     Output variables
      real(dp),      intent(out) :: lambda(neig), Q(n,neig)
C     Local variables
      integer              :: kmax, mmin, mmax, it, i, j, k, m, m2,
     &                        offset, gm_it, gm_it_tot, nc
      real(dp)             :: eps, dotp, norm, rteta, nrm2
      real(dp),    pointer :: t(:), t_tilde(:), V(:,:), Va(:,:),
     &                        Vb(:,:), VAV(:), VAVtmp(:), VAUX(:,:),
     &                        teta(:), S(:),
     &                        vi(:), via(:), vib(:),
     &                        vmb(:), s1(:), u(:), p(:),
     &                        ua(:), r(:), si(:), siplus1(:), Z(:,:),
     &                        xx(:), work(:)
      real(dp)             :: dvdot, residual_norm
      real(dp),  parameter :: ONE=1.0, ZERO=0.0
      integer              :: reorthogonalize, nreort
#ifdef MPI
      integer              :: MPIerror
#endif

      t   => jd_v0
      kmax = jd_kmax
      mmin = jd_mmin 
      mmax = jd_mmax 
      eps  = jd_eps

      k         = 0
      m         = 0
      gm_it_tot = 0
      nullify( t_tilde, V, Va, Vb, VAUX, VAV,
     &         teta, S, u, p, ua, r, Z, xx, work )
      m2 = ((mmax+1)**2+(mmax+1))/2
      nc = A%nc
      call re_alloc( t_tilde, 1, n, 't_tilde', 'jd_solver' )
      call re_alloc( V, 1, n, 1, mmax+1, 'V', 'jd_solver' )
      call re_alloc( Va, 1, n, 1, mmax+1, 'Va', 'jd_solver' )
      call re_alloc( Vb, 1, n, 1, mmax+1, 'Vb', 'jd_solver' )
      call re_alloc( VAUX, 1, n, 1, mmax+1, 'VAUX', 'jd_solver' )
      call re_alloc( VAV, 1, m2, 'VAV', 'jd_solver' )
      call re_alloc( VAVtmp, 1, m2, 'VAVtmp', 'jd_solver' )
      call re_alloc( teta, 1, mmax+1, 'teta', 'jd_solver' )
      call re_alloc( S, 1, (mmax+1)**2, 'S', 'jd_solver' )
      call re_alloc( u, 1, n, 'u', 'jd_solver' )
      call re_alloc( p, 1, n, 'p', 'jd_solver' )
      call re_alloc( ua, 1, n, 'ua', 'jd_solver' )
      call re_alloc( r, 1, n, 'r', 'jd_solver' )
      call re_alloc( Z, 1, n, 1, kmax, 'Z', 'jd_solver' )
      call re_alloc( xx, 1, nc, 'xx', 'jd_solver' )
      call re_alloc( work, 1, 3*mmax, 'work', 'jd_solver' )

      ! (1)
      t_tilde = 0.0
      V       = 0.0
      Va      = 0.0
      Vb      = 0.0

      it = 1
      reorthogonalize = 0
      nreort          = 0
      ! (2)
      do while ( k .lt. kmax .and. gm_it_tot.le.200*kmax )
        write(23,*) ' '
        write(23,*) '  NEW JD_ITER:', it, 'ITSOL_IT:', gm_it_tot,
     &             'JD_dim:', m+1, '  EV found:', k, ' nreort:', nreort
        flush(23)


        ! (3) - (5)
        if (m.gt.0)
     &    call GramSchmidt_B( N, M, Vb, V, t, work, work(mmax+1) )

#ifdef RGT
        do i= 1, k
          norm = dvdot( n, Z(1,i), t )
          write(23,*) i, k, 't _|_ Z?', norm
        enddo
        do i= 1, m
          norm = dvdot( n, Vb(1,i), t )
          write(23,*) i, m, 't _|_ VB?', norm
        enddo
#endif

        ! (6) 
        m = m + 1
        call dsmv( n, nc, B, t, xx, t_tilde )
        norm = dvdot( n, t_tilde, t )
        norm = 1.0/sqrt(norm)

        ! (7)
        call dscal( N, norm, t, 1 )             ! t     = t/norm
        call dcopy( N, t, 1, V(1,m), 1 )        ! Vm    = t

        call dsmv( n, nc, A, t, xx, Va(1,m) )   ! VAm   = A*Vm 

        call dscal( N, norm, t_tilde, 1 )       ! t_til = t_til/norm
        call dcopy( N, t_tilde, 1, Vb(1,m), 1 ) ! VBm   = t_til

        ! (8) (9) (10)
        offset = (m-1)*m/2
        call ddmtv( n, m, V, Va(1,m), VAV(offset+1), work )

        ! (11)
        call dcopy( offset+m, VAV, 1, VAVtmp, 1 )
        call compute_eigenpairs( m, VAVtmp, teta, S, work )
        write(23,*) '      CANDIDATE EIGENVALUES:', teta(1:m)
        flush(23)

        ! (12)
        call DGEMV( 'N', n, m, ONE, V, n, S, 1, ZERO, u, 1 )
        call DGEMV( 'N', n, m, ONE, Vb, n, S, 1, ZERO, p, 1 )
        call DGEMV( 'N', n, m, ONE, Va, n, S, 1, ZERO, ua, 1 )
        do i= 1, n
          r(i) = ua(i) - teta(1)*p(i)
        enddo

        ! (13) 
        do while ( residual_norm( n, r, u) .le. eps )
          write(23,*) '      GOOD EIGENVALUE(', k+1, ')=', teta(1)
          flush(23)

          ! (14)
          call dcopy( n, u, 1, Q(1,k+1), 1 )
          call dcopy( n, p, 1, Z(1,k+1), 1 )
          lambda(k+1)   = teta(1)
          k             = k + 1

          ! (15)
          if ( k .eq. kmax ) goto 100

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

            ! (19) (20)
            offset = 0
            do i = 1, m
              do j= 1, i-1
                offset = offset + 1
                VAV(offset) = 0.0
              enddo
              offset = offset + 1
              VAV(offset) = teta(i+1)
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

!            reorthogonalize = reorthogonalize + 1
          endif
          ! (21)
          call dcopy( n, V,  1, u, 1 )
          call dcopy( n, Vb, 1, p, 1 )
          do i= 1, n
            r(i) = Va(i,1) - teta(1)*p(i)
          enddo

#ifdef RGT
          do i= 1, m
            do j= 1, k
              norm = dvdot( n, Z(1,j), V(1,i) )
              write(23,*) i, j, 'V _|_ Z?', norm,
     &        dvdot( n, Q(1,j), Vb(1,i) )
            enddo
            do j = 1, i
              norm = dvdot( n, V(1,i), Vb(1,j) )
              write(23,*) i, j, 'V^t _|_ Vb?', norm,
     &        dvdot( n, V(1,j), Vb(1,i) )
            enddo
          enddo
#endif
        ! (22)  
        enddo


        ! (23)
        if ( m .gt. mmax ) then
#ifdef RGT
          write(23,*) '  Restart Krylov Space'
          flush(23)
#endif
          reorthogonalize = reorthogonalize + 1

          ! (24)
          VAV = 0.0
          ! (25) - (27)
          call DGEMM( 'N', 'N', n, mmin-1, m, ONE, V, n, S(m+1),
     &                m, ZERO, VAUX, n )
          call dcopy( n, u, 1, V, 1 )
          call dcopy( n*(mmin-1), VAUX, 1, V(1,2), 1 )

          call DGEMM( 'N', 'N', n, mmin-1, m, ONE, Va, n, S(m+1),
     &                m, ZERO, VAUX, n )
          call dcopy( n, ua, 1, Va, 1 )
          call dcopy( n*(mmin-1), VAUX, 1, Va(1,2), 1 )

          call DGEMM( 'N', 'N', n, mmin-1, m, ONE, Vb, n, S(m+1),
     &                m, ZERO, VAUX, n )
          call dcopy( n, p, 1, Vb, 1 )
          call dcopy( n*(mmin-1), VAUX, 1, Vb(1,2), 1 )

          ! (28)
          m      = mmin
          offset = 0
          do i = 1, m
            do j= 1, i-1
              offset = offset + 1
              VAV(offset) = 0.0
            enddo
            offset = offset + 1
            VAV(offset) = teta(i+1)
          enddo

#ifdef RGT
          do i= 1, m
            do j= 1, k
              norm = dvdot( n, Z(1,j), V(1,i) )
              write(23,*) i, j, 'V _|_ Z?', norm,
     &        dvdot( n, Q(1,j), Vb(1,i) )
            enddo
            do j = 1, i
              norm = dvdot( n, V(1,i), Vb(1,j) )
              write(23,*) i, j, 'V^t _|_ Vb?', norm,
     &        dvdot( n, V(1,j), Vb(1,i) )
            enddo
          enddo
#endif
        ! (29)
        endif

        if (reorthogonalize.ge.jd_reorthogonalize) then
          nreort = nreort + 1
          write(23,*) '  Reorthogonalization of V & Vb'
          flush(23)
          do i= 1, m
            do j= 1, k
              norm = -dvdot( n, Z(1,j), V(1,i) )
              call DAXPY( n, norm, Q(1,j), 1, V(1,i), 1 )
            enddo

            if (i.gt.1)
     &        call GramSchmidt_B( N, i-1, Vb, V, V(1,i), work,
     &                            work(mmax+1) )

            ! (6) 
            call dsmv( n, nc, B, V(1,i), xx, t_tilde )
            norm = dvdot( n, t_tilde, V(1,i) )
            norm = 1.0/sqrt(norm)

            ! (7)
            call dscal( N, norm, V(1,i), 1 )           ! t     = t/norm
            call dsmv( n, nc, A, V(1,i), xx, Va(1,i) ) ! VAm   = A*Vm 
            call dscal( N, norm, t_tilde, 1 )          ! t_til = t_til/norm
            call dcopy( N, t_tilde, 1, Vb(1,i), 1 )    ! VBm   = t_til

            ! (8) (9) (10)
            offset = (i-1)*i/2
            call ddmtv( n, i, V, Va(1,i), VAV(offset+1), work )
          enddo

#ifdef RGT
          do i= 1, m
            do j= 1, k
              norm = dvdot( n, Z(1,j), V(1,i) )
              write(23,*) i, j, 'RESTART V _|_ Z?', norm,
     &        dvdot( n, Q(1,j), Vb(1,i) )
            enddo
            do j = 1, i
              norm = dvdot( n, V(1,i), Vb(1,j) )
              write(23,*) i, j, 'RESTART V _|_ Vb?', norm,
     &         dvdot( n, V(1,j), Vb(1,i) )
            enddo
          enddo
#endif
          reorthogonalize = 0
        endif
        ! (30)
        rteta = teta(1)
        call dcopy( n, p, 1, Z(1,k+1), 1 )
        call dcopy( n, u, 1, Q(1,k+1), 1 )

        ! (31)
        t  = 0.0
        call medirTiempo( 1, 7 )
        call jd_correction( n, k+1, Z, A, B, rteta, r, gm_it, t )
        call medirTiempo( 2, 7 )

#ifdef RGT
        do i= 1, k+1
          norm = dvdot( n, Z(1,i), t )
          write(23,*) i, k+1, 't _|_ Z?', norm
        enddo
#endif

        gm_it_tot = gm_it_tot + gm_it
        it = it + 1

      ! (32)  
      end do


100   continue

!!!      call MPI_Finalize( it )
      STOP

      call de_alloc( work, 'work', 'jd_solver' )
      call de_alloc( xx, 'xx', 'jd_solver' )
      call de_alloc( Z, 'Z', 'jd_solver' )
      call de_alloc( r, 'r', 'jd_solver' )
      call de_alloc( ua, 'ua', 'jd_solver' )
      call de_alloc( p, 'p', 'jd_solver' )
      call de_alloc( u, 'u', 'jd_solver' )
      call de_alloc( S, 'S', 'jd_solver' )
      call de_alloc( teta, 'teta', 'jd_solver' )
      call de_alloc( VAVtmp, 'VAVtmp', 'jd_solver' )
      call de_alloc( VAV, 'VAV', 'jd_solver' )
      call de_alloc( Vb, 'Vb', 'jd_solver' )
      call de_alloc( Va, 'Va', 'jd_solver' )
      call de_alloc( V, 'V', 'jd_solver' )
      call de_alloc( t_tilde, 't_tilde', 'jd_solver' )
      end subroutine jd_solver

      subroutine compute_eigenpairs( m, VAV, W, VR, work )
      use sys,       only : die
      implicit none
C     Input parameters
      integer,    intent(in) :: m
      real(dp)               :: VAV((m*m+m)/2)
C     Output parameters
      real(dp)               :: W(m), VR(m*m), work(3*M)
C     Local Variables
      integer                :: info

      call medirTiempo( 1, 12 )
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
      call DSPEV( 'V', 'U', m, VAV, W, VR, m, work, info )

      if (info.ne.0) call die( 'ERROR in compute_eigenpairs->DSPEV' )
      call medirTiempo( 2, 12 )

      end subroutine compute_eigenpairs


      subroutine jd_correction( n, ne, Z, A, B, teta, r, gm_it, t )
      use sys,       only : die
      implicit none
C     Input parameters
      integer               :: n, ne
      real(dp)              :: Z(n,ne)
      type(matrix_t)        :: A, B
      real(dp)              :: teta, r(n)
C     Output parameters
      real(dp)              :: t(n)
      integer               :: gm_it
C     Local variables
      integer               :: info
      real(dp),     pointer :: Z_hat(:,:), M(:,:), r_hat(:), gama(:),
     &                         alpha(:), Maux(:,:), work(:)
      integer,      pointer :: ipiv(:)
      real*8,      parameter :: ONE=1.0, ZERO=0.0, M_ONE=-1.0

C     Allocate local memory
      call re_alloc( Z_hat, 1, n, 1, ne, 'Z_hat', 'jd_correction' )
      call re_alloc( M, 1, ne, 1, ne, 'M', 'jd_correction' )
      call re_alloc( ipiv, 1, ne, 'ipiv', 'jd_correction' )
      call re_alloc( r_hat, 1, n, 'r_hat', 'jd_correction' )
      call re_alloc( gama, 1, ne, 'gama', 'jd_correction' )
      call re_alloc( alpha, 1, ne, 'alpha', 'jd_correction' )
      call re_alloc( Maux, 1, ne, 1, ne, 'Maux', 'jd_correction' )
      call re_alloc( work, 1, ne, 'work', 'jd_correction' )

      ! (1)
      call precond_solver_mat( n, ne, PC, Z, Z_hat )

      ! (2)
      call ddmtm( n, ne, Z, Z_hat, M, Maux )

      ! (3)
C     GETRF computes an LU factorization of a general M-by-N matrix A
C     using partial pivoting with row interchanges.
      call DGETRF( ne, ne, M, ne, ipiv, info )

      if ( info .ne. 0 ) call die( 'ERROR in DGETRF' )


      ! (4)
      ! (5)
      call precond_solver_vec( n, PC, r, r_hat )

      ! (6)
      call ddmtv( n, ne, Z, r_hat, gama, work )

      ! (7) (8)
C     ZGETRS solves a system of linear equations
C     A * X = B,  A**T * X = B,  or  A**H * X = B
C     with a general N-by-N matrix A using the LU factorization 
C     computed by DGETRF.
      call DGETRS( 'N', ne, 1, M, ne, ipiv, gama, ne, info )
      if (info.ne.0) call die( 'ERROR in jd_correction->DGETRS' )
      alpha = gama

      ! (9) r_hat = Z*alpha - r_hat
      call DGEMV( 'N', n, ne, ONE, Z_hat, n, alpha, 1, M_ONE, r_hat, 1 )

      ! (10)
      call medirTiempo( 1, 8 )
!      call com_gmres( n, ne, A, B, PC, teta, r_hat, Z, Z_hat,
!     &                M, ipiv, t, gm_it )
      call com_minres( n, ne, A, B, PC, teta, r_hat, Z, Z_hat,
     &                M, ipiv, t, gm_it )
      call medirTiempo( 2, 8 )
!      call com_cg( n, ne, A, B, PC, teta, r_hat, Z, Z_hat,
!     &                M, ipiv, t, gm_it )

C     Free local memory
      call de_alloc( work, 'work', 'jd_correction' )
      call de_alloc( Maux, 'Maux', 'jd_correction' )
      call de_alloc( alpha, 'alpha', 'jd_correction' )
      call de_alloc( gama, 'gama', 'jd_correction' )
      call de_alloc( r_hat, 'r_hat', 'jd_correction' )
      call de_alloc( ipiv, 'ipiv', 'jd_correction' )
      call de_alloc( M, 'M', 'jd_correction' )
      call de_alloc( Z_hat, 'Z_hat', 'jd_correction' )

      end subroutine jd_correction


      end module jacobi_davidson
