      MODULE GMRES
      use precision
      use alloc,    only : re_alloc, de_alloc
      use precond,  only : precond_t, precond_solver_vec
      use matrixop
      implicit none

      real(dp)            :: gmres_eps     = 1.0e-8
      integer             :: gmres_dim     = 257
      integer             :: gmres_maxiter = 500
      real(dp), parameter :: gmres_epsmac  = 1.0d-14


      CONTAINS

      subroutine com_gmres( n, ne, A, B, PC, teta, bb, Z, Z_hat,
     &                      M, ipiv, xx, iters )
C     Input variables
      integer,         intent(in) :: n, ne, ipiv(ne)
      type(matrix_t),  intent(in) :: A, B
      type(precond_t), intent(in) :: PC
      real(dp),        intent(in) :: teta, bb(n), Z(n,ne), Z_hat(n,ne),
     &                               M(ne,ne)
C     Output variables
      real(dp),       intent(out) :: xx(n)
      integer,        intent(out) :: iters
C     Local variables
      logical                     :: conv, FIN2
      integer                     :: ierr, ii, jj, jj1, kk, idx, dim2
      real(dp)                    :: raux, gamma, stopcri
      real(dp),           pointer :: kryl(:,:), rs(:), hh(:),
     &                               cc(:), ss(:)

      ierr  = 0
      iters = 0
      conv  = .FALSE.
      xx    = 0.0

!     raux = ||bb||
      raux = SQRT(dvdot( n, bb, bb ))

      if (raux .eq. 0.0d0) return

      stopcri = gmres_eps * raux
      dim2    = gmres_dim*(gmres_dim+3)/2
      call re_alloc( kryl, 1, n, 1, gmres_dim+1, 'kryl', 'com_gmres' )
      call re_alloc( rs, 1, gmres_dim+1, 'rs', 'com_gmres' )
      call re_alloc( hh, 1, dim2, 'hh', 'com_gmres' )
      call re_alloc( cc, 1, gmres_dim, 'cc', 'com_gmres' )
      call re_alloc( ss, 1, gmres_dim, 'ss', 'com_gmres' )
!------------------------------------------------- MAIN LOOP
      DO WHILE (.NOT. conv)
        call gmres_MV( n, ne, A, B, PC, teta, xx, Z, Z_hat, M, ipiv,
     &                 kryl(:,1) )

        kryl(1:n,1) = bb - kryl(1:n,1)
C       raux = ||kryl(*,1)||
        raux = SQRT(dvdot( n, kryl(1:n,1), kryl(1:n,1)))

        if (raux.le.stopcri) then
C         The initial guess is the solution
          conv = .TRUE.
        else
!         Initialize 1-st term of the rhs of hessenberg system
          rs(1) = raux
!         Ortonormalize kryl(*,1)
          raux = 1.0/raux
          kryl(1:n,1) = kryl(1:n,1)*raux
        endif

        jj   = 0
        idx  = -1
        FIN2 = conv

C       Inner loop. Restarted each kryldim iterations.
        DO WHILE (.NOT. FIN2)
          iters = iters + 1
          jj    = jj + 1
          jj1   = jj + 1
          idx   = idx + jj

          call gmres_MV( n, ne, A, B, PC, teta, kryl(:,jj), Z, Z_hat,
     &                   M, ipiv, kryl(:,jj1) )
C         Modified Gram-Schmidt
C         For i= 1, j
C           H(i,j) = <v_i, v_j1>
C           v_j1   = v_j1 - H(i,j) * v_i
          DO ii= 1, jj
            raux = dvdot( n, kryl(1:n,ii), kryl(1:n,jj1))
            hh(idx+ii) = raux
            kryl(1:n,jj1) = kryl(1:n,jj1) - raux*kryl(1:n,ii)
          ENDDO

C         H(jj1,jj) = ||kryl(*,jj1)||
          raux = SQRT(dvdot( n, kryl(1:n,jj1), kryl(1:n,jj1)))
          hh(idx+jj1) = raux

          if (raux.eq.0.0d0) then
            ierr  = 2
            FIN2  = .TRUE.
            conv  = .TRUE.
            idx   = idx - jj
            jj    = jj - 1
            iters = iters - 1
          else
C           Ortonormalize kryl(*,jj1)
            raux = 1.0/raux
            kryl(1:n,jj1) = kryl(1:n,jj1)*raux

C           Update factorization of H. Perform previous 
C           transformations on jj-th column of H
            do ii= 1, jj-1
              kk   = ii + 1
              raux = hh(idx+ii)

              hh(idx+ii) =  cc(ii)*raux + ss(ii)*hh(idx+kk)
              hh(idx+kk) = -ss(ii)*raux + cc(ii)*hh(idx+kk)
            enddo

            gamma = hh(idx+jj)*hh(idx+jj) + hh(idx+jj1)*hh(idx+jj1)
            gamma = SQRT(gamma)
     
C           if gamma is zero then take any small
C           value will affect only residual estimate
            if (gamma.eq.0.0) gamma = gmres_epsmac

C           Get next plane rotation
            cc(jj) = hh(idx+jj)/gamma
            ss(jj) = hh(idx+jj1)/gamma

            hh(idx+jj) = cc(jj)*hh(idx+jj) + ss(jj)*hh(idx+jj1)

C           Update the rhs of the LS problem
            rs(jj1) = -ss(jj)*rs(jj)
            rs(jj)  =  cc(jj)*rs(jj)

C           Convergence Test
            raux = ABS( rs(jj1) )

            if (raux.le.stopcri) then
              write(23,*) '         GMRES GOOD...', iters, raux, stopcri
              call pxfflush(23)
                conv = .TRUE.
                FIN2 = .TRUE.
            else
              if (iters.ge.gmres_maxiter) then
              write(23,*) '         GMRES MAXITER', iters
              call pxfflush(23)
                conv = .TRUE.
                FIN2 = .TRUE.
                ierr = 1
              else if (jj .ge. gmres_dim) then
                FIN2 = .TRUE.
              endif
            endif
          endif
        END DO   ! Krylov Loop End
        if (jj .gt. 0) then
C         Compute y => Solve upper triangular system
          do ii= jj, 2, -1
            rs(ii) = rs(ii)/hh(idx+ii)
            raux   = rs(ii)

            do kk= 1, ii-1
              rs(kk) = rs(kk) - hh(idx+kk) * raux
            enddo

            idx = idx - ii
          enddo
          rs(1) = rs(1) / hh(1)

C         Linear combination of kryl(*,jj)'s to get the solution.
          do ii= 1, jj
            raux = rs(ii)
            xx = xx + raux*kryl(1:n,ii)
          enddo
        endif
      END DO 

      call de_alloc( ss, 'ss', 'com_gmres' )
      call de_alloc( cc, 'cc', 'com_gmres' )
      call de_alloc( hh, 'hh', 'com_gmres' )
      call de_alloc( rs, 'rs', 'com_gmres' )
      call de_alloc( kryl, 'kryl', 'com_gmres' )
      end subroutine com_gmres

      subroutine gmres_mv( n, ne, A, B, PC, teta, v, Z, Z_hat,
     &                     M, ipiv, zz )
      use sys,       only : die
      implicit none
C     Input variables
      integer,         intent(in) :: n, ne, ipiv(ne)
      type(matrix_t),  intent(in) :: A, B
      type(precond_t), intent(in) :: PC
      real(dp),        intent(in) :: teta, v(n), Z(n,ne), Z_hat(n,ne),
     &                               M(ne,ne)
C     Output variables
      real(dp),       intent(out) :: zz(n)
C     Local variables
      integer                     :: info
      real(dp),           pointer :: y(:), y_hat(:), gama(:)

      call re_alloc( y, 1, n, 'y', 'gmres_mv' )
      call re_alloc( y_hat, 1, n, 'y_hat', 'gmres_mv' )
      call re_alloc( gama, 1, ne, 'gama', 'gmres_mv' )

      ! (11)
      call dsmlamv( n, A, teta, B, v, y )

      ! (12)
      call precond_solver_vec( n, PC, y, y_hat )

      ! (13)
      call ddmtv( n, ne, Z, y, gama )

      ! (14) (15)
C     DGETRS solves a system of linear equations
C         A * X = B,  A**T * X = B,  or  A**H * X = B
C     with a general N-by-N matrix A using the LU factorization
C     computed by DGETRF.
      call DGETRS ( 'N', ne, 1, M, ne, ipiv, gama, ne, info )
      if (info.ne.0) call die( 'ERROR in gmres_mv->DGETRS' )

      ! (16)
      call ddvlmv( n, ne, y_hat, Z_hat, gama, zz )

      call de_alloc( gama, 'gama', 'gmres_mv' )
      call de_alloc( y_hat, 'y_hat', 'gmres_mv' )
      call de_alloc( y, 'y', 'gmres_mv' )

      end subroutine gmres_mv

      END MODULE GMRES
