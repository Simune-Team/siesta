      module slepc_diagon
#ifdef _SLEPC_
      use precision,    only : dp
      use sys,          only : die
      use alloc,        only : re_alloc, de_alloc
      implicit none

      PUBLIC :: slepc_diag
      PRIVATE
#include "finclude/slepc.h"
#include "finclude/slepceps.h"

C---------------         SLEPC parameters         ---------------
      PetscScalar :: sigma   ! The value of the shift 
      PetscInt    :: nev     ! number of eigenvalues to compute
      PetscInt    :: mpd     ! The maximum dimension allowed for the
                             ! projected problem
      PetscInt    :: ncv     ! The maximum dimension of the subspace to
                             ! be used by the solver
      PetscReal   :: tol     ! The convergence tolerance
      PetscInt    :: maxit   ! Maximum number of iterations to use 
      EPSType     :: epsty   ! Type odf solver: dsitrlanczos
                             !                  lanczos
                             !                  arpack

      contains

      subroutine slepc_diag( nuotot, nuo, nnz, nspin, neig, fixspin,
     &                       temp, qtot, qs, occtol, numh, listh, H, S,
     &                       ef, efs, Entropy, Dnew, Enew )
      use domain_decom, only : slepc_perm
      use fdf
      implicit none
#include "mpif.h"
#include "finclude/petsc.h"
#include "finclude/petscmat.h"
C     Input variables
      integer,   intent(in) :: nuotot, nuo, nspin
      logical,   intent(in) :: fixspin
      real(dp),  intent(in) :: temp, qtot, qs(nspin), occtol
      integer,   intent(in) :: numh(nuo), listh(nnz)
      real(dp),  intent(in) :: H(nnz), S(nnz)
C     Output variables
      integer,  intent(out) :: neig
      real(dp), intent(out) :: ef, efs(nspin), Entropy
      real(dp), intent(out) :: Dnew(nnz), Enew(nnz)
C     Local variables
      real(dp),     pointer :: w(:), Z(:,:), qo(:)

      integer               :: i, ierr, nnz
      integer,      pointer :: col(:)
      Mat                   :: A, B
      logical               :: getD, getPSI

      call write_debug( '    PRE SLEPC_DIAG' )
      call timer( 'slepc_diag', 1 )

!     Number of eigenvalues to compute
      neig = (qtot*6)/5
      neig = min(neig,nuotot)

!     Initialize matrix structures for SLEPC
      nullify( col )
      call re_alloc( col, 1, nnz, 'col', 'slepc_diagon' )
      do i= 1, nnz
        col(i) = slepc_perm(listh(i))-1
      enddo
      call create_slepc_matrix( nuotot, nuo, nnz, numh, col, H, A )
      call create_slepc_matrix( nuotot, nuo, nnz, numh, col, S, B )
      call de_alloc( col, 'col', 'slepc_diagon' )

!     Set slepc parameters
      sigma = -1.0
      nev   = neig
      mpd   = min(200,nuotot)
      ncv   = min(mpd+nev,nuotot)
      tol   = 1.0E-8
      maxit = 10*nuotot
      epsty = 'dsitrlanczos'

      sigma = fdf_double(  'slepc.sigma', sigma )
      nev   = fdf_integer( 'slepc.nev',   nev )
      mpd   = fdf_integer( 'slepc.mpd',   mpd )
      ncv   = fdf_integer( 'slepc.ncv',   ncv )
      maxit = fdf_integer( 'slepc.maxit', maxit )
      tol   = fdf_double(  'slepc.tol',   tol )
      epsty = fdf_string(  'slepc.epsty', epsty )

      nev   = min(nev,nuotot)
      mpd   = min(mpd,nuotot)
      ncv   = min(ncv,nev+mpd)
      neig  = nev

      write(23,*) 'sigma:', sigma
      write(23,*) 'nev  :', nev
      write(23,*) 'mpd  :', mpd
      write(23,*) 'ncv  :', ncv
      write(23,*) 'maxit:', maxit
      write(23,*) 'tol  :', tol
      write(23,*) 'epsty:', epsty
      flush(23)

!     Allocate memory
      nullify( w, Z, qo )
      call re_alloc( w, 1, neig, 'w', 'slepc_diagon' )
      call re_alloc( Z, 1, nuo, 1, neig, 'Z', 'slepc_diagon' )
      call re_alloc( qo, 1, neig, 'qo', 'slepc_diagon' )

      getD   = .true.
      getPSI = .false.

      call slepc_solver( nuo, nnz, nspin, neig, getD, getPSI, fixspin,
     &                   temp, qtot, qs, occtol, numh, listh, A, B, ef,
     &                   efs, Entropy, w, Z, qo, Dnew, Enew )

      call de_alloc( qo, 'qo', 'slepc_diagon' )
      call de_alloc( Z,  'Z',  'slepc_diagon' )
      call de_alloc( w,  'w',  'slepc_diagon' )

      call MatDestroy( B, ierr )
      call MatDestroy( A, ierr )

      call timer( 'slepc_diag', 2 )
      call write_debug( '    POS SLEPC_DIAG' )
      end subroutine slepc_diag

      subroutine slepc_solver( n, nnz, nspin, neig, getD, getPSI,
     &                         fixspin, temp, qtot, qs, occtol, numh,
     &                         listh, A, B, ef, efs, Entropy, w, Z, qo,
     &                         Dnew, Enew )
      use m_fermid,     only : fermid, fermispin
      use domain_decom, only : dd_ncolum, dd_cperm
      implicit none
#include "mpif.h"
#include "finclude/petsc.h"
#include "finclude/petscmat.h"
#include "finclude/slepc.h"
#include "finclude/slepceps.h"
#include "finclude/slepcst.h"
#include "finclude/petscvec.h"
C     Input variables
      integer,        intent(in) :: n, nnz, nspin, neig, numh(n),
     &                              listh(nnz)
      logical,        intent(in) :: getD, getPSI, fixspin
      real(dp)                   :: temp, qtot, qs(nspin), occtol
      Mat                        :: A, B
C     Output variables
      real(dp),      intent(out) :: ef, efs(nspin), Entropy
      real(dp),      intent(out), target :: Z(n,neig)
      real(dp),      intent(out) :: w(neig), qo(neig),
     &                              Dnew(nnz), Enew(nnz)
C     Local variables
      integer,         parameter :: nk = 1
      real(dp),        parameter :: wk(1) = (/ 1.0_dp /)
      logical                    :: notfound
      integer                    :: io, ie, j, jo, k, uneig, nc
      integer,           pointer :: col(:)
      real(dp)                   :: qe, ee, qei, eei
      real(dp),          pointer :: ZZ(:,:)
      PetscInt    :: nconv, i, ierr
      PetscScalar :: kr, ki
      PetscScalar, pointer :: array(:)
      EPS :: SOLVER
      ST  :: SpecT
      KSP :: ksp
      PC  :: precond 
      Vec :: x
      EPSProblemType mytype
      EPSConvergedReason reason
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C     Create the eigensolver and display info
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      call EPSCreate( PETSC_COMM_WORLD, SOLVER, ierr )
      if (ierr.ne.0) call die( 'EPSCreate error' )

      call EPSSetOperators( SOLVER, A, B, ierr )
      if (ierr.ne.0) call die( 'EPSSetOperators error' )

      call EPSSetProblemType( SOLVER, EPS_GHEP, ierr )
      if (ierr.ne.0) call die( 'EPSSetProblemType error' )

      call EPSSetWhichEigenpairs( SOLVER, EPS_SMALLEST_REAL, ierr )
      if (ierr.ne.0) call die( 'EPSSetWhichEigenpairs error' )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C     Set the method and initial shift
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      call EPSSetType( SOLVER, epsty, ierr)
!      call EPSSetType( SOLVER, 'lanczos', ierr)
!      call EPSSetType( SOLVER, 'arpack', ierr)
      if (ierr.ne.0) call die( 'EPSSetType error' )

      call EPSGetST( SOLVER, SpecT, ierr )
      if (ierr.ne.0) call die( 'EPSGetST error' )

      call STSetType( SpecT, STSINV, ierr )
      if (ierr.ne.0) call die( 'STSetType error' )

      call STSetShift( SpecT, sigma, ierr ) 
      if (ierr.ne.0) call die( 'STSetShift error' )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C     Choose a solver
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      call STGetKSP( SpecT, ksp, ierr )
      if (ierr.ne.0) call die( 'STGetKSP error' )

      call KSPGetPC( ksp, precond, ierr )
      if (ierr.ne.0) call die( 'KSPGetPC error' )

      call KSPSetType( ksp, KSPPREONLY, ierr )
      if (ierr.ne.0) call die( 'KSPSetType error' )

      call PCSetType( precond, PCLU, ierr ) 
      if (ierr.ne.0) call die( 'PCSetType error' )

      call PCFactorSetMatSolverPackage( precond, MAT_SOLVER_MUMPS,
     &                                  ierr )
      if (ierr.ne.0) call die( 'PCFactorSetMatSolverPackage error' )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C     Set the parameters and stoping criterium
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      call EPSSetDimensions( SOLVER, nev, ncv, mpd, ierr )
      if (ierr.ne.0) call die( 'EPSSetDimensions error' )

      call EPSSetTolerances( SOLVER, tol, maxit, ierr )
      if (ierr.ne.0) call die( 'EPSSetTolerances error' )

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      call EPSSetFromOptions( SOLVER, ierr )
      if (ierr.ne.0) call die( 'EPSSetFromOptions error' )

!      call EPSGetDimensions( SOLVER, nev, ncv, mpd, ierr )
!      call EPSGetProblemType( SOLVER, mytype, ierr )
!      call EPSGetTolerances( SOLVER, tol, maxit, ierr )
!      call EPSGetType( SOLVER, epsty, ierr )

      call medirTiempo( 1, 8 )
      call timer( 'EPSSolve', 1 )
      write(23,*) 'PRE EPSSolve'
      flush(23)
      call EPSSolve( SOLVER, ierr )
      write(23,*) 'POS EPSSolve'
      flush(23)
      call timer( 'EPSSolve', 2 )
      if (ierr.ne.0) call die( 'EPSSolve error' )
      call medirTiempo( 2, 8 )
      call showTiempo( 23, 8,  'Tiempo EPSSolve:' )

      call EPSGetConvergedReason( SOLVER, reason, ierr )

      select case( reason )
      case(EPS_CONVERGED_TOL)
        write(23,*) 'EPS_CONVERGED_TOL'
      case(EPS_CONVERGED_ITERATING)
        write(23,*) 'EPS_CONVERGED_ITERATING'
      case(EPS_DIVERGED_ITS)
        write(23,*) 'EPS_DIVERGED_ITS'
      case(EPS_DIVERGED_BREAKDOWN)
        write(23,*) 'EPS_DIVERGED_BREAKDOWN'
      case(EPS_DIVERGED_NONSYMMETRIC)
        write(23,*) 'EPS_DIVERGED_NONSYMMETRIC'
      end select
      flush(23)



! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Recover the solution
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      call EPSGetConverged( SOLVER, nconv, ierr )
      write(23,*) 'nconv', nconv, ' of', neig, '  ierr', ierr
      flush(23)
      if (ierr.ne.0 .or. nconv.lt.nev)
     &  call die( 'EPSGetConverged error' )


      call MatGetVecs( A, PETSC_NULL_OBJECT, x, ierr )
      if (ierr.ne.0) call die( 'MatGetVecs error' )

      allocate( array(n) )

      write(23,*) '=== AUTOVALORES de salida ==='
      flush(23)
      do i= 0, nev-1
        call VecPlaceArray( x, array, ierr )
        if (ierr.ne.0) call die( 'VecPlaceArray error' )

        call EPSGetEigenpair( SOLVER, i, kr, ki, x,
     &                        PETSC_NULL_OBJECT, ierr )
        if (ierr.ne.0) call die( 'EPSGetEigenpair error' )

        w(i+1) = PetscRealPart(kr)
        do j= 1, n
          Z(j,i+1) = PetscRealPart(array(j))
        enddo

        call VecResetArray( x, ierr )
      enddo
      do i= 1, min(nev-5,5)
        write(23,*) 'EV(', i+1, '):', W(i+1)   ! , Z(:,i+1)
        flush(23)
      enddo
      do i= max(nev-5,1), nev
        write(23,*) 'EV(', i+1, '):', W(i+1)   ! , Z(:,i+1)
        flush(23)
      enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Free work space
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      deallocate( array )
      call EPSDestroy( SOLVER, ierr )
      write(23,*) ' call EPSDestroy '
      flush(23)

C Check if we are donecol
      if (.not.getD) return

C Find new Fermi energy and occupation weights
      if (fixspin) then
        call fermispin( nspin, nspin, nk, wk, nev, nev, 
     &                  w, temp, qs, qo, efs, Entropy )
      else
        call fermid( nspin, nspin, nk, wk, nev, nev, 
     &               w, temp, qtot, qo, ef, Entropy )
      endif

C     Count the number of usefull eigenvalues/vectors
      uneig    = 1
      notfound = .true.
      do while(notfound)
        if (qo(uneig).gt.occtol) then
          if (uneig+1.gt.nev) then
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
      nc = dd_ncolum
      nullify(ZZ)
      call re_alloc( ZZ, 1, nc, 1, uneig, 'ZZ', 'slepc_diagon' )
      call reduce_matrix( n, nc, uneig, Z, ZZ )
#else
      ZZ => Z
#endif

      call re_alloc( col, 1, nnz, 'col', 'slepc_diagon' )
      do j= 1, nnz
        col(j) = dd_cperm(listh(j))
      enddo

      Dnew = 0.0
      Enew = 0.0
      do ie = 1, uneig
        qe = qo(ie)
        ee = qe*w(ie)
        j = 1
        do io = 1, n
          qei = qe*ZZ(io,ie)
          eei = ee*ZZ(io,ie)
          do k = 1, numh(io)
            jo = col(j)
            Dnew(j) = Dnew(j) + qei*ZZ(jo,ie)
            Enew(j) = Enew(j) + eei*ZZ(jo,ie)
            j = j + 1
          enddo
        enddo
      enddo

      call de_alloc( col, 'col', 'slepc_diagon' )
#ifdef MPI
      call de_alloc( ZZ, 'ZZ', 'slepc_diagon' )
#endif
      end subroutine slepc_solver

      subroutine create_slepc_matrix( nn, nr, nnz, numh, col, H, A )
      implicit none
#include "mpif.h"
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
C     Input variables
      integer,   intent(in) :: nn, nr, nnz
      integer,   intent(in) :: numh(0:nr-1), col(nnz)
      real(dp),  intent(in) :: H(nnz)
C     Output variables
      Mat,      intent(out) :: A

C     Local parameters
!      integer               :: i, j, k, col, Istart, Iend, ierr
      integer                :: i, j, k, maxnz
      PetscInt               :: n
      PetscInt               :: nrows
      PetscInt               :: row
      PetscInt               :: clm
      PetscInt               :: Istart
      PetscInt               :: Iend
      PetscInt               :: ierr
      PetscScalar,   pointer :: val(:)
      PetscInt,      pointer :: ourlens(:)
      PetscInt,      pointer :: offlens(:)

      call timer( 'crea_matri', 1 )
      call medirTiempo( 1, 9 )
      write(23,*) '     Creating SLEPC matrix:', nn, nr, nnz
      flush(23)


      call MatCreate( PETSC_COMM_WORLD, A, ierr )
      if (ierr.ne.0) call die( 'Matcreate error' )

      n     = nn
      nrows = nr
      call MatSetSizes( A, nrows, nrows, n, n, ierr )
      if (ierr.ne.0) call die( 'MatSetSizes error' )

      call MatSetType( A, MATMPIAIJ, ierr )
      if (ierr.ne.0) call die( 'MatSetType error' )

      call MatSetFromOptions( A, ierr )
      if (ierr.ne.0) call die( 'MatSetFromOptions error' )

      call MatGetOwnershipRange( A, Istart, Iend, ierr )
      write(23,*) '     matrix range:', Istart, Iend
      flush(23)

      allocate(ourlens(0:nr-1))
      allocate(offlens(0:nr-1))
      k     = 0
      do i= 0, nr-1
        maxnz = max(maxnz,numh(i))
        offlens(i) = 0
        ourlens(i) = 0
        do j= 1, numh(i)
          k = k + 1
          if (col(k).lt.Istart .or. col(k).ge.Iend) then
            offlens(i) = offlens(i) + 1
          else
            ourlens(i) = ourlens(i) + 1
          endif
        enddo
      enddo

      call MatMPIAIJSetPreallocation( A, 0, ourlens, 0, offlens, ierr )

      deallocate(offlens)
      deallocate(ourlens)

      maxnz = 0
      do i= 0, nr-1
        maxnz = max(maxnz,numh(i))
      enddo

      allocate( val(maxnz) )
      k = 0
      row = Istart
      call medirTiempo( 1, 10 )
      do i= 0, nr-1
        do j= 1, numh(i)
          val(j) = H(j+k)
        enddo
        call MatSetValues( A, 1, row, numh(i), col(k+1), val,
     &                     INSERT_VALUES, ierr )
        k = k + numh(i)
        row = row + 1
      enddo
      call medirTiempo( 2, 10 )
      deallocate( val )

      call medirTiempo( 1, 11 )
      call MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY, ierr )
      call MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY, ierr )
      call medirTiempo( 2, 11 )

      call medirTiempo( 2, 9 )
      call showTiempo( 23, 9,  'Total matrix     :' )
      call showTiempo( 23, 10, '  Insert row      :' )
      call showTiempo( 23, 11, '  Assembly matrix :' )
      call timer( 'crea_matri', 2 )

      end subroutine create_slepc_matrix

#endif    /* _SLEPC_ */

      end module slepc_diagon      
