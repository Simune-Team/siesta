module m_minim

use densematrix,  only : psi
use parallel,     only : BlockSize, Node, Nodes
use parallelsubs, only : GetNodeOrbs, GlobalToLocalOrb, WhichNodeOrb
use precision,    only : dp
use sys,          only : die
use m_timer,      only : timer_start, timer_stop
#ifdef MPI
use mpi_siesta,   only : mpi_integer, mpi_double_precision, mpi_comm_world, mpi_double_precision, mpi_sum, mpi_status_size
#endif

implicit none

!**** PRIVATE ***********************************!

private

real(dp), parameter :: Pi=3.141592653589793238462643383279502884197_dp

logical, save :: FirstCall=.true.      ! First time minim_cg is called?
logical, save :: LineSearchFit=.false. ! Use fitting method for line search?

integer, save :: N_occ_loc  ! Num. of (L)WFs (local)
integer, save :: desc1(1:9) ! Descriptor for operator matrices in orbital basis
integer, save :: desc2(1:9) ! Descriptor for operator matrices in (L)WF basis
integer, save :: desc3(1:9) ! Descriptor for (L)WF coeffs. matrix

!**** PUBLIC ************************************!

public :: minim

real(dp), allocatable, public, save :: chi(:,:) ! (L)WF coeffs. matrix

#ifdef MPI
integer, public, save :: BlockSize_chi ! ScaLAPACK blocking factor for distributing the (L)WF coeffs. matrix
#endif

!************************************************!

contains

subroutine minim(calc_Escf,iscf,nbasis,nspin,nuotot,nhmax,numh,listhptr,listh,Dscf,eta,qs,h,s)
  implicit none

  !**** INPUT ***********************************!

  logical, intent(in) :: calc_Escf ! Calculate the energy-density matrix from the existing coeffs.?

  integer, intent(in) :: iscf            ! SCF iteration num.
  integer, intent(in) :: nbasis          ! Num. of basis orbitals
  integer, intent(in) :: nspin           ! Num. of spins
  integer, intent(in) :: nuotot          ! Num. of orbitals in unit cell (global)
  integer, intent(in) :: nhmax           ! First dimension of listh and H
  integer, intent(in) :: numh(nbasis)    ! Num. of nonzero elements of each row of H
  integer, intent(in) :: listhptr(nhmax) ! Pointer to start of row in listh
  integer, intent(in) :: listh(nbasis)   ! List of nonzero elements of each row of H

  real(dp), intent(in) :: qs(2)                    ! Num. of electrons per spin
  real(dp), intent(in) :: eta(2)                   ! Chemical potential for Kim functional
  real(dp), intent(in), optional :: h(nhmax,nspin) ! Hamiltonian matrix (sparse)
  real(dp), intent(in), optional :: s(nhmax)       ! Overlap matrix (sparse)

  !**** OUTPUT **********************************!

  real(dp), intent(out) :: Dscf(nhmax,nspin) ! Density (or energy-density) matrix (sparse)

  !**** LOCAL ***********************************!

  integer :: nuo ! Num. of orbitals in unit cell (local)

  integer :: ispin, io, jo, j, k, ind, lwork

  real(dp), allocatable :: Haux(:,:) ! Hamiltonian matrix (dense)
  real(dp), allocatable :: Saux(:,:) ! Overlap matrix (dense)
  real(dp), allocatable :: Daux(:,:) ! Density (or energy-density) matrix (dense)

  !**********************************************!

  call timer_start('minim_mod')

  if (nspin/=1) call die('ERROR: only spin-unpolarized calculations allowed with minimization (for now)!')

#ifdef MPI
  call GetNodeOrbs(nuotot,Node,Nodes,nuo)
#else
  nuo=nuotot
#endif

  allocate(Daux(nuotot,nuo))

  if (.not. calc_Escf) then

    allocate(Haux(nuotot,nuo))
    allocate(Saux(nuotot,nuo))

    ! Convert the Hamiltonian and overlap matrices from sparse to dense
    do ispin=1,nspin
      Haux=0.0d0
      Saux=0.0d0
      do io=1,nuo
        do j=1,numh(io)
          ind=listhptr(io)+j
          jo=listh(ind)
          Haux(jo,io)=Haux(jo,io)+H(ind,ispin)
          Saux(jo,io)=Saux(jo,io)+S(ind)
        end do
      end do
    end do

    ! Shift the eingevalue spectrum w.r.t. the chemical potential reference
    Haux=Haux-eta(1)*Saux

  end if

  ! Call the routine to perform the energy minimization
  call minim_cg(calc_Escf,iscf,nuo,nuotot,nint(0.5_dp*qs(1)),eta(1),Daux,psi,nspin,Haux,Saux)

  if (.not. calc_Escf) then
    deallocate(Saux)
    deallocate(Haux)
  end if

  ! Convert the density (or energy-density) matrix from dense to sparse
  Dscf=0.0_dp
  do ispin=1,nspin
    do io=1,nuo
      do j=1,numh(io)
        ind=listhptr(io)+j
        jo=listh(ind)
        Dscf(ind,ispin)=Daux(jo,io)
      end do
    end do
  end do

  deallocate(Daux)

  call timer_stop('minim_mod')

  end subroutine minim

!================================================!
! Minimize the Kim functional by conjugate       !
! gradients                                      !
!================================================!
subroutine minim_cg(calc_Escf,iscf,Hp_dim_loc,Hp_dim,N_occ,eta,Daux,psi,nspin,Hp,Sp)
  implicit none

  !**** INPUT ***********************************!

  logical, intent(in) :: calc_Escf ! Calculate the energy-density matrix from the existing coeffs.?

  integer, intent(in) :: iscf       ! SCF iteration num.
  integer, intent(in) :: Hp_dim_loc ! Num. of orbitals (local)
  integer, intent(in) :: Hp_dim     ! Num. of orbitals (global)
  integer, intent(in) :: N_occ      ! Num. of (L)WFs (global)
  integer, intent(in) :: nspin      ! Num. of spins

  real(dp), intent(in) :: eta                                ! Chemical potential for Kim functional
  real(dp), intent(in) :: psi(1:Hp_dim,1:Hp_dim_loc,1:nspin) ! Eigenvectors from diagonalization
  real(dp), intent(in), optional :: Hp(:,:)                  ! Hamiltonian matrix in orbital basis
  real(dp), intent(in), optional :: Sp(:,:)                  ! Overlap matrix in orbital basis

  !**** OUTPUT **********************************!

  real(dp), intent(out) :: Daux(:,:) ! Density (or energy-denity) matrix in orbital basis

  !**** LOCAL ***********************************!

  logical :: conv, ls_conv

  integer :: n_step_max ! Max. steps of cg minimization
  integer :: icg, info
  integer :: i, j, k, l, m
  integer, save :: ictxt
#ifdef MPI
  integer :: N_occ_loc_i, mpi_status(mpi_status_size)
#endif

  real(dp) :: cg_tol     ! Convergence tolerance of cg minimization
  real(dp) :: step_init  ! Initial trial step length for line search
  real(dp) :: coeff(0:4) ! Coeffs. of the quartic equation
  real(dp) :: rn, rn2, E_omg, E_omg_old
  real(dp) :: lambda, lambda_n, lambda_d, E_diff, TrS
  real(dp), allocatable, save :: H(:,:)
  real(dp), allocatable, save :: S(:,:)
  real(dp), allocatable :: Hd(:,:)
  real(dp), allocatable :: Sd(:,:)
  real(dp), allocatable :: Hdd(:,:)
  real(dp), allocatable :: Sdd(:,:)
  real(dp), allocatable :: g(:,:)
  real(dp), allocatable :: g_p(:,:)
  real(dp), allocatable :: d(:,:)
  real(dp), allocatable :: work1(:,:)
  real(dp), allocatable :: work2(:,:)
  real(dp), allocatable, save :: work3(:,:)
  real(dp), allocatable :: work4(:,:)
  real(dp), allocatable :: work5(:,:)
  real(dp), external :: ddot
#ifdef MPI
  real(dp) :: lambda_n_tot, lambda_d_tot
  real(dp), external :: pdlatra
#endif

  !**********************************************!

  call timer_start('minim_cg')

  ! If this is the first time the minimization module is called, several things need to be done
  ! (detailed below)  
  if (FirstCall) then

#ifdef MPI
    if (N_occ<Nodes) call die('ERROR: Number of nodes is larger than number of occupied states!')

    ! Calculate the ScaLAPACK blocking factor for distributing the (L)WF coeffs. matrix
    ! Probably a good idea in future to let the user override this
    BlockSize_chi=int(N_occ/Nodes)
    if (BlockSize_chi==0) BlockSize_chi=1
    k=N_occ
    N_occ_loc=0
    conv=.false.
    do i=1,100
      do j=0,Nodes-1
        l=BlockSize_chi
        if (k<BlockSize_chi) l=k
        k=k-l
        if (Node==j) N_occ_loc=N_occ_loc+l
        if (k==0) then
          conv=.true.
          exit
        end if
      end do
      if (conv) exit
    end do
    !write(1200+Node,*) Node, Nodes, BlockSize_chi, N_occ, N_occ_loc, Hp_dim, Hp_dim_loc

    ! Initialize the BLACS process grid
    call blacs_get(0,0,ictxt)
    call blacs_gridinit(ictxt,'C',1,Nodes)

    ! Initialize the matrix descriptors
    call descinit(desc1,Hp_dim,Hp_dim,BlockSize,BlockSize,0,0,ictxt,Hp_dim,info)
    if (info/=0) call die('ERROR: desc1 setup has failed in minim!')
    call descinit(desc2,N_occ,N_occ,BlockSize_chi,BlockSize_chi,0,0,ictxt,N_occ,info)
    if (info/=0) call die('ERROR: desc2 setup has failed in minim!')
    call descinit(desc3,Hp_dim,N_occ,BlockSize,BlockSize_chi,0,0,ictxt,Hp_dim,info)
    if (info/=0) call die('ERROR: desc3 setup has failed in minim!')
#else
    N_occ_loc=N_occ
#endif

    allocate(chi(1:Hp_dim,1:N_occ_loc))

    ! If this is the first SCF step, then we need to initialize the (L)WF coeffs. matrix with random
    ! numbers between -0.5 and 0.5 (normalize at the end to avoid instabilities)
    if (iscf==1) then
      call rand_init
      do i=1,N_occ_loc
        do j=1,Hp_dim
          call random_number(rn)
          call random_number(rn2)
          chi(j,i)=sign(rn,rn2-0.5_dp)
        end do
      end do
      chi=1.0d-2*chi/sqrt(real(Hp_dim,dp))
    ! If this is *not* the first SCF step, but it *is* the first time the minimization routine is
    ! called, then we have already calculated the eigenfuctions from a diagonaliazation call--in
    ! this case, we take the lowest N_occ eigenfunctions as our initial guess, but we must take care
    ! to distribute them properly amongst the MPI processes for parallel runs
    else
#ifdef MPI
      ! i: Receiving node
      ! j: Local orbital num. on receiving node
      ! k: Global orbital num.
      ! l: Sending node
      ! m: Local orbital num. on sending node
      k=0
      do i=0,Nodes-1
        N_occ_loc_i=N_occ_loc
        call mpi_bcast(N_occ_loc_i,1,mpi_integer,i,mpi_comm_world,info)
        do j=1,N_occ_loc_i
          k=k+1
          call WhichNodeOrb(k,Nodes,l)
          call GlobalToLocalOrb(k,l,Nodes,m)
          if (Node==l) then
            if (Node==i) then
              chi(1:Hp_dim,j)=psi(1:Hp_dim,m,1)
            else
              call mpi_send(psi(1,m,1),Hp_dim,mpi_double_precision,i,k,mpi_comm_world,info)
            end if
          else if (Node==i) then
            call mpi_recv(chi(1,j),Hp_dim,mpi_double_precision,l,k,mpi_comm_world,mpi_status,info)
          end if
        end do
      end do
#else
      chi=psi(1:Hp_dim,1:N_occ_loc,1)
#endif
    end if

    FirstCall=.false.
  end if

  if (calc_Escf) then

    allocate(work1(1:N_occ,1:N_occ_loc))
    allocate(work2(1:Hp_dim,1:N_occ_loc))

    ! Calculate the energy-density matrix: E=C*[(2*I-S)*(H+eta*S)]*C^T
#ifdef MPI
    call pdsymm('L','U',N_occ,N_occ,1.0_dp,work3,1,1,desc2,H+eta*S,1,1,desc2,0.0_dp,work1,1,1,desc2)
#else
    call dsymm('L','U',N_occ,N_occ,1.0_dp,work3,N_occ,H+eta*S,N_occ,0.0_dp,work1,N_occ)
#endif
    call calc_densmat(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,work1,chi,Daux,work2)

    deallocate(work2)
    deallocate(work1)
    deallocate(work3)
    deallocate(S)
    deallocate(H)

    call timer_stop('minim_cg')

    return

  end if

  if (.not. allocated(H)) allocate(H(1:N_occ,1:N_occ_loc))
  if (.not. allocated(S)) allocate(S(1:N_occ,1:N_occ_loc))
  if (.not. allocated(work3)) allocate(work3(1:N_occ,1:N_occ_loc))

  if (.not. LineSearchFit) then
    allocate(Hd(1:N_occ,1:N_occ_loc))
    allocate(Sd(1:N_occ,1:N_occ_loc))
    allocate(Hdd(1:N_occ,1:N_occ_loc))
    allocate(Sdd(1:N_occ,1:N_occ_loc))
  end if
  allocate(g(1:Hp_dim,1:N_occ_loc))
  allocate(g_p(1:Hp_dim,1:N_occ_loc))
  allocate(d(1:Hp_dim,1:N_occ_loc))
  allocate(work1(1:Hp_dim,1:N_occ_loc))
  allocate(work2(1:Hp_dim,1:N_occ_loc))
  if (LineSearchFit) then
    allocate(work4(1:Hp_dim,1:N_occ_loc))
    allocate(work5(1:Hp_dim,1:N_occ_loc))
  end if

  n_step_max=10000
  cg_tol=1.0d-8
  if (LineSearchFit) step_init=1.0d-1

  ! First we calculate the energy and gradient for our initial guess, with the following steps:
  ! -calculate the Hamiltonian in (L)WFs basis: H=C^T*h*C
  call calc_A(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Hp,chi,H,work1)
  ! -calculate the overlap matrix in (L)WFs basis: S=C^T*s*C
  call calc_A(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Sp,chi,S,work2)
  ! -calculate the energy if we are using the fitting procedure: E=2*Tr(H)-Tr(H*S)
  if (LineSearchFit) then
#ifdef MPI
    call calc_E(N_occ_loc,N_occ,H,S,E_omg,work3)
#else
    call calc_E(N_occ_loc,N_occ,H,S,E_omg)
#endif
  end if
  ! -calculate the gradient: G=2*(2*h*C-s*C*H-h*C*S)
  !  (note that we *reuse* h*C and s*C contained in work1 and work2 from the previous calls to
  !   calc_A)
  call calc_grad(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,H,S,g,work1,work2)
  ! -if we are not using the fitting procedure, calculate the additional matrices:
  !  Hd=G^T*h*C
  !  Sd=G^T*s*C
  !  Hdd=G^T*h*G
  !  Sdd=G^T*s*G
  !  (again, h*C and s*C have already been calculated, although h*G and s*G have not)
  !  and, finally, the coeffs. of the quartic line search equation in the direction G
  !  (the energy at C is given by the zeroth-order coeff. c(0))
  if (.not. LineSearchFit) then
#ifdef MPI
    call pdgemm('T','N',N_occ,N_occ,Hp_dim,1.0_dp,g,1,1,desc3,work1,1,1,desc3,0.0_dp,Hd,1,1,desc2)
    call pdgemm('T','N',N_occ,N_occ,Hp_dim,1.0_dp,g,1,1,desc3,work2,1,1,desc3,0.0_dp,Sd,1,1,desc2)
#else
    call dgemm('T','N',N_occ,N_occ,Hp_dim,1.0_dp,g,Hp_dim,work1,Hp_dim,0.0_dp,Hd,N_occ)
    call dgemm('T','N',N_occ,N_occ,Hp_dim,1.0_dp,g,Hp_dim,work2,Hp_dim,0.0_dp,Sd,N_occ)
#endif
    call calc_A(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Hp,g,Hdd,work1)
    call calc_A(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Sp,g,Sdd,work2)
    call calc_coeff(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,H,S,Hd,Sd,Hdd,Sdd,coeff,work3)
    E_omg=coeff(0)
  end if

  ! This is the main loop of the cg algorithm. We perform a series of line minimizations, with the
  ! gradient g at each new step being modified to obtain the search direction d
  conv=.false.
  d=0.0_dp
  icg=0
  do i=1,n_step_max
    lambda=0.0_dp
    do j=1,Hp_dim*N_occ-1
      d=g+lambda*d
      g_p=g
      E_omg_old=E_omg
      ! Call the routine to perform the line search (either fitting or analytical)
      if (LineSearchFit) then
        call line_search_fit(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Hp,Sp,H,S,d,E_omg,g,ls_conv,step_init,work1,work2,work3,work4,work5,&
                             coeff)
      else
        call line_search_exact(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Hp,Sp,H,S,Hd,Sd,Hdd,Sdd,d,E_omg,work1,work2,work3,coeff,icg)
        ls_conv=.true.
      end if
      icg=icg+1
      E_diff=2.0_dp*abs((E_omg-E_omg_old)/(E_omg+E_omg_old))
      if (Node==0) print*, i, j, E_omg, E_diff
      if (E_diff<=cg_tol) then
        conv=.true.
        exit
      end if
      if (.not. LineSearchFit) call calc_grad(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,H,S,g,work1,work2)
      if (ls_conv) then
        lambda_n=ddot(Hp_dim*N_occ_loc,g,1,g-g_p,1)
        lambda_d=ddot(Hp_dim*N_occ_loc,g_p,1,g_p,1)
#ifdef MPI
        call mpi_allreduce(lambda_n,lambda_n_tot,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
        call mpi_allreduce(lambda_d,lambda_d_tot,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
        lambda=lambda_n_tot/lambda_d_tot
#else
        lambda=lambda_n/lambda_d
#endif
      else
        exit
      end if
    end do
    if (conv) exit
  end do
  if (i>n_step_max) then
    if (Node==0) print*, '#WARNING: Geometry optimization failed to converge!'
  end if

  ! Calculate the density matrix: D=C*(2*I-S)*C^T
#ifdef MPI
  call pdlaset('A',N_occ,N_occ,0.0_dp,2.0_dp,work3,1,1,desc2)
  work3=work3-S
#else
  work3=-S
  do i=1,N_occ
    work3(i,i)=2.0_dp+work3(i,i)
  end do
#endif
  call calc_densmat(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,work3,chi,Daux,work2)

  ! Calculate the trace of S to make sure we are occupying the right number of eigenstates in our
  ! solution
#ifdef MPI
  TrS=pdlatra(N_occ,S,1,1,desc2)
#else
  TrS=0.0_dp
  do i=1,N_occ
    TrS=TrS+S(i,i)
  end do
#endif
  if (Node==0) then
    print*, 'minim: icg              = ', icg
    print*, 'minim: Tr(S)            = ', TrS
  end if

  if (LineSearchFit) then
    deallocate(work5)
    deallocate(work4)
  end if
  deallocate(work2)
  deallocate(work1)
  deallocate(d)
  deallocate(g_p)
  deallocate(g)
  if (.not. LineSearchFit) then
    deallocate(Sdd)
    deallocate(Hdd)
    deallocate(Sd)
    deallocate(Hd)
  end if

  call timer_stop('minim_cg')

end subroutine minim_cg

!================================================!
! Perform line search by exact calculation       !
!================================================!
subroutine line_search_exact(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Hp,Sp,H,S,Hd,Sd,Hdd,Sdd,d,E_omg,work1,work2,work3,coeff,icg)
  implicit none

  !**** INPUT ***********************************!

  integer, intent(in) :: icg        ! CG iteration num.
  integer, intent(in) :: Hp_dim_loc ! Num. of orbitals (local)
  integer, intent(in) :: Hp_dim     ! Num. of orbitals (global)
  integer, intent(in) :: N_occ_loc  ! Num. of (L)WFs (local)
  integer, intent(in) :: N_occ      ! Num. of (L)WFs (global)

  real(dp), intent(in) :: Hp(:,:) ! Hamiltonian matrix in orbital basis
  real(dp), intent(in) :: Sp(:,:) ! Overlap matrix in orbital basis
  real(dp), intent(in) :: d(:,:)  ! Gradient of Kim functional along line

  !**** INOUT ***********************************!

  real(dp), intent(inout) :: E_omg      ! Energy of Kim functional
  real(dp), intent(inout) :: coeff(0:4) ! Coeffs. of the quartic equation
  real(dp), intent(inout) :: H(:,:)     ! Hamiltonian matrix in (L)WF basis
  real(dp), intent(inout) :: S(:,:)     ! Overlap matrix in (L)WF basis
  real(dp), intent(inout) :: Hd(:,:)
  real(dp), intent(inout) :: Sd(:,:)
  real(dp), intent(inout) :: Hdd(:,:)
  real(dp), intent(inout) :: Sdd(:,:)
  real(dp), intent(inout) :: work1(:,:)
  real(dp), intent(inout) :: work2(:,:)
  real(dp), intent(inout) :: work3(:,:)

  !**** LOCAL ***********************************!

  logical :: fail ! Did we fail to find a minimum?

  real(dp) :: x_min ! Position of minimum

  !**********************************************!

  call timer_start('minim_line_search_exact')

  ! If this is not the first cg step, we have to recalculate Hd, Sd, Hdd, Sdd, and the coeffs.
  if (icg>0) then
#ifdef MPI
    call pdgemm('T','N',N_occ,N_occ,Hp_dim,1.0_dp,d,1,1,desc3,work1,1,1,desc3,0.0_dp,Hd,1,1,desc2)
    call pdgemm('T','N',N_occ,N_occ,Hp_dim,1.0_dp,d,1,1,desc3,work2,1,1,desc3,0.0_dp,Sd,1,1,desc2)
#else
    call dgemm('T','N',N_occ,N_occ,Hp_dim,1.0_dp,d,Hp_dim,work1,Hp_dim,0.0_dp,Hd,N_occ)
    call dgemm('T','N',N_occ,N_occ,Hp_dim,1.0_dp,d,Hp_dim,work2,Hp_dim,0.0_dp,Sd,N_occ)
#endif
    call calc_A(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Hp,d,Hdd,work1)
    call calc_A(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Sp,d,Sdd,work2)
    call calc_coeff(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,H,S,Hd,Sd,Hdd,Sdd,coeff,work3)
  end if

  ! Using the coeffs. calculated anlytically, we can find the minimum of the functional in the
  ! search direction, and calculate the energy at that minimum
  call solve_quartic(coeff(0:4),x_min,fail)
  E_omg=coeff(4)*x_min**4+&
        coeff(3)*x_min**3+&
        coeff(2)*x_min**2+&
        coeff(1)*x_min+&
        coeff(0)
  ! In certain regions of the coeffs. space the line search gives no minimum. If the initial guess
  ! for the coeffs. is not too large, this should never occur in the Ordejon-Mauri functional.
  if (fail) call die('ERROR: no minimum in line search!')

  ! Move to the minimum and recalculate H and S at this point
  chi=chi+x_min*d
  call calc_A(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Hp,chi,H,work1)
  call calc_A(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Sp,chi,S,work2)

  call timer_stop('minim_line_search_exact')

end subroutine line_search_exact

!================================================!
! Perform line search by quartic fitting         !
!================================================!
subroutine line_search_fit(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Hp,Sp,H,S,d,E_omg,g,ls_conv,step_init,work1,work2,work3,chi_step,&
                           g_step,coeff)
  implicit none

  !**** INPUT ***********************************!

  integer, intent(in) :: Hp_dim_loc ! Num. of orbitals (local)
  integer, intent(in) :: Hp_dim     ! Num. of orbitals (global)
  integer, intent(in) :: N_occ_loc  ! Num. of (L)WFs (local)
  integer, intent(in) :: N_occ      ! Num. of (L)WFs (global)

  real(dp), intent(in) :: Hp(:,:) ! Hamiltonian matrix in orbital basis
  real(dp), intent(in) :: Sp(:,:) ! Overlap matrix in orbital basis
  real(dp), intent(in) :: d(:,:)  ! Gradient of Kim functional along line

  !**** OUTPUT **********************************!

  logical, intent(out) :: ls_conv ! Was the line search successful?

  !**** INOUT ***********************************!

  real(dp), intent(inout) :: step_init     ! Initial step length
  real(dp), intent(inout) :: coeff(0:4)    ! Coeffs. of the quartic equation
  real(dp), intent(inout) :: E_omg         ! Energy of Kim functional
  real(dp), intent(inout) :: H(:,:)        ! Hamiltonian matrix in (L)WF basis
  real(dp), intent(inout) :: S(:,:)        ! Overlap matrix in (L)WF basis
  real(dp), intent(inout) :: g(:,:)        ! Gradient of Kim functional
  real(dp), intent(inout) :: work1(:,:)
  real(dp), intent(inout) :: work2(:,:)
  real(dp), intent(inout) :: work3(:,:)
  real(dp), intent(inout) :: chi_step(:,:)
  real(dp), intent(inout) :: g_step(:,:)

  !**** LOCAL ***********************************!

  logical :: fail

  integer :: ls_max, info
  integer :: i

  real(dp) :: ls_tol
  real(dp) :: ls(1:3,1:4)
  real(dp), external :: ddot
#ifdef MPI
  real(dp) :: ddot_val
#endif

  !******************************************************************************!

  call timer_start('minim_line_search_fit')

  ls_max=10
  ls_tol=1.0d-0

  ! first point
  ls(1,1)=0.0_dp
  ls(2,1)=E_omg
#ifdef MPI
  ddot_val=ddot(Hp_dim*N_occ_loc,g,1,d,1)
  call mpi_allreduce(ddot_val,ls(3,1),1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
#else
  ls(3,1)=ddot(Hp_dim*N_occ_loc,g,1,d,1)
#endif

  ! second point
  ls(1,2)=-step_init
  chi_step=chi+ls(1,2)*d
  call calc_A(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Hp,chi_step,H,work1)
  call calc_A(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Sp,chi_step,S,work2)
#ifdef MPI
  call calc_E(N_occ_loc,N_occ,H,S,ls(2,2),work3)
#else
  call calc_E(N_occ_loc,N_occ,H,S,ls(2,2))
#endif
  call calc_grad(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,H,S,g_step,work1,work2)
#ifdef MPI
  ddot_val=ddot(Hp_dim*N_occ_loc,g_step,1,d,1)
  call mpi_allreduce(ddot_val,ls(3,2),1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
#else
  ls(3,2)=ddot(Hp_dim*N_occ_loc,g_step,1,d,1)
#endif

  ! third point
  ls(1,3)=-2.0_dp*step_init
  chi_step=chi+ls(1,3)*d
  call calc_A(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Hp,chi_step,H,work1)
  call calc_A(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Sp,chi_step,S,work2)
#ifdef MPI
  call calc_E(N_occ_loc,N_occ,H,S,ls(2,3),work3)
#else
  call calc_E(N_occ_loc,N_occ,H,S,ls(2,3))
#endif

  ! fourth point
  do i=1,ls_max
    call fit_quartic(ls(1,1:3),ls(2,1:3),ls(3,1:2),coeff(0:4))
    call solve_quartic(coeff(0:4),ls(1,4),fail)
    if (fail) call die('ERROR: no minimum in line search!')
    chi_step=chi+ls(1,4)*d
    call calc_A(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Hp,chi_step,H,work1)
    call calc_A(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Sp,chi_step,S,work2)
#ifdef MPI
    call calc_E(N_occ_loc,N_occ,H,S,ls(2,4),work3)
#else
    call calc_E(N_occ_loc,N_occ,H,S,ls(2,4))
#endif
    call calc_grad(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,H,S,g_step,work1,work2)
#ifdef MPI
    ddot_val=ddot(Hp_dim*N_occ_loc,g_step,1,d,1)
    call mpi_allreduce(ddot_val,ls(3,4),1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
#else
    ls(3,4)=ddot(Hp_dim*N_occ_loc,g_step,1,d,1)
#endif
    if (abs(ls(3,4))<ls_tol) then
      chi=chi_step
      g=g_step
      E_omg=ls(2,4)
      step_init=sqrt(step_init*abs(ls(1,4)))
      ls_conv=.true.
      exit
    end if
    if (i==1) then
      ls(1:3,3)=ls(1:3,4)
    else
      if (all(ls(2,1)>=ls(2,2:3))) then
        ls(1:3,1)=ls(1:3,4)
      else if (ls(2,2)>=ls(2,3)) then
        ls(1:3,2)=ls(1:3,4)
      else
        ls(1:3,3)=ls(1:3,4)
      end if
    end if
  end do
  if (i>ls_max) then
    if (Node==0) print*, '#WARNING: Line search failed to converge!'
    chi=chi_step
    g=g_step
    E_omg=ls(2,4)
    ls_conv=.false.
  end if

  call timer_stop('minim_line_search_fit')

end subroutine line_search_fit

!================================================!
! Calculate the coeffs. of the quartic line      !
! search equation from three energy points and   !
! two gradient points                            !
!================================================!
subroutine fit_quartic(x,y,g,c)
  implicit none

  !**** INPUT ***********************************!

  real(dp), intent(in) :: x(1:3) ! Three x-points {x_i}
  real(dp), intent(in) :: y(1:3) ! y(x_i) at the three points
  real(dp), intent(in) :: g(1:2) ! (dy/dx)|x_i at the three points

  !**** OUTPUT **********************************!

  real(dp), intent(out) :: c(0:4) ! Coeffs. of the quartic equation

  !**********************************************!

  !call timer_start('minim_solve_quartic')

  ! The following expressions for the coeffs. were produced automatically using Maple 12
  c(4)=(x(3)**3*x(2)*g(1)-3*x(1)*x(2)**2*y(1)+3*y(3)*x(1)*x(2)**2+x(1)**2*x(2)**2*g(1)+x(3)*x(2)**3*&
       g(1)+2*x(1)**2*x(3)**2*g(2)-3*x(2)*x(3)**2*y(1)+3*y(2)*x(1)**2*x(2)-x(3)**3*x(1)*g(1)+x(3)**3*&
       x(2)*g(2)-x(2)**2*x(3)**2*g(2)-x(1)**2*x(2)**2*g(2)-2*x(2)**2*x(3)**2*g(1)+3*x(2)*x(3)**2*&
       y(2)+x(1)**2*x(3)**2*g(1)-x(1)*x(2)**3*g(1)-3*x(1)*x(3)**2*y(1)-x(3)**3*x(1)*g(2)+3*x(1)*&
       x(3)**2*y(2)+x(2)*g(2)*x(1)**3-3*y(3)*x(1)**2*x(2)-x(3)*g(2)*x(1)**3+6*x(1)*x(3)*x(2)*y(1)+2*&
       x(1)*x(3)*x(2)**2*g(2)+x(1)*x(3)**2*x(2)*g(1)-x(1)*x(3)**2*x(2)*g(2)-2*x(1)**2*x(3)*x(2)*g(1)-&
       x(1)**2*x(3)*x(2)*g(2)+x(1)*x(3)*x(2)**2*g(1)-6*x(1)*x(3)*x(2)*y(2)+2*x(3)**3*y(1)-2*x(3)**3*&
       y(2)+x(2)**3*y(1)+y(3)*x(1)**3-y(3)*x(2)**3-y(2)*x(1)**3)/(-2*x(3)**3*x(1)**4+x(3)**4*x(1)**3-&
       x(1)**2*x(2)**5-x(3)**4*x(2)**3-x(2)**5*x(3)**2-3*x(1)**4*x(2)**3+2*x(3)**3*x(2)**4+x(1)**5*&
       x(3)**2+3*x(2)**4*x(1)**3+4*x(3)**3*x(2)*x(1)**3-4*x(3)**3*x(1)*x(2)**3+2*x(1)*x(3)*x(2)**5+4*&
       x(1)**4*x(3)*x(2)**2+8*x(1)**2*x(3)**2*x(2)**3+x(1)**4*x(3)**2*x(2)-x(1)*x(3)**2*x(2)**4-2*&
       x(1)**5*x(3)*x(2)-4*x(1)**2*x(3)*x(2)**4+x(1)**5*x(2)**2-8*x(2)**2*x(3)**2*x(1)**3-3*x(3)**4*&
       x(1)**2*x(2)+3*x(3)**4*x(1)*x(2)**2)

  c(3)=-(-x(1)*g(1)+2*c(4)*x(1)**4-x(1)*g(2)+4*x(1)*c(4)*x(2)**3+x(2)*g(1)-4*x(2)*c(4)*x(1)**3+x(2)*&
       g(2)-2*c(4)*x(2)**4+2*y(1)-2*y(2))/(x(1)**3+3*x(1)*x(2)**2-3*x(2)*x(1)**2-x(2)**3)

  c(2)=-(-y(2)+c(4)*x(2)**4+c(3)*x(2)**3+x(2)*g(1)-4*x(2)*c(4)*x(1)**3-3*x(2)*c(3)*x(1)**2+y(1)+3*&
       c(4)*x(1)**4+2*c(3)*x(1)**3-x(1)*g(1))/(x(1)**2-2*x(1)*x(2)+x(2)**2)

  c(1)=g(1)-4*c(4)*x(1)**3-3*c(3)*x(1)**2-2*c(2)*x(1)

  c(0)=y(1)-c(4)*x(1)**4-c(3)*x(1)**3-c(2)*x(1)**2-c(1)*x(1)

  !if (Node==0) print*, 'f(x)=',c(4),'*x**4+',c(3),'*x**3+',c(2),'*x**2+',c(1),'*x+',c(0)

  !call timer_stop('minim_fit_quartic')

end subroutine fit_quartic

!================================================!
! Find the minimum for the quartic line search   !
! equation                                       !
!================================================!
subroutine solve_quartic(c,x_min,fail)
  implicit none

  !**** INPUT ***********************************!

  real(dp), intent(in) :: c(0:4) ! Coeffs. of the quartic equation

  !**** OUTPUT **********************************!

  logical, intent(out) :: fail ! Did we fail to find a minimum?

  real(dp), intent(out) :: x_min ! Position of minimum

  !**** LOCAL ***********************************!

  integer :: i, x_order(1:3)

  real(dp) :: t(1:3), z(1:3), a, b, d, Q, R, theta, S, U

  !**********************************************!

  !call timer_start('minim_solve_quartic')

  fail=.false.

  !if (c(4)<0.0_dp) then
  !  if (Node==0) print*, '#WARNING: Function is unbounded!'
  !  !stop
  !end if

  ! In order to find the minimum of the quartic equation, we have to solve a cubic equation. The
  ! following method is taken from Numerical Recipes
  a=3.0_dp*c(3)/(4.0_dp*c(4))
  b=2.0_dp*c(2)/(4.0_dp*c(4))
  if ((abs(b)>=1.0d11) .or. (abs(c(4))<=1.0d-11)) then
    !if (Node==0) print*, '#WARNING: Function is quadratic!'
    x_min=-0.5_dp*c(1)/c(2)
    return
  end if
  d=c(1)/(4.0_dp*c(4))

  Q=(a**2-3.0_dp*b)/9.0_dp
  R=(2.0_dp*a**3-9.0_dp*a*b+27.0_dp*d)/54.0_dp
  if (R**2<Q**3) then
    theta=acos(R/sqrt(Q**3))
    t(1)=-2.0_dp*sqrt(Q)*cos(theta/3.0_dp)-a/3.0_dp
    t(2)=-2.0_dp*sqrt(Q)*cos((theta+2.0_dp*Pi)/3.0_dp)-a/3.0_dp
    t(3)=-2.0_dp*sqrt(Q)*cos((theta-2.0_dp*Pi)/3.0_dp)-a/3.0_dp
    z(1:3)=c(4)*t(1:3)**4+c(3)*t(1:3)**3+c(2)*t(1:3)**2+c(1)*t(1:3)+c(0)
    if (c(4)>0.0_dp) then
      if (all(z(1)>=z(2:3))) then
        x_order(1:3)=(/1,2,3/)
      else if (z(2)>z(3)) then
        x_order(1:3)=(/2,3,1/)
      else
        x_order(1:3)=(/3,1,2/)
      end if
      if ((0.0_dp<=t(x_order(1))) .and. (t(x_order(2))<=t(x_order(1)))) then
        x_min=t(x_order(2))
      else
        x_min=t(x_order(3))
      end if
    else
      if (all(z(1)<=z(2:3))) then
        x_min=t(1)
      else if (z(2)<z(3)) then
        x_min=t(2)
      else
        x_min=t(3)
      end if
    end if
  else
    S=-sign(1.0_dp,R)*(abs(R)+sqrt(R**2-Q**3))**(1.0_dp/3.0_dp)
    if (S==0.0_dp) then
      U=0.0_dp
    else
      U=Q/S
    end if
    x_min=(S+U)-(a/3.0_dp)
    if (c(4)<0.0_dp) fail=.true.
  end if

  !call timer_stop('minim_solve_quartic')

end subroutine solve_quartic

!================================================!
! Calculate the gradient of the Kim functional:  !
! G=2*(2*h*C-s*C*H-h*C*S)                        !
!================================================!
subroutine calc_grad(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,H,S,grad,work1,work2)
  implicit none

  !**** INPUT ***********************************!

  integer, intent(in) :: Hp_dim_loc ! Num. of orbitals (local)
  integer, intent(in) :: Hp_dim     ! Num. of orbitals (global)
  integer, intent(in) :: N_occ_loc  ! Num. of (L)WFs (local)
  integer, intent(in) :: N_occ      ! Num. of (L)WFs (global)

  real(dp), intent(in) :: H(:,:) ! Hamiltonian matrix in (L)WF basis
  real(dp), intent(in) :: S(:,:) ! Overlap matrix in (L)WF basis

  !**** INOUT ***********************************!

  real(dp), intent(inout) :: grad(:,:)  ! Gradient of Kim functional
  real(dp), intent(inout) :: work1(:,:)
  real(dp), intent(inout) :: work2(:,:)

  !**********************************************!

  call timer_start('minim_calc_grad')

  grad=4.0_dp*work1

#ifdef MPI
  call pdsymm('R','U',Hp_dim,N_occ,-2.0_dp,S,1,1,desc2,work1,1,1,desc3,1.0_dp,grad,1,1,desc3)
  call pdsymm('R','U',Hp_dim,N_occ,-2.0_dp,H,1,1,desc2,work2,1,1,desc3,1.0_dp,grad,1,1,desc3)
#else
  call dsymm('R','U',Hp_dim,N_occ,-2.0_dp,S,N_occ,work1,Hp_dim,1.0_dp,grad,Hp_dim)
  call dsymm('R','U',Hp_dim,N_occ,-2.0_dp,H,N_occ,work2,Hp_dim,1.0_dp,grad,Hp_dim)
#endif

  call timer_stop('minim_calc_grad')

end subroutine calc_grad

!================================================!
! Calculate the energy of the Kim functional:    !
! E=2*Tr(H)-Tr(H*S)                              !
!================================================!
subroutine calc_E(N_occ_loc,N_occ,H,S,E,SH)
  implicit none

  !**** INPUT ***********************************!

  integer, intent(in) :: N_occ_loc ! Num. of (L)WFs (local)
  integer, intent(in) :: N_occ     ! Num. of (L)WFs (global)

  real(dp), intent(in) :: H(:,:) ! Hamiltonian matrix in (L)WF basis
  real(dp), intent(in) :: S(:,:) ! Overlap matrix in (L)WF basis

  !**** INOUT ***********************************!

  real(dp), intent(inout) :: E ! Energy of Kim functional
  real(dp), optional, intent(inout) :: SH(:,:)

  !**** LOCAL ***********************************!

  integer :: i, j

#ifdef MPI
  real(dp) :: TrH, TrSH
  real(dp), external :: pdlatra
#endif

  !**********************************************!

  call timer_start('minim_calc_E')

#ifdef MPI
  TrH=pdlatra(N_occ,H,1,1,desc2)
  call pdsymm('L','U',N_occ,N_occ,1.0_dp,S,1,1,desc2,H,1,1,desc2,0.0_dp,SH,1,1,desc2)
  TrSH=pdlatra(N_occ,SH,1,1,desc2)
  E=2.0_dp*TrH-TrSH
#else
  E=0.0_dp
  do i=1,N_occ
    do j=1,i-1
      E=E-S(i,j)*H(i,j)
    end do
    E=E+2.0_dp*H(i,i)-S(i,i)*H(i,i)
    do j=i+1,N_occ
      E=E-S(i,j)*H(i,j)
    end do
  end do
#endif

  call timer_stop('minim_calc_E')

end subroutine calc_E

!================================================!
! Calculate operator matrix in (L)WF basis       !
!================================================!
subroutine calc_A(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,Ap,C,A,ApC)
  implicit none

  !**** INPUT ***********************************!

  integer, intent(in) :: Hp_dim_loc   ! Num. of orbitals (local)
  integer, intent(in) :: Hp_dim       ! Num. of orbitals (global)
  integer, intent(in) :: N_occ_loc    ! Num. of (L)WFs (local)
  integer, intent(in) :: N_occ        ! Num. of (L)WFs (global)

  real(dp), intent(in) :: Ap(:,:)     ! Operator matrix in orbital basis
  real(dp), intent(in) :: C(:,:)      ! (L)WF coefficients in orbitals basis

  !**** INOUT ***********************************!

  real(dp), intent(inout) :: A(:,:)   ! Operator matrix in (L)WF basis
  real(dp), intent(inout) :: ApC(:,:) ! Work matrix

  !**********************************************!

  call timer_start('minim_calc_A')

#ifdef MPI
  call pdsymm('L','U',Hp_dim,N_occ,       1.0_dp,Ap,1,1,desc1,C,  1,1,desc3,0.0_dp,ApC,1,1,desc3)
  call pdgemm('T','N',N_occ, N_occ,Hp_dim,1.0_dp,C, 1,1,desc3,ApC,1,1,desc3,0.0_dp,A,  1,1,desc2)
#else
  call dsymm('L','U',Hp_dim,N_occ,       1.0_dp,Ap,Hp_dim,C,  Hp_dim,0.0_dp,ApC,Hp_dim)
  call dgemm('T','N',N_occ, N_occ,Hp_dim,1.0_dp,C, Hp_dim,ApC,Hp_dim,0.0_dp,A,  N_occ)
#endif

  call timer_stop('minim_calc_A')

end subroutine calc_A

!================================================!
! Calculate operator matrix in orbital basis     !
!================================================!
subroutine calc_densmat(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,A,C,Ap,CA)
  implicit none

  !**** INPUT ***********************************!

  integer, intent(in) :: Hp_dim_loc ! Num. of orbitals (local)
  integer, intent(in) :: Hp_dim     ! Num. of orbitals (global)
  integer, intent(in) :: N_occ_loc  ! Num. of (L)WFs (local)
  integer, intent(in) :: N_occ      ! Num. of (L)WFs (global)

  real(dp), intent(in) :: A(:,:) ! Operator matrix in (L)WF basis
  real(dp), intent(in) :: C(:,:) ! (L)WF coefficients in orbitals basis

  !**** INOUT ***********************************!

  real(dp), intent(inout) :: Ap(:,:) ! Operator matrix in orbital basis
  real(dp), intent(inout) :: CA(:,:)

  !**********************************************!

  call timer_start('minim_calc_densmat')

#ifdef MPI
  call pdsymm('R','U',Hp_dim,N_occ,       1.0_dp,A, 1,1,desc2,C,1,1,desc3,0.0_dp,CA,1,1,desc3)
  call pdgemm('N','T',Hp_dim,Hp_dim,N_occ,2.0_dp,CA,1,1,desc3,C,1,1,desc3,0.0_dp,Ap,1,1,desc1)
#else
  call dsymm('R','U',Hp_dim,N_occ,       1.0_dp,A, N_occ, C,Hp_dim,0.0_dp,CA,Hp_dim)
  call dgemm('N','T',Hp_dim,Hp_dim,N_occ,2.0_dp,CA,Hp_dim,C,Hp_dim,0.0_dp,Ap,Hp_dim)
#endif

  call timer_stop('minim_calc_densmat')

end subroutine calc_densmat

!================================================!
! Calculate coeffs. of the quartic line search   !
! equation using analytical expressions          !
!================================================!
subroutine calc_coeff(Hp_dim_loc,Hp_dim,N_occ_loc,N_occ,H,S,Hd,Sd,Hdd,Sdd,coeff,SdT)
  implicit none

  !**** INPUT ***********************************!

  integer, intent(in) :: Hp_dim_loc ! Num. of orbitals (local)
  integer, intent(in) :: Hp_dim     ! Num. of orbitals (global)
  integer, intent(in) :: N_occ_loc  ! Num. of (L)WFs (local)
  integer, intent(in) :: N_occ      ! Num. of (L)WFs (global)

  real(dp), intent(in) :: H(:,:) ! Hamiltonian matrix in (L)WF basis
  real(dp), intent(in) :: S(:,:) ! Overlap matrix in (L)WF basis
  real(dp), intent(in) :: Hd(:,:)
  real(dp), intent(in) :: Sd(:,:)
  real(dp), intent(in) :: Hdd(:,:)
  real(dp), intent(in) :: Sdd(:,:)

  !**** INOUT ***********************************!

  real(dp), intent(inout) :: coeff(0:4) ! Coeffs. of the quartic equation
  real(dp), intent(inout) :: SdT(:,:)

  !**** LOCAL ***********************************!

#ifdef MPI
  integer :: info
#else
  integer :: i, j
#endif

  real(dp) :: Tr_loc
  real(dp) :: TrH, TrHS
  real(dp) :: TrHd, TrHdS, TrHSd
  real(dp) :: TrHdd, TrHddS, TrHSdd, TrHdSd, TrHdSdT
  real(dp) :: TrHddSd, TrHdSdd
  real(dp) :: TrHddSdd
  real(dp), external :: ddot, pdlatra

  !**********************************************!

  call timer_start('minim_calc_coeff')

#ifdef MPI
  TrH=pdlatra(N_occ,H,1,1,desc2)
  Tr_loc=ddot(N_occ*N_occ_loc,H,1,S,1)
  call mpi_allreduce(Tr_loc,TrHS,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)

  TrHd=pdlatra(N_occ,Hd,1,1,desc2)
  Tr_loc=ddot(N_occ*N_occ_loc,Hd,1,S,1)
  call mpi_allreduce(Tr_loc,TrHdS,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
  Tr_loc=ddot(N_occ*N_occ_loc,H,1,Sd,1)
  call mpi_allreduce(Tr_loc,TrHSd,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)

  TrHdd=pdlatra(N_occ,Hdd,1,1,desc2)
  Tr_loc=ddot(N_occ*N_occ_loc,Hdd,1,S,1)
  call mpi_allreduce(Tr_loc,TrHddS,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
  Tr_loc=ddot(N_occ*N_occ_loc,H,1,Sdd,1)
  call mpi_allreduce(Tr_loc,TrHSdd,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
  call pdtran(N_occ,N_occ,1.0_dp,Sd,1,1,desc2,0.0_dp,SdT,1,1,desc2)
  Tr_loc=ddot(N_occ*N_occ_loc,Hd,1,SdT,1)
  call mpi_allreduce(Tr_loc,TrHdSd,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
  Tr_loc=ddot(N_occ*N_occ_loc,Hd,1,Sd,1)
  call mpi_allreduce(Tr_loc,TrHdSdT,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)

  Tr_loc=ddot(N_occ*N_occ_loc,Hdd,1,Sd,1)
  call mpi_allreduce(Tr_loc,TrHddSd,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
  Tr_loc=ddot(N_occ*N_occ_loc,Hd,1,Sdd,1)
  call mpi_allreduce(Tr_loc,TrHdSdd,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)

  Tr_loc=ddot(N_occ*N_occ_loc,Hdd,1,Sdd,1)
  call mpi_allreduce(Tr_loc,TrHddSdd,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
#else
  TrH=0.0_dp
  do i=1,N_occ
    TrH=TrH+H(i,i)
  end do
  TrHS=ddot(N_occ*N_occ,H,1,S,1)

  TrHd=0.0_dp
  do i=1,N_occ
    TrHd=TrHd+Hd(i,i)
  end do
  TrHdS=ddot(N_occ*N_occ,Hd,1,S,1)
  TrHSd=ddot(N_occ*N_occ,H,1,Sd,1)

  TrHdd=0.0_dp
  do i=1,N_occ
    TrHdd=TrHdd+Hdd(i,i)
  end do
  TrHddS=ddot(N_occ*N_occ,Hdd,1,S,1)
  TrHSdd=ddot(N_occ*N_occ,H,1,Sdd,1)
  do i=1,N_occ
    do j=1,N_occ
      SdT(j,i)=Sd(i,j)
    end do
  end do
  TrHdSd=ddot(N_occ*N_occ,Hd,1,SdT,1)
  TrHdSdT=ddot(N_occ*N_occ,Hd,1,Sd,1)

  TrHddSd=ddot(N_occ*N_occ,Hdd,1,Sd,1)
  TrHdSdd=ddot(N_occ*N_occ,Hd,1,Sdd,1)

  TrHddSdd=ddot(N_occ*N_occ,Hdd,1,Sdd,1)
#endif

  coeff(0)=2.0_dp*TrH-TrHS
  coeff(1)=2.0_dp*(2.0_dp*TrHd-TrHdS-TrHSd)
  coeff(2)=2.0_dp*(TrHdd-TrHdSd-TrHdSdT)-TrHddS-TrHSdd
  coeff(3)=-2.0_dp*(TrHddSd+TrHdSdd)
  coeff(4)=-TrHddSdd

  call timer_stop('minim_calc_coeff')

end subroutine calc_coeff

!================================================!
! Random number generator                        !
! -initialize with:                              !
!  call rand_init()                              !
! -generate new number with:                     !
!  call random_numer(rn)                         !
!  where where rn is a real(dp) variable         !
!================================================!
subroutine rand_init
  implicit none

  !**** LOCAL ***********************************!

  character(10) :: system_time

  integer :: i, rand_size
  integer, allocatable :: rand_seed(:)

  real(dp) :: rtime, rn

  !**********************************************!

  call random_seed(size=rand_size)
  allocate(rand_seed(1:rand_size))
  call date_and_time(time=system_time)
  read (system_time,*) rtime
  rand_seed=(Node+1)*int(rtime*1000.0_dp)
  call random_seed(put=rand_seed)
  deallocate(rand_seed)

  do i=1,10000
    call random_number(rn)
  end do

end subroutine rand_init

end module m_minim
