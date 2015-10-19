!
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
!
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

module m_tbt_tri_scat

  use precision, only : dp
  use units, only : Pi
  use m_region

  use class_zTriMat

  use m_ts_tri_scat, only : GF_Gamma_GF
  use m_ts_tri_common, only : GFGGF_needed_worksize

  use m_ts_electype

  implicit none

  private

  public :: A_DOS   ! Spectral function density of states
  public :: GF_DOS  ! Green's function density of states
  public :: A_Gamma ! Calculate the transmission from spectral function . Gamma
  public :: A_Gamma_Block ! Calculate the transmission from spectral function . Gamma (in block form)
  public :: TT_eigen ! Eigenvalue calculation of the transmission eigenvalues
  public :: GF_Gamma ! Calculate the transmission from Green function . Gamma (same-lead contribution)
  public :: GF_T, GF_T_solve
  
  public :: insert_Self_Energy
  public :: insert_Self_Energy_Dev

  ! From ts_tri_scat
  public :: GF_Gamma_GF
  public :: GFGGF_needed_worksize
#ifdef NCDF_4
  public :: orb_current
#endif

  ! Used for BLAS calls (local variables)
  complex(dp), parameter :: z0  = dcmplx( 0._dp, 0._dp)
  complex(dp), parameter :: z1  = dcmplx( 1._dp, 0._dp)
  complex(dp), parameter :: zm1 = dcmplx(-1._dp, 0._dp)
  complex(dp), parameter :: zi  = dcmplx( 0._dp, 1._dp)

contains

  ! A simple routine to calculate the DOS
  ! from a partially calculated GF
  ! When entering this routine Gf_tri
  ! should contain:
  ! all GF_nn
  ! all Yn/Bn-1 and all Xn/Cn+1
  ! This lets us calculate all entries
  subroutine GF_DOS(r,Gf_tri,S_1D,DOS,nwork,work)
    use intrinsic_missing, only : SFIND
    use class_Sparsity
    use class_zSpData1D

    type(tRgn), intent(in) :: r
    type(zTriMat), intent(inout) :: Gf_tri
    type(zSpData1D), intent(inout) :: S_1D
    real(dp), intent(out) :: DOS(r%n)
    integer, intent(in) :: nwork
    complex(dp), intent(inout), target :: work(nwork)

    type(Sparsity), pointer :: sp
    complex(dp), pointer :: S(:), Gf(:), Mnn(:), XY(:)
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    integer :: off1, off2, n, in
    integer :: jo, ii, i, j, no_o, no_i, ind, np

#ifdef TBTRANS_TIMING
    call timer('GF-DOS',1)
#endif

    S  => val(S_1D)
    sp => spar(S_1D)
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col)
    
    ! Initialize DOS to 0
!$OMP parallel workshare default(shared)
    DOS(:) = 0._dp
!$OMP end parallel workshare

    off2 = 0
    np = parts(Gf_tri)
    do n = 1 , np

       no_o = nrows_g(Gf_tri,n)

       do in = max(1,n-1) , min(n+1,np)

          no_i = nrows_g(Gf_tri,in)

          if ( in < n ) then
             off1 = off2 - no_i
          else if ( n < in ) then
             off1 = off2 + no_o
          else
             off1 = off2
          end if

          if ( in == n ) then
             ! Retrieve the central part of the
             ! matrix
             Gf => val(Gf_tri,n,n)

          else

             XY  => val(Gf_tri,in,n)
             Mnn => val(Gf_tri,n,n)

             Gf  => work(1:no_o*no_i)

             ! We need to calculate the 
             ! Mnm1n/Mnp1n Green's function
#ifdef USE_GEMM3M
             call zgemm3m( &
#else
             call zgemm( &
#endif
                  'N','N',no_i,no_o,no_o, &
                  zm1, XY,no_i, Mnn,no_o,z0, Gf,no_i)

          end if

!$OMP parallel do default(shared), private(j,ii,jo,ind,i)
          do j = 1 , no_o
             ii = (j-1) * no_i
             jo = r%r(off2+j)
             ! get the equivalent one in the
             ! overlap matrix
             ! REMEMBER, S is transposed!
             ! Hence we do not need conjg :)
             do i = 1 , no_i
                ind = SFIND(l_col(l_ptr(jo)+1:l_ptr(jo)+ncol(jo)),r%r(off1+i))
                if ( ind == 0 ) cycle
                ind = l_ptr(jo) + ind
                DOS(off2+j) = DOS(off2+j) - dimag( Gf(ii+i) * S(ind) )
             end do
          end do
!$OMP end parallel do

       end do

       ! Update the offset
       off2 = off2 + no_o

    end do

!$OMP parallel workshare default(shared)
    DOS(:) = DOS(:) / Pi
!$OMP end parallel workshare

#ifdef TBTRANS_TIMING
    call timer('GF-DOS',2)
#endif

  end subroutine GF_DOS

#ifdef NCDF_4
#ifdef NOT_WORKING
  ! A simple routine to calculate the DOS
  ! from a partially calculated GF
  ! When entering this routine Gf_tri
  ! should contain:
  ! all GF_nn
  ! all Yn/Bn-1 and all Xn/Cn+1
  ! This lets us calculate all entries
  subroutine GF_DOS_proj(r,Gf_tri,S_1D,N_mol,mols,DOS,bGfk,nwork,work)
    use class_Sparsity
    use class_zSpData1D

    use m_tbt_proj

    type(tRgn), intent(in) :: r
    type(zTriMat), intent(inout) :: Gf_tri
    type(zSpData1D), intent(inout) :: S_1D
    integer, intent(in) :: N_mol
    type(tProjMol), intent(in) :: mols(N_mol)
    real(dp), intent(out) :: DOS(r%n)
    complex(dp), intent(out) :: bGfk(:)
    integer, intent(in) :: nwork
    complex(dp), intent(inout), target :: work(nwork)

    type(Sparsity), pointer :: sp
    complex(dp), pointer :: S(:), Gf(:), Mnn(:), XY(:)
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    integer :: off1, off2, n, in
    integer :: jo, ii, i, j, no_o, no_i, ind, np, iD

    ! For looping the molecule projections
    integer :: Ns, Nl, Nm_dos
    integer :: im, ip, idx, no, i_o
    integer :: step_o

    S  => val(S_1D)
    sp => spar(S_1D)
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col)
    
    ! Initialize DOS to 0
!$OMP parallel workshare default(shared)
    DOS(:) = 0._dp
!$OMP end parallel workshare

    off2 = 0
    np = parts(Gf_tri)

    ! Find maximum size of molecule orbitals
    no = maxval(mols(:)%orb%n)
    ! Calculate the size of the calculated bra at each index
    no_i = 0
    Nm_dos = 0
    do im = 1 , N_mol
       if ( .not. mols(im)%DOS ) cycle
       Nm_dos = Nm_dos + 1
       no_i = no_i + size(mols(im)%proj) * no
    end do
    ! Find maximum work size needed to retain the Gf
    no_o = 0
    do n = 1 , np - 1
       no_o = max(no_o,nrows_g(Gf_tri,n)*nrows_g(Gf_tri,n+1))
    end do
    ! Get the starting position of the projection matrices
    idx = no_o + 1
    ! Calculate maximum number of state orbitals we can
    ! accomodate simultaneously
    max_p = (nwork - no_o) / no_i
    if ( max_p < 1 ) then
       call die('Work size for projection of Gf not sufficient. &
            Try and use fewer projections, or simply do not calculate &
            the DOS projection.')
    end if

    do n = 1 , np

       no_o = nrows_g(Gf_tri,n)

       ! Calculate the step size for the projection
       ! on this column
       step_o = min(max_p,no_o)

       ! Loop over smaller group of columns in this block-column
       do i_o = 1 , no_o, step_o

       im = 0
       do i = 1 , N_mol
          if ( .not. mols(i)%DOS ) cycle
          no = mols(i)%orb%n
          Ns = size(mols(i)%proj)
          ! step calculated DOS for molecule
          im = im + 1
!$OMP parallel do default(shared), private(j,ip,ii), collapse(2)
          do j = 1 , step_o
             ! Calculate the projection matrix on these column
             ! indices
             do ip = 1 , Ns
                ! We have all molecules
                ii = idx + (((im-1)*step_o+j-1)*Ns+ip-1) * no + 1
                call proj_state_bra(mols(i),mols(i)%proj(ip), &
                     i_o+j, zwork(ii:ii+no-1) )
             end do
          end do
!$OMP end parallel do
       end do

       do in = max(1,n-1) , min(n+1,np)

          no_i = nrows_g(Gf_tri,in)

          if ( in < n ) then
             off1 = off2 - no_i
          else if ( n < in ) then
             off1 = off2 + no_o
          else
             off1 = off2
          end if

          if ( in == n ) then
             ! Retrieve the central part of the
             ! matrix
             Gf => val(Gf_tri,n,n)
             ! re-point
             Gf => Gf((i_o-1)*no_o+1:)

          else

             XY  => val(Gf_tri,in,n)
             Mnn => val(Gf_tri,n,n)
             ! re-point
             Mnn => Mnn((i_o-1)*no_o+1:)

             Gf  => work(1:no_o*no_i)

             ! We need to calculate the 
             ! Mnm1n/Mnp1n Green's function
#ifdef USE_GEMM3M
             call zgemm3m( &
#else
             call zgemm( &
#endif
                  'N','N',no_i,step_o,no_o, &
                  zm1, XY,no_i, Mnn,no_o,z0, Gf,no_i)

          end if

!$OMP parallel do default(shared), private(j,ii,jo,ind,i,ip,im,iD)
          do j = 1 , step_o
             ii = (j-1) * no_i
             iD = off2 + j
             jo = r%r(iD)
             ! get the equivalent one in the
             ! overlap matrix
             ! REMEMBER, S is transposed!
             ! Hence we do not need conjg :)
             do i = 1 , no_i
                ind = SFIND(l_col(l_ptr(jo)+1:l_ptr(jo)+ncol(jo)),r%r(off1+i))
                if ( ind == 0 ) cycle
                ind = l_ptr(jo) + ind
                DOS(iD) = DOS(iD) - dimag( Gf(ii+i) * S(ind) )
             end do
          end do
!$OMP end parallel do

       end do

       ! Update the offset
       off2 = off2 + no_o

    end do

!$OMP parallel workshare default(shared)
    DOS(:) = DOS(:) / Pi
!$OMP end parallel workshare

  contains

    subroutine calc_state_Gf(N_mol,mols,Gf,step_o,zw,bGfk)
      
      iG = 0
       im = 0
       do i = 1 , N_mol
          if ( .not. mols(i)%DOS ) cycle
          no = mols(i)%orb%n
          Ns = size(mols(i)%proj)
          ! step calculated DOS for molecule
          im = im + 1
          do ip = 1 , Ns
             ! We have all molecules
             ii = idx + (((im-1)*step_o+j-1)*Ns+ip-1) * no + 1
             iG = iG + 1
             bGfk(iG) = bGfk(iG) + zw(
          end do
       end do
    end subroutine calc_state_Gf

  end subroutine GF_DOS_PROJ

#endif
#endif

  ! A simple routine to calculate the DOS
  ! from a full calculated spectral function
  subroutine A_DOS(r,A_tri,S_1D,DOS)
    use intrinsic_missing, only : SFIND
    use class_Sparsity
    use class_zSpData1D

    type(tRgn), intent(in) :: r
    type(zTriMat), intent(inout) :: A_tri
    type(zSpData1D), intent(inout) :: S_1D
    real(dp), intent(out) :: DOS(r%n)

    type(Sparsity), pointer :: sp
    complex(dp), pointer :: S(:), A(:)
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    integer :: off1, off2, n, in
    integer :: jo, ii, i, j, no_o, no_i, ind, np

#ifdef TBTRANS_TIMING
    call timer('A-DOS',1)
#endif

    S  => val(S_1D)
    sp => spar(S_1D)
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col)

    ! Initialize DOS to 0
!$OMP parallel workshare default(shared)
    DOS(:) = 0._dp
!$OMP end parallel workshare

    off2 = 0
    np = parts(A_tri)
    do n = 1 , np

       no_o = nrows_g(A_tri,n)

       do in = max(1,n-1) , min(n+1,np)

          A => val(A_tri,in,n)

          no_i = nrows_g(A_tri,in)

          if ( in < n ) then
             off1 = off2 - no_i
          else if ( n < in ) then
             off1 = off2 + no_o
          else
             off1 = off2
          end if

!$OMP parallel do default(shared), private(j,ii,jo,ind,i)
          do j = 1 , no_o
             ii = (j-1) * no_i
             jo = r%r(off2+j)
             ! get the equivalent one in the
             ! overlap matrix
             ! REMEMBER, S is transposed!
             ! Hence we are doing it correctly
             do i = 1 , no_i
                ind = SFIND(l_col(l_ptr(jo)+1:l_ptr(jo)+ncol(jo)),r%r(off1+i))
                if ( ind == 0 ) cycle
                ind = l_ptr(jo) + ind
                DOS(off2+j) = DOS(off2+j) + dreal( A(ii+i) * S(ind) )
             end do
          end do
!$OMP end parallel do

       end do

       ! Update the offset
       off2 = off2 + no_o

    end do

    ! The spectral function has a factor two

!$OMP parallel workshare default(shared)
    DOS(:) = DOS(:) / (2._dp * Pi)
!$OMP end parallel workshare


#ifdef TBTRANS_TIMING
    call timer('A-DOS',2)
#endif

  end subroutine A_DOS

  
  ! The simplest routine to do the transport calculation
  ! It takes the spectral function and multiplies it with
  ! the scattering matrix of the down-projected self-energy
  ! and calculates the transmission.
  ! We do this by taking advantage of the transposed scattering
  ! matrix: \Gamma
  subroutine A_Gamma(A_tri,El,T)

    type(zTriMat), intent(inout) :: A_tri ! Spectral function
    type(Elec), intent(in) :: El
    real(dp), intent(out) :: T

    ! Here we need a double loop
    integer :: no
    integer :: i_Elec, ii, isN, in, A_i
    integer :: j_Elec, jj, jsN, jn, A_j
    integer :: o
    integer, pointer :: crows(:)
    complex(dp), pointer :: A(:)

    ! External BLAS routine
    complex(dp), external :: zdotu
    
#ifdef TBTRANS_TIMING
    call timer('A-Gamma',1)
#endif

    ! Get data from tri-diagonal matrix
    crows => cum_rows(A_tri)

    no = El%inDpvt%n

    ! This code is based on the down-folded self-energies
    ! which are determined by the col region
    T = 0._dp

    ! Loop columns
    i_Elec = 1
    do while ( i_Elec <= no ) 

       ! We start by creating a region of consecutive memory.
       call consecutive_index(A_tri,El,i_Elec,in,ii)
       isN = nrows_g(A_tri,in)

       ! Get starting placement of column in the current block
       ! of the spectral function (zero based)
       if ( in == 1 ) then
          A_i = El%inDpvt%r(i_Elec) - 1
       else
          A_i = El%inDpvt%r(i_Elec) - crows(in-1) - 1
       end if

       if ( ii == no ) then
          
          ! The easy calculation, note that ii == no, only
          ! if the entire electrode sits in one block
          A => val(A_tri,in,in)
          do o = 0 , no - 1
             T = T - aimag(zdotu(no,A((A_i+o)*isN+A_i+1),1,El%Gamma(o*no+1),1))
          end do
          
          ! Quick break of loop
          exit

       end if

       ! Loop rows
       j_Elec = 1
       do while ( j_Elec <= no ) 

          ! We start by creating a region of consecutive memory.
          call consecutive_index(A_tri,El,j_Elec,jn,jj)
          jsN = nrows_g(A_tri,jn)

          ! Get the block with the spectral function
          A => val(A_tri,jn,in)

          if ( jn == 1 ) then
             A_j = El%inDpvt%r(j_Elec)
          else
             A_j = El%inDpvt%r(j_Elec) - crows(jn-1)
          end if

          do o = 0 , ii - 1
             T = T - aimag(zdotu(jj,A((A_i+o)*jsN+A_j),1, &
                  El%Gamma((i_Elec-1+o)*no+j_Elec),1))
          end do

          j_Elec = j_Elec + jj

       end do
       
       i_Elec = i_Elec + ii

    end do

#ifdef TBTRANS_TIMING
    call timer('A-Gamma',2)
#endif
    
  end subroutine A_Gamma

  
  ! On entry A_tri is the spectral function
  ! on return the first El%o_inD%n x El%o_inD%n will be the
  ! G.Gamma.Gf.El%Gamma matrix
  ! This will enable eigenvalue calculators and possibly
  ! speed up the calculation of the transmission.
  subroutine A_Gamma_Block(A_tri,El,T,nwork,work)

    use intrinsic_missing, only : transpose, trace
    
    type(zTriMat), intent(inout) :: A_tri ! Spectral function
    type(Elec), intent(in) :: El
    real(dp), intent(out) :: T
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)

    ! Here we need a double loop
    integer :: no
    integer :: i_Elec, ii, isN, in, A_i
    integer :: j_Elec, jj, jsN, jn, A_j
    integer, pointer :: crows(:)
    complex(dp), pointer :: A(:)
    complex(dp) :: z
    
#ifdef TBTRANS_TIMING
    call timer('A-Block-Gamma',1)
#endif

    ! Get data from tri-diagonal matrix
    crows => cum_rows(A_tri)

    no = El%inDpvt%n
    if ( no ** 2 > nwork ) then
       call die('A_Gamma_Block: Insufficient work-size')
    end if

    ! "sadly" Gamma is saved in transposed form, hence
    ! we transpose, and return it to original form, when returning
    call transpose(no,El%Gamma)

    ! This code is based on the down-folded self-energies
    ! which are determined by the col region

    ! Loop columns
    i_Elec = 1
    ! The first column calculation initializes the result
    z = z0
    do while ( i_Elec <= no ) 

       ! We start by creating a region of consecutive memory.
       call consecutive_index(A_tri,El,i_Elec,in,ii)
       isN = nrows_g(A_tri,in)

       ! Get starting placement of column in the current block
       ! of the spectral function (zero based)
       if ( in == 1 ) then
          A_i = El%inDpvt%r(i_Elec) - 1
       else
          A_i = El%inDpvt%r(i_Elec) - crows(in-1) - 1
       end if

       if ( ii == no ) then
          ! The easy calculation, note that ii == no, only
          ! if the entire electrode sits in one block
          A => val(A_tri,in,in)
          
#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
              'N','N',no,no,no, zi, A(A_i*(isN+1)+1), isN, &
              El%Gamma(1), no, z0, work(1), no)

          ! Quick break of loop
          exit

       end if

       ! Loop rows
       j_Elec = 1
       do while ( j_Elec <= no ) 

          ! We start by creating a region of consecutive memory.
          call consecutive_index(A_tri,El,j_Elec,jn,jj)
          jsN = nrows_g(A_tri,jn)

          ! Get the block with the spectral function
          A => val(A_tri,jn,in)

          if ( jn == 1 ) then
             A_j = El%inDpvt%r(j_Elec)
          else
             A_j = El%inDpvt%r(j_Elec) - crows(jn-1)
          end if

#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
              'N','N',jj,no,ii, zi, A(A_i*jsN + A_j), jsN, &
              El%Gamma(i_Elec), no, z, work(j_Elec), no)

          j_Elec = j_Elec + jj

       end do
       
       i_Elec = i_Elec + ii
       ! Now we have already filled the first entries, sum...
       z = z1

    end do

    ! Calculate transmission
    T = dreal(trace(no,work))
    
    ! Now we have the square matrix product
    !   tt = G \Gamma_1 G^\dagger \Gamma_El

    call transpose(no,El%Gamma)

#ifdef TBTRANS_TIMING
    call timer('A-Block-Gamma',2)
#endif
    
  end subroutine A_Gamma_Block


  subroutine TT_eigen(n,tt,nwork,work,eig)
    integer, intent(in) :: n
    complex(dp), intent(inout) :: tt(n*n)
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)
    complex(dp), intent(inout) :: eig(n)

    real(dp) :: rwork(n*2)
    complex(dp) :: z
    integer :: i, j

#ifdef TBTRANS_TIMING
    call timer('TT-eig',1)
#endif

    ! To remove any singular values we add a 1e-3 to the diagonal
    do i = 1 , n
       tt((i-1)*n+i) = tt((i-1)*n+i) + 1.e-3_dp
    end do
    call zgeev('N','N',n,tt,n,eig,work(1),1,work(1),1, &
         work,nwork,rwork,i)
    if ( i /= 0 ) then
       print *,i
       call die('TT_eigen: Could not calculate eigenvalues.')
    end if

    ! Sort the eigenvalues, and simultaneously shift them back
    eig(1) = eig(1) - 1.e-3_dp
    do i = 2 , n
       eig(i) = eig(i) - 1.e-3_dp
       do j = 1 , i - 1
          if ( dreal(eig(j)) < dreal(eig(i)) ) then
             z = eig(j)
             eig(j) = eig(i)
             eig(i) = z
          end if
       end do
    end do

#ifdef TBTRANS_TIMING
    call timer('TT-eig',2)
#endif
    
  end subroutine TT_eigen
  
  subroutine GF_Gamma(Gfcol,El,T)

    use m_ts_trimat_invert, only : TriMat_Bias_idxs

    type(zTriMat), intent(inout) :: Gfcol
    type(Elec), intent(inout) :: El
    real(dp), intent(out) :: T

    complex(dp), pointer :: Gf(:)
    complex(dp), pointer :: z(:)

    integer :: no, np
    integer :: i, ii, i_Elec
    integer, allocatable :: cumsum(:)
    integer :: sN, n, nb
    ! BLAS routines
    complex(dp), external :: zdotu, zdotc

#ifdef TBTRANS_TIMING
    call timer('Gf-Gamma',1)
#endif

    no = El%inDpvt%n
    np = parts(Gfcol)
    allocate(cumsum(np))
    cumsum(1) = 0
    do n = 2 , np
       cumsum(n) = cumsum(n-1) + nrows_g(Gfcol,n-1)
    end do

    ! This code is based on the down-folded self-energies
    ! which are determined by the col region

    ! Point to the matrices
    z => val(Gfcol,all=.true.)

    T = 0._dp

    i_Elec = 1
    do while ( i_Elec <= no ) 

       ! We start by creating a region of consecutive memory.
       call consecutive_index(Gfcol,El,i_Elec,n,nb)
       sN = nrows_g(Gfcol,n)

       ! get placement of the diagonal block in the column
       call TriMat_Bias_idxs(Gfcol,no,n,i,ii)

       i = i + El%inDpvt%r(i_Elec) - cumsum(n) - 1
       Gf => z(i:ii)

#ifdef TBT_T_G_GAMMA_OLD
       
       ! Number of columns that we want to do product of
       ii = 1
       do i = 1 , no
          T = T - zdotu(nb,Gf(ii),1,El%Gamma(i_Elec+(i-1)*no),1) ! G \Gamma
          ii = ii + sN
       end do
       ! Note that Tr[G^\dagger \Gamma] = Tr[ \Gamma G^\dagger ] =
       !    Tr[(G \Gamma)^\dagger]
       ! Hence the below calculation shouldn't be necessary
       ii = (i_Elec - 1) * no + 1
       do i = 1 , nb
          T = T + zdotc(no,Gf(i),sN,El%Gamma(ii),1) ! G^\dagger \Gamma
          ii = ii + no
       end do

#else
       
       ! Note that Tr[G^\dagger \Gamma] = Tr[ \Gamma G^\dagger ] =
       !    Tr[(G \Gamma)^\dagger]
       ! Hence we only calculate one of the contributions and double
       ! it after
       ii = 1
       do i = 1 , no
          T = T - zdotu(nb,Gf(ii),1,El%Gamma(i_Elec+(i-1)*no),1) ! G \Gamma
          ii = ii + sN
       end do
       
#endif

       i_Elec = i_Elec + nb

    end do

    ! Now we have:
    !   T = G \Gamma - G^\dagger \Gamma
#ifndef TBT_T_G_GAMMA_OLD
    T = T * 2._dp
#endif
    deallocate(cumsum)

#ifdef TBTRANS_TIMING
    call timer('Gf-Gamma',2)
#endif
    
  end subroutine GF_Gamma

  subroutine GF_T(Gfcol,El,T_Gf,T_self,nzwork,zwork)

    use intrinsic_missing, only: TRACE
    use m_ts_trimat_invert, only : TriMat_Bias_idxs

    type(zTriMat), intent(inout) :: Gfcol
    type(Elec), intent(inout) :: El
    real(dp), intent(out) :: T_Gf, T_self
    integer, intent(in) :: nzwork
    complex(dp), intent(out) :: zwork(nzwork)

    complex(dp), pointer :: z(:)

    integer :: no, np
    integer :: i, ii, i_Elec
    integer, allocatable :: cumsum(:)
    integer :: sN, n
    integer :: nb
    ! BLAS routines
    complex(dp), external :: zdotu

#ifdef TBTRANS_TIMING
    call timer('Gf-T',1)
#endif

    no = El%inDpvt%n
    np = parts(Gfcol)
    allocate(cumsum(np))
    cumsum(1) = 0
    do n = 2 , np
       cumsum(n) = cumsum(n-1) + nrows_g(Gfcol,n-1)
    end do

#ifndef TS_NOCHECKS
    if ( no**2 > nzwork ) call die('GF_T: no**2 < nzwork')

    ! First we check that we can use the first elements
    ! of Gfcol as temporary storage
    call TriMat_Bias_idxs(Gfcol,no,1,i,ii)
    if ( i < no ** 2 ) then
       write(*,'(a)') 'Remove TBT.T.Gf from your fdf file. &
            &It is not possible in your current setup.'
       call die('GF_T: Size of temporary array not possible.')
    end if

#endif

    ! This code is based on the down-folded self-energies
    ! which are determined by the col region

    ! Point to the matrices
    z => val(Gfcol,all=.true.)

    i_Elec = 1
    do while ( i_Elec <= no ) 

       ! We start by creating a region of consecutive memory.
       call consecutive_index(Gfcol,El,i_Elec,n,nb)
       sN = nrows_g(Gfcol,n)
       
       ! get placement of the diagonal block in the column
       call TriMat_Bias_idxs(Gfcol,no,n,i,ii)
       
       i = i + El%inDpvt%r(i_Elec) - cumsum(n) - 1

       ! Calculate the G \Gamma
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','T',nb,no,no, z1, z(i),sN, &
            El%Gamma(1), no, z0, zwork(i_Elec), no)
       
       i_Elec = i_Elec + nb

    end do

    ! Note that Tr[G^\dagger \Gamma] = Tr[ \Gamma G^\dagger ] =
    !    Tr[(G \Gamma)^\dagger]
    ! Now we have:
    !    = G \Gamma
    T_Gf = - real( TRACE(no,zwork) ,dp) * 2._dp


    ! Now we need to correct for the current electrode
    
    ! Now we can calculate the spectral function for this
    ! electrode where it lives
    i_Elec = 1
    do while ( i_Elec <= no ) 
       
       ! We start by creating a region of consecutive memory.
       call consecutive_index(Gfcol,El,i_Elec,n,nb)
       sN = nrows_g(Gfcol,n)
       
       ! get placement of the diagonal block in the column
       call TriMat_Bias_idxs(Gfcol,no,n,i,ii)
       
       i = i + El%inDpvt%r(i_Elec) - cumsum(n) - 1

       ii = ( i_Elec-1 ) * no + 1
       ! Calculate the G \Gamma G^\dagger
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','C',no,nb,no, z1, zwork(1),no, &
            z(i), sN, z0, z(ii), no)
       
       i_Elec = i_Elec + nb

    end do

    ! The remaining calculation is very easy as we simply
    ! need to do the sum of the trace
    T_self = - real( zdotu(no*no,z(1),1,El%Gamma(1),1) , dp)

    deallocate(cumsum)

#ifdef TBTRANS_TIMING
    call timer('Gf-T',2)
#endif
    
  end subroutine GF_T

  subroutine GF_T_solve(N_Elec,T,has_all)
    integer, intent(in) :: N_Elec
    real(dp), intent(inout) :: T(N_Elec+1,N_Elec)
    ! Whether or not we have calculated all transmissions
    logical, intent(in) :: has_all

    ! System of equations
    real(dp), allocatable :: A(:,:), B(:), work(:)
    real(dp) :: S(N_Elec)
    integer :: lwork

    integer :: i, d, j, Nt, Ne

    ! For one electrodes, we simply return immediately
    if ( N_Elec == 1 ) return

    ! The simple case is when we have 2 electrodes
    if ( N_Elec == 2 ) then

       T(2,1) = T(N_Elec+1,1) - T(1,1)
       if ( has_all ) then
          T(1,2) = T(N_Elec+1,2) - T(2,2)
       else
          T(1,2) = T(2,1)
       end if

       return

    end if

    if ( .not. has_all ) then
       call die('GF_T_solve: Can not separate transmissions.')
    end if

    ! Setup the real system of equations
    Nt = 0
    do i = 1 , N_Elec - 1
       Nt = Nt + i
    end do

    ! Determine work-size (we have too much here)
    lwork = 4 * N_Elec + max(2*N_Elec,Nt)
    
    ! Allocate
    allocate(A(N_Elec,Nt),B(Nt),work(lwork))

    ! Setup RHS and check that we should be able to calculate
    ! the exact transmission using the diagonal part of the
    ! Green function
    Ne = 0
    do i = 1 , N_Elec
       B(i) = T(N_Elec+1,i) - T(i,i)
       if ( i > 1 ) then
        do j = 1 , i - 1
           if ( abs(B(i) - B(j)) < 1.e-4_dp ) then
              Ne = Ne + 1
           end if
        end do
       end if
    end do
    if ( Ne > 0 ) then
       write(*,'(a,i0,a)')'tbt: Transmission from diagonal requires unique &
            &transmission coefficients. You have ',Ne,' similar coefficients.'
       write(*,'(a)')'tbt: This happens with symmetric junctions as &
            &the system is underdetermined.'
       call die('TBtrans cannot calculate the actual transmission using only &
            &the diagonal as two coefficients are the same.')
    end if
    
    
    ! Setup system of equations
    A = 0._dp
    d = 0
    Ne = N_Elec - 1
    do i = 1 , N_Elec

       ! Each partition can be written as
       ! a sequence of 1'nes, and then 
       ! a diagonal matrix below.
       ! Fill rows of A
       A(i,d+1:d+Ne) = 1._dp
       ! Add a diagonal below the current electrode
       do j = 1 , Ne
          A(i+j,d+j) = 1._dp
       end do

       d = d + Ne

       Ne = Ne - 1

    end do

    ! The last two elements of A is also one, always
    ! This corresponds to fill of 1 and diagonal of size 1 below.
    A(N_Elec-1:N_Elec,Nt) = 1._dp
!!$    do i = 1 , N_Elec
!!$       print'(1000(tr1,f2.0))',A(i,:)
!!$    end do

    ! Solve the system of equations.
    ! This is not a fully determined one, but we _know_ it is
    ! Hence, the linear least squares should provide enough 
    ! accuracy.
    ! -1._dp signals machine precision for singular values
    ! However, there shouldn't be any singular values
    ! unless the transmission is zero... ???
    call dgelss(N_Elec,Nt,1,A,N_Elec,B,Nt,S,-1._dp,d,work,lwork,i)
    if ( i /= 0 .or. d /= N_Elec ) then
       call die('GF_T_solve: Could not solve using least-squares.')
    end if

    ! d is the rank of the matrix

    ! Now we have the actual transmissions,
    ! copy over
    d = 0
    do i = 1 , N_Elec - 1
       do j = i + 1 , N_Elec
          d = d + 1
          T(j,i) = B(d)
          T(i,j) = B(d)
       end do
    end do

    deallocate(A,B,work)

  end subroutine GF_T_solve

  subroutine consecutive_index(Tri,El,current,p,n)
    type(zTriMat), intent(inout) :: Tri
    type(Elec), intent(in) :: El
    integer, intent(in) :: current
    integer, intent(out) :: p, n

    ! Local variables
    integer :: idx_Elec, i
    
    idx_Elec = El%inDpvt%r(current)
    p = which_part(Tri,idx_Elec)

    n = 1
    do while ( current + n <= El%inDpvt%n )
       i = El%inDpvt%r(current+n)
       ! In case it is not consecutive
       if ( i - idx_Elec /= n ) exit
       ! In case the block changes, then
       ! we cut the block size here.
       if ( p /= which_part(Tri,i) ) exit
       n = n + 1
    end do

  end subroutine consecutive_index


#ifdef NCDF_4
  subroutine orb_current(spH,A_tri,r,orb_J)

    use class_Sparsity
    use class_zSpData1D
    use class_dSpData1D
    use intrinsic_missing, only : SFIND

    type(zSpData1D), intent(inout) :: spH
    type(zTriMat), intent(inout) :: A_tri
    ! The region that specifies the size of orb_J
    type(tRgn), intent(in) :: r
    type(dSpData1D), intent(inout) :: orb_J

    type(Sparsity), pointer :: i_sp
    integer, pointer :: i_ncol(:), i_ptr(:), i_col(:)
    type(Sparsity), pointer :: sp
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)

    complex(dp), pointer :: H(:)
    complex(dp), pointer :: A(:)
    real(dp), pointer :: J(:)
    integer :: iu, io, ind, iind, idx, ju, jo

#ifdef TBTRANS_TIMING
    call timer('orb-current',1)
#endif

    sp => spar(spH)
    H  => val (spH)
    call attach(sp, n_col=l_ncol, list_ptr=l_ptr, list_col=l_col)

    i_sp => spar(orb_J)
    J    => val (orb_J)
    call attach(i_sp, n_col=i_ncol, list_ptr=i_ptr, list_col=i_col)

    ! We do not initialize J as every entry is overwritten

    A => val(A_tri)

!$OMP parallel do default(shared), private(iu,io,ju,jo,iind,ind,idx)
    do iu = 1, r%n
       io = r%r(iu)

#ifndef TS_NOCHECKS
       if ( i_ncol(io) == 0 ) call die('orb_current: J has zero columns &
            &for at least one row')
#endif

       ! Loop on entries here...
       do ju = 1 , r%n
          ! We search the transposed sparse J matrix as 
          ! nnzs(H) > nnzs(J), always.
          ! J(iind) = J(jo,io)
          jo = r%r(ju)
          iind = SFIND(i_col(i_ptr(jo)+1:i_ptr(jo)+i_ncol(jo)),io)
          if ( iind == 0 ) cycle
          iind = i_ptr(jo) + iind

          ! Check if the orbital exists in the region
          ! We are dealing with a UC sparsity pattern.
          ! this wil ALWAYS be non-zero, note that 
          ! the device region is a subset of the full sparsity
          ! pattern
          ind = l_ptr(io) + &
               SFIND(l_col(l_ptr(io)+1:l_ptr(io)+l_ncol(io)),jo)

          ! H(ind) = H(io,jo) ^ * = H(jo,io)

          ! Notice that H and S are transposed
          jo  = index(A_tri,iu,ju)
          idx = index(A_tri,ju,iu)

          ! Jji = - Im(Hji A_ij - Hji^* A_ji)
          ! I think we need a factor 1/2, but as the units
          ! are more or less never used, I refrain from
          ! dividing by two!
          ! Currently we calculate it using the intrinsic
          ! overlap matrix, however bond-currents are 
          ! not well defined for non-orthogonal basis sets.
          ! We should not shift the Hamiltonian with respect 
          ! to the overlap matrix.
          ! This can easily be seen using the Loewdin basis.
          J(iind) = - aimag( H(ind) * A(jo) - dconjg( H(ind) ) * A(idx) )
          
       end do
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('orb-current',2)
#endif

  end subroutine orb_current
#endif

  subroutine insert_Self_energy(n1,n2,M,r,El,off1,off2)

    ! The sizes of the matrix
    integer, intent(in) :: n1, n2
    complex(dp), intent(inout) :: M(n1,n2)
    ! the region which describes the current segment of insertion
    type(tRgn), intent(in) :: r
    ! Electrodes...
    type(Elec), intent(inout) :: El
    ! The offsets of the matrix
    integer, intent(in) :: off1, off2

    ! local variables
    integer :: j, je, i, ie, no

#ifdef TBTRANS_TIMING
    call timer('insert-SE',1)
#endif

    El%idx_o = El%idx_o - 1
    no = TotUsedOrbs(El)

    ! We are dealing with the intrinsic electrode
    ! self energy
    ! Here we have two options,
    ! Bulk) We are dealing with a bulk electrode
    ! not bulk) A non-bulk electrode

    if ( El%Bulk ) then
!$OMP parallel do default(shared), private(j,je,i,ie)
       do j = 1 , n2
          je = r%r(off2+j) - El%idx_o
          if ( 1 <= je .and. je <= no ) then
          je = (je - 1) * no
          do i = 1 , n1
             ie = r%r(off1+i) - El%idx_o
             if ( ie < 1 ) cycle
             if ( no < ie ) cycle
             
             M(i,j) = El%Sigma(je + ie)
             
          end do
          end if
       end do
!$OMP end parallel do
    else
!$OMP parallel do default(shared), private(j,je,i,ie)
       do j = 1 , n2
          je = r%r(off2+j) - El%idx_o
          if ( 1 <= je .and. je <= no ) then
          je = (je - 1) * no
          do i = 1 , n1
             ie = r%r(off1+i) - El%idx_o
             if ( ie < 1 ) cycle
             if ( no < ie ) cycle
             
             M(i,j) = M(i,j) - El%Sigma(je + ie)
             
          end do
          end if
       end do
!$OMP end parallel do
    end if
    
    El%idx_o = El%idx_o + 1

#ifdef TBTRANS_TIMING
    call timer('insert-SE',2)
#endif

  end subroutine insert_Self_energy


  subroutine insert_Self_energy_Dev(Gfinv_tri,Gfinv,r,El)

    type(zTriMat), intent(inout) :: GFinv_tri
    complex(dp), intent(inout) :: Gfinv(:)
    ! the region which describes the current segment of insertion
    type(tRgn), intent(in) :: r
    type(Elec), intent(in) :: El

    ! local variables
    integer :: j, je, i, ii, idx, no

#ifdef TBTRANS_TIMING
    call timer('insert-SED',1)
#endif

    no = El%o_inD%n

    ! A down-folded self-energy, this
    ! is always considered to be "non-bulk" as
    ! we have it downfolded.

!$OMP parallel do default(shared), private(j,ii,je,i,idx)
    do j = 1 , no
       ii = (j-1)*no
       ! grab the index in the full tri-diagonal matrix
       je = El%inDpvt%r(j)
       do i = 1 , no

          idx = index(GFinv_tri,El%inDpvt%r(i),je)
          
          Gfinv(idx) = Gfinv(idx) - El%Sigma(ii+i)
          
       end do
    end do
!$OMP end parallel do

#ifdef TBTRANS_TIMING
    call timer('insert-SED',2)
#endif

  end subroutine insert_Self_energy_Dev

end module m_tbt_tri_scat
