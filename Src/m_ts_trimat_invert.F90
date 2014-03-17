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

! This particular solution method relies on solving the GF
! with the tri-diagonalization routine.
! This will leverage memory usage and also the execution time.
! It is customized to deal with calculating the inverted matrix
! only in a certain region.

module m_ts_trimat_invert

  ! We use the general inversion module to obtain the
  ! same routines etc.

  use precision, only : dp
  use class_zTrimat
  use m_trimat_invert, only : attach2piv
  use m_trimat_invert, only : calc_Mnn_inv
  use m_trimat_invert, only : calc_Xn_div_Cn_p1, calc_Yn_div_Bn_m1
  use m_trimat_invert, only : Xn_div_Cn_p1, Yn_div_Bn_m1
  

  implicit none

  private

  ! Current size of the pivoting arrays
  integer, save          :: Npiv = 0
  ! The pivoting array
  integer, save, pointer :: ipiv(:) => null()
  ! determines whether it is attached, or created
  logical, save          :: attached = .false.

  ! Used for BLAS calls (local variables)
  complex(dp), parameter :: z0  = dcmplx( 0._dp, 0._dp)
  complex(dp), parameter :: z1  = dcmplx( 1._dp, 0._dp)
  complex(dp), parameter :: zm1 = dcmplx(-1._dp, 0._dp)

  public :: invert_BiasTriMat_prep
  public :: invert_BiasTriMat_col
  public :: init_BiasTriMat_inversion
  public :: clear_BiasTriMat_inversion
  public :: TriMat_Bias_idxs

contains

  subroutine invert_BiasTriMat_prep(M,Minv,no_BufL, N_Elec,Elecs, has_El)
    use m_mat_invert
    use m_ts_electype

    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in) :: no_BufL
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    logical, intent(in) :: has_El(N_Elec)

    complex(dp), pointer :: Mpinv(:)

    integer :: lsPart, lePart
    integer :: sN, sNm1, sNp1, n
    integer :: iEl, idx_no, off, sCol, eCol
    logical, allocatable :: Mnn_parts(:)
    logical :: piv_initialized

    if ( parts(M) /= parts(Minv) ) then
       call die('Could not calculate the inverse on non equal sized &
            &matrices')
    end if
    if ( parts(M) == 1 ) then
       call die('This matrix is not tri-diagonal')
    end if
    piv_initialized = .true.
    do n = 1 , parts(M) 
       if ( Npiv < nrows_g(M,n) ) piv_initialized = .false.
    end do
    if ( .not. piv_initialized ) then
       call die('Pivoting array for inverting matrix not set.')
    end if

    call timer('V_TM_Pinv',1)

    lsPart = 1
    !if ( present(sPart) ) lsPart = sPart
    lePart = parts(M)
    !if ( present(ePart) ) lePart = ePart

    ! Calculate all Xn/Cn+1
    do n = lePart - 1 , lsPart , -1 
       Mpinv => val(Minv,n+1,n+1)
       sNp1 = nrows_g(M,n+1)
       call calc_Xn_div_Cn_p1(M,Minv, n, Mpinv, sNp1**2 )
    end do
    ! Calculate all Yn/Bn-1
    do n = 2 , lePart
       Mpinv => val(Minv,n-1,n-1)
       sNm1 = nrows_g(M,n-1)
       call calc_Yn_div_Bn_m1(M,Minv, n, Mpinv, sNm1**2 )
    end do

    ! We calculate all the required Mnn
    ! Here it is permissable to overwrite the old A
    off = 0
    do n = lsPart , lePart
       do iEl = 1 , N_Elec
          if ( .not. has_El(iEl) ) cycle
          idx_no = Elecs(iEl)%idx_no - no_BufL
          sNm1 = idx_no
          sN   = nrows_g(M,n)
          sNp1 = idx_no + TotUsedOrbs(Elecs(iEl)) - 1
          if ( which_part(M,sNm1) <= n .and. &
               n <= which_part(M,sNp1) ) then
             ! get number of columns that belongs to
             ! the electrode in the 'n' diagonal part
             ! this means we only calculate "what is needed"
             sCol = max(sNm1 - off ,  1)
             eCol = min(sNp1 - off , sN)
             call calc_Mnn_inv_cols(M,Minv,n,sCol,eCol)
             exit ! the electrode loop
          end if
       end do
       off = off + nrows_g(M,n)
    end do

    call timer('V_TM_Pinv',2)

  end subroutine invert_BiasTriMat_prep

  subroutine invert_BiasTriMat_col(M,Minv,no_BufL,El)

    use m_ts_electype

    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in) :: no_BufL
    type(Elec), intent(in) :: El

    complex(dp), pointer :: fullMinv(:)
    complex(dp), pointer :: Mpinv(:), Mp(:)
    complex(dp), pointer :: Xn(:), Yn(:)
    complex(dp), pointer :: z(:), fz(:)

    integer :: nr, np, no, ip
    integer :: idx_no
    integer :: sPart, ePart
    integer :: sColF, eColF, sIdxF, eIdxF
    integer :: sColT, eColT, sIdxT, eIdxT
    integer :: sN, sNc, sNm1, sNp1, n, s
    integer :: off, off_inv, i
    logical :: piv_initialized

    ! In this routine M should have been processed through invert_PrepTriMat
    ! So M contains all *needed* inv(Mnn) and all Xn/Cn+1 and Yn/Bn-1.
    ! So we will save the result in Minv

    if ( parts(M) /= parts(Minv) ) then
       call die('Could not calculate the inverse on non equal sized &
            &matrices')
    end if
    if ( parts(M) == 1 ) then
       call die('This matrix is not tri-diagonal')
    end if
    piv_initialized = .true.
    do n = 1 , parts(M) 
       if ( Npiv < nrows_g(M,n) ) piv_initialized = .false.
    end do
    if ( .not. piv_initialized ) then
       call die('Pivoting array for inverting matrix not set.')
    end if

    call timer('V_TM_inv',1)

    nr = nrows_g(M)
    np = parts(M)

    no = TotUsedOrbs(El)
    idx_no = El%idx_no - no_BufL

    sPart = which_part(M,idx_no)
    ePart = which_part(M,idx_no+no-1)
    if ( sPart < 1 ) call die('Error in the Bias inversion')
    if ( ePart - sPart + 1 > 2 ) call die('Error in trimat partition')
    if ( ePart > parts(M) ) call die('Error in the Bias inversion')

    ! Point to the matrices
    z => val(Minv)

    ! First we need to copy over the Mnn with the electrode part!

    ! with this offset we can calculate the column offset for
    ! the current part
    off = 0
    do n = 1 , sPart - 1
       off = off + nrows_g(M,n)
    end do
    idx_no = idx_no - off
    if ( idx_no <= 0 ) call die('Error in electrode setup')
    do n = sPart , ePart

       ! current count of orbitals in the tri-diagonal segment
       if ( n > 1  ) sNm1 = nrows_g(M,n-1)
                     sN   = nrows_g(M,n  )
       if ( n < np ) sNp1 = nrows_g(M,n+1)
       
       ! placement of the already inverted matrix
       Mp => val(M,n,n)

       ! get number of columns that belongs to
       ! the electrode in the 'n' diagonal part
       sColF = max(idx_no          ,  1)
       eColF = min(idx_no + no - 1 , sN)
       if ( eColF < sColF ) &
            call die('Here: Something went wrong')
       sIdxF = (sColF-1) * sN + 1
       eIdxF =  eColF    * sN

       ! get placement of the diagonal block in the column
       call TriMat_Bias_idxs(M,no,n,sIdxT,eIdxT)
       Mpinv => z(sIdxT:eIdxT)

       ! get the placement in the inversed column
       if ( 1 <= idx_no ) then
          ! we are taking the first part of the inversed matrix
          sColT = 1
          eColT = min(eColF-sColF+1,no)
       else
          sColT = -idx_no + 2 ! we have to pass zero
          eColT = no
       end if
       sIdxT = (sColT-1) * sN + 1
       eIdxT =  eColT    * sN

       if ( eIdxT - sIdxT /= eIdxF - sIdxF ) & 
            call die('Error in determining column')
       Mpinv(sIdxT:eIdxT) = Mp(sIdxF:eIdxF)

       if ( sPart == ePart ) cycle ! we have everything! :)

       ! We need to calculate the remaining 
       ! inverted matrix (they currently reside 
       ! with the neighbouring cells)

       ! first calculate the missing number of columns
       sNc = no - (eColF - sColF + 1)

       if ( n == sPart ) then
          ! we miss the right part of Mnn
          ! we thus need the
          ! Mmn = -Ym+1/Bm * Mm+1n

          Yn => Yn_div_Bn_m1(M,n+1)

          ! placement of the inverted matrix
          Mp => val(M,n+1,n+1)

          call zgemm('N','N',sN,sNc,sNp1, &
               zm1, Yn, sN, Mp(1),sNp1,z0, Mpinv(eIdxT+1),sN)

       else if ( n == ePart ) then
          ! we miss the left part of Mnn
          ! we thus need the
          ! Mmn = -Xm-1/Cm * Mm-1n

          Xn => Xn_div_Cn_p1(M,n-1)

          ! placement of the inverted matrix
          Mp => val(M,n-1,n-1)
          s = sNm1 ** 2

          call zgemm('N','N',sN,sNc,sNm1, &
               zm1, Xn, sN, Mp(s-sNm1*sNc+1),sNm1,z0, Mpinv(1),sN)

       end if

       ! update offset on rows
       off = off + sN
       idx_no = idx_no - sN

    end do

    ! We now have inv(Mnn) in the correct place.

    ! We now calculate:
    !  Mmn = -Ym+1/Bm * Mm+1n, for m<n
    do n = sPart - 1 , 1 , - 1
       
       sN   = nrows_g(M,n)
       sNp1 = nrows_g(M,n+1)

       ! get Ym+1/Bm
       Yn => Yn_div_Bn_m1(M,n+1)

       ! Get Mm+1n
       call TriMat_Bias_idxs(M,no,n+1,sIdxF,eIdxF)
       Mp    => z(sIdxF:eIdxF)
       ! Get Mmn
       call TriMat_Bias_idxs(M,no,n,sIdxF,eIdxF)
       Mpinv => z(sIdxF:eIdxF)

       call zgemm('N','N',sN,no,sNp1, &
            zm1, Yn, sN, Mp(1),sNp1,z0, Mpinv(1),sN)
       
    end do

    ! We now calculate:
    !  Mmn = -Xm-1/Cm * Mm-1n, for m>n
    do n = ePart + 1 , np 

       sNm1 = nrows_g(M,n-1)
       sN   = nrows_g(M,n)

       ! get Xm-1/Cm
       Xn => Xn_div_Cn_p1(M,n-1)
       
       ! Get Mm-1n
       call TriMat_Bias_idxs(M,no,n-1,sIdxF,eIdxF)
       Mp    => z(sIdxF:eIdxF)
       ! Get Mmn
       call TriMat_Bias_idxs(M,no,n,sIdxF,eIdxF)
       Mpinv => z(sIdxF:eIdxF)
       
       call zgemm('N','N',sN,no,sNm1, &
            zm1, Xn, sN, Mp(1),sNm1,z0, Mpinv(1),sN)
       
    end do

    call timer('V_TM_inv',2)

  end subroutine invert_BiasTriMat_col

  ! We will partition the system by:
  ! 1. nrows_g(tri,1) x no
  ! 2. nrows_g(tri,2) x no
  ! ...
  subroutine TriMat_Bias_idxs(M,no,p,sIdx,eIdx)
    type(zTriMat), intent(in) :: M
    ! no is the number of orbitals we wish to take out
    ! p is the part that we wish to point to
    integer, intent(in) :: no, p
    integer, intent(out) :: sIdx, eIdx
    integer :: cum

    ! the size of the partition
    cum = 0
    ! we are requesting the first column,
    ! hence we order the matrix in from the
    ! beginning...
    do eIdx = 1 , p - 1
       cum = cum + nrows_g(M,eIdx)
    end do
    ! This is the number of elements already occupied
    sIdx = no * cum + 1
    eIdx = sIdx + no * nrows_g(M,p) - 1

  end subroutine TriMat_Bias_idxs

  subroutine calc_Mnn_inv_cols(M,Minv,n,sCol,eCol)
    use m_mat_invert
    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in) :: n, sCol, eCol
    ! Local variables
    complex(dp), pointer :: Mp(:), Mpinv(:)
    complex(dp), pointer :: Xn(:), Yn(:), Cn(:), Bn(:)
    integer :: sNm1, sN, sNp1, ierr, i, j

    if ( 1 < n )        sNm1 = nrows_g(M,n-1)
                        sN   = nrows_g(M,n)
    if ( n < parts(M) ) sNp1 = nrows_g(M,n+1)

    if ( sCol == 1 .and. eCol == sN ) then
       call calc_Mnn_inv(M,Minv,n)
       return
    end if
    
    ! Retrieve Ann
    Mp => val(M,n,n)
    if ( n == 1 ) then
       ! First we calculate M11^-1
       ! Retrieve the X1/C2 array
       Xn => Xn_div_Cn_p1(Minv,n)
       ! The C2 array
       Cn => val(M,n,n+1)
       ! Calculate: A1 - X1
       call zgemm('N','N',sN,sN,sNp1, &
            zm1, Cn,sN, Xn,sNp1,z1, Mp,sN)

    else if ( n == parts(M) ) then

       ! Retrieve the Yn/Bn-1 array
       Yn => Yn_div_Bn_m1(Minv,n)
       ! The Bn-1 array
       Bn => val(M,n,n-1)
       ! Calculate: An - Yn
       call zgemm('N','N',sN,sN,sNm1, &
            zm1, Bn,sN, Yn,sNm1,z1, Mp,sN)

    else
       ! Retrieve the Xn/Cn+1 array
       Xn => Xn_div_Cn_p1(Minv,n)
       ! The Cn+1 array
       Cn => val(M,n,n+1)
       ! Calculate: An - Xn
       call zgemm('N','N',sN,sN,sNp1, &
            zm1, Cn,sN, Xn,sNp1,z1, Mp,sN)
       ! Retrieve the Yn/Bn-1 array
       Yn => Yn_div_Bn_m1(Minv,n)
       ! The Bn-1 array
       Bn => val(M,n,n-1)
       ! Calculate: An - Xn - Yn
       call zgemm('N','N',sN,sN,sNm1, &
            zm1, Bn,sN, Yn,sNm1,z1, Mp,sN)
  
    end if

    ! Retrieve the position in the inverted matrix
    Mpinv => val(Minv,n,n)
    Mpinv((sCol-1)*sN+1:eCol*sN) = dcmplx(0._dp,0._dp)
    do j = sCol - 1 , eCol - 1
       Mpinv(j * sN + j + 1) = dcmplx(1._dp,0._dp)
    end do
    j = eCol - sCol + 1
    call zgesv(sN,j,Mp,sN,ipiv,Mpinv((sCol-1)*sN+1),sN,ierr)
    if ( ierr /= 0 ) then
       call die('Error in inverting the partial bias block')
    end if

  end subroutine calc_Mnn_inv_cols

  ! We initialize the pivoting array for rotating the inversion
  subroutine init_BiasTriMat_inversion(M)
    use alloc, only : re_alloc
    type(zTriMat), intent(in) :: M
    integer :: i

    call clear_BiasTriMat_inversion()
    Npiv = 0
    do i = 1 , parts(M)
       if ( nrows_g(M,i) > Npiv ) then
          Npiv = nrows_g(M,i)
       end if
    end do
    call attach2piv(Npiv,ipiv,i)
    attached = .true.

    if ( i /= 0 ) then
       attached = .false.
       ! Allocate space for the pivoting array
       call re_alloc(ipiv,1, Npiv, &
            name="TriMat_piv",routine='TriMatInversion')
    end if

  end subroutine init_BiasTriMat_inversion

  subroutine clear_BiasTriMat_inversion()
    use alloc, only: de_alloc
    if ( Npiv == 0 ) return

    Npiv = 0
    if ( attached ) then
       nullify(ipiv)
       attached = .false.
    else
       ! Deallocate the pivoting array
       call de_alloc(ipiv, &
            name="TriMat_piv",routine='TriMatInversion')
    end if

  end subroutine clear_BiasTriMat_inversion

end module m_ts_trimat_invert
    
