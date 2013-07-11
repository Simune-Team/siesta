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

  implicit none

  private

  ! Current size of the pivoting arrays
  integer, save          :: Npiv = 0
  ! The pivoting array
  integer, save, pointer :: ipiv(:) => null()

  ! Used for BLAS calls (local variables)
  complex(dp), parameter :: z0  = dcmplx( 0._dp, 0._dp)
  complex(dp), parameter :: z1  = dcmplx( 1._dp, 0._dp)
  complex(dp), parameter :: zm1 = dcmplx(-1._dp, 0._dp)

  public :: invert_BiasTriMat
  public :: init_BiasTriMat_inversion
  public :: clear_BiasTriMat_inversion
  public :: TriMat_Bias_idxs

contains

  subroutine invert_BiasTriMat(UpdateDMCR,M,Minv,no)
    use m_mat_invert
    use intrinsic_missing, only : EYE

    ! If we only need the middle part of the Gf^-1
    logical, intent(in) :: UpdateDMCR
    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in) :: no
    complex(dp), pointer :: fullMinv(:)
    complex(dp), pointer :: Mpinv(:), Mp(:), XYn(:)
    complex(dp), pointer :: Xn(:), Yn(:), Cn(:), Bn(:)
    complex(dp), pointer :: z(:), fz(:)

    integer :: nr, np, ip
    integer :: sNB, sN, sNC
    integer :: sPart, ePart
    integer :: GfsPart, GfePart
    integer :: sIdx, eIdx
    integer :: sNm1, sNp1, n
    integer :: i, ierr
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

    call timer('V_TM_inv',1)

    nr = nrows_g(M)
    np = parts(M)

    ! We can still utilize the trimat-inversion
    if ( no > 0 ) then
       ! Left
       sPart = which_part(M,1)
       ePart = which_part(M,no)
       if ( sPart /= 1 ) call die('Error in the Bias inversion')
       if ( ePart >  2 ) call die('Error in the Bias inversion')
    else
       ! Right
       sPart = which_part(M,nr+no+1)
       ePart = which_part(M,nr)
       if ( sPart <  np-1 ) call die('Error in the Bias inversion')
       if ( ePart /= np   ) call die('Error in the Bias inversion')
    end if

    if ( UpdateDMCR ) then
       ! If we only require part of the full Gf
       ! column
       n = 0
       do ip = 1 , np
          n = n + nrows_g(M,ip)
          GfsPart = ip
          if ( n > abs(no) ) exit
       end do
       n = 0
       do ip = np , 1 , -1
          n = n + nrows_g(M,ip)
          GfePart = ip
          if ( n > abs(no) ) exit
       end do
    else
       ! we need everything
       GfsPart = 1
       GfePart = np
    end if

#ifdef TRANSIESTA_DEBUG 
    print '(a,i0,tr1,i0)','Inverting in parts ',GfsPart, GfePart
#endif


    ! if there is not enough room for the column, error out
    if ( UpdateDMCR ) then
       if ( (nr-nrows_g(M,1) ) * abs(no) > elements(M) .or. &
            (nr-nrows_g(M,np)) * abs(no) > elements(M) ) then
          call die('The matrix column can not be calculated, &
               &the tri-diagonal parts cannot contain the full &
               &column.')
       end if
    else
       if ( nr * abs(no) > elements(M) ) then
          call die('The matrix column can not be calculated, &
               &the tri-diagonal parts cannot contain the full &
               &column.')
       end if
    end if

    ! Notice that we know that the work array we use will never be
    ! overwritten as we have "removed" the diagonal entries 
    ! in Minv

    ! Calculate all Xn/Cn+1
    do n = np - 1 , sPart , -1 
       sNp1 = nrows_g(Minv,n+1) ** 2
       z => work_array(Minv,no,sNp1)
       call calc_Xn_div_Cn_p1(M,Minv, n, no, z,sNp1)
    end do
    ! Calculate all Yn/Bn-1
    do n = 2 , ePart
       sNm1 = nrows_g(Minv,n-1) ** 2
       z => work_array(Minv,no,sNm1)
       call calc_Yn_div_Bn_m1(M,Minv, n, no, z, sNm1)
    end do

    ! Calculate the diagonal entries (this will let us get rid of
    ! the Y/X arrays which we do not need)

    ! Lets note that the first two parts MUST (they are forced to be)
    ! equal or exceed no.
    ! Otherwise the full \Sigma_L/R can not be contained in the tri-matrix
    ! This means that it is safe to calculate the first calculated positions

    if ( no > 0 ) then
       n = sPart
    else
       n = ePart
    end if

    ! We calculate the first diagonal matrix
      
    if ( 1 < n ) sNm1 = nrows_g(M,n-1)
                 sN   = nrows_g(M,n)
    if ( n < np) sNp1 = nrows_g(M,n+1)
    if ( sN > abs(no) ) call die('This system we cannot handle')

    Mpinv => work_array(Minv,no,sN**2)
    Mp    => val(M,n,n)
    
    if ( n == 1 ) then
       ! First we calculate M11^-1
       ! Retrieve the X1/C2 array
       Xn => Xn_div_Cn_p1(Minv,n,no)
       ! The C2 array
       Cn => val(M,n,n+1)
       ! Calculate: A1 - X1
       call zgemm('N','N',sN,sN,sNp1, &
            zm1, Cn,sN, Xn,sNp1,z1, Mp,sN)
       
    else if ( n == np ) then

       ! Retrieve the Yn/Bn-1 array
       Yn => Yn_div_Bn_m1(Minv,n,no)
       ! The Bn-1 array
       Bn => val(M,n,n-1)
       ! Calculate: An - Yn
       call zgemm('N','N',sN,sN,sNm1, &
            zm1, Bn,sN, Yn,sNm1,z1, Mp,sN)

    else

       call die('Something is totally wrong')

    end if

    call mat_invert(Mp,Mpinv,sN)

    fullMinv => val(M)

    ! How many columns do we still need to calculate...?
    sNC = abs(no) - sN

    if ( sNC > 0 ) then
       ! We need to calculate the second column in the
       ! diagonal inverse part

       ! For this we need the Y2 and X2 array ( or Yn-1, Xn-1 )

       ! First calculate the respective Y2 and X2 arrays

       if ( no > 0 ) then
          n = 2
       else
          n = np - 1
       end if
       
       sNm1 = nrows_g(M,n-1)
       sN   = nrows_g(M,n)
       sNp1 = nrows_g(M,n+1)

       ! Retrieve the Ann array
       z => val(M,n,n)
       ! retrieve a work array for retaining the values
       Mp => work_array(Minv,no,sN**2)
       ! copy over Ann
       Mp(:) = z(:)
       
       ! Retrieve the Xn/Cn+1 array
       Xn => Xn_div_Cn_p1(Minv,n,no)
       ! The Cn+1 array
       Cn => val(M,n,n+1)
       ! Calculate: An - Xn
       call zgemm('N','N',sN,sN,sNp1, &
            zm1, Cn,sN, Xn,sNp1,z1, Mp,sN)
       ! Retrieve the Yn/Bn-1 array
       Yn => Yn_div_Bn_m1(Minv,n,no)
       ! The Bn-1 array
       Bn => val(M,n,n-1)
       ! Calculate: An - Xn - Yn
       call zgemm('N','N',sN,sN,sNm1, &
            zm1, Bn,sN, Yn,sNm1,z1, Mp,sN)
       
       ! Put the Mnn in the correct place
       call TriMat_Bias_idxs(M,no,n,sIdx,eIdx)
       if ( no > 0 ) then
          ! Move the start index into the Mnn region
          sIdx = sIdx + nrows_g(M,n) * nrows_g(M,n-1)
       else
          ! Move the end index into the Mnn region
          eIdx = eIdx - nrows_g(M,n) * nrows_g(M,n+1)
       end if

       ! Now we don't need the Bn and Cn's
       ! Hence we can do the remainder of the inversions
       ! without taken notice of them

       ! placement of the inverted matrix
       z => val(M)
       Mpinv => z(sIdx:eIdx)
       Mpinv(:) = dcmplx(0._dp,0._dp)
       if ( no > 0 ) then
          do i = 1 , sNC
             Mpinv((i-1)*sN+i) = dcmplx(1._dp,0._dp)
          end do
       else
          do i = sNC , 1 , -1
             Mpinv(i*sN-(sNC - i)) = dcmplx(1._dp,0._dp)
          end do
       end if

       ! Calculate Mnn^-1
       call zgesv(sN,sNC,Mp,sN,ipiv,Mpinv,sN,ierr)
       if ( ierr /= 0 ) call die('Error on inverting Mnn^-1')

       ! Calculate the remaining parts in the first two sections
       if ( no > 0 ) then

          ! obtain the indices for the inverted matrix above
          call TriMat_Bias_idxs(M,no,n-1,sIdx,eIdx)
          fz => fullMinv(sIdx:eIdx) ! fz is the part x no array

          ! obtain the indices for the current inverted matrix
          call TriMat_Bias_idxs(M,no,n,sIdx,eIdx)
          z  => fullMinv(sIdx:eIdx) ! z is the part x no array

          ! Obtain the X2/C3 matrix
          XYn => Xn_div_Cn_p1(Minv,n-1,no)

          ! First calculate the M21 region
          call zgemm('N','N',sN,sNm1,sNm1, &
               zm1, XYn,sN, fz(1),sNm1,z0, z(1),sN)

          ! Obtain the Y2/B1 matrix
          XYn => Yn_div_Bn_m1(Minv,n,no)
          ! Index move of the 1 row
          sIdx = 1 + sNm1 ** 2
          ! Index move of the 2 row
          eIdx = 1 + sN * sNm1

          ! ... calculate the M12 region
          call zgemm('N','N',sNm1,sNC,sN, &
               zm1, XYn,sNm1, z(eIdx),sN,z0, fz(sIdx),sNm1)

       else

          ! obtain the indices for the inverted matrix below
          call TriMat_Bias_idxs(M,no,n,sIdx,eIdx)
          fz  => fullMinv(sIdx:eIdx) ! fz is the part x no array

          ! obtain the indices for the current inverted matrix
          call TriMat_Bias_idxs(M,no,n+1,sIdx,eIdx)
          z => fullMinv(sIdx:eIdx) ! z is the part x no array

          ! Obtain the Xn/Cn+1 matrix
          XYn => Xn_div_Cn_p1(Minv,n,no)

          ! First calculate the Mn+1n region
          call zgemm('N','N',sNp1,sNC,sN, &
               zm1, XYn,sNp1, fz(1),sN,z0, z(1),sNp1)

          ! Obtain the Y2/B1 matrix
          XYn => Yn_div_Bn_m1(Minv,n+1,no)
          ! Index move of the n row
          sIdx = 1 + sN * sNC
          ! Index move of the n+1 row
          eIdx = 1 + sNp1 * sNC

          ! ... calculate the Mnn+1 region
          call zgemm('N','N',sN,sNp1,sNp1, &
               zm1, XYn,sN, z(eIdx),sNp1,z0, fz(sIdx),sN)

       end if

    end if


    ! Now we are ready to do the calculation of the full column

    sNB = abs(no)

    if ( no > 0 ) then
       
       do ip = ePart + 1 , GfePart
      
          sN   = nrows_g(M,ip-1)
          sNp1 = nrows_g(M,ip)
          XYn => Xn_div_Cn_p1(Minv,ip-1,no)

          ! obtain the indices for the inverted matrix above
          call TriMat_Bias_idxs(M,no,ip-1,sIdx,eIdx)
          fz => fullMinv(sIdx:eIdx) ! fz is the part x no array

          ! obtain the indices for the current inverted matrix
          call TriMat_Bias_idxs(M,no,ip,sIdx,eIdx)
          z  => fullMinv(sIdx:eIdx) ! z is the part x no array

          call zgemm('N','N',sNp1,sNB,sN, &
               zm1, XYn,sNp1, fz(1),sN,z0, z(1),sNp1)

       end do

    else

       do ip = sPart - 1 , GfsPart , -1

          sNm1 = nrows_g(M,ip)
          sN   = nrows_g(M,ip+1)
          XYn => Yn_div_Bn_m1(Minv,ip+1,no)

          ! obtain the indices for the inverted matrix below
          call TriMat_Bias_idxs(M,no,ip+1,sIdx,eIdx)
          fz => fullMinv(sIdx:eIdx) ! fz is the part x no array

          ! obtain the indices for the current inverted matrix
          call TriMat_Bias_idxs(M,no,ip,sIdx,eIdx)
          z  => fullMinv(sIdx:eIdx) ! z is the part x no array

          call zgemm('N','N',sNm1,sNB,sN, &
               zm1, XYn,sNm1, fz(1),sN,z0, z(1),sNm1)

       end do

    end if

    call timer('V_TM_inv',2)

  end subroutine invert_BiasTriMat

  subroutine TriMat_Bias_idxs(M,no,p,sIdx,eIdx)
    type(zTriMat), intent(in) :: M
    ! no is the number of orbitals we wish to take out
    ! p is the part that we wish to point to
    integer, intent(in) :: no, p
    integer, intent(out) :: sIdx, eIdx
    sIdx = TriMat_Bias_idx(M,no,p)
    eIdx = sIdx + abs(no) * nrows_g(M,p) - 1
  end subroutine TriMat_Bias_idxs


  ! We will partition the system by:
  ! 1. nrows_g(tri,1) x no
  ! 2. nrows_g(tri,2) x no
  ! ...
  function TriMat_Bias_idx(M,no,p) result(IDX)
    type(zTriMat), intent(in) :: M
    ! no is the number of orbitals we wish to take out
    ! p is the part that we wish to point to
    integer, intent(in) :: no, p
    integer :: IDX

    integer :: ip, lno, cum

    ! the size of the partition
    lno = abs(no)

    if ( no > 0 ) then

       cum = 0
       ! we are requesting the first column,
       ! hence we order the matrix in from the
       ! beginning...
       do ip = 1 , p - 1
          cum = cum + nrows_g(M,ip)
       end do
       ! This is the number of elements already occupied
       idx = lno * cum

    else
       
       cum = 0
       ! We are requesting the last column
       do ip = p, parts(M)
          cum = cum + nrows_g(M,ip)
       end do
       ! This is the number of elements already free
       idx = elements(M) - lno * cum
       
    end if

    idx = idx + 1

  end function TriMat_Bias_idx


  ! In this case we do not need the diagonal part, so we have
  ! this going from 
  function Xn_div_Cn_p1(M,n,no) result(Xn)
    type(zTriMat), intent(in) :: M
    integer, intent(in) :: n, no
    complex(dp), pointer :: Xn(:), z(:)
    integer :: idx

    ! The last Xn is 0
    if ( n == parts(M) ) return

    z => val(M)
    idx = idx_Xn_div_Cn_p1(M,n,no)
    Xn => z(idx:idx-1+nrows_g(M,n)*nrows_g(M,n+1))
    if ( idx -1+nrows_g(M,n)*nrows_g(M,n+1) > size(z) ) &
         call die('aeostuhsaehu')

  end function Xn_div_Cn_p1

  function Yn_div_Bn_m1(M,n,no) result(Yn)
    type(zTriMat), intent(in) :: M
    integer, intent(in) :: n, no
    complex(dp), pointer :: Yn(:), z(:)
    integer :: idx

    ! The first Yn is 0
    if ( n == 1 ) return

    z => val(M)
    idx = idx_Yn_div_Bn_m1(M,n,no)
    Yn => z(idx:idx-1+nrows_g(M,n)*nrows_g(M,n-1))

    if ( idx -1+nrows_g(M,n)*nrows_g(M,n-1) > size(z) ) &
         call die('aeostuhsaehu')

  end function Yn_div_Bn_m1

  ! In this case we do not need the diagonal part, so we have
  ! this going from 
  function idx_XYn(M,n,no) result(idx)
    type(zTriMat), intent(in) :: M
    integer, intent(in) :: n, no
    integer :: idx, i, lp

    if ( no > 0 ) then

       ! The last part needed for calculating the entire column
       lp = which_part(M,no)
       
       ! In this case we have the Xn taken up the full back
       ! so direct sequence
       idx = elements(M) + 1
       do i = parts(M) - 1 , max(lp,n) + 1 , -1
          idx = idx - nrows_g(M,i) * nrows_g(M,i+1)
       end do
       if ( lp > n ) then ! then lp must also be larger than 1
          ! We need to mingle the Yn and Xn arrays
          do i = lp , n + 1 , -1 
             ! this is Yn
             idx = idx - nrows_g(M,i) * nrows_g(M,i-1)
             ! this is Xn
             idx = idx - nrows_g(M,i) * nrows_g(M,i+1)
          end do
       end if
       ! We are now at the edge of the position, if the Yn exist
       ! add that first
       if ( lp > 1 .and. lp == n ) then
          idx = idx - nrows_g(M,lp) * nrows_g(M,lp-1)
       end if

       ! return the idx of the first element of Xn
       idx = idx - nrows_g(M,n) * nrows_g(M,n+1)

    else
       
       lp = which_part(M,nrows_g(M) - abs(no) + 1)

       ! In this case we have the Yn taking up the full start
       idx = 1
       do i = 2 , min(lp,n) - 1
          idx = idx + nrows_g(M,i) * nrows_g(M,i-1)
       end do
       if ( lp < n ) then ! then lp must also be smaller than parts(M)
          ! We need to mingle the Yn and Xn arrays
          do i = lp , n - 1
             ! this is Xn
             idx = idx + nrows_g(M,i) * nrows_g(M,i+1)
             ! this is Yn
             idx = idx + nrows_g(M,i) * nrows_g(M,i-1)
          end do
       end if

    end if

  end function idx_XYn

  ! In this case we do not need the diagonal part, so we have
  ! this going from 
  function idx_Xn_div_Cn_p1(M,n,no) result(idx)
    type(zTriMat), intent(in) :: M
    integer, intent(in) :: n, no
    integer :: idx

    ! retrieve the index for the equivalent
    idx = idx_XYn(M,n,no)

    if ( no < 0 ) then
       idx = idx + nrows_g(M,n) * nrows_g(M,n-1)
    end if

  end function idx_Xn_div_Cn_p1

  ! In this case we do not need the diagonal part, so we have
  ! this going from 
  function idx_Yn_div_Bn_m1(M,n,no) result(idx)
    type(zTriMat), intent(in) :: M
    integer, intent(in) :: n, no
    integer :: idx

    ! retrieve the index for the equivalent
    idx = idx_XYn(M,n,no)

    if ( no > 0 ) then
       idx = idx + nrows_g(M,n) * nrows_g(M,n+1)
    end if

  end function idx_Yn_div_Bn_m1

  function work_array(M,no,s) result(z)
    type(zTriMat), intent(in) :: M
    integer, intent(in) :: no, s
    complex(dp), pointer :: z(:), fz(:)
    integer :: els
    fz => val(M)
    if ( no > 0 ) then
       z => fz(1:s)
    else
       els = elements(M)
       z => fz(els-s+1:els)
    end if
  end function work_array


! ***************** Direct copies of the routines in m_trimat_invert *************************!
! We have only changed the Yn_div_Bn_m1 and Xn_div_Cn_p1 routines, i.e. we could do with
! passing the functions as an interface, but...

  ! We will calculate the Xn/Cn+1 component of the 
  ! tri-diagonal inversion algorithm.
  ! The Xn/Cn+1 will be saved in the Minv n,n-1 (as that has
  ! the same size).
  subroutine calc_Xn_div_Cn_p1(M,Minv,n,no,zwork,nz)
    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in) :: n, no, nz
    complex(dp), intent(inout) :: zwork(nz)
    ! Local variables
    complex(dp), pointer :: ztmp(:), Xn(:), Cnp2(:)
    integer :: sN, sNp1, sNp1SQ, sNp2, ierr

    if ( n < 1 .or. parts(M) <= n .or. parts(M) /= parts(Minv) ) then
       call die('Could not calculate Xn on these matrices')
    end if
    ! Collect all matrix sizes for this step...
    sN     = nrows_g(M,n)
    sNp1   = nrows_g(M,n+1)
    sNp1SQ = sNp1 ** 2
    if ( nz < sNp1SQ ) then
       call die('Work array in Xn calculation not sufficiently &
            &big.')
    end if

    ! Copy over the Bn array
    ztmp  => val(M   ,n+1,n)
    ! This is where the inverted matrix will be located 
    Xn    => Xn_div_Cn_p1(Minv,n,no)
    Xn(:) =  ztmp(:)

    ! Copy over the An+1 array
    ztmp            => val(M,n+1,n+1)
    zwork(1:sNp1SQ) =  ztmp(:)

    ! If we should calculate X_N-1 then X_N == 0
    if ( n < parts(M) - 1 ) then
       ! Size...
       sNp2 =  nrows_g(M,n+2)
       ! Retrieve the Xn+1/Cn+2 array
       ztmp => Xn_div_Cn_p1(Minv,n+1,no)
       ! Retrieve the Cn+2 array
       Cnp2 => val(M,n+1,n+2)
       ! Calculate: An+1 - Xn+1
       call zgemm('N','N',sNp1,sNp1,sNp2, &
            zm1, Cnp2,sNp1, ztmp,sNp2,z1, zwork,sNp1)
    end if

    ! Calculate Xn/Cn+1
    call zgesv(sNp1,sN,zwork,sNp1,ipiv,Xn,sNp1,ierr)
    if ( ierr /= 0 ) call die('Error on inverting Xn/Cn+1')

  end subroutine calc_Xn_div_Cn_p1

  ! We will calculate the Yn/Bn-1 component of the 
  ! tri-diagonal inversion algorithm.
  ! The Yn/Bn-1 will be saved in the Minv n-1,n (as that has
  ! the same size).
  subroutine calc_Yn_div_Bn_m1(M,Minv,n,no, zwork,nz)
    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in) :: n, no, nz
    complex(dp), intent(inout) :: zwork(nz)
    ! Local variables
    complex(dp), pointer :: ztmp(:), Yn(:), Bnm2(:)
    integer :: sN, sNm1, sNm1SQ, sNm2, ierr

    if ( n < 2 .or. parts(M) < n .or. parts(M) /= parts(Minv) ) then
       call die('Could not calculate Xn on these matrices')
    end if
    ! Collect all matrix sizes for this step...
    sN     = nrows_g(M,n)
    sNm1   = nrows_g(M,n-1)
    sNm1SQ = sNm1 ** 2
    if ( nz < sNm1SQ ) then
       call die('Work array in Xn calculation not sufficiently &
            &big.')
    end if

    ! Copy over the Cn array
    ztmp  => val(M   ,n-1,n)
    ! This is where the inverted matrix will be located 
    Yn    => Yn_div_Bn_m1(Minv,n,no)
    Yn(:) =  ztmp(:)

    ! Copy over the An-1 array
    ztmp            => val(M,n-1,n-1)
    zwork(1:sNm1SQ) =  ztmp(:)

    ! If we should calculate Y_2 then Y_1 == 0
    if ( 2 < n ) then
       ! Size...
       sNm2 =  nrows_g(M,n-2)
       ! Retrieve the Yn-1/Bn-2 array
       ztmp => Yn_div_Bn_m1(Minv,n-1,no)
       ! Retrieve the Bn-2 array
       Bnm2 => val(M,n-1,n-2)
       ! Calculate: An-1 - Yn-1
       call zgemm('N','N',sNm1,sNm1,sNm2, &
            zm1, Bnm2,sNm1, ztmp,sNm2,z1, zwork,sNm1)
    end if

    ! Calculate Yn/Bn-1
    call zgesv(sNm1,sN,zwork,sNm1,ipiv,Yn,sNm1,ierr)
    if ( ierr /= 0 ) call die('Error on inverting Yn/Bn-1')

  end subroutine calc_Yn_div_Bn_m1

  ! We initialize the pivoting array for rotating the inversion
  subroutine init_BiasTriMat_inversion(M)
    use alloc, only : re_alloc
    type(zTriMat), intent(in) :: M
    integer :: i
    Npiv = 0
    do i = 1 , parts(M)
       if ( nrows_g(M,i) > Npiv ) then
          Npiv = nrows_g(M,i)
       end if
    end do

    ! Allocate space for the pivoting array
    call re_alloc(ipiv,1, Npiv, &
         name="TriMat_piv",routine='TriMatInversion')

  end subroutine init_BiasTriMat_inversion

  subroutine clear_BiasTriMat_inversion()
    use alloc, only: de_alloc
    Npiv = 0
    ! Deallocate the pivoting array
    call de_alloc(ipiv, &
         name="TriMat_piv",routine='TriMatInversion')
  end subroutine clear_BiasTriMat_inversion

end module m_ts_trimat_invert
    
