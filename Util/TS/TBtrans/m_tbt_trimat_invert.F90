! Module for creating the column of tbtrans
! in the tri-diagonal matrix.

module m_tbt_trimat_invert

  use m_ts_trimat_invert, only : invert_BiasTriMat_prep
  use m_ts_trimat_invert, only : init_BiasTriMat_inversion
  use m_ts_trimat_invert, only : clear_BiasTriMat_inversion
  use m_ts_trimat_invert, only : TriMat_Bias_idxs

  implicit none

contains

  subroutine invert_BiasTriMat_col(M,Minv,r,r_col)

    use precision, only : dp

    use class_zTriMat
    use m_trimat_invert, only : Xn_div_Cn_p1, Yn_div_Bn_m1

    use m_region

    type(zTriMat), intent(inout) :: M, Minv
    type(tRgn), intent(in) :: r
    type(tRgn), intent(in) :: r_col

    complex(dp), pointer :: Mpinv(:), Mp(:)
    complex(dp), pointer :: XY(:)
    complex(dp), pointer :: z(:)

    ! Used for BLAS calls (local variables)
    complex(dp), parameter :: z0  = dcmplx( 0._dp, 0._dp)
    complex(dp), parameter :: zm1 = dcmplx(-1._dp, 0._dp)

    integer :: nr, np
    integer :: sPart, ePart
    integer :: i, ii, i_Elec, idx_Elec
    integer :: sIdxF, eIdxF, sIdxT, eIdxT
    integer :: sN, n, in, s, sNo
    integer, allocatable :: cumsum(:)
    type(tRgn) :: rB

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

    call timer('V_TM_inv',1)

    nr = nrows_g(M)
    np = parts(M)
    allocate(cumsum(np))
    cumsum(1) = 0
    do n = 2 , np
       cumsum(n) = cumsum(n-1) + nrows_g(M,n-1)
    end do

    ! This code is based on the down-folded self-energies
    ! which are determined by the col region

    sPart = huge(1)
    ePart = 0
    do n = 1 , r_col%n
       s = rgn_pivot(r,r_col%r(n))
       sPart = min(sPart,which_part(M,s))
       ePart = max(ePart,which_part(M,s))
    end do
    if ( sPart < 1 ) call die('Error in the Bias inversion, sPart')
    if ( ePart - sPart + 1 > 2 ) call die('Error in trimat partition')
    if ( ePart > np ) call die('Error in the Bias inversion, ePart')

    ! Point to the matrices
    z => val(Minv)

    ! CHECK
    ! This requires that the o_inD is sorted
    ! according to the device region.
    ! Check m_tbt_regions to assert this!

    i_Elec = 1
    rB%n = 1
    do while ( i_Elec <= r_col%n ) 

       idx_Elec = rgn_pivot(r,r_col%r(i_Elec))

       ! We start by copying over the Mnn in blocks

       ! We start by creating a region of consecutive
       ! memory.
       n = which_part(M,idx_Elec)
       sN = nrows_g(M,n)
       ii = 0
       do
          i = rgn_pivot(r,r_col%r(i_Elec+ii))
          ! In case it is not consecutive
          if ( i - idx_Elec /= ii ) exit
          ! In case the block changes, then
          ! we cut the block size here.
          if ( n /= which_part(M,i) ) exit
          ii = ii + 1
          if ( i_Elec + ii > r_col%n ) exit
       end do
       ! The consecutive memory block is this size 'ii'
       call rgn_list(rB,ii,r_col%r(i_Elec:i_Elec+ii-1))

       ! Copy over this portion of the Mnn array
       
       ! Figure out which part we have Mnn in
       i = rgn_pivot(r,rB%r(1))
       n = which_part(M,i)

       ! get placement of the diagonal block in the column
       call TriMat_Bias_idxs(Minv,r_col%n,n,sIdxT,eIdxT)
       Mpinv => z(sIdxT:eIdxT)

       ! Correct the copied elements
       ! Figure out the placement in the copied to array
       ! First we calculate the starting index of the block
       sIdxT = ( i_Elec -   1 ) * sN + 1
       eIdxT = ( i_Elec + rB%n - 1) * sN

       ! *** Now we have the matrix that we can save the 
       !     block Mnn in

       ! We now need to figure out the placement of the 
       ! Mnn part that we have already calculated.
       Mp => val(M,n,n)
       i = rgn_pivot(r,rB%r(1))
       sIdxF = (i-cumsum(n)-1) * sN + 1
       i = rgn_pivot(r,rB%r(rB%n))
       eIdxF = (i-cumsum(n)) * sN

       ! Check that we have something correct...
!print *,trim(El%name),sIdxT,eIdxT,sIdxF,eIdxF
!print *,trim(El%name),eIdxT-sIdxT,eIdxF-sIdxF

       if ( eIdxT - sIdxT /= eIdxF - sIdxF ) & 
            call die('Error in determining column')

       ! Copy over diagonal block
!$OMP parallel workshare default(shared)
       Mpinv(sIdxT:eIdxT) = Mp(sIdxF:eIdxF)
!$OMP end parallel workshare

       ! Prepare for next segment
       Mp => Mp(sIdxF:)

       ! Calculate the off-diagonal Green's function in the regions
       ! of interest
       do in = n - 1 , sPart , - 1

!print *,'calculating Mn-1n',in,n, i_Elec
          ! Number of orbitals in the other segment
          sNo = nrows_g(M,in)
          ! grab the Yn matrix to perform the 
          ! Gf off-diagonal calculation.
          XY => Yn_div_Bn_m1(M,in+1)

          ! placement of the inverted matrix
          call TriMat_Bias_idxs(Minv,r_col%n,in,sIdxT,eIdxT)
          Mpinv => z(sIdxT:eIdxT)
          sIdxT = ( i_Elec - 1 ) * sNo + 1

#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
               'N','N',sNo,rB%n,sN, &
               zm1, XY, sNo, Mp(1),sN,z0, Mpinv(sIdxT),sNo)

          Mp => Mpinv(sIdxT:)
          sN = sNo
          
       end do
          
       ! Reset arrays to just before off-diagonal
       sN = nrows_g(M,n)
       Mp => val(M,n,n)
       Mp => Mp(sIdxF:)

       ! Calculate the off-diagonal Green's function in the regions
       ! of interest
       do in = n + 1 , ePart

!print *,'calculating Mn+1n',in,n, i_Elec
          ! Number of orbitals in the other segment
          sNo = nrows_g(M,in)
          ! grab the Xn matrix to perform the 
          ! Gf off-diagonal calculation.
          XY => Xn_div_Cn_p1(M,in-1)

          ! placement of the inverted matrix
          call TriMat_Bias_idxs(Minv,r_col%n,in,sIdxT,eIdxT)
          Mpinv => z(sIdxT:eIdxT)
          sIdxT = ( i_Elec - 1 ) * sNo + 1

#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
               'N','N',sNo,rB%n,sN, &
               zm1, XY, sNo, Mp(1),sN,z0, Mpinv(sIdxT),sNo)

          Mp => Mpinv(sIdxT:)
          sN = sNo
          
       end do
       
       ! Update current segment of the electrode copied entries.
       i_Elec = i_Elec + rB%n

    end do

    ! Now we need to calculate the remaining column

    ! We now calculate:
    !  Mmn = -Ym+1/Bm * Mm+1n, for m<n
    do n = sPart - 1 , 1 , - 1
       
       sN  = nrows_g(M,n)
       sNo = nrows_g(M,n+1)

       ! get Ym+1/Bm
       XY => Yn_div_Bn_m1(M,n+1)

       ! Get Mm+1n
       call TriMat_Bias_idxs(Minv,r_col%n,n+1,sIdxF,eIdxF)
       Mp    => z(sIdxF:eIdxF)
       ! Get Mmn
       call TriMat_Bias_idxs(Minv,r_col%n,n,sIdxF,eIdxF)
       Mpinv => z(sIdxF:eIdxF)

#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',sN,r_col%n,sNo, &
            zm1, XY, sN, Mp(1),sNo,z0, Mpinv(1),sN)
       
    end do

    ! We now calculate:
    !  Mmn = -Xm-1/Cm * Mm-1n, for m>n
    do n = ePart + 1 , np

       sNo = nrows_g(M,n-1)
       sN  = nrows_g(M,n)

       ! get Xm-1/Cm
       XY => Xn_div_Cn_p1(M,n-1)
       
       ! Get Mm-1n
       call TriMat_Bias_idxs(Minv,r_col%n,n-1,sIdxF,eIdxF)
       Mp    => z(sIdxF:eIdxF)
       ! Get Mmn
       call TriMat_Bias_idxs(Minv,r_col%n,n,sIdxF,eIdxF)
       Mpinv => z(sIdxF:eIdxF)
       
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'N','N',sN,r_col%n,sNo, &
            zm1, XY, sN, Mp(1),sNo,z0, Mpinv(1),sN)
       
    end do

    call rgn_delete(rB)

    ! At this point the total 
    ! inverted column is placed at the end of
    ! the tri-mat inversion.

    call timer('V_TM_inv',2)

  end subroutine invert_BiasTriMat_col

end module m_tbt_trimat_invert
