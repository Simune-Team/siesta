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
  use m_trimat_invert

  use precision, only : dp

  implicit none

  private :: dp

  ! Used for BLAS calls (local variables)
  complex(dp), private, parameter :: z0  = dcmplx( 0._dp, 0._dp)
  complex(dp), private, parameter :: z1  = dcmplx( 1._dp, 0._dp)
  complex(dp), private, parameter :: zm1 = dcmplx(-1._dp, 0._dp)

contains

  subroutine invert_TriMat_Bias(UpdateDMCR,M,Minv,no)
    use class_zTriMat
    
    ! If we only need the middle part of the Gf^-1
    logical, intent(in) :: UpdateDMCR
    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in) :: no
    complex(dp), pointer :: fullMinv(:)
    complex(dp), pointer :: Mpinv(:), Mp(:), XYn(:)
    complex(dp), pointer :: z(:), fz(:)

    integer :: nr, np , ip
    integer :: sNB, sN, sNC
    integer :: sPart, ePart
    integer :: GfsPart, GfePart
    integer :: sIdx, eIdx, prevIdx
    integer :: sNm1, sNp1, nf, n
    integer :: i

    nr = nrows_g(M)
    np = parts(M)

    ! We can still utilize the trimat-inversion
    if ( no > 0 ) then
       ! Left
       sPart = which_part(M,1)
       ePart = which_part(M,no)
    else
       ! Right
       sPart = which_part(M,nr+no+1)
       ePart = which_part(M,nr)
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

    ! ***** Dev notice *****
    ! this routine can be improved to not do too many calculations
    ! simply by calculating all Xn/Cn+1 and Yn/Bn-1 and
    ! explicitly only calculating the matrix-elements required
    !  ******

    ! We calculate the inverted matrix in the parts requested
    call invert_TriMat(M,Minv, &
         sPart = sPart, ePart = ePart)

    ! all needed quantities should still be available
    ! We have Xn/Cn+1 and Yn/Bn-1 in the needed parts

    ! Lets note that the first two parts MUST (they are forced to be)
    ! equal or exceed no.
    ! Otherwise the full \Sigma_L/R can not be contained in the tri-matrix
    ! This means that it is safe to correct the first calculated positions

    ! obtain the full inverted array values
    fullMinv => val(Minv)

    if ( no > 0 .and. abs(no) > nrows_g(M,1) ) then

       ! If we are in the left region we move about the 2,1 and 2,2
       ! note that 1,1 and 1,2 are placed correctly in the array
       ! however, 1,2 can be too large 
       
#ifdef TRANSIESTA_DEBUG
       write(*,*) 'Moving left region 2'
#endif

       ! Copy over M_{2,1}^{-1}, otherwise we will overwrite it
       Mpinv => val(Minv,2,1)
       z     => val(M,2,1)
       z(:)  =  Mpinv(:)

       ! obtain the actual indices for the second values
       call TriMat_Bias_idxs(M,no,1,sIdx,eIdx)
       ! This is the placement for the M_{1,2}^-1 array
       Mpinv => fullMinv(sIdx+nrows_g(M,1)**2:eIdx)

       ! The M_{1,2}^-1 array
       z => val(Minv,1,2)

       ! obtain the actual indices for the second values
       call TriMat_Bias_idxs(M,no,2,sIdx,eIdx)
       Mp => fullMinv(sIdx:eIdx)

       ! retain the values in M_{2,1}^{-1}
       XYn => val(M,2,1)
       
       ! we know that nrows_g(1) + nrows_g(2) >= no so the 
       ! remaining filling must be smaller than or equal to no
       ! hence, size(Mpinv) <= size(XYn)
       do i = 1 , size(Mpinv)
          ! copy over M_{1,2}^-1
          Mpinv(i) = z(i)
          ! copy over M_{2,1}^-1
          Mp(i)    = XYn(i)
       end do
       ! now M_{1,2}^{-1} lies in the correct place
       ! retain the values in M_{2,2}^{-1}
       z => val(Minv,2,2)
       sN  = size(Mpinv)
       sNC = size(XYn)
       do i = 1 , size(XYn) - size(Mpinv)
          ! copy over M_{2,1}^-1
          Mp(sN+i)  = XYn(sN+i)
          ! copy over M_{2,2}^-1
          Mp(sNC+i) = z(i)
       end do
       ! now M_{2,1}^{-1} lies in the correct place
       do i = size(XYn) - size(Mpinv) + 1, size(Mp)
          ! copy over M_{2,2}^-1
          Mp(sNC+i) = z(i)
       end do
       ! now M_{2,2}^{-1} lies in the correct place

    else if ( no < 0 .and. abs(no) > nrows_g(M,np) ) then

       ! If we are in the right region we move about the np-1,np and np-1,np-1
       ! note that np,np and np,np-1 are placed correctly in the array
       ! however, np,np-1 can be too large 

#ifdef TRANSIESTA_DEBUG
       write(*,*) 'Moving right region np-1'
#endif

       ! Copy over M_{np-1,np}^{-1}, otherwise we will overwrite it
       Mpinv => val(Minv,np-1,np)
       z     => val(M,np-1,np)
       z(:)  =  Mpinv(:)

       ! obtain the actual indices for the second values
       call TriMat_Bias_idxs(M,no,np,sIdx,eIdx)
       ! This is the placement for the M_{np,np-1}^-1 array
       Mpinv => fullMinv(eIdx-nrows_g(M,np)**2:sIdx:-1)

       ! The M_{np,np-1}^-1 array
       fz => val(Minv,np,np-1)
       z  => fz(size(fz):1:-1)

       ! obtain the actual indices for the second values
       call TriMat_Bias_idxs(M,no,np-1,sIdx,eIdx)
       Mp => fullMinv(eIdx:sIdx:-1)

       ! retain the values in M_{np-1,np}^{-1}
       fz => val(M,np-1,np)
       XYn => fz(size(fz):1:-1)
       
       ! we know that nrows_g(np-1) + nrows_g(np) >= no so the 
       ! remaining filling must be smaller than or equal to no
       ! hence, size(Mpinv) <= size(XYn)
       do i = 1 , size(Mpinv)
          ! copy over M_{np,np-1}^-1 
          Mpinv(i) = z(i)
          ! copy over M_{np-1,np}^-1
          Mp(i)    = XYn(i)
       end do
       ! now M_{np-1,np}^{-1} lies in the correct place
       ! retain the values in M_{np-1,np-1}^{-1}
       fz => val(Minv,np-1,np-1)
       z  => fz(size(fz):1:-1)
       sN  = size(Mpinv)
       sNC = size(XYn)
       do i = 1, size(XYn) - size(Mpinv)
          ! copy over M_{np-1,np}^-1
          Mp(sN+i)  = XYn(sN+i)
          ! copy over M_{np-1,np-1}^-1
          Mp(sNC+i) = z(i)
       end do
       ! now M_{np-1,np}^{-1} lies in the correct place
       do i = size(XYn) - size(Mpinv) + 1, size(Mp)
          ! copy over M_{np-1,np-1}^-1
          Mp(sNC+i) = z(i)
       end do
       ! now M_{np-1,np-1}^{-1} lies in the correct place
      
    end if

    ! we have now corrected the inverted matrix to only contain the first 
    !   no * nrows_g(tri,1+2)
    ! states. The rest is not needed

    ! We now only need to calculate the remaning part
    ! first we copy over the Xn/Cn+1 or the Yn/Bn-1 arrays
    if ( no > 0 ) then

       do ip = ePart , min(np - 1,GfePart)
          z  => Xn_div_Cn_p1(Minv,ip)
          fz => Xn_div_Cn_p1(M,ip)
          fz(:) = z(:)
       end do

    else

       do ip = sPart , max(GfsPart,2) , -1
          z  => Yn_div_Bn_m1(Minv,ip)
          fz => Yn_div_Bn_m1(M,ip)
          fz(:) = z(:)
       end do

    end if
       
    ! Now we are ready to do the calculation of the column

    if ( no > 0 ) then

       sNB = nrows_g(M,1)
       
       do ip = ePart + 1 , GfePart
      
          sN   = nrows_g(M,ip-1)
          sNp1 = nrows_g(M,ip)
          XYn => Xn_div_Cn_p1(M,ip-1)

          ! obtain the indices for the inverted matrix above
          call TriMat_Bias_idxs(M,no,ip-1,sIdx,eIdx)
          fz => fullMinv(sIdx:eIdx) ! fz is the part x no array

          ! obtain the indices for the current inverted matrix
          call TriMat_Bias_idxs(M,no,ip,sIdx,eIdx)
          z  => fullMinv(sIdx:eIdx) ! z is the part x no array

          sNC = sNB
          
          ! First we calculate the easy thing (the boundary inverted matrix)
          call zgemm('N','N',sNp1,sNC,sN, &
               zm1, XYn,sNp1, fz(1),sN,z0, z(1),sNp1)

          ! the second dimension changes to only the part we miss 
          sNC = abs(no) - sNB

          if ( sNC > 0 ) then

             ! Calculate the remaining size of the array,
             ! i.e. Gf[...,nrows_g(M,1)+1:no]
             ! calculate the offset in the full matrices for the placement
             n  = sNp1 * sNC + 1
             nf = sN   * sNC + 1

             ! this calculates the column next to part 1
             call zgemm('N','N',sNp1,sNC,sN, &
                  zm1, XYn,sNp1, fz(nf),sN,z0, z(n),sNp1)
          end if

       end do

    else

       sNB = nrows_g(M,np)

       do ip = sPart - 1 , GfsPart , -1

          sNm1 = nrows_g(M,ip)
          sN   = nrows_g(M,ip+1)
          XYn => Yn_div_Bn_m1(M,ip+1)

          ! obtain the indices for the inverted matrix below
          call TriMat_Bias_idxs(M,no,ip+1,sIdx,eIdx)
          fz => fullMinv(sIdx:eIdx) ! fz is the part x no array

          ! obtain the indices for the current inverted matrix
          call TriMat_Bias_idxs(M,no,ip,sIdx,eIdx)
          z  => fullMinv(sIdx:eIdx) ! z is the part x no array

          ! Here we calculate the Gf[...,-no+1:-sNB]
          sNC = abs(no) - sNB
         
          if ( sNC > 0 ) then
             call zgemm('N','N',sNm1,sNC,sN, &
                  zm1, XYn,sNm1, fz(1),sN,z0, z(1),sNm1)
          end if

          ! Calculate offset for the matrix
          n  = sNm1 * sNC + 1
          nf = sN   * sNC + 1

          ! the number of rows in the boundary 
          sNC = sNB

          ! First we calculate the easy thing (the boundary inverted matrix)
          call zgemm('N','N',sNm1,sNC,sN, &
               zm1, XYn,sNm1, fz(nf),sN,z0, z(n),sNm1)

       end do

    end if
    
  end subroutine invert_TriMat_Bias

  subroutine TriMat_Bias_idxs(M,no,p,sIdx,eIdx)
    use class_zTriMat
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
    use class_zTriMat
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

end module m_ts_trimat_invert
    
