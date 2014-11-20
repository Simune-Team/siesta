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

module m_ts_tri_scat

  use precision, only : dp
  use m_ts_method

  implicit none

  private

  public :: GF_Gamma_GF
  public :: GFGGF_needed_worksize
  public :: ts_needed_mem
  public :: has_full_part
  public :: insert_Self_Energies

contains

  ! The problem of this routine is that we wish not to
  ! overwrite the old half-inverted matrix, that would mean
  ! that we need to do the full calculation for each electrode
  ! (which can be quite time-consuming!)
  subroutine GF_Gamma_GF(Gf_tri, El, no, calc_parts, nwork, work)

    use alloc, only : re_alloc, de_alloc

    use class_zTriMat
    use m_ts_trimat_invert, only : TriMat_Bias_idxs
    use m_ts_electype

    implicit none

! *********************
! * INPUT variables   *
! *********************
    ! The Green's function column
    type(zTriMat), intent(inout) :: Gf_tri
    type(Elec), intent(in) :: El ! contains: i (Sigma - Sigma^dagger) ^T
    integer, intent(in) :: no ! The dimension of i (Sigma - Sigma^dagger) ^T
    logical, intent(in) :: calc_parts(:)

! *********************
! * OUTPUT variables  *
! *********************
    integer, intent(in) :: nwork
    complex(dp), intent(inout) :: work(nwork)

    complex(dp), parameter :: z0  = dcmplx( 0._dp, 0._dp)
    complex(dp), parameter :: z1  = dcmplx( 1._dp, 0._dp)
    complex(dp), parameter :: zi  = dcmplx( 0._dp, 1._dp)

    ! local variables
    complex(dp), pointer :: fGf(:), Gf(:), GGG(:)
    integer :: nr, np
    integer :: sIdx, eIdx
    integer :: ip, cp, n
    integer :: sN, sNc

    integer :: lsPart, lePart
    integer :: BsPart, BePart

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE GFGammaGF' )
#endif

    call timer("GFGGF",1)

    ! tri-diagonal parts information
    nr = nrows_g(Gf_tri)
    np = parts(Gf_tri)

    ! Which parts are needed
    ! In this case we need to calculate till the end
    do n = 1 , np
       if ( calc_parts(n) ) then
          lsPart = max(1,n-1)
          exit
       end if
    end do
    do n = np , 1 , -1
       if ( calc_parts(n) ) then
          lePart = min(n+1,np)
          exit
       end if
    end do

    ! Capture the full elements
    fGf => val(Gf_tri)

    ip = 0
    do n = lsPart , lePart
       ip = max(ip,no * nrows_g(Gf_tri,n))
    end do

    if ( nwork < ip ) then
       print *,nwork,ip
       call die('Work size not big enough')
    end if

    do n = lsPart , lePart

       ! Calculate the \Gamma Gf^\dagger n,1
       sN = nrows_g(Gf_tri,n)

       ! correct to the quantities that is available
       BsPart = max(n-1,lsPart)
       BePart = min(n+1,lePart)

       if ( .not. calc_parts(n) ) cycle

       call TriMat_Bias_idxs(Gf_tri,no,n,sIdx,eIdx)
       ! obtain the Gf in the respective column
       Gf => fGf(sIdx:eIdx)
#ifdef USE_GEMM3M
       call zgemm3m( &
#else
       call zgemm( &
#endif
            'T','C',no,sN,no, zi, El%Gamma, no, &
            Gf, sN, z0, work, no)
       
       ! Now we are ready to perform the multiplication
       ! for the requested region

#ifdef TRANSIESTA_DEBUG
       write(*,'(a,2(tr1,i0),a,2(tr1,i0))')'GfGGf at:',BsPart,n,' --',BePart,n
#endif
       
       ! this will populate in ascending column major order
       do cp = BsPart , BePart

          ! skip unneeded elements
          if ( .not. calc_parts(cp) ) cycle

          sNc = nrows_g(Gf_tri,cp)

          ! Retrieve Gf block
          call TriMat_Bias_idxs(Gf_tri,no,cp,sIdx,eIdx)
          Gf => fGf(sIdx:eIdx)

          ! retrieve the GGG block
          GGG => val(Gf_tri,cp,n)

          ! We need only do the product in the closest
          ! regions (we don't have information anywhere else)
#ifdef USE_GEMM3M
          call zgemm3m( &
#else
          call zgemm( &
#endif
               'N','N', sNc, sN, no, z1, &
               Gf, sNc, work, no, z0, GGG, sNc)
          
       end do
       
    end do
       
    call timer('GFGGF',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS GFGammaGF' )
#endif

  end subroutine GF_Gamma_GF

  subroutine GFGGF_needed_worksize(N_tri_part, tri_parts, &
       no_max, padding, worksize)
    integer, intent(in) :: N_tri_part
    integer, intent(in) :: tri_parts(N_tri_part)
    integer, intent(in) :: no_max
    integer, intent(out) :: padding, worksize

    integer :: els, n, tn, io
    integer :: cur_n

    ! We just need to find the maximum overlap of
    ! two regions.
    ! This will give the "pushed" number of elements
    ! that is required to prevent overwriting two
    ! quantities.
    ! At minimum this will most likely be the size
    ! of the last two parts due to that being calculated
    ! last.

    els = tri_parts(N_tri_part)**2
    worksize = tri_parts(N_tri_part)
    do n = 1 , N_tri_part - 1
       els = els + tri_parts(n)*( tri_parts(n) + 2 * tri_parts(n+1) )
       worksize = max(worksize,tri_parts(n))
    end do
    worksize = worksize * no_max

    ! subtract total column size
    ! to get the first matrix element of the current processing
    ! block (with an index shift of 1, so actually previous element
    ! of what is needed)
    tn = els - sum(tri_parts(1:n) * no_max)

    cur_n = 0
    io = 1
    do n = 1 , N_tri_part

       if ( 1 < n ) &
            cur_n = cur_n + tri_parts(n-1) * tri_parts(n)
       cur_n = cur_n + tri_parts(n) ** 2
       if ( n < N_tri_part ) &
            cur_n = cur_n + tri_parts(n) * tri_parts(n+1)
       
       if ( cur_n > tn ) then
          ! we have an overlap, calculate overlap
          ! and correct tn
          ! With ">" we do not need to correct the tn initialization
          ! of element - 1 as noted above
          padding = cur_n - tn
          ! We correct the starting index of tn
          tn = tn + padding
       end if
       
       ! we need to retain the column block
       ! for the next block...
       ! in that way we can still multiply the previous
       ! block with the current block.
       if ( n > 1 ) then
          tn = tn + tri_parts(n-1) * no_max
       end if

    end do
    tn = tn + tri_parts(N_tri_part) * no_max

    ! the padding must be the excess size we have appended to the matrix
    padding = tn - els

  end subroutine GFGGF_needed_worksize

  subroutine ts_needed_mem(IsVolt, N_Elec, Elecs, N_tri_part, tri_parts, worksize)
    use m_ts_electype
    logical, intent(in) :: IsVolt
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(in) :: N_tri_part
    integer, intent(in) :: tri_parts(N_tri_part)
    integer, intent(out) :: worksize

    integer :: pad, n, no

    ! find at which point they will cross...
    worksize = tri_parts(N_tri_part)**2
    do n = 1 , N_tri_part - 1
       worksize = worksize + &
            tri_parts(n)*( tri_parts(n) + 2 * tri_parts(n+1) )
    end do
    ! for the two arrays
    worksize = worksize * 2

    if ( IsVolt ) then
       no = maxval(TotUsedOrbs(Elecs(:)))
       call GFGGF_needed_worksize(N_tri_part, tri_parts, &
            no, pad, n)
       worksize = worksize + pad + n
    end if
    
  end subroutine ts_needed_mem

  function has_full_part(N_tri_part,tri_parts, &
       part,io1,io2) result(has)
    integer, intent(in) :: N_tri_part, tri_parts(N_tri_part), part, io1, io2
    logical :: has
    integer :: i, io

    io = 1
    do i = 1 , part - 1
       io = io + tri_parts(i)
    end do
    
    has = io1 <= io .and. &
         io + tri_parts(i) - 1 <= io2

  end function has_full_part

  ! Generic routine for inserting the self-energies in the 
  ! tri-diagonal matrices
  subroutine insert_Self_Energies(no_u, Gfinv_tri, Gfinv, El)
    use m_ts_electype
    use class_zTriMat
    integer, intent(in) :: no_u
    type(zTriMat), intent(inout) :: GFinv_tri
    complex(dp), intent(inout) :: Gfinv(:)
    type(Elec), intent(in) :: El

    integer :: no, off, i, j, ii, idx
    
    no = TotUsedOrbs(El)
    off = El%idx_o - orb_offset(El%idx_o) - 1

    if ( El%Bulk ) then
!$OMP do private(j,i,idx,ii)
       do j = 0 , no - 1
          ii = j * no
          do i = 1 , no
             idx = index(GFinv_tri,i+off,j+1+off)
             Gfinv(idx) = El%Sigma(ii+i)
          end do
       end do
!$OMP end do nowait
    else
!$OMP do private(j,i,idx,ii)
       do j = 0 , no - 1
          ii = j * no
          do i = 1 , no
             idx = index(GFinv_tri,i+off,j+1+off)
             Gfinv(idx) = Gfinv(idx) - El%Sigma(ii+i)
          end do
       end do
!$OMP end do nowait
    end if
    
  end subroutine insert_Self_Energies

end module m_ts_tri_scat
