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
! Nick Papior Andersen, 2015, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

module m_ts_tri_common

  use precision, only : dp, i8b

  implicit none

  private

  public :: GFGGF_needed_worksize
  public :: needed_mem
  public :: nnzs_tri, nnzs_tri_i8b

contains

  subroutine GFGGF_needed_worksize(N_tri, tri, &
       N_Elec, Elecs, padding, worksize)

    use m_ts_electype

    integer, intent(in) :: N_tri
    integer, intent(in) :: tri(N_tri)
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(out) :: padding, worksize

    ! els "could" be very large
    integer(i8b) :: els
    integer :: n, tn, io
    integer :: cur_n, no_max

    ! calculate the maximum electrode size
    no_max = maxval(TotUsedOrbs(Elecs))

    ! We just need to find the maximum overlap of
    ! two regions.
    ! This will give the "pushed" number of elements
    ! that is required to prevent overwriting two
    ! quantities.
    ! At minimum this will most likely be the size
    ! of the last two parts due to that being calculated
    ! last.

    els = nnzs_tri(N_tri,tri)
    worksize = maxval(tri(:)) * no_max

    ! subtract total column size
    ! to get the first matrix element of the current processing
    ! block (with an index shift of 1, so actually previous element
    ! of what is needed)
    tn = els - sum(tri(:)) * no_max

    cur_n = 0
    io = 1
    do n = 1 , N_tri

       if ( 1 < n ) &
            cur_n = cur_n + tri(n-1) * tri(n)
       cur_n = cur_n + tri(n) ** 2
       if ( n < N_tri ) &
            cur_n = cur_n + tri(n) * tri(n+1)
       
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
          tn = tn + tri(n-1) * no_max
       end if

    end do
    tn = tn + tri(N_tri) * no_max

    ! the padding must be the excess size we have appended to the matrix
    padding = tn - els

  end subroutine GFGGF_needed_worksize

  subroutine needed_mem(IsVolt, N_Elec, Elecs, N_tri, tri, worksize)

    use m_ts_electype

    logical, intent(in) :: IsVolt
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(in) :: N_tri
    integer, intent(in) :: tri(N_tri)
    integer, intent(out) :: worksize

    integer :: pad, n

    ! find at which point they will cross...
    ! worksize for ONE array
    worksize = nnzs_tri(N_tri,tri)

    if ( IsVolt ) then
       call GFGGF_needed_worksize(N_tri, tri, &
            N_Elec, Elecs, pad, n)
       worksize = worksize + pad + n
    end if
    
  end subroutine needed_mem

  function nnzs_tri(N_tri,tri) result(elem)
    integer, intent(in) :: N_tri, tri(N_tri)
    integer :: elem, i
    
    elem = tri(N_tri)**2
    do i = 1 , N_tri - 1
       elem = elem + tri(i)*( tri(i) + 2 * tri(i+1) )
    end do
    
  end function nnzs_tri

  function nnzs_tri_i8b(N_tri,tri) result(elem)
    integer, intent(in) :: N_tri, tri(N_tri)
    integer(i8b) :: elem
    integer :: i
    elem = tri(N_tri)**2
    do i = 1 , N_tri - 1
       elem = elem + tri(i)*( tri(i) + 2 * tri(i+1) )
    end do
    
  end function nnzs_tri_i8b
  
end module m_ts_tri_common
