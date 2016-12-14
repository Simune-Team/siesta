!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
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

  public :: ts_pivot_tri_sort
  public :: ts_pivot_tri_sort_El

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
#ifdef TBTRANS
    ! The work size needed for the electrodes
    ! depends on the size of the connecting region
    ! Hence we need to check the connecting region
    ! size
    do io = 1 , N_Elec
       no_max = max(no_max,Elecs(io)%o_inD%n)
    end do
#endif

    ! We just need to find the maximum overlap of
    ! two regions.
    ! This will give the "pushed" number of elements
    ! that is required to prevent overwriting two
    ! quantities.
    ! At minimum this will most likely be the size
    ! of the last two parts due to that being calculated
    ! last.

    els = nnzs_tri_i8b(N_tri,tri)
    worksize = maxval(tri(:)) * no_max

    ! subtract total column size
    ! to get the first matrix element of the current processing
    ! block (with an index shift of 1, so actually previous element
    ! of what is needed)
    tn = els - sum(tri(:)) * no_max
    ! Initialize the padding
    padding = 0

    cur_n = 0
    do n = 1 , N_tri

       ! Calculate the last element needed to be written for
       ! the triple product using the n'th block of
       ! the Green function
       cur_n = cur_n + tri(n) ** 2
       if ( n < N_tri ) then
          ! We fill both the n+1,n and n,n+1 blocks
          ! using the n'th block of the Green function column
          cur_n = cur_n + 2 * tri(n) * tri(n+1)
       end if

       if ( cur_n > tn ) then
          ! we have an overlap, calculate overlap
          ! and correct tn
          io = cur_n - tn
          ! We _have_ to extend the padding if previous blocks
          ! have added padding
          padding = padding + io
          ! We correct the starting index of the column matrix
          tn = tn + io
       end if
       
       ! update the placement of the retained column matrix
       ! We are now allowed to overwrite column[tri(n)]
       tn = tn + tri(n) * no_max

    end do

    if ( padding < tri(N_tri) * no_max ) then

       ! Something seems to have gone wrong in the padding-query
       ! estimation (most likely a bug)

       call die('bug in padding estimation')

    end if

  end subroutine GFGGF_needed_worksize

  subroutine needed_mem(IsVolt, N_Elec, Elecs, N_tri, tri, worksize)

    use m_ts_electype

    logical, intent(in) :: IsVolt
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(in) :: N_tri
    integer, intent(in) :: tri(N_tri)
    integer, intent(out) :: worksize

#ifdef TRANSIESTA_GFGGF_COLUMN
    integer :: pad, n
#endif
    
    ! find at which point they will cross...
    ! worksize for ONE array
    worksize = nnzs_tri(N_tri,tri)

#ifdef TRANSIESTA_GFGGF_COLUMN
    if ( IsVolt ) then
       call GFGGF_needed_worksize(N_tri, tri, &
            N_Elec, Elecs, pad, n)
       worksize = worksize + pad + n
    end if
#endif
    
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


  ! This routine sorts the orbitals/atoms within each
  ! block to create the most consecutive blocks.
  subroutine ts_pivot_tri_sort(pvt,tri)

    use intrinsic_missing, only: SORT_QUICK
    use m_region

    ! The pivoting region
    ! We know that any pivoting of orbitals
    ! within each region will not result in any
    ! changes to the tri-diagonal matrix
    ! Hence we can swap indices within each tri-diagonal
    ! part, at will.
    type(tRgn), intent(inout) :: pvt
    ! The tri-diagonal parts
    type(tRgn), intent(in) :: tri

    ! Local variables
    integer :: i, n, off

    ! First check that the sizes are the same
    if ( pvt%n /= sum(tri) ) then
       call die('Pivoting sort cannot be performed &
            &with non-equal sizes of pivoting table &
            &and block tri-diagonal matrix.')
    end if

    ! The beauty of the tri-diagonal matrix
    ! is that we can sort each of the
    ! blocks to make the memory more consecutive.
    ! This will also partly ensure electrodes
    ! are more or less streamlined
    off = 0
    do i = 1 , tri%n

       n = tri%r(i)
       call sort_quick(n, pvt%r(off+1:))

       off = off + n
       
    end do

  end subroutine ts_pivot_tri_sort


  ! This routine sorts the orbitals/atoms within each
  ! block to create the most consecutive blocks.
  ! It also sorts each block to create the most consecutive
  ! electrode regions as possible.
  subroutine ts_pivot_tri_sort_el(pvt,N_Elec,Elecs,tri)

    use m_region
    use m_ts_electype

    ! The pivoting region
    ! We know that any pivoting of orbitals
    ! within each region will not result in any
    ! changes to the tri-diagonal matrix
    ! Hence we can swap indices within each tri-diagonal
    ! part, at will.
    type(tRgn), intent(inout) :: pvt
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    ! The tri-diagonal parts
    type(tRgn), intent(in) :: tri

    ! Local variables
    type(tRgn) :: r, rp, rsub
    integer :: iE, i, n1, n2, n, off

    ! This sorts each block and checks for correct size as well
    call ts_pivot_tri_sort(pvt,tri)

    do iE = 1 , N_Elec

#ifdef TBTRANS
       call rgn_assoc(r,Elecs(iE)%o_inD)
#else
       call rgn_range(r,Elecs(iE)%idx_o, &
            Elecs(iE)%idx_o + TotUsedOrbs(Elecs(iE))-1)
#endif

       ! Figure out the parts of the electrode
       ! in the tri-diagonal matrix
       n1 = huge(1)
       n2 = 0
       do i = 1 , r%n
          n = which_part(tri,rgn_pivot(pvt,r%r(i)))
          n1 = min(n1,n)
          n2 = max(n,n2)
       end do

       ! Calculate offset in the tri-diagonal subspace
       off = 0
       do n = 1 , n1 - 1
          off = off + tri%r(n)
       end do


       ! For each part in the tri-diagonal
       ! matrix, we push the parts
       ! to the back/front
       ! depending on whether n < tri%n/2 / tri%n/2 < n
       do n = n1 , n2 ! maximum: n2-n1 == 1


       ! Create copy of tri-diagonal block
       i = tri%r(n)
       call rgn_list(rsub,i,pvt%r(off+1:off+i))

       ! Create push-list for new pivoting table
       call rgn_init(rp,i)
       rp%n = 0 ! init for zero elements

       if ( n <= tri%n / 2 ) then
          ! push to back of tri-diagonal part

          ! First add all electrode elements
          do i = 1 , r%n
             if ( in_rgn(rsub,r%r(i)) ) then
                if ( .not. rgn_push(rp,r%r(i)) ) then
                   call die('Error in programming(A) 1')
                end if
             end if
          end do

          ! Last add the remaning elements in the block
          do i = 1 , rsub%n
             if ( .not. in_rgn(rp,rsub%r(i)) ) then
                if ( .not. rgn_push(rp,rsub%r(i)) ) then
                   call die('Error in programming(A) 2')
                end if
             end if
          end do
          
       else
          ! push to front of tri-diagonal part

          ! First add the elements not being the electrode
          ! in the block
          do i = 1 , rsub%n
             if ( .not. in_rgn(r,rsub%r(i)) ) then
                if ( .not. rgn_push(rp,rsub%r(i)) ) then
                   call die('Error in programming(B) 1')
                end if
             end if
          end do

          ! Last add all electrode elements
          do i = 1 , r%n
             if ( in_rgn(rsub,r%r(i)) ) then
                if ( .not. rgn_push(rp,r%r(i)) ) then
                   call die('Error in programming(B) 2')
                end if
             end if
          end do

       end if

       ! Check that we have populated the full tri-diagonal
       ! block
       if ( rp%n /= rsub%n ) then
          call die('Error in programming 3')
       end if
       
       ! Copy back the pivoting table 
       pvt%r(off+1:off+tri%r(n)) = rp%r(:)

       ! Temporary clean
       call rgn_delete(rsub,rp)

       ! Update offset for next part
       off = off + tri%r(n)

       end do
       
    end do

#ifndef TBTRANS
    call rgn_delete(r)
#endif
       
  contains

    function which_part(tri,i) result(n)
      type(tRgn), intent(in) :: tri
      integer, intent(in) :: i
      integer :: n
      integer :: j, t

      n = 0
      
      t = 0
      do j = 1 , tri%n
         t = t + tri%r(j)
         if ( i <= t ) then
            n = j
            return
         end if
      end do
      
    end function which_part
    
  end subroutine ts_pivot_tri_sort_el

end module m_ts_tri_common
