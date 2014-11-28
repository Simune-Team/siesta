module m_supercell

! Calculates the supercell factors by examining all atoms.
! 1. We find the smallest vector that connects any atom ia1 to
!    the neighbouring cell atom ia2 [in direction xyz]
! 2. Loop on supercell factors until it is not overlapping
! 3. Save supercell factor

! This will correctly capture no supercell factors in
! slaps with only UC connections. This supercell size will typically
! be one less than the usual 1+2*ceiling(rmaxh/vnorm(ucell(:,xyz)))

  ! TODO : create a fully parallel version (each node can take na_u/Nodes
  ! atoms)
  
  use precision, only : dp
  
  implicit none
  
  public :: exact_sc_size
  private
  
contains

  subroutine exact_sc_size(rmaxh,ucell,na_u,xa,nsc)
    use intrinsic_missing, only : VNORM
    ! Calculates the exact size required to encompass all super-cells
    ! This is a brute force calculation which is over the top, yet
    ! it is effictive in reducing the super-cell to the minimum basis
    real(dp), intent(in) :: rmaxh
    real(dp), intent(in) :: ucell(3,3)
    integer,  intent(in) :: na_u
    real(dp), intent(in) :: xa(3,na_u)
    integer,  intent(inout) :: nsc(3)

    ! Local variables
    integer :: xyz, ia, ja, tnsc(3), idiag
    real(dp) :: recell(3,3), v1(3), v2(3), ucdiag(3)
    logical :: on_boundary

    call reclat(ucell,recell,0) ! do not add 2 Pi

    ! Initialize new nsc 
    tnsc(:) = 1

    ! We loop over all atoms and find the two atoms
    ! that will be connected by the least 
    ! length across the first cell boundary
    do ia = 1 , na_u
       do ja = ia , na_u
          on_boundary = (ia == ja)

          ! Vector between atoms
          v1(:) = xa(:,ja) - xa(:,ia)

          do xyz = 1 , 3

             ! from i -> J
             v2(:) = v1(:) + ucell(:,xyz)
             ! update the nsc quantity based on this vector
             call update_nsc(rmaxh,ucell(:,xyz),v2,on_boundary,tnsc(xyz))

             ! from j -> I
             v2(:) = -v1(:) + ucell(:,xyz)
             ! update the nsc quantity based on this vector
             call update_nsc(rmaxh,ucell(:,xyz),v2,on_boundary,tnsc(xyz))

             ! Possibly track the diagonal path when
             ! having skewed unit-cells. Specifically if 
             ! |ucell(:,1)|,|ucell(:,2)| > |ucell(:,1) + ucell(:,2)|
             ! Thus we should track the diagonal of each direction

             if ( xyz == 3 ) cycle
             do idiag = xyz + 1 , 3
                
                ! get the diagonal unit-cell direction
                ucdiag(:) = ucell(:,xyz) + ucell(:,idiag)

                ! from i -> J
                v2(:) = v1(:) + ucdiag(:)
                ! update the nsc quantity based on this vector (note that
                ! the diagonal now takes two into account
                call update_nsc(rmaxh,ucdiag,v2,on_boundary,tnsc(xyz))
                call update_nsc(rmaxh,ucdiag,v2,on_boundary,tnsc(idiag))

                ! from j -> I
                v2(:) = -v1(:) + ucdiag(:)
                ! update the nsc quantity based on this vector (note that
                ! the diagonal now takes two into account
                call update_nsc(rmaxh,ucdiag,v2,on_boundary,tnsc(xyz))
                call update_nsc(rmaxh,ucdiag,v2,on_boundary,tnsc(idiag))

             end do
                
          end do

          ! Do the last diagonal unit-cell direction (a+b+c)
          ucdiag(:) = ucell(:,1) + ucell(:,2) + ucell(:,3)
          
          ! from i -> J
          v2(:) = v1(:) + ucdiag(:)
          ! update the nsc quantity based on this vector (note that
          ! the diagonal now takes two into account
          call update_nsc(rmaxh,ucdiag,v2,on_boundary,tnsc(1))
          call update_nsc(rmaxh,ucdiag,v2,on_boundary,tnsc(2))
          call update_nsc(rmaxh,ucdiag,v2,on_boundary,tnsc(3))

          ! from j -> I
          v2(:) = -v1(:) + ucdiag(:)
          ! update the nsc quantity based on this vector (note that
          ! the diagonal now takes two into account
          call update_nsc(rmaxh,ucdiag,v2,on_boundary,tnsc(1))
          call update_nsc(rmaxh,ucdiag,v2,on_boundary,tnsc(2))
          call update_nsc(rmaxh,ucdiag,v2,on_boundary,tnsc(3))

       end do
    end do

    ! Copy over calculated super-cell factors
    nsc(:) = tnsc(:)

  contains

    subroutine update_nsc(rmaxh,vcell,v,on_boundary,nsc)
      real(dp), intent(in) :: rmaxh
      real(dp), intent(in) :: vcell(3)
      real(dp), intent(inout) :: v(3)
      logical, intent(in) :: on_boundary
      integer, intent(inout) :: nsc
      real(dp), parameter :: EPS = 1.e-8_dp
      real(dp) :: vl
      integer :: n

      vl = VNORM(v)

      ! we do not have any connections in the xyz direction
      ! Hence we can entirely skip that direction
      ! This will be the case for slap calculations
      if ( vl > rmaxh + EPS ) return

      n = 1
      ! Loop and find the minimum supercell factor
      do while ( vl <= rmaxh + EPS )
         ! increment supercell atom
         n = n + 1
         v(:) = v(:) + vcell(:)
         vl = vnorm(v)
      end do

      ! As vl now is larger than rmaxh we can subtract
      ! one from the number of supercells, truncate
      ! at 1 as we do have connects across the unit-cell
      ! This will typically reduce nsc by one from the
      ! naive auxillary cell as the naive calculation
      ! must take into account the two atoms lying on the
      ! boundary. But this is direct interaction! :)
      if ( .not. on_boundary ) then
         ! TODO, I think this should ALWAYS be n-1
         ! Yet a case example of unit cell:
         ! Unit cell vectors (Ang):
         !        1.414210   -2.449490    0.000000
         !        1.414210    2.449490    0.000000
         !        0.000000    0.000000    6.928200
         ! yields the wrong supercell indices if using n-1
         ! I think maybe the hsparse takes too many non-zero elements?
         ! Or are there some additional corrections for the
         ! orbital range I do not know of?
         n = max(1,n - 1)
      end if

      ! UC + n cells on both sides
      n = 1 + 2 * n
      if ( n > nsc ) nsc = n

    end subroutine update_nsc
    
  end subroutine exact_sc_size
  
end module m_supercell
