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

! Here we implement routines that enables the reduction of the bandwidth of 
! a sparse matrix

! Currently we have implemented 3 different algorithms
!  1) Cuthill-Mckee
!  2) Cuthill-Mckee with z-position priority
!  3) Papior algorithm (intuitive algoritm)

! In all cases the reverse algorithm (reduces fill-in) can be 
! performed by rotating the input matrix by 180 degrees

module m_bandwidth

  implicit none

  public 

  integer, parameter :: BW_CUTHILL_MCKEE = 1
  integer, parameter :: BW_CUTHILL_MCKEE_Z_PRIORITY = 2
  integer, parameter :: BW_PAPIOR = 3

  ! This integer is added to create the reverse algorithm
  integer, parameter :: BW_REVERSE = 1000

contains

  subroutine BandWidth_pivoting(method,na_u,a_mm,R,na_L,na_R,priority)
    integer, intent(in) :: method, na_u
    integer, intent(inout) :: a_mm(na_u,na_u)
    integer, intent(out) :: R(na_u)
    integer, intent(in), optional :: na_L, na_R
    integer, intent(inout), optional :: priority(na_u)
    integer :: i, ia, iac, lL, lR, tmp, band_algo

    lL = 0
    if ( present(na_L) ) lL = na_L
    lR = 0
    if ( present(na_R) ) lR = na_R

    band_algo = method

    if ( method >= BW_REVERSE ) then
       band_algo = method - BW_REVERSE
       ! Reversing means that the ending atoms are switched
       if ( present(na_L) ) lR = na_L
       if ( present(na_R) ) lL = na_R

       ! We need to rotate the a_mm array by 180 degrees.
       ! This will enable us to do the reverse algorithm without 
       ! altering the routines
       a_mm = CSHIFT(a_mm,SHIFT=na_u,DIM=1)
       a_mm = CSHIFT(a_mm,SHIFT=na_u,DIM=2)

       ! Reverse the priority
       if ( present(priority) ) then
          do i = 1 , na_u
             priority(na_u+1-i) = priority(i)
          end do
       end if

    end if

    ! Initialize the pivoting array
    do i = 1 , na_u
       R(i) = i
    end do

    select case ( band_algo )
    case ( BW_CUTHILL_MCKEE )
       call Cuthill_Mckee(na_u,lL,lR,a_mm,R)
    case ( BW_CUTHILL_MCKEE_Z_PRIORITY )
       call Cuthill_Mckee(na_u,lL,lR,a_mm,R, priority = priority)
    case ( BW_PAPIOR )
       call Papior(na_u,lL,lR,a_mm,R, priority = priority)
    case default
       call die('Could not recognize the bandwidth reduction algorithm')
    end select

    if ( method >= BW_REVERSE ) then

       ! Revert to correct atomic indices
       do ia = lL + 1 , na_u - lR
          R(ia) = na_u + 1 - R(ia)
       end do

       ! transfer to the reversed Cuthill-Mckee algorithm
       ! We could do this and the above in one go, but that would require a 
       ! correction if mod(io,2) == 1
       i = na_u - lL - lR
       do ia = lL + 1 , lL + i / 2
          iac = R(ia)
          R(ia) = R(na_u+1-ia)
          R(na_u+1-ia) = iac
       end do

    end if
    
  end subroutine BandWidth_pivoting

  subroutine Cuthill_Mckee(na_u,na_L,na_R,a_mm,R,priority)
    ! We return the pivoting indexes for the atoms
    use parallel, only : IONode
    use alloc
    ! number of atoms in the cell, also the last orbitals of each
    ! atom.
    integer, intent(in) :: na_u
    ! The atoms which need not be accounted for (na_BufL/R + na_L/R)
    integer, intent(in) :: na_L, na_R
    integer, intent(in) :: a_mm(na_u,na_u)
    integer, intent(inout) :: R(na_u)
    integer, intent(in), optional :: priority(na_u)

    integer, pointer :: Q(:)
    integer :: na_S, na_E, na
    integer :: ia, nQ, iQ, i_ca, ca
    integer :: degree

    na = na_u - na_L - na_R
    na_S = na_L + 1
    na_E = na_u - na_R

    nullify(Q)
    call re_alloc(Q,1,na)

    ! we can now perform the Cuthill-Mckee algorithm
    ! Notice that we force the last atom in the electrode
    ! to be in the result array (thus we force
    ! the construction of the Cuthill-Mckee to start
    ! with the electrode)

    ! Prepare counters for the algorithm
    i_ca = na_L

    ! the counter for the number of saved atoms
    ca = 0

    ! start the algorithm 
    do while ( ca < na )

       ! initialize the current node
       ! remark that the first iteration fixes
       ! the node building from the last atom
       ! specified
       nQ = 0
       if ( i_ca > 0 ) then
          ! In case we want a pure Cuthill-Mckee
          ! we shall not construct the queue the first
          ! time
          call add_Queue(na_u,a_mm,na_S,na_E,i_ca,Q,nQ, &
               priority=priority)
       end if

       ! loop the queed items
       iQ = 0
       ! if the number of connected nodes is zero
       ! we skip (the node is already added)
       do while ( iQ < nQ ) 
          iQ = iQ + 1 
          if ( ca == 0 ) then
             ! Force the addition to the result list
             ca = ca + 1
             R(na_L+ca) = Q(iQ)
             call add_Queue(na_u,a_mm,na_S,na_E,Q(iQ),Q,nQ, &
                  priority=priority)
          else if ( any( R(na_L+1:na_L+ca) == Q(iQ) ) ) then
             ! Just delete the entry (we do not need to worry about
             ! it anymore)
             Q(iQ:nQ-1) = Q(iQ+1:nQ)
             iQ = iQ - 1
             nQ = nQ - 1
          else
             ! We are allowed to add the queued item to the result list
             ca = ca + 1
             R(na_L+ca) = Q(iQ)
             ! update the queue-list
             call add_Queue(na_u,a_mm,na_S,na_E,Q(iQ),Q,nQ, &
                  priority=priority)
          end if
       end do

       ! Do a quick exit if possible
       if ( ca == na ) cycle

       ! From here on it is the actual Cuthill-Mckee algorithm
       ! Now we need to select the one with the lowest degree

       i_ca = 0
       degree = huge(1)
       do ia = na_S , na_E
          iQ = 0

          if ( ca == 0 ) then
             ! If no items has been added
             iQ = 1
          else if ( .not. any( R(na_L+1:na_L+ca) == ia ) ) then
             ! Check it hasn't been processed 
             iQ = 1
          end if

          if ( iQ == 1 ) then
             iQ = sum(a_mm(na_S:na_E,ia))
             if ( degree > iQ ) then
                
                ! We save this as the new node
                degree = iQ
                i_ca = ia
             else if ( degree == iQ .and. &
                  present(priority) .and. i_ca > 0 ) then
                if ( priority(i_ca) < priority(ia) ) then
                   degree = iQ
                   i_ca = ia
                end if

             end if
          end if
       end do
       
       ! Append to the result list
       if ( i_ca > 0 ) then
          ca = ca + 1
          R(ca) = i_ca
       end if
       
    end do

    call de_alloc(Q)
    
  contains 
    
    subroutine add_Queue(na_u,a_mm,na_S,na_E,ia,Q,nQ, priority)
      integer, intent(in) :: na_u, a_mm(na_u,na_u)
      integer, intent(in) :: na_S, na_E, ia
      integer, intent(inout) :: Q(na_E - na_S + 1)
      integer, intent(inout) :: nQ
      integer, intent(in), optional :: priority(na_u)

      integer :: i, iQ, ii, lnQ, j
      integer :: degree(nQ+1:na_u)
      logical :: append, switch

      ! We prepend to the queue
      iQ = nQ
      do i = na_S , na_E
         ! if they are connected and not the same atom
         if ( a_mm(i,ia) == 1 .and. i /= ia ) then
            append = .false.
            if ( nQ == 0 ) then
               append = .true.
            else if ( .not. any( Q(1:nQ) == i ) ) then
               append = .true.
            end if
            if ( append ) then
               ! We have a connected node
               ! save its degree
               iQ = iQ + 1 
               Q(iQ) = i
               degree(iQ) = sum(a_mm(na_S:na_E,i))
            end if
         end if
      end do
      ! capture the connected nodes
      lnQ = iQ

      ! Sort the queued items (original)
      do i = nQ + 1 , lnQ
         ! we sort by the first find the one to switch with
         do ii = lnQ , i + 1 , -1
            switch = (degree(ii) < degree(i))

            if ( present(priority) ) then
               if ( degree(ii) == degree(i) ) &
                    switch = ( priority(ii) > priority(i) )
            end if

            if ( switch ) then
               iQ = degree(ii)
               degree(ii) = degree(i)
               degree(i) = iQ
               iQ = Q(ii)
               Q(ii) = Q(i)
               Q(i) = iQ
            end if
         end do

      end do

      ! Save the new queued count
      nQ = lnQ

    end subroutine add_Queue

  end subroutine Cuthill_Mckee

  subroutine Papior(na_u,na_L,na_R,a_mm,R,priority)
    ! We return the pivoting indexes for the atoms
    use parallel, only : IONode
    use alloc
    ! number of atoms in the cell, also the last orbitals of each
    ! atom.
    integer, intent(in) :: na_u
    ! The atoms which need not be accounted for (na_BufL/R + na_L/R)
    integer, intent(in) :: na_L, na_R
    integer, intent(in) :: a_mm(na_u,na_u)
    integer, intent(inout) :: R(na_u)
    integer, intent(in) :: priority(na_u)

    integer, pointer :: Q(:), pri_idx(:)
    integer :: na_S, na_E, na
    integer :: ia, nQ, iQ, i_ca, ca, ja
    integer :: degree

    na = na_u - na_L - na_R
    na_S = na_L + 1
    na_E = na_u - na_R

    nullify(Q,pri_idx)
    call re_alloc(Q,1,na)
    call re_alloc(pri_idx,1,na)

    ! Notice that we force the last atom in the electrode
    ! to be in the result array (thus we force
    ! the construction of the Papior to start
    ! with the electrode)

    ! Prepare counters for the algorithm
    i_ca = na_L

    ! the counter for the number of saved atoms
    ca = 0

    ! start the algorithm 
    do while ( ca < na )

       ! initialize the current node
       ! remark that the first iteration fixes
       ! the node building from the last atom
       ! specified
       nQ = 0
       if ( i_ca > 0 ) then
          ! In case we want a pure Papior
          ! we shall not construct the queue the first
          ! time
          call add_Queue(na_u,a_mm,na_S,na_E,i_ca,Q,nQ, &
               priority=priority)
       end if

       ! loop the queed items
       iQ = 0
       ! if the number of connected nodes is zero
       ! we skip (the node is already added)
       do while ( iQ < nQ ) 
          iQ = iQ + 1 
          if ( ca == 0 ) then
             ! Force the addition to the result list
             ca = ca + 1
             R(na_L+ca) = Q(iQ)
             call add_Queue(na_u,a_mm,na_S,na_E,Q(iQ),Q,nQ, &
                  priority=priority)
          else if ( any( R(na_L+1:na_L+ca) == Q(iQ) ) ) then
             ! Just delete the entry (we do not need to worry about
             ! it anymore)
             Q(iQ:nQ-1) = Q(iQ+1:nQ)
             iQ = iQ - 1
             nQ = nQ - 1
          else
             ! We are allowed to add the queued item to the result list
             ca = ca + 1
             R(na_L+ca) = Q(iQ)
             ! update the queue-list
             call add_Queue(na_u,a_mm,na_S,na_E,Q(iQ),Q,nQ, &
                  priority=priority)
          end if
       end do

       ! Do a quick exit if possible
       if ( ca == na ) cycle

       ! From here on it is the actual Papior algorithm
       ! Now we need to select the one with the lowest degree

       i_ca = 0
       degree = huge(1)
       do ia = na_S , na_E
          iQ = 0

          if ( ca == 0 ) then
             ! If no items has been added
             iQ = 1
          else if ( .not. any( R(na_L+1:na_L+ca) == ia ) ) then
             ! Check it hasn't been processed 
             iQ = 1
          end if

          if ( iQ == 1 ) then
             do ja = ia , na_E ! from i ensures that we get at least a zero
                if ( a_mm(ja,ia) == 1 ) iQ = ia - ja
             end do
             if ( degree > iQ ) then
                
                ! We save this as the new node
                degree = iQ
                i_ca = ia
             else if ( degree == iQ .and. i_ca > 0 ) then
                if ( priority(i_ca) < priority(ia) ) then
                   degree = iQ
                   i_ca = ia
                end if

             end if
          end if
       end do
       
       ! Append to the result list
       if ( i_ca > 0 ) then
          ca = ca + 1
          R(ca) = i_ca
       end if
       
    end do

    call de_alloc(Q)
    call de_alloc(pri_idx)
    
  contains 
    
    ! Instead of the Papior we use a different
    ! sorting algorithm
    ! We take that the degree is the distance from the current 
    ! segment to the farthest segment, and then we compare that the
    ! farthest segment is the smallest
    subroutine add_Queue(na_u,a_mm,na_S,na_E,ia,Q,nQ, priority)
      integer, intent(in) :: na_u, a_mm(na_u,na_u)
      integer, intent(in) :: na_S, na_E, ia
      integer, intent(inout) :: Q(na_E - na_S + 1)
      integer, intent(inout) :: nQ
      integer, intent(in) :: priority(na_u)

      integer :: i, iQ, ii, lnQ, j
      integer :: degree(nQ+1:na_u)
      logical :: append, switch

      ! We prepend to the queue
      iQ = nQ
      do i = na_S , na_E
         ! if they are connected and not the same atom
         if ( a_mm(i,ia) == 1 .and. i /= ia ) then
            append = .false.
            if ( nQ == 0 ) then
               append = .true.
            else if ( .not. any( Q(1:nQ) == i ) ) then
               append = .true.
            end if
            if ( append ) then
               ! We have a connected node
               ! save its degree
               iQ = iQ + 1 
               Q(iQ) = i
               do j = i , na_E ! from i ensures that we get at least a zero
                  if ( a_mm(j,i) == 1 ) degree(iQ) = i - j
               end do
            end if
         end if
      end do
      ! capture the connected nodes
      lnQ = iQ

      ! Sort the queued items
      do i = nQ + 1 , lnQ
         do ii = i + 1 , lnQ
            switch = (degree(ii) < degree(i))

            if ( degree(ii) == degree(i) ) &
                 switch = ( priority(ii) > priority(i) )
               
            if ( switch ) then
               iQ = degree(ii)
               degree(ii) = degree(i)
               degree(i) = iQ
               iQ = Q(ii)
               Q(ii) = Q(i)
               Q(i) = iQ
            end if
         end do
      end do
      
      ! Save the new queued count
      nQ = lnQ

    end subroutine add_Queue
    
  end subroutine Papior

end module m_bandwidth
