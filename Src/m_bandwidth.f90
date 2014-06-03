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
!  2) Cuthill-Mckee with transport direction priority
!  3) Papior algorithm (intuitive algoritm)

! In all cases the reverse algorithm (reduces fill-in) can be 
! performed by rotating the input matrix by 180 degrees

module m_bandwidth

  implicit none

  public 

  integer, parameter :: BW_CUTHILL_MCKEE = 1
  integer, parameter :: BW_CUTHILL_MCKEE_PRIORITY = 2
  integer, parameter :: BW_PAPIOR = 3

  ! This integer is added to create the reverse algorithm
  integer, parameter :: BW_REVERSE = 1000

contains

  subroutine BandWidth_pivoting(method,na_u,a_mm,R,priority)
    integer, intent(in) :: method, na_u
    integer, intent(inout) :: a_mm(na_u,na_u)
    integer, intent(out) :: R(na_u)
    integer, intent(inout), optional :: priority(na_u)
    integer :: i, ia, iac, tmp, band_algo

    band_algo = method

    if ( method >= BW_REVERSE ) then
       band_algo = method - BW_REVERSE

       ! We need to rotate the a_mm array by 180 degrees.
       ! This will enable us to do the reverse algorithm without 
       ! altering the routines
       do i = 1 , na_u
          do ia = 1 , na_u / 2
             tmp = a_mm(ia,i)
             a_mm(ia,i) = a_mm(na_u-ia+1,i)
             a_mm(na_u-ia+1,i) = tmp
          end do
       end do
       do ia = 1 , na_u / 2
          do i = 1 , na_u
             tmp = a_mm(i,ia)
             a_mm(i,ia) = a_mm(i,na_u-ia+1)
             a_mm(i,na_u-ia+1) = tmp
          end do
       end do

       ! Reverse the priority
       if ( present(priority) ) then
          do i = 1 , na_u / 2
             tmp = priority(na_u+1-i)
             priority(na_u+1-i) = priority(i)
             priority(i) = tmp
          end do
       end if

    end if

    ! Initialize the pivoting array
    do i = 1 , na_u
       R(i) = i
    end do

    select case ( band_algo )
    case ( BW_CUTHILL_MCKEE )
       call Cuthill_Mckee(na_u,a_mm,R)
    case ( BW_CUTHILL_MCKEE_PRIORITY )
       call Cuthill_Mckee(na_u,a_mm,R, priority = priority)
    case ( BW_PAPIOR )
       call Papior(na_u,a_mm,R, priority = priority)
    case default
       call die('Could not recognize the bandwidth reduction algorithm')
    end select

    if ( method >= BW_REVERSE ) then

       ! Revert to correct atomic indices
       do ia = 1 , na_u
          R(ia) = na_u + 1 - R(ia)
       end do

       ! transfer to the reversed Cuthill-Mckee algorithm
       ! We could do this and the above in one go, but that would require a 
       ! correction if mod(io,2) == 1
       do ia = 1 , na_u / 2
          iac = R(ia)
          R(ia) = R(na_u+1-ia)
          R(na_u+1-ia) = iac
       end do

    end if
    
  end subroutine BandWidth_pivoting

  subroutine Cuthill_Mckee(na_u,a_mm,R,priority)
    ! We return the pivoting indexes for the atoms
    use parallel, only : IONode
    use alloc
    ! number of atoms in the cell, also the last orbitals of each
    ! atom.
    integer, intent(in) :: na_u
    integer, intent(in) :: a_mm(na_u,na_u)
    integer, intent(inout) :: R(na_u)
    integer, intent(in), optional :: priority(na_u)

    integer, pointer :: Q(:)
    integer :: ia, nQ, iQ, i_ca, ca
    integer :: degree

    ! 1. initialize queue array
    nullify(Q)
    call re_alloc(Q,1,na_u)

    ! we can now perform the Cuthill-Mckee algorithm
    ! Notice that we force the last atom in the electrode
    ! to be in the result array (thus we force
    ! the construction of the Cuthill-Mckee to start
    ! with the electrode)

    ! Prepare counters for the algorithm
    i_ca = 0

    ! the counter for the number of saved atoms
    ca = 0

    ! start the algorithm 
    do while ( ca < na_u )

       ! initialize the current node
       ! remark that the first iteration fixes
       ! the node building from the last atom
       ! specified
       nQ = 0
       if ( i_ca > 0 ) then
          ! In case we want a pure Cuthill-Mckee
          ! we shall not construct the queue the first
          ! time
          call add_Queue(na_u,a_mm,i_ca,Q,nQ, &
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
             R(ca) = Q(iQ)
             call add_Queue(na_u,a_mm,Q(iQ),Q,nQ, &
                  priority = priority)
          else if ( any( R(1:ca) == Q(iQ) ) ) then
             ! Just delete the entry (we do not need to worry about
             ! it anymore)
             Q(iQ:nQ-1) = Q(iQ+1:nQ)
             iQ = iQ - 1
             nQ = nQ - 1
          else
             ! We are allowed to add the queued item to the result list
             ca = ca + 1
             R(ca) = Q(iQ)
             ! update the queue-list
             call add_Queue(na_u,a_mm,Q(iQ),Q,nQ, &
                  priority=priority)
          end if
       end do

       ! Do a quick exit if possible
       if ( ca == na_u ) cycle

       ! From here on it is the actual Cuthill-Mckee algorithm
       ! Now we need to select the one with the lowest degree

       ! 2. select lowest degree
       i_ca = 0
       degree = huge(1)
       do ia = 1 , na_u

          ! check that it has not been added to the queue
          iQ = 0
          if ( ca == 0 ) then
             ! If no items has been added
             iQ = 1
          else if ( .not. any( R(1:ca) == ia ) ) then
             ! Check it hasn't been processed 
             iQ = 1
          end if

          if ( iQ == 1 ) then
             ! calculate the degree of the index
             iQ = sum(a_mm(:,ia))
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
    
    subroutine add_Queue(na_u,a_mm,ia,Q,nQ, priority)
      integer, intent(in) :: na_u, a_mm(na_u,na_u)
      integer, intent(in) :: ia
      integer, intent(inout) :: Q(na_u)
      integer, intent(inout) :: nQ
      integer, intent(in), optional :: priority(na_u)

      integer :: i, iQ, ii, lnQ, j
      integer :: degree(nQ+1:na_u)
      logical :: append, switch

      ! 3. extend the queue for the current parent (ia)
      ! We prepend to the queue
      iQ = nQ
      do i = 1 , na_u
         ! check that it is not the same atom
         if ( i /= ia ) cycle
         ! check that it is connected to the parent
         if ( a_mm(i,ia) /= 1 ) cycle

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
            degree(iQ) = sum(a_mm(:,i))
         end if
      end do

      ! capture the total connected nodes
      lnQ = iQ

      ! Sort the queued items (original)
      do i = nQ + 1 , lnQ - 1
         ! we sort by the first find the one to switch with
         do ii = lnQ , i + 1 , -1
            ! we switch to have the lowest degree first
            switch = degree(ii) < degree(i)

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

  subroutine Papior(na_u,a_mm,R,priority)
    ! We return the pivoting indexes for the atoms
    use parallel, only : IONode
    use alloc
    ! number of atoms in the cell, also the last orbitals of each
    ! atom.
    integer, intent(in) :: na_u
    integer, intent(in) :: a_mm(na_u,na_u)
    integer, intent(inout) :: R(na_u)
    integer, intent(in) :: priority(na_u)

    integer, pointer :: Q(:), pri_idx(:)
    integer :: ia, nQ, iQ, i_ca, ca, ja
    integer :: degree

    nullify(Q,pri_idx)
    call re_alloc(Q,1,na_u)
    call re_alloc(pri_idx,1,na_u)

    ! Notice that we force the last atom in the electrode
    ! to be in the result array (thus we force
    ! the construction of the Papior to start
    ! with the electrode)

    ! Prepare counters for the algorithm
    i_ca = 0

    ! the counter for the number of saved atoms
    ca = 0

    ! start the algorithm 
    do while ( ca < na_u )

       ! initialize the current node
       ! remark that the first iteration fixes
       ! the node building from the last atom
       ! specified
       nQ = 0
       if ( i_ca > 0 ) then
          ! In case we want a pure Papior
          ! we shall not construct the queue the first
          ! time
          call add_Queue(na_u,a_mm,i_ca,Q,nQ, &
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
             R(ca) = Q(iQ)
             call add_Queue(na_u,a_mm,Q(iQ),Q,nQ, &
                  priority=priority)
          else if ( any( R(1:ca) == Q(iQ) ) ) then
             ! Just delete the entry (we do not need to worry about
             ! it anymore)
             Q(iQ:nQ-1) = Q(iQ+1:nQ)
             iQ = iQ - 1
             nQ = nQ - 1
          else
             ! We are allowed to add the queued item to the result list
             ca = ca + 1
             R(ca) = Q(iQ)
             ! update the queue-list
             call add_Queue(na_u,a_mm,Q(iQ),Q,nQ, &
                  priority=priority)
          end if
       end do

       ! Do a quick exit if possible
       if ( ca == na_u ) cycle

       ! From here on it is the actual Papior algorithm
       ! Now we need to select the one with the lowest degree

       i_ca = 0
       degree = huge(1)
       do ia = 1 , na_u
          iQ = 0

          if ( ca == 0 ) then
             ! If no items has been added
             iQ = 1
          else if ( .not. any( R(1:ca) == ia ) ) then
             ! Check it hasn't been processed 
             iQ = 1
          end if

          if ( iQ == 1 ) then
             do ja = ia , na_u ! from i ensures that we get at least a zero
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
    subroutine add_Queue(na_u,a_mm,ia,Q,nQ, priority)
      integer, intent(in) :: na_u, a_mm(na_u,na_u)
      integer, intent(in) :: ia
      integer, intent(inout) :: Q(na_u)
      integer, intent(inout) :: nQ
      integer, intent(in) :: priority(na_u)

      integer :: i, iQ, ii, lnQ, j
      integer :: degree(nQ+1:na_u)
      logical :: append, switch

      ! We prepend to the queue
      iQ = nQ
      do i = 1 , na_u
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
               do j = i , na_u ! from i ensures that we get at least a zero
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
