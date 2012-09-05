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
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!
module m_hs_matrix
!
! Routines that are used for the distribution of the Hamiltonian 
! and scattering matrices into a full size matrix instead of sparse matrices.
! It has several options for creating different kinds of matrices.
! This module is a serial version. It requires that the Node has the full
! sparse matrix available.
!
! Specifically it can be used to remove the z-direction connections.
! Also the inner-cell distances are in the xij arrays. Through routine calls
! these inner cell distances can be removed.
!
! The reason for having this is future use of routines which can utilize the
! Hamiltonian created in different ways.
! Say you want to test calculate a transmission from a Hamiltonian which is
! created by SIESTA. For that purpose you need to remove the cell connection 
! in the z-direction.
!
! The usage of this module is highly encouraged in future utilities where
! the need for the Hamiltonian and/or overlap matrix is needed.
! With this in mind the code can further be optimized for speed.
! However, for normal sizes it is quite fast.
! 
! Also for testing against Inelastica the use of inner-cell distances is not
! used. Therefore the Hamiltonians can not be numerically compared.
! If this is to be enforced later on. The option is there.
! 
! The use of this module is a straight forward call:
! 
!   call set_HS_matrix(Gamma,ucell,na_u,no_u,no_s,maxnh, &
!       xij,numh,listhptr,listh,indxuo,H,S, &
!       k,Hk,Sk, &
!       xa,iaorb, &
!       RemZConnection,RemUCellDistances,RemNFirstOrbitals,RemNLastOrbitals)
! 
! xa, iaorb, RemZConnection, RemUCellDistances, RemNFirstOrbitals, RemNLastOrbitals
! are all optional arguments. xa and iaorb are required if RemZConnection or RemUCellDistances
! are true.
! Notice, that any calls to the optional arguments MUST be with keywords! Otherwise
! the program will end!
! TODO : CONSIDERATION : iaorb could be replaced with lasto, however, iaorb is easier to use (perhaps both could be optioned?)
!
! Gamma denotes whether it is a Gamma calculation (if true, it will not
! add k-phases, no matter if k /= \Gamma-point.
! na_u,no_u,no_s,maxnh,xij,numh,listhptr,listh,indxuo,H,S are all variables
! needed in the definition of the entire H and S matrices in the sparse format.
!
! k is the k-point that will be created for the Hamiltonian.
! Hk and Sk are returned for the user.
!
! The RemZConnection can be used to remove any matrix elements <i|H|j> where i and j
! are connections in the next unit cell in the Z-direction.
! This is usefull if one wishes to see the matrix as it would look while doing 
! transmission calculations.
!
! RemUCellDistances can be set to true so that inner-cell distances are removed.
! Several people on the SIESTA mailing list have "complained/asked questions" about 
! the difference in the Hamiltonians which are not always constructed similarly. 
! However, inner cell differences can be neglected. This is merely to get the same
! matrix as is created in Inelastica, for example.
!
! RemNFirstOrbitals can be used to fully remove that many states from the start of the
! Hamiltonian. This is used when there are buffer regions which should be disregarded.
! RemNLastOrbitals is the same, albeit in the end of the Hamiltonian.
!
! If the sizes of the incoming pointers Hk and Sk do not match the above
! they will be reallocated to the correct size.
!
! Furthermore there are the routines:
!   call matrix_rem_left_right(no_tot,Hk,Sk,no_L,no_R)
!   and 
!   call matrix_symmetrize(no_tot,Hk,Sk,Ef)
! Which are used to remove connections of left/right regions
! and used to symmetrize and shift the Hamiltonian Ef, respectively.
! 
! NOTICE that a call to matrix_symmetrize is almost always needed!
 
  implicit none

  private

  interface set_HS_matrix
     module procedure set_HS_matrix_1d
     module procedure set_HS_matrix_2d
  end interface set_HS_matrix

  interface matrix_rem_left_right
     module procedure matrix_rem_left_right_1d
     module procedure matrix_rem_left_right_2d
  end interface matrix_rem_left_right

  interface matrix_symmetrize
     module procedure matrix_symmetrize_1d
     module procedure matrix_symmetrize_2d
  end interface matrix_symmetrize

  public :: set_HS_matrix
  public :: matrix_rem_left_right
  public :: matrix_symmetrize

contains

!*****************
! Setting the Hamiltonian for a specific k-point.
! It requires that the Node has the full sparse matrix available
!*****************
  subroutine set_HS_matrix_1d(Gamma,ucell,na_u,no_u,no_s,maxnh, &
       xij,numh,listhptr,listh,indxuo,H,S, &
       k,Hk,Sk, &
       DUMMY, & ! Ensures that the programmer makes EXPLICIT keywork passing
       xa,iaorb, &
       RemZConnection,RemUCellDistances,RemNFirstOrbitals,RemNLastOrbitals)
    use precision, only : dp
    use sys,       only : die 
! ***********************
! * INPUT variables     *
! ***********************
    logical, intent(in)               :: Gamma ! Is it a Gamma Calculation?
    real(dp), intent(in)              :: ucell(3,3) ! The unit cell of system
    integer, intent(in)               :: na_u ! Unit cell atoms
    integer, intent(in)               :: no_u ! Unit cell orbitals
    integer, intent(in)               :: no_s ! Supercell orbitals
    integer, intent(in)               :: maxnh ! Hamiltonian size
    real(dp), intent(in)              :: xij(3,maxnh) ! differences with unitcell, differences with unitcell
    integer, intent(in)               :: numh(no_u),listhptr(no_u)
    integer, intent(in)               :: listh(maxnh),indxuo(no_s)
    real(dp), intent(in)              :: H(maxnh) ! Hamiltonian
    real(dp), intent(in)              :: S(maxnh) ! Overlap
    real(dp), intent(in)              :: k(3) ! k-point in [1/Bohr]
! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), pointer, intent(out) :: Hk(:),Sk(:)

! ***********************
! * OPTIONAL variables  *
! ***********************
    logical, intent(in), optional :: DUMMY ! Do not supply this, it merely requires the coder
!                                          ! to use the keyworded arguments!
    real(dp), intent(in),optional :: xa(3,na_u) ! Atomic coordinates (needed for RemZConnection & RemUCellDistances)
    integer, intent(in), optional :: iaorb(no_u) ! The equivalent atomic index for a given orbital (needed for RemUCellDistances)
    logical, intent(in), optional :: RemZConnection, RemUCellDistances
    integer, intent(in), optional :: RemNFirstOrbitals, RemNLastOrbitals

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: recell(3,3) ! The reciprocal unit cell
    real(dp) :: xo(3), xc
    real(dp), allocatable :: xuo(:)
    real(dp) :: kxij
    complex(dp) :: cphase
    integer :: no_tot
    integer :: i,j,iu,iuo,juo,iind,ind
    logical :: l_RemZConnection, l_RemUCellDistances
    integer :: l_RemNFirstOrbitals, l_RemNLastOrbitals 

    if ( present(DUMMY) ) &
         call die("You must specify the keyworded arguments &
         &for set_HS_matrix")

    ! Option collecting
    l_RemZConnection = .false.
    if ( present(RemZConnection) ) &
         l_RemZConnection = RemZConnection
    if (l_RemZConnection .and. .not. present(xa)) &
         call die("You need xa in set_HS_matrix when removing &
         &the z-connection.")
    if (l_RemZConnection .and. .not. present(iaorb)) &
         call die("You need iaorb in set_HS_matrix when removing &
         &the z-connection.")
    l_RemUCellDistances = .false.
    if ( present(RemUCellDistances) ) &
         l_RemUCellDistances = RemUCellDistances
    if (l_RemUCellDistances .and. .not. present(xa)) &
         call die("You need xa in set_HS_matrix when removing &
         &unit cell distances.")
    if (l_RemUCellDistances .and. .not. present(iaorb)) &
         call die("You need iaorb in set_HS_matrix when removing &
         &unit cell distances.")

    ! Make l_RemNFirstOrbitals contain the number of orbitals
    ! to be removed from the Hamiltonian in the BEGINNING
    l_RemNFirstOrbitals = 0
    if ( present(RemNFirstOrbitals) ) &
         l_RemNFirstOrbitals = RemNFirstOrbitals
    ! Make l_RemNLastOrbitals contain the number of orbitals
    ! to be removed from the Hamiltonian in the END
    l_RemNLastOrbitals = 0
    if ( present(RemNLastOrbitals) ) &
         l_RemNLastOrbitals = RemNLastOrbitals
    no_tot = no_u - (l_RemNLastOrbitals + l_RemNFirstOrbitals)


    if ( associated(Hk) .and. associated(Sk) ) then
       if ( size(Hk) /= no_tot*no_tot ) then
          call memory('D','Z',size(Hk)+size(Sk),'set_HS')
          deallocate(Hk,Sk)
          nullify(Hk,Sk)
          allocate(Hk(no_tot*no_tot),Sk(no_tot*no_tot))
          call memory('A','Z',2*no_tot*no_tot,'set_HS')
       end if
    else
       ! No need to nullify...
       allocate(Hk(no_tot*no_tot),Sk(no_tot*no_tot))
       call memory('A','Z',2*no_tot*no_tot,'set_HS')
    end if
    
    if ( l_RemZConnection ) then
       ! Prepare the cell to calculate the index of the atom
       call reclat(ucell,recell,0) ! Without 2*Pi
       
       ! Find the actual coordinates of the orbitals in the form of the sparse matrices
       ! Notice that this array is without the removed orbitals
       allocate(xuo(no_tot))
       call memory('A','D',no_tot,'consHS')

       do iuo = 1 , no_tot
          i = iaorb(iuo + l_RemNFirstOrbitals)
          xuo(iuo) = &
               xa(1,i) * recell(1,3) + &
               xa(2,i) * recell(2,3) + &
               xa(3,i) * recell(3,3)
       end do !io in uc
    end if


!
! Setup H,S for this k-point:
!
    do i = 1,no_tot*no_tot
       Hk(i) = dcmplx(0.d0,0.d0)
       Sk(i) = dcmplx(0.d0,0.d0)
    end do

    xo(:) = 0.0_dp

    setup_HS: if (.not.Gamma ) then

       do iuo = 1 , no_tot
          iu = iuo + l_RemNFirstOrbitals
          iind = listhptr(iu)
          do j = 1,numh(iu)
             ind = iind + j
             juo = indxuo(listh(ind)) - l_RemNFirstOrbitals

             ! Cycle if we are not in the middle region
             if ( juo < 1 .or. no_tot < juo ) cycle

             ! If we have a removal of the z-direction do the following
             if ( l_RemZConnection ) then
                xc = - (xuo(juo) - xuo(iuo))
                xc = xc + xij(1,ind) * recell(1,3) + &
                          xij(2,ind) * recell(2,3) + &
                          xij(3,ind) * recell(3,3)
                if ( nint(xc) /= 0 ) cycle
             end if

             ! We also wish to remove the connection in
             ! in the inner cell
             if ( l_RemUCellDistances ) then
                xo(1) = xa(1,iaorb(juo + l_RemNFirstOrbitals)) &
                     -  xa(1,iaorb(iu))
                xo(2) = xa(2,iaorb(juo + l_RemNFirstOrbitals)) &
                     -  xa(2,iaorb(iu))
                xo(3) = xa(3,iaorb(juo + l_RemNFirstOrbitals)) &
                     -  xa(3,iaorb(iu))
             end if

             kxij = &
                  k(1) * (xij(1,ind) - xo(1)) + &
                  k(2) * (xij(2,ind) - xo(2)) + &
                  k(3) * (xij(3,ind) - xo(3))
             cphase = exp(dcmplx(0d0,1d0)*kxij)
             i = iuo+(juo-1)*no_tot
             Hk(i) = Hk(i)+H(ind)*cphase
             Sk(i) = Sk(i)+S(ind)*cphase
          end do
       end do

    else setup_HS
       ! It is not a Gamma calculation, thus we do not have any
       ! neighbouring cells etc.
       do iuo = 1 , no_tot
          iu = iuo + l_RemNFirstOrbitals
          iind = listhptr(iu)
          do j = 1,numh(iu)
             ind = iind + j
             juo = listh(ind) - l_RemNFirstOrbitals

             ! Cycle if we are not in the middle region
             if ( juo < 1 .or. no_tot < juo ) cycle

             ! If we have a removal of the z-direction do the following
             if ( l_RemZConnection ) then
                xc = - (xuo(juo) - xuo(iuo))
                xc = xc + xij(1,ind) * recell(1,3) + &
                          xij(2,ind) * recell(2,3) + &
                          xij(3,ind) * recell(3,3)
                if ( nint(xc) /= 0 ) cycle
             end if

             i = iuo+(juo-1)*no_tot
             Hk(i) = Hk(i)+H(ind)
             Sk(i) = Sk(i)+S(ind)
          end do
       end do

    end if setup_HS

    if ( l_RemZConnection ) then
       call memory('D','D',no_tot,'consHS')
       deallocate(xuo)
    end if

  end subroutine set_HS_matrix_1d

!*****************
! Setting the Hamiltonian for a specific k-point.
! It requires that the Node has the full sparse matrix available
!*****************
  subroutine set_HS_matrix_2d(Gamma,ucell,na_u,no_u,no_s,maxnh, &
       xij,numh,listhptr,listh,indxuo,H,S, &
       k,Hk,Sk, &
       DUMMY, & ! Ensures that the programmer makes EXPLICIT keywork passing
       xa,iaorb, &
       RemZConnection,RemUCellDistances,RemNFirstOrbitals,RemNLastOrbitals)
    use precision, only : dp
    use sys,       only : die 
! ***********************
! * INPUT variables     *
! ***********************
    logical, intent(in)               :: Gamma ! Is it a Gamma Calculation?
    real(dp), intent(in)              :: ucell(3,3) ! The unit cell of system
    integer, intent(in)               :: na_u ! Unit cell atoms
    integer, intent(in)               :: no_u ! Unit cell orbitals
    integer, intent(in)               :: no_s ! Total orbitals
    integer, intent(in)               :: maxnh ! Hamiltonian size
    real(dp), intent(in)              :: xij(3,maxnh) ! differences with unitcell, differences with unitcell
    integer, intent(in)               :: numh(no_u),listhptr(no_u)
    integer, intent(in)               :: listh(maxnh),indxuo(no_s)
    real(dp), intent(in)              :: H(maxnh) ! Hamiltonian
    real(dp), intent(in)              :: S(maxnh) ! Overlap
    real(dp), intent(in)              :: k(3) ! k-point in [1/Bohr]
! ***********************
! * OUTPUT variables    *
! ***********************
    complex(dp), pointer, intent(out) :: Hk(:,:),Sk(:,:)
! ***********************
! * OPTIONAL variables  *
! ***********************
    logical, intent(in), optional :: DUMMY ! Do not supply this, it merely requires the coder
!                                          ! to use the keyworded arguments!
    real(dp), intent(in),optional :: xa(3,na_u) ! Atomic coordinates (needed for RemZConnection & RemUCellDistances)
    integer, intent(in), optional :: iaorb(no_u) ! The equivalent atomic index for a given orbital (needed for RemUCellDistances)
    logical, intent(in), optional :: RemZConnection, RemUCellDistances
    integer, intent(in), optional :: RemNFirstOrbitals, RemNLastOrbitals

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: recell(3,3)
    real(dp) :: xo(3), xc
    real(dp), allocatable :: xuo(:)
    integer :: no_tot
    real(dp) :: kxij
    complex(dp) :: cphase
    integer :: i,j,iuo,iu,juo,iind,ind
    logical :: l_RemZConnection, l_RemUCellDistances
    integer :: l_RemNFirstOrbitals, l_RemNLastOrbitals 

    if ( present(DUMMY) ) &
         call die("You must specify the keyworded arguments &
         &for set_HS_matrix")

    ! Option collecting
    l_RemZConnection = .false.
    if ( present(RemZConnection) ) &
         l_RemZConnection = RemZConnection
    if (l_RemZConnection .and. .not. present(xa)) &
         call die("You need xa in set_HS_matrix when removing &
         &the z-connection.")
    if (l_RemZConnection .and. .not. present(iaorb)) &
         call die("You need iaorb in set_HS_matrix when removing &
         &the z-connection.")
    l_RemUCellDistances = .false.
    if ( present(RemUCellDistances) ) &
         l_RemUCellDistances = RemUCellDistances
    if (l_RemUCellDistances .and. .not. present(xa)) &
         call die("You need xa in set_HS_matrix when removing &
         &unit cell distances.")
    if (l_RemUCellDistances .and. .not. present(iaorb)) &
         call die("You need iaorb in set_HS_matrix when removing &
         &unit cell distances.")


    ! Make l_RemNFirstOrbitals contain the number of orbitals
    ! to be removed from the Hamiltonian in the BEGINNING
    l_RemNFirstOrbitals = 0
    if ( present(RemNFirstOrbitals) ) &
         l_RemNFirstOrbitals = RemNFirstOrbitals
    ! Make l_RemNLastOrbitals contain the number of orbitals
    ! to be removed from the Hamiltonian in the END
    l_RemNLastOrbitals = 0
    if ( present(RemNLastOrbitals) ) &
         l_RemNLastOrbitals = RemNLastOrbitals
    no_tot = no_u - (l_RemNLastOrbitals + l_RemNFirstOrbitals)


    if ( associated(Hk) .and. associated(Sk) ) then
       if ( size(Hk) /= no_tot*no_tot ) then
          call memory('D','Z',size(Hk)+size(Sk),'set_HS')
          deallocate(Hk,Sk)
          nullify(Hk,Sk)
          allocate(Hk(no_tot,no_tot),Sk(no_tot,no_tot))
          call memory('A','Z',2*no_tot*no_tot,'set_HS')
       end if
    else
       ! No need to nullify
       allocate(Hk(no_tot,no_tot),Sk(no_tot,no_tot))
       call memory('A','Z',2*no_tot*no_tot,'set_HS')
    end if

    if ( l_RemZConnection ) then
       ! Prepare the cell to calculate the index of the atom
       call reclat(ucell,recell,0) ! Without 2*Pi
       
       ! Find the actual coordinates of the orbitals in the form of the sparse matrices
       ! Notice that this array is without the removed orbitals
       allocate(xuo(no_tot))
       call memory('A','D',no_tot,'consHS')

       do iuo = 1 , no_tot
          i = iaorb(iuo + l_RemNFirstOrbitals)
          xuo(iuo) = &
               xa(1,i) * recell(1,3) + &
               xa(2,i) * recell(2,3) + &
               xa(3,i) * recell(3,3)
       end do !io in uc
    end if


!
! Setup H,S for this k-point:
!
    do juo = 1,no_tot
       do iuo = 1,no_tot
          Hk(iuo,juo) = dcmplx(0.d0,0.d0)
          Sk(iuo,juo) = dcmplx(0.d0,0.d0)
       end do
    end do

    xo(:) = 0.0_dp

    setup_HS: if (.not.Gamma ) then

       do iuo = 1,no_tot
          iu = iuo + l_RemNFirstOrbitals
          iind = listhptr(iu)
          do j = 1,numh(iu)
             ind = iind + j
             juo = indxuo(listh(ind)) - l_RemNFirstOrbitals

             ! Cycle if we are not in the middle region
             if ( juo < 1 .or. no_tot < juo ) cycle

             ! If we have a removal of the z-direction do the following
             if ( l_RemZConnection ) then
                xc = - (xuo(juo) - xuo(iuo))
                xc = xc + xij(1,ind) * recell(1,3) + &
                          xij(2,ind) * recell(2,3) + &
                          xij(3,ind) * recell(3,3)
                if ( nint(xc) /= 0 ) cycle
             end if

             ! We also wish to remove the connection in
             ! in the inner cell
             if ( l_RemUCellDistances ) then
                xo(1) = xa(1,iaorb(juo + l_RemNFirstOrbitals)) &
                     -  xa(1,iaorb(iu))
                xo(2) = xa(2,iaorb(juo + l_RemNFirstOrbitals)) &
                     -  xa(2,iaorb(iu))
                xo(3) = xa(3,iaorb(juo + l_RemNFirstOrbitals)) &
                     -  xa(3,iaorb(iu))
             end if

             kxij = &
                  k(1) * (xij(1,ind) - xo(1)) + &
                  k(2) * (xij(2,ind) - xo(2)) + &
                  k(3) * (xij(3,ind) - xo(3))
             cphase = exp(dcmplx(0d0,1d0)*kxij)
             Hk(iuo,juo) = Hk(iuo,juo)+H(ind)*cphase
             Sk(iuo,juo) = Sk(iuo,juo)+S(ind)*cphase
          end do
       end do

    else setup_HS

       do iuo = 1,no_tot
          iu = iuo + l_RemNFirstOrbitals
          iind = listhptr(iu)
          do j = 1,numh(iu)
             ind = iind + j
             juo = listh(ind) - l_RemNFirstOrbitals

             ! Cycle if we are not in the middle region
             if ( juo < 1 .or. no_tot < juo ) cycle

             ! If we have a removal of the z-direction do the following
             if ( l_RemZConnection ) then
                xc = - (xuo(juo) - xuo(iuo))
                xc = xc + xij(1,ind) * recell(1,3) + &
                          xij(2,ind) * recell(2,3) + &
                          xij(3,ind) * recell(3,3)
                if ( nint(xc) /= 0 ) cycle
             end if

             Hk(iuo,juo) = Hk(iuo,juo)+H(ind)
             Sk(iuo,juo) = Sk(iuo,juo)+S(ind)
          end do
       end do

    end if setup_HS

    if ( l_RemZConnection ) then
       call memory('D','D',no_tot,'consHS')
       deallocate(xuo)
    end if

  end subroutine set_HS_matrix_2d


  ! Routine for removing left right overlaps of certain regions.
  ! Is used to fully remove the connection between left and right
  ! states in the Hamiltonian
  subroutine matrix_rem_left_right_1d(no_tot,Hk,Sk,no_L,no_R)
    use precision, only : dp

! **************************
! * INPUT variables        *
! **************************
    integer, intent(in) :: no_tot, no_L, no_R

! **************************
! * OUTPUT variables       *
! **************************
    complex(dp), intent(inout) :: Hk(no_tot*no_tot),Sk(no_tot*no_tot)

! **************************
! * LOCAL variables        *
! **************************
    integer :: i,j

    ! If nothing is to be removed return immidiately...
    if ( no_L == 0 .or. no_R == 0 ) return

    do j = no_tot - no_R + 1 , no_tot
       do i = 1 , no_L
          Hk(i+(j-1)*no_tot) = dcmplx(0.d0,0.d0)
          Sk(i+(j-1)*no_tot) = dcmplx(0.d0,0.d0)
          Hk(j+(i-1)*no_tot) = dcmplx(0.d0,0.d0)
          Sk(j+(i-1)*no_tot) = dcmplx(0.d0,0.d0)
       end do
    end do

  end subroutine matrix_rem_left_right_1d

  ! Routine for removing left right overlaps of certain regions.
  ! Is used to fully remove the connection between left and right
  ! states in the Hamiltonian
  subroutine matrix_rem_left_right_2d(no_tot,Hk,Sk,no_L,no_R)
    use precision, only : dp

! **************************
! * INPUT variables        *
! **************************
    integer, intent(in) :: no_tot, no_L, no_R

! **************************
! * OUTPUT variables       *
! **************************
    complex(dp), intent(inout) :: Hk(no_tot,no_tot),Sk(no_tot,no_tot)

! **************************
! * LOCAL variables        *
! **************************
    integer :: i,j

    ! If nothing is to be removed return immidiately...
    if ( no_L == 0 .or. no_R == 0 ) return

    do j = no_tot - no_R + 1 , no_tot
       do i = 1 , no_L
          Hk(i,j) = dcmplx(0.d0,0.d0)
          Sk(i,j) = dcmplx(0.d0,0.d0)
          Hk(j,i) = dcmplx(0.d0,0.d0)
          Sk(j,i) = dcmplx(0.d0,0.d0)
       end do
    end do

  end subroutine matrix_rem_left_right_2d



  subroutine matrix_symmetrize_1d(no_tot,Hk,Sk,Ef)
    use precision, only : dp

! **************************
! * INPUT variables        *
! **************************
    integer, intent(in) :: no_tot
    real(dp), intent(in) :: Ef

! **************************
! * OUTPUT variables       *
! **************************
    complex(dp), intent(inout) :: Hk(no_tot*no_tot),Sk(no_tot*no_tot)
    
! **************************
! * LOCAL variables        *
! **************************
    integer :: i,j,iuo,juo

    do i = 1,no_tot
       do j = 1,i-1
          iuo = i + no_tot*(j-1)
          juo = j + no_tot*(i-1)

          Sk(juo) = 0.5d0*( Sk(juo) + dconjg(Sk(iuo)) )
          Sk(iuo) =  dconjg(Sk(juo))
          
          Hk(juo) = 0.5d0*( Hk(juo) + dconjg(Hk(iuo)) ) &
               - Ef*Sk(juo)
          Hk(iuo) =  dconjg(Hk(juo))
          
       end do
       iuo = i + no_tot*(i-1)
       Sk(iuo)=Sk(iuo) - dcmplx(0d0,1d0)*dimag(Sk(iuo))
       
       Hk(iuo)=Hk(iuo) - dcmplx(0d0,1d0)*dimag(Hk(iuo)) &
            - Ef*Sk(iuo) 
    end do

  end subroutine matrix_symmetrize_1d



  ! Routine for symmetrizing and shifting the matrix Ef
  subroutine matrix_symmetrize_2d(no_tot,Hk,Sk,Ef)
    use precision, only : dp

! **************************
! * INPUT variables        *
! **************************
    integer, intent(in) :: no_tot
    real(dp), intent(in) :: Ef

! **************************
! * OUTPUT variables       *
! **************************
    complex(dp), intent(inout) :: Hk(no_tot,no_tot),Sk(no_tot,no_tot)
    
! **************************
! * LOCAL variables        *
! **************************
    integer :: iuo,juo

    do iuo = 1,no_tot
       do juo = 1,iuo-1
          
          Sk(juo,iuo) = 0.5d0*( Sk(juo,iuo) + dconjg(Sk(iuo,juo)) )
          Sk(iuo,juo) =  dconjg(Sk(juo,iuo))
          
          Hk(juo,iuo) = 0.5d0*( Hk(juo,iuo) + dconjg(Hk(iuo,juo)) ) &
               - Ef*Sk(juo,iuo)
          Hk(iuo,juo) =  dconjg(Hk(juo,iuo))
          
       end do
       
       Sk(iuo,iuo)=Sk(iuo,iuo) - dcmplx(0d0,1d0)*dimag(Sk(iuo,iuo))
       
       Hk(iuo,iuo)=Hk(iuo,iuo) - dcmplx(0d0,1d0)*dimag(Hk(iuo,iuo)) &
            - Ef*Sk(iuo,iuo) 
    end do

  end subroutine matrix_symmetrize_2d

end module m_hs_matrix
  
