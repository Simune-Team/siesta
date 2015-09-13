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
! Nick Papior Andersen, 2014, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! This particular solution method relies on solving the GF
! with the tri-diagonalization routine.
! This will leverage memory usage and also the execution time.

module m_tbt_tri_init

  use precision, only : dp, i8b
  use m_region

  implicit none

  public :: tbt_tri_init
  public :: tbt_tri_print_opti

  type(tRgn), allocatable, target :: ElTri(:)
  type(tRgn) :: DevTri

  public :: ElTri, DevTri
  public :: fold_elements

  private
  
contains

  subroutine tbt_tri_init_elec( dit , sp )

    use parallel, only : Node, Nodes
    use class_OrbitalDistribution
    use class_Sparsity
    use create_Sparsity_Union

    use m_ts_rgn2trimat
    use m_ts_electype
#ifdef MPI
    use mpi_siesta
#endif
#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

    use m_sparsity_handling
    use m_tbt_options, only : N_Elec, Elecs, BTD_method
    use m_tbt_regions

    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp

    type(Sparsity) :: tmpSp1, tmpSp2
    integer :: i, iEl

#ifdef MPI
    integer :: MPIerror
#endif

    call timer('tri-init-elec',1)

    ! This works as creating a new sparsity deletes the previous
    ! and as it is referenced several times it will not be actually
    ! deleted...
    allocate(ElTri(N_Elec))
    do i = 1 + Node , N_Elec , Nodes

       ! Retain region
       call Sp_retain_rgn(dit,sp,r_oElpD(i),tmpSp2)

       ! Add the self-energy of the electrode (in its original position)
       call crtSparsity_Union_region(dit,tmpSp2, r_oEl_alone(i),tmpSp1)
       call delete(tmpSp2)

#ifdef TRANSIESTA_DEBUG
       open(file='ELEC_'//trim(Elecs(i)%name)//'_SP',unit=1400,form='formatted')
       call sp_to_file(1400,tmpSp1)
       close(1400)
#endif

       ! Create tri-diagonal parts for this electrode
       ! IF parts == 0 will create new partition
       call ts_rgn2TriMat(N_Elec, Elecs, .false., &
            dit, tmpSp1, r_oElpD(i), ElTri(i)%n, ElTri(i)%r, &
            BTD_method, last_eq = Elecs(i)%o_inD%n , par = .false. )
       call delete(tmpSp1)

    end do

    ! The i'th processor has the following electrodes
    do iEl = 1 , N_Elec

       ! The node having this electrode is
       i = iEl - 1
       do while ( i >= Nodes )
          i = i - Nodes
       end do
          
       ! Set the name 
       ElTri(iEl)%name = '[TRI] '//trim(Elecs(iEl)%name)
       
#ifdef MPI
       call MPI_Bcast(ElTri(iEl)%n,1,MPI_Integer, &
            i,MPI_Comm_World,MPIerror)
       if ( Node /= i ) then
          allocate(ElTri(iEl)%r(ElTri(iEl)%n))
       end if
       call MPI_Bcast(ElTri(iEl)%r,ElTri(iEl)%n,MPI_Integer, &
            i,MPI_Comm_World,MPIerror)
#endif
       
    end do

    call timer('tri-init-elec',2)
    
  end subroutine tbt_tri_init_elec

  subroutine tbt_tri_init( dit , sp , proj )

    use parallel, only : IONode
    use class_OrbitalDistribution
    use class_Sparsity
    use create_Sparsity_Union

    use m_ts_rgn2trimat
    use m_ts_tri_common, only : nnzs_tri_i8b
    use m_ts_electype
#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

    use m_sparsity_handling
    use m_tbt_options, only : N_Elec, Elecs, BTD_method
    use m_tbt_regions

    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp
    ! An array of additional projection regions
    ! which determines the projection of a molecule
    ! onto seperate regions
    type(tRgn), intent(in), optional :: proj(:)

    type(Sparsity) :: tmpSp1, tmpSp2
    integer :: i
    integer(i8b) :: els

    call timer('tri-init',1)

    if ( IONode ) then
       write(*,'(/,a)')'tbtrans: Creating electrode tri-diagonal matrix blocks'
    end if

    ! Copy over sparsity pattern
    tmpSp1 = sp
    
    do i = 1 , N_Elec

       ! Add the self-energy of the electrode in the projected position
       ! of the "device" region.
       if ( mod(i,2) == 1 ) then
          call crtSparsity_Union_region(dit,tmpSp1,Elecs(i)%o_inD,tmpSp2)
       else
          call crtSparsity_Union_region(dit,tmpSp2,Elecs(i)%o_inD,tmpSp1)
       end if

    end do
    if ( mod(N_Elec,2) == 1 ) then
       tmpSp1 = tmpSp2
    end if
    call delete(tmpSp2)

    ! We have now already added the projected position of the
    ! self-energies. Create the tri-diagonal matrices
    ! for the electrode down-folding regions
    call tbt_tri_init_elec( dit , tmpSp1 )

    if ( IONode ) then
       write(*,'(a)')'tbtrans: Creating device tri-diagonal matrix blocks'
    end if

    if ( present(proj) ) then
       do i = 1 , size(proj)

          ! Add the self-energy of the electrode in the projected position
          ! of the "device" region.
          call crtSparsity_Union_region(dit,tmpSp1,proj(i),tmpSp2)
          tmpSp1 = tmpSp2
          
       end do
       call delete(tmpSp2)
    end if

    ! Create the device region sparsity pattern by removing everything
    ! else....
    call Sp_retain_rgn(dit,tmpSp1,r_oDev,tmpSp2)
    call delete(tmpSp1)

#ifdef TRANSIESTA_DEBUG
    open(file='DEV_FULL_SP',unit=1400,form='formatted')
    call sp_to_file(1400,tmpSp2)
    close(1400)
#endif

    call rgn_delete(DevTri)

    ! Create tri-diagonal parts for this one...
    call ts_rgn2TriMat(N_Elec, Elecs, .true., &
       dit, tmpSp2, r_oDev, DevTri%n, DevTri%r, &
       BTD_method, last_eq = 0, par = .true. )
    call delete(tmpSp2) ! clean up

    DevTri%name = '[TRI] device region'

    if ( IONode ) then
       ! Print out stuff
       call rgn_print(DevTri, seq_max = 10 , collapse = .false.)
       ! Print out memory estimate
       els = nnzs_tri_i8b(DevTri%n,DevTri%r)
       write(*,'(a,i0)') 'tbtrans: Matrix elements in BTD: ', els

       write(*,'(/,a)') 'tbtrans: Electrodes tri-diagonal matrices'
       do i = 1 , N_Elec
          call rgn_print(ElTri(i), seq_max = 10 , collapse = .false.)
       end do
    end if

    call timer('tri-init',2)

  end subroutine tbt_tri_init

  subroutine tbt_tri_print_opti(na_u,lasto,r_oDev,N_Elec)
    use parallel, only : IONode
    use m_region
    use m_verbosity, only : verbosity
    integer, intent(in) :: na_u, lasto(0:na_u)
    type(tRgn), intent(in) :: r_oDev
    integer, intent(in) :: N_Elec

    integer :: cum_sum, li, i, off
    type(tRgn) :: ro, ra

    if ( .not. IONode ) return
    if ( verbosity < 2 ) return
    if ( N_Elec > 2 ) return

    ! In case we have more than two tri-mat regions we can advice
    ! the user to a minimal tri-mat matrix

    ! In fact we should do this analysis for a sparsity pattern
    ! without any self-energies added, this would correctly
    ! get the size of the couplings and can shrink the 
    ! transmission region based on the mininum connections
    ! the self-energies will follow.

    if ( DevTri%n > 2 ) then
       
       ! We take the minimal region in the middle
       li = 0
       cum_sum = huge(1)
       do i = 2 , DevTri%n - 1
          if ( DevTri%r(i) < cum_sum ) then
             li = i
             cum_sum = DevTri%r(i)
          end if
       end do 

       if ( DevTri%r(1) + DevTri%r(2) < cum_sum ) then
          li = 0
          cum_sum = DevTri%r(1) + DevTri%r(2)
          call rgn_list(ro,cum_sum,r_oDev%r)
       end if

       if ( DevTri%r(DevTri%n-1) + DevTri%r(DevTri%n) < cum_sum ) then
          li = 0
          cum_sum = DevTri%r(DevTri%n-1) + DevTri%r(DevTri%n)
          call rgn_list(ro,cum_sum,r_oDev%r(r_oDev%n-cum_sum+1:))
       end if

       if ( li > 0 ) then
          off = 1
          do i = 1 , li - 1
             off = off + DevTri%r(i)
          end do
          call rgn_list(ro,cum_sum,r_oDev%r(off:))
       end if

       call rgn_Orb2Atom(ro,na_u,lasto,ra)
       call rgn_delete(ro)

       ! Sort transmission region atoms
       call rgn_sort(ra)

       write(*,*) ''
       write(*,'(a)') 'tbtrans: Suggested atoms for fastest transmission calculation:'

       ra%name = '[A]-Fast transmission'
       call rgn_print(ra, seq_max = 12)
       write(*,*) ''

       ! Clean-up
       call rgn_delete(ra)

    end if

  end subroutine tbt_tri_print_opti

  function fold_elements(N_tri,tri) result(elem)
    integer, intent(in) :: N_tri, tri(N_tri)
    integer :: elem, i, tmp

    elem = 0
    tmp = 0
    do i = 1 , N_tri - 1
       tmp = tri(i)**2
       tmp = tmp + tri(i)*( tri(i) + 2 * tri(i+1) )
       tmp = tmp + tri(i+1) ** 2
       elem = max(elem,tmp)
    end do

  end function fold_elements

end module m_tbt_tri_init
