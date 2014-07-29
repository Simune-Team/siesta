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

module m_ts_tri_init

  use precision, only : dp

  implicit none

  ! arrays for containing the tri-diagonal matrix part sizes
  integer, pointer, save :: tri_parts(:) => null()
  integer, save :: N_tri_part = 0

  public :: ts_tri_init
  public :: N_tri_part, tri_parts

  private
  
contains

  subroutine ts_tri_init()

    use alloc, only : re_alloc, de_alloc
    use parallel, only : IONode
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World, MPI_Comm_Self
#endif

    use class_OrbitalDistribution
    use class_Sparsity
    use create_Sparsity_Union

    use m_ts_electype
    use m_ts_options, only : N_Elec, Elecs, IsVolt
    use m_ts_sparse, only : ts_sp_uc

    use m_ts_method

    use m_ts_Sparsity2TriMat

    use m_ts_tri_scat, only : GFGGF_needed_worksize

    type(OrbitalDistribution) :: dit
    type(Sparsity) :: tmpSp1, tmpSp2

    integer :: idx, no
    integer :: i, els, no_u_TS
    integer :: padding, worksize
    
    no_u_TS = nrows_g(ts_sp_uc) - no_Buf

    ! In order to ensure that the electrodes are in the
    ! tri-diagonal sparsity pattern, we can easily create
    ! the full sparsity pattern with the electrodes included
    ! and then recreate the tri-diagonal sparsity pattern
    ! This is probably the crudest way of doing it.
#ifdef MPI
    call newDistribution(nrows_g(ts_sp_uc),MPI_Comm_Self,dit, &
         name='TranSIESTA UC distribution')
#else    
    call newDistribution(nrows_g(ts_sp_uc),-1,dit, &
         name='TranSIESTA UC distribution')
#endif

    ! This works as creating a new sparsity deletes the previous
    ! and as it is referenced several times it will not be actually
    ! deleted...
    tmpSp2 = ts_sp_uc
    do i = 1 , N_Elec

       idx = Elecs(i)%idx_o
       no = TotUsedOrbs(Elecs(i))

       ! we first create the super-set sparsity
       tmpSp1 = tmpSp2
       call crtSparsity_Union(dit,tmpSp1, &
            idx,idx,no,no, tmpSp2)
    end do
    call delete(tmpSp1)

    ! This will create and even out the parts
    if ( associated(tri_parts) ) &
         call de_alloc(tri_parts, &
         routine='tsSp2TM', name='n_part')

    N_tri_part = 0
    nullify(tri_parts)
    if ( IONode ) &
         write(*,'(/,a)') 'transiesta: Determining an optimal tri-matrix...'
    call ts_Sparsity2TriMat(dit,tmpSp2,N_tri_part,tri_parts)
    call delete(tmpSp2)
    call delete(dit)

    if ( N_tri_part < 3 ) then
       call die('Erroneous transiesta update sparsity pattern. &
            &Check with the developers')
    end if

    ! Calculate size of the tri-diagonal matrix
    els = tri_parts(N_tri_part)**2
    do i = 1 , N_tri_part - 1
       els = els + tri_parts(i)*( tri_parts(i) + 2 * tri_parts(i+1) )
    end do

    if ( IONode ) then
       write(*,'(a)') 'transiesta: Established a near-optimal partition &
            &for the tri-diagonal matrix.'
       write(*,'(a)') 'transiesta: Size of the partitions:'
       do i = 1 , N_tri_part
          write(*,'(tr15,i3,'':'',tr1,i6)') i,tri_parts(i)
       end do
       write(*,'(a,i0,a,i0)') 'transiesta: Matrix elements in tri / full: ', &
            els,' / ',nrows_g(ts_sp_uc)**2
    end if

    if ( .not. IsVolt ) return

    ! Get the padding for the array to hold the entire column
    call GFGGF_needed_worksize(N_tri_part,tri_parts, &
         N_Elec, Elecs, padding, worksize)
    if ( IONode ) then
       write(*,'(a,i0)') 'transiesta: Padding + work size: ',padding + worksize
    end if

    ! Check that we can contain the full column
    if ( maxval(TotUsedOrbs(Elecs)) * no_u_TS > els ) then
       call die('Transiesta tri-diagonal partitioning is too good to &
            &perform GFGGF calculation. Cannot sustain the full column. &
            &Use the sparse solutionmethod instead')
    end if
    
  end subroutine ts_tri_init

end module m_ts_tri_init
