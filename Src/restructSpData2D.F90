module m_restruct_SpData2D

  implicit none

  public :: restruct_dSpData2D
  
contains
  
  subroutine restruct_dSpData2D(SpMin, sp_out, SpMout, dim2, show_warning)

    use class_dSpData2D
    use class_dData2D
    use class_Sparsity
    use class_OrbitalDistribution

    type(dSpData2D), intent(in) :: SpMin
    type(Sparsity), intent(in) :: sp_out
    ! Note!!  inout is essential to avoid memory leaks...
    type(dSpData2D), intent(inout) :: SpMout
    ! The size of the second dimension in SpMout
    ! This is only used if passed, otherwise it defaults to size(SpMin, 2)
    integer, intent(in), optional :: dim2
    ! Changes the sparsity pattern of a matrix, re-using
    ! as much information as possible. Newly appearing elements are
    ! set to zero. Presumably, elements which dissappear from
    ! the pattern were very close to zero to begin with.
    ! (that could be checked)
    logical, intent(in), optional :: show_warning

    ! Limitations:
    !       This assumes that a prior fixing of the number of columns has
    !       been corrected to account for changing supercell information etc.
    !       I.e. this routine assumes the supercell indices are coherent between
    !       the two sparse patterns.

    integer, parameter :: dp = selected_real_kind(10,100)

    integer :: dim2_in, dim2_min, ldim2
    integer :: i, in, ind, j, k, size_in, size_out
    integer :: maxval_j_out, maxval_j_in, max_col

    type(dData2D)  :: a2d_out

    real(dp), pointer              :: a_in(:,:), a_out(:,:)
    integer, dimension(:), pointer :: n_col_in, n_col_out
    integer, dimension(:), pointer :: listptr_in, listptr_out
    integer, dimension(:), pointer :: list_in, list_out

    real(dp), allocatable  :: aux(:)
    type(OrbitalDistribution) :: dit

    logical :: lshow_warning

    ! For non-block-cyclic distributions, we are assuming
    ! that the distribution pattern is the same (it might not
    ! be for spatial- or interaction-based distributions).
    ! This test is just a partial sanity check
    ! One should compare nl2g(i) for i 1..nrows

    if (nrows(SpMin) /= nrows(sp_out)) then
       call die("Incompatible nrows in SpMatrices")
    endif

    lshow_warning = .true.
    if ( present(show_warning) ) lshow_warning = show_warning

    n_col_in => n_col(SpMin)
    n_col_out => n_col(sp_out)
    listptr_in => list_ptr(SpMin)
    listptr_out => list_ptr(sp_out)
    list_in => list_col(SpMin)
    list_out => list_col(sp_out)

    size_in  = nnzs(SpMin)
    size_out = nnzs(sp_out)

    ! We need to check the maximum supercell...
    maxval_j_in  = (maxval(list_in(1:size_in))-1)/nrows_g(SpMin)
    maxval_j_in  = ( maxval_j_in + 1 ) * nrows_g(SpMin)
    maxval_j_out = (maxval(list_out(1:size_out))-1)/nrows_g(sp_out)
    maxval_j_out = ( maxval_j_out + 1 ) * nrows_g(sp_out)

    if ( maxval_j_in /= maxval_j_out .and. lshow_warning ) then
       ! Print out the different values
       write(*,'(a,tr2,i0,a,tr2,i0)') &
            'WARNING: Connected supercells may have changed: [in]/[out]', &
            maxval_j_in, '/',maxval_j_out
    end if

    ! Maximum "column" index
    max_col = max(maxval_j_in, maxval_j_out)
    allocate(aux(1:max_col))

    a_in => val(SpMin)
    ! This enables the use of the same density matrix in
    ! two different cases
    dim2_in = size(a_in,dim=2)
    if ( present(dim2) ) then
       ldim2 = dim2
    else
       ldim2 = dim2_in
    end if

    ! This is the actual number of components.
    dim2_min = min(dim2_in, ldim2)
    call newdData2D(a2d_out,size_out,ldim2,"(new in restruct)")
    a_out => val(a2d_out)

    ! Use two loops, the outer one to cover the extra dimension
    ! The alternative could be prone to cache misses
    ! (there will be also cache misses due to the sparse de-referencing,
    ! though.
    aux(:) = 0.0_dp

    do k = 1, dim2_min
       do i = 1, nrows(SpMin)
          do ind = listptr_in(i) + 1, listptr_in(i) + n_col_in(i)
            aux(list_in(ind)) = a_in(ind,k)
          end do
          do ind = listptr_out(i) + 1, listptr_out(i) + n_col_out(i)
            a_out(ind,k) = aux(list_out(ind))
          end do
          do ind = listptr_in(i) + 1, listptr_in(i) + n_col_in(i)
            aux(list_in(ind)) = 0._dp
          end do
       end do
    end do

    deallocate(aux)

    ! In cases where SpMout is equal to SpMin it is vital that
    ! we have a pure reference to the object
    dit = dist(SpMin)

    call newdSpData2D(sp_out,a2d_out,dit, &
         SpMout,name="Re-structured SpM")

    call delete(dit)
    call delete(a2d_out)

  contains
    
    subroutine die(str)
      character(len=*), optional :: str
      if (present(str)) then
         print *, trim(str)
      endif
      stop
    end subroutine die

  end subroutine restruct_dSpData2D

end module m_restruct_SpData2D
