module m_restruct_SpData2D

  implicit none

  public :: restructdSpData2D
  
contains
  
  subroutine restructdSpData2D(SpMin, sp_out, SpMout, dim2)

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

    ! Limitations:
    !       It is assumed that the number of rows and columns of the
    !       underlying matrix is the same. The number of columns might
    !       conceivably change, but if so, are we guaranteed that the
    !       smaller-numbered columns are the same in both patterns?
    !       (i.e., when the auxiliary supercell changes, does column
    !       number 32 refer to the same atom as before, even if column
    !       235 refers to a completely new image atom?)
    ! 

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

    ! For non-block-cyclic distributions, we are assuming
    ! that the distribution pattern is the same (it might not
    ! be for spatial- or interaction-based distributions).
    ! This test is just a partial sanity check
    ! One should compare nl2g(i) for i 1..nrows

    if (nrows(SpMin) /= nrows(sp_out)) then
       call die("Incompatible nrows in SpMatrices")
    endif


    n_col_in => n_col(SpMin)
    n_col_out => n_col(sp_out)
    listptr_in => list_ptr(SpMin)
    listptr_out => list_ptr(sp_out)
    list_in => list_col(SpMin)
    list_out => list_col(sp_out)
    !print *, &
    !     maxval(abs(n_col_in-n_col_out)), &
    !     maxval(abs(list_in-list_out)),node,nnzs(SpMin),nnzs(sp_out)

    size_in  = nnzs(SpMin)
    size_out = nnzs(sp_out)

    ! We need to check the maximum supercell...
    maxval_j_in  = (maxval(list_in(1:size_in))-1)/nrows_g(SpMin)
    maxval_j_in  = ( maxval_j_in + 1 ) * nrows_g(SpMin)
    maxval_j_out = (maxval(list_out(1:size_out))-1)/nrows_g(sp_out)
    maxval_j_out = ( maxval_j_out + 1 ) * nrows_g(sp_out)

    if ( maxval_j_in /= maxval_j_out ) then
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

    do k = 1, dim2_min
       do i=1,nrows(SpMin)
          aux(:) = 0.0_dp
          do in=1,n_col_in(i)
             ind = listptr_in(i) + in
             j = list_in(ind)
             aux(j) = a_in(ind,k)
          enddo
          do in=1,n_col_out(i)
             ind = listptr_out(i) + in
             j = list_out(ind)
             a_out(ind,k) = aux(j)
          enddo
       enddo
    enddo

    deallocate(aux)

    call newdSpData2D(sp_out,a2d_out,dist(SpMin), &
         SpMout,name="Re-structured SpM")

    call delete(a2d_out)

  contains
    
    subroutine die(str)
      character(len=*), optional :: str
      if (present(str)) then
         print *, trim(str)
      endif
      stop
    end subroutine die

  end subroutine restructdSpData2D

end module m_restruct_SpData2D
