 subroutine restructdSpData2D(SpMin,sp_out,SpMout)

   use class_dSpData2D
   use class_dData2D
   use class_Sparsity
   use class_OrbitalDistribution

   type(dSpData2D), intent(in)    :: SpMin
   type(Sparsity), intent(in)    :: sp_out
   ! Note!!  inout is essential to avoid memory leaks...
   type(dSpData2D), intent(inout) :: SpMout

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

      integer :: i, in, ind, j, k, dim2, size_in, size_out
      integer :: maxval_j_old, maxval_j_out, max_col

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

      size_in = nnzs(SpMin)
      size_out = nnzs(sp_out)

      maxval_j_in = maxval(list_in(1:size_in))
      maxval_j_out = maxval(list_out(1:size_out))

      ! Maximum "column" index
      max_col = max(maxval_j_in, maxval_j_out)
      allocate(aux(1:max_col))

      a_in  => val(SpMin)
      dim2 = size(a_in, dim=2)

      call newdData2D(a2d_out,size_out,dim2,"(new in restruct)")

      a_out => val(a2d_out)

      ! Use two loops, the outer one to cover the extra dimension
      ! The alternative could be prone to cache misses
      ! (there will be also cache misses due to the sparse de-referencing,
      ! though.

      do k = 1, dim2
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

      CONTAINS
        subroutine die(str)
          character(len=*), optional :: str
          if (present(str)) then
             print *, trim(str)
          endif
          stop
        end subroutine die

    end subroutine restructdSpData2D
