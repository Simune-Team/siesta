      module m_sparse

      use precision,    only:  dp
      use alloc,        only : re_alloc, de_alloc

      public :: same_sparsity, change_sparsity
      private

      CONTAINS

      subroutine change_sparsity (nrows,ncol_old,listptr_old,list_old,
     $                                  ncol_new,listptr_new,list_new,
     $                                  m_old, m_new,
     $                            work)

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
      ! If present, work is supposed to be a work pointer.
      ! It is the responsibility of the caller to nullify it before
      ! entry and to deallocate it when it is no longer needed.

      integer, intent(in) :: nrows
      integer, intent(in) :: ncol_old(nrows), ncol_new(nrows)
      integer, intent(in) :: listptr_old(nrows), listptr_new(nrows)
      integer, intent(in) :: list_old(:), list_new(:)
      real(dp), intent(in)  :: m_old(:)
      real(dp), intent(out) :: m_new(:)

      real(dp), pointer, optional :: work(:)

      integer :: size_old, size_new, i, in, ind, j
      integer :: maxval_j_old, maxval_j_new, max_col

      real(dp), pointer :: aux(:)

      size_old = size(m_old)
      size_new = size(m_new)

      maxval_j_old = maxval(list_old(:))
      maxval_j_new = maxval(list_new(:))

      max_col = max(maxval_j_old, maxval_j_new)

      if (present(work)) then
         call re_alloc(work,1,max_col)
         aux => work
      else
         nullify(aux)
         call re_alloc(aux,1,max_col)
      endif

      aux(1:max_col) = 0.0_dp
      do i=1,nrows
         do in=1,ncol_old(i)
            ind = listptr_old(i) + in
            j = list_old(ind)
            aux(j) = m_old(ind)
         enddo
         do in=1,ncol_new(i)
            ind = listptr_new(i) + in
            j = list_new(ind)
            m_new(ind) = aux(j)
         enddo
         aux(1:max_col) = 0.0_dp
      enddo

      if (.not. present(work)) then
         call de_alloc(aux)
      endif

      end subroutine change_sparsity

      function same_sparsity (nrows,ncol_old,list_old,
     $                              ncol_new,list_new)
     $         result(same_sp)

      ! Checks whether two sparsity patterns are the same

      integer, intent(in) :: nrows
      integer, intent(in) :: ncol_old(nrows), ncol_new(nrows)
      integer, intent(in) :: list_old(:), list_new(:)

      integer :: size_old, size_new, i

      same_sp = .true.
      if (nrows == 0) RETURN

      size_old = size(list_old)
      size_new = size(list_new)

      if (size_old /= size_new) then
         same_sp = .false.
         RETURN
      endif

      do i = 1,nrows
         if (ncol_old(i) .ne. ncol_new(i)) then
            same_sp = .false.
            RETURN
         endif
      enddo

      do i = 1,size_old
         if (list_old(i) .ne. list_new(i))  then
            same_sp = .false.
            RETURN
         endif
      enddo

      end function same_sparsity
      end module m_sparse
