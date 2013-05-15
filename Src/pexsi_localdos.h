!
! Header for the interface block for SI LDOS routine
!
subroutine f_ppexsi_localdos_interface(&
		nrows,&
		nnz,&
		nnzLocal,&
		numColLocal,&
		colptrLocal,&
		rowindLocal,&
		HnzvalLocal,&
		isSIdentity,&
		SnzvalLocal,&
		Energy,&
		broadening,&
		ordering,&
		mpi_Comm,&
		localDOSnzvalLocal,&
		info)

   integer, intent(in) :: nrows, nnz, nnzLocal, numColLocal
   integer, intent(in) :: colptrLocal(:), rowindLocal(:)

   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: HnzvalLocal(:)
   integer, intent(in)                            :: isSIdentity
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: SnzvalLocal(:)

   ! Reference energy and broadening
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: energy
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: broadening

   ! Ordering choice
   !   0   : PARMETIS
   !   1   : METIS_AT_PLUS_A
   !   2   : MMD_AT_PLUS_A
   integer, intent(in)                           :: ordering
   
   integer, intent(in)                           :: mpi_comm
   
   ! Partial DM containing the LDOS info
   real(SELECTED_REAL_KIND(10,100)), intent(out) :: localDOSnzvalLocal(:)

   integer, intent(out)                          :: info

end subroutine f_ppexsi_localdos_interface
