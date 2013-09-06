module m_pexsi_interface

! Interface blocks for PEXSI library, using the ISO_C_BINDING
! mechanism, included in the F2003 standard and implemented
! by most compilers.

! Valid for PEXSI library versions >= 0.4.6

! Alberto Garcia, September 2013

! Note that the current use of the FORTRAN() macro in the original C++
! source implies that the interface specification needs to specify an
! explicit C name with an extra underscore in the "BIND" section. This
! could be changed in future versions of the library.

!
! Dummy C routines compatible with this interface are available
! in Src/dummy_pexsi/dummy_pexsi.c
!
!--------------------------------------------------------------------
!
! Technical notes:
!
! The C-interop mechanism guarantees that the data types are consistent,
! by using the kind names in the ISO_C_BINDING module (i.e., C_INT, C_DOUBLE)
!
! Array arguments are *required* by the standard to be explicit shape
! (e.g. a(3)), or assumed size (e.g. a(*)). This avoids problems with
! the assumed shape specification (e.g. a(:)), which would need to
! pass extra information in general. This permits the use of pointers
! and array sections as actual arguments. The compiler would
! presumably make on-the-fly copies in the case of non-contiguous data
! (a situation that should be avoided on performance grounds).
!
!---------------------------
! Note that not all arguments are properly documented in this file.
!
!-----------------------------------------------------------------------------
!
! Interface block for PEXSI solver
!
interface  
subroutine f_ppexsi_solve_interface(&
        nrows,&
        nnz,&
        nnzLocal,&
        numColLocal,&
        colptrLocal,&
        rowindLocal,&
        HnzvalLocal,&
        isSIdentity,&
        SnzvalLocal,&
        temperature,&
        numElectronExact,&
        mu0,&
        muMin0,&
        muMax0,&
        gap,&
        deltaE,&
        numPole,&
        muMaxIter,&
        PEXSINumElectronTolerance,&
        ordering,&
        npPerPole,&
	npSymbFact,&
        comm_global,&
        DMnzvalLocal,&
        EDMnzvalLocal,&
        FDMnzvalLocal,&
        muPEXSI,&
        numElectron,&
        muMinPEXSI,&
        muMaxPEXSI,&
        muIter,&
        muList,&
        numElectronList,&
        numElectronDrvList,&
        info) &
   BIND(C, Name="f_ppexsi_solve_interface_")

   USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT,C_DOUBLE


   integer(C_INT), intent(in)   :: nrows, nnz, nnzLocal, numColLocal
   integer(C_INT), intent(in)   :: colptrLocal(*), rowindLocal(*)
   real(C_DOUBLE), intent(in)   :: HnzvalLocal(*)
   integer(C_INT), intent(in)   :: isSIdentity
   real(C_DOUBLE), intent(in)   :: SnzvalLocal(*)
   real(C_DOUBLE), intent(in)   :: temperature
   real(C_DOUBLE), intent(in)   :: numElectronExact
   real(C_DOUBLE), intent(in)   :: mu0, &         
                                   muMin0, &
                                   muMax0
   real(C_DOUBLE), intent(in)   :: gap, & 
                                   deltaE
   integer(C_INT), intent(in)   :: numPole
   
   ! Maximum number of allowed iterations
   integer(C_INT), intent(in)   :: muMaxIter
   
   real(C_DOUBLE), intent(in)   :: PEXSInumElectronTolerance

   ! Ordering 
   !   0   : PARMETIS
   !   1   : METIS_AT_PLUS_A
   !   2   : MMD_AT_PLUS_A
   integer(C_INT), intent(in)   :: ordering
   
   integer(C_INT), intent(in)   :: npPerPole, &
                                   comm_global 
   
   ! Number of processors used for symbolic factorization
   ! (Maximum: npPerPole)
   ! Only relevant if PARMETIS/PT-SCOTCH is used.
   integer(C_INT), intent(in)   :: npSymbFact

   real(C_DOUBLE), intent(out)  :: DMnzvalLocal(*),&
                                   EDMnzvalLocal(*), &
                                   FDMnzvalLocal(*)
   real(C_DOUBLE), intent(out)  :: muPEXSI
   real(C_DOUBLE), intent(out)  :: numElectron
   real(C_DOUBLE), intent(out)  :: muMinPEXSI, &
                                   muMaxPEXSI

   ! Variables related to mu history
   
   ! Actual number of iterations performed
   integer(C_INT), intent(out)  :: muIter

   ! List of values of mu, N_e, d(N_e)/d_mu
   real(C_DOUBLE), intent(out)  :: muList(muMaxIter)
   real(C_DOUBLE), intent(out)  :: numElectronList(muMaxIter)
   real(C_DOUBLE), intent(out)  :: numElectronDrvList(muMaxIter)

   integer(C_INT), intent(out)  :: info

end subroutine f_ppexsi_solve_interface
end interface

!
! Inertia Counts
!
interface
subroutine f_ppexsi_inertiacount_interface( &
      nrows,&
      nnz,&
      nnzLocal,&
      numColLocal,&
      colptrLocal,&
      rowindLocal,&
      HnzvalLocal,&
      isSIdentity,&
      SnzvalLocal,&
      temperature,&
      numElectronExact,&
      muMin0,&
      muMax0,&
      numPole,&
      inertiaMaxIter,&
      inertiaNumElectronTolerance,&
      ordering,&
      npPerPole,&
      npSymbFact,&
      comm_global,&
      muMinInertia,&
      muMaxInertia,&
      muLowerEdge,&
      muUpperEdge,&
      inertiaIter,&
      shiftList,&
      inertiaList,&
      info) &
   BIND(C,Name="f_ppexsi_inertiacount_interface_")

   USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT,C_DOUBLE


   integer(C_INT), intent(in) :: nrows, nnz, nnzLocal, numColLocal
   integer(C_INT), intent(in) :: colptrLocal(*), rowindLocal(*)
   real(C_DOUBLE), intent(in) :: HnzvalLocal(*)
   integer(C_INT), intent(in) :: isSIdentity
   real(C_DOUBLE), intent(in) :: SnzvalLocal(*)
   real(C_DOUBLE), intent(in) :: temperature
   real(C_DOUBLE), intent(in) :: numElectronExact
   real(C_DOUBLE), intent(in) :: muMin0, muMax0
   integer(C_INT), intent(in) :: numPole
   integer(C_INT), intent(in) :: inertiaMaxIter
   real(C_DOUBLE), intent(in) :: inertiaNumElectronTolerance
   
   ! Ordering 
   !   0   : PARMETIS
   !   1   : METIS_AT_PLUS_A
   !   2   : MMD_AT_PLUS_A
   integer(C_INT), intent(in)  :: ordering
   
   integer(C_INT), intent(in)  :: npPerPole, comm_global
   
   ! Number of processors used for symbolic factorization
   ! (Maximum: npPerPole)
   ! Only relevant if PARMETIS/PT-SCOTCH is used.

   integer(C_INT), intent(in)  :: npSymbFact

   real(C_DOUBLE), intent(out) :: muMinInertia, muMaxInertia
   real(C_DOUBLE), intent(out) :: muLowerEdge, muUpperEdge
   integer(C_INT), intent(out) :: inertiaIter
   real(C_DOUBLE), intent(out) :: shiftList(numPole), inertiaList(numPole)
   integer(C_INT), intent(out) :: info
end subroutine f_ppexsi_inertiacount_interface
end interface

!
! Interface block for selected-inversion LDOS routine
!
interface

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
		npSymbFact,&
		mpi_Comm,&
		localDOSnzvalLocal,&
		info)  &
   BIND(C, Name="f_ppexsi_localdos_interface_")

   USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT,C_DOUBLE

   integer(C_INT), intent(in)   :: nrows, nnz, nnzLocal, numColLocal
   integer(C_INT), intent(in)   :: colptrLocal(*), rowindLocal(*)

   real(C_DOUBLE), intent(in)   :: HnzvalLocal(*)
   integer(C_INT), intent(in)   :: isSIdentity
   real(C_DOUBLE), intent(in)   :: SnzvalLocal(*)

   ! Reference energy and broadening
   real(C_DOUBLE), intent(in)   :: energy
   real(C_DOUBLE), intent(in)   :: broadening

   ! Ordering choice
   !   0   : PARMETIS
   !   1   : METIS_AT_PLUS_A
   !   2   : MMD_AT_PLUS_A
   integer(C_INT), intent(in)   :: ordering

   ! Number of processors used for symbolic factorization
   ! (Maximum: number of procs in mpi_comm)
   integer(C_INT), intent(in)   :: npSymbFact
   
   integer(C_INT), intent(in)   :: mpi_comm
   
   ! Partial DM containing the LDOS info
   real(C_DOUBLE), intent(out)  :: localDOSnzvalLocal(*)

   integer(C_INT), intent(out)  :: info

end subroutine f_ppexsi_localdos_interface
end interface
   
end module m_pexsi_interface
