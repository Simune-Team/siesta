!
! Header for the interface block for PEXSI
!
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
      comm_global,&
      muMinInertia,&
      muMaxInertia,&
      muLowerEdge,&
      muUpperEdge,&
      inertiaIter,&
      shiftList,&
      inertiaList,&
      info)


   integer, intent(in) :: nrows, nnz, nnzLocal, numColLocal
   integer, intent(in) :: colptrLocal(:), rowindLocal(:)
   real(SELECTED_REAL_KIND(10,100)), intent(in) :: HnzvalLocal(:)
   integer, intent(in) :: isSIdentity
   real(SELECTED_REAL_KIND(10,100)), intent(in) :: SnzvalLocal(:)
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: temperature
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: numElectronExact
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: muMin0, muMax0
   integer, intent(in)                            :: numPole
   integer, intent(in)                            :: inertiaMaxIter
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: inertiaNumElectronTolerance
   
   ! Ordering 
   !   0   : PARMETIS
   !   1   : METIS_AT_PLUS_A
   !   2   : MMD_AT_PLUS_A
   integer, intent(in)                           :: ordering
   
   integer, intent(in)                           :: npPerPole, comm_global
   
   real(SELECTED_REAL_KIND(10,100)), intent(out) :: muMinInertia, muMaxInertia
   real(SELECTED_REAL_KIND(10,100)), intent(out) :: muLowerEdge, muUpperEdge
   integer, intent(out)                          :: inertiaIter
   real(SELECTED_REAL_KIND(10,100)), intent(out) :: shiftList(numPole), inertiaList(numPole)
   integer, intent(out)                          :: info
   
