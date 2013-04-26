!
! Header for the interface block for PEXSI
!
  
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
        numElectronDrvList)


   integer, intent(in) :: nrows, nnz, nnzLocal, numColLocal
   integer, intent(in) :: colptrLocal(:), rowindLocal(:)
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: HnzvalLocal(:)
   integer, intent(in)                            :: isSIdentity
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: SnzvalLocal(:)
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: temperature
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: numElectronExact
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: mu0, &
                                                     muMin0, &
                                                     muMax0
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: gap, & 
                                                     deltaE
   integer, intent(in)                            :: numPole
   
   ! Maximum number of allowed iterations
   integer, intent(in)                            :: muMaxIter
   
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: PEXSInumElectronTolerance

   ! Ordering 
   !   0   : PARMETIS
   !   1   : METIS_AT_PLUS_A
   !   2   : MMD_AT_PLUS_A
   integer, intent(in)                           :: ordering
   
   integer, intent(in)                           :: npPerPole, &
                                                    comm_global 
   
   real(SELECTED_REAL_KIND(10,100)), intent(out) :: DMnzvalLocal(:),&
                                                    EDMnzvalLocal(:), &
                                                    FDMnzvalLocal(:)
   real(SELECTED_REAL_KIND(10,100)), intent(out) :: muPEXSI
   real(SELECTED_REAL_KIND(10,100)), intent(out) :: numElectron
   real(SELECTED_REAL_KIND(10,100)), intent(out) :: muMinPEXSI, &
                                                    muMaxPEXSI

   ! Variables related to mu history
   
   ! Actual number of iterations performed
   integer, intent(out)                          :: muIter

   ! List of values of mu, N_e, d(N_e)/d_mu
   real(SELECTED_REAL_KIND(10,100)), intent(out) :: muList(muMaxIter)
   real(SELECTED_REAL_KIND(10,100)), intent(out) :: numElectronList(muMaxIter)
   real(SELECTED_REAL_KIND(10,100)),intent(out)  :: numElectronDrvList(muMaxIter)

