***************
*** 16,27 ****
        use atomlist,           only: no_u, rmaxkb, amass, lasto, qtot,
       &                              iza, rmaxo, zvaltot, superc,
       &                              initatomlists, no_l
-       use m_fdf_global,       only: fdf_global_get
        use sys,                only: die, bye
        use xcmod,              only: setXC
        use molecularmechanics, only: inittwobody
        use metaforce,          only: initmeta
-       use m_mpi_utils,        only: broadcast
        use alloc,              only: re_alloc, alloc_report
        use phonon,             only: phonon_num_disps, phonon_setup
        use parallelsubs,       only: getnodeorbs
--- 16,26 ----
        use atomlist,           only: no_u, rmaxkb, amass, lasto, qtot,
       &                              iza, rmaxo, zvaltot, superc,
       &                              initatomlists, no_l
+       use fdf,                only: fdf_parallel
        use sys,                only: die, bye
        use xcmod,              only: setXC
        use molecularmechanics, only: inittwobody
        use metaforce,          only: initmeta
        use alloc,              only: re_alloc, alloc_report
        use phonon,             only: phonon_num_disps, phonon_setup
        use parallelsubs,       only: getnodeorbs
***************
*** 75,80 ****
        call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
        call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
        call debugMpiOn( )
  #endif
  
        IOnode = (Node .eq. 0)
--- 74,84 ----
        call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
        call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
        call debugMpiOn( )
+ 
+       if (.not. fdf_parallel()) then
+         call die('siesta_init: ERROR: FDF module doesn''t have ' //
+                  'parallel support')
+       endif
  #endif
  
        IOnode = (Node .eq. 0)
***************
*** 129,135 ****
        call siesta_cml_init()
  
  C     Set allocation report level
-       call fdf_global_get( level, 'alloc_report_level', 0 )
        call alloc_report( level=level, file=trim(slabel)//'.alloc' )
  
  C     Initialise exchange-correlation functional information
--- 133,139 ----
        call siesta_cml_init()
  
  C     Set allocation report level
+       call fdf_get( level, 'alloc_report_level', 0 )
        call alloc_report( level=level, file=trim(slabel)//'.alloc' )
  
  C     Initialise exchange-correlation functional information
***************
*** 154,168 ****
        call initatomlists( )    ! Sets iza
  
  C     early exit if only checking the structure
-       call fdf_global_get(struct_only,'Output-Structure-Only',.false.)
        if (IONode) then
-          call write_struct( ucell, na_u, isa, iza, xa )
-          if (lUseZmatrix) then
-             call write_canonical_ucell_and_Zmatrix()
-          endif
        endif
        if (struct_only) then
-          call bye("End of structure processing")
        endif
        
  C     End of Initial Structure Processing
--- 158,172 ----
        call initatomlists( )    ! Sets iza
  
  C     early exit if only checking the structure
+       call fdf_get(struct_only,'Output-Structure-Only',.false.)
        if (IONode) then
+         call write_struct( ucell, na_u, isa, iza, xa )
+         if (lUseZmatrix) then
+           call write_canonical_ucell_and_Zmatrix()
+         endif
        endif
        if (struct_only) then
+         call bye("End of structure processing")
        endif
        
  C     End of Initial Structure Processing
***************
*** 249,255 ****
        numhold(:)       = 0
  
  C     Get number of eigenstates that need to be calculated
-       call fdf_global_get(neigwanted,'NumberOfEigenStates',no_u)
  
  C     Check number of eigenstates - cannot be larger than number of
  C     basis functions or smaller than number of occupied states + 1
--- 253,259 ----
        numhold(:)       = 0
  
  C     Get number of eigenstates that need to be calculated
+       call fdf_get(neigwanted,'NumberOfEigenStates',no_u)
  
  C     Check number of eigenstates - cannot be larger than number of
  C     basis functions or smaller than number of occupied states + 1
