***************
*** 21,27 ****
        use fdf
        use sys,            only : die
        use mesh,           only : xdsp, nsm, nsp
-       use parsing
        use m_iorho,        only : write_rho
        use m_forhar,       only : forhar
        use xcmod,          only : nXCfunc, XCauth
--- 21,26 ----
        use fdf
        use sys,            only : die
        use mesh,           only : xdsp, nsm, nsp
        use m_iorho,        only : write_rho
        use m_forhar,       only : forhar
        use xcmod,          only : nXCfunc, XCauth
***************
*** 33,41 ****
       &                           field, isefld, npcc, uharrs, scell,
       &                           shape, qspiral, g2mesh, nm, nml,
       &                           bcell, nmpl
- #ifdef MPI
-       use mpi_siesta
- #endif
        implicit none
        integer, intent(in)     :: nspin, norb, iaorb(norb), iphorb(norb),
       &                           nuo, nuotot, nua, na, isa(na),
--- 32,37 ----
       &                           field, isefld, npcc, uharrs, scell,
       &                           shape, qspiral, g2mesh, nm, nml,
       &                           bcell, nmpl
        implicit none
        integer, intent(in)     :: nspin, norb, iaorb(norb), iphorb(norb),
       &                           nuo, nuotot, nua, na, isa(na),
***************
*** 52,60 ****
        logical                 :: samesh, samexa
        logical, external       :: leqi
        real(dp), external      :: volcel, ddot
- #ifdef MPI
-       integer                 :: MPIerror
- #endif
        real(grid_p)            :: dummy_Drho(1,1), dummy_Vaux(1),
       &                           dummy_Vscf(1)
        logical, save           :: frstme = .true.   ! Keeps state
--- 48,53 ----
        logical                 :: samesh, samexa
        logical, external       :: leqi
        real(dp), external      :: volcel, ddot
        real(grid_p)            :: dummy_Drho(1,1), dummy_Vaux(1),
       &                           dummy_Vscf(1)
        logical, save           :: frstme = .true.   ! Keeps state
***************
*** 69,82 ****
  
        if (frstme) then
  C       nsm lives in module mesh in meshsubs
-         if (Node.eq.0) then
-           nsm = fdf_integer( 'MeshSubDivisions', 2 )
-           nsm = max( nsm, 1 )
-         endif
- 
- #ifdef MPI
-         call MPI_Bcast( nsm,1,MPI_integer,0,MPI_Comm_World,MPIerror )
- #endif
  
  C       Set mesh sub-division variables & perform one off allocation
          nsp = nsm*nsm*nsm
--- 62,69 ----
  
        if (frstme) then
  C       nsm lives in module mesh in meshsubs
+         nsm = fdf_integer( 'MeshSubDivisions', 2 )
+         nsm = max( nsm, 1 )
  
  C       Set mesh sub-division variables & perform one off allocation
          nsp = nsm*nsm*nsm
***************
*** 442,451 ****
        use precision
        use parallel,      only  : Node, Nodes
        use atmfuncs,      only  : rcut, rcore
-       use fdf
        use sys,           only  : die
        use mesh,          only  : xdsp, nsm, nsp
-       use parsing
        use m_iorho,       only  : write_rho
        use m_forhar,      only  : forhar
        use xcmod,         only  : nXCfunc, XCauth
--- 429,436 ----
        use precision
        use parallel,      only  : Node, Nodes
        use atmfuncs,      only  : rcut, rcore
        use sys,           only  : die
        use mesh,          only  : xdsp, nsm, nsp
        use m_iorho,       only  : write_rho
        use m_forhar,      only  : forhar
        use xcmod,         only  : nXCfunc, XCauth
