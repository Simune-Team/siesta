***************
*** 24,32 ****
        use parallel,  only: Node
        use fdf,       only: fdf_boolean
        use MeshComm,  only: setnewMeshLim
- #ifdef MPI
-       use mpi_siesta
- #endif
        implicit none
  C Passed arguments
        integer, intent(in)  :: nmpl, norb, iaorb(norb), iphorb(norb),
--- 24,29 ----
        use parallel,  only: Node
        use fdf,       only: fdf_boolean
        use MeshComm,  only: setnewMeshLim
        implicit none
  C Passed arguments
        integer, intent(in)  :: nmpl, norb, iaorb(norb), iphorb(norb),
***************
*** 35,43 ****
        integer              :: ia, io, iop, ip, ip0, iphi, is, isp, n,
       &                        nlist, nliste
  
- #ifdef MPI
-       integer              :: MPIerror
- #endif
        logical              :: within
        logical,       save  :: firsttime = .true.
        real(dp)             :: dxsp(3,nsp), grphi(3), phip, r2o,
--- 32,37 ----
        integer              :: ia, io, iop, ip, ip0, iphi, is, isp, n,
       &                        nlist, nliste
  
        logical              :: within
        logical,       save  :: firsttime = .true.
        real(dp)             :: dxsp(3,nsp), grphi(3), phip, r2o,
***************
*** 47,59 ****
  !------------------------------------------------------------------------- BEGIN
  C     On first call set the logical DirectPhi
        if (firsttime) then
-         if (Node.eq.0) then
-           DirectPhi  = fdf_boolean( 'DirectPhi', .false. )
-         endif
- #ifdef MPI
-         call MPI_Bcast( DirectPhi, 1, MPI_logical, 0, MPI_Comm_World,
-      &                  MPIerror )
- #endif
          firsttime = .false.
        endif
  
--- 41,47 ----
  !------------------------------------------------------------------------- BEGIN
  C     On first call set the logical DirectPhi
        if (firsttime) then
+         DirectPhi  = fdf_boolean( 'DirectPhi', .false. )
          firsttime = .false.
        endif
  
