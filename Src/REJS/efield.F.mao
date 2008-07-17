***************
*** 59,71 ****
        use parallelsubs, only : HowManyMeshPerNode
        use atmfuncs,     only : rcut
        use fdf
-       use parsing
        use siesta_cml
        use sys,          only : die
        use mesh,         only : meshLim
- #ifdef MPI
-       use mpi_siesta
- #endif
  
        implicit          none
  
--- 59,67 ----
        use parallelsubs, only : HowManyMeshPerNode
        use atmfuncs,     only : rcut
        use fdf
        use siesta_cml
        use sys,          only : die
        use mesh,         only : meshLim
  
        implicit          none
  
***************
*** 81,101 ****
  C Internal variables
        logical           found, isfield, orthog
        logical,  save :: frstme = .true.
-       character         eunits*10, shape*8, line*130, names*20
        integer           i0(3), i1, i2, i3, ia, imesh, int(1),
       $                  is, iu, ix,
-      .                  j1, j2, j3, last, lc(0:1), 
-      .                  nbcell, ni, nn, nr, nv, meshl(3),
       .                  Yoffset, Zoffset, i30, i20
- #ifdef MPI
-       integer           MPIerror, npl
- #endif
        real(dp)          b1xb2(3), bcell(3,3), cfactor, dplane(3),
       .                  e(3), e0(3), eb1, eb2, eb3,
       .                  f(3), rc, rcell(3,3), v0,
       .                  xfrac, xmax(3), xmean, xmin(3)
        save              e, f, isfield, i0, v0
  
  C Find and store the electric field only the first time
        if (frstme) then
          frstme = .false.
--- 77,97 ----
  C Internal variables
        logical           found, isfield, orthog
        logical,  save :: frstme = .true.
+       character         eunits*10, shape*8
        integer           i0(3), i1, i2, i3, ia, imesh, int(1),
       $                  is, iu, ix,
+      .                  j1, j2, j3, last,
+      .                  nbcell, meshl(3),
       .                  Yoffset, Zoffset, i30, i20
        real(dp)          b1xb2(3), bcell(3,3), cfactor, dplane(3),
       .                  e(3), e0(3), eb1, eb2, eb3,
       .                  f(3), rc, rcell(3,3), v0,
       .                  xfrac, xmax(3), xmean, xmin(3)
        save              e, f, isfield, i0, v0
  
+       type(block_fdf)            :: bfdf
+       type(parsed_line), pointer :: pline
+ 
  C Find and store the electric field only the first time
        if (frstme) then
          frstme = .false.
***************
*** 104,130 ****
  
          if (ionode) then
  C Read the electric field block from the fdf input file
-           found = fdf_block('ExternalElectricField',iu)
            if (found) then
-             read(iu,'(a)') line
-             last = index(line,'#') - 1
-             if (last .le. 0) last = len(line)
-             call parse( line(1:last), nn, lc, names, nv, e,
-      &                ni, int, nr, e0 )
-             eunits = names(lc(0)+1:lc(1))
              cfactor = fdf_convfac(eunits,'Ry/Bohr/e')
              do ix = 1,3
                if (e(ix) .ne. 0.0_dp) isfield = .true.
                e(ix) = e(ix) * cfactor
                e0(ix) = e(ix)
              enddo
            endif
          endif
- #ifdef MPI
-         call MPI_Bcast(e,3,MPI_double_precision,0,MPI_Comm_World,
-      &    MPIerror)
-         call MPI_Bcast(isfield,1,MPI_logical,0,MPI_Comm_World,MPIerror)
- #endif
  
  C Check that the field is orthogonal to the bulk directions
          if (isfield) then
--- 100,120 ----
  
          if (ionode) then
  C Read the electric field block from the fdf input file
+           found = fdf_block('ExternalElectricField',bfdf)
            if (found) then
+             if (.not. fdf_bline(bfdf,pline)) then
+               call die('efield: ExternalElectricField block is empty')
+             endif
+             eunits = fdf_bnames(pline,1)
              cfactor = fdf_convfac(eunits,'Ry/Bohr/e')
              do ix = 1,3
+               e(ix) = fdf_bvalues(pline,ix)
                if (e(ix) .ne. 0.0_dp) isfield = .true.
                e(ix) = e(ix) * cfactor
                e0(ix) = e(ix)
              enddo
            endif
          endif
  
  C Check that the field is orthogonal to the bulk directions
          if (isfield) then
