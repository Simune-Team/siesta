! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine grdsam(nspin, norb, iaorb, iphorb, 
     .                  nuo, nuotot, nua, na, isa, xa, indxua,
     .                  cell, ntm, ifa, istr, maxnd,
     .                  maxnh, numh, listhptr, listh, Dscf, Datm, Hmat,
     .                  Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .                  Exc, Dxc, dipol, fa, stress, fal, stressl )

c ***************************************************************************
c Final call to dhscf of the scf cycle + grid-cell sampling
c
c   After a first call to DHSCF, there are more calls for rigidly shifted
c   atomic coordinates, sampling the small cell defined by the grid, to get
c   average energy, forces, stress, and dipole. It does it with fix density 
c   matrix. It can be regarded as a (discrete sampling) symmetrization
c   to restore the homogeneity of space, which was lost with the 
c   grid summations. 
c
c   The first call to DHSCF is the usual one (the point in which 
c   selfconsistency has been performed). It corresponds to the 
c   shift (0,0,0) in the grid cell.
c
c   The sampling is done on the additional points given in block
c   GridCellSampling, which are given in fractional coordinates of
c   the small grid cell. 
c
c   Even if the sum is large the symmetrization is never complete 
c   since the density matrix is the one obtained for one position of
c   the grid (no new selfconsistency for each shift).
c      
c Written by E. Artacho. January 1998. 
c *********************** INPUT TOWARDS DHSCF *******************************
c integer nspin         : number of spins considered 
c integer norb          : total number of basis orbitals
c integer iaorb(norb)   : atom to which each orbital belongs
c integer iphorb(norb)  : orbital index (within atom) of each orbital
c integer nuo           : number of orbitals in a unit cell (local)
c integer nuotot        : number of orbitals in a unit cell (global)
c integer nua           : number of atoms in unit cell
c integer na            : number of atoms in supercell
c integer isa(na)       : species indexes
c real*8  xa(3,na)      : atomic positions
c real*8  cell(3,3)     : unit cell vectors
c integer ifa           : scf contrib to forces calculated or not
c integer istr          : scf contrib to stress calculated or not
c integer maxnd         : first dimension of dscf
c integer maxnh         : first dimension of listh, hmat
c real*8  dscf(maxnd,nspin): scf DM elements
c real*8  datm(norb)    : Harris DM diagonal elements
c **** DHSCF INPUT OR OUTPUT (DEPENDING ON WHETHER MESH IS CALCULATED)*******
c integer ntm(3)        : number of mesh divisions of each cell vector
c **** DHSCF INPUT OR OUTPUT (DEPENDING ON ARGUMENT ILH) ********************
c integer numh(nuo)     : number of nonzero H elements for each row
c integer listhptr(nuo) : pointer to each row (-1) of H
c integer listh(maxnh)  : nonzero-H-element column indexes for each row
c real*8  Hmat(maxnh,nspin): Hamiltonian matrix in sparse form
c ************************* DHSCF OUTPUT ******i*****************************
c real*8  enaatm : integral of vna * rhoatm
c real*8  enascf : integral of vna * rhoscf
c real*8  uatm   : Harris hartree electron-interaction energy
c real*8  uscf   : scf hartree electron-interaction energy
c real*8  duscf  : electrostatic (hartree) energy of rhoscf-rhoatm density
c real*8  duext  : interaction energy with external electric field
c real*8  exc    : scf XC energy
c real*8  dxc    : scf double-counting correction to exc
c real*8  dipol(3): electric dipole (in a.u.)
c *********************** DHSCF INPUT AND OUTPUT ****************************
c real*8  fa(3,na)    : atomic forces
c real*8  stress(3,3) : stress tensor
c real*8  fal(3,na)   : atomic forces (local to Node)
c real*8  stressl(3,3): stress tensor (local to Node)
c ***************************************************************************

C
C  Modules
C
      use precision
      use fdf
      use parallel,  only : Node
      use sys,       only : die
      use m_dhscf,   only : dhscf
      use alloc,     only : re_alloc, de_alloc

      implicit          none

      integer           maxnh, maxnd, na, norb, nspin, nuo, nuotot, 
     .                  nua, iaorb(norb), ifa, indxua(na), 
     .                  iphorb(norb), isa(na), istr, listh(maxnh),
     .                  listhptr(nuo), numh(nuo), ntm(3)

      real(dp)          cell(3,3), datm(norb), dipol(3), 
     .                  dscf(maxnd,nspin), DUscf, DUext, Dxc, 
     .                  enaatm, enascf, exc, fa(3,nua),
     .                  Hmat(maxnh,nspin), stress(3,3),
     .                  uatm, uscf, xa(3,na), fal(3,nua), stressl(3,3)

C Internal variables and arrays
 
      logical           samesh
      integer           nr, ipt, maxpt, ia,
     .                  ix, iv
      real(dp)          reals(4), avexc,
     .                  strold(3,3), avdipo(3), avstre(3,3), avdxc,
     .                  avenaa, avenas, avuatm, avuscf, avdusc, anpt,
     .                  strold2(3,3), avstre2(3,3)

      parameter         (maxpt = 100)

      real(dp), pointer, save :: pt(:,:), avfa(:,:), faold(:,:),
     &                           xanew(:,:), avfa2(:,:), faold2(:,:),
     &                           dpt(:,:)

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline


      !Saved variables:
      logical, save :: blread = .false.
      logical, save :: samgrd = .false.

      real(dp), save :: lstcll(3,3) = 0.743978657912656D50
      integer, save  :: lstntm(3) = (/10000, 1, 1/) 
      integer, save  :: npt

c ----------------------------------------------------------------------------

C Allocate local memory
      call re_alloc( pt, 1, 3, 1, maxpt+1, 'pt', 'grdsam' )
      call re_alloc( dpt, 1, 3, 1, maxpt+1, 'dpt', 'grdsam' )

      call re_alloc( avfa, 1, 3, 1, na, 'avfa', 'grdsam' )
      call re_alloc( avfa2, 1, 3, 1, na, 'avfa2', 'grdsam' )
      call re_alloc( faold, 1, 3, 1, na, 'faold', 'grdsam' )
      call re_alloc( faold2, 1, 3, 1, na, 'faold2', 'grdsam' )
      call re_alloc( xanew, 1, 3, 1, na, 'xanew', 'grdsam' )

C Check whether GridCellSampling block has been looked for ------------------

      if ( .not. blread ) then

C Look for block and read it if found ---------------------------------------
 
        samgrd = fdf_block('GridCellSampling',bfdf)

        if (samgrd) then
          if (Node.eq.0) then
            write(6,'(a)') 'grdsam: Reading %block GridCellSampling'
          endif
          ipt = 1
          do while(fdf_bline(bfdf,pline) .and. (ipt .le. maxpt+1))
            nr = fdf_bnreals(pline)
            if (nr .ge. 3) then
              pt(1,ipt) = fdf_breals(pline,1)
              pt(2,ipt) = fdf_breals(pline,2)
              pt(3,ipt) = fdf_breals(pline,3)
              ipt = ipt + 1
            else
              write(6,'(a,/,a)')
     &          'grdsam: syntax ERROR in %block GridCellSampling:',
     &          trim(fdf_getline(bfdf%mark))
            endif
          enddo
          if (ipt .gt. maxpt+1) then
            call die('grdsam: ERROR: Number of points' //
     $               ' in block larger than MAXPT')
          endif
          npt = ipt - 1
        endif
        blread = .true.
      endif

c if sampling, store the forces and stresses prior to dhscf -----------------

      if (samgrd) then
        faold(1:3,1:nua) = fa(1:3,1:nua)
        faold2(1:3,1:nua) = fal(1:3,1:nua)
        strold(1:3,1:3) = stress(1:3,1:3)
        strold2(1:3,1:3) = stressl(1:3,1:3)
        if (Node.eq.0) then
          write(6,'(a,i3)') 'grdsam: Grid-cell sampling, point ', 0
        endif
      endif


c first (normal) call to DHSCF, ihmat=1 --------------------------------------

      call dhscf( nspin, norb, iaorb, iphorb, nuo, nuotot, nua,
     .            na, isa, xa, indxua, ntm,
     .            ifa, istr, 1, ' ', ' ', ' ', ' ', ' ', ' ',
     .            maxnd, numh, listhptr, listh, Dscf, Datm,
     .            maxnh, Hmat, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .            Exc, Dxc, dipol, stress, fal, stressl)

c Sampling ------------------------------------------------------------------

      if (samgrd) then

c find if mesh is new -------------------------------------------------------

        samesh = .true.
        do iv = 1, 3
          if ( ntm(iv) .ne. lstntm(iv) ) samesh = .false.
          lstntm(iv) = ntm(iv)
          do ix = 1, 3
            if ( cell(ix,iv) .ne. lstcll(ix,iv) ) samesh = .false.
            lstcll(ix,iv) = cell(ix,iv)
          enddo
        enddo

c generate displacements (in Bohr cartesian) if first time or new mesh ------

        if (.not. samesh) then
          if (Node.eq.0) then
            write(6,'(a)') 
     .       'grdsam: Generating displacements for grid-cell sampling'
          endif
          do ipt = 1, npt
            do ix = 1, 3
              dpt(ix,ipt) = 0.0d0
              do iv = 1, 3
                dpt(ix,ipt) = dpt(ix,ipt) + 
     .                          cell(ix,iv)*pt(iv,ipt)/dble(ntm(iv))
              enddo
            enddo
          enddo
        endif

c initialize averages with output of regular DHSCF run -----------------------

        avfa(1:3,1:nua) = fa(1:3,1:nua)
        avfa2(1:3,1:nua) = fal(1:3,1:nua)
        avstre(1:3,1:3) = stress(1:3,1:3)
        avstre2(1:3,1:3) = stressl(1:3,1:3)
        avdipo(1:3) = dipol(1:3)
        avenaa = enaatm
        avenas = enascf
        avuatm = uatm
        avuscf = uscf
        avdusc = duscf
        avexc  = exc
        avdxc  = dxc

c loop sampling on displacements making averages ----------------------------

        do ipt = 1, npt

          do ia = 1, na
            do ix = 1, 3
              xanew(ix,ia) = xa(ix,ia) + dpt(ix,ipt)
            enddo
          enddo
          fa(1:3,1:nua) = faold(1:3,1:nua)
          fal(1:3,1:nua) = faold2(1:3,1:nua)
          stress(1:3,1:3) = strold(1:3,1:3)
          stressl(1:3,1:3) = strold2(1:3,1:3)
        
          if (Node.eq.0) then
          write(6,'(a,i3)') 'grdsam: Grid-cell sampling, point ', ipt
          endif

          call dhscf( nspin, norb, iaorb, iphorb, nuo, nuotot, nua,
     .                na, isa, xanew, indxua, ntm,
     .                ifa, istr, 0, ' ', ' ', ' ', ' ', ' ', ' ',
     .                maxnd, numh, listhptr, listh, Dscf, Datm,
     .                maxnh, Hmat, Enaatm, Enascf, Uatm, Uscf, DUscf,
     .                DUext, Exc, Dxc, dipol, stress, fal, stressl)

          avfa(1:3,1:nua) = avfa(1:3,1:nua) + fa(1:3,1:nua)
          avfa2(1:3,1:nua) = avfa2(1:3,1:nua) + fal(1:3,1:nua)
          avstre(1:3,1:3) = avstre(1:3,1:3) + stress(1:3,1:3)
          avstre2(1:3,1:3) = avstre2(1:3,1:3) + stressl(1:3,1:3)
          avdipo(1:3) = avdipo(1:3) + dipol(1:3)
          avenaa = avenaa + enaatm
          avenas = avenas + enascf
          avuatm = avuatm + uatm
          avuscf = avuscf + uscf
          avdusc = avdusc + duscf
          avexc  = avexc  + exc
          avdxc  = avdxc  + dxc
        enddo

c final averages ------------------------------------------------------------

        anpt = dble(npt+1)

        fa(1:3,1:nua) = avfa(1:3,1:nua)/anpt
        fal(1:3,1:nua) = avfa2(1:3,1:nua)/anpt
        stress(1:3,1:3) = avstre(1:3,1:3)/anpt
        stressl(1:3,1:3) = avstre2(1:3,1:3)/anpt
        dipol(1:3) = avdipo(1:3)/anpt
        enaatm = avenaa/anpt
        enascf = avenas/anpt
        uatm   = avuatm/anpt
        uscf   = avuscf/anpt
        duscf  = avdusc/anpt
        exc    = avexc/anpt
        dxc    = avdxc/anpt

      endif

C Deallocate local memory
      call de_alloc( xanew, 'xanew', 'grdsam' )
      call de_alloc( faold2, 'faold2', 'grdsam' )
      call de_alloc( faold, 'faold', 'grdsam' )
      call de_alloc( avfa2, 'avfa2', 'grdsam' )
      call de_alloc( avfa, 'avfa', 'grdsam' )

      return
      end

