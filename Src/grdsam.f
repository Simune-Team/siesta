c $Id: grdsam.f,v 1.5 1999/02/26 10:06:49 wdpgaara Exp $

      subroutine grdsam(nspin, maxorb, norb, iaorb, iphorb, indxuo,
     .                  nua, na, isa, xa, indxua,
     .                  cell, mscell, g2max, ntm, 
     .                  ilh, ifa, istr,
     .                  maxno, numh, listh, Dscf, Datm, Hmat,
     .                  Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .                  Exc, Dxc, dipol, fa, stress, ierror )

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
c integer maxorb        : second dimension of dscf and hmat
c integer norb          : total number of basis orbitals
c integer iaorb(norb)   : atom to which each orbital belongs
c integer iphorb(norb)  : orbital index (within atom) of each orbital
c integer indxuo(norb)  : index of equivalent atom in unit cell
c integer nua           : number of atoms in unit cell
c integer na            : number of atoms in supercell
c integer isa(na)       : species indexes
c real*8  xa(3,na)      : atomic positions
c real*8  cell(3,3)     : unit cell vectors
c integer mscell(3,3)   : supercell vectors in units of ucell
c integer ilh           : numh and listh: input or output?
c integer ifa           : scf contrib to forces calculated or not
c integer istr          : scf contrib to stress calculated or not
c integer maxno         : first dimension of listh, dscf, hmat
c real*8  dscf(maxno,maxorb,nspin): scf DM elements
c real*8  datm(norb)        : Harris DM diagonal elements
c **** DHSCF INPUT OR OUTPUT (DEPENDING ON WHETHER MESH IS CALCULATED)*******
c integer ntm(3)            : number of mesh divisions of each cell vector
c **** DHSCF INPUT OR OUTPUT (DEPENDING ON ARGUMENT ILH) ********************
c integer numh(norb)        : number of nonzero H elements for each row
c integer listh(maxno,norb) : nonzero-H-element column indexes for each row
c real*8  hmat(maxno,maxorb,nspin): Hamiltonian matrix in sparse form
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
c integer ierror : ierror=0 => no error occurred
c                  ierror=1 => bad internal dimensions. Recompile 
c *********************** DHSCF INPUT AND OUTPUT ****************************
c real*8  g2max       : effective planewave cutoff in Ry 
c real*8  fa(3,na)    : atomic forces
c real*8  stress(3,3) : stress tensor
c ***************************************************************************

      implicit          none

      integer           maxno, maxorb, na, norb, nspin, nua,
     .                  iaorb(norb), ierror, ifa, ilh,
     .                  indxua(na), indxuo(norb),
     .                  iphorb(norb), isa(na), istr, listh(maxno,norb),
     .                  mscell(3,3), numh(norb), ntm(3)

      double precision  cell(3,3), datm(norb), dipol(3), 
     .                  dscf(maxno,maxorb,nspin), DUscf, DUext, Dxc, 
     .                  enaatm, enascf, exc, fa(3,nua), g2max,
     .                  hmat(maxno,maxorb,nspin), stress(3,3),
     .                  uatm, uscf, xa(3,na)

      external          parse

c ----------------------------------------------------------------------------

c Internal variables and arrays
 
      character         line*130, names*80
      logical           samgrd, blread, samesh
      integer           ni, nn, nr, nv, npt, ipt, maxpt, maxat, iu, ia,
     .                  ix, iv, integs(4), lastc, lc(0:3), lstntm(3)
      double precision  reals(4), values(4), lstcll(3,3), avexc,
     .                  strold(3,3), avdipo(3), avstre(3,3), avdxc,
     .                  avenaa, avenas, avuatm, avuscf, avdusc, anpt

      parameter         (maxpt = 100, maxat = 5000)
      double precision  pt(3,maxpt), dpt(3,maxpt), 
     .                  avfa(3,maxat), faold(3,maxat), xanew(3,maxat)

      save              blread, samgrd, lstcll, lstntm, npt, dpt
    

c enable FDF input/output

      include 'fdf/fdfdefs.h'


      data              blread      /.false./,
     .                  samgrd      /.false./,
     .                  lstcll(1,1) / 0.743978657912656D50 /,
     .                  lstntm      / 10000, 1, 1 /

c ----------------------------------------------------------------------------

c check whether GridCellSampling block has been looked for ------------------

      if ( .not. blread ) then

c look for block and read it if found ---------------------------------------
 
        samgrd = fdf_block('GridCellSampling',iu)

        if (samgrd) then
          write(6,'(a)') 'grdsam: Reading %block GridCellSampling'
          do ipt = 1, maxpt + 1
            read(iu,'(a)', end=50) line
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )
            if (nr .ge. 3) then
              pt(1,ipt) = reals(1)
              pt(2,ipt) = reals(2)
              pt(3,ipt) = reals(3)
            else
              npt = ipt - 1
              goto 50
            endif
          enddo
          write(6,'(/,a)') 
     .     'grdsam: ERROR: Number of points in block larger than MAXPT'
          stop 
     .     'grdsam: ERROR: Number of points in block larger than MAXPT'
   50     continue
        endif
        blread = .true.
      endif

c if sampling, store the forces and stresses prior to dhscf -----------------

      if (samgrd) then
        if (maxat .lt. na) then
          write(6,"(/,a,i5)") 
     .    'grdsam: ERROR: Parameter MAXAT should be increased to ', na
          stop 'grdsam: ERROR: Insufficient MAXAT'
        endif
        do ia = 1, nua
          do ix = 1, 3
            faold(ix,ia) = fa(ix,ia)
          enddo
        enddo
        do iv = 1, 3
          do ix = 1, 3
            strold(ix,iv) = stress(ix,iv)
          enddo
        enddo
        write(6,'(a,i3)') 'grdsam: Grid-cell sampling, point ', 0
      endif


c first (normal) call to DHSCF, ihmat=1 --------------------------------------

      call dhscf( nspin, maxorb, norb, iaorb, iphorb, indxuo,
     .            na, isa, xa, indxua, cell, mscell, g2max, ntm,
     .            ilh, ifa, istr, 1, ' ', ' ', ' ', ' ',
     .            maxno, numh, listh, Dscf, Datm,
     .            maxno, numh, listh, Hmat,  
     .            Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .            Exc, Dxc, dipol, fa, stress, ierror )
      if (ierror .ne. 0) return


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
          write(6,'(a)') 
     .       'grdsam: Generating displacements for grid-cell sampling'
          do ipt = 1, npt
            do ix = 1, 3
              dpt(ix,ipt) = 0.d0
              do iv = 1, 3
                dpt(ix,ipt) = dpt(ix,ipt) + 
     .                          cell(ix,iv)*pt(iv,ipt)/dble(ntm(iv))
              enddo
            enddo
          enddo
        endif

c initialize averages with output of regular DHSCF run -----------------------

        do ia = 1, nua
          do ix = 1, 3
            avfa(ix,ia) = fa(ix,ia)
          enddo
        enddo
        do iv = 1, 3
          do ix = 1, 3
            avstre(ix,iv) = stress(ix,iv)
          enddo
        enddo
        do ix = 1, 3
          avdipo(ix) = dipol(ix)
        enddo
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
          do ia = 1, nua
            do ix = 1, 3
              xanew(ix,ia) = xa(ix,ia) + dpt(ix,ipt)
            enddo
          enddo
          do iv = 1, 3
            do ix = 1, 3
              stress(ix,iv) = strold(ix,iv)
            enddo
          enddo
        
          write(6,'(a,i3)') 'grdsam: Grid-cell sampling, point ', ipt

          call dhscf( nspin, maxorb, norb, iaorb, iphorb, indxuo,
     .                na, isa, xanew, indxua, cell, mscell, g2max, ntm,
     .                ilh, ifa, istr, 0, ' ', ' ', ' ', ' ',
     .                maxno, numh, listh, Dscf, Datm,
     .                maxno, numh, listh, Hmat, 
     .                Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .                Exc, Dxc, dipol, fa, stress, ierror )
          if (ierror .ne. 0) return

          do ia = 1, nua
            do ix = 1, 3
              avfa(ix,ia) = avfa(ix,ia) + fa(ix,ia)
            enddo
          enddo
          do iv = 1, 3
            do ix = 1, 3
              avstre(ix,iv) = avstre(ix,iv) + stress(ix,iv)
            enddo
          enddo
          do ix = 1, 3
            avdipo(ix) = avdipo(ix) + dipol(ix)
          enddo
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

        do ia = 1, nua
          do ix = 1, 3
            fa(ix,ia) = avfa(ix,ia)/anpt
          enddo
        enddo
        do iv = 1, 3
          do ix = 1, 3
            stress(ix,iv) = avstre(ix,iv)/anpt
          enddo
        enddo
        do ix = 1, 3
          dipol(ix) = avdipo(ix)/anpt
        enddo
        enaatm = avenaa/anpt
        enascf = avenas/anpt
        uatm   = avuatm/anpt
        uscf   = avuscf/anpt
        duscf  = avdusc/anpt
        exc    = avexc/anpt
        dxc    = avdxc/anpt

      endif

      return
      end

