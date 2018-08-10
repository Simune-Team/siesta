! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

program unfold

! Reads the .fdf, .ion, .psf and .WFSX files of a SIESTA calculation and generates
! unfolded bands.
! Ref: "Band unfolding made simple", S.Garcia-Mayo and J.M.Soler, Draft to be published
! S.Garcia-Mayo and J.M.Soler, Aug.2018

  use precision,    only: dp, sp
  use alloc,        only: re_alloc
  use chemical,     only: species_label
  use atm_types,    only: species_info, nspecies, species
  use atmfuncs,     only: lofio, mofio, nofis, rcut, rphiatm, zetafio
  use atom_options, only: get_atom_options
  use atom,         only: setup_atom_tables
  use basis_specs,  only: read_basis_specs
  use basis_types,  only: nsp, basis_specs_transfer, write_basis_specs
  use basis_io,     only: read_ion_ascii
  use sys,          only: die
  use fdf,          only: block_fdf, fdf_bintegers, fdf_bline, fdf_block, &
                          fdf_bmatch, fdf_bnames, fdf_bnvalues, fdf_bvalues, &
                          fdf_convfac, fdf_get, fdf_init, parsed_line
  use m_get_kpoints_scale, &
                    only: get_kpoints_scale
  use m_radfft,     only: radfft
  use cellsubs,     only: reclat, volcel
  use spher_harm,   only: lofilm, rlylm
  use siesta_geom,  only: ucell, xa, isa   ! unit cell, atomic coords and species

  implicit none

  ! Internal parameters
  integer, parameter :: nr = 4096       ! number of radial points for basis orbitals
  integer, parameter :: maxl = 5        ! max angular momentum
  integer, parameter :: maxlines = 100  ! max number of unfolded band lines
  integer, parameter :: maxig = 10      ! max index of refolding G vectors
  real(dp),parameter :: g2c = 300_dp    ! default mesh cutoff (ry)
  real(dp),parameter :: qc = 25_dp      ! cutoff for FT of atomic orbitals (bohr^-1)

  ! Internal variables
  integer          :: i, i1, i2, i3, i123(3), ia, iao, ierr, ig, ik, ikx(3), &
                      ikq(8), iline, io, iostat, iq, iqx(3), ir, irq, &
                      isp, ispin, iispin, iu, iw, &
                      j, je(0:1), jk, jlm, jspin(maxlines), jw, &
                      kdsc(3,3), kscell(3,3), l, lastq(0:maxlines), lmax, &
                      m, maux(3,3,2), maxorb, maxw, mk, &
                      ml(3,3), mr(3,3), mw, &
                      na, ne, ng, nk, nktot, nkx(3), nlines, &
                      nlm, no, nq, nqline, nrq, nspin, nw, proj(3,3), t
  real(dp)         :: alat, c0, dkcell(3,3), de, dek, diqx(3), dq, dqline(3), &
                      dqx(3), dr, drq, dscell(3,3), emax, emin, &
                      gcut, gq(3,8), gnew(3), gnorm, grad, gylm(3,maxl*maxl), &
                      k0(3), kcell(3,3), kmax, knew(3), kq(3), pi, &
                      q0(3), qmod, qcell(3,3), qg(3), qline(3), qmax, qx(3), &
                      r, rc, rcell(3,3), refoldCell(3,3), refoldBcell(3,3), rmax, rq, &
                      scell(3,3), vol, wkq(8), wq, xmax, ylm(maxl*maxl)
  complex(dp)      :: ck, ii, phi, psik, ukg(8)
  logical          :: ccq(8), found, gamma
  character(len=50):: eunit, fname, formatstr,  iostr, isstr, numstr, slabel
  type(block_fdf)  :: bfdf

  integer :: nx
  real(dp):: dg, drmin, dx, g0(3), gmod, gvec(3), normphiq, normphir, x(3), wr, xmod

  ! Allocatable arrays and pointers
  integer,          allocatable:: cnfigfio(:), iaorb(:), ie(:,:,:,:), &
                                  indk(:,:,:), iphorb(:)
  real(sp),         allocatable:: psi(:,:,:,:,:)
  real(dp),         allocatable:: dos(:,:,:), dossum(:,:), ek(:,:,:), &
                                  k(:,:), &
                                  phir(:,:,:), phiq(:,:,:), &
                                  we(:,:,:,:), wk(:)
  real(dp),         pointer    :: g(:,:)=>null(), q(:,:)=>null()
  logical,          allocatable:: cc(:,:,:)
  character(len=20),allocatable:: symfio(:), labelfis(:)
  type(species_info),  pointer :: spp
  type(parsed_line),   pointer :: pline

!--------------------

  ! Read atomic species
  call fdf_init('stdin','unfold_fdf.log')
  call get_atom_options()
  call read_xc_info()
  call read_basis_specs()
  call basis_specs_transfer()
  call setup_atom_tables(nsp)

  ! Allocate arrays for atomic orbitals
  nspecies = nsp
  allocate(species(nsp))
  maxorb = 0               ! max. number of orbitals in any atom
  do isp=1,nsp
    spp => species(isp)
    spp%label = species_label(isp)
    spp%read_from_file = .true.
    call read_ion_ascii(spp)
    maxorb = max(maxorb,nofis(isp))
  enddo
  allocate(phir(0:nr,nsp,maxorb))
  allocate(phiq(0:nr,nsp,maxorb))

  ! Find cutoffs of atomic orbitals
  lmax = 0
  rmax = 0
  do isp = 1,nsp            ! species index
    do io = 1,nofis(isp)    ! orbital index within species
      rmax = max(rmax,rcut(isp,io))
      lmax = max(lmax,lofio(isp,io))
    enddo ! io
  enddo ! isp
  print*, 'lmax,rmax=',lmax,rmax

  ! Find atomic orbitals in real space
  print*,'unfold: reading atomic orbitals'
  pi = acos(-1._dp)
  nq = nr
  dq = qc/nr                ! radial interval in reciprocal space
  dr = pi/qc                ! radial interval in real space
  rc = nr*dr                ! radial cutoff in real space
  do isp = 1,nsp            ! species index
    do io = 1,nofis(isp)    ! orbital index within species
!      write(6,*) "# Orbital (#, l, z, m, rc):", &
!        lofio(isp,io), zetafio(isp,io), mofio(isp,io), rcut(isp,io)
      do ir=0,nr
        r = dr*ir
        call rphiatm(isp,io,r,phir(ir,isp,io),grad)
      enddo
    enddo ! io
  enddo ! isp
  print*,'unfold: nr,dr,rc=',nr,dr,rc
  print*,'unfold: nq,dq,qc=',nq,dq,qc
  print'(a,/,(2i4,f12.6))','unfold: isp,io,rc=', &
    ((isp,io,rcut(isp,io),io=1,nofis(isp)),isp=1,nsp)

  ! Fourier transform atomic orbitals
  print*,'unfold: Fourier-transforming atomic orbitals'
  do isp = 1,nsp
    do io = 1,nofis(isp)
      call radfft(lofio(isp,io),nr,rc,phir(:,isp,io),phiq(:,isp,io))
    enddo
  enddo

  ! Write radial dependence of atomic orbitals in real and reciprocal space
  do isp = 1,nsp
    do io = 1,nofis(isp)
      write(isstr,*) isp
      write(iostr,*) io
      fname = 'species'//trim(adjustl(isstr))//'orbital'//trim(adjustl(iostr))//'.r'
      open(21,file=fname,status='unknown',form='formatted',action='write')
      write(21,'(2f12.6)') (ir*dr,phir(ir,isp,io),ir=0,nr)
      close(21)
      fname = 'species'//trim(adjustl(isstr))//'orbital'//trim(adjustl(iostr))//'.q'
      open(21,file=fname,status='unknown',form='formatted',action='write')
      close(21)
    enddo
  enddo

  ! Find norm of atomic orbitals in real space, as a check
!  gcut = 0.5*sqrt( fdf_get('MeshCutoff',g2c,'ry') )
!  xmax = rmax+1
!  dx = pi/gcut
!  nx = ceiling(xmax/dx)
!  print*,'unfold: dx,nx=',dx,nx
!  do isp = 1,nsp            ! species index
!    do io = 1,nofis(isp)    ! orbital index within species
!      l = lofio(isp,io)
!      m = mofio(isp,io)
!      jlm = l*(l+1) + m+1                 !! ilm(l,m)
!      normphir = 0
!      do i3 = -nx,nx
!      do i2 = -nx,nx
!      do i1 = -nx,nx
!        x = (/i1,i2,i3/)*dx
!        xmod = sqrt(sum(x**2))+1.e-15
!        ir = ceiling(xmod/dr)
!        wr = (ir*dr-xmod)/dr
!        if (ir<=nr) then
!          phi = phir(ir-1,isp,io)*(1-wr) + phir(ir,isp,io)*wr
!          call rlylm(l,x/xmod,ylm,gylm)
!          phi = phi*ylm(jlm)
!          normphir = normphir + abs(phi)**2*dx**3
!        endif
!      enddo
!      enddo
!      enddo
!      print*,'ulfold: isp,io,norm(phi(r))=',isp,io,normphir
!    enddo
!  enddo

  ! Find norm of atomic orbitals in reciprocal space, as a check
  dg = 2*pi/(2*xmax)
!  ng = ceiling(gcut/dg)
!  print*,'unfold: dg,ng=',dg,ng
!  do isp = 1,nsp            ! species index
!    do io = 1,nofis(isp)    ! orbital index within species
!      l = lofio(isp,io)
!      m = mofio(isp,io)
!      jlm = l*(l+1) + m+1                 !! ilm(l,m)
!      normphiq = 0
!      do i3 = -ng,ng
!      do i2 = -ng,ng
!      do i1 = -ng,ng
!        gvec = (/i1,i2,i3/)*dg
!        gmod = sqrt(sum(gvec**2))+1.e-15
!        iq = ceiling(gmod/dq)
!        wq = (iq*dq-gmod)/dq
!        if (iq<=nq) then
!          phi = phiq(iq-1,isp,io)*(1-wq) + phiq(iq,isp,io)*wq
!          call rlylm(l,gvec/gmod,ylm,gylm)
!          phi = phi*ylm(jlm)
!          normphiq = normphiq + abs(phi)**2*dg**3
!        endif
!      enddo
!      enddo
!      enddo
!      print*,'ulfold: isp,io,norm(phi(q))=',isp,io,normphiq
!    enddo
!  enddo

  ! Find unit cell and initialize atomic coords
  print'(a,/,(3f12.6))','unfold: reading system geometry'
  alat = fdf_get('LatticeConstant',1._dp,'bohr')
  call coor(na,ucell)        ! atomic coordinates xa stored in module siesta_geom
  vol = volcel(ucell)        ! unit cell volume
  call reclat(ucell,rcell,1) ! rcell = reciprocal cell vectors
  print'(a,f12.6,/,(3f21.15))','unfold: vol, unit_cell (bohr) =',vol,ucell
  print'(a,/,(2i6,3f12.6))','unfold: ia,ispecies,xa (bohr) =', &
    (ia,isa(ia),xa(:,ia),ia=1,na)
  print'(a,/,(3f21.15))','unfold: reciprocal_cell (1/bohr) =',rcell

  ! Find k-points in FBZ
  print*,'unfold: generating k-points'
  kscell = 0
  kmax = 0
  call kgridinit( ucell, kscell, k0, kmax, nk )  ! read kscell from fdf file
  print*,'unfold: from kgridinit nk=',nk
  allocate( k(3,nk), wk(nk) )
  call kgrid( ucell, kscell, k0, nk, k, wk )     ! generate k vectors

  ! Find cell vectors of k grid (see kgrid.F)
  call idiag( 3, kscell, kdsc, ml, mr, maux )    ! kdsc = diagonal supercell
  proj = 0
  forall(j=1:3) proj(j,j)=sign(1,kdsc(j,j))
  kdsc = matmul(kdsc,proj)
  mr = matmul(mr,proj)
  dscell = matmul(ucell,real(kscell,dp))         ! reciprocal k-grid vectors
  dscell = matmul(dscell,mr)                     ! same but 'diagonal'
  call reclat(dscell,dkcell,1)                   ! dkcell = k-grid unit vectors
  call reclat(ucell,kcell,1)                     ! kcell = reciprocal lattice vecs
  k0 = matmul(dkcell,k0)                         ! k-grid origin
  nktot = nint(volcel(kcell)/volcel(dkcell))     ! total num. of k vecs in FBZ
  print'(a,i6,/,(3f12.6))','unfold: nk, dkcell =',nk,dkcell
!  print'(a,/,(i6,3f12.6,3x,3f9.3))','unfolding: k, indk =', &
!    (ik,k(:,ik),matmul((k(:,ik)-k0),dscell)/(2*pi),ik=1,nk)

  ! Find index of all vectors in FBZ
  print*,'unfold: indexing k vectors, nk=',nk
  forall(j=1:3) nkx(j)=kdsc(j,j)
  allocate( indk(nkx(1),nkx(2),nkx(3)), cc(nkx(1),nkx(2),nkx(3)) )
  indk = 0
  do ik = 1,nk
    ikx = nint(matmul(k(:,ik)-k0,dscell)/(2*pi))
    ikx = modulo(ikx,nkx)+1
    indk(ikx(1),ikx(2),ikx(3)) = ik
    cc(ikx(1),ikx(2),ikx(3)) = .false.
    ikx = nint(matmul(-k(:,ik)-k0,dscell)/(2*pi))
    ikx = modulo(ikx,nkx)+1
    if (indk(ikx(1),ikx(2),ikx(3))==0) then
      indk(ikx(1),ikx(2),ikx(3)) = ik
      cc(ikx(1),ikx(2),ikx(3)) = .true.
    endif
  enddo
  if (any(indk==0)) &
    call die('unfold ERROR: some k-points in FBZ not found')

!  print'(a,/,(3i4,i6,6f12.6))','unfold: ikx,indk,k,kfrac =', &
!    (((i1,i2,i3,indk(i1,i2,i3),k(:,indk(i1,i2,i3)), &
!    matmul(k(:,indk(i1,i2,i3)),ucell)/(2*pi),i1=1,nkx(1)),i2=1,nkx(2)),i3=1,nkx(3))

  ! Read band energies and wavefunctions at k points of the FBZ
  print*,'unfold: reading wavefunctions'
  slabel = fdf_get('SystemLabel','unknown')
  fname = trim(slabel)//'.fullBZ.WFSX'
  iu = 11
  open(iu, file=fname, form='unformatted', status='old' )
  read(iu) mk, gamma    ! number of k-points in FBZ, gamma-point only?
  if (mk/=nk) &
    call die('unfold ERROR: unexpected nk in .WFSX file')
  read(iu) nspin        ! number of spin components
  read(iu) no           ! number of atomic orbitals in unit cell
  maxw = no             ! max number of bands
  if (gamma) then
    allocate(psi(1,no,maxw,nk,nspin))
  else
    allocate(psi(2,no,maxw,nk,nspin))
  endif
  allocate( ek(maxw,nk,nspin) )
  allocate( iaorb(no), labelfis(no), iphorb(no), cnfigfio(no), symfio(no) )
  read(iu) (iaorb(io),labelfis(io),iphorb(io),cnfigfio(io),symfio(io),io=1,no)
  do ik = 1,nk
    do ispin = 1,nspin
      read(iu) jk, knew  ! k-point coordinates
      if (jk/=ik) &
        call die('unfold ERROR: unexpected k-point index')
      if (sqrt(sum((knew-k(:,ik))**2))>1.e-12_dp) &
        call die('unfold ERROR: unexpected k vector')
      read(iu) iispin                
      if (iispin/=ispin) &
        call die('unfold ERROR: unexpected spin index')
      read(iu) nw           ! number of wavefunctions (bands)
      if (ik==1 .and. ispin==1) then
        mw = nw
      else if (nw/=mw) then
          call die('unfold ERROR: number of bands not constant')
      endif
      do iw = 1,nw
        read(iu) jw
        if (jw/=iw) &
          call die('unfold ERROR: unexpected band index')
        read(iu) ek(iw,ik,ispin)                  ! band energy
        read(iu) (psi(:,io,iw,ik,ispin),io=1,no)  ! wavefunction coeffs.
!        print'(a,2i4,15e12.3)','unfold: iw,ik,psi=',iw,ik,psi(1,:,iw,ik,1)
!        print'(a,2i4,e15.6)','unfold: iw,ik,norm(psi)=',iw,ik,sum(psi(:,:,iw,ik,1)**2)
      enddo
    enddo
  enddo
  close (iu)

!  print'(a,/,(i6,3f12.6,3x,3f9.3))','unfolding: k, indk =', &
!    (ik,k(:,ik),matmul((k(:,ik)-k0),dscell)/(2*pi),ik=1,nk)
!  print'(a,/,(2i6,e15.6))','unfold: ik,iw,ek =',((ik,iw,ek(iw,ik,1),iw=1,nw),ik=1,nk)
  print*,'unfold: ekmin =', minval(ek),' eV'

  ! Read fdf block UnfoldedBandLines
  print*,'unfold: reading UnfoldedBandLines block'
  call get_kpoints_scale('BandLinesScale',qcell,ierr)
  if (ierr/=0) &
    call die('unfold: ERROR calling get_kpoints_scale')
  found = fdf_block('UnfoldedBandLines',bfdf)
  if (.not.found) &
    call die('unfold ERROR: fdf block UnfoldedBandLines not found')
  if (.not.(fdf_bline(bfdf,pline).and.fdf_bmatch(pline,'ivvs'))) &
    call die('unfold ERROR: wrong format in fdf block UnfoldedBandLines')
  ne    = fdf_bintegers(pline,1)
  eunit = fdf_bnames(pline,1)          ! energy unit
  emin  = fdf_bvalues(pline,1,after=1)*fdf_convfac(eunit,'eV')
  emax  = fdf_bvalues(pline,2,after=1)*fdf_convfac(eunit,'eV')
  nlines = 0
  nq = 0
  lastq(0) = 0
  do while( fdf_bline(bfdf,pline) )
    if (fdf_bmatch(pline,'ivvv').or.fdf_bmatch(pline,'ivvvs')) then
      nqline = fdf_bintegers(pline,1)
!      print*,'unfold: nqline=',nqline
      qline = qcell(:,1)*fdf_bvalues(pline,1,after=1) &
            + qcell(:,2)*fdf_bvalues(pline,2,after=1) &
            + qcell(:,3)*fdf_bvalues(pline,3,after=1)
      if (nqline==1) then
        nlines = nlines+1
        if (nlines>maxlines) &
          call die('unfold ERROR: parameter maxlines too small')
        nq = nq+1
        call re_alloc(q,1,3,1,nq)
        q(:,nq) = qline
      else
        call re_alloc(q,1,3,1,nq+nqline)
        do iq = 1,nqline
          q(:,nq+iq) = q(:,nq) + (qline-q(:,nq))*iq/nqline
        enddo
        nq = nq+nqline
      endif 
      lastq(nlines) = nq
    else
      call die('unfold ERROR: wrong format in fdf block UnfoldedBandLines')
    endif
  enddo

! Read fdf block RefoldingLatticeVectors
  print*,'unfold: reading RefoldingLatticeVectors block'
  gcut = 0.5*sqrt( fdf_get('MeshCutoff',g2c,'ry') )
  if (fdf_block('RefoldingLatticeVectors',bfdf)) then
    print*,'unfold: block RefoldingLatticeVectors found'
    do j = 1,3
      if (.not.(fdf_bline(bfdf,pline).and.fdf_bmatch(pline,'vvv'))) &
        call die('unfold ERROR: wrong format in fdf block RefoldingLatticeVectors')
      do i = 1,3
        refoldCell(i,j)  = fdf_bvalues(pline,i)*alat
      enddo
    enddo
    call reclat(refoldCell,refoldBcell,1)
    ng = 0
    do i3 = -maxig,maxig
    do i2 = -maxig,maxig
    do i1 = -maxig,maxig
      gnew = refoldBcell(:,1)*i1 + refoldBcell(:,2)*i2 + refoldBcell(:,3)*i3
      gnorm = sqrt(sum(gnew**2))
      if (gnorm<gcut) then
        ng = ng+1
        call re_alloc(g,1,3,1,ng)
        g(:,ng) = gnew
      endif
    enddo
    enddo
    enddo
  else
    ng = 1
    call re_alloc(g,1,3,1,ng)
    g(:,1) = 0
  endif
  print'(a,/,(3f12.6))','unfold: refoldCell=',refoldCell
  print*,'unfold: alat,gcut,ng=',alat,gcut,ng

!  print'(a,/,(i6,3f12.6))','unfold: iq,q=',(iq,q(:,iq),iq=1,nq)
!  print*,'unfold: nlines,lastq=',nlines,lastq(0:nlines)

  ! Find index and weight of band energies
  allocate( ie(0:1,nw,nk,nspin), we(0:1,nw,nk,nspin) )
  de = (emax-emin)/ne
  do ispin = 1,nspin
    do ik = 1,nk
      do iw = 1,nw
        ie(0,iw,ik,ispin) = floor((ek(iw,ik,ispin)-emin)/de)
        ie(1,iw,ik,ispin) = ie(0,iw,ik,ispin)+1
        dek = emin + ie(1,iw,ik,ispin)*de - ek(iw,ik,ispin)
        we(0,iw,ik,ispin) = dek/de
        we(1,iw,ik,ispin) = 1-we(0,iw,ik,ispin)
      enddo
    enddo
  enddo

  ! Loop on unfolded band lines
  print*,'unfold: main loop'
  ii = cmplx(0,1)
  nlm = (lmax+1)**2
  c0 = (2*pi)**1.5_dp / vol
  do iline = 1,nlines
    allocate( dos(lastq(iline-1)+1:lastq(iline),0:ne,nspin) )
    dos = 0
    do iq = lastq(iline-1)+1,lastq(iline)
!      print*,'unfold: q=',q(:,iq)
      do ig = 1,ng
        qg = q(:,iq)+g(:,ig) ! note: this is a refolding g, not a superlattice one
        qx = matmul(qg-k0,dscell)/(2*pi)  ! qg in mesh coords: qg=k0+matmul(dkcell,qx)
        iqx = floor(qx-1.e-12)            ! k-mesh point in mesh coords
        diqx = iqx+1-qx
        j = 0
        do i3 = 0,1
        do i2 = 0,1
        do i1 = 0,1
          j = j+1
          i123 = (/i1,i2,i3/)
          ikx = modulo(iqx+i123,nkx)+1
          ikq(j) = indk(ikx(1),ikx(2),ikx(3))
          wkq(j) = product(i123-(2*i123-1)*diqx(:))
          ccq(j) = cc(ikx(1),ikx(2),ikx(3))
          kq = k(:,ikq(j))
          if (ccq(j)) kq = -kq
          gq(:,j) = qg-kq
        enddo
        enddo
        enddo
        qmod = sqrt(sum(qg**2))+1.e-15
        call rlylm( lmax, qg/qmod, ylm, gylm )
        do ispin = 1,nspin
          do iw = 1,nw
            ukg = 0
            do io = 1,no
              ia = iaorb(io)
              isp = isa(ia)
              iao = iphorb(io)
              l = lofio(isp,iao)
              m = mofio(isp,iao)
              jlm = l*(l+1) + m+1                 !! ilm(l,m)
              qmod = sqrt(sum(qg**2))
              irq = floor(qmod/dq)
              drq = (irq+1)*dq-qmod
              wq = drq/dq
              phi = phiq(irq,isp,iao)*wq + phiq(irq+1,isp,iao)*(1-wq)
              phi = (-ii)**l * ylm(jlm) * phi
              do j = 1,8
                if (gamma) then
                  ck = psi(1,io,iw,ikq(j),ispin)
                else
                  ck = cmplx(psi(1,io,iw,ikq(j),ispin), &
                             psi(2,io,iw,ikq(j),ispin))
                endif
                if (ccq(j)) ck = conjg(ck)
                psik = c0*ck*phi*exp(-ii*sum(gq(:,j)*xa(:,ia)))
                ukg(j) = ukg(j) + psik
              enddo ! j
            enddo ! io
            do j = 1,8
              je(:) = ie(:,iw,ikq(j),ispin)
              if (je(0)>=0 .and. je(1)<=ne) then
                dos(iq,je(:),ispin) = dos(iq,je(:),ispin) &
                              + vol*we(:,iw,ikq(j),ispin)*wkq(j)*abs(ukg(j))**2
              endif
            enddo ! j
          enddo ! iw
        enddo ! ispin
      enddo ! ig
    enddo ! iq

    do ispin = 1,nspin
      iu = 12
      write(numstr,*) iline
      fname = 'unfoldedBandLine'//adjustl(numstr)
      if (nspin==2 .and. ispin==1) then
        fname = trim(fname)//'spinUp.out'
      elseif (nspin==2 .and. ispin==2) then
        fname = trim(fname)//'spinDown.out'
      else
        fname = trim(fname)//'.out'
      endif
      open(iu,file=fname,status='unknown',form='formatted',action='write')
      write(iu,*) lastq(iline)-lastq(iline-1),ne+1,emin,emax
      do iq = lastq(iline-1)+1,lastq(iline)
        write(iu,*) q(:,iq)
        do j = 0,ne
          write(iu,'(f12.6)') dos(iq,j,ispin)
        enddo
      enddo ! iq
      close(iu)
    enddo
    deallocate(dos)
  enddo ! iline

end program unfold

