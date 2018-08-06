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
  integer, parameter :: nr=1024         ! number of radial points for basis orbitals
  integer, parameter :: maxlines = 100  ! max number of unfolded band lines

  ! Internal variables
  integer          :: i1, i2, i3, i123(3), ia, iao, ierr, ik, ikx(3), &
                      ikq(8), iline, io, iostat, iq, iqx(3), ir, irq, &
                      isp, ispin, iispin, iu, iw, &
                      j, je(0:1), jk, jlm, jspin(maxlines), jw, &
                      kdsc(3,3), kscell(3,3), l, lmax, &
                      m, maux(3,3,2), maxorb, maxw, mk, &
                      ml(3,3), mr(3,3), mw, &
                      na, ne(maxlines), nk, nktot, nkx(3), nlines, &
                      nlm, no, nq(maxlines), nspin, nw, proj(3,3), t
  real(dp)         :: c0, dkcell(3,3), de, deline, dq, dqline(3), dqx(3), &
                      dr, drq,  dscell(3,3), emax(maxlines), emin(maxlines), &
                      g(3,8), grad, gradv(3), k0(3), kcell(3,3), kmax, kq(3,8), &
                      phi, pi, q(3), qmod, qcell(3,3), qmax(3,maxlines), qx(3), &
                      r, rc, scell(3,3), vol, wkq(8), wq
  complex(dp)      :: ii, psik, psiq
  logical          :: ccq(8), found, gamma
  character(len=50):: eunit, fname, formatstr,  slabel
  type(block_fdf)  :: bfdf

  ! Allocatable arrays and pointers
  integer,          allocatable:: cnfigfio(:), iaorb(:), ie(:,:,:,:), &
                                  indk(:,:,:), iphorb(:)
  real(sp),         allocatable:: psi(:,:,:,:,:)
  real(dp),         allocatable:: dos(:,:,:), dossum(:,:), ek(:,:,:), &
                                  k(:,:), gylm(:,:), &
                                  phir(:,:,:), phiq(:,:,:), &
                                  we(:,:,:,:), wk(:), ylm(:,:)
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

  ! Find atomic orbitals in real space
  print*,'unfold: reading atomic orbitals'
  lmax = 0
  do isp = 1,nsp            ! species index
    do io = 1,nofis(isp)    ! orbital index within species
!      write(6,*) "# Orbital (#, l, z, m, rc):", &
!        lofio(isp,io), zetafio(isp,io), mofio(isp,io), rcut(isp,io)
      lmax = max(lmax,lofio(isp,io))
      rc = rcut(isp,io)
      dr = rc/nr
      do ir=0,nr
        r = dr*ir
        call rphiatm(isp,io,r,phir(ir,isp,io),grad)
      enddo
    enddo ! io
  enddo ! isp
!  print*, 'lmax=',lmax

  ! Fourier transform atomic orbitals
  print*,'unfold: Fourier-transforming atomic orbitals'
  do isp = 1,nsp
    do io = 1,nofis(isp)
      call radfft(lofio(isp,io),nr,rcut(isp,io),phir(:,isp,io),phiq(:,isp,io))
    enddo
  enddo

  ! Find unit cell and initialize atomic coords
  print'(a,/,(3f12.6))','unfold: reading system geometry'
  call coor(na,ucell)
  vol = volcel(ucell)
!  print'(a,/,(3f12.6))','unfold: unit_cell (bohr) =',ucell
!  print'(a,/,(2i6,3f12.6))','unfold: ia,ispecies,xa (bohr) =', &
!    (ia,isa(ia),xa(:,ia),ia=1,na)

  ! Find k-points in FBZ
  print*,'unfold: generating k-points'
  kscell = 0
  kmax = 0
  call kgridinit( ucell, kscell, k0, kmax, nk )  ! read kscell from fdf file
  allocate( k(3,nk), wk(nk) )
  call kgrid( ucell, kscell, k0, nk, k, wk )     ! generate k vectors

  ! Find cell vectors of k grid (see kgrid.F)
  pi = acos(-1._dp)
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
!  print'(a,i6,/,(3f12.6))','unfold: nktot, dkcell =',nktot,dkcell
!  print'(a,/,(3f12.6))','unfold: dkcell =',dkcell
!  print'(a,/,(i6,3f12.6,3x,3f9.3))','unfolding: k, indk =', &
!    (ik,k(:,ik),matmul((k(:,ik)-k0),dscell)/(2*pi),ik=1,nk)

  ! Find index of all vectors in FBZ
  forall(j=1:3) nkx(j)=kdsc(j,j)
  allocate( indk(nkx(1),nkx(2),nkx(3)), cc(nkx(1),nkx(2),nkx(3)) )
  do ik = 1,nk
    ikx = nint(matmul(k(:,ik)-k0,dscell)/(2*pi))
    ikx = modulo(ikx,nkx)+1
    indk(ikx(1),ikx(2),ikx(3)) = ik
    cc(ikx(1),ikx(2),ikx(3)) = .false.
    ikx = nint(matmul(-k(:,ik)-k0,dscell)/(2*pi))
    ikx = modulo(ikx,nkx)+1
    indk(ikx(1),ikx(2),ikx(3)) = ik
    cc(ikx(1),ikx(2),ikx(3)) = .true.
  enddo
  if (any(indk==0)) &
    call die('unfold ERROR: some k-points in FBZ not found')
!  print'(a,/,(3i4,i6,l6))','unfold: indk =', &
!    (((i1,i2,i3,indk(i1,i2,i3),cc(i1,i2,i3),i1=1,nkx(1)),i2=1,nkx(2)),i3=1,nkx(3))

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
      read(iu) jk, k(:,ik)  ! k-point coordinates
      if (jk/=ik) &
        call die('unfold ERROR: unexpected k-point index')
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
      enddo
    enddo
  enddo
  close (iu)

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
  nlines = 0
  do while( fdf_bline(bfdf,pline) )
    if ( fdf_bmatch(pline,'iivvvvvs') .or. fdf_bmatch(pline,'iivvvvvsv') ) then
      nlines = nlines+1
      if (nlines>maxlines) &
        call die('unfold ERROR: parameter maxlines too small')
      nq(nlines) = fdf_bintegers(pline,1)
      ne(nlines) = fdf_bintegers(pline,2)
      qmax(:,nlines) = qcell(:,1)*fdf_bvalues(pline,1,after=2) &
                     + qcell(:,2)*fdf_bvalues(pline,2,after=2) &
                     + qcell(:,3)*fdf_bvalues(pline,3,after=2)
      eunit        = fdf_bnames(pline,1)          ! energy unit
      emin(nlines) = fdf_bvalues(pline,1,after=5)*fdf_convfac(eunit,'eV')
      emax(nlines) = fdf_bvalues(pline,2,after=5)*fdf_convfac(eunit,'eV')
      if (fdf_bnvalues(pline,after=7)>0) then
        jspin(nlines) = fdf_bvalues(pline,1,after=7)
      else
        jspin(nlines) = 0
      endif
!      print*,'iline,nq,ne,qmax,emin,emax=', &
!        nlines,nq(nlines),ne(nlines),qmax(:,nlines),emin(nlines),emax(nlines)
    else
      call die('unfold ERROR: wrong format in fdf block UnfoldedBandLines')
    endif
  enddo

  ! Spherical harmonics for q vectors
  print*,'unfold: finding spherical harmonics'
  nlm = (lmax+1)**2
  allocate(ylm(nlm,nk),gylm(3,nlm))
  do iline = 1,nlines
    qmod = sqrt(sum(qmax(:,iline)**2))
    call rlylm( lmax, qmax(:,iline)/qmod, ylm(:,iline), gylm )
  enddo

  ! Loop on unfolded band lines
  allocate( ie(0:1,nw,nk,nspin), we(0:1,nw,nk,nspin) )
  ii = cmplx(0,1)
  c0 = (2*pi)**1.5_dp / vol
  do iline = 1,nlines
    allocate( dos(0:nq(iline),0:ne(iline),nspin) )
    dos = 0
    iu = 12
    if (iline<10) then
      formatstr = '(A16,I1,A4)'
    else
      formatstr = '(A16,I2,A4)'
    endif
    write(fname,formatstr) 'unfoldedBandLine',iline,'.out'
    print*, fname
    open(iu,file=fname,status='unknown',form='formatted',action='write')
!    write(iu,*) nq(iline), ne(iline), jspin(iline)
    dqline(:) = qmax(:,iline)/nq(iline)
    deline = (emax(iline)-emin(iline))/ne(iline)
    do iw = 1,nw
      do ispin = 1,nspin
        do ik = 1,nk
          ie(0,iw,ik,ispin) = floor((ek(iw,ik,ispin)-emin(iline))/deline)
          ie(1,iw,ik,ispin) = ie(0,iw,ik,ispin)+1
          de = emin(iline) + ie(1,iw,ik,ispin)*deline - ek(iw,ik,ispin)
          we(0,iw,ik,ispin) = de/deline
          we(1,iw,ik,ispin) = 1-we(0,iw,ik,ispin)
        enddo
      enddo
      do iq = 0,nq(iline)
        q = qmax(:,iline)*iq/nq(iline)
        qx = matmul(q-k0,dscell)/(2*pi)
        iqx = floor(qx)
        dqx = iqx+1-qx
        j = 0
        do i3 = 0,1
        do i2 = 0,1
        do i1 = 0,1
          j = j+1
          i123 = (/i1,i2,i3/)
          ikx = modulo(iqx+i123,nkx)+1
          ikq(j) = indk(ikx(1),ikx(2),ikx(3))
          wkq(j) = product(i123-(2*i123-1)*dqx(:))
          ccq(j) = cc(ikx(1),ikx(2),ikx(3))
          kq(:,j) = matmul(dkcell,real(ikx-1,dp))
          g(:,j) = q(:) - kq(:,j)   
        enddo
        enddo
        enddo
        do io = 1,no
          ia = iaorb(io)
          isp = isa(ia)
          iao = iphorb(io)
          l = lofio(isp,iao)
          m = mofio(isp,iao)
          jlm = l*(l+1) + m+1                 !! ilm(l,m)
          qmod = sqrt(sum(q**2))
          dq = pi/rcut(isp,iao)
          irq = floor(qmod/dq)      ! irq, qmod, dq
          drq = (irq+1)*dq-qmod
          wq = drq/dq                ! drq, dq, wq
          phi = phiq(irq,isp,io)*wq + phiq(irq+1,isp,io)*(1-wq)
          phi = phi*ylm(jlm,iline)
          do j = 1,8
            do ispin = 1,nspin
              if (size(psi,1)==2) then
                psik = cmplx(psi(1,io,iw,ikq(j),ispin), &
                             psi(2,io,iw,ikq(j),ispin))
                if (ccq(j)) psik = conjg(psik)
              else
                psik = psi(1,io,iw,ikq(j),ispin)
              endif
              psiq = c0*phi*psik*exp(-ii*sum(g(:,j)*xa(:,ia)))
              je(:) = ie(:,iw,ikq(j),ispin)
              if (je(0)>=0 .and. je(1)<=ne(iline)) then
                dos(iq,je(:),ispin) = dos(iq,je(:),ispin) &
                            + we(:,iw,ikq(j),ispin)*wkq(j)*abs(psiq)**2
!                print'(a,6e15.6)','unfold: we,wkq,psi2,dos=', &
!                  we(:,iw,ikq(j),ispin),wkq(j),abs(psiq)**2,dos(iq,je(:),ispin)
              endif
            enddo ! ispin
          enddo ! j
        enddo ! io
      enddo ! iq
    enddo ! iw
    ispin = jspin(iline)
    do iq = 0,nq(iline)
      if (nspin==2) then
        if (ispin==0) then
          write(iu,*) dos(iq,:,1) + dos(iq,:,2)
        else
          write(iu,*) dos(iq,:,ispin)
        endif
      else
        if (ispin==0) then
          write(iu,*) 2*dos(iq,:,1)
        else
          write(iu,*) dos(iq,:,1)
        endif
      endif
    enddo ! iq
    close(iu)
    deallocate(dos)
  enddo ! iline

end program unfold

