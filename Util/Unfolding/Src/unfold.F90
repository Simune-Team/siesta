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

  use alloc,        only: alloc_report, de_alloc, re_alloc
  use atmfuncs,     only: lofio, mofio, nofis, rcut, rphiatm, zetafio
  use basis_types,  only: nsp
  use basis_types,  only: basis_parameters, initialize
  use basis_io,     only: read_basis_ascii
  use cellsubs,     only: reclat, volcel
  use fdf,          only: block_fdf, fdf_bintegers, fdf_bline, fdf_block, &
                          fdf_bmatch, fdf_bnames, fdf_bnvalues, fdf_bvalues, &
                          fdf_convfac, fdf_get, fdf_init, fdf_parallel, parsed_line
  use hsx_m,        only: hsx_t, read_hsx_file
!  use m_diag,       only: diag_init
  use m_get_kpoints_scale, &
                    only: get_kpoints_scale
  use m_mpi_utils,  only: broadcast
#ifdef MPI
  use m_diag_option,only: ParallelOverK, diag_serial=>Serial
  use mpi_siesta,   only: MPI_Comm_Rank, MPI_Comm_Size, MPI_Comm_World, &
                          MPI_double_precision, MPI_Init, MPI_Finalize, &
                          MPI_Reduce, MPI_Sum
#endif
  use m_radfft,     only: radfft
  use m_timer,      only: timer_init, timer_report, timer_start, timer_stop
  use parallel,     only: Nodes, Node
  use precision,    only: dp, sp
  use spher_harm,   only: lofilm, rlylm
  use siesta_geom,  only: ucell, xa, isa   ! unit cell, atomic coords and species
  use sys,          only: die


  implicit none

  ! Internal parameters
  character(len=*),parameter:: myName = 'unfold '
  integer, parameter :: nr = 4096       ! number of radial points for basis orbitals
  integer, parameter :: maxl = 5        ! max angular momentum
  integer, parameter :: maxlines = 100  ! max number of unfolded band lines
  integer, parameter :: maxig = 10      ! max index of refolding G vectors
  real(dp),parameter :: g2c = 300_dp    ! default mesh cutoff (ry)
  real(dp),parameter :: qc = 50_dp      ! cutoff for FT of atomic orbitals (bohr^-1)
  logical, parameter :: writeOrbitals = .false.  ! write atomic orbital files?
  real(dp),parameter :: tolSuperCell = 1.e-6 ! tolerance for comparing cell and supercell
  integer, parameter :: allocReportLevelDefault = 2 ! default allocation report level

  ! Internal variables
  integer          :: i, i1, i2, i3, ia, iao, ib, ie, ierr, ig, ij, &
                      iline, io, ios, iostat, iou, iq, iq1, iq2, iqNode, iqx(3), &
                      ir, irq, iscf, isp, ispin, iu, j, je, jk, jlm, jo, jos, jou, &
                      kdsc(3,3), kscell(3,3), l, lastq(0:maxlines), level, ll, lmax, &
                      m, maxorb, myNode, na, nbands, ne, ng, nh, nlines, nlm, nNodes, &
                      nos, nou, nq, nqline, nrq, nspin, ntmp, nw, t, z
  real(dp)         :: alat, c0, cellRatio(3,3), ddos, de, dek, dq, dqline(3), &
                      dqx(3), dr, drq, dscell(3,3), emax, emin, &
                      gcut, gnew(3), gnorm, gq(3), grad, gylm(3,maxl*maxl), &
                      kq(3), kxij, pi, &
                      qmod, qcell(3,3), qg(3), qline(3), qmax, qx(3), &
                      r, rc, rcell(3,3), refoldCell(3,3), refoldBcell(3,3), rmax, rq, &
                      scell(3,3), threshold, vol, we, wq, xmax, ylm(maxl*maxl)
  complex(dp)      :: ck, ii, phase, phi, psik, ukg
  logical          :: found, gamma, notSuperCell
  character(len=50):: eunit, fname, formatstr,  iostr, isstr, numstr, slabel
  character(len=20):: labelfis, symfio
  character(len=200):: line
  type(block_fdf)  :: bfdf
  type(hsx_t)      :: hsx

  integer :: nx
  real(dp):: dg, drmin, dx, g0(3), gmod, gvec(3), normphiq, normphir, x(3), wr, xmod

  ! Allocatable arrays and pointers
  real(dp),         pointer:: dos(:,:,:)=>null(), dossum(:,:)=>null(), eb(:,:)=>null(), &
                              phir(:,:,:)=>null(), phiq(:,:,:)=>null(), &
                              tmp1(:)=>null(), tmp2(:)=>null()
  real(dp),         pointer:: g(:,:)=>null(), q(:,:)=>null()
  complex(dp),      pointer:: h(:,:)=>null(), psi(:,:,:)=>null(), s(:,:)=>null()
  logical,          pointer:: cc(:,:,:)=>null()
  type(parsed_line),pointer:: pline=>null()

!--------------------

! Initialize MPI
#ifdef MPI
  call MPI_Init( ierr )
  call MPI_Comm_Rank( MPI_Comm_World, Node, ierr )
  call MPI_Comm_Size( MPI_Comm_World, Nodes, ierr )
  myNode = Node
  nNodes = Nodes
#else
  myNode = 0
  nNodes = 1
#endif

  ! Initialize input
#ifdef MPI
  if (.not.fdf_parallel()) &
    call die('unfold ERROR: FDF has no parallel support')
#endif
  call fdf_init( fileOutput='unfold.fdflog', unitInput=5 )

  ! Initialize timer
  threshold = fdf_get('TimerReportThreshold', 0._dp)
  call timer_report( file='unfold.times', threshold=threshold )
  call timer_init()
  call timer_start(myName//'init')

  ! Set allocation report parameters
  level = fdf_get( 'AllocReportLevel', AllocReportLevelDefault )
  threshold = fdf_get( 'AllocReportThreshold', 0._dp )
  call alloc_report( level, file='unfold.alloc', printNow=.false., threshold=threshold )

  ! Read atomic basis
  if (myNode==0) then
    print*,'unfold: reading atomic orbitals'
    allocate(basis_parameters(nsp))
    do isp=1,nsp
      call initialize(basis_parameters(isp))
    enddo
    call read_basis_ascii(nsp)
  endif
  call broadcast(nsp)
  call broadcast_basis()

  ! Allocate arrays for atomic orbitals
  maxorb = 0               ! max. number of orbitals in any atom
  do isp=1,nsp
    maxorb = max(maxorb,nofis(isp))
  enddo
  call re_alloc( phir, 0,nr, 1,nsp, 1,maxorb, myName//'phir' )
  call re_alloc( phiq, 0,nr, 1,nsp, 1,maxorb, myName//'phiq' )

  ! Find cutoffs of atomic orbitals
  lmax = 0
  rmax = 0
  do isp = 1,nsp            ! species index
    do io = 1,nofis(isp)      ! orbital index within species
      rmax = max(rmax,rcut(isp,io))
      lmax = max(lmax,lofio(isp,io))
    enddo ! io
  enddo ! isp
  if (myNode==0) print*, 'lmax,rmax=',lmax,rmax

  ! Find atomic orbitals in real space
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
  if (myNode==0) then
    print*,'unfold: nr,dr,rc=',nr,dr,rc
    print*,'unfold: nq,dq,qc=',nq,dq,qc
    print'(a,/,(2i4,f12.6))','unfold: isp,io,rc=', &
      ((isp,io,rcut(isp,io),io=1,nofis(isp)),isp=1,nsp)
  endif

  ! Fourier transform atomic orbitals
  if (myNode==0) print*,'unfold: Fourier-transforming atomic orbitals'
  do isp = 1,nsp
    do io = 1,nofis(isp)
      call radfft(lofio(isp,io),nr,rc,phir(:,isp,io),phiq(:,isp,io))
    enddo
  enddo

  ! Write radial dependence of atomic orbitals in real and reciprocal space
  if (myNode==0 .and. writeOrbitals) then
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
  endif

  ! Find norm of atomic orbitals in real space, as a check
!  gcut = 0.5*sqrt( fdf_get('MeshCutoff',g2c,'ry') )
!  xmax = rmax+1
!  dx = pi/gcut
!  nx = ceiling(xmax/dx)
!  if (myNode==0) print*,'unfold: dx,nx=',dx,nx
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
!      if (myNode==0) print*,'ulfold: isp,io,norm(phi(r))=',isp,io,normphir
!    enddo
!  enddo

  ! Find norm of atomic orbitals in reciprocal space, as a check
!  dg = 2*pi/(2*xmax)
!  ng = ceiling(gcut/dg)
!  if (myNode==0) print*,'unfold: dg,ng=',dg,ng
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
!      if (myNode==0) print*,'ulfold: isp,io,norm(phi(q))=',isp,io,normphiq
!    enddo
!  enddo

  ! Find unit cell and initialize atomic coords
  if (myNode==0) print'(a,/,(3f12.6))','unfold: reading system geometry'
  alat = fdf_get('LatticeConstant',1._dp,'bohr')
  call coor(na,ucell)        ! atomic coordinates xa stored in module siesta_geom
  vol = volcel(ucell)        ! unit cell volume
  call reclat(ucell,rcell,1) ! rcell = reciprocal cell vectors
  if (myNode==0) then
    print'(a,f12.6,/,(3f21.15))','unfold: vol, unit_cell (bohr) =',vol,ucell
    print'(a,/,(2i6,3f12.6))','unfold: ia,ispecies,xa (bohr) =', &
      (ia,isa(ia),xa(:,ia),ia=1,na)
    print'(a,/,(3f21.15))','unfold: reciprocal_cell (1/bohr) =',rcell
  endif

  ! Read HSX file with H, S, and xij matrices
  slabel = fdf_get('SystemLabel','unknown')
  fname = trim(slabel)//'.HSX'
  if (myNode==0) then
    call read_hsx_file(hsx,fname)
    if (hsx%nspecies/=nsp) call die('unfold ERROR: hsx%nspecies/=nsp')
    if (hsx%na_u/=na) call die('unfold ERROR: hsx%na_u/=na')
  endif
  call broadcast(hsx%nspecies)
  call broadcast(hsx%na_u)
  call broadcast(hsx%no_u)
  call broadcast(hsx%no_s)
  call broadcast(hsx%nspin)
  call broadcast(hsx%nh)
  nou = hsx%no_u
  nos = hsx%no_s
  nh  = hsx%nh
  nspin = hsx%nspin
  if (myNode/=0) then
    call re_alloc(hsx%iaorb,    1,nou,          myName//'hsx%iaorb')
    call re_alloc(hsx%iphorb,   1,nou,          myName//'hsx%iphorb')
    call re_alloc(hsx%numh,     1,nou,          myName//'hsx%numh')
    call re_alloc(hsx%listhptr, 1,nou,          myName//'hsx%listhptr')
    call re_alloc(hsx%listh,    1,nh,           myName//'hsx%listh')
    call re_alloc(hsx%indxuo,   1,nos,          myName//'hsx%indxuo')
    call re_alloc(hsx%hamilt,   1,nh,  1,nspin, myName//'hsx%hamilt')
    call re_alloc(hsx%Sover,    1,nh,           myName//'hsx%Sover')
    call re_alloc(hsx%xij, 1,3, 1,nh,           myName//'hsx%xij')
  endif
  call broadcast(hsx%iaorb)
  call broadcast(hsx%iphorb)
  call broadcast(hsx%numh)
  call broadcast(hsx%listhptr)
  call broadcast(hsx%listh)
  call broadcast(hsx%indxuo)
  call broadcast(hsx%hamilt)
  call broadcast(hsx%Sover)
  call broadcast(hsx%xij)

  ! Read fdf block UnfoldedBandLines
  if (myNode==0) print*,'unfold: reading UnfoldedBandLines block'
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
!      if (myNode==0) print*,'unfold: nqline=',nqline
      qline = qcell(:,1)*fdf_bvalues(pline,1,after=1) &
            + qcell(:,2)*fdf_bvalues(pline,2,after=1) &
            + qcell(:,3)*fdf_bvalues(pline,3,after=1)
      if (nqline==1) then
        nlines = nlines+1
        if (nlines>maxlines) &
          call die('unfold ERROR: parameter maxlines too small')
        nq = nq+1
        call re_alloc( q, 1,3, 1,nq, myName//'q' )
        q(:,nq) = qline
      else
        call re_alloc( q, 1,3, 1,nq+nqline, myName//'q' )
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

!  print'(a,/,(i6,3f12.6))','unfold: iq,q=',(iq,q(:,iq),iq=1,nq)
!  print*,'unfold: nlines,lastq=',nlines,lastq(0:nlines)

  ! Read fdf block RefoldingLatticeVectors
  if (myNode==0) print*,'unfold: reading RefoldingLatticeVectors block'
  gcut = 0.5*sqrt( fdf_get('MeshCutoff',g2c,'ry') )
  if (fdf_block('RefoldingLatticeVectors',bfdf)) then
    if (myNode==0) print*,'unfold: block RefoldingLatticeVectors found'
    do j = 1,3
      if (.not.(fdf_bline(bfdf,pline).and.fdf_bmatch(pline,'vvv'))) &
        call die('unfold ERROR: wrong format in fdf block RefoldingLatticeVectors')
      do i = 1,3
        refoldCell(i,j) = fdf_bvalues(pline,i)*alat
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
        call re_alloc( g, 1,3, 1,ng, myName//'g', copy=.true. )
        g(:,ng) = gnew
      endif
    enddo
    enddo
    enddo
  else
    refoldCell = 0
    refoldBcell = 0
    ng = 1
    call re_alloc( g, 1,3, 1,ng, myName//'g' )
    g(:,1) = 0
  endif
  if (myNode==0) print'(a,/,(3f12.6))','unfold: refoldCell=',refoldCell
  if (myNode==0) print*,'unfold: alat,gcut,ng=',alat,gcut,ng

  ! Find if simulation cell is a supercell of refoldCell
  cellRatio = matmul(refoldBcell,ucell)/(2*pi)
  notSuperCell = .not.all(abs(cellRatio-nint(cellRatio))<tolSuperCell)
  if (myNode==0) print'(a,l,a,/,(3f12.6))', &
    'unfold: notSuperCell= ',notSuperCell,'   cellRatio=',cellRatio
  
  ! Set some constants
  ii = cmplx(0,1,dp)
  nlm = (lmax+1)**2
  c0 = (2*pi)**1.5_dp / vol
  de = (emax-emin)/ne
  ParallelOverK = .true.
  diag_serial = .true.
  nbands = nou
  iscf = 1
  call re_alloc( h,   1,nou, 1,nou,          myName//'h'   )
  call re_alloc( s,   1,nou, 1,nou,          myName//'s'   )
  call re_alloc( psi, 1,nou, 1,nou, 1,nspin, myName//'psi' ) 
  call re_alloc( eb,  1,nou,        1,nspin, myName//'eb'  )

  ! Main loops on unfolded band lines and q vectors along them
  call timer_stop(myName//'init')
  call timer_start(myName//'main loop')
  if (myNode==0) print*,'unfold: main loop'
  do iline = 1,nlines
    iq1 = lastq(iline-1)+1
    iq2 = lastq(iline)
    call re_alloc( dos, iq1,iq2, 0,ne, 1,nspin, myName//'dos', &
                   copy=.false., shrink=.true. )
    dos = 0
    do iq = iq1,iq2
      iqNode = mod((iq-1),nNodes)
      if (myNode==iqNode) then
!        print*,'unfold: q=',q(:,iq)
        do ig = 1,ng
          qg = q(:,iq)+g(:,ig) ! note: this is a refolding g, not a superlattice G
          if (ig==1 .or. notSuperCell) then
            qx = matmul(qg,ucell)/(2*pi)   ! qg in mesh coords: qg=matmul(rcell,qx)
            kq = matmul(rcell,qx-nint(qx)) ! q+g translated to FBZ of SC, i.e. kq=k+g-G
          endif
          gq = qg-kq                       ! supercell G vector 
          do ispin = 1,nspin

            if (ig==1 .or. notSuperCell) then
              call timer_start(myName//'diag')
              h = 0
              s = 0
              do iou = 1,nou
                do j = 1,hsx%numh(iou)
                  ij = hsx%listhptr(iou) + j
                  jos = hsx%listh(ij)
                  jou = hsx%indxuo(jos)
                  kxij = sum(kq*hsx%xij(:,ij))
                  phase = exp(ii*kxij)
                  s(iou,jou) = s(iou,jou) + hsx%Sover(ij)*phase
                  h(iou,jou) = h(iou,jou) + hsx%hamilt(ij,ispin)*phase
                enddo
              enddo
              call cdiag(h,s,nou,nou,nou,eb(:,ispin),psi(:,:,ispin), &
                         nbands,iscf,ierr,-1)
              if (ierr/=0) print*,'unfold: ERROR in cdiag'
              eb(:,ispin) = eb(:,ispin)*fdf_convfac('ry','ev') ! from Ry to eV
              call timer_stop(myName//'diag')
            endif

            call timer_start(myName//'g sum')
            qmod = sqrt(sum(qg**2))+1.e-15
            call rlylm( lmax, qg/qmod, ylm, gylm )
            do ib = 1,nbands
              je = floor((eb(ib,ispin)-emin)/de)
              dek = emin + (je+1)*de - eb(ib,ispin)
              we = dek/de
              ukg = 0
              do io = 1,nou
                ia = hsx%iaorb(io)
                isp = isa(ia)
                iao = hsx%iphorb(io)
                l = lofio(isp,iao)
                m = mofio(isp,iao)
                jlm = l*(l+1) + m+1                 !! ilm(l,m)
                irq = floor(qmod/dq)
                drq = (irq+1)*dq-qmod
                wq = drq/dq
                phi = phiq(irq,isp,iao)*wq + phiq(irq+1,isp,iao)*(1-wq)
                phi = (-ii)**l * ylm(jlm) * phi
                psik = c0*psi(io,ib,ispin)*phi*exp(-ii*sum(gq*xa(:,ia)))
                ukg = ukg + psik
              enddo ! io
              ddos = vol*abs(ukg)**2
              if (je>=0) dos(iq,je,ispin) = dos(iq,je,ispin) + ddos*we
              if (je<ne) dos(iq,je+1,ispin) = dos(iq,je+1,ispin) + ddos*(1-we)
            enddo ! ib
            call timer_stop(myName//'g sum')
          enddo ! ispin
        enddo ! ig
      endif ! (myNode==iqNode)
    enddo ! iq
#ifdef MPI
    ntmp = (iq2-iq1+1)*(ne+1)*nspin
    call re_alloc( tmp1, 1,ntmp, myName//'tmp1', copy=.false., shrink=.true. )
    call re_alloc( tmp2, 1,ntmp, myName//'tmp2', copy=.false., shrink=.true. )
    tmp1 = reshape(dos,(/ntmp/))
    tmp2 = 0
    call MPI_reduce(tmp1,tmp2,ntmp,MPI_double_precision,MPI_sum,0, &
                    MPI_COMM_WORLD,ierr)
    dos(iq1:iq2,0:ne,:) = reshape(tmp2,(/iq2-iq1+1,ne+1,nspin/))
#endif

    if (myNode==0) then
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
        do iq = iq1,iq2
          write(iu,*) q(:,iq)
          do j = 0,ne
            write(iu,'(f12.6)') dos(iq,j,ispin)
!            write(iu,'(3i4,f12.6)') iq,j,ispin,dos(iq,j,ispin)
          enddo
        enddo ! iq
        close(iu)
      enddo
    endif ! (myNode==0)
  enddo ! iline
  call timer_stop(myName//'main loop')

  ! Write allocation report
  call alloc_report( printNow=.true. )

  ! Write timer report
  call timer_report( printNow=.true. )

#ifdef MPI
  call MPI_finalize(ierr)
#endif

end program unfold

