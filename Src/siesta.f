      Program SIESTA

C ********************************************************************
C SIESTA Density Functional LCAO program package.
C Copyright by E.Artacho, P.Ordejon, D.Sanchez-Portal and 
C   J.M.Soler, 1996, 1997.
C The use of this program is allowed for non-for-profit research only.
C Copy or disemination of all or part of this package is not
C permitted without prior and explicit authorization by the authors.
C Send comments/suggestions/bug-reports to:
C   Jose M. Soler: jose.soler@uam.es
C   Pablo Ordejon: ordejon@condmat2.ciencias.uniovi.es
C ********************************************************************
C This program uses the flexible data format (FDF) package,
C copyrighted by A.Garcia and J.M.Soler, 1996, 1997
C ********************************************************************
C This program uses exchange-correlation routines copyrighted
C by C.Balbas, J.L.Martins and J.M.Soler, 1996, 1997.
C ********************************************************************
C Routine idiag copyrighted by J.Moreno and J.M.Soler, 1990.
C ********************************************************************
C Routine pulayx copyrighted by In-Ho Lee and P.Ordejon, 1997.
C ********************************************************************
C This program uses routines from 'Numerical Recipes', copyrighted
C by W.H.Press, A.A.Teukolsky, W.T.Veterlig and B.P.Flannery.
C Some of these routines have been modified or their name changed.
C ********************************************************************
C The atomic part uses old routines from unknown authors. If you wrote
C any of these routines, please let us know, to include the apropriate
C acknowledgements.
C ********************************************************************
C This module was written by P.Ordejon and J.M.Soler, 1996, 1997.
C ********************************************************************
C    6  10        20        30        40        50        60        7072

      implicit none

      include 'siesta.h'
      include 'fdf/fdfdefs.h'

      integer
     .  amax, auxdim, fincoor,
     .  i, i1, i2, ia, iaKB(maxkb), ianneal, iaorb(maxo),
     .  idyn, ierror, ifa, ifinal, ik, ikb, ilv, in,
     .  indxua(maxa), indxuo(maxo),
     .  inicoor, integers, io, ioa, ioptlwf, iord, iphKB(maxkb), 
     .  iphorb(maxo), iquench, is, isa(maxa), isel, iscf, 
     .  isolve, ispin, istp, istart, 
     .  istep, istr, iua, iunit, iuo, iv, ix, iza(maxa), izs(maxs), 
     .  j, ja, jamin, jna(maxna), jo, jx, kbmax, kbamax, kscell(3,3),
     .  lastc, lastkb(0:maxa), lasto(0:maxa), lc(0:1), linmin,
     .  listh(maxno,maxo), listhold(maxno,maxo), 
     .  lmax, lmaxs(maxs), lmxkbs(maxs), maxsav, mullipop, mua,
     .  na, namax, nbcell, nbk, ncells, ncgmax,
     .  ni, nkbs(maxs), nkpnt, nmove,
     .  nn, nnia, no, nokb, nomax, nos(maxs), nr, ns, nscf, nspin,
     .  nua, numel, numh(maxo), numhold(maxo), nuo, nv, 
     .  nzls(0:maxl,maxs), omax, one, osmax, smax, spnmax, zetmax

      double precision
     .  amass(maxa), Ang, aux(2,maxo), bcell(3,3), bk(3,maxbk), bulkm, 
     .  cell(3,3), cfa(3,maxa), cgaux(3*maxa*2), cgcntr(0:20), 
     .  contrf(maxzet,0:maxl,maxs), cstress(3,3),
     .  Datm(maxo), Debye, dipol(3), dDmax, 
     .  dDtol, DEna,  Dold(maxno,maxo,maxspn), Dscf(maxno,maxo,maxspn), 
     .  Dscfsave(maxno,maxo,maxspn), dt, DUext, DUscf, Dxc, dxmax,
     .  e1, e2, ebk(maxuo,maxbk,maxspn), Ecorrec, ef, Eharrs, Eions,
     .  Ekin, Ekinion, Ena, Enaatm, Enascf, 
     .  Enl, Entrop, eo(maxuo,maxspn,maxk), 
     .  Eold(maxno,maxo,maxspn), Escf(maxno,maxo,maxspn), 
     .  eta, etol, Etot, eV, Exc, E0,
     .  fa(3,maxa), factor, fdf_convfac,
     .  fmax, fmean, FreeE, fres, ftem, ftol, ftot(3),
     .  g2max, g2maxsave, getot,
     .  H(maxno,maxo,maxspn), H0(maxno,maxo),
     .  kBar, kcutof, kdispl(3), Kelvin,
     .  kn, kpoint(3,maxk), kweight(maxk), kpr, 
     .  mn, mpr, Pint, Pmol, Psol
      double precision
     .  qa(maxa), qo(maxuo,maxspn,maxk), qos(maxos,maxs),
     .  qspin, qsol, qtot,
     .  rc, rckb(maxkb), rcls(maxzet,0:maxl,maxs), rco(maxo),
     .  rcoor, rcut, reals(2), rijmin, rh,
     .  rmax, rmaxkb, rmaxo, rmaxv, rmin, r2ij(maxna), r2min,
     .  S(maxno,maxo), smass(maxs), ssign, stress(3,3),
     .  taurelax, temp, tempinit, tempion, tiny, tp, tt,
     .  Uatm, ucell(3,3), uion, Uscf,
     .  va(3,maxa), values(2), vcell(3,3), virial, vn,
     .  volcel, volume, vpr, we, wmix, wo,
     .  xa(3,maxa), xij(3,maxna), xijo(3,maxno,maxo), xmax, xmin

      real*8
     .  aux1(dimaux), aux2(dimaux)

      logical
     .  default, first, found, gamma, inspn, itest, last,
     .  negl, overflow, pulfile, relaxd,
     .  same, savehs, savevh, savevt, savrho,
     .  usesavelwf, usesavedm, usesavexv, writedim

      character
     .  filevh*25, filevt*25, filrho*25,
     .  line*150, names*80, paste*25,
     .  slabel*20, sname*150, shape*10, atm_label(maxs)*20

      external
     .  anneal, chkdim, conjgr, fixed,
     .  dhscf, diagon, dnaefs, extrapol,
     .  fdf_convfac, hsparse, initatom, iodm, ioxv, kgrid, kinefsm,
     .  mixing, mulliken, naefs, neighb, nlefsm, nose, npr, 
     .  ordern, overfsm, parse, paste, prmem, rcut, redata, shaper,
     .  timer, uion, verlet2, volcel, xijorb
c
      data
     .  e1, e2 / 1.d0, -1.d0 /
     .  kcutof, kscell, kdispl / 0.d0, 9*0, 3*0.d0 /
     .  nkpnt, (kpoint(i,1),i=1,3), kweight(1) / 1, 3*0.d0, 1.d0 /
     .  relaxd /.false./
     .  tiny   /1.d-15/
     .  amax, auxdim, lmax, kbmax, kbamax, namax, nbk, nomax /8*1/
     .  nuo, omax, osmax, smax, spnmax, zetmax /6*1/
     .  overflow / .false. /
c    .  pi / 3.1415926d0 /
c
c---------------------------------------------------------------------
c
c      Print version information
c
      call prversion
c
C Start time counter ..................................................
      call timer( 'siesta', 0 )
      call timer( 'siesta', 1 )
C ..................

C Print initial date and time .........................................
      call prdate( 'siesta' )
C ....................

C Factors of conversion to internal units (Bohr,Ry) ...................
      Ang    = 1.d0 / 0.529177d0
      eV     = 1.d0 / 13.60580d0
      kBar   = 1.d0 / 1.471d5
      Kelvin = eV / 11604.45d0
      Debye  = 0.393430d0
C ..................

C Print array sizes ...................................................
      call prmem( 0, 'siesta', 'Dold',     'd', maxno*maxo*maxspn )
      call prmem( 0, 'siesta', 'Dscf',     'd', maxno*maxo*maxspn )
      call prmem( 0, 'siesta', 'Dscfsave', 'd', maxno*maxo*maxspn )
      call prmem( 0, 'siesta', 'Eold',     'd', maxno*maxo*maxspn )
      call prmem( 0, 'siesta', 'Escf',     'd', maxno*maxo*maxspn )
      call prmem( 0, 'siesta', 'H',        'd', maxno*maxo*maxspn )
      call prmem( 0, 'siesta', 'H0',       'd', maxno*maxo        )
      call prmem( 0, 'siesta', 'listh',    'i', maxno*maxo        )
      call prmem( 0, 'siesta', 'listhold', 'i', maxno*maxo        )
      call prmem( 0, 'siesta', 'S',        'd', maxno*maxo*maxspn )
      call prmem( 0, 'siesta', 'xijo',     'd', 3*maxno*maxo      )
      call prmem( 0, 'siesta', ' ',        ' ', 0                 )
C .....................

C Read simulation data ................................................
      call redata(maxa, maxl, maxs, maxspn, maxzet,
     .            amax, lmax, smax, spnmax, zetmax, overflow,
     .            slabel, sname,
     .            na, ns, nspin, isa, izs, lmaxs, lmxkbs, 
     .            nzls, rcls, contrf, atm_label,
     .            cell, ucell, indxua, xa, smass, g2maxsave, negl, 
     .            nscf, dDtol, wmix, isolve,
     .            temp, ncgmax, ftol, eta, etol, rcoor, ioptlwf,
     .            idyn, istart, ifinal, nmove, ianneal, iquench,
     .            dt,  dxmax, tt, tp, mn, mpr, bulkm, taurelax,
     .            writedim, usesavelwf, usesavedm, mullipop,
     .            inspn, maxsav, pulfile, tempinit)

      g2max = g2maxsave
      if (overflow) goto 444
C ................

C Find some switches ..................................................
      default   = fdf_boolean('UseSaveData',.false.)
      usesavexv = fdf_boolean('MD.UseSaveXV',default)
      savehs    = fdf_boolean('SaveHS',.false.)
C .....................

C Read cell shape and atomic positions from a former run ..............
      if (usesavexv) then
        call ioxv( 'read', cell, vcell, na, isa, iza, xa, va, found )
        if (.not.found) write(6,'(/,a)')
     .    'siesta: WARNING: XV file not found'
      endif
C ..................

C Initialize pseudopotentials, atomic orbitals and atom lists .........
      call initatom(ns, na, izs, smass, lmxkbs, 
     .              lmaxs, maxl, lmax,
     .              nzls, maxos, maxzet, zetmax, maxo, maxkb,
     .              rcls, contrf, isa, atm_label,
     .              overflow, omax, osmax, kbmax,
     .              nkbs, nos, qos,
     .              no, nokb, qtot, rmaxv, rmaxo, rmaxkb,
     .              lasto, lastkb, iza, amass,
     .              iaorb, iphorb, Datm, qa, iaKB, iphKB)
      if (overflow) goto 444
C ..................

C Check dimensions for routine matel ..................................
      call matel0( ns, lmaxs, lmxkbs, maxl, nzls )
C ..................

C Find cell volume ....................................................
      volume = volcel( cell )
C ..................

C Automatic cell generation ...........................................
      if (volume .lt. 1.d-8) then
        do iv = 1,3
          do ix = 1,3
            cell(ix,iv) = 0.d0
            ucell(ix,iv) = 0.d0
          enddo
        enddo
        do ix = 1,3
          xmin =  1.d30
          xmax = -1.d30
          do ia = 1,na
            is = isa(ia)
            rc = rcut(is,0)
            xmin = min( xmin, xa(ix,ia)-rc )
            xmax = max( xmax, xa(ix,ia)+rc )
          enddo
C         Use a 10% margin for atomic movements
          cell(ix,ix) = 1.10d0 * (xmax - xmin)
          ucell(ix,ix) = cell(ix,ix)
        enddo
        volume = volcel( cell )
        write(6,'(/,a,3(/,a,3f12.6))')
     .    'siesta: Automatic unit cell vectors (Ang):',
     .    ('siesta:', (cell(ix,iv)/Ang,ix=1,3), iv =1,3)
      endif
C ..................

C Find system shape ...................................................
      call shaper( cell, na, isa, xa, shape, nbcell, bcell )
      write(6,'(/,2a)') 'siesta: System type = ', shape
C ......................

C Find indxuo .......................................
      ncells = nint( volcel(cell) / volcel(ucell) )
      nua = na / ncells
      nuo = no / ncells
      mua = 0
      do ia = 1,na
        iua = indxua(ia)
        iuo = lasto(iua-1)
        do io = lasto(ia-1)+1,lasto(ia)
          iuo = iuo + 1
          indxuo(io) = iuo
        enddo
        mua = max( nua, iua )
      enddo
      if (mua .ne. nua) stop 'siesta: ERROR: incorrect indxua'
      if (nuo .gt. maxuo) overflow = .true.
C ..................

C Some printout for debugging .......................
*     write(6,'(/,a)') 'siesta: indxua, indxuo ='
*     do ia = 1,na
*       write(6,'(i6,3x,20i4)')
*    .    indxua(ia), (indxuo(io),io=lasto(ia-1)+1,lasto(ia))
*     enddo
C ..................

C Find kbmax and kbamax .............................
      kbmax = lastkb(na)
      kbamax = 0
      do ia = 1,na
        kbamax = max( kbamax, lastkb(ia)-lastkb(ia-1) )
      enddo
C ......................

C Find rco and rckb .............................
      do ia = 1,na
        is = isa(ia)
        do io = lasto(ia-1)+1,lasto(ia)
          ioa = iphorb(io)
          rco(io) = rcut(is,ioa)
        enddo
        do ikb = lastkb(ia-1)+1,lastkb(ia)
          ioa = iphKB(ikb)
          rckb(ikb) = rcut(is,ioa)
        enddo
      enddo
C ......................

C Find k-grid for Brillouin zone integration ..........................
      nkpnt = maxk
      call kgrid( ucell, kscell, kdispl,
     .            kcutof, nkpnt, kpoint, kweight )
      if (nkpnt .gt. maxk) overflow = .true.
C ......................

c Find number of band k-points ........................................
      nbk = 0
      call bands( no, nspin, maxno, maxo, 0,
     .            numh, listh, H, S, xijo, maxuo, indxuo,
     .            Ef, .false., nbk, bk, ebk )
      if (nbk .gt. maxbk) overflow = .true.
C ......................

C Find if only gamma point is used ....................................
      if (nkpnt.eq.1 .and. abs(kpoint(1,1)).lt.tiny .and.
     .                     abs(kpoint(2,1)).lt.tiny .and.
     .                     abs(kpoint(3,1)).lt.tiny) then
        gamma = .true.
      else
        gamma = .false.
      endif
C ....................

C Print k-points ......................................................
      if (.not.gamma .and. nkpnt.le.maxk) then
        write(6,'(/,a)')
     .   'siesta: k-point coordinates (Bohr**-1) and weights:'
        write(6,'(a,i4,3f12.6,3x,f12.6)')
     .    ('siesta: ', ik, (kpoint(ix,ik),ix=1,3), kweight(ik),
     .    ik=1,nkpnt)
        write(6,'(a,f12.6,a)')
     .    'siesta: k-mesh cutoff =', kcutof/Ang, ' Ang'
      endif
C ....................

C Initialize atomic velocities to zero ................................
      if (.not. usesavexv) then
        do ia = 1,na
          do ix = 1,3
            va(ix,ia) = 0.0
          enddo
        enddo
      endif
C ..................

C Begin of coordinate relaxation iteration ============================
C Notice that this loop is not indented
      if (idyn .eq. 0) then
        inicoor = 0
        fincoor = nmove
      else 
        inicoor = istart
        fincoor = ifinal
      endif

C Build initial velocities according to Maxwell-Bolzmann distribution....
      if (idyn .ne. 0 .and. (.not. usesavexv)) 
     .               call vmb(na,tempinit,amass,va)
C ...

      istp = 0
      do istep = inicoor,fincoor
      call timer( 'IterMD', 1 )
      istp = istp + 1
      write(6,'(/,a)') 'siesta:    ==============================='
      if (idyn .eq. 0) 
     . write(6,'(a,i6)') 'siesta:        Begin CG move = ',istep
      if (idyn .ne. 0) 
     . write(6,'(a,i6)') 'siesta:        Begin MD step = ',istep
      write(6,'(a)')   'siesta:    ==============================='

C Print atomic coordinates ............................................
      call outcoor(cell, xa, isa, na, atm_label, ' ')
C ...................

C Print unit cell for Parrinello-Rahman (with or w/o Nose).............
      if((idyn.eq.3).or.(idyn.eq.4)) call outcell(cell)
C ...................

C Set grid cutoff and volume ..........................................
      g2max = g2maxsave
      volume = volcel( cell )
C ...................



C Initialize neighb subroutine ........................................
      ia = 0
      isel = 0
      rmax = max( 2.d0*rmaxv, 2.d0*rmaxo, rmaxo+rmaxkb )
      if (negl) then
        rh = 2.d0*rmaxo
      else
        rh = 2.d0*rmaxo + 2.0*rmaxkb
      endif
      nnia = maxna
      call neighb( cell, rmax, na, xa, ia, isel,
     .             nnia, jna, xij, r2ij )
      namax = 0
      do ia = 1,na
        nnia = 0
        call neighb( cell, rmax, na, xa, ia, isel,
     .               nnia, jna, xij, r2ij )
        namax = max( namax, nnia )
      enddo
      if (namax .gt. maxna) then
C       Increase namax with safety margin when atoms move
        namax = namax + 0.10 * namax + 10
        overflow = .true.
      endif
      if (overflow) goto 444
C ..................

C Check if any two atoms are unreasonably close .......................
      rijmin = fdf_physical( 'WarningMinimumAtomicDistance',
     .                        1.d0, 'Bohr' )
      do ia = 1,na
        r2min = 1.d30
        jamin = 0
        nnia = maxna
        call neighb( cell, rmax, na, xa, ia, isel,
     .               nnia, jna, xij, r2ij )
        do j = 1,nnia
          ja = jna(j)
          if ( r2ij(j).lt.r2min .and. ja.ge.ia ) then
C           Check that it is not the same atom
            if ( ja.ne.ia .or. r2ij(j).gt.1.d-12 ) then
              r2min = r2ij(j)
              jamin = ja
            endif
          endif
        enddo
        rmin = sqrt( r2min )
        if ( rmin .lt. rijmin ) write(6,'(a,2i6,a,f12.6,a)')
     .    'siesta: WARNING: Atoms', ia, jamin, ' too close: rij =',
     .     rmin/Ang, ' Ang'
      enddo
C ..................

C List of nonzero hamiltonian matrix elements .........................
      nomax = maxno
      call hsparse( negl, cell, na, isa, xa,
     .              lasto, lastkb, iphorb, iphKB,
     .              nomax, numh, listh )
      if (nomax .gt. maxno) then
C       Increase nomax with safety margin when atoms move
        nomax = nomax + 0.10 * nomax + 40
        overflow = .true.
      endif
      if (overflow) goto 444
      if (writedim) goto 111
C ..................

C Some printout for debugging ........................................
*     write(6,'(/,a)') 'siesta: connected orbitals'
*     do io = 1,no
*       write(6,'(i6,4x,15i4)') 
*    .    io, (listh(j,io),j=1,numh(io))
*     enddo
*     write(6,*) ' '
C ..................

C Find vectors between orbital centers ................................
      if (.not.gamma .or. savehs) then
        nomax = maxno
        call xijorb( negl, cell, nua, na, xa,
     .               lasto, lastkb, rco, rckb,
     .               nomax, numh, listh, xijo )
      endif
C ..................

C Initialize density matrix ...........................................
C set density matrix for first step
      if (istp .eq. 1) 
     .   call initdm(Datm, Dscf, Dold, Escf, lasto, maxa,
     .               maxno, maxo, maxspn, na, no, nspin,
     .               numh, numhold, listh, listhold, iaorb,
     .               found, inspn, usesavedm)

C Extrapolate density matrix between steps
      itest = .false.
      iord = 1
      if (idyn .eq. 0) iord = 0
C  If DM has just been read from disk, 
C  call extrapol with istep = 2 and iord = 0
C  to make it update the structure of DM, if needed
      if (found .and. istp .eq. 1) then
        istp = 2
        iord = 0
        itest = .true.
      endif
      call extrapol(istp, iord, nspin, no, maxo, maxno,  numh, listh,
     .              aux, numhold, listhold, Dscfsave, Dscf)
C  If DM have just been read, restore istp
      if (itest) istp = 1
      itest = .false.
C ..................

C Check for Pulay auxiliary matrices sizes ...................................
111   continue
      auxdim = 1
      if (pulfile) then
        one = 1
        if (dimaux .ne. one)  then
          overflow = .true.
          goto 444
        endif
      else
        numel = 0
        do ispin = 1,nspin
        do io = 1,no
        do in = 1,numh(io)
          numel = numel+1
        enddo
        enddo
        enddo
        numel = numel * maxsav
        call chkdime(dimaux,numel,overflow,auxdim)
C       Increase auxdim with safety margin when atoms move
        auxdim = auxdim + 0.1 * auxdim + 10
        if (overflow) goto 444
      endif
      if (writedim) goto 444
C ....................

C Start of SCF iteration _____________________________________________
      first = .true.
      last  = .false.
      do iscf = 1, nscf
        if (iscf .eq. nscf) last = .true.
        call timer( 'IterSCF', 1 )
        
C Find overlap matrix ...............................................
        if (first) then
           call overfsm(na, no, cell, xa, rmaxo, maxo,
     .                 maxna, maxno, maxno, lasto, iphorb, isa, volume,
     .                 numh, listh, numh, listh, nspin, Escf, 
     .                 jna, xij, r2ij,
     .                 fa, stress, S)
        endif
C ..................

C Initialize Hamiltonian ........................................
        do ispin = 1,nspin
          do io = 1,no
            do j = 1,numh(io)
              H(j,io,ispin) = 0.0d0
            enddo
          enddo
        enddo
C ..................

C Initialize forces and stress ...................
        if (first.or.last) then
          do ia = 1,na
            do ix = 1,3
              fa(ix,ia) = 0.d0
            enddo
          enddo
          do ix = 1,3
            do jx = 1,3
              stress(jx,ix) = 0.d0
            enddo
          enddo
        endif
C ..................

C Self-energy of isolated ions ........................................
        if (first) then
          Eions = 0.d0
          do ia = 1,na
            is = isa(ia)
            Eions = Eions + uion(is)
          enddo
        endif
C ..................

C Neutral-atom: energy, forces and stress ............................
C First time for energy, last time for forces
        if (first.or.last) then
          call naefs(na, ns, cell, xa, rmaxv,
     .               maxna, isa, izs, volume,
     .               jna, xij, r2ij,
     .               Ena, fa, stress)
          call dnaefs(na, ns, cell, xa, rmaxv,
     .               maxna, isa, izs, volume,
     .               jna, xij, r2ij,
     .               DEna, fa, stress)
          Ena = Ena + DEna
        endif
C ..................

C Normalize density matrix to exact charge ...........................
        qsol = 0.d0
        do ispin = 1,nspin
          do io = 1,no
            do in = 1,numh(io)
	    qsol = qsol + Dscf(in,io,ispin) * s(in,io)
            enddo
          enddo
        enddo
        if (.not.first .and.
     .       abs(qsol/qtot-1.d0).gt.1.d-2) write(6,'(a,2f15.6)')
     .      'siesta: WARNING: Qtot, Tr[D*S] =', qtot, qsol
        do ispin = 1,nspin
          do io = 1,no
            do in = 1,numh(io)
              Dscf(in,io,ispin) = Dscf(in,io,ispin) * qtot/qsol
            enddo
          enddo
        enddo
C ..................

C Kinetic: energy, forces, stress and matrix elements .................
        if (first.or.last) then
          call kinefsm(na, no, cell, xa, rmaxo, maxo,
     .                 maxna, maxno, maxno, lasto, iphorb, isa, volume,
     .                 numh, listh, numh, listh, nspin, Dscf, 
     .                 jna, xij, r2ij,
     .                 Ekin, fa, stress, H)
        endif
C ..................

C Non-local-pseudop: energy, forces, stress and matrix elements .......
        if (first.or.last) then
          call nlefsm(cell, na, isa, xa, maxo, maxno, maxno,
     .                lasto, lastkb, iphorb, iphKB,
     .                numh, listh, numh, listh, nspin, Dscf, 
     .                Enl, fa, stress, H)
        endif
C ..................

C Save or get partial Hamiltonian (non-SCF part) ......................
        if (first.or.last) then
          do io = 1,no
            do j = 1,numh(io)
              H0(j,io) = H(j,io,1)
            enddo
          enddo
        else
          do ispin = 1,nspin
            do io = 1,no
              do j = 1,numh(io)
                H(j,io,ispin) = H0(j,io)
              enddo
            enddo
          enddo          
        endif
C ..................

C Non-SCF part of total energy .......................................
        if (first.or.last) then
          E0 = -Eions + Ena + Ekin + Enl
        else
          E0 = -Eions + Ena
          do ispin = 1,nspin
            do io = 1,no
              do j = 1,numh(io)
                E0 = E0 + H0(j,io) * Dscf(j,io,ispin)
              enddo
            enddo
          enddo
        endif
C ..................

C Add SCF contribution to energy and matrix elements ..................
        ilv = 0
        if (last) then
          ifa  = 1
          istr = 1
        else
          ifa  = 0
          istr = 0
        endif
        call dhscf( nspin, maxo, no, iaorb, iphorb,
     .              na, isa, xa, cell, g2max,
     .              ilv, ifa, istr, ' ', ' ', ' ',
     .              maxno, numh, listh, Dscf, Datm,
     .              maxno, numh, listh, H,
     .              Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, Exc, Dxc,
     .              dipol, fa, stress, ierror )
            
        if (istp.eq.1 .and. iscf.eq.1) write(6,'(/,a,f10.3,a)')
     .    'siesta: dhscf mesh cutoff =', g2max, ' Ry'

        if (ierror.ne.0) then
          if (ierror.eq.1) then
            write(6,'(/,a,/)')
     .        'siesta: BAD DIMENSIONS in DHSCF. RECOMPILE.'
          else
            write(6,'(/,a,I3,A)') 'siesta: ERROR', ierror, ' IN DHSCF.'
          endif
          goto 999
        endif
C ..................

C Orthonormalization forces ...........................................
        if (last) then
          call overfsm(na, no, cell, xa, rmaxo, maxo,
     .                 maxna, maxno, maxno, lasto, iphorb, isa, volume,
     .                 numh, listh, numh, listh, nspin, Escf, 
     .                 jna, xij, r2ij,
     .                 fa, stress, S)
        endif
C ..................

C Impose supercell translational symmetry (posibly broken by mesh) ....
        if (ncells .gt. 1) then
          do ia = nua+1,na
            iua = indxua(ia)
            do ix = 1,3
              fa(ix,iua) = fa(ix,iua) + fa(ix,ia)
            enddo
          enddo
          do ia = 1,nua
            do ix = 1,3
              fa(ix,ia) = fa(ix,ia) / ncells
            enddo
          enddo
          do ia = nua+1,na
            iua = indxua(ia)
            do ix = 1,3
              fa(ix,ia) = fa(ix,iua) 
            enddo
          enddo
        endif
C ..................

C Find entropy ........................................................
        if (isolve .eq. 0) then
          Entrop = 0.d0
          if (istp.gt.1 .or. iscf.gt.1) then
            do ik = 1,nkpnt
              do io = 1,nuo
                if (nspin .eq. 1) then
                  wo = qo(io,1,ik) / kweight(ik) * 0.5d0
                  we = 1.d0 - wo
                  wo = max( wo, tiny )
                  we = max( we, tiny )
                  Entrop = Entrop - 2.d0 * kweight(ik) *
     .                             ( wo*log(wo) + we*log(we) )
                else
                  do ispin = 1,2
                    wo = qo(io,ispin,ik) / kweight(ik)
                    we = 1.d0 - wo
                    wo = max( wo, tiny )
                    we = max( we, tiny )
                    Entrop = Entrop - kweight(ik) *
     .                               ( wo*log(wo) + we*log(we) )
                  enddo
                endif
              enddo
            enddo
            Entrop = Entrop * no / nuo
          endif
        endif
C ..................

C Save present density matrix ........................................
      do ispin = 1,nspin
        do io = 1,no
          do jo = 1,numh(io)
            Dold(jo,io,ispin) = Dscf(jo,io,ispin)
            Eold(jo,io,ispin) = Escf(jo,io,ispin)
          enddo
        enddo
      enddo
C ..................

C Save hamiltonian and overlap matrices ............................
      if (savehs) then
        call iohs( 'write', no, nspin, maxo, maxno,
     .             numh, listh, H, S, qtot, temp, xijo )
      endif
C ..................

C Solve eigenvalue problem .........................................
        if (isolve .eq. 0) then
          call diagon(no, nspin, maxspn, maxo, maxno, maxno,
     .                numh, listh, numh, listh, H, S,
     .                qtot, temp, e1, e2,
     .                xijo, maxuo, indxuo, nkpnt, kpoint, kweight,
     .                eo, qo, Dscf, Escf, ef)
          Ecorrec = 0.d0
C         Save density matrix on disk
          call iodm( 'write', maxno, maxo, no, nspin,
     .               numh, listh, Dscf, found )
        elseif (isolve .eq. 1) then
          call ordern(usesavelwf,ioptlwf,ns,na,no,lasto,maxl,lmaxs,
     .                nzls,isa,qa,rcoor,rh,cell,xa,iscf,istp,ncgmax,
     .                etol,eta,qtot,maxno,numh,listh,h,s,
     .                Dscf,Escf,Ecorrec)
        endif
C ..................

C Harris-functional energy ............................................
        Eharrs = 0.d0
        do ispin = 1,nspin
          do io = 1,no
            do j = 1,numh(io)
              Eharrs = Eharrs + H(j,io,ispin) *
     .                         ( Dscf(j,io,ispin) - Dold(j,io,ispin) )
            enddo
          enddo
        enddo
C ..................

C Mix input and output energy-density and density matrices ............
C  Pulay mixing implementation (ihlee) ......
        if(last) goto 994
        if (wmix .eq. 0.d0) wmix = 0.5
        if(wmix .le. 0.d0) wmix = -wmix
        if(maxsav.ge.2) then
          call pulayx(pulfile, iscf, no, maxo, maxno, numh, listh,
     .         nspin, maxsav, wmix, aux1, aux2, dimaux,
     .         Dscf, Dold, dDmax)

c          do ispin = 1,nspin
c            write(6,'(/,a)') 'siesta: Mulliken Analisys:'
c            if (nspin .eq. 2) then
c              if(ispin .eq. 1) write(6,'(/,a)') '        Spin UP '
c              if(ispin .eq. 2) write(6,'(/,a)') '        Spin DOWN '
c            endif
c           write(6,*) mullipop, na, no, maxno, numh, listh
c           call mulliken(1, na, no, maxno, numh, listh, S,
c     .                    Dscf(1,1,ispin), lasto)
c           enddo
        endif
C ........
C Linear mixing .....
994     continue
        if (last .or. maxsav .lt. 2) then

          do ispin = 1,nspin
            call mixing (iscf, no, maxo, maxno, numh, listh, wmix,
     .                   Escf(1,1,ispin), Eold(1,1,ispin), dDmax)
            call mixing (iscf, no, maxo, maxno, numh, listh, wmix,
     .                   Dscf(1,1,ispin), Dold(1,1,ispin), dDmax)

c            write(6,'(/,a)') 'siesta: Mulliken Analisys:'
c            if (nspin .eq. 2) then
c              if(ispin .eq. 1) write(6,'(/,a)') '        Spin UP '
c              if(ispin .eq. 2) write(6,'(/,a)') '        Spin DOWN '
c            endif
c            call mulliken(1, na, no, maxno, numh, listh, S,
c     .                Dscf(1,1,ispin), lasto)
          enddo

        endif
C ...................

C Print energies ......................................................
        DEna = Enascf - Enaatm
        Etot = E0 + DEna + DUscf + DUext + Exc + Ecorrec
        Eharrs = Etot + Eharrs
        FreeE  = Etot - Temp * Entrop
        if (istp.eq.1 .and. first) write(6,'(/,a,/,(a,f15.6))')
     .   'siesta: Program''s energy decomposition (eV):',
     .   'siesta: Eions   =', Eions/eV,
     .   'siesta: Ena     =', Ena/eV,
     .   'siesta: Ekin    =', Ekin/eV,
     .   'siesta: Enl     =', Enl/eV,
     .   'siesta: DEna    =', DEna/eV,
     .   'siesta: DUscf   =', DUscf/eV,
     .   'siesta: DUext   =', DUext/eV,
     .   'siesta: Exc     =', Exc/eV,
     .   'siesta: eta*DQ  =', Ecorrec/eV,
     .   'siesta: Eharris =', Eharrs/eV,
     .   'siesta: Etot    =', Etot/eV,
     .   'siesta: FreeEng =', FreeE/eV
C ...................

C Print total energy and density matrix error .........................
        if (isolve .eq. 0) then
          if (iscf .eq. 1) write(6,'(/,a12,3a14,2a8)')
     .    'siesta: iscf', '   Eharris(eV)', 
     .    '      E_KS(eV)', '   FreeEng(eV)', 
     .    '   dDmax', '  Ef(eV)'
          write(6,'(a8,i4,3f14.4,2f8.4)')
     .    'siesta: ',iscf, Eharrs/eV, Etot/eV, FreeE/eV, dDmax, Ef/eV
        endif
        if (isolve .eq. 1) then
          write(6,'(/,a15,i4)') 'siesta: iscf = ',iscf
          write(6,'(a14,f15.4,a13,f15.4,a10,f7.4/)') 
     .    'Eharris(eV) = ',Eharrs/eV,
     .    '  E_KS(eV) = ',Etot/eV,'  dDmax = ',dDmax
        endif
C ...................

C If last iteration, exit SCF loop ....................................
        if (last) then
          do ispin = 1,nspin
            do io = 1,no
              do jo = 1,numh(io)
                Dscf(jo,io,ispin) = Dold(jo,io,ispin)
                Escf(jo,io,ispin) = Eold(jo,io,ispin)
              enddo
            enddo
          enddo
*         call iodm( 'write', maxno, maxo, no, nspin,
*    .               numh, listh, Dscf, found )
          call timer( 'IterSCF', 2 )
          goto 50
        endif
C ...................

C If converged, make last iteration to find forces ....................
        if (dDmax.lt.dDtol) last = .true.
C ...................

        call timer( 'IterSCF', 2 )
        if (istep.eq.inicoor .and. first) call timer( 'IterSCF', 3 )
        first = .false.
      enddo
   50 continue
C End of SCF iteration_________________________________________________

C Write final Kohn-Sham Energy ........................................
      write(6,"(/a,f14.4)") 'siesta: E_KS(eV) = ', Etot/eV

C Write atomic forces .................................................
      fmax = 0.d0
      fres = 0.d0
      do ix = 1,3
        ftot(ix) = 0.d0
        do ia = 1,na
          ftem = fa(ix,ia)
          ftot(ix) = ftot(ix) + ftem
          fres = fres + ftem*ftem
          fmax = max( fmax, dabs(ftem) )
        enddo
      enddo
      fres = dsqrt( fres / (3.d0*na) )
      write(6,'(/,a)') 'siesta: Atomic forces (eV/Ang):'
      write(6,'(i4,3f12.6)') (ia,(fa(ix,ia)*Ang/eV,ix=1,3),ia=1,na)
      write(6,'(40(1h-),/,a4,3f12.6)') 'Tot',(ftot(ix)*Ang/eV,ix=1,3)
      write(6,'(40(1h-),/,a4, f12.6)') 'Max',fmax*Ang/eV
      write(6,'(a4,f12.6,a)')'Res',fres*Ang/eV,
     .                   '    [ sqrt( Sum f_i^2 / 3N ) ]'

C Write stress tensor for Parrinello-Rahman (with or w/o Nose).........
      if((idyn.eq.3).or.(idyn.eq.4)) then
         write(6,'(/,a,3(/,a,3f12.6))')
     .     'siesta: Stress tensor (eV/Ang**3):',
     .     ('     ',(stress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)

C Pressure (only for the solid)
      Psol = - (( stress(1,1) + stress(2,2) + stress(3,3) ) / 3.0)
      write(6,'(/,a,f20.8,a)')
     .  'siesta: Pressure:', Psol/kBar, '  kBar'

      endif
C ...................

C  Mulliken population analysis .......................................
      if (mullipop .ne. 0) then
        write(6,'(/,a)') 'siesta: Mulliken analysis:'
        do ispin = 1,nspin
          if (nspin .eq. 2) then
            if(ispin .eq. 1) write(6,'(a)') 'siesta: Spin UP '
            if(ispin .eq. 2) write(6,'(a)') 'siesta: Spin DOWN '
          endif
          call mulliken(mullipop, na, no, maxno, numh, listh, S,
     .                  Dscf(1,1,ispin), lasto)
        enddo
      endif
C ...................

C Impose constraints to atomic movements by changing forces ...........
      call fixed( cell, stress, na, isa, amass, xa, fa, cstress, cfa )
C ...................

C Move atoms ..........................................................
      if (idyn .eq. 0) then
        if (nmove .ne. 0) then
          if (istp .eq. 1) then
            relaxd = .false.
            cgcntr(0) = 0
            linmin = 1
          endif
          call conjgr( 3*na, xa, cfa, dxmax, ftol, cgcntr, cgaux )
          if (nint(cgcntr(1)) .ne. linmin) then
            write(6,'(a,i4,a,f10.4)')
     .       'siesta: Finished line minimization ', linmin,
     .       '.  Mean atomic displacement =', cgcntr(18)/sqrt(dble(na))
            linmin = nint(cgcntr(1))
          endif
          if (nint(cgcntr(0)) .eq. 0) then
            relaxd = .true.
C           Exit coordinate relaxation loop
            call timer( 'IterMD', 2 )
            goto 60
          endif
        endif
      endif

      Ekinion  = 0.d0
      vn       = 0.d0
      vpr      = 0.d0
      kn       = 0.d0
      kpr      = 0.d0

      iunit = 2

      if (idyn .eq. 1) then
        call verlet2(istp, iunit, iquench, na, cfa, dt, 
     .              amass, va, xa, Ekinion, tempion)
      
      elseif (idyn .eq. 2) then
        call nose(istp, iunit, na, cfa, tt, dt, amass, mn, 
     .           va, xa, Ekinion, kn, vn, tempion)

      elseif (idyn .eq. 3) then
        call pr(istp, iunit, iquench, na, cfa, cstress, tp, dt, amass, 
     .          mpr, va, xa, cell, Ekinion, kpr, vpr, tempion, Pint)
        write(6,'(/,a,f12.3,a)')
     .    'siesta: E_kin PR =', kpr/Kelvin, ' K'
      
      elseif (idyn .eq. 4) then
        call npr(istp, iunit, na, cfa, cstress, tp, tt, dt, amass, mn, 
     .           mpr, va, xa, cell, Ekinion, kn, kpr, vn, vpr, tempion, 
     .           Pint)

      elseif (idyn .eq. 5) then
        call anneal(istp, iunit, ianneal, taurelax, bulkm,
     .              na, cfa, cstress, tp, tt, dt, amass,
     .              va, xa, cell, Ekinion, tempion, Pint)
      endif

      if (idyn .gt. 0) then
        write(6,'(/,a,f12.3,a)')
     .    'siesta: Temp_ion =', tempion, ' K'
      endif

C ...................

C Save atomic positions and velocities ...............................
      call ioxv( 'write', cell, vcell, na, isa, iza, xa, va, found )
C ...................

      getot = Etot + Ekinion + kn + kpr + vn + vpr
      call timer( 'IterMD', 2 )
      enddo
   60 continue
C End of coordinate-relaxation loop ===================================


C Print atomic coordinates (and also unit cell for ParrRah.)
      if (nmove .ne. 0) then
        if (relaxd) 
     .            call outcoor(cell, xa, isa, na, atm_label, 'Relaxed')
        if (.not.relaxd) 
     .  call outcoor(cell, xa, isa, na, atm_label, 'Final (unrelaxed)')
        if((idyn.eq.3).or.(idyn.eq.4)) call outcell(cell)
      endif


C Print coordinates in xmol format in a separate file

      if(fdf_boolean('WriteCoorXmol',.false.)) 
     .   call coxmol(iza, xa, na, slabel)

C Print coordinates in cerius format in a separate file

      if(fdf_boolean('WriteCoorCerius',.false.))
     .   call coceri(iza, xa, cell, na, sname, slabel)


c Find and print bands
      nbk = 0
      call bands( no, nspin, maxno, maxo, maxbk,
     .            numh, listh, H, S, xijo, maxuo, indxuo,
     .            Ef, .true., nbk, bk, ebk )
      if (nbk.gt.0 .and. nbk.le.maxbk) then
        write(6,'(/,a,/,a4,a12)')
     .   'siesta: Band k vectors (Bohr**-1):', 'ik', 'k'
        do ik = 1,nbk
          write(6,'(i4,3f12.6)') ik, (bk(ix,ik),ix=1,3)
        enddo
        write(6,'(/,a,/,a4,a3,a7)')
     .   'siesta: Band energies (eV):', 'ik', 'is', 'eps'
        do ispin = 1,nspin
          do ik = 1,nbk
            write(6,'(i4,i3,10f7.2)')
     .        ik, ispin, (ebk(io,ik,ispin)/eV,io=1,min(10,nuo))
            if (nuo.gt.10) write(6,'(7x,10f7.2)')
     .          (ebk(io,ik,ispin)/eV,io=11,nuo)
          enddo
        enddo
      endif

C Print eigenvalues
      if (isolve.eq.0 .and. nuo.lt.1000) then
        write(6,'(/,a,/,2a3,a7)')
     .   'siesta: Eigenvalues (eV):', 'ik', 'is', 'eps'
        do ik = 1,nkpnt
          do ispin = 1,nspin
            write(6,'(2i3,10f7.2)')
     .        ik, ispin, (eo(io,ispin,ik)/eV,io=1,min(10,nuo))
            if (nuo.gt.10) write(6,'(6x,10f7.2)')
     .          (eo(io,ispin,ik)/eV,io=11,nuo)
          enddo
        enddo
        write(6,'(a,f15.6,a)') 'siesta: Fermi energy =', ef/eV, ' eV'
      endif

C Print program's energy decomposition
      write(6,'(/,a,/,(a,f15.6))')
     .   'siesta: Program''s energy decomposition (eV):',
     .   'siesta:-Eions   =', (-Eions)/eV,
     .   'siesta: Ena     =', Ena/eV,
     .   'siesta: Ekin    =', Ekin/eV,
     .   'siesta: Enl     =', Enl/eV,
     .   'siesta: DEna    =', DEna/eV,
     .   'siesta: DUscf   =', DUscf/eV,
     .   'siesta: DUext   =', DUext/eV,
     .   'siesta: Exc     =', Exc/eV,
     .   'siesta: eta*DQ  =', Ecorrec/eV,
     .   'siesta: Ekinion =', Ekinion/eV,
     .   'siesta: Eharris =', (Eharrs+Ekinion)/eV,
     .   'siesta: Etot    =', (Etot+Ekinion)/eV,
     .   'siesta: FreeEng =', (FreeE+Ekinion)/eV

C Print standard energy decomposition
      write(6,'(/,a)') 'siesta: Final energy (eV):'
      write(6,'(a,a15,f15.6)')
     .  'siesta: ',      'Kinetic =', Ekin/eV,
     .  'siesta: ',      'Hartree =', Uscf/eV,
     .  'siesta: ',   'Ext. field =', DUext/eV,
     .  'siesta: ',  'Exch.-corr. =', Exc/eV,
     .  'siesta: ', 'Ion-electron =', (Enascf+Enl+DUscf-Uscf-Uatm)/eV,
     .  'siesta: ',      'Ion-ion =', (Ena+Uatm-Enaatm-Eions)/eV,
     .  'siesta: ',      'Ekinion =', Ekinion/eV,
     .  'siesta: ',        'Total =', (Etot+Ekinion)/eV

C Print atomic forces
      fmax = 0.d0
      do ix = 1,3
        ftot(ix) = 0.d0
        do ia = 1,na
          fmax = max( fmax, dabs(fa(ix,ia)) )
          ftot(ix) = ftot(ix) + fa(ix,ia)
        enddo
      enddo
      if (fmax .gt. ftol) then
        write(6,'(/,a)') 'siesta: Atomic forces (eV/Ang):'
        write(6,'(a,i4,3f12.6)')
     .    ('siesta: ', ia,(fa(ix,ia)*Ang/eV,ix=1,3),ia=1,na)
        write(6,'(a,40(1h-),/,a,a4,3f12.6)')
     .    'siesta: ','siesta: ','Tot',(ftot(ix)*Ang/eV,ix=1,3)
      endif

C Print constrained atomic forces
      same = .true.
      do ia = 1,na
        do ix = 1,3
          if (cfa(ix,ia) .ne. fa(ix,ia)) same = .false.
        enddo
      enddo
      if (.not.same) then
        fmax = 0.d0
        do ix = 1,3
          ftot(ix) = 0.d0
          do ia = 1,na
            fmax = max( fmax, dabs(cfa(ix,ia)) )
            ftot(ix) = ftot(ix) + cfa(ix,ia)
          enddo
        enddo
        if (fmax .gt. ftol) then
          write(6,'(/,a)') 'siesta: Constrained forces (eV/Ang):'
          write(6,'(a,i4,3f12.6)')
     .      ('siesta: ', ia,(cfa(ix,ia)*Ang/eV,ix=1,3),ia=1,na)
          write(6,'(a,40(1h-),/,a,a4,3f12.6)')
     .      'siesta: ','siesta: ','Tot',(ftot(ix)*Ang/eV,ix=1,3)
        endif
      endif

C Print stress tensor
      write(6,'(/,a,3(/,a,3f12.6))')
     .  'siesta: Stress tensor (eV/Ang**3):',
     .  ('siesta: ',(stress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)

C Print constrained stress tensor
      same = .true.
      do ix = 1,3
        do jx = 1,3
          if (cstress(jx,ix) .ne. stress(jx,ix)) same = .false.
        enddo
      enddo
      if (.not.same)
     .  write(6,'(/,a,3(/,a,3f12.6))')
     .    'siesta: Constrained stress tensor (eV/Ang**3):',
     .    ('siesta: ',(cstress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)

C Find pressure
      virial = 0.d0
      do ix = 1,3
        fmean = 0.d0
        do ia = 1,na
          fmean = fmean + fa(ix,ia) / na
        enddo
        do ia = 1,na
          virial = virial + xa(ix,ia) * (fa(ix,ia) - fmean)
        enddo
      enddo
      Psol = - (( stress(1,1) + stress(2,2) + stress(3,3) ) / 3.0)
      Pmol = Psol - virial / volume / 3.0
      write(6,'(/,a,f12.6,a)')
     .  'siesta: Cell volume =', volume/Ang**3, ' Ang**3'
      write(6,'(/,a,/,a,2a20,a,3(/,a,2f20.8,a))')
     .  'siesta: Pressure:',
     .  'siesta: ','Solid',        'Molecule',      '  Units',
     .  'siesta: ', Psol,           Pmol,           '  Ry/Bohr**3',
     .  'siesta: ', Psol*Ang**3/eV, Pmol*Ang**3/eV, '  eV/Ang**3',
     .  'siesta: ', Psol/kBar,      Pmol/kBar,      '  kBar'

c Print spin polarization
      if (nspin .eq. 2) then
        i1 = 1
        i2 = 2
        qspin = 0.d0
        do ispin = 1,nspin
          if (ispin .eq. 1) ssign = +1
          if (ispin .eq. 2) ssign = -1
          do io = 1,no
            do j = 1,numh(io)
              jo = listh(j,io)
              qspin = qspin + ssign * Dscf(j,io,ispin) * S(j,io)
            enddo
          enddo
        enddo
        write(6,'(/,a,f12.6)')
     .   'siesta: Total spin polarization (Qup-Qdown) =', qspin
      endif

c Print electric dipole
      if (shape .ne. 'bulk') then
        write(6,'(/,a,3f12.6)')
     .    'siesta: Electric dipole (a.u.)  =', dipol
        write(6,'(a,3f12.6)')
     .    'siesta: Electric dipole (Debye) =', 
     .    (dipol(ix)/Debye,ix=1,3)
      endif

c Save electron density and potential
      savrho = fdf_boolean( 'SaveRho',                    .false. )
      savevh = fdf_boolean( 'SaveElectrostaticPotential', .false. )
      savevt = fdf_boolean( 'SaveTotalPotential',         .false. )
      if (savrho .or. savevh .or. savevt) then
        filrho = ' '
        filevh = ' '
        filevt = ' '
        if (savrho) filrho = paste( slabel, '.RHO' )
        if (savevh) filevh = paste( slabel, '.VH'  )
        if (savevt) filevt = paste( slabel, '.VT'  )
        do ispin = 1,nspin
          do io = 1,no
            do j = 1,numh(io)
              H(j,io,ispin) = H0(j,io)
            enddo
          enddo
        enddo          
        call dhscf( nspin, maxo, no, iaorb, iphorb,
     .              na, isa, xa, cell, g2max,
     .              ilv, 0, 0, filrho, filevh, filevt,
     .              maxno, numh, listh, Dscf, Datm,
     .              maxno, numh, listh, H,
     .              Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, Exc, Dxc,
     .              dipol, fa, stress, ierror )
      endif

c Find local density of states
      if ( fdf_block('LocalDensityOfStates',iunit) ) then

c       Find the desired energy range
        read(iunit,'(a)') line
        lastc = index(line,'#') - 1
        if (lastc .eq. -1) lastc = len(line)
        call parse( line(1:lastc), 
     .              nn, lc, names, 
     .              nv, values,
     .              ni, integers,
     .              nr, reals )
        factor = fdf_convfac( names(1:lc(1)), 'Ry' )
        e1 = values(1) * factor
        e2 = values(2) * factor

c       Find the density matrix for states between e1 and e2
        if (isolve .eq. 0) then
          call diagon(no, nspin, maxspn, maxo, maxno, maxno,
     .                numh, listh, numh, listh, H, S, 
     .                qtot, temp, e1, e2,
     .                xijo, maxuo, indxuo, nkpnt, kpoint, kweight,
     .                eo, qo, Dscf, Escf, ef)
        else
          write(6,*)
     .     'siesta: ERROR: LDOS implemented only with diagon'
          goto 70
        endif
     
c       Find the LDOS in the real space mesh
        filrho = paste( slabel, '.LDOS' )
        do ispin = 1,nspin
          do io = 1,no
            do j = 1,numh(io)
              H(j,io,ispin) = H0(j,io)
            enddo
          enddo
        enddo          
        call dhscf( nspin, maxo, no, iaorb, iphorb,
     .              na, isa, xa, cell, g2max,
     .              ilv, 0, 0, filrho, ' ', ' ',
     .              maxno, numh, listh, Dscf, Datm,
     .              maxno, numh, listh, H,
     .              Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, Exc, Dxc,
     .              dipol, fa, stress, ierror )
      endif
   70 continue

c Stop time counter
      call timer( 'siesta', 2 )
      call timer( 'all', 3 )

C Print final date and time
      call prdate( 'siesta' )

c Go to exit point
      goto 999

C Bad-dimensions exception handling
  444 continue

      open( 1, file='siesta.h', status='unknown' )
      write(1,'(a)') 'c Dimension parameters for siesta:'
      write(1,'(a)')
     .'c  maxa   : Maximum number of atoms'
      write(1,'(a)')
     .'c  maxbk : Maximum number of band k points'
      write(1,'(a)')
     .'c  maxk   : Maximum number of k points'
      write(1,'(a)') 
     .'c  maxkb  : Maximum total number of Kleinman-Bylander projectors'
      write(1,'(a)') 
     .'c  maxkba : Maximum number of KB projectors of one atom'
      write(1,'(a)') 
     .'c  maxl   : Maximum angular momentum'
      write(1,'(a)') 
     .'c  maxna  : Maximum number of neighbours of any atom'
      write(1,'(a,/,a)') 
     .'c  maxno  : Max basis orbitals connected to any given one,',
     .'c           either directly or through a KB projector'
      write(1,'(a)') 
     .'c  maxo   : Maximum total number of basis orbitals'
      write(1,'(a)') 
     .'c  maxos  : Maximum number of basis orbitals of any atom'
      write(1,'(a)') 
     .'c  maxs   : Maximum number of species'
      write(1,'(a)') 
     .'c  maxspn : Maximum number of spin polarizations (1 or 2)'
      write(1,'(a)') 
     .'c  maxuo  : Maximum number of basis orbitals in unit cell'
      write(1,'(a)')
     .'c  maxzet : Maximum number of orbital zetas'
      write(1,'(a)')
     .'c  dimaux : Dimension of Pulay auxiliary vectors'

      write(1,'(6x,a)')'integer maxa, maxbk, maxk, maxkb, maxkba, maxl'
      write(1,'(6x,a)')'integer maxna, maxno, maxo, maxos'
      write(1,'(6x,a)')'integer maxs, maxspn, maxuo, maxzet,dimaux'
      write(1,'(6x,a,i12,a)') 'parameter ( maxa   =', amax,      ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxbk  =', max(nbk,1),' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxk   =', nkpnt,     ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxkb  =', kbmax,     ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxkba =', kbamax,    ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxl   =', lmax,      ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxna  =', namax,     ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxno  =', nomax,     ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxo   =', omax,      ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxos  =', osmax,     ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxs   =', smax,      ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxspn =', spnmax,    ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxuo  =', max(nuo,1),' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxzet =', zetmax,    ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( dimaux =', auxdim,    ' )'
      close(1)
      if (writedim) write(6,'(/,a)')
     .   'siesta: Minimum dimensions dumped to siesta.h'
      if (overflow) write(6,'(/,a,/)')
     .   'siesta: BAD DIMENSIONS. RECOMPILE.'

c Exit point
  999 continue

c Stop memory counter
      call prmem( 0, ' ', ' ', ' ', 0 )

      if (overflow) stop 'siesta: BAD DIMENSIONS. RECOMPILE.'
      end


