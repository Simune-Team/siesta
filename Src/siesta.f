C $Id: siesta.f,v 1.67 1999/05/26 09:46:03 emilio Exp $

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
C This module was written by P.Ordejon and J.M.Soler, 1996, 1997, 1998.
C ********************************************************************
C    6  10        20        30        40        50        60        7072

      implicit none

      include 'siesta.h'
      include 'fdf/fdfdefs.h'

      integer
     .  fincoor,
     .  i, ia, ia1, ia2, iadispl, iaKB(maxkb), ianneal, iaorb(maxo),
     .  idyn, ierror, ifa, ifinal, ihmat, iiscf, ik, ikb, ilv, in,
     .  indxua(maxa), indxuo(maxo),
     .  inicoor, integers, io, ioa, ioptlwf, iord, iphKB(maxkb), 
     .  iphorb(maxo), iquench, is, isa(maxa), isel, iscf, 
     .  isolve, ispin, istp, istpsave, istart, istep, istr,
     .  iu, iunit,  iv, ix, ixdispl, iza(maxa),  
     .  j, ja, jamin, jna(maxna), jo, jx, kscell(3,3),
     .  lastc, lastkb(0:maxa), lasto(0:maxa), lc(0:1), 
     .  listh(maxno,maxo), listhold(maxno,maxo), 
     .  maxsav, mscell(3,3), mullipop,
     .  na, nkick, nnamax, nauxpul, nbcell, nbk, ncells, ncgmax,
     .  ni, nkpnt, nmove,
     .  nn, nnia, no, nokb, nnomax, nr,
     .  ns, nsc(3), nscf, nspin,
     .  ntm(3), nua, numh(maxo), numhold(maxo), nuo, nuokb, nv, 
     .  nxij, one, pmax

      double precision
     .  amass(maxua), Ang, aux(2,maxo), auxpul(maxpul,2),
     .  bcell(3,3), beta, bk(3,maxbk), bulkm, 
     .  charnet, cfa(3,maxua), cfmax, cftem, const, cstress(3,3),
     .  Datm(maxo), Debye, dipol(3), dDmax, dDtol, DEna,
     .  Dold(maxno,maxuo,maxspn), Dscf(maxno,maxuo,maxspn), 
     .  Dscfsave(maxno,maxuo,maxspn), dt, DUext, DUscf, Dxc, 
     .  dx, dxmax, e1, e2, ebk(maxuo,maxspn,maxbk), 
     .  Ecorrec, ef, efs(maxspn),
     .  Eharrs, Eions, Ekin, Ekinion, Emad, Ena, Enaatm, Enascf, 
     .  Enl, Entrop, eo(maxuo,maxspn,maxk), 
     .  Eold(maxno,maxuo,maxspn), Escf(maxno,maxuo,maxspn), 
     .  eta, etol, Etot, eV, Exc, E0,
     .  fa(3,maxua), factor, fdf_convfac,
     .  fmax, fmean, FreeE, fres, ftem, ftol, ftot(3),
     .  g2cut, g2max, getot,
     .  H(maxno,maxuo,maxspn), H0(maxno,maxuo),
     .  kBar, kcutof, kdispl(3), Kelvin,
     .  kn, kpoint(3,maxk), kweight(maxk), kpr, 
     .  mn, mpr, Pint, Pmol, Psol
      double precision
     .  qa(maxa), qaux, qo(maxuo,maxspn,maxk), 
     .  qspin(4), qs(maxspn), qsol, qtot,
     .  rc, rckb(maxkb), rco(maxo),
     .  rcoor, rcoorcp, rcut, reals(2), rijmin, 
     .  rmax, rmaxh, rmaxkb, rmaxo, rmaxv, rmin, r2ij(maxna), r2min,
     .  S(maxno,maxuo), scell(3,3), 
     .  stot, stress(3,3), strtol, svec(3),
     .  taurelax, temp, tempinit, tempion, tiny, tp, ts, tt,
     .  Uatm, ucell(3,3), uion, Uscf,
     .  va(3,maxua), values(2), vcell(3,3), virial, vn,
     .  volcel, volume, vpr, we, wmix, wmixkick, wo,
     .  xa(3,maxa), xij(3,maxna), xijo(3,maxxij), xmax, xmin

      logical
     .  chebef, default, dumpcharge, first, fixspin, found, foundxv, 
     .  gamma, inspn, itest, last, lastst, mix, mmix, negl, noeta, 
     .  outlng, overflow, pulfile, relaxd, same, savehs, savevh, savevt, 
     .  savdrh,  savrho, usesavecg, usesavelwf, usesavedm, usesavexv, 
     .  writbk, writmd, writpx, writb, writec, writedim, writef, 
     .  writei, writek, writic, varcel

      character
     .  filevh*25, filevt*25, fildrh*25, filrho*25,
     .  line*150, names*80, paste*25,
     .  slabel*20, sname*150, shape*10

      external
     .  anneal, chkdim, cgvc, fixed,
     .  dhscf, diagon, dnaefs, extrapol,
     .  fdf_convfac, hsparse, initatom, iodm, ioxv, kgrid, kinefsm,
     .  mulliken, naefs, neighb, nlefsm, nose, npr, 
     .  ordern, overfsm, parse, paste, prmem, pulayx, rcut, redata, 
     .  shaper, spnvec, superc, superx,
     .  timer, uion, verlet2, volcel, xijorb,
     .  lomaxfis, lmxkbfis, nztfl, ioeig, iofa, iokp, iomd


      data
     .  e1, e2 / 1.d0, -1.d0 /
     .  kcutof, kscell, kdispl / 0.d0, 9*0, 3*0.d0 /
     .  nkpnt, (kpoint(i,1),i=1,3), kweight(1) / 1, 3*0.d0, 1.d0 /
     .  relaxd /.false./
     .  tiny /1.d-15/
     .  one  /1/
     .  na, no, nokb, nua, nuo, nuokb, nnamax, nnomax /8*1/
     .  nauxpul, nbk, ns, nspin, nxij /5*1/
     .  ncells, nsc, mscell / 1,   1,1,1,   1,0,0, 0,1,0, 0,0,1 /
     .  overflow / .false. /
c---------------------------------------------------------------------

C Print version information ...........................................
      call prversion
C ..................

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
      kBar   = 1.d0 / 1.47108d5
      Kelvin = eV / 11604.45d0
      Debye  = 0.393430d0
C ..................

C Initialize some variables
      Eharrs = 0.0d0
      Eions = 0.0d0
      FreeE = 0.0d0

C Print array sizes ...................................................
      call prmem( 0, 'siesta', 'auxpul',   'd', maxpul*2           )
      call prmem( 0, 'siesta', 'Dold',     'd', maxno*maxuo*maxspn )
      call prmem( 0, 'siesta', 'Dscf',     'd', maxno*maxuo*maxspn )
      call prmem( 0, 'siesta', 'Dscfsave', 'd', maxno*maxuo*maxspn )
      call prmem( 0, 'siesta', 'Eold',     'd', maxno*maxuo*maxspn )
      call prmem( 0, 'siesta', 'ebk',      'd', maxuo*maxspn*maxbk )
      call prmem( 0, 'siesta', 'eo',       'd', maxuo*maxspn*maxk  )
      call prmem( 0, 'siesta', 'Escf',     'd', maxno*maxuo*maxspn )
      call prmem( 0, 'siesta', 'H',        'd', maxno*maxuo*maxspn )
      call prmem( 0, 'siesta', 'H0',       'd', maxno*maxuo        )
      call prmem( 0, 'siesta', 'listh',    'i', maxno*maxo         )
      call prmem( 0, 'siesta', 'listhold', 'i', maxno*maxo         )
      call prmem( 0, 'siesta', 'S',        'd', maxno*maxuo        )
      call prmem( 0, 'siesta', 'xijo',     'd', 3*maxxij           )
      call prmem( 0, 'siesta', ' ',        ' ', 0                  )
C .....................

C Read simulation data ................................................
      call redata(maxa, maxspn, nua, ns, nspin, overflow, 
     .            slabel, sname, isa, ucell, xa, outlng, g2cut,
     .            charnet, negl, nscf, dDtol, mix, wmix, isolve,
     .            temp, fixspin, ts, ncgmax, ftol, strtol, eta, 
     .            etol, rcoor, 
     .            ioptlwf, chebef, noeta, rcoorcp, beta, pmax,
     .            idyn, istart, ifinal, nmove, ianneal, iquench,
     .            dt, ia1, ia2, dx, dxmax, tt, tp, mn, mpr, 
     .            bulkm, taurelax,
     .            writedim, usesavelwf, usesavedm, usesavecg,
     .            mullipop, inspn, maxsav, nkick, wmixkick, 
     .            pulfile, tempinit, dumpcharge, varcel)
      na = nua
      if (overflow) goto 444
C ................

C Find some switches ..................................................
      writek    = fdf_boolean('WriteKpoints'    , outlng )
      writef    = fdf_boolean('WriteForces'     , outlng )
      writei    = fdf_boolean('WriteEigenvalues', outlng )
      writb     = fdf_boolean('WriteBands'      , outlng )
      writbk    = fdf_boolean('WriteKbands'     , outlng )
      writec    = fdf_boolean('WriteCoorStep'   , outlng )
      writic    = fdf_boolean('WriteCoorInitial', .true. )
      writmd    = fdf_boolean('WriteMDhistory'  , .false.)
      writpx    = fdf_boolean('WriteMDXmol'     , .not. writec)
      default   = fdf_boolean('UseSaveData'     , .false.)
      usesavexv = fdf_boolean('MD.UseSaveXV'    , default)
      savehs    = fdf_boolean('SaveHS'          , .false.)

      rijmin    = fdf_physical( 'WarningMinimumAtomicDistance',
     .                                              1.d0, 'Bohr' )

C .....................

C Read cell shape and atomic positions from a former run ..............
      foundxv = .false.
      if (usesavexv) then
        call ioxv( 'read', ucell, vcell, nua, isa, iza, xa, va, found )
        if (.not.found) write(6,'(/,a)')
     .    'siesta: WARNING: XV file not found'
        foundxv = found
      endif
C ..................

C Dump initial coordinates to output ..................................
      if ( writic ) then
        write(6,'(/a)') 'siesta: Atomic coordinates (Bohr) and species'
        write(6,"('siesta: ',2x,3f10.5,i3,3x,i4)")
     .           ( (xa(ix,ia), ix=1,3), isa(ia), ia, ia=1, na)
      endif
C ..................

C Initialize pseudopotentials, atomic orbitals and atom lists .........
      call initatom(ns, nua, maxuo, maxkb, isa, overflow,
     .              nuo, nuokb, qtot, rmaxv, rmaxo, rmaxkb,
     .              lasto, lastkb, iza, amass,
     .              iaorb, iphorb, Datm, qa, iaKB, iphKB)
      no = nuo
      nokb = nuokb
      qtot = qtot - charnet

C  calculate spin populations for fixed spin case...
      if (fixspin) then
        if (nspin .ne. 2) then
          write(6,'(a)')
     .   'siesta: ERROR: You can only fix the spin of the system'
          write(6,'(a)')
     . '  siesta:        for collinear spin polarized calculations.'
          stop
        endif
        do i=1,2
           qs(i) = (qtot + (3-2*i)*ts) / 2.0d0
cag      qs(1) = (qtot + ts) / 2.0d0
cag      qs(2) = (qtot - ts) / 2.0d0
        enddo
      endif
C ...
       
      if (overflow) goto 444
C ..................

C Find rco and rckb .............................
      do ia = 1,nua
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

C Find maximum interaction range ......................................
      if (negl) then
        rmaxh = 2.d0*rmaxo
      else
        rmaxh = 2.d0*rmaxo + 2.0*rmaxkb
      endif
C ......................

C Automatic cell generation ...........................................
      if (volcel(ucell) .lt. 1.d-8) then
        do iv = 1,3
          do ix = 1,3
            ucell(ix,iv) = 0.d0
            scell(ix,iv) = 0.d0
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
          ucell(ix,ix) = 1.10d0 * (xmax - xmin)
          scell(ix,ix) = ucell(ix,ix)
        enddo
C build cubic cell if system is charged
        if (charnet .ne. 0.0d0) then
          xmax = -1.d30
          do ix = 1,3
            if (ucell(ix,ix) .gt. xmax) xmax = ucell(ix,ix)
          enddo
          do ix = 1,3
            ucell(ix,ix) = xmax
            scell(ix,ix) = xmax
          enddo
        endif
C
        volume = volcel( ucell )
        write(6,'(/,a,3(/,a,3f12.6))')
     .    'siesta: Automatic unit cell vectors (Ang):',
     .    ('siesta:', (ucell(ix,iv)/Ang,ix=1,3), iv =1,3)
      endif
C ......................

C Find system shape ...................................................
      call shaper( ucell, nua, isa, xa, shape, nbcell, bcell )
      write(6,'(/,2a)') 'siesta: System type = ', shape
C ......................

C Madelung correction for charged systems .............................
      if (charnet .ne. 0.0d0) then
        call madelung(ucell, shape, charnet, Emad)
      endif
C ......................

C Check dimensions for routine matel ..................................
      call matel0( ns )
C ......................

C Find k-grid for Brillouin zone integration ..........................
      nkpnt = maxk
      call kgrid( ucell, kscell, kdispl,
     .            kcutof, nkpnt, kpoint, kweight )
      if (nkpnt .gt. maxk) overflow = .true.
C ......................

C Find number of band k-points ........................................
      nbk = 0
      call bands( no, nspin, maxspn, maxuo, maxno, 0,
     .            numh, listh, H, S, Ef, xijo, indxuo,
     .            .false., nbk, bk, ebk )
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
      if (nbk .gt. 0) gamma = .false.
C ....................

C Print k-points ......................................................
      if (.not.gamma .and. nkpnt.le.maxk) then
        if ( writek ) then
          write(6,'(/,a)')
     .     'siesta: k-point coordinates (Bohr**-1) and weights:'
          write(6,'(a,i4,3f12.6,3x,f12.6)')
     .      ('siesta: ', ik, (kpoint(ix,ik),ix=1,3), kweight(ik),
     .      ik=1,nkpnt)
        else
          call iokp( nkpnt, kpoint, kweight )
        endif
        write(6,'(/a,i6)')
     .    'siesta: k-grid: Number of k-points =', nkpnt
        write(6,'(a,f10.3,a)')
     .    'siesta: k-grid: Cutoff             =', kcutof/Ang, ' Ang'
        write(6,'(a)')
     .    'siesta: k-grid: Supercell and displacements'
        write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .                             (kscell(i,1),i=1,3), kdispl(1)
        write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .                             (kscell(i,2),i=1,3), kdispl(2)
        write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .                             (kscell(i,3),i=1,3), kdispl(3)
      endif
C ....................

C Find auxiliary supercell (required only for k sampling) ............
      call superc( gamma, rmaxh, ucell, nua, nuo, nuokb,
     .             maxa, maxo, maxkb,
     .             overflow, xa, isa, iza, lasto, lastkb,
     .             iaorb, iphorb, rco, iakb, iphkb, rckb,
     .             scell, mscell, nsc, ncells,
     .             na, no, nokb, indxua, indxuo )
      if (.not.gamma) write(6,'(/,a,i6,a,i6,a,i6,a,i8)')
     .  'siesta: Internal auxiliary supercell:',
     .   nsc(1), '  x', nsc(2), '  x', nsc(3), '  =', ncells
      if (overflow) goto 444
C ..................

C Some printout for debugging .......................
*     write(6,'(/,a,/,(3f12.6))') 'siesta: scell =', scell
*     write(6,'(/,a,/,(3i6,2i8,3f12.6))')
*    .  'siesta: ia,isa,iza,lasto,lastkb,xa =',
*    .  (ia,isa(ia),iza(ia),lasto(ia),lastkb(ia),
*    .   (xa(ix,ia),ix=1,3),ia=1,na)
*     write(6,'(/,a,3i6)') 'siesta: na, no, nokb =', na, no, nokb
*     write(6,'(/,a)') 'siesta: indxua, indxuo ='
*     do ia = 1,na
*       write(6,'(i6,3x,20i4)')
*    .    indxua(ia), (indxuo(io),io=lasto(ia-1)+1,lasto(ia))
*     enddo
C ..................

C Initialize atomic velocities to zero ................................
      if (.not. foundxv) then
        do ia = 1,nua
          do ix = 1,3
            va(ix,ia) = 0.0
          enddo
        enddo
        do i = 1,3
          do j = 1,3
            vcell(i,j)=0.0
          enddo
        enddo
      endif
C ..................

C Begin of coordinate relaxation iteration ============================
C Notice that this loop is not indented
      if (idyn .eq. 0) then
        inicoor = 0
        fincoor = nmove
      else if (idyn .ge. 1 .and. idyn .le. 5) then
        inicoor = istart
        fincoor = ifinal
      else if (idyn .eq. 6) then
        inicoor = 0
        fincoor = (ia2-ia1+1)*3*2
      else
        stop 'siesta: wrong idyn'
      endif

C Build initial velocities according to Maxwell-Bolzmann distribution....
      if (idyn .ne. 0 .and. idyn .ne. 6 .and. (.not. foundxv)) 
     .               call vmb(nua,tempinit,amass,va)
C ..................

      istp = 0
      do istep = inicoor,fincoor
      call timer( 'IterMD', 1 )
      istp = istp + 1
      write(6,'(/,a)') 'siesta:    ==============================='
      if (idyn .eq. 0) 
     . write(6,'(a,i6)') 'siesta:        Begin CG move = ',istep
      if (idyn .ne. 0) 
     . write(6,'(a,i6)') 'siesta:        Begin MD step = ',istep
      if (idyn .eq. 6)  then
       write(6,'(a,i6)') 'siesta:        Begin FC step = ',istep
       if (istep .eq. 0) then
         write(6,'(a)') 'siesta:        Undisplaced coordinates'
       else
         iadispl = (istep-mod(istep-1,6))/6+ia1
         write(6,'(a,i6)') 'siesta:        displace atom = ', iadispl
         ix = mod(istep-1,6)+1
         ixdispl = (ix - mod(ix-1,2) +1)/2
         write(6,'(a,i6)') 'siesta:        in direction  = ', ixdispl
         dx=-dx
         write(6,'(a,f8.4,a)') 'siesta:        by            = ',
     .                      dx, ' Bohr'
c displace atom by dx...
         xa(ixdispl,iadispl)=xa(ixdispl,iadispl)+dx
       endif
      endif
      write(6,'(a)')   'siesta:    ==============================='

C Print atomic coordinates ............................................
      call outcoor( ucell, xa, isa, nua, ' ', writec )
C ...................


C Actualize things if variable cell ...................................

      if ( varcel .and. (istep.ne.inicoor) .and. (.not.gamma) ) then

C k-grid 
        nkpnt = maxk
        call kgrid( ucell, kscell, kdispl,
     .              kcutof, nkpnt, kpoint, kweight )
        if (nkpnt .gt. maxk) overflow = .true.
 
C Print k-points 
        if (nkpnt.le.maxk) then
          if ( writek ) then
            write(6,'(/,a)')
     .       'siesta: k-point coordinates (Bohr**-1) and weights:'
            write(6,'(a,i4,3f12.6,3x,f12.6)')
     .        ('siesta: ', ik, (kpoint(ix,ik),ix=1,3), kweight(ik),
     .        ik=1,nkpnt)
          else
            call iokp( nkpnt, kpoint, kweight )
          endif
          write(6,'(/a,i6)')
     .      'siesta: k-grid: Number of k-points =', nkpnt
          write(6,'(a,f10.3,a)')
     .      'siesta: k-grid: Cutoff             =', kcutof/Ang, ' Ang'
          write(6,'(a)')
     .      'siesta: k-grid: Supercell and displacements'
          write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .                               (kscell(i,1),i=1,3), kdispl(1)
          write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .                               (kscell(i,2),i=1,3), kdispl(2)
          write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .                               (kscell(i,3),i=1,3), kdispl(3)
        endif
 
C auxiliary supercell
        call superc( gamma, rmaxh, ucell, nua, nuo, nuokb,
     .               maxa, maxo, maxkb,
     .               overflow, xa, isa, iza, lasto, lastkb,
     .               iaorb, iphorb, rco, iakb, iphkb, rckb,
     .               scell, mscell, nsc, ncells,
     .               na, no, nokb, indxua, indxuo )
        if (.not.gamma) write(6,'(/,a,i6,a,i6,a,i6,a,i8)')
     .    'siesta: Internal auxiliary supercell:',
     .     nsc(1), '  x', nsc(2), '  x', nsc(3), '  =', ncells
        if (overflow) goto 444

C Madelung correction for charged systems .............................
        if (charnet .ne. 0.0d0) then
          call madelung(ucell, shape, charnet, Emad)
        endif

      endif
C End variable cell actualization


C Print unit cell for variable cell ...................................
      if ( varcel ) call outcell(ucell)
C ...................

C Find unit cell volume ...............................................
      volume = volcel( ucell )
C ...................

C Initialize neighb subroutine ........................................
      ia = 0
      isel = 0
      rmax = max( 2.d0*rmaxv, 2.d0*rmaxo, rmaxo+rmaxkb )
      nnia = maxna
      call neighb( scell, rmax, na, xa, ia, isel,
     .             nnia, jna, xij, r2ij )
      nnamax = 0
      do ia = 1,na
        nnia = 0
        call neighb( scell, rmax, na, xa, ia, isel,
     .               nnia, jna, xij, r2ij )
        nnamax = max( nnamax, nnia )
      enddo
      if (nnamax .gt. maxna) then
C       Increase nnamax with safety margin when atoms move
        nnamax = nnamax + 0.10 * nnamax + 10
        overflow = .true.
      else
        nnamax = maxna
      endif
      if (overflow) goto 444
C ..................

C Check if any two atoms are unreasonably close .......................
      do ia = 1,na
        r2min = 1.d30
        jamin = 0
        nnia = maxna
        call neighb( scell, rmax, na, xa, ia, isel,
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
      nnomax = maxno
      call hsparse( negl, scell, nsc, na, isa, xa,
     .              lasto, lastkb, iphorb, iphKB,
     .              nnomax, numh, listh )
      if (nnomax .gt. maxno) then
C       Increase nnomax with safety margin when atoms move
        nnomax = nnomax + 0.10 * nnomax + 40
        overflow = .true.
      else
        nnomax = maxno
      endif
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
      if (.not.gamma) then
        nxij = nnomax * nuo
        if (nxij .le. maxxij) then
          call xijorb( negl, scell, nua, na, xa,
     .                 lasto, lastkb, rco, rckb,
     .                 maxno, numh, listh, xijo )
        else
          overflow = .true.
        endif
      else
        nxij = 1
      endif
      if (overflow) goto 444
C ..................

C Initialize density matrix ...........................................
C set density matrix for first step
      if (istp .eq. 1) 
     .   call initdm(Datm, Dscf, Dold, Escf, lasto, maxa,
     .               maxno, maxuo, maxspn, nua, nuo, nspin,
     .               numh, numhold, listh, listhold, iaorb,
     .               found, inspn, usesavedm)

      if (istp .ne. 1 .and. idyn .eq. 6)
     .   call initdm(Datm, Dscf, Dold, Escf, lasto, maxa,
     .               maxno, maxuo, maxspn, nua, nuo, nspin,
     .               numh, numhold, listh, listhold, iaorb,
     .               found, inspn, usesavedm)

C Extrapolate density matrix between steps
      itest = .false.
      istpsave = 0
      iord = 1
      if (idyn .eq. 0) iord = 0
      if (idyn .eq. 6) iord = 0
C  If DM has just been read from disk, 
C  call extrapol with istep = 2 and iord = 0
C  to make it update the structure of DM, if needed
      if (found .and. ((istp .eq. 1) .or. (idyn .eq. 6))) then
        istpsave = istp
        istp = 2
        iord = 0
        itest = .true.
      endif
      call extrapol(istp, iord, nspin, nuo, maxuo, maxno,  numh, listh,
     .              aux, numhold, listhold, Dscfsave, Dscf)
C  If DM have just been read, restore istp
      if (itest) istp = istpsave
      itest = .false.
C ..................

C Check for Pulay auxiliary matrices sizes ...................................
111   continue
      if (pulfile .or. maxsav .le. 0) then
        if (maxpul .ne. one)  then
          nauxpul = 1
          overflow = .true.
          goto 444
        endif
      else
        nauxpul = 0
        do io = 1,nuo
          nauxpul = nauxpul + numh(io)
        enddo
        nauxpul = nauxpul * nspin * maxsav
        call chkdime(maxpul,nauxpul,overflow,nauxpul)
C       Increase nauxpul with safety margin when atoms move
        nauxpul = nauxpul + 0.1 * nauxpul + 10
        if (overflow) goto 444
      endif
      if (writedim) goto 444
C ....................

C Start of SCF iteration _____________________________________________
      first = .true.
      last  = .false.
      if (wmix .le. 0.d0) then
        write(6,'(/,a,f15.8)')
     .   'siesta: WARNING: Mixing weight for SCF loop =', wmix
        last = .true.
      endif

      do iscf = 1, nscf
        if (iscf .eq. nscf) last = .true.
        call timer( 'IterSCF', 1 )

C Find overlap matrix ...............................................
        if (first) then
           call overfsm(nua, na, no, scell, xa, indxua, rmaxo, maxuo,
     .                 maxna, maxno, maxno, lasto, iphorb, isa, 
     .                 numh, listh, numh, listh, min(nspin,2), Escf, 
     .                 jna, xij, r2ij,
     .                 fa, stress, S)
        endif
C ..................

C Normalize density matrix to exact charge ...........................
        qsol = 0.d0
        do ispin = 1,min(nspin,2)
          do io = 1,nuo
            do in = 1,numh(io)
              qsol = qsol + Dscf(in,io,ispin) * s(in,io)
            enddo
          enddo
        enddo
        if (.not.first .and.
     .       abs(qsol/qtot-1.d0).gt.1.d-2) write(6,'(a,2f15.6)')
     .      'siesta: WARNING: Qtot, Tr[D*S] =', qtot, qsol
        do ispin = 1,nspin
          do io = 1,nuo
            do in = 1,numh(io)
              Dscf(in,io,ispin) = Dscf(in,io,ispin) * qtot/qsol
              Escf(in,io,ispin) = Escf(in,io,ispin) * qtot/qsol
            enddo
          enddo
        enddo
C ..................

C Initialize Hamiltonian ........................................
        do ispin = 1,nspin
          do io = 1,nuo
            do j = 1,numh(io)
              H(j,io,ispin) = 0.0d0
            enddo
          enddo
        enddo
C ..................

C Initialize forces and stress ...................
        if (first.or.last) then
          do ia = 1,nua
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
          do ia = 1,nua
            is = isa(ia)
            Eions = Eions + uion(is)
          enddo
        endif
C ..................

C Neutral-atom: energy, forces and stress ............................
C First time for energy, last time for forces
        if (first.or.last) then
          call naefs(nua, na, scell, xa, indxua, rmaxv,
     .               maxna, isa, jna, xij, r2ij,
     .               Ena, fa, stress)
          call dnaefs(nua, na, scell, xa, indxua, rmaxv,
     .               maxna, isa, jna, xij, r2ij,
     .               DEna, fa, stress) 
          Ena = Ena + DEna
        endif
C ..................

C Kinetic: energy, forces, stress and matrix elements .................
        if (first.or.last) then
          call kinefsm(nua, na, no, scell, xa, indxua, rmaxo, maxuo,
     .                 maxna, maxno, maxno, lasto, iphorb, isa, 
     .                 numh, listh, numh, listh, min(nspin,2), Dscf, 
     .                 jna, xij, r2ij,
     .                 Ekin, fa, stress, H) 
        endif
C ..................

C Non-local-pseudop: energy, forces, stress and matrix elements .......
        if (first.or.last) then
          call nlefsm(scell, nua, na, isa, xa, indxua,
     .                maxuo, maxno, maxno,
     .                lasto, lastkb, iphorb, iphKB,
     .                numh, listh, numh, listh, min(nspin,2), Dscf, 
     .                Enl, fa, stress, H) 
        endif
C ..................

C Save or get partial Hamiltonian (non-SCF part) ......................
        if (first.or.last) then
          do io = 1,nuo
            do j = 1,numh(io)
              H0(j,io) = H(j,io,1)
            enddo
          enddo
        else
          do ispin = 1,nspin
            do io = 1,nuo
              do j = 1,numh(io)
                if (ispin .le. 2) then
                  H(j,io,ispin) = H0(j,io)
                else
                  H(j,io,ispin) = 0.d0
                endif
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
          do ispin = 1,min(nspin,2)
            do io = 1,nuo
              do j = 1,numh(io)
                E0 = E0 + H0(j,io) * Dscf(j,io,ispin)
              enddo
            enddo
          enddo
        endif
C ..................

C Add SCF contribution to energy and matrix elements ..................
        ilv = 0
        g2max = g2cut
        if (last) then
c         Last call to dhscf and grid-cell sampling if requested
          ifa  = 1
          istr = 1
          call grdsam( nspin, maxuo, no, iaorb, iphorb, indxuo,
     .                 nua, na, isa, xa, indxua,
     .                 ucell, mscell, g2max, ntm,
     .                 ilv, ifa, istr,
     .                 maxno, numh, listh, Dscf, Datm, H,
     .                 Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .                 Exc, Dxc, dipol, fa, stress, ierror ) 
        else
          ifa  = 0
          istr = 0
          ihmat = 1
          call dhscf( nspin, maxuo, no, iaorb, iphorb, indxuo,
     .                na, isa, xa, indxua, ucell, mscell, g2max, ntm,
     .                ilv, ifa, istr, ihmat, ' ', ' ', ' ', ' ',
     .                maxno, numh, listh, Dscf, Datm,
     .                maxno, numh, listh, H,
     .                Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .                Exc, Dxc, dipol, fa, stress, ierror )
        endif
            
*       if (istp.eq.1 .and. iscf.eq.1) write(6,'(/,a,f10.3,a)')
*    .    'siesta: dhscf mesh cutoff =', g2max, ' Ry'

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
           call overfsm(nua, na, no, scell, xa, indxua, rmaxo, maxuo,
     .                 maxna, maxno, maxno, lasto, iphorb, isa, 
     .                 numh, listh, numh, listh, min(nspin,2), Escf, 
     .                 jna, xij, r2ij,
     .                 fa, stress, S) 
        endif
C ..................

C Find entropy ........................................................
        Entrop = 0.d0
        if (isolve .eq. 0) then
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
          endif
        endif
C ..................

C Save present density matrix ........................................
        do ispin = 1,nspin
          do io = 1,nuo
            do jo = 1,numh(io)
              Dold(jo,io,ispin) = Dscf(jo,io,ispin)
              Eold(jo,io,ispin) = Escf(jo,io,ispin)
            enddo
          enddo
        enddo
C ..................

C Save hamiltonian and overlap matrices ............................
        if (savehs) then
          call iohs( 'write', gamma, nuo, nspin, maxuo, maxno,
     .               numh, listh, H, S, qtot, temp, xijo )
        endif
C ..................

C Solve eigenvalue problem .........................................
        if (isolve .eq. 0) then
          call diagon(no, nspin, maxspn, maxuo, maxno, maxno,
     .                numh, listh, numh, listh, H, S,
     .                qtot, fixspin, qs, temp, e1, e2,
     .                gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     .                eo, qo, Dscf, Escf, ef, efs)
          Ecorrec = 0.d0
        elseif (isolve .eq. 1) then
          call ordern(usesavelwf,ioptlwf,na,no,lasto,
     .                isa,qa,rcoor,rmaxh,ucell,xa,iscf,istp,ncgmax,
     .                etol,eta,qtot,maxno,numh,listh,h,s,
     .                chebef, noeta, rcoorcp, beta, pmax,
     .                Dscf,Escf,Ecorrec)
        else
          write(6,'(/,a,I3,A)') 'siesta: ERROR: wrong solution method'
          stop
        endif
C ..................

C Harris-functional energy ............................................
        Eharrs = 0.0d0
        do ispin = 1,nspin
C         const factor takes into account that there are two nondiagonal
C         elements in non-collinear spin density matrix, stored as
C         ispin=1 => D11; ispin=2 => D22, ispin=3 => Real(D12);
C         ispin=4 => Imag(D12)
          const = 1.d0
          if (ispin .gt. 2) const = 2.d0
          do io = 1,nuo
            do j = 1,numh(io)
              Eharrs = Eharrs + H(j,io,ispin) * const * 
     .                         ( Dscf(j,io,ispin) - Dold(j,io,ispin) )
            enddo
          enddo
        enddo
C ..................

C Mix input and output energy-density and density matrices ............
C Following line for using and saving the density matrix without mix ..
        if (wmix.ne.0.d0) then
C Pulay mixing
          mmix  = mix
          iiscf = iscf
          if (maxsav .le. 0) then
            iiscf = 1
            if (iscf .ne. 1) mmix = .true.
          endif
          call pulayx( pulfile, iiscf, mmix, nuo, maxuo, maxno, numh,
     .                 nspin, maxsav, wmix, nkick, wmixkick,
     .                 auxpul(1,1), auxpul(1,2), maxpul,
     .                 Dscf, Dold, dDmax)
        endif
C ...................

C Save density matrix on disk, after mixing ...........................
        if (idyn .eq. 6) then
          if (istp .eq.1)
     .    call iodm( 'write', maxno, maxuo, nuo, nspin,
     .               numh, listh, Dscf, found )
        else
          call iodm( 'write', maxno, maxuo, nuo, nspin,
     .               numh, listh, Dscf, found )
        endif
C ...................

C Print energies ......................................................
        DEna = Enascf - Enaatm
        Etot = E0 + DEna + DUscf + DUext + Exc + Ecorrec + Emad
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
     .   'siesta: Emadel  =', Emad/eV,
     .   'siesta: Eharris =', Eharrs/eV,
     .   'siesta: Etot    =', Etot/eV,
     .   'siesta: FreeEng =', FreeE/eV
C ...................

C Print total energy and density matrix error .........................
        if (isolve .eq. 0) then
          if (fixspin) then
            if (iscf .eq. 1)  write(6,'(/,a12,3a14,a8,a7,a11)')
     .      'siesta: iscf', '   Eharris(eV)', 
     .      '      E_KS(eV)', '   FreeEng(eV)', 
     .      '   dDmax', '  Ef_up', '  Ef_dn(eV)'
            write(6,'(a8,i4,3f14.4,f8.4,2f9.4)')
     .      'siesta: ',iscf, Eharrs/eV, Etot/eV, FreeE/eV, dDmax, 
     .                 (Efs(i)/eV,i=1,2)

          else
            if (iscf .eq. 1) write(6,'(/,a12,3a14,2a8)')
     .      'siesta: iscf', '   Eharris(eV)', 
     .      '      E_KS(eV)', '   FreeEng(eV)', 
     .      '   dDmax', '  Ef(eV)'
            write(6,'(a8,i4,3f14.4,2f8.4)')
     .      'siesta: ',iscf, Eharrs/eV, Etot/eV, FreeE/eV, dDmax, Ef/eV
          endif
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
            do io = 1,nuo
              do jo = 1,numh(io)
                Dscf(jo,io,ispin) = Dold(jo,io,ispin)
                Escf(jo,io,ispin) = Eold(jo,io,ispin)
              enddo
            enddo
          enddo
          if (dumpcharge) then
            call plcharge( nspin, nuo, no, na, maxo, maxa, maxno,
     .                     scell, rmaxo, xa, 
     .                     lasto, isa, iphorb, Datm,
     .                     numh, listh, indxuo )
          endif
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

C Impose constraints to atomic movements by changing forces ...........
      call fixed( ucell, stress, nua, isa, amass, xa, fa, cstress, cfa )
C ...................

C Write atomic forces .................................................
      fmax = 0.d0
      cfmax = 0.d0
      fres = 0.d0
      do ix = 1,3
        ftot(ix) = 0.d0
        do ia = 1,nua
          ftem = fa(ix,ia)
          cftem = cfa(ix,ia)
          ftot(ix) = ftot(ix) + ftem
          fres = fres + ftem*ftem
          fmax = max( fmax, dabs(ftem) )
          cfmax = max( cfmax, dabs(cftem) )
        enddo
      enddo
      fres = dsqrt( fres / (3.d0*nua) )
      write(6,'(/,a)') 'siesta: Atomic forces (eV/Ang):'
      if ( writef ) then
        write(6,'(i4,3f12.6)') (ia,(fa(ix,ia)*Ang/eV,ix=1,3),ia=1,nua)
      else
        call iofa( nua, fa )
      endif
      write(6,'(40(1h-),/,a4,3f12.6)') 'Tot',(ftot(ix)*Ang/eV,ix=1,3)
      write(6,'(40(1h-),/,a4, f12.6)') 'Max',fmax*Ang/eV
      write(6,'(a4,f12.6,a)')'Res',fres*Ang/eV,
     .                   '    sqrt( Sum f_i^2 / 3N )'
      write(6,'(40(1h-),/,a4, f12.6,a)') 'Max',cfmax*Ang/eV, 
     .                   '    constrained'

C Write Force Constant matrix if FC calculation ...

      if (idyn .eq. 6) then
        call ofc(fa,dx,nua)
      endif

C Write stress tensor for any variable cell dynamics ..................
      if ( varcel ) then
         write(6,'(/,a,3(/,a,3f12.6))')
     .     'siesta: Stress tensor (eV/Ang**3):',
     .     ('     ',(stress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)

C       Pressure (only for the solid)
        Psol = - (( stress(1,1) + stress(2,2) + stress(3,3) ) / 3.0)
        write(6,'(/,a,f20.8,a)')
     .    'siesta: Pressure:', Psol/kBar, '  kBar'

      endif
C ...................

C Mulliken population analysis .......................................
      call mulliken( mullipop, nspin, nua, nuo, maxuo, maxno,
     .               numh, listh, S, Dscf, isa, lasto, iaorb, iphorb )

C ...................

C Move atoms ..........................................................
      if (idyn .eq. 0) then
        if (nmove .ne. 0) then
          call cgvc( nua, xa, cfa, ucell, cstress, volume, dxmax,
     .               tp, ftol, strtol, varcel, relaxd, usesavecg )
          if (relaxd) then
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
        call verlet2(istp, iunit, iquench, nua, cfa, dt, 
     .               amass, va, xa, Ekinion, tempion)
      
      elseif (idyn .eq. 2) then
        call nose(istp, iunit, nua, cfa, tt, dt, amass, mn, 
     .            va, xa, Ekinion, kn, vn, tempion)

      elseif (idyn .eq. 3) then
        call pr(istp, iunit, iquench, nua, cfa, cstress, tp, dt, amass,
     .          mpr, va, xa, vcell, ucell, Ekinion, kpr, vpr, tempion,
     .          Pint)
        write(6,'(/,a,f12.3,a)')
     .    'siesta: E_kin PR =', kpr/Kelvin, ' K'
      
      elseif (idyn .eq. 4) then
        call npr(istp, iunit, nua, cfa, cstress, tp, tt, dt, amass, mn,
     .           mpr, va, xa, vcell, ucell, Ekinion, kn, kpr, vn, vpr,
     .           tempion, Pint)

      elseif (idyn .eq. 5) then
        call anneal(istp, iunit, ianneal, taurelax, bulkm,
     .              nua, cfa, cstress, tp, tt, dt, amass,
     .              va, xa, ucell, Ekinion, tempion, Pint)
      endif

      if (idyn .gt. 0) then
        write(6,'(/,a,f12.3,a)')
     .    'siesta: Temp_ion =', tempion, ' K'
      endif
C ...................

C Change supercell and expand atomic positions to supercell ..........
      call superx( ucell, nsc, nua, maxa, xa, scell )
C ...................

C Save last atomic positions and velocities ..........................
      call ioxv( 'write', ucell, vcell, nua, isa, iza, xa, va, found )
C ...................

      getot = Etot + Ekinion + kn + kpr + vn + vpr

c restore original coordinates after FC displacements ...
      if (idyn .eq. 6 .and. istep .ne. 0) then
        xa(ixdispl,iadispl)=xa(ixdispl,iadispl)-dx
      endif
c ...

C Save atomic positions and velocities accumulatively ................
      if (writmd) call iomd( nua, isa, iza, xa, va, ucell, vcell, 
     .         varcel, istep, inicoor, fincoor, tempion, Etot, getot)

C Accumulate coor in Xmol file for animation .........................
      lastst = fincoor .le. istep
      if (writpx) call pixmol(iza, xa, nua, slabel, lastst)
C ...................

      call timer( 'IterMD', 2 )
      enddo
   60 continue
C End of coordinate-relaxation loop ==================================

C Print atomic coordinates (and also unit cell for ParrRah.)
      if (nmove .ne. 0) then
        if (relaxd) 
     .    call outcoor(ucell, xa, isa, nua,'Relaxed', .true. )
        if (.not.relaxd) 
     .    call outcoor(ucell, xa, isa, nua,
     .                 'Final (unrelaxed)', .true. )
        if ( varcel ) call outcell(ucell)
      endif


C Print coordinates in xmol format in a separate file

      if(fdf_boolean('WriteCoorXmol',.false.)) 
     .   call coxmol(iza, xa, nua, slabel)

C Print coordinates in cerius format in a separate file

      if(fdf_boolean('WriteCoorCerius',.false.))
     .   call coceri(iza, xa, ucell, nua, sname, slabel)


c Find and print bands
      nbk = 0
      call bands( no, nspin, maxspn, maxuo, maxno, maxbk,
     .            numh, listh, H, S, Ef, xijo, indxuo,
     .            .true., nbk, bk, ebk )
      if (nbk.gt.0 .and. nbk.le.maxbk ) then
        if ( writbk ) then
          write(6,'(/,a,/,a4,a12)')
     .     'siesta: Band k vectors (Bohr**-1):', 'ik', 'k'
          do ik = 1,nbk
            write(6,'(i4,3f12.6)') ik, (bk(ix,ik),ix=1,3)
          enddo
        endif
        
        if ( writb ) then
          write(6,'(/,a,/,a4,a3,a7)')
     .     'siesta: Band energies (eV):', 'ik', 'is', 'eps'
          do ispin = 1,min(nspin,2)
            do ik = 1,nbk
              write(6,'(i4,i3,10f7.2)')
     .          ik, ispin, (ebk(io,ispin,ik)/eV,io=1,min(10,nuo))
              if (nuo.gt.10) write(6,'(7x,10f7.2)')
     .            (ebk(io,ispin,ik)/eV,io=11,nuo)
            enddo
          enddo
        endif
      endif


C Print eigenvalues
      if (isolve.eq.0 .and. writei) then
        if (nspin .le. 2) then
          write(6,'(/,a,/,a4,a3,a7)')
     .     'siesta: Eigenvalues (eV):', 'ik', 'is', 'eps'
          do ik = 1,nkpnt
            do ispin = 1,nspin
              write(6,'(i4,i3,10f7.2)')
     .          ik, ispin, (eo(io,ispin,ik)/eV,io=1,min(10,nuo))
              if (nuo.gt.10) write(6,'(7x,10f7.2)')
     .          (eo(io,ispin,ik)/eV,io=11,nuo)
            enddo
          enddo
        else
          write(6,'(/,a)') 'siesta: Eigenvalues (eV):'
          do ik = 1,nkpnt
            write(6,'(a,i6)') 'ik =', ik
            write(6,'(10f7.2)')
     .        ((eo(io,ispin,ik)/eV,io=1,nuo),ispin=1,2)
          enddo
        endif
        write(6,'(a,f15.6,a)') 'siesta: Fermi energy =', ef/eV, ' eV'
      endif

      if (isolve.eq.0) 
     .     call ioeig(eo, ef, nuo, nspin, nkpnt, maxuo, maxspn, maxk)

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
     .   'siesta: Emadel  =', Emad/eV,
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
        do ia = 1,nua
          fmax = max( fmax, dabs(fa(ix,ia)) )
          ftot(ix) = ftot(ix) + fa(ix,ia)
        enddo
      enddo
      if (fmax .gt. ftol) then
        write(6,'(/,a)') 'siesta: Atomic forces (eV/Ang):'
        write(6,'(a,i4,3f12.6)')
     .    ('siesta: ', ia,(fa(ix,ia)*Ang/eV,ix=1,3),ia=1,nua)
        write(6,'(a,40(1h-),/,a,a4,3f12.6)')
     .    'siesta: ','siesta: ','Tot',(ftot(ix)*Ang/eV,ix=1,3)
      endif

C Print constrained atomic forces
      same = .true.
      do ia = 1,nua
        do ix = 1,3
          if (cfa(ix,ia) .ne. fa(ix,ia)) same = .false.
        enddo
      enddo
      if (.not.same) then
        fmax = 0.d0
        do ix = 1,3
          ftot(ix) = 0.d0
          do ia = 1,nua
            fmax = max( fmax, dabs(cfa(ix,ia)) )
            ftot(ix) = ftot(ix) + cfa(ix,ia)
          enddo
        enddo
        if (fmax .gt. ftol) then
          write(6,'(/,a)') 'siesta: Constrained forces (eV/Ang):'
          write(6,'(a,i4,3f12.6)')
     .      ('siesta: ', ia,(cfa(ix,ia)*Ang/eV,ix=1,3),ia=1,nua)
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
        do ia = 1,nua
          fmean = fmean + fa(ix,ia) / nua
        enddo
        do ia = 1,nua
          virial = virial + xa(ix,ia) * (fa(ix,ia) - fmean)
        enddo
      enddo
      Psol = - (( stress(1,1) + stress(2,2) + stress(3,3) ) / 3.0)
      Pmol = Psol - virial / volume / 3.0
      write(6,'(/,a,f18.6,a)')
     .  'siesta: Cell volume =', volume/Ang**3, ' Ang**3'
      write(6,'(/,a,/,a,2a20,a,3(/,a,2f20.8,a))')
     .  'siesta: Pressure:',
     .  'siesta: ','Solid',        'Molecule',      '  Units',
     .  'siesta: ', Psol,           Pmol,           '  Ry/Bohr**3',
     .  'siesta: ', Psol*Ang**3/eV, Pmol*Ang**3/eV, '  eV/Ang**3',
     .  'siesta: ', Psol/kBar,      Pmol/kBar,      '  kBar'

c Print spin polarization
      if (nspin .ge. 2) then
        do ispin = 1,nspin
          qspin(ispin) = 0.d0
          do io = 1,nuo
            do j = 1,numh(io)
              jo = listh(j,io)
              qspin(ispin) = qspin(ispin) +
     .                       Dscf(j,io,ispin) * S(j,io)
            enddo
          enddo
        enddo
        if (nspin .eq. 2) then
          write(6,'(/,a,f12.6)')
     .     'siesta: Total spin polarization (Qup-Qdown) =', 
     .     qspin(1) - qspin(2)
        elseif (nspin .eq. 4) then
          call spnvec( nspin, qspin, qaux, stot, svec )
          write(6,'(/,a,f12.6)')
     .     'siesta: Total spin polarization (Qup-Qdown) =', stot
          write(6,'(a,3f12.6)') 'siesta: Spin vector =', svec
        endif
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
      savdrh = fdf_boolean( 'SaveDeltaRho',               .false. )
      savevh = fdf_boolean( 'SaveElectrostaticPotential', .false. )
      savevt = fdf_boolean( 'SaveTotalPotential',         .false. )
      if (savrho .or. savdrh .or. savevh .or. savevt) then
        filrho = ' '
        fildrh = ' '
        filevh = ' '
        filevt = ' '
        if (savrho) filrho = paste( slabel, '.RHO' )
        if (savdrh) fildrh = paste( slabel, '.DRHO' )
        if (savevh) filevh = paste( slabel, '.VH'  )
        if (savevt) filevt = paste( slabel, '.VT'  )
        g2max = g2cut
        call dhscf( nspin, maxuo, no, iaorb, iphorb, indxuo,
     .              na, isa, xa, indxua, ucell, mscell, g2max, ntm,
     .              ilv, 0, 0, 0, filrho, fildrh, filevh, filevt,
     .              maxno, numh, listh, Dscf, Datm,
     .              maxno, numh, listh, H,
     .              Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, Exc, Dxc,
     .              dipol, fa, stress, ierror )
      endif

c Find local density of states
      if ( fdf_block('LocalDensityOfStates',iu) ) then

c       Find the desired energy range
        read(iu,'(a)') line
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
          call diagon( no, nspin, maxspn, maxuo, maxno, maxno,
     .                 numh, listh, numh, listh, H, S, 
     .                 qtot, fixspin, qs, temp, e1, e2,
     .                 gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     .                 eo, qo, Dscf, Escf, ef, efs )
        else
          write(6,*)
     .     'siesta: ERROR: LDOS implemented only with diagon'
          goto 70
        endif
     
c       Find the LDOS in the real space mesh
        filrho = paste( slabel, '.LDOS' )
        g2max = g2cut
        call dhscf( nspin, maxuo, no, iaorb, iphorb, indxuo,
     .              na, isa, xa, indxua, ucell, mscell, g2max, ntm,
     .              ilv, 0, 0, 0, filrho, ' ', ' ', ' ',
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
     .'c  maxa   : Maximum total number of atoms'
      write(1,'(a)')
     .'c  maxbk : Maximum number of band k points'
      write(1,'(a)')
     .'c  maxk   : Maximum number of k points'
      write(1,'(a)') 
     .'c  maxkb  : Maximum total number of Kleinman-Bylander projectors'
      write(1,'(a)') 
     .'c  maxna  : Maximum number of neighbours of any atom'
      write(1,'(a,/,a)') 
     .'c  maxno  : Max basis orbitals connected to any given one,',
     .'c           either directly or through a KB projector'
      write(1,'(a)') 
     .'c  maxo   : Maximum total number of basis orbitals'
      write(1,'(a)')
     .'c  maxpul : Dimension of Pulay auxiliary vectors'
      write(1,'(a)') 
     .'c  maxspn : Maximum number of spin polarizations (1 or 2)'
      write(1,'(a)') 
     .'c  maxua  : Maximum number of atoms in unit cell'
      write(1,'(a)') 
     .'c  maxuo  : Maximum number of basis orbitals in unit cell'
      write(1,'(a)')
     .'c  maxxij : Dimension of array xijo'

      write(1,'(6x,a)') 'integer maxa, maxbk, maxk, maxkb'
      write(1,'(6x,a)') 'integer maxna, maxno, maxo, maxpul'
      write(1,'(6x,a)') 'integer maxspn, maxua, maxuo'
      write(1,'(6x,a)') 'integer maxxij'
      write(1,'(6x,a,i12,a)') 'parameter ( maxa   =', na,        ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxbk  =', max(nbk,1),' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxk   =', nkpnt,     ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxkb  =', nokb,      ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxna  =', nnamax,    ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxno  =', nnomax,    ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxo   =', no,        ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxpul =', nauxpul,   ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxspn =', nspin,     ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxua  =', nua,       ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxuo  =', nuo,       ' )'
      write(1,'(6x,a,i12,a)') 'parameter ( maxxij =', nxij,      ' )'
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


