      subroutine redata(maxa, maxl, maxs, maxspn, maxzet,
     .                  amax, lmax, smax, spnmax, zetmax, overflow,
     .                  slabel, sname,
     .                  na, ns, nspin, isa, izs, lmaxs, lmxkbs, 
     .                  nzls, rcls, contrf, atm_label,
     .                  cell, ucell, indxua, xa, smass, g2max, negl, 
     .                  nscf, dDtol, wmix, isolve,
     .                  temp, ncgmax, ftol, eta, etol, rcoor, ioptlwf,
     .                  idyn, istart, ifinal, nmove, ianneal, iquench,
     .                  dt,  dxmax, tt, tp, mn, mpr, bulkm, taurelax, 
     .                  writedim, usesavelwf, usesavedm, mullipop,
     .                  inspn, maxsav, pulfile, tempinit)
C *********************************************************************
C Subroutine to read the data for the SIESTA program
C
C     It uses the FDF (Flexible Data Fromat) package 
C     of J.M.Soler and A.Garcia
C
C     This routine is compatible with fdf Version 0.6.1.5
C 
C Writen by P.Ordejon, December'96
C ***************************** INPUT *********************************
C integer maxa             : Maximum number of atoms
C integer maxl             : Maximum angular momentum of basis orbitals
C integer maxs             : Maximum number of species
C integer maxspn           : Maximum number of spins (1 or 2)
C integer maxzet           : Maximum number of Zetas of any atom
C **************************** OUTPUT *********************************
C integer amax             : Actual maxa for this system
C integer lmax             : Actual maxl for this system
C integer smax             : Actual maxs for this system
C integer spnmax           : Actual maxspn for this system
C integer zetmax           : Actual maxzet for this system
C logical overflow         : True = Some of the dimensions is too small
C character*20 slabel      : System Label (to name output files)
C character*150 sname      : System Name
C integer na               : Number of atoms
C integer ns               : Number of species
C integer nspin            : Spin polarization
C integer isa(maxa)        : Species index of each atom
C integer izs(maxs)        : Atomic Number of each species
C integer lmaxs(maxs)      : Angular momentum cutoff for the orbitals 
C                            of each species
C integer lmxkbs(maxs)     : Angular momentum cutoff for the KB proj.
C                            of each species
C integer nzls(0:maxl,maxs): Number of zetas for each l of each species
C real*8 rcls(maxzet,0:maxl,maxs) : Cutoff R of atomic orbitals (Bohr)
C real*8 contrf(maxzet,0:maxl,maxs) : Scaling factor for orbitals
C real*8 cell(3,3)         : (Super) lattice vectors CELL(ixyz,ivector)
C                            (in Bohr)
C real*8 ucell(3,3)        : Unit-cell vectors (Bohr)
C integer indxua(maxa)     : Index of equivalent atom in unit cell
C real*8 xa(3,maxa)        : Atomic coordinates (Bohr)
C real*8 smass(maxs)       : Atmic masses of each species (amu)
C real*8 g2max             : PW cutoff energy (Ry)
C logical negl             : True = Neglect interactions between
C                            non-overlaping orbitals (coming from
C                            KB projectors)
C integer nscf             : Maximum number of SCF cycles per time step
C real*8 dDtol             : Maximum Density Matrix tolerance
C real*8 wmix              : Amount of output DM for new DM
C integer isolve           : Method of solution.  0 = Diagonalization
C                                                 1 = Order-N
C real*8 temp              : Temperature for Fermi smearing (Ry)
C integer ncgmax           : Maximum number of CG steps for 
C                            band structure energy minimization
C real*8 etol              : Relative tolerance in CG minimization
C                            of band structure energy
C real*8 eta               : Fermi level parameter of Kim functional
C real*8 rcoor             : Cutoff radius of LWF's (Bohr)
C integer ioptlwf          : Option to build LWF's according to:
C                             0 = Read blindly from disk
C                             1 = Functional of Kim et al.
C                             2 = Functional of Ordejon-Mauri
C integer idyn             : Atomic dynamics option:
C                             0 = CG geometry optimization
C                             1 = Standard MD run (Verlet)
C                             2 = Nose thermostat MD
C                             3 = Parrinello-Rahman MD
C                             4 = Nose thermostat + Parrinello-Rahman MD
C                             5 = Annealing MD
C integer istart           : Initial time step for MD
C integer ifinal           : Final time step for MD
C integer nmove            : Number of CG steps in CG optimization
C real*8 ftol              : Maximum force for CG structure optimization
C integer ianneal          : Annealing option for idyn = 5
C                             1 = Temperature 
C                             2 = Pressure
C                             3 = Temperature and Pressure
C integer iquench          : Quench option: 0 = No;  1 = Yes
C real*8 dt                : Length of time step (fs)
C real*8 dxmax             : Maximum atomic displacement in one CG move
C real*8 tt                : Target temperature (Kelvin)
C real*8 tp                : Target Pressure (Ry/Bohr**3)
C real*8 mn                : Mass of Nose variable (Ry/fs**2)
C real*8 mpr               : Mass of Parrinello-R. variable (Ry/fs**2)
C real*8 bulkm             : Estimate of bulk modulus (Ry/Bohr**3)
C real*8 taurelax          : Annealing time to reach targer T and P (fs)
C logical writedim         : True = dump minimum dimensions to siesta.h
C                              and stop.
C logical usesavelwf       : True = try to use continuation LWF files 
C                              from disk
C logical usesavedm        : True = try to use continuation DM files 
C                              from disk
C integer mullipop         : Option for Mulliken Pop. analisys
C logical inspn            : Spin initialization for spin-polarized
C                              .true.  -> Antiferro
C                              .false. -> Ferro
C integer maxsav           : Number of density-matrices stored for Pulay
C                            mixing. Pulay mixing is done every maxsav 
C                            iterations, the rest is linear mixing.
C                              .lt.2 => linear mixing only
C                              .ge.2 => pulay mixing
C logical pulfile          : Use file (.true.) or memory (.false.)
C                            to store Pulay miximg intermediate vectors
C                            Default: .false.
C real*8 tempinit          : Initial temperature (Kelvin) of the MD simulation
C **********************************************************************

      implicit none

      integer
     .  maxa, maxl, maxs, maxspn, maxzet

      character
     .  slabel*20, sname*150, atm_label(*)*20

      integer
     .  amax, lmax, smax, spnmax, zetmax

      integer
     .  ianneal, idyn, ifinal, ioptlwf,
     .  iquench, isolve, istart, maxsav,
     .  mullipop, na, ncgmax, nmove, ns, nscf, nspin

      integer 
     .  indxua(maxa), isa(maxa), izs(maxs), lmaxs(maxs), 
     .  lmxkbs(maxs), nzls(0:maxl,maxs)

      double precision
     .  bulkm, cell(3,3), contrf(maxzet,0:maxl,maxs), 
     .  dcell(3,3), dDtol, dt, dxmax, eta, etol, ftol, g2max,
     .  mn, mpr, rcoor, rcls(maxzet,0:maxl,maxs), smass(maxs),
     .  taurelax, temp, tempinit, tp, tt, ucell(3,3), wmix, 
     .  xa(3,maxa)

      logical
     .  inspn, negl, overflow, pulfile,
     .  usesavelwf, usesavedm, writedim

C Internal parameters ................................................
C na_diag      : maximun number of atoms with diagon as default method
C g2max_defect : Mesh cutoff default, in Ry
C temp_defect  : Electronic temperature default, in Ry
      integer na_diag
      double precision g2max_defect, temp_defect
      parameter ( na_diag      = 100      )
      parameter ( g2max_defect = 50.d0    )
      parameter ( temp_defect  = 1.900d-3 )
C ................

C  Internal variables .................................................
      character
     .  annop*22, dyntyp*22, 
     .  filein*20, fileout*20, 
     .  method*6, line*150, names_parse*150,
     .  lwfopt*13

      character
     .  annop_defect*22, dyntyp_defect*22, 
     .  lwfopt_defect*13,  slabel_defect*150, 
     .  sname_defect*20, atomic_label*20,symbol*2

      integer 
     .  ifinal_defect, istart_defect, maxsv_defect, mpop_defect,
     .  na_defect, ncgmax_defect, nmove_defect,
     .  ns_defect, nscf_defect

      integer 
     .  i, i1, i2, i3, ia, ic, il, is, iunit, ix, izeta, izsr, 
     .  j, js, l, lmaxsr, lmxkbsr, maux(3,3,2), ml(3,3), mr(3,3),
     .  ncells, ni, nn, nr, nscell(3,3), nsd(3,3), nua, nv, nzlsr,
     .  integers_parse(20),lc(0:20),lastc

      double precision
     .  alat, atmass

      double precision
     .  bulkm_defect, 
     .  dDtol_defect, dt_defect, dxmax_defect,
     .  eta_defect, etol_defect, ftol_defect,
     .  mn_defect, mpr_defect, rcoor_defect, 
     .  taurelax_defect, ti_defect, tp_defect, 
     .  tt_defect, wmix_defect, values_parse(20),
     .  reals_parse(20), volcel

      double precision
     .  alp, blp, clp, alplp, betlp, gamlp, pi, xxx

      logical
     .  leqi, qnch, sppol,chemblock 

      logical
     .  inspn_defect, negl_defect, pul_defect,
     .  qnch_defect, sppol_defect, 
     .  usdm_defect, uslwf_defect, wd_defect

C ................

C Define FDF calls ....................................................
      include 'fdf/fdfdefs.h'
C ................

      data pi / 3.1415926d0 /

      overflow = .false.
      amax = 1
      lmax = 0
      smax = 1
      spnmax = 1
      zetmax = 1

C Print Welcome and Presentation .......................................
      write(6,'(/a)') 
     . '                           ***********************       '
      write(6,'(a)') 
     . '                           *  WELCOME TO SIESTA  *       '
      write(6,'(a)')
     . '                           ***********************       '
C ..................

C Dump data file to output file .......................................
      write(6,'(/,a,18(1h*),a,28(1h*))')
     .  'redata: ', ' Dump of input data file '
   10 continue
        read(5,'(a80)',end=20) line
        do ic = len(line),1,-1
          if (line(ic:ic) .ne. ' ') goto 15
        enddo
        ic = 0
   15   if (ic .gt. 0) write(6,'(a)') line(1:ic)
      goto 10
   20 continue
      write(6,'(a,18(1h*),a,29(1h*),/)')
     .  'redata: ', ' End of input data file '
C ..................

C Read data from FDF file..............................................

C Set up fdf ...
      filein = 'stdin'
      fileout = 'out.fdf'
      call fdf_init(filein,fileout)
C ...

      write(6,'(/,a,18(1h*),a,30(1h*))')
     .  'redata: ', ' Simulation parameters '
      write(6,'(a)')  'redata:'
      write(6,'(a)')  'redata:  The following are some of the '//
     .                         'parameters of the simulation.'
      write(6,'(a)')  'redata:  A complete list of the parameters '//
     .                         'used, including defect values,'
      write(6,'(a,a)')'redata:  can be found in file ',fileout
      write(6,'(a)')  'redata:'

C Defile Name of the system ...
      sname_defect = ' '
      sname = fdf_string('SystemName',sname_defect)
      write(6,'(a,a)') 
     . 'redata: System Name                      = ',sname
C ...

C Defile System Label (short name to label files) ...
      slabel_defect  = 'siesta'
      slabel = fdf_string('SystemLabel',slabel_defect)
      write(6,'(a,a)') 
     . 'redata: System Label                     = ',slabel
C ...

C Read Number of Atoms ...
      na_defect = 0
      na = fdf_integer('NumberOfAtoms',na_defect)
      write(6,'(a,i5)') 
     . 'redata: Number of Atoms                  = ',na
      if (na .le. 0) then
        write(6,100) 
        write(6,101) 
        write(6,'(a)') 
     . 'redata: ERROR: Number of atoms must be larger than zero.'
        write(6,102)
        stop
      endif
C     check if dimension of number of atomos is large enough
      call chkdime(maxa, na, overflow, amax)
C ...

C Defile Number of species ...
      ns_defect = 0
      ns = fdf_integer('NumberOfSpecies',ns_defect)
      write(6,'(a,3x,i2)') 
     . 'redata: Number of Atomic Species         = ',ns
      if (ns .le. 0) then
        write(6,100) 
        write(6,101) 
        write(6,'(a)') 
     . 'redata: ERROR: Number of species must be larger than zero.'
        write(6,102)
        stop
      endif
C     check if dimension of number of species is large enough
      call chkdime(maxs, ns, overflow, smax)
C ...

C Dump minimum sizes to siesta.h and stop
      wd_defect = .false.
      writedim  = fdf_boolean('WriteSiestaDim',wd_defect)
      write(6,'(a,4x,l1)') 
     . 'redata: Write minimum siesta.h dimensions  = ',writedim
C ...

C Perform Mulliken Population Analisys
      mpop_defect = 0
      mullipop = fdf_integer('WriteMullikenPop',mpop_defect)
      if (mullipop .eq. 0) then
      write(6,'(a,4x,l1)') 
     . 'redata: Write Mulliken Pop. = NO'
      elseif (mullipop .eq. 1) then
      write(6,'(a,4x,l1)') 
     . 'redata: Write Mulliken Pop. = Atomic and Orbital charges'
      elseif (mullipop .eq. 2) then
      write(6,'(a,4x,l1)') 
     . 'redata: Write Mulliken Pop. = Atomic and Orbital charges'
      write(6,'(a,4x,l1)') 
     . '                              plus Atomic Overlap Populations'
      elseif (mullipop .eq. 3) then
      write(6,'(a,4x,l1)') 
     . 'redata: Write Mulliken Pop. = Atomic and Orbital charges'
      write(6,'(a,4x,l1)') 
     . '                              plus Atomic Overlap Populations'
      write(6,'(a,4x,l1)') 
     . '                              plus Overlap Overlap Populations'
      else
        stop 'redata: Wrong value for WriteMullikenPop'
      endif

C Next lines added by D. Sanchez-Portal. 8-4-1997.
C Define the chemical especies of the atoms.

      chemblock=fdf_block('PAO_basis_and_PS_lmax',i)
      if (fdf_block('Chemical_species_label',
     .         iunit)) then 
        do js=1,ns
          read(iunit,'(a)') line
          lastc = index(line,'#') - 1
          if(lastc.eq.-1) lastc=150
          call parse( line(1:lastc), 
     .         nn, lc, names_parse, 
     .            nv, values_parse,
     .            ni, integers_parse,
     .            nr, reals_parse )
          if (nn.gt.1) then 
             write(6,'(a,i3)') 
     .       'redata: More than one label for species',js
              write(6,'(a,/,(2x,a))')
     .       'label =', (names_parse(lc(i-1)+1:lc(i)),i=1,nn)
              stop 
          endif 
          if (ni.ne.2) then 
             write(6,'(a)') 
     .         'redata: ERROR: number of integer read in '
             write(6,'(a)')
     .        '%block  Chemical_species_label, ni=',ni
              stop 
         endif

          is=integers_parse(1) 
          izsr=integers_parse(2)
        
          if (is .ne. js)
     .    stop 'redata: Unexp. index is. Order species consecutively.'
          izs(is) = izsr

C  Floating orbitals have negative atomic numbers.
          if(izsr.gt.0) then 
            smass(is) = atmass(izsr)
          else
            smass(is) = 0.0d0
          endif 
C  At this  moment the atomic_label is not read......    
           atomic_label=names_parse(1:lc(1))
           atm_label(is)=atomic_label
           write(6,'(a,i3,3x,a)')
     .      'redata: label for specie ',is,atomic_label
        enddo  
      else
       if(.not.chemblock)then 
         write(6,100)
          write(6,101)
          write(6,'(a)')
     .   'redata: ERROR: The chemical species must be specified.'
          write(6,102)
          stop
        endif
      endif 
       

C Define Basis Set and Pseudopotentials...
  
      chemblock=fdf_block('Chemical_species_label',i)
      if ( fdf_block('PAO_basis_and_PS_lmax',iunit) ) then
        do js = 1,ns
          read(iunit,*) is, izsr, lmaxsr, lmxkbsr
          if(izsr.eq.-100) lmxkbsr=0
C  Check if there are previously read atomic numbers.
          if(chemblock) then
             if(izs(js).ne.izsr) then 
               write(6,'(2a)') 'redata: ERROR:',
     .  ' Atomic numbers specified in Chemical_species_label block'
               write(6,'(3a)') 'redata: ERROR:',
     .  ' are NOT consistent with those specified in',
     .  ' PAO_basis_and_PS_lmax block'
             stop 
             endif 
          endif 

          if (is .ne. js)
     .    stop 'redata: Unexp. index is. Order species consecutively.'
          call chkdime(maxl, lmaxsr, overflow, lmax)
          if (.not.overflow) then
            izs(is)=izsr
            lmaxs(is) = lmaxsr
            lmxkbs(is) = lmxkbsr
          endif
          do il = 0,lmaxsr
            read(iunit,*) l, nzlsr
            if (l .ne. il) stop 'redata: Unexpected index l'
            call chkdime(maxzet, nzlsr, overflow, zetmax)
            if (.not.overflow) then
              nzls(l,is) = nzlsr
              read(iunit,*) (rcls(izeta,l,is),izeta=1,nzls(l,is))
              read(iunit,*) (contrf(izeta,l,is),izeta=1,nzls(l,is))
            else
              read(iunit,*)
              read(iunit,*)
            endif
          enddo
        enddo
        if(.not.chemblock) then 
          do js=1,ns
            if(izs(js).ne.-100) then 
                atm_label(js)=symbol(abs(izs(js))) 
            else
                atm_label(js)='Bessel'
            endif 

C  Floating orbitals have negative atomic numbers.
            if(izs(js).gt.0) then 
              smass(js) = atmass(izs(js))
            else
              smass(js) = 0.0d0
            endif 
          enddo 
        endif 

      else
C    Next lines modified by D. Sanchez-Portal. 8-4-1997.
C    Now  it is not necessary to specify the basis set. 
C    The program is able to provided you a 'default' basis-set.
C    Anyway, the performance of this basis should be tested.
C    There is no guarantee that this is the 'best' basis!!!!

      do js=1,ns 
        lmaxs(js)=-1
        lmxkbs(js)=-1
          
        do il=0,maxl
          nzls(il,js)=0     
          do izeta=1,maxzet
             rcls(izeta,il,js)=0.0d0 
             contrf(izeta,il,js)=1.0d0
          enddo 
        enddo 
      enddo 
      endif
C ...

C Read atomic masses from fdf if a change is wanted 

      call remass(smass,maxs,ns)

C Lattice constant
      alat = fdf_physical('LatticeConstant',0.d0,'Bohr')
      if (alat .eq. 0.d0) then
        write(6,'(a,3(/,a))') 
     .   'redata:',
     .   'redata: WARNING: No valid lattice constant specified',
     .   'redata: WARNING: Cell will be generated automatically',
     .   'redata:'
      endif

C Lattice vectors

      if ( fdf_block('LatticeParameters',iunit) .and.
     .     fdf_block('LatticeVectors',iunit) ) then
         write(6,'(2a)')'redata: ERROR: Lattice vectors doubly ',
     .     'specified: by LatticeVectors and by LatticeParameters.' 
         stop 'redata: ERROR: Double input for lattice vectors'
      endif

      if ( fdf_block('LatticeParameters',iunit) ) then
         read(iunit,*) alp, blp, clp, alplp, betlp, gamlp
         write(6,'(a)')
     .    'redata: Lattice Parameters (units of Lattice Constant) ='
         write(6,'(a,3f10.5,3f9.3)')
     .    'redata: ',alp,blp,clp,alplp,betlp,gamlp
         alplp = alplp * pi/180.d0
         betlp = betlp * pi/180.d0
         gamlp = gamlp * pi/180.d0
         cell(1,1) = alp
         cell(2,1) = 0.d0
         cell(3,1) = 0.d0
         cell(1,2) = blp * cos(gamlp)
         cell(2,2) = blp * sin(gamlp)
         cell(3,2) = 0.d0
         cell(1,3) = clp * cos(betlp)
         xxx = (cos(alplp) - cos(betlp)*cos(gamlp))/sin(gamlp)
         cell(2,3) = clp * xxx
         cell(3,3) = clp * sqrt(sin(betlp)*sin(betlp) - xxx*xxx)
      elseif ( fdf_block('LatticeVectors',iunit) ) then
        do i = 1,3
          read(iunit,*) (cell(j,i), j=1,3)
        enddo
      else
        do i = 1,3
          do j  = 1,3
            cell(i,j) = 0.d0
          enddo
          cell(i,i) = 1.d0
        enddo
      endif

      if (alat .ne. 0.d0) then
        write(6,'(a,f10.4,a)') 
     .   'redata: Lattice Constant                 = ',alat,'  Bohr'
        write(6,'(a)') 
     .   'redata: Lattice vectors (in units of Lattice Constant) ='
        do i = 1,3
          write(6,'(a,3f10.5)')
     .   '        ',(cell(j,i), j=1,3)
        enddo
      endif
C ...

C Multiply cell vectors by by lattice constant ........................
      do i = 1,3
        do ix = 1,3
          cell(ix,i) = alat * cell(ix,i)
        enddo
      enddo
      if (alat .ne. 0.d0) then
        write(6,'(a)') 
     .   'redata: Lattice vectors (in Bohr) ='
        do i = 1,3
          write(6,'(a,3f10.5)') '        ',(cell(j,i), j=1,3)
        enddo
      endif
C ..................

C Spin Polarization ...
      sppol_defect = .false.
      sppol = fdf_boolean('SpinPolarized',sppol_defect)
      if (sppol) then
        nspin = 2
      else 
        nspin = 1
      endif
      write(6,'(a,4x,i1)') 
     . 'redata: Spin Polarization                = ',nspin
      call chkdime(maxspn, nspin, overflow, spnmax)
C ...
 
C Planewave cutoff of the real space mesh ...
      g2max = fdf_physical('MeshCutoff',g2max_defect,'Ry')
      write(6,'(a,f10.4,a)') 
     . 'redata: Mesh Cutoff                      = ',g2max,'  Ry'
C ...
       
C SCF Loop parameters ...
C     Maximum number of SCF iterations
      nscf_defect = 50
      nscf = fdf_integer('MaxSCFIterations',nscf_defect)
      write(6,'(a,i5)') 
     . 'redata: Max. number of SCF Iter          = ',nscf

C     Pulay mixing, numer of iterations for one Pulay mixing (maxsav)
      maxsv_defect = 0
      maxsav = fdf_integer('DM.NumberPulay',maxsv_defect)
      if(maxsav .gt. 1) then
        write(6,'(a,i5,a)') 
     .   'redata: One Pulay mixing every           = ',maxsav,
     .   ' iterations'
      else
         write(6,'(a)')'redata: Mixing is linear'
      endif

C     Use disk or memory to store intermediate Pulay miximg vectors
C     (pulfile)
      pul_defect = .false.
      pulfile  = fdf_boolean('DM.PulayOnFile',pul_defect)
      write(6,'(a,4x,l1)')
     .  'redata: Write Pulay info on disk?       = ',pulfile
C ...


C     Density Matrix Mixing  (proportion of output DM in new input DM)
      wmix_defect = 0.25d0
      wmix = fdf_double('DM.MixingWeight',wmix_defect)
      write(6,'(a,f10.4,a)') 
     . 'redata: New DM Mixing Weight             = ',wmix

C     Density Matrix Tolerance for achieving Self-Consistency
      dDtol_defect = 1.d-4
      dDtol = fdf_double('DM.Tolerance',dDtol_defect)
      write(6,'(a,f12.6,a)') 
     . 'redata: DM Tolerance for SCF             = ',dDtol

C     Initial spin density: Maximum polarization, Ferro (false), AF (true)
      if(sppol) then
        inspn_defect = .false.
        inspn = fdf_boolean('DM.InitSpinAF',inspn_defect)
        write(6,'(a,4x,l1)')
     .   'redata: Antiferro initial spin density   = ',inspn
      endif
C ...

C Use continumation DM files
*     usdm_defect = .true.
      usdm_defect = fdf_boolean('UseSaveData',.false.)
      usesavedm  = fdf_boolean('DM.UseSaveDM',usdm_defect)
      write(6,'(a,4x,l1)') 
     . 'redata: Use continuation files for DM    = ',
     . usesavedm
C ...

C Neglect Interactions between non-overlapping orbitals ...
      negl_defect = .false.
      negl  = fdf_boolean('NeglNonOverlapInt',negl_defect)
      write(6,'(a,4x,l1)') 
     . 'redata: Neglect nonoverlap interactions  = ',negl
C ...

C Method to Solve LDA Hamitonian ...
      if (na .le. na_diag) then
        method = fdf_string('SolutionMethod','diagon')
      else
        method = fdf_string('SolutionMethod','ordern')
      endif
      if (leqi(method,'diagon')) then
        isolve = 0
        write(6,'(a,a)') 
     .   'redata: Method of Calculation            = ',
     .   'Diagonalization'
      else if (leqi(method,'ordern')) then
        isolve = 1
        write(6,'(a,a)') 
     .   'redata: Method of Calculation            = ',
     .   'Order-N'
        if (nspin .eq. 2) then
         write(6,100) 
         write(6,101) 
         write(6,'(a)') 
     .   'redata:    You chose the Order-N solution option'
         write(6,'(a)') 
     .   'redata:    together with nspin=2.  This is not  '
         write(6,'(a)') 
     .   'redata:    allowed in this version of siesta    '
         write(6,102)
         stop
        endif
      else
        write(6,100) 
        write(6,101) 
        write(6,'(a)') 
     .  'redata:    The method of solution must be either'
        write(6,'(a)') 
     .  'redata:           OrderN or Diagon'
        write(6,102)
        stop
      endif
C ...

C Electronic temperature for Fermi Smearing ...
      temp = fdf_physical('ElectronicTemperature',temp_defect,'Ry')
      if (isolve .eq. 0) then
        write(6,'(a,f10.4,a)') 
     .  'redata: Electronic Temperature           = ',temp,'  Ry'
      endif
C ...

C Order-N solution parameters ...
C     Maximum number of CG minimization iterations
      ncgmax_defect = 1000
      ncgmax = fdf_integer('ON.MaxNumIter',ncgmax_defect)

C     Relative tolerance in total band structure energy
      etol_defect = 1.d-8
      etol = fdf_double('ON.etol',etol_defect)

C     Fermi level parameter
      eta_defect = 0.d0
      eta = fdf_physical('ON.eta',eta_defect,'Ry')

C     Cutoff radius for Localized Wave Functions
      rcoor_defect = 9.5d0
      rcoor = fdf_physical('On.RcLWF',rcoor_defect,'Bohr')

C     Use continumation LWF files
*     uslwf_defect = .true.
      uslwf_defect = fdf_boolean('UseSaveData',.false.)
      usesavelwf  = fdf_boolean('ON.UseSaveLWF',uslwf_defect)

C     Option on how to build LWF's (disk or functionals)
      lwfopt_defect = 'kim'
      lwfopt  = fdf_string('ON.functional',lwfopt_defect)
      if (leqi(lwfopt,'files')) then
        ioptlwf = 0
      else if (leqi(lwfopt,'kim')) then
        ioptlwf = 1
      else if (leqi(lwfopt,'ordejon-mauri')) then
        ioptlwf = 2
      else
        write(6,'(a)') 'redata: wrong ON.funcional option'
        stop
      endif
C ...


      if (isolve .eq. 1) then
        write(6,'(a,i5)') 
     .  'redata: Maximum number of iterations     = ',ncgmax
        write(6,'(a,d12.2)') 
     .  'redata: Relative tolerance               = ',etol
        write(6,'(a,f10.4,a)') 
     .  'redata: Eta (Fermi level parameter)      = ',eta,'  Ry'
        write(6,'(a,f8.2,a)') 
     .  'redata: Radius of LWFs                   = ',rcoor,
     .  '    Bohr'
        write(6,'(a,4x,l1)') 
     .  'redata: Use continuation files for LWF = ',
     .  usesavelwf
        write(6,'(a,a)') 
     .  'redata: Method to build LWFs = ',lwfopt
      endif
C ...

C Dynamics parameters ...
C     Type of dynamics
      dyntyp_defect = 'verlet'
      dyntyp = fdf_string('MD.TypeOfRun',dyntyp_defect)
      if (leqi(dyntyp,'cg')) then
        idyn = 0
        write(6,'(a,a)') 
     .   'redata: Dynamics option                  = ',
     .   'CG coord. optimization'
      else if (leqi(dyntyp,'verlet')) then
        idyn = 1
        write(6,'(a,a)') 
     .   'redata: Dynamics option                  = ',
     .   'Verlet MD run'
      else if (leqi(dyntyp,'nose')) then
        idyn = 2
        write(6,'(a,a)') 
     .   'redata: Dynamics option                  = ',
     .   'Nose termostat MD run'
      else if (leqi(dyntyp,'parrinellorahman')) then
        idyn = 3
        write(6,'(a,a)') 
     .   'redata: Dynamics option                  = ',
     .   'Parrinello-Rahman MD run'
      else if (leqi(dyntyp,'noseparrinellorahman')) then
        idyn = 4
        write(6,'(a,a)') 
     .   'redata: Dynamics option                  = ',
     .   'Nose-Parrinello-Rahman MD run'
      else if (leqi(dyntyp,'anneal')) then
        idyn = 5
        write(6,'(a,a)') 
     .   'redata: Dynamics option                  = ',
     .   'Annealing MD run'
      else
        write(6,100) 
        write(6,101) 
        write(6,'(a)') 'redata:  Wrong Dynamics Option Selected       '
        write(6,'(a)') 'redata:  You must choose one of the following:'
        write(6,'(a)') 'redata:                                       '
        write(6,'(a)') 'redata:      - CG                             '
        write(6,'(a)') 'redata:      - Verlet                         '
        write(6,'(a)') 'redata:      - Nose                           '
        write(6,'(a)') 'redata:      - Parrinello-Rahman              '
        write(6,'(a)') 'redata:      - Nose-Parrinello-Rahman         '
        write(6,'(a)') 'redata:      - Anneal                         '
        write(6,102)
        stop
      endif
    
C     Maximum number of steps in CG coordinate optimization
      nmove_defect = 0
      nmove = fdf_integer('MD.NumCGsteps',nmove_defect)

C     Maximum atomic displacement in one CG step
      dxmax_defect = 0.2d0
      dxmax = fdf_physical('MD.MaxCGDispl',dxmax_defect,'Bohr')

C     Tolerance in the maximum atomic force
      ftol_defect = 0.01d0 
      ftol = fdf_physical('MD.MaxForceTol',ftol_defect,'Ry/Bohr')

      if (idyn .eq. 0) then
        write(6,'(a,i5)') 
     .  'redata: Maximum number of CG moves       = ',nmove
        write(6,'(a,f10.4,a)') 
     .  'redata: Max atomic displ per move        = ',dxmax,'  Bohr'
        write(6,'(a,f10.4,a)') 
     .  'redata: Force tolerance                  = ',ftol,'  Ry/Bohr'
      endif
  
C     Initial time step for MD
      istart_defect = 1
      istart = fdf_integer('MD.InitialTimeStep',istart_defect)

C     Final time step for MD
      ifinal_defect = 1
      ifinal = fdf_integer('MD.FinalTimeStep',ifinal_defect)

C     Length of time step for MD
      dt_defect = 1.d0
      dt = fdf_physical('MD.LengthTimeStep',dt_defect,'fs')

C     Quench Option
      qnch_defect = .false.
      qnch = fdf_boolean('MD.Quench',qnch_defect)
      if (qnch .and. (idyn .eq. 2 .or. idyn .eq. 4)) then
        write(6,100) 
        write(6,101) 
        write(6,'(a)') 
     .  'redata: ERROR: You cannot quench and use a Nose'
        write(6,'(a)') 
     .  'redata: ERROR: thermostat simultaneously.'
        write(6,102)
        stop
      endif
      iquench = 0
      if (qnch) iquench = 1

C     Initial Temperature of MD simulation
C     (draws random velocities from the Maxwell-Boltzmann distribition
C      at the given temperature)
      ti_defect = 0.d0
      tempinit = fdf_physical('MD.InitialTemperature',ti_defect,'K')

      if (idyn .ge. 1 .and. idyn .le. 5) then
        write(6,'(a,i5)') 
     .  'redata: Initial MD time step             = ',istart
        write(6,'(a,i5)') 
     .  'redata:   Final MD time step             = ',ifinal
        write(6,'(a,f10.4,a)') 
     .  'redata: Length of MD time step           = ',dt,'  fs'
        write(6,'(a,f10.4,a)') 
     .  'redata: Length of MD time step           = ',dt,'  fs'
        write(6,'(a,f10.4,a)') 
     .  'redata: Initial Temperature of MD run    = ',tempinit,'  K'
        if (idyn .ne. 5) then
          write(6,'(a,4x,l1)') 
     .    'redata: Perform a MD quench              = ',qnch
        endif
      endif

C     Target Temperature
      tt_defect = 0.d0
      tt = fdf_physical('MD.TargetTemperature',tt_defect,'K')
      
C     Target Pressure
      tp_defect = 0.d0
      tp = fdf_physical('MD.TargetPressure',tp_defect,'Ry/Bohr**3')

C     Mass of Nose variable
      mn_defect = 1.d2
      mn = fdf_physical('MD.NoseMass',mn_defect,'Ry*fs**2')

C     Mass of Parrinello-Rahman variables
      mpr_defect = 1.d2
      mpr = fdf_physical('MD.ParrinelloRahmanMass',
     .                    mpr_defect,'Ry*fs**2')

      if (idyn .eq. 2 .or. idyn .eq. 4) then
        write(6,'(a,f10.4,a)') 
     .  'redata: Nose mass                        = ',mn,'  Ry/fs**2'
      endif

      if (idyn .eq. 3 .or. idyn .eq. 4) then
        write(6,'(a,f10.4,a)') 
     .  'redata: Parrinello-Rahman mass           = ',mpr,'  Ry/fs**2'
      endif

C     Annealing option
      annop_defect = 'TemperatureAndPressure'
      annop = fdf_string('MD.AnnealOption',annop_defect)
      if (idyn .eq. 5) then
        if (leqi(annop,'Temperature')) then
          ianneal = 1
          write(6,'(a,a)') 
     .     'redata: Annealing Option                 = ',
     .     'Temperature'
        else if (leqi(annop,'Pressure')) then
          ianneal = 2
          write(6,'(a,a)') 
     .     'redata: Annealing Option                 = ',
     .     'Pressure'
        else if (leqi(annop,'TemperatureAndPressure')) then
          ianneal = 3
          write(6,'(a,a)') 
     .     'redata: Annealing Option                 = ',
     .     'Temperature and Pressure'
        else
          write(6,100) 
          write(6,101) 
          write(6,'(a)') 
     .    'redata:           You have chosen annealing MD, and you '
          write(6,'(a)') 
     .    'redata:           must use one of the following options:'
          write(6,'(a)') 
     .    'redata:           - Temperature                         '
          write(6,'(a)') 
     .    'redata:           - Pressure                            '
          write(6,'(a)') 
     .    'redata:           - TemperatureAndPressure              '
          write(6,102)
          stop
        endif
      endif

      if (idyn .eq. 2 .or. idyn .eq. 4 .or. 
     .   (idyn .eq. 5 .and. (ianneal .eq. 1 .or. ianneal .eq. 3))) then
        write(6,'(a,f10.4,a)') 
     .  'redata: Target Temperature               = ',tt,'  Kelvin'
      endif

      if (idyn .eq. 3 .or. idyn .eq. 4 .or. 
     .   (idyn .eq. 5 .and. (ianneal .eq. 2 .or. ianneal .eq. 3))) then
        write(6,'(a,f10.4,a)') 
     .  'redata: Target Pressure                  = ',tt,'  Ry/Bohr**3'
      endif

C     Relaxation Time for Annealing
      taurelax_defect = 1.d2
      taurelax = fdf_physical('MD.TauRelax',taurelax_defect,'fs')
      if (idyn .eq. 5) then
        write(6,'(a,f10.4,a)') 
     .  'redata: Annealing Relaxation Time        = ',
     .   taurelax,'  fs'
      endif
        
C     Estimated Bulk modulus (for Pressure annealing)
      bulkm_defect = 1.d2
      bulkm = fdf_double('MD.BulkModulus',bulkm_defect)
      if (idyn .eq. 5 .and. (ianneal .eq. 2 .or. ianneal .eq. 3)) then
        write(6,'(a,f10.4,a)') 
     .  'redata: Approx. Bulk Modulus             = ',
     .   bulkm,'  Ry/Bohr**3'
      endif
C ...

C Atomic Coordinates ..................................
      call recoor(overflow, cell, alat, xa, isa, na)
      
C Super cell ..........................................
      if ( fdf_block('SuperCell',iunit) ) then
        if ( alat .eq. 0.d0 ) then
          write(6,*) 'redata: ERROR: LatticeConstant required',
     .                      ' to define SuperCell'
          stop 'redata: ERROR: LatticeConstant required'
        endif
        do i = 1,3
          read(iunit,*) (nscell(j,i),j=1,3)
        enddo
      else
        do i = 1,3
          do j = 1,3
            nscell(j,i) = 0
          enddo
          nscell(i,i) = 1
        enddo
      endif
      do i = 1,3
        do ix = 1,3
          ucell(ix,i) = cell(ix,i)
        enddo
      enddo
      do i = 1,3
        do ix = 1,3
          cell(ix,i) = ucell(ix,1) * nscell(1,i) +
     .                 ucell(ix,2) * nscell(2,i) +
     .                 ucell(ix,3) * nscell(3,i)
        enddo
      enddo
      if ( volcel(cell) .lt. 1.d-12 ) then
        ncells = 1
      else
        ncells = nint( volcel(cell) / volcel(ucell) )
      endif
      nua = na
      na = nua * ncells
      call chkdime(maxa, na, overflow, amax)
      if (ncells.gt.1) then
        write(6,'(a,/,(a,3f12.6))')
     .    'redata: Total-cell (supercell) vectors (Bohr)',
     .    ('redata:',(cell(ix,i),ix=1,3),i=1,3)
        write(6,'(a,i6)') 'redata: Number of unit cells  =', ncells
        write(6,'(a,i6)') 'redata: Total number of atoms =', na
      endif
      if (.not.overflow) then
*       do i = 1,3
*         do j = 1,3
*           if (j.ne.i .and. nscell(j,i).ne.0) then
*             write(6,*) 'redata: ERROR: Non-diagonal supercells',
*    .                   ' are not implemented yet'
*             stop 'redata: ERROR: Non-diagonal supercell.'
*           endif
*         enddo
*       enddo
        call idiag( 3, nscell, nsd, ml, mr, maux )
        do i = 1,3
          do ix = 1,3
            dcell(ix,i) = ( cell(ix,1) * mr(1,i) +
     .                      cell(ix,2) * mr(2,i) +
     .                      cell(ix,3) * mr(3,i) ) / nsd(i,i)
          enddo
        enddo
        na = 0
        do ia = 1,nua
          na = na + 1
          indxua(na) = ia
        enddo
        do i3 = 0,nsd(3,3)-1
        do i2 = 0,nsd(2,2)-1
        do i1 = 0,nsd(1,1)-1
          if (i1.ne.0 .or. i2.ne.0 .or. i3.ne.0) then
            do ia = 1,nua
              na = na + 1
              isa(na) = isa(ia)
              indxua(na) = ia
              do ix = 1,3
                xa(ix,na) = xa(ix,ia) + dcell(ix,1) * i1 +
     .                                  dcell(ix,2) * i2 +
     .                                  dcell(ix,3) * i3
              enddo
            enddo
          endif
        enddo
        enddo
        enddo
      endif
C ..................

      write(6,102)
      write(6,*) ' '

100   format(/,'redata: ',71(1h*))
101   format('redata:                  INPUT ERROR')
102   format('redata: ',71(1h*))
103   format('redata: ',i4,2x,3f10.5,i3) 
      end

