C $Id: redata.f,v 1.45 1999/05/25 13:56:03 ordejon Exp $

      subroutine redata(maxa, maxspn, na, ns, nspin, overflow,
     .                  slabel, sname, isa, cell, xa, outlng, g2max,
     .                  charnet, negl, nscf, dDtol, mix, wmix, isolve, 
     .                  temp, fixspin, ts, ncgmax, ftol, strtol, eta, 
     .                  etol, rcoor, 
     .                  ioptlwf, chebef, noeta, rcoorcp, beta, pmax,
     .                  idyn, istart, ifinal, nmove, ianneal, iquench,
     .                  dt, ia1, ia2, dx, dxmax, tt, tp, mn, mpr, 
     .                  bulkm, taurelax, 
     .                  writedim, usesavelwf, usesavedm, usesavecg,
     .                  mullipop, inspn, maxsav, nkick, wmixkick, 
     .                  pulfile, tempinit, dumpcharge, varcel )
C *********************************************************************
C Subroutine to read the data for the SIESTA program
C
C     It uses the FDF (Flexible Data Fromat) package 
C     of J.M.Soler and A.Garcia
C
C Writen by P.Ordejon, December'96
C ***************************** INPUT *********************************
C integer maxa             : Maximum number of atoms
C integer maxspn           : Max. number of spin components (1,2 or 4)
C **************************** OUTPUT *********************************
C integer na               : Number of atoms
C integer ns               : Number of species
C integer nspin            : Spin polarization
C logical overflow         : True = Some of the dimensions is too small
C character*20 slabel      : System Label (to name output files)
C character*59 sname       : System Name
C integer isa(maxa)        : Species index of each atom
C real*8 cell(3,3)         : (Super) lattice vectors CELL(ixyz,ivector)
C                            (in Bohr)
C real*8 xa(3,maxa)        : Atomic coordinates (Bohr)
C real*8 charnet           : Net charge (in units of |e|)
C logical outlng           : Long (true) or short (false) output
C real*8 g2max             : PW cutoff energy (Ry)
C logical negl             : True = Neglect interactions between
C                            non-overlaping orbitals (coming from
C                            KB projectors)
C integer nscf             : Maximum number of SCF cycles per time step
C real*8 dDtol             : Maximum Density Matrix tolerance
C logical mix              : Perform mix in first SCF step
C real*8 wmix              : Amount of output DM for new DM
C integer isolve           : Method of solution.  0 = Diagonalization
C                                                 1 = Order-N
C real*8 temp              : Temperature for Fermi smearing (Ry)
C logical fixspin          : Fix the spin of the system?
C real*8  ts               : Total spin of the system
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
C logical chebef          : Compute the chemical potential 
C logical noeta            : Use computed Chem.pot. instead of eta
C real*8 rcoorcp           : Cutoff (Bohr) to compute the chem.pot.
C real*8 beta              : Inverse temperature to compute chem.pot.
C integer pmax             : Order of Chebi expansion for chem.pot.
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
C real*8 strtol            : Maximum stress for CG structure optimization
C integer ianneal          : Annealing option for idyn = 5
C                             1 = Temperature 
C                             2 = Pressure
C                             3 = Temperature and Pressure
C integer iquench          : Quench option: 0 = No;  1 = Yes
C real*8 dt                : Length of time step (fs)
C real*8 dx                : Atomic displacement for Force Constants
C                             calculation
C integer ia1              : First atom to displace for force constants
C integer ia2              : Last atom to displace for force constants
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
C logical usesavecg        : True = try to use continuation CG files
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
C integer nkick            : Perform a linear mixing eack nkick scf cycles
C real*8 wmixkick          : Mixing parameter for linear mixing each nkick scf
C                            cycles
C logical pulfile          : Use file (.true.) or memory (.false.)
C                            to store Pulay miximg intermediate vectors
C                            Default: .false.
C real*8 tempinit          : Initial temperature (Kelvin) of the MD simulation
C logical dumpcharge       : True: Dump information to plot charge contours
C                            by the external DENCHAR application program.
C logical varcel           : variable shape for CG optimization or dynamics
C **********************************************************************

      implicit none

      integer
     .  maxa, maxspn

      character
     .  slabel*20, sname*59

      integer
     .  ia1, ia2, ianneal, idyn, ifinal, ioptlwf,
     .  iquench, isa(maxa), isolve, istart, maxsav,
     .  mullipop, na, ncgmax, nkick, nmove, ns, nscf, nspin,
     .  pmax

      double precision
     .  beta, bulkm, cell(3,3), charnet,
     .  dDtol, dt, dx, dxmax, eta, etol, ftol, g2max,
     .  mn, mpr, rcoor, rcoorcp, strtol,
     .  taurelax, temp, tempinit, tcp, tp, ts, tt, wmix, wmixkick,
     .  xa(3,maxa)

      logical
     .  chebef, dumpcharge, fixspin, inspn, mix, negl, noeta, overflow, 
     .  outlng, pulfile, usesavecg, usesavelwf, usesavedm, varcel, 
     .  writedim

C Internal parameters ................................................
C na_diag      : maximum number of atoms with diagon as default method
C g2max_default : Mesh cutoff default, in Ry
C temp_default  : Electronic temperature default, in Ry
      integer na_diag
      double precision g2max_default, temp_default
      parameter ( na_diag       = 100      )
      parameter ( g2max_default = 50.d0    )
      parameter ( temp_default  = 1.900d-3 )
C ................

C  Internal variables .................................................
      character
     .  annop*22, dyntyp*22, 
     .  filein*20, fileout*20, 
     .  method*6, line*150, 
     .  lwfopt*13

      character
     .  annop_default*22, dyntyp_default*22, 
     .  lwfopt_default*13,  slabel_default*59, 
     .  sname_default*20

      integer 
     .  ia1_default, ia2_default,
     .  ifinal_default, istart_default, maxsv_default, mpop_default,
     .  na_default, ncgmax_default, nkick_default, nmove_default,
     .  ns_default, nscf_default, pmax_default

      integer 
     .  i, ic, ix, mscell(3,3), ncells

      integer length, lun

      double precision
     .  alat, ucell(3,3), volcel

      double precision
     .  bulkm_default, cnet_default,
     .  dDtol_default, dt_default, dx_default, dxmax_default,
     .  eta_default, etol_default, ftol_default,
     .  mn_default, mpr_default, rccp_default, rcoor_default, 
     .  taurelax_default, tcp_default, ti_default, tp_default, 
     .  tt_default, wmix_default, wmixk_default

      logical
     .  leqi, noncol, qnch, sppol

      logical
     .  chebef_default, dc_default, fs_default, inspn_default, 
     .  negl_default, mix_default, noeta_default, pul_default, 
     .  qnch_default, usdm_default, uslwf_default,  wd_default, 
     .  uscg_default

C ................

C Define FDF calls ....................................................
      include 'fdf/fdfdefs.h'
C ................

      overflow = .false.
      na = 1
      ns = 1
      nspin = 1

C Print Welcome and Presentation .......................................
      write(6,'(/a)') 
     . '                           ***********************       '
      write(6,'(a)') 
     . '                           *  WELCOME TO SIESTA  *       '
      write(6,'(a)')
     . '                           ***********************       '
C ..................

C Dump data file to output file .......................................
C and generate scratch file for FDF to read from
C
      write(6,'(/,a,18(1h*),a,28(1h*))')
     .  'redata: ', ' Dump of input data file '

      call io_assign(lun)
      open(lun,file='FDF_STDIN',form='formatted',status='unknown')
      rewind(lun)

 10   continue
         read(5,err=20,end=20,fmt='(a)') line
         call chrlen(line,0,length)
         write(lun,'(a)') line(1:length)
         if (length .ne. 0) write(6,'(a)') line(1:length)
         goto 10
 20   continue
      call io_close(lun)

      write(6,'(a,18(1h*),a,29(1h*))')
     .  'redata: ', ' End of input data file '
C ..................

C Read data from FDF file..............................................

C Set up fdf ...
      filein = 'FDF_STDIN'
      fileout = 'out.fdf'
      call fdf_init(filein,fileout)
C ...

      write(6,'(/,a,18(1h*),a,30(1h*))')
     .  'redata: ', ' Simulation parameters '
      write(6,'(a)')  'redata:'
      write(6,'(a)')  'redata: The following are some of the '//
     .                         'parameters of the simulation.'
      write(6,'(a)')  'redata: A complete list of the parameters '//
     .                         'used, including defect values,'
      write(6,'(a,a)')'redata: can be found in file ',fileout
      write(6,'(a)')  'redata:'

C Defile Name of the system ...
      sname_default = ' '
      sname = fdf_string('SystemName',sname_default)
      write(6,'(a,71(1h-))') 'redata: '
      write(6,'(a,a)') 
     . 'redata: System Name: ',sname
      write(6,'(a,71(1h-))') 'redata: '
C ...

C Defile System Label (short name to label files) ...
      slabel_default  = 'siesta'
      slabel = fdf_string('SystemLabel',slabel_default)
      write(6,'(a,4x,a)') 
     . 'redata: System Label                     = ',slabel
C ...

C Type of output
      outlng = fdf_boolean( 'LongOutput', .false. )
      write(6,'(a,4x,l1)')
     . 'redata: Long output                      = ',outlng

C Read Number of Atoms ...
      na_default = 0
      na = fdf_integer('NumberOfAtoms',na_default)
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
      call chkdime(maxa, na, overflow, na)
C ...

C Defile Number of species ...
      ns_default = 0
      ns = fdf_integer('NumberOfSpecies',ns_default)
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

C Dump minimum sizes to siesta.h and stop
      wd_default = .false.
      writedim  = fdf_boolean('WriteSiestaDim',wd_default)
      write(6,'(a,4x,l1)') 
     . 'redata: Write minimum siesta dimensions  = ',writedim
C ...

C Dump information to plot charge contours
C by the external DENCHAR application program.
      dc_default = .false.
      dumpcharge  = fdf_boolean('WriteDenchar',dc_default)
      write(6,'(a,4x,l1)') 
     . 'redata: Dump information for DENCHAR     = ',dumpcharge
C ...

C Perform Mulliken Population Analisys
      mpop_default = 0
      if ( outlng ) mpop_default = 1
      mullipop = fdf_integer('WriteMullikenPop',mpop_default)
      if (mullipop .eq. 0) then
      write(6,'(a)') 
     . 'redata: Write Mulliken Pop.              =     NO'
      elseif (mullipop .eq. 1) then
      write(6,'(a,a)') 
     . 'redata: Write Mulliken Pop.              =     ',
     . 'Atomic and Orbital charges'
      elseif (mullipop .eq. 2) then
      write(6,'(a,a/45x,a)') 
     . 'redata: Write Mulliken Pop.              =     ',
     . 'Atomic and Orbital charges','plus Atomic Overlap Populations'
      elseif (mullipop .eq. 3) then
      write(6,'(a,a/45x,a/45x,a)') 
     . 'redata: Write Mulliken Pop.              =     ',
     . 'Atomic and Orbital charges','plus Atomic Overlap Populations',
     . 'plus Oorbital Overlap Populations'
      else
        stop 'redata: Wrong value for WriteMullikenPop'
      endif
C ...

C Lattice constant and lattice vectors .........
      call redcel( alat, ucell, cell, mscell )
      if (alat .ne. 0.d0) then
        write(6,'(a,f10.4,a)')
     . 'redata: Lattice Constant                 = ',alat,'  Bohr'
          write(6,'(a,/,(a,1x,3f12.6))')
     .     'redata: Lattice vectors (in units of Lattice Constant)',
     .     ('redata:',(cell(ix,i)/alat,ix=1,3),i=1,3)
      endif
      if (volcel(ucell) .lt. 1.d-8) then
        ncells = 1
      else
        ncells = nint( volcel(cell) / volcel(ucell) )
        if (ncells.gt.1) then
          write(6,'(a,/,(8x,3f12.6))')
     .      'redata: Total-cell (supercell) vectors (Bohr)', cell
          write(6,'(a,i6)') 'redata: Number of unit cells  =', ncells
          write(6,'(a,i6)') 'redata: Total number of atoms =', na*ncells
        else
          write(6,'(a,/,(a,1x,3f12.6))')
     .     'redata: Lattice vectors (in Bohr)', 
     .     ('redata:',(cell(ix,i),ix=1,3),i=1,3)
        endif
      endif
C ..................

C Spin Polarization ...
      sppol  = fdf_boolean('SpinPolarized',.false.)
      noncol = fdf_boolean('NonCollinearSpin',.false.)
      if (noncol) then
        nspin = 4
      elseif (sppol) then
        nspin = 2
      else 
        nspin = 1
      endif
      write(6,'(a,4x,i1)') 
     . 'redata: Number of spin components        = ',nspin
      call chkdime(maxspn, nspin, overflow, nspin)
C ...
 
C Planewave cutoff of the real space mesh ...
      g2max = fdf_physical('MeshCutoff',g2max_default,'Ry')
      write(6,'(a,f10.4,a)') 
     . 'redata: Mesh Cutoff                      = ',g2max,'  Ry'
C ...

C Net charge in the cell ...
      cnet_default = 0.0d0
      charnet=fdf_double('NetCharge',cnet_default)
      write(6,'(a,f10.4,a)') 
     . 'redata: Net charge of the system         = ',charnet,' |e|'
C ...
       
C SCF Loop parameters ...
C     Maximum number of SCF iterations
      nscf_default = 50
      nscf = fdf_integer('MaxSCFIterations',nscf_default)
      write(6,'(a,i5)') 
     . 'redata: Max. number of SCF Iter          = ',nscf

C     Pulay mixing, numer of iterations for one Pulay mixing (maxsav)
      maxsv_default = 0
      maxsav = fdf_integer('DM.NumberPulay',maxsv_default)
      if(maxsav .gt. 1) then
        write(6,'(a,i5,a)') 
     .   'redata: Pulay mixing using               = ',maxsav,
     .   ' iterations'
      else
         write(6,'(a)')'redata: Mixing is linear'
      endif
 
C     Mix density matrix on first SCF step
C     (mix)
      mix_default = .false.
      mix  = fdf_boolean('DM.MixSCF1',mix_default)
      write(6,'(a,4x,l1)')
     .  'redata: Mix DM in first SCF step ?       = ',mix

C     Use disk or memory to store intermediate Pulay miximg vectors
C     (pulfile)
      pul_default = .false.
      pulfile  = fdf_boolean('DM.PulayOnFile',pul_default)
      write(6,'(a,4x,l1)')
     . 'redata: Write Pulay info on disk         = ',pulfile
C ...


C     Density Matrix Mixing  (proportion of output DM in new input DM)
      wmix_default = 0.25d0
      wmix = fdf_double('DM.MixingWeight',wmix_default)
      write(6,'(a,f10.4,a)') 
     . 'redata: New DM Mixing Weight             = ',wmix

C     Perform linear mixing each nkick SCF iterations (to kick system
C     when it is pinned in a poorly convergent SCF loop)
      nkick_default = 0
      nkick = fdf_integer('DM.NumberKick',nkick_default)
      if(nkick .ge. 1) then
        write(6,'(a,i5,a)')
     .   'redata: Kick with linear mixing every    = ',nkick,
     .   ' iterations'
      else
         write(6,'(a)')'redata: No kicks to SCF'
      endif
 
C     Density Matrix Mixing each nkick SCF iterations
      wmixk_default = 0.50d0
      wmixkick = fdf_double('DM.KickMixingWeight',wmixk_default)
      write(6,'(a,f10.4,a)')
     . 'redata: DM Mixing Weight for Kicks       = ',wmixkick

C     Density Matrix Tolerance for achieving Self-Consistency
      dDtol_default = 1.d-4
      dDtol = fdf_double('DM.Tolerance',dDtol_default)
      write(6,'(a,f12.6,a)') 
     . 'redata: DM Tolerance for SCF             = ',dDtol

C     Initial spin density: Maximum polarization, Ferro (false), AF (true)
      if(sppol) then
        inspn_default = .false.
        inspn = fdf_boolean('DM.InitSpinAF',inspn_default)
        write(6,'(a,4x,l1)')
     .   'redata: Antiferro initial spin density   = ',inspn
      endif
C ...

C Use continumation DM files
*     usdm_default = .true.
      usdm_default = fdf_boolean('UseSaveData',.false.)
      usesavedm  = fdf_boolean('DM.UseSaveDM',usdm_default)
      write(6,'(a,4x,l1)') 
     . 'redata: Use continuation files for DM    = ',
     . usesavedm
C ...

C Neglect Interactions between non-overlapping orbitals ...
      negl_default = .false.
      negl  = fdf_boolean('NeglNonOverlapInt',negl_default)
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
     .   '    Diagonalization'
      else if (leqi(method,'ordern')) then
        isolve = 1
        write(6,'(a,a)') 
     .   'redata: Method of Calculation            = ',
     .   '    Order-N'
        if (nspin .gt. 1) then
         write(6,100) 
         write(6,101) 
         write(6,'(a)') 
     .   'redata:    You chose the Order-N solution option'
         write(6,'(a)') 
     .   'redata:    together with nspin>1.  This is not  '
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
      temp = fdf_physical('ElectronicTemperature',temp_default,'Ry')
      if (isolve .eq. 0) then
        write(6,'(a,f10.4,a)') 
     .  'redata: Electronic Temperature           = ',temp,'  Ry'
      endif
C ...

C Fix the spin of the system to a given value 
      fs_default = .false.
      fixspin  = fdf_boolean('FixSpin',fs_default)
      if (fixspin .and. nspin .ne. 2) then
        write(6,'(a)') 
     . 'redata: ERROR: You can only fix the spin of the system' 
        write(6,'(a)') 
     . 'redata:        for collinear spin polarized calculations.'
        write(6,102)
        stop
      endif
      write(6,'(a,4x,l1)') 
     . 'redata: Fix the spin of the system       = ',fixspin 
C ...

C Value of the Spin of the system (only used if fixspin = TRUE
      if (fixspin) then
        ts = fdf_double('TotalSpin',0.0d0)
        write(6,'(a,f10.4)') 
     .   'redata: Value of the Spin of the System  = ',ts
      else
        ts = 0.0
      endif
C ...
        

C Order-N solution parameters ...
C     Maximum number of CG minimization iterations
      ncgmax_default = 1000
      ncgmax = fdf_integer('ON.MaxNumIter',ncgmax_default)

C     Relative tolerance in total band structure energy
      etol_default = 1.d-8
      etol = fdf_double('ON.etol',etol_default)

C     Fermi level parameter
      eta_default = 0.d0
      eta = fdf_physical('ON.eta',eta_default,'Ry')

C     Cutoff radius for Localized Wave Functions
      rcoor_default = 9.5d0
      rcoor = fdf_physical('On.RcLWF',rcoor_default,'Bohr')

C     Use continumation LWF files
*     uslwf_default = .true.
      uslwf_default = fdf_boolean('UseSaveData',.false.)
      usesavelwf  = fdf_boolean('ON.UseSaveLWF',uslwf_default)

C     Option on how to build LWF's (disk or functionals)
      lwfopt_default = 'kim'
      lwfopt  = fdf_string('ON.functional',lwfopt_default)
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

C     Option to calculate the Chemical potential in O(N)
      chebef_default = .false.
      chebef = fdf_boolean('ON.ChemicalPotential',chebef_default)
      
C     Option to use the Chemical Potential calculated instead
C     of the eta variable of the input
      noeta_default = .false.
      noeta = fdf_boolean('ON.ChemicalPotentialUse',noeta_default)

      if (noeta) chebef = .true.

C     Cutoff radius to calculate the Chemical Potential by projection
      rccp_default = 9.5d0
      rcoorcp = fdf_physical('ON.ChemicalPotentialRc',
     .                        rccp_default,'Bohr')

C     Temperature of the Fermi distribution to calculate the
C     Chemical potential by projection
      tcp_default = 0.05
      tcp = fdf_physical('ON.ChemicalPotentialTemperature',
     .                   tcp_default,'Ry')
      beta = 1.0d0/tcp

C     Order of the Chebishev expansion to calculate the Chemical
C     potential
      pmax_default = 100
      pmax = fdf_integer('ON.ChemicalPotentialOrder',pmax_default)
C ...


      if (isolve .eq. 1) then
        write(6,'(a,i5)') 
     .  'redata: Maximum number of iterations     = ',ncgmax
        write(6,'(a,d12.2)') 
     .  'redata: Relative tolerance               = ',etol
        write(6,'(a,f10.4,a)') 
     .  'redata: Eta (Fermi level parameter)      = ',eta,'  Ry'
        write(6,'(a,f10.4,a)') 
     .  'redata: Radius of LWFs                   = ',rcoor,
     .  '  Bohr'
        write(6,'(a,4x,l1)') 
     .  'redata: Use continuation files for LWF   = ',
     .  usesavelwf
        write(6,'(a,a)') 
     .  'redata: Method to build LWFs             =     ',lwfopt

        if (chebef) then
        write(6,'(a,l1)')
     .  'redata: Compute Chemical Potential       =     ',chebef
        write(6,'(a)')
     .  'redata: Use the calculated Chemical ..'
        write(6,'(a,l1)')
     .  'redata: ..Potential instead of eta       =     ',noeta
        write(6,'(a,f10.4,a)') 
     .  'redata: Radius to compute the Chem. Pot. = ',rcoorcp,
     .  '  Bohr'
        write(6,'(a)')
     .  'redata: Temp. for Fermi distribution ..'
        write(6,'(a,f10.4,a)') 
     .  'redata: .. to compute the Chem. Pot.     = ',tcp,
     .  '    Ry'
        write(6,'(a,i5)') 
     .  'redata: Order of the Chebishev expansion = ',pmax
        endif
        
      endif
C ...

C Dynamics parameters ...
      varcel = fdf_boolean('MD.VariableCell', .false. )
C     Type of dynamics 
      dyntyp_default = 'verlet'   
      dyntyp = fdf_string('MD.TypeOfRun',dyntyp_default)
      if (leqi(dyntyp,'cg')) then
        idyn = 0
        write(6,'(a,a)') 
     .   'redata: Dynamics option                  = ',
     .   '    CG coord. optimization'
        write(6,'(a,4x,l1)')
     .   'redata: Variable cell                    = ', varcel
        uscg_default = fdf_boolean('UseSaveData',.false.)
        usesavecg  = fdf_boolean('MD.UseSaveCG',uscg_default)
        write(6,'(a,4x,l1)')
     .   'redata: Use continuation files for CG    = ',
     .   usesavecg
      else if (leqi(dyntyp,'verlet')) then
        idyn = 1
        write(6,'(a,a)') 
     .   'redata: Dynamics option                  = ',
     .   '    Verlet MD run'
      else if (leqi(dyntyp,'nose')) then
        idyn = 2
        write(6,'(a,a)') 
     .   'redata: Dynamics option                  = ',
     .   '    Nose termostat MD run'
      else if (leqi(dyntyp,'parrinellorahman')) then
        idyn = 3
        write(6,'(a,a)') 
     .   'redata: Dynamics option                  = ',
     .   '    Parrinello-Rahman MD run'
      else if (leqi(dyntyp,'noseparrinellorahman')) then
        idyn = 4
        write(6,'(a,a)') 
     .   'redata: Dynamics option                  = ',
     .   '    Nose-Parrinello-Rahman MD run'
      else if (leqi(dyntyp,'anneal')) then
        idyn = 5
        write(6,'(a,a)') 
     .   'redata: Dynamics option                  = ',
     .   '    Annealing MD run'
      else if (leqi(dyntyp,'fc')) then
        idyn = 6
        write(6,'(a,a)')
     .   'redata: Dynamics option                  = ',
     .   '    Force Constants Matrix Calculation'
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
        write(6,'(a)') 'redata:      - FC                             '
        write(6,102)
        stop
      endif 


C     Maximum number of steps in CG coordinate optimization
      nmove_default = 0
      nmove = fdf_integer('MD.NumCGsteps',nmove_default)

C     Maximum atomic displacement in one CG step
      dxmax_default = 0.2d0
      dxmax = fdf_physical('MD.MaxCGDispl',dxmax_default,'Bohr')

C     Tolerance in the maximum atomic force (def 0.04 eV/Ang)
      ftol_default = 0.00155574d0 
      ftol = fdf_physical('MD.MaxForceTol',ftol_default,'Ry/Bohr')

C     Tolerance in the maximum residual stress (var cell) def = 1 GPa 
      strtol = fdf_physical('MD.MaxStressTol',6.79773d-5,'Ry/Bohr**3')

      if (idyn .eq. 0) then
        write(6,'(a,i5)') 
     .  'redata: Maximum number of CG moves       = ',nmove
        write(6,'(a,f10.4,a)') 
     .  'redata: Max atomic displ per move        = ',dxmax,'  Bohr'
        write(6,'(a,f10.4,a)') 
     .  'redata: Force tolerance                  = ',ftol,'  Ry/Bohr'
        if ( varcel ) then
           strtol = dabs(strtol)
           write(6,'(a,f10.4,a)')
     .  'redata: Stress tolerance                 = ',strtol/6.79773d-5,
     .                                              '  GPa'
        endif

      endif
  
C     Initial time step for MD
      istart_default = 1
      istart = fdf_integer('MD.InitialTimeStep',istart_default)

C     Final time step for MD
      ifinal_default = 1
      ifinal = fdf_integer('MD.FinalTimeStep',ifinal_default)

C     Length of time step for MD
      dt_default = 1.d0
      dt = fdf_physical('MD.LengthTimeStep',dt_default,'fs')

C     Quench Option
      qnch_default = .false.
      qnch = fdf_boolean('MD.Quench',qnch_default)
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
      ti_default = 0.d0
      tempinit = fdf_physical('MD.InitialTemperature',ti_default,'K')

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
      tt_default = 0.d0
      tt = fdf_physical('MD.TargetTemperature',tt_default,'K')
      
C     Target Pressure
      tp_default = 0.d0
      tp = fdf_physical('MD.TargetPressure',tp_default,'Ry/Bohr**3')

C     Mass of Nose variable
      mn_default = 1.d2
      mn = fdf_physical('MD.NoseMass',mn_default,'Ry*fs**2')

C     Mass of Parrinello-Rahman variables
      mpr_default = 1.d2
      mpr = fdf_physical('MD.ParrinelloRahmanMass',
     .                    mpr_default,'Ry*fs**2')

      if (idyn .eq. 2 .or. idyn .eq. 4) then
        write(6,'(a,f10.4,a)') 
     .  'redata: Nose mass                        = ',mn,'  Ry/fs**2'
      endif

      if (idyn .eq. 3 .or. idyn .eq. 4) then
        write(6,'(a,f10.4,a)') 
     .  'redata: Parrinello-Rahman mass           = ',mpr,'  Ry/fs**2'
      endif

C     Annealing option
      annop_default = 'TemperatureAndPressure'
      annop = fdf_string('MD.AnnealOption',annop_default)
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
     .  'redata: Target Pressure                  = ',tp,'  Ry/Bohr**3'
      endif

C     Relaxation Time for Annealing
      taurelax_default = 1.d2
      taurelax = fdf_physical('MD.TauRelax',taurelax_default,'fs')
      if (idyn .eq. 5) then
        write(6,'(a,f10.4,a)') 
     .  'redata: Annealing Relaxation Time        = ',
     .   taurelax,'  fs'
      endif
        
C     Estimated Bulk modulus (for Pressure annealing)
      bulkm_default = 1.d2
      bulkm = fdf_physical('MD.BulkModulus',bulkm_default,'Ry/Bohr**3')
      if (idyn .eq. 5 .and. (ianneal .eq. 2 .or. ianneal .eq. 3)) then
        write(6,'(a,f10.4,a)') 
     .  'redata: Approx. Bulk Modulus             = ',
     .   bulkm,'  Ry/Bohr**3'
      endif

C     Atomic displacement for force constant calculation
      dx_default = 0.04d0
      dx = fdf_physical('MD.FCDispl',dx_default,'Bohr')

C     First and last atoms to displace for calculation of force constants
      ia1_default = 1
      ia1 = fdf_integer('MD.FCfirst',ia1_default)
      ia2_default = na
      ia2 = fdf_integer('MD.FClast',ia2_default)

      if (idyn .eq. 6) then
        write(6,'(a,f10.4,a)')
     .  'redata: Atomic displ for force constans  = ',dx,'  Bohr'
        write(6,'(a,i8)')
     .  'redata: First atom to move               = ',ia1
        write(6,'(a,i8)')
     .  'redata: Last atom to move                = ',ia1
      endif
C ...

C Variable cell shape? Depending on input and type of dynamics
      varcel = varcel .or. (idyn.eq.3) .or. (idyn.eq.4) 
     .                .or. (idyn.eq.5 .and. ianneal.ne.1)
      varcel = varcel .and. (idyn.ne.1) .and. (idyn.ne.2) 
     .                .and. (idyn.ne.6)
     .                .and. (.not. (idyn.eq.5 .and. ianneal.ne.1) )

C Atomic Coordinates ..................................
      call recoor( maxa, na, isa, xa )
      call chkdime(maxa, na, overflow, na)
      
      write(6,102)

100   format(/,'redata: ',71(1h*))
101   format('redata:                  INPUT ERROR')
102   format('redata: ',71(1h*))
103   format('redata: ',i4,2x,3f10.5,i3) 
      end

