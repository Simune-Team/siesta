      subroutine redata(na, ncells, ns, nspin, isa, xa, outlng, g2max,
     .                  charnet, negl, nscf, dDtol, mix, wmix, isolve, 
     .                  temp, fixspin, ts, ncgmax, ftol, strtol, eta, 
     .                  etol, rcoor, 
     .                  ioptlwf, chebef, noeta, rcoorcp, beta, pmax,
     .                  idyn, istart, ifinal, nmove, ianneal, iquench,
     .                  dt, ia1, ia2, dx, dxmax, tt, tp, mn, mpr, 
     .                  bulkm, taurelax, 
     .                  usesavelwf, usesavedm, usesavecg,
     .                  mullipop, inspn, maxsav, nkick, wmixkick, 
     .                  pulfile, tempinit, dumpcharge, varcel )
C *********************************************************************
C Subroutine to read the data for the SIESTA program
C
C     It uses the FDF (Flexible Data Fromat) package 
C     of J.M.Soler and A.Garcia
C
C Writen by P.Ordejon, December'96
C Modified for introduction of dynamic memory in SIESTA by JDG Sept 99
C ***************************** INPUT *********************************
C integer na               : Number of atoms
C integer ncells           : Number of unit cells in supercell
C integer nspin            : Spin polarization
C **************************** OUTPUT *********************************
C integer ns               : Number of species
C integer isa(na)          : Species index of each atom
C real*8 xa(3,na)          : Atomic coordinates (Bohr)
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

C
C  Modules
C
      use precision
      use fdf
#ifdef MPI
      use mpi
#endif
      use diagmemory

      implicit none

      integer
     .  ia1, ia2, ianneal, idyn, ifinal, ioptlwf,
     .  iquench, na, ncells, isa(na), isolve, istart, maxsav,
     .  mullipop, ncgmax, nkick, nmove, ns, nscf, nspin,
     .  pmax

      double precision
     .  beta, bulkm, charnet,
     .  dDtol, dt, dx, dxmax, eta, etol, ftol, g2max,
     .  mn, mpr, rcoor, rcoorcp, strtol,
     .  taurelax, temp, tempinit, tcp, tp, ts, tt, wmix, wmixkick,
     .  xa(3,na)

      logical
     .  chebef, dumpcharge, fixspin, inspn, mix, negl, noeta, 
     .  outlng, pulfile, usesavecg, usesavelwf, usesavedm, varcel

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
     .  method*6, lwfopt*13

      character
     .  annop_default*22, dyntyp_default*22, 
     .  lwfopt_default*13

      integer 
     .  ia1_default, ia2_default,
     .  ifinal_default, istart_default, maxsv_default, mpop_default,
     .  ncgmax_default, nkick_default, nmove_default,
     .  ns_default, nscf_default, pmax_default

      integer 
     .  Node

#ifdef MPI
      integer
     .  MPIerror
#endif

      double precision
     .  bulkm_default, cnet_default,
     .  dDtol_default, dt_default, dx_default, dxmax_default,
     .  eta_default, etol_default, ftol_default,
     .  mn_default, mpr_default, rccp_default, rcoor_default, 
     .  taurelax_default, tcp_default, ti_default, tp_default, 
     .  tt_default, wmix_default, wmixk_default

      logical
     .  leqi, qnch

      logical
     .  chebef_default, dc_default, fs_default, inspn_default, 
     .  negl_default, mix_default, noeta_default, pul_default, 
     .  qnch_default, usdm_default, uslwf_default, uscg_default

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

      ns = 1

C Type of output
      if (Node.eq.0) then
        outlng = fdf_boolean( 'LongOutput', .false. )
        write(6,'(a,4x,l1)')
     .   'redata: Long output                      = ',outlng
      endif
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(outlng,1,MPI_logical,0,MPI_Comm_World,
     .  MPIerror)
#endif

C Defile Number of species ...
      ns_default = 0
      if (Node.eq.0) then
        ns = fdf_integer('NumberOfSpecies',ns_default)
        write(6,'(a,3x,i2)') 
     .   'redata: Number of Atomic Species         = ',ns
      endif
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(ns,1,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
#endif
      if (ns .le. 0) then
        if (Node.eq.0) then
          write(6,100) 
          write(6,101) 
          write(6,'(a)') 
     .   'redata: ERROR: Number of species must be larger than zero.'
          write(6,102)
        endif
        stop
      endif

C ...

C Dump information to plot charge contours
C by the external DENCHAR application program.
      dc_default = .false.
      if (Node.eq.0) then
        dumpcharge  = fdf_boolean('WriteDenchar',dc_default)
        write(6,'(a,4x,l1)') 
     .   'redata: Dump information for DENCHAR     = ',dumpcharge
      endif
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(dumpcharge,1,MPI_logical,0,MPI_Comm_World,
     .  MPIerror)
#endif
C ...

C Perform Mulliken Population Analisys
      mpop_default = 0
      if ( outlng ) mpop_default = 1
      if (Node.eq.0) then
        mullipop = fdf_integer('WriteMullikenPop',mpop_default)
      endif
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(mullipop,1,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
#endif
      if (Node.eq.0) then
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
      else
        if (mullipop.lt.0.or.mullipop.gt.3) then
          stop
        endif
      endif
C ..................

C Planewave cutoff of the real space mesh ...
      if (Node.eq.0) then
        g2max = fdf_physical('MeshCutoff',g2max_default,'Ry')
        write(6,'(a,f10.4,a)') 
     .   'redata: Mesh Cutoff                      = ',g2max,'  Ry'
      endif
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(g2max,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
C ...

C Net charge in the cell ...
      cnet_default = 0.0d0
      if (Node.eq.0) then
        charnet=fdf_double('NetCharge',cnet_default)
        write(6,'(a,f10.4,a)') 
     .   'redata: Net charge of the system         = ',charnet,' |e|'
      endif
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(charnet,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
C ...
       
C SCF Loop parameters ...
C     Maximum number of SCF iterations
      nscf_default = 50
      if (Node.eq.0) then
        nscf = fdf_integer('MaxSCFIterations',nscf_default)
        write(6,'(a,i5)') 
     .   'redata: Max. number of SCF Iter          = ',nscf
      endif
#ifdef MPI
      call MPI_Bcast(nscf,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

C     Pulay mixing, numer of iterations for one Pulay mixing (maxsav)
      maxsv_default = 0
      if (Node.eq.0) then
        maxsav = fdf_integer('DM.NumberPulay',maxsv_default)
        if (maxsav .gt. 1) then
          write(6,'(a,i5,a)') 
     .     'redata: One Pulay mixing every           = ',maxsav,
     .     ' iterations'
        else
           write(6,'(a)')'redata: Mixing is linear'
        endif
      endif
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(maxsav,1,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
#endif
 
C     Mix density matrix on first SCF step
C     (mix)
      mix_default = .false.
      if (Node.eq.0) then
        mix  = fdf_boolean('DM.MixSCF1',mix_default)
        write(6,'(a,4x,l1)')
     .    'redata: Mix DM in first SCF step ?       = ',mix
      endif
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(mix,1,MPI_logical,0,MPI_Comm_World,
     .  MPIerror)
#endif

C     Use disk or memory to store intermediate Pulay miximg vectors
C     (pulfile)
      pul_default = .false.
      if (Node.eq.0) then
        pulfile  = fdf_boolean('DM.PulayOnFile',pul_default)
        write(6,'(a,4x,l1)')
     .    'redata: Write Pulay info on disk?       = ',pulfile
      endif
#ifdef MPI
      call MPI_Bcast(pulfile,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
C ...


C     Density Matrix Mixing  (proportion of output DM in new input DM)
      wmix_default = 0.25d0
      if (Node.eq.0) then
        wmix = fdf_double('DM.MixingWeight',wmix_default)
        write(6,'(a,f10.4,a)') 
     .   'redata: New DM Mixing Weight             = ',wmix
      endif
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(wmix,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

C     Perform linear mixing each nkick SCF iterations (to kick system
C     when it is pinned in a poorly convergent SCF loop)
      nkick_default = 0
      if (Node.eq.0) then
        nkick = fdf_integer('DM.NumberKick',nkick_default)
        if (nkick .ge. 1) then
          write(6,'(a,i5,a)')
     .     'redata: Kick with linear mixing every    = ',nkick,
     .     ' iterations'
        else
           write(6,'(a)')'redata: No kicks to SCF'
        endif
      endif
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(nkick,1,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
#endif
 
C     Density Matrix Mixing each nkick SCF iterations
      wmixk_default = 0.50d0
      if (Node.eq.0) then
        wmixkick = fdf_double('DM.KickMixingWeight',wmixk_default)
        write(6,'(a,f10.4,a)')
     .   'redata: DM Mixing Weight for Kicks       = ',wmixkick
      endif
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(wmixkick,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

C     Density Matrix Tolerance for achieving Self-Consistency
      dDtol_default = 1.d-4
      if (Node.eq.0) then
        dDtol = fdf_double('DM.Tolerance',dDtol_default)
        write(6,'(a,f12.6,a)') 
     .   'redata: DM Tolerance for SCF             = ',dDtol
      endif
#ifdef MPI
      call MPI_Bcast(dDtol,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

C     Initial spin density: Maximum polarization, Ferro (false), AF (true)
      if (nspin.eq.2) then
        inspn_default = .false.
        if (Node.eq.0) then
          inspn = fdf_boolean('DM.InitSpinAF',inspn_default)
          write(6,'(a,4x,l1)')
     .     'redata: Antiferro initial spin density   = ',inspn
        endif
#ifdef MPI
        call MPI_Bcast(inspn,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
      endif
C ...

C Use continumation DM files
*     usdm_default = .true.
      if (Node.eq.0) then
        usdm_default = fdf_boolean('UseSaveData',.false.)
        usesavedm  = fdf_boolean('DM.UseSaveDM',usdm_default)
        write(6,'(a,4x,l1)') 
     .   'redata: Use continuation files for DM    = ',
     .   usesavedm
      endif
#ifdef MPI
      call MPI_Bcast(usesavedm,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
C ...

C Neglect Interactions between non-overlapping orbitals ...
      negl_default = .false.
      if (Node.eq.0) then
        negl  = fdf_boolean('NeglNonOverlapInt',negl_default)
        write(6,'(a,4x,l1)') 
     .   'redata: Neglect nonoverlap interactions  = ',negl
      endif
#ifdef MPI
      call MPI_Bcast(negl,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
C ...

C Method to Solve LDA Hamitonian ...
      if (Node.eq.0) then
        if (na .le. na_diag) then
          method = fdf_string('SolutionMethod','diagon')
        else
          method = fdf_string('SolutionMethod','ordern')
        endif
      endif
#ifdef MPI
      call MPI_Bcast(method,6,MPI_character,0,MPI_Comm_World,MPIerror)
#endif
      if (leqi(method,'diagon')) then
        isolve = 0
        if (Node.eq.0) then
          write(6,'(a,a)') 
     .     'redata: Method of Calculation            = ',
     .     'Diagonalization'
        endif
      else if (leqi(method,'ordern')) then
        isolve = 1
        if (Node.eq.0) then
          write(6,'(a,a)') 
     .     'redata: Method of Calculation            = ',
     .     'Order-N'
        endif
        if (nspin .gt. 1) then
          if (Node.eq.0) then
           write(6,100) 
           write(6,101) 
           write(6,'(a)') 
     .     'redata:    You chose the Order-N solution option'
           write(6,'(a)') 
     .     'redata:    together with nspin>1.  This is not  '
           write(6,'(a)') 
     .     'redata:    allowed in this version of siesta    '
           write(6,102)
         endif
         stop
        endif
      else
        if (Node.eq.0) then
          write(6,100) 
          write(6,101) 
          write(6,'(a)') 
     .    'redata:    The method of solution must be either'
          write(6,'(a)') 
     .    'redata:           OrderN or Diagon'
          write(6,102)
        endif
        stop
      endif

C Memory scaling factor for rdiag/cdiag - cannot be less than 1.0
      if (Node.eq.0) then
        MemoryFactor = fdf_double('DiagMemory', 2.0d0 )
        MemoryFactor = max(MemoryFactor,1.0d0)
      endif
#ifdef MPI
      call MPI_Bcast(MemoryFactor,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

C Memory increasing allowed if diagonalisation fails to converge?
      if (Node.eq.0) then
        TryMemoryIncrease = fdf_boolean('TryMemoryIncrease', .true. )
      endif
#ifdef MPI
      call MPI_Bcast(TryMemoryIncrease,1,MPI_logical,0,MPI_Comm_World,
     .  MPIerror)
#endif
C ...

C Electronic temperature for Fermi Smearing ...
      if (Node.eq.0) then
        temp = fdf_physical('ElectronicTemperature',temp_default,'Ry')
        if (isolve .eq. 0) then
          write(6,'(a,f10.4,a)') 
     .    'redata: Electronic Temperature           = ',temp,'  Ry'
        endif
      endif
#ifdef MPI
      call MPI_Bcast(temp,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
C ...

C Fix the spin of the system to a given value 
      fs_default = .false.
      if (Node.eq.0) then
        fixspin  = fdf_boolean('FixSpin',fs_default)
        if (fixspin .and. nspin .ne. 2) then
          write(6,'(a)') 
     .   'redata: ERROR: You can only fix the spin of the system' 
          write(6,'(a)') 
     .   'redata:        for collinear spin polarized calculations.'
          write(6,102)
        else
          write(6,'(a,4x,l1)') 
     .     'redata: Fix the spin of the system       = ',fixspin 
        endif
      endif
#ifdef MPI
      call MPI_Bcast(fixspin,1,MPI_logical,0,MPI_Comm_World,
     .  MPIerror)
#endif
      if (fixspin .and. nspin .ne. 2) then
        stop
      endif
C ...

C Value of the Spin of the system (only used if fixspin = TRUE
      if (fixspin) then
        if (Node.eq.0) then
          ts = fdf_double('TotalSpin',0.0d0)
          write(6,'(a,f10.4)') 
     .     'redata: Value of the Spin of the System  = ',ts
        endif
#ifdef MPI
        call MPI_Bcast(ts,1,DAT_double,0,MPI_Comm_World,
     .    MPIerror)
#endif
      else
        ts = 0.0d0
      endif
C ...
        

C Order-N solution parameters ...
C     Maximum number of CG minimization iterations
      ncgmax_default = 1000
      if (Node.eq.0) then
        ncgmax = fdf_integer('ON.MaxNumIter',ncgmax_default)
      endif
#ifdef MPI
      call MPI_Bcast(ncgmax,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

C     Relative tolerance in total band structure energy
      etol_default = 1.d-8
      if (Node.eq.0) then
        etol = fdf_double('ON.etol',etol_default)
      endif
#ifdef MPI
      call MPI_Bcast(etol,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

C     Fermi level parameter
      eta_default = 0.d0
      if (Node.eq.0) then
        eta = fdf_physical('ON.eta',eta_default,'Ry')
      endif
#ifdef MPI
      call MPI_Bcast(eta,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

C     Cutoff radius for Localized Wave Functions
      rcoor_default = 9.5d0
      if (Node.eq.0) then
        rcoor = fdf_physical('On.RcLWF',rcoor_default,'Bohr')
      endif
#ifdef MPI
      call MPI_Bcast(rcoor,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

C     Use continumation LWF files
*     uslwf_default = .true.
      if (Node.eq.0) then
        uslwf_default = fdf_boolean('UseSaveData',.false.)
        usesavelwf  = fdf_boolean('ON.UseSaveLWF',uslwf_default)
      endif
#ifdef MPI
      call MPI_Bcast(usesavelwf,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

C     Option on how to build LWF's (disk or functionals)
      lwfopt_default = 'kim'
      if (Node.eq.0) then
        lwfopt  = fdf_string('ON.functional',lwfopt_default)
      endif
#ifdef MPI
      call MPI_Bcast(lwfopt,13,MPI_character,0,MPI_Comm_World,MPIerror)
#endif
      if (leqi(lwfopt,'files')) then
        ioptlwf = 0
      else if (leqi(lwfopt,'kim')) then
        ioptlwf = 1
      else if (leqi(lwfopt,'ordejon-mauri')) then
        ioptlwf = 2
      else
        if (Node.eq.0) then
          write(6,'(a)') 'redata: wrong ON.funcional option'
        endif
        stop
      endif

C     Option to calculate the Chemical potential in O(N)
      chebef_default = .false.
      noeta_default = .false.
      if (Node.eq.0) then
        chebef = fdf_boolean('ON.ChemicalPotential',chebef_default)
      
C     Option to use the Chemical Potential calculated instead
C     of the eta variable of the input
        noeta = fdf_boolean('ON.ChemicalPotentialUse',noeta_default)
      endif
#ifdef MPI
      call MPI_Bcast(chebef,1,MPI_logical,0,MPI_Comm_World,
     .  MPIerror)
      call MPI_Bcast(noeta,1,MPI_logical,0,MPI_Comm_World,
     .  MPIerror)
#endif

      if (noeta) chebef = .true.

C     Cutoff radius to calculate the Chemical Potential by projection
      rccp_default = 9.5d0
      if (Node.eq.0) then
        rcoorcp = fdf_physical('ON.ChemicalPotentialRc',
     .                        rccp_default,'Bohr')
      endif
#ifdef MPI
      call MPI_Bcast(rcoorcp,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

C     Temperature of the Fermi distribution to calculate the
C     Chemical potential by projection
      tcp_default = 0.05
      if (Node.eq.0) then
        tcp = fdf_physical('ON.ChemicalPotentialTemperature',
     .                   tcp_default,'Ry')
      endif
#ifdef MPI
      call MPI_Bcast(tcp,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
      beta = 1.0d0/tcp

C     Order of the Chebishev expansion to calculate the Chemical
C     potential
      pmax_default = 100
      if (Node.eq.0) then
        pmax = fdf_integer('ON.ChemicalPotentialOrder',pmax_default)
      endif
#ifdef MPI
      call MPI_Bcast(pmax,1,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
#endif
C ...


      if (isolve .eq. 1.and.Node.eq.0) then
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
      if (Node.eq.0) then
        varcel = fdf_boolean('MD.VariableCell', .false. )
      endif
#ifdef MPI
      call MPI_Bcast(varcel,1,MPI_logical,0,MPI_Comm_World,
     .  MPIerror)
#endif
C     Type of dynamics 
      dyntyp_default = 'verlet'   
      if (Node.eq.0) then
        dyntyp = fdf_string('MD.TypeOfRun',dyntyp_default)
      endif
#ifdef MPI
      call MPI_Bcast(dyntyp,22,MPI_character,0,MPI_Comm_World,
     .  MPIerror)
#endif
      if (leqi(dyntyp,'cg')) then
        idyn = 0
        if (Node.eq.0) then
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    CG coord. optimization'
          write(6,'(a,4x,l1)')
     .     'redata: Variable cell                    = ', varcel
          uscg_default = fdf_boolean('UseSaveData',.false.)
          usesavecg  = fdf_boolean('MD.UseSaveCG',uscg_default)
          write(6,'(a,4x,l1)')
     .     'redata: Use continuation files for CG    = ',
     .     usesavecg
        endif
#ifdef MPI
        call MPI_Bcast(usesavecg,1,MPI_logical,0,MPI_Comm_World,
     .    MPIerror)
#endif
      else if (leqi(dyntyp,'verlet')) then
        idyn = 1
        if (Node.eq.0) then
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    Verlet MD run'
        endif
      else if (leqi(dyntyp,'nose')) then
        idyn = 2
        if (Node.eq.0) then
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    Nose termostat MD run'
        endif
      else if (leqi(dyntyp,'parrinellorahman')) then
        idyn = 3
        if (Node.eq.0) then
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    Parrinello-Rahman MD run'
        endif
      else if (leqi(dyntyp,'noseparrinellorahman')) then
        idyn = 4
        if (Node.eq.0) then
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    Nose-Parrinello-Rahman MD run'
        endif
      else if (leqi(dyntyp,'anneal')) then
        idyn = 5
        if (Node.eq.0) then
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    Annealing MD run'
        endif
      else if (leqi(dyntyp,'fc')) then
        idyn = 6
        if (Node.eq.0) then
          write(6,'(a,a)')
     .     'redata: Dynamics option                  = ',
     .     '    Force Constants Matrix Calculation'
        endif
      else
        if (Node.eq.0) then
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
        endif
        stop
      endif 


C     Maximum number of steps in CG coordinate optimization
      nmove_default = 0
      if (Node.eq.0) then
        nmove = fdf_integer('MD.NumCGsteps',nmove_default)
      endif
#ifdef MPI
      call MPI_Bcast(nmove,1,MPI_integer,0,MPI_Comm_World,
     .    MPIerror)
#endif

C     Maximum atomic displacement in one CG step
      dxmax_default = 0.2d0
      if (Node.eq.0) then
        dxmax = fdf_physical('MD.MaxCGDispl',dxmax_default,'Bohr')
      endif
#ifdef MPI
      call MPI_Bcast(dxmax,1,DAT_double,0,MPI_Comm_World,
     .    MPIerror)
#endif

C     Tolerance in the maximum atomic force (def 0.04 eV/Ang)
      ftol_default = 0.00155574d0 
      if (Node.eq.0) then
        ftol = fdf_physical('MD.MaxForceTol',ftol_default,'Ry/Bohr')
      endif
#ifdef MPI
      call MPI_Bcast(ftol,1,DAT_double,0,MPI_Comm_World,
     .    MPIerror)
#endif

C     Tolerance in the maximum residual stress (var cell) def = 1 GPa 
      if (Node.eq.0) then
        strtol = fdf_physical('MD.MaxStressTol',6.79773d-5,'Ry/Bohr**3')
        strtol = dabs(strtol)
      endif
#ifdef MPI
      call MPI_Bcast(strtol,1,DAT_double,0,MPI_Comm_World,
     .    MPIerror)
#endif

      if (idyn .eq. 0.and.Node.eq.0) then
        write(6,'(a,i5)') 
     .  'redata: Maximum number of CG moves       = ',nmove
        write(6,'(a,f10.4,a)') 
     .  'redata: Max atomic displ per move        = ',dxmax,'  Bohr'
        write(6,'(a,f10.4,a)') 
     .  'redata: Force tolerance                  = ',ftol,'  Ry/Bohr'
        if ( varcel ) then
           write(6,'(a,f10.4,a)')
     .  'redata: Stress tolerance                 = ',strtol/6.79773d-5,
     .                                              '  GPa'
        endif

      endif
  
C     Initial time step for MD
      istart_default = 1
      if (Node.eq.0) then
        istart = fdf_integer('MD.InitialTimeStep',istart_default)
      endif
#ifdef MPI
      call MPI_Bcast(istart,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

C     Final time step for MD
      ifinal_default = 1
      if (Node.eq.0) then
        ifinal = fdf_integer('MD.FinalTimeStep',ifinal_default)
      endif
#ifdef MPI
      call MPI_Bcast(ifinal,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

C     Length of time step for MD
      dt_default = 1.d0
      if (Node.eq.0) then
        dt = fdf_physical('MD.LengthTimeStep',dt_default,'fs')
      endif
#ifdef MPI
      call MPI_Bcast(dt,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

C     Quench Option
      qnch_default = .false.
      if (Node.eq.0) then
        qnch = fdf_boolean('MD.Quench',qnch_default)
      endif
#ifdef MPI
      call MPI_Bcast(qnch,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

      if (qnch .and. (idyn .eq. 2 .or. idyn .eq. 4)) then
        if (Node.eq.0) then
          write(6,100) 
          write(6,101) 
          write(6,'(a)') 
     .    'redata: ERROR: You cannot quench and use a Nose'
          write(6,'(a)') 
     .    'redata: ERROR: thermostat simultaneously.'
          write(6,102)
        endif
        stop
      endif
      iquench = 0
      if (qnch) iquench = 1

C     Initial Temperature of MD simulation
C     (draws random velocities from the Maxwell-Boltzmann distribition
C      at the given temperature)
      ti_default = 0.d0
      if (Node.eq.0) then
        tempinit = fdf_physical('MD.InitialTemperature',ti_default,'K')
      endif
#ifdef MPI
      call MPI_Bcast(tempinit,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

      if (Node.eq.0) then
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
      endif

C     Target Temperature
      tt_default = 0.d0
      if (Node.eq.0) then
        tt = fdf_physical('MD.TargetTemperature',tt_default,'K')
      endif
#ifdef MPI
      call MPI_Bcast(tt,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
      
C     Target Pressure
      tp_default = 0.d0
      if (Node.eq.0) then
        tp = fdf_physical('MD.TargetPressure',tp_default,'Ry/Bohr**3')
      endif
#ifdef MPI
      call MPI_Bcast(tp,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

C     Mass of Nose variable
      mn_default = 1.d2
      if (Node.eq.0) then
        mn = fdf_physical('MD.NoseMass',mn_default,'Ry*fs**2')
      endif
#ifdef MPI
      call MPI_Bcast(mn,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

C     Mass of Parrinello-Rahman variables
      mpr_default = 1.d2
      if (Node.eq.0) then
        mpr = fdf_physical('MD.ParrinelloRahmanMass',
     .                    mpr_default,'Ry*fs**2')
      endif
#ifdef MPI
      call MPI_Bcast(mpr,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

      if (Node.eq.0) then
      if (idyn .eq. 2 .or. idyn .eq. 4) then
        write(6,'(a,f10.4,a)') 
     .  'redata: Nose mass                        = ',mn,'  Ry/fs**2'
      endif

      if (idyn .eq. 3 .or. idyn .eq. 4) then
        write(6,'(a,f10.4,a)') 
     .  'redata: Parrinello-Rahman mass           = ',mpr,'  Ry/fs**2'
      endif
      endif

C     Annealing option
      annop_default = 'TemperatureAndPressure'
      if (Node.eq.0) then
        annop = fdf_string('MD.AnnealOption',annop_default)
      endif
#ifdef MPI
      call MPI_Bcast(annop,22,MPI_character,0,MPI_Comm_World,MPIerror)
#endif
      ianneal = 0
      if (idyn .eq. 5) then
        if (leqi(annop,'Temperature')) then
          ianneal = 1
          if (Node.eq.0) then
            write(6,'(a,a)') 
     .       'redata: Annealing Option                 = ',
     .       'Temperature'
          endif
        else if (leqi(annop,'Pressure')) then
          ianneal = 2
          if (Node.eq.0) then
            write(6,'(a,a)') 
     .       'redata: Annealing Option                 = ',
     .       'Pressure'
          endif
        else if (leqi(annop,'TemperatureAndPressure')) then
          ianneal = 3
          if (Node.eq.0) then
            write(6,'(a,a)') 
     .       'redata: Annealing Option                 = ',
     .       'Temperature and Pressure'
          endif
        else
          if (Node.eq.0) then
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
          endif
          stop
        endif
      endif

      if (Node.eq.0) then
        if (idyn .eq. 2 .or. idyn .eq. 4 .or. 
     .   (idyn .eq. 5 .and. (ianneal .eq. 1 .or. ianneal .eq. 3))) then
          write(6,'(a,f10.4,a)') 
     .    'redata: Target Temperature               = ',tt,'  Kelvin'
      endif

      if (idyn .eq. 3 .or. idyn .eq. 4 .or. 
     .   (idyn .eq. 5 .and. (ianneal .eq. 2 .or. ianneal .eq. 3))) then
        write(6,'(a,f10.4,a)') 
     .  'redata: Target Pressure                  = ',tt,'  Ry/Bohr**3'
      endif
      endif

C     Relaxation Time for Annealing
      taurelax_default = 1.d2
      if (Node.eq.0) then
      taurelax = fdf_physical('MD.TauRelax',taurelax_default,'fs')
      if (idyn .eq. 5) then
        write(6,'(a,f10.4,a)') 
     .  'redata: Annealing Relaxation Time        = ',
     .   taurelax,'  fs'
      endif
      endif
#ifdef MPI
      call MPI_Bcast(taurelax,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
        
C     Estimated Bulk modulus (for Pressure annealing)
      bulkm_default = 1.d2
      if (Node.eq.0) then
        bulkm = fdf_double('MD.BulkModulus',bulkm_default)
        if (idyn .eq. 5 .and. (ianneal.eq.2.or.ianneal.eq.3)) then
          write(6,'(a,f10.4,a)') 
     .    'redata: Approx. Bulk Modulus             = ',
     .     bulkm,'  Ry/Bohr**3'
        endif
      endif
#ifdef MPI
      call MPI_Bcast(bulkm,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

C     Atomic displacement for force constant calculation
      dx_default = 0.04d0
      if (Node.eq.0) then
        dx = fdf_physical('MD.FCDispl',dx_default,'Bohr')
      endif
#ifdef MPI
      call MPI_Bcast(dx,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

C     First and last atoms to displace for calculation of force constants
      ia1_default = 1
      ia2_default = na
      if (Node.eq.0) then
        ia1 = fdf_integer('MD.FCfirst',ia1_default)
        ia2 = fdf_integer('MD.FClast',ia2_default)
      endif
#ifdef MPI
      call MPI_Bcast(ia1,1,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
      call MPI_Bcast(ia2,1,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
#endif

      if (idyn .eq. 6.and.Node.eq.0) then
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
      if (Node.eq.0) then
        call recoor( na*ncells, na, isa, xa )
      endif

#ifdef MPI
      call MPI_Bcast(na,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(isa,na,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(xa(1,1),3*na,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
      
      if (Node.eq.0) then
        write(6,102)
      endif

100   format(/,'redata: ',71(1h*))
101   format('redata:                  INPUT ERROR')
102   format('redata: ',71(1h*))
103   format('redata: ',i4,2x,3f10.5,i3) 

      return
      end
