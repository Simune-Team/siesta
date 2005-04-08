      subroutine redata(na, ns, nspin, outlng, g2max,
     .                  charnet, negl, nscf, dDtol, mix, wmix, isolve, 
     .                  temp, fixspin, ts, ncgmax, ftol, strtol, eta, 
     .                  etol, rcoor, 
     .                  ioptlwf, chebef, noeta, rcoorcp, beta, pmax,
     .                  idyn, istart, ifinal, nmove, ianneal, iquench,
     .                  dt, ia1, ia2, dx, dxmax, tt, tp, mn, mpr, 
     .                  bulkm, taurelax, 
     .                  usesavelwf, usesavedm, usesavecg,
     .                  mullipop, inspn, maxsav, nkick, wmixkick, 
     .                  pulfile, tempinit, dumpcharge, varcel, 
     .                  harrisfun, occtol )
C *********************************************************************
C Subroutine to read the data for the SIESTA program
C
C     It uses the FDF (Flexible Data Format) package 
C     of J.M.Soler and A.Garcia
C
C Writen by P.Ordejon, December'96
C Modified for introduction of dynamic memory in SIESTA by JDG Sept 99
C ***************************** INPUT *********************************
C integer na               : Number of atoms
C integer nspin            : Spin polarization
C **************************** OUTPUT *********************************
C integer ns               : Number of species
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
C real*8 eta(2)            : Fermi level parameter of Kim functional
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
C                             6 = Force constants
C                             7 = Forces for PHONON program
C                             8 = Force evaluation
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
C integer mullipop         : Option for Mulliken Pop. analysis
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
C     (This is now obsolete: info will appear in .RHO file)
C logical varcel           : variable shape for CG optimization or dynamics
C logical harrisfun        : swith that indicates if harris functional will
C                            be used or not
C real*8  occtol           : Occupancy threshold for DM build
C **********************************************************************

C
C  Modules
C
      use precision, only : dp
      use parallel,  only : IOnode,node,  Nodes
      use fdf
      use sys
#ifdef MPI
      use mpi_siesta, only : MPI_Comm_World, MPI_logical, MPI_integer,
     .                       MPI_character
#ifdef NODAT
      use mpi_siesta, only: mpi_double_precision
#else
      use mpi_siesta, only: DAT_double
#endif
#endif
      use diagmemory, only: memoryfactor

      use siesta_cml

      implicit none

! Subroutine arguments:
      integer, intent(in)  :: na, nspin

      integer, intent(out) :: 
     .  ia1, ia2, ianneal, idyn, ifinal, ioptlwf,
     .  iquench, isolve, istart, maxsav,
     .  mullipop, ncgmax, nkick, nmove, ns, nscf,
     .  pmax
      real(dp), intent(out) ::
     .  beta, bulkm, charnet,
     .  dDtol, dt, dx, dxmax, eta(2), etol, ftol, g2max,
     .  mn, mpr, occtol, rcoor, rcoorcp, strtol,
     .  taurelax, temp, tempinit, tp, ts, tt, wmix, wmixkick
      logical, intent(out) ::
     .  chebef, dumpcharge, fixspin, inspn, mix, negl, noeta, 
     .  outlng, pulfile, usesavecg, usesavelwf, usesavedm, varcel,
     .  harrisfun

C Internal parameters ................................................
C na_diag      : maximum number of atoms with diagon as default method
C g2max_default : Mesh cutoff default, in Ry
C temp_default  : Electronic temperature default, in Ry
      integer na_diag
      real(dp) g2max_default, temp_default
      parameter ( na_diag       = 100         )
      parameter ( g2max_default = 100.e0_dp   )
      parameter ( temp_default  = 1.900e-3_dp )
C ................

C  Internal variables .................................................
      real(dp) :: tcp

      character
     .  annop*22, dyntyp*22, 
     .  method*6, lwfopt*13

#ifdef MPI
      integer
     .  MPIerror
#endif

      logical
     .  DaC, leqi, qnch

! Template to read parameter:
!
!      param = default_value
!      if (ionode) then
!        param = fdf_integer('ParamName', param)
!        write(6,'(a,i)'), 'redata: ParamName           = ',param
!      endif
!      if (cml_p) call cmlAddParameter(xf=mainXML, name='ParamName',
!     .                 value=param, dictref='siesta:param')
!#ifdef MPI
!      call MPI_Bcast(param, 1, MPI_integer, 0, MPI_Comm_World, MPIerror)
!#endif
!
! Adjust as appropriate for logical/real parameters. May be more complicated
! where dependencies exist.

C Start of code
      if (cml_p) then
        call cmlStartParameterList(mainXML, title='Input Parameters')
      endif

C Type of output
      if (ionode) then
        outlng = fdf_boolean( 'LongOutput', .false. )
        write(6,'(a,4x,l1)')
     .   'redata: Long output                      = ',outlng
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, name='LongOutput',
     .       value=outlng, dictRef='siesta:verbosity')
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(outlng,1,MPI_logical,0,MPI_Comm_World,
     .  MPIerror)
#endif

C Define Number of species ...
      ns = 0
      if (ionode) then
        ns = fdf_integer('NumberOfSpecies', ns)
        write(6,'(a,3x,i2)') 
     .   'redata: Number of Atomic Species         = ',ns
        if (ns .le. 0) call die('redata: ERROR:'// 
     .       'Number of species must be larger than zero.')
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='NumberOfSpecies', value=ns, dictRef='siesta:ns')
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(ns,1,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
#endif

C ...

C Dump information to plot charge contours
C by the external DENCHAR application program.
      dumpcharge = .false.
      if (ionode) then
        dumpcharge  = fdf_boolean('WriteDenchar', dumpcharge)
        write(6,'(a)') 
     .   'redata: Charge density info will appear in .RHO file'
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='WriteDenChar', value=dumpcharge)
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(dumpcharge,1,MPI_logical,0,MPI_Comm_World,
     .  MPIerror)
#endif
C ...

C Perform Mulliken Population Analysis
      mullipop = 0
      if (outlng) mullipop = 1
      if (ionode) then
        mullipop = fdf_integer('WriteMullikenPop', mullipop)
        select case (mullipop)
          case(0)
            write(6,'(a)') 
     .           'redata: Write Mulliken Pop.              =     NO'
          case(1)
            write(6,'(a,a)') 
     .           'redata: Write Mulliken Pop.              =     ',
     .           'Atomic and Orbital charges'
          case(2)
            write(6,'(a,a/45x,a)') 
     .           'redata: Write Mulliken Pop.              =     ',
     .           'Atomic and Orbital charges',
     .           'plus Atomic Overlap Populations'
          case(3)
            write(6,'(a,a/45x,a/45x,a)') 
     .           'redata: Write Mulliken Pop.              =     ',
     .           'Atomic and Orbital charges',
     .           'plus Atomic Overlap Populations',
     .           'plus Orbital Overlap Populations'
          case default
            call die('redata: Invalid value for WriteMullikenPop')
        end select
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='WriteMullikenPop', value=mullipop)
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(mullipop,1,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
#endif
C ...

C Planewave cutoff of the real space mesh ...
      if (ionode) then
        g2max = g2max_default
        g2max = fdf_physical('MeshCutoff',g2max,'Ry')
        write(6,'(a,f10.4,a)') 
     .   'redata: Mesh Cutoff                      = ',g2max,'  Ry'
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, name='MeshCutOff',
     .     value=g2max, dictRef='siesta:g2max', units='Ry')
C Global broadcast of values
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(g2max,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(g2max,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif
C ...

C Net charge in the cell ...
      charnet = 0.0_dp
      if (ionode) then
        charnet = fdf_double('NetCharge',charnet)
        write(6,'(a,f10.4,a)') 
     .   'redata: Net charge of the system         = ',charnet,' |e|'
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, name='NetCharge',
     .   value=charnet, dictRef='siesta:NetCharge', units='|e|')
C Global broadcast of values
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(charnet,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(charnet,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif
C ...
       
C SCF Loop parameters ...
C     Maximum number of SCF iterations
      nscf = 50
      if (ionode) then
        nscf = fdf_integer('MaxSCFIterations',nscf)
        write(6,'(a,i5)') 
     .   'redata: Max. number of SCF Iter          = ',nscf
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='MaxSCFIterations', value=nscf, dictRef='siesta:maxscf')
#ifdef MPI
      call MPI_Bcast(nscf,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

C Pulay mixing, numer of iterations for one Pulay mixing (maxsav)
      maxsav = 0
      if (ionode) then
        maxsav = fdf_integer('DM.NumberPulay', maxsav)
        if (maxsav .gt. 1) then
          write(6,'(a,i5,a)') 
     .     'redata: One Pulay mixing every           = ',maxsav,
     .     ' iterations'
        else
           write(6,'(a)')'redata: Mixing is linear'
        endif
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, name='DM.NumberPulay',
     .     value=maxsav, dictRef='siesta:maxsav')
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(maxsav,1,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
#endif
C ...

C     Mix density matrix on first SCF step
C     (mix)
      mix = .false.
      if (ionode) then
        mix  = fdf_boolean('DM.MixSCF1',mix)
        write(6,'(a,4x,l1)')
     .    'redata: Mix DM in first SCF step ?       = ',mix
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, name='DM.MixSCF1',
     .     value=mix, dictRef='siesta:mix')
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(mix,1,MPI_logical,0,MPI_Comm_World,
     .  MPIerror)
#endif

C Use disk or memory to store intermediate Pulay miximg vectors
C (pulfile)
      pulfile = .false.
      if (ionode) then
        pulfile  = fdf_boolean('DM.PulayOnFile',pulfile)
        if (pulfile.and.Nodes>1) then
          call die('redata: Cannot use DM.PulayOnFile=.true.'//
     .         'when running in parallel')
        endif
        write(6,'(a,4x,l1)')
     .    'redata: Write Pulay info on disk?        = ',pulfile
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, name='DM.PulayOnFile',
     .     value=pulfile, dictRef='siesta:pulfile')
#ifdef MPI
      call MPI_Bcast(pulfile,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
C ...

C Density Matrix Mixing  (proportion of output DM in new input DM)
      wmix = 0.25_dp
      if (ionode) then
        wmix = fdf_double('DM.MixingWeight',wmix)
        write(6,'(a,f10.4,a)') 
     .   'redata: New DM Mixing Weight             = ',wmix
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='DM.MixingWeight', value=wmix, dictRef='siesta:wmix')
C Global broadcast of values
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(wmix,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(wmix,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif

C Density Matrix occupancy tolerance
      occtol = 1.0e-12_dp
      if (ionode) then
        occtol = fdf_double('DM.OccupancyTolerance',occtol)
        write(6,'(a,f14.12)') 
     .   'redata: New DM Occupancy tolerance       = ',occtol
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='DM.OccupancyTolerance', value=occtol, 
     .     dictRef='siesta:occtol')
C Global broadcast of values
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(occtol,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(occtol,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif

C Perform linear mixing each nkick SCF iterations (to kick system
C when it is pinned in a poorly convergent SCF loop)
      nkick = 0
      if (ionode) then
        nkick = fdf_integer('DM.NumberKick',nkick)
        if (nkick .ge. 1) then
          write(6,'(a,i5,a)')
     .     'redata: Kick with linear mixing every    = ',nkick,
     .     ' iterations'
        else
           write(6,'(a)')'redata: No kicks to SCF'
        endif
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='DM.NumberKick', value=nkick, dictRef='siesta:nkick')
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(nkick,1,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
#endif
 
C     Density Matrix Mixing each nkick SCF iterations
      wmixkick = 0.50_dp
      if (ionode) then
        wmixkick = fdf_double('DM.KickMixingWeight',wmixkick)
        write(6,'(a,f10.4,a)')
     .   'redata: DM Mixing Weight for Kicks       = ',wmixkick
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='DM.KickMixingWeight', value=wmixkick, 
     .     dictRef='siesta:wmixkick')
C Global broadcast of values
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(wmixkick,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(wmixkick,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif

C     Density Matrix Tolerance for achieving Self-Consistency
      dDtol = 1.e-4_dp
      if (ionode) then
        dDtol = fdf_double('DM.Tolerance',dDtol)
        write(6,'(a,f12.6,a)') 
     .   'redata: DM Tolerance for SCF             = ',dDtol
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, name='DM.Tolerance', 
     .     value=dDtol, dictRef='siesta:dDtol')
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(dDtol,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(dDtol,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif

C     Initial spin density: Maximum polarization, Ferro (false), AF (true)
      if (nspin.eq.2) then
        inspn = .false.
        if (ionode) then
          inspn = fdf_boolean('DM.InitSpinAF',inspn)
          write(6,'(a,4x,l1)')
     .     'redata: Antiferro initial spin density   = ',inspn
        endif
        if (cml_p) call cmlAddParameter(xf=mainXML, 
     .       name='DM.InitSpinAF', value=inspn, dictRef='siesta:inspn')
#ifdef MPI
        call MPI_Bcast(inspn,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
      endif
C ...

C Use continuation DM files
      usesavedm = .false.
      if (ionode) then
        usesavedm = fdf_boolean('UseSaveData',usesavedm)
        usesavedm  = fdf_boolean('DM.UseSaveDM',usesavedm)
        write(6,'(a,4x,l1)') 
     .   'redata: Use continuation files for DM    = ',
     .   usesavedm
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, name='DM.UseSaveDM', 
     .     value=usesavedm, dictRef='siesta:usesavedm')
#ifdef MPI
      call MPI_Bcast(usesavedm,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
C ...

C Neglect Interactions between non-overlapping orbitals ...
      negl = .false.
      if (ionode) then
        negl  = fdf_boolean('NeglNonOverlapInt',negl)
        write(6,'(a,4x,l1)') 
     .   'redata: Neglect nonoverlap interactions  = ',negl
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='NeglNonOverlapInt', value=negl, dictRef='siesta:negl')
#ifdef MPI
      call MPI_Bcast(negl,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
C ...

C Method to Solve LDA Hamiltonian ...
      if (ionode) then
        if (na .le. na_diag) then
          method = fdf_string('SolutionMethod','diagon')
        else
          method = fdf_string('SolutionMethod','ordern')
        endif
        if (leqi(method,'diagon')) then
          isolve = 0
          write(6,'(a,4x,a)') 
     .         'redata: Method of Calculation            = ',
     .         'Diagonalization'
          DaC=.false.
          if (Nodes.eq.1) then
            DaC  = fdf_boolean('Diag.DivideAndConquer', DaC)
          endif
          write(6,'(a,4x,l1)') 
     .         'redata: Divide and Conquer               = ',DaC
        else if (leqi(method,'ordern')) then
          isolve = 1
          if (ionode) then
            write(6,'(a,4x,a)') 
     .           'redata: Method of Calculation            = ',
     .           'Order-N'
          endif
          if (nspin .gt. 2) then
            call die('redata: You chose the Order-N solution option '//
     .           'together with nspin>2.  This is not allowed in '//
     .           'this version of siesta')
          endif
        else
          call die('redata: The method of solution must be either '//
     .         'OrderN or Diagon')
        endif
      endif
      if (cml_p) then
        call cmlAddParameter(xf=mainXML, name='SolutionMethod',
     .       value=method, dictRef='siesta:SCFmethod')
        call cmlAddParameter(xf=mainXML, name='Diag.DivideAndConquer',
     .       value=DaC, dictRef='siesta:DaC')
      endif
#ifdef MPI
      call MPI_Bcast(method,6,MPI_character,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(isolve,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(DaC,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
C ...

C Memory scaling factor for rdiag/cdiag - cannot be less than 1.0
      MemoryFactor = 1.0_dp
      if (ionode) then
        MemoryFactor = fdf_double('Diag.Memory', MemoryFactor )
        MemoryFactor = max(MemoryFactor,1.0_dp)
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='Diag.Memory', value=MemoryFactor, 
     .     dictRef='siesta:MemoryFactor')
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(MemoryFactor,1,MPI_double_precision,0,
     .  MPI_Comm_World,MPIerror)
#else
      call MPI_Bcast(MemoryFactor,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif
C ...

C Electronic temperature for Fermi Smearing ...
      temp = temp_default
      if (ionode) then
        temp = fdf_physical('ElectronicTemperature',temp,'Ry')
        if (isolve .eq. 0) then
          write(6,'(a,f10.4,a)') 
     .    'redata: Electronic Temperature           = ',temp,'  Ry'
        endif
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='ElectronicTemperature', value=temp, 
     .     dictRef='siesta:etemp', units='Ry')
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(temp,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(temp,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif
C ...

C Fix the spin of the system to a given value ; and
C value of the Spin of the system (only used if fixspin = TRUE
      fixspin = .false.
      ts = 0.0_dp
      if (ionode) then
        fixspin  = fdf_boolean('FixSpin',fixspin)
        if (fixspin) then
          if (nspin .ne. 2) then
            call die('redata: ERROR: You can only fix the spin of '//
     .          'the system for collinear spin polarized calculations.')
          endif
          ts = fdf_double('TotalSpin',ts)
          write(6,'(a,4x,l1)') 
     .         'redata: Fix the spin of the system       = ',fixspin 
          write(6,'(a,f10.4)') 
     .     'redata: Value of the Spin of the System  = ',ts
        endif
      endif
      if (cml_p) then 
        call cmlAddParameter(xf=mainXML, name='FixSpin',
     .       value=fixspin, dictref='siesta:fixspin')
        call cmlAddParameter(xf=mainXML, name='TotalSpin',
     .       value=ts, dictref='siesta:ts')
      endif
#ifdef MPI
      call MPI_Bcast(fixspin,1,MPI_logical,0,MPI_Comm_World,
     .  MPIerror)
#ifdef NODAT
        call MPI_Bcast(ts,1,MPI_double_precision,0,MPI_Comm_World,
     .    MPIerror)
#else
        call MPI_Bcast(ts,1,DAT_double,0,MPI_Comm_World,
     .    MPIerror)
#endif
#endif
C ...

C Order-N solution parameters ...
C     Maximum number of CG minimization iterations
      ncgmax = 1000
      if (ionode) then
        ncgmax = fdf_integer('ON.MaxNumIter',ncgmax)
      endif
#ifdef MPI
      call MPI_Bcast(ncgmax,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
C ...

C     Relative tolerance in total band structure energy
      etol = 1.e-8_dp
      if (ionode) then
        etol = fdf_double('ON.etol',etol)
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(etol,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(etol,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif
C ...

C     Fermi level parameter
      eta = 0.0_dp
      if (ionode) then
        eta(1) = fdf_physical('ON.eta',eta(1),'Ry')
        eta(2) = eta(1)
        eta(1) = fdf_physical('ON.eta_alpha',eta(1),'Ry')
        eta(2) = fdf_physical('ON.eta_beta',eta(2),'Ry')
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(eta,2,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(eta,2,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif
C ...

C     Cutoff radius for Localized Wave Functions
      rcoor = 9.5_dp
      if (ionode) then
        rcoor = fdf_physical('On.RcLWF',rcoor,'Bohr')
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(rcoor,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(rcoor,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif
C ...

C     Use continumation LWF files
      usesavelwf = .false.
      if (ionode) then
        usesavelwf = fdf_boolean('UseSaveData',usesavelwf)
        usesavelwf  = fdf_boolean('ON.UseSaveLWF',usesavelwf)
      endif
#ifdef MPI
      call MPI_Bcast(usesavelwf,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
C ...

C     Option on how to build LWF's (disk or functionals)
      lwfopt = 'kim'
      if (ionode) then
        lwfopt  = fdf_string('ON.functional',lwfopt)
        if (leqi(lwfopt,'files')) then
          ioptlwf = 0
        else if (leqi(lwfopt,'kim')) then
          ioptlwf = 1
        else if (leqi(lwfopt,'ordejon-mauri')) then
          ioptlwf = 2
        else
          call die('redata: wrong ON.funcional option')
        endif
      endif
#ifdef MPI
      call MPI_Bcast(lwfopt,13,MPI_character,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(ioptlwf,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
C ...

C     Option to calculate the Chemical potential in O(N)
      chebef = .false.
      noeta = .false.
      if (ionode) then
C     Option to use the Chemical Potential calculated instead
C     of the eta variable of the input
        noeta = fdf_boolean('ON.ChemicalPotentialUse',noeta)
        ! if so, we must (obviously) calculate the chemical potential
        if (noeta) then
          chebef=.true.
        else
          ! otherwise, we may still want to calculate it but not use it.
          chebef = fdf_boolean('ON.ChemicalPotential',chebef)
        endif
      endif
#ifdef MPI
      call MPI_Bcast(chebef,1,MPI_logical,0,MPI_Comm_World,
     .  MPIerror)
      call MPI_Bcast(noeta,1,MPI_logical,0,MPI_Comm_World,
     .  MPIerror)
#endif
C ...

C     Cutoff radius to calculate the Chemical Potential by projection
      rcoorcp = 9.5_dp
      if (ionode) then
        rcoorcp = fdf_physical('ON.ChemicalPotentialRc',
     .                        rcoorcp,'Bohr')
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(rcoorcp,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(rcoorcp,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif
C ...

C     Temperature of the Fermi distribution to calculate the
C     Chemical potential by projection
      tcp = 0.05_dp
      if (ionode) then
        tcp = fdf_physical('ON.ChemicalPotentialTemperature',
     .                   tcp,'Ry')
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(tcp,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(tcp,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif
      beta = 1.0_dp/tcp

C     Order of the Chebishev expansion to calculate the Chemical
C     potential
      pmax = 100
      if (ionode) then
        pmax = fdf_integer('ON.ChemicalPotentialOrder',pmax)
      endif
#ifdef MPI
      call MPI_Bcast(pmax,1,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
#endif
C ...


      if (isolve==1) then
        if (ionode) then
          write(6,'(a,i5)') 
     .       'redata: Maximum number of iterations     = ',ncgmax
          write(6,'(a,d12.2)') 
     .       'redata: Relative tolerance               = ',etol
          if (nspin.eq.2) then
            write(6,'(a,f10.4,a)') 
     .       'redata: Eta (Fermi level) Alpha spin     = ',eta(1),'  Ry'
            write(6,'(a,f10.4,a)') 
     .       'redata: Eta (Fermi level) Beta spin      = ',eta(2),'  Ry'
          else
            write(6,'(a,f10.4,a)') 
     .       'redata: Eta (Fermi level parameter)      = ',eta(1),'  Ry'
          endif
          write(6,'(a,f10.4,a)') 
     .       'redata: Radius of LWFs                   = ',rcoor,
     .       '  Bohr'
          write(6,'(a,4x,l1)') 
     .       'redata: Use continuation files for LWF   = ',
     .         usesavelwf
          write(6,'(a,a)') 
     .       'redata: Method to build LWFs             =     ',lwfopt
          if (chebef) then
            write(6,'(a,l1)')
     .       'redata: Compute Chemical Potential       =     ',chebef
            write(6,'(a)')
     .       'redata: Use the calculated Chemical ..'
            write(6,'(a,l1)')
     .       'redata: ..Potential instead of eta       =     ',noeta
            write(6,'(a,f10.4,a)') 
     .       'redata: Radius to compute the Chem. Pot. = ',rcoorcp,
     .       '  Bohr'
            write(6,'(a)')
     .       'redata: Temp. for Fermi distribution ..'
            write(6,'(a,f10.4,a)') 
     .       'redata: .. to compute the Chem. Pot.     = ',tcp,
     .       '    Ry'
            write(6,'(a,i5)') 
     .       'redata: Order of the Chebishev expansion = ',pmax
          endif
        endif
        if (cml_p) then
          call cmlAddParameter(xf=mainXML, name='ON.MaxNumIter',
     .         value=ncgmax, dictref='siesta:ncgmax')
          call cmlAddParameter(xf=mainXML, name='ON.etol',
     .         value=etol, dictref='siesta:etol')
          if (nspin==2) then
            call cmlAddParameter(xf=mainXML, name='ON.eta_alpha',
     .           value=eta(1), dictref='siesta:eta1', units='Ry')
            call cmlAddParameter(xf=mainXML, name='ON.eta_beta',
     .           value=eta(2), dictref='siesta:eta2', units='Ry')
          else
            call cmlAddParameter(xf=mainXML, name='ON.eta',
     .           value=eta(1), dictref='siesta:eta', units='Ry')
          endif
          call cmlAddParameter(xf=mainXML, name='On.RcLWF',
     .         value=rcoor, dictref='siesta:rcoor', units='Bohr')
          call cmlAddParameter(xf=mainXML, name='On.UseSaveLWF',
     .         value=usesavelwf, dictref='siesta:usesavelwf')
          call cmlAddParameter(xf=mainXML, name='ON.functional',
     .         value=lwfopt, dictref='siesta:lwfopt')
          if (chebef) then
            call cmlAddParameter(xf=mainXML, 
     .           name='ON.ChemicalPotential', value=chebef, 
     .           dictref='siesta:chebef')
            call cmlAddParameter(xf=mainXML, 
     .           name='ON.ChemicalPotentialUse', value=noeta, 
     .           dictref='siesta:noeta')
            call cmlAddParameter(xf=mainXML,
     .           name='ON.ChemicalPotentialRc', value=rcoorcp, 
     .           dictref='siesta:rcoorcp', units='Bohr')
            call cmlAddParameter(xf=mainXML,
     .           name='ON.ChemicalPotentialTemperature',
     .           value=tcp, dictref='siesta:tcp', units='Ry')
            call cmlAddParameter(xf=mainXML,
     .           name='ON.ChemicalPotentialOrder', value=pmax,
     .         dictref='siesta:pmax')
          endif !chebef
        endif !cml_p
      endif !Order-N
C ...

C Dynamics parameters ...
      varcel = .false.
      if (ionode) then
        varcel = fdf_boolean('MD.VariableCell', varcel )
      endif
! NB reset below ...
C ...

C     Type of dynamics 
      dyntyp = 'verlet'   
      if (ionode) then
        dyntyp = fdf_string('MD.TypeOfRun',dyntyp)
      endif
C ...

      if (ionode) then
! Two ways to identify a phonon calculation ...
        if (leqi(dyntyp,'Phonon').or.fdf_defined('MD.ATforPhonon')) then
        idyn = 7
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    PHONON forces calculation'
          if (cml_p) call cmlAddParameter(xf=mainXML, 
     .         name='MD.TypeOfRun', value='Phonon')
        else if (leqi(dyntyp,'cg')) then
          idyn = 0
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    CG coord. optimization'
          write(6,'(a,4x,l1)')
     .     'redata: Variable cell                    = ', varcel
          usesavecg = .false.
          usesavecg = fdf_boolean('UseSaveData', usesavecg)
          usesavecg = fdf_boolean('MD.UseSaveCG', usesavecg)
          write(6,'(a,4x,l1)')
     .     'redata: Use continuation files for CG    = ',
     .     usesavecg
          if (cml_p) then
            call cmlAddParameter(xf=mainXML, 
     .           name='MD.TypeOfRun', value='CG')
            call cmlAddParameter(xf=mainXML,
     .           name='MD.UseSaveCG', value=usesavecg)
          endif
        else if (leqi(dyntyp,'verlet')) then
          idyn = 1
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    Verlet MD run'
          if (cml_p) call cmlAddParameter(xf=mainXML,
     .         name='MD.TypeOfRun', value='Verlet')
        else if (leqi(dyntyp,'nose')) then
          idyn = 2
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    Nose thermostat MD run'
          if (cml_p) call cmlAddParameter(xf=mainXML,
     .         name='MD.TypeOfRun', value='Nose')
        else if (leqi(dyntyp,'parrinellorahman')) then
          idyn = 3
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    Parrinello-Rahman MD run'
          if (cml_p) call cmlAddParameter(xf=mainXML,
     .         name='MD.TypeOfRun', value='Parrinello-Rahman')
        else if (leqi(dyntyp,'noseparrinellorahman')) then
          idyn = 4
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    Nose-Parrinello-Rahman MD run'
          if (cml_p) call cmlAddParameter(xf=mainXML,
     .         name='MD.TypeOfRun', value='Nose-Parrinello-Rahman')
        else if (leqi(dyntyp,'anneal')) then
          idyn = 5
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    Annealing MD run'
          if (cml_p) call cmlAddParameter(xf=mainXML,
     .         name='MD.TypeOfRun', value='Annealing')
        else if (leqi(dyntyp,'fc')) then
          idyn = 6
          write(6,'(a,a)')
     .     'redata: Dynamics option                  = ',
     .     '    Force Constants Matrix Calculation'
          if (cml_p) call cmlAddParameter(xf=mainXML,
     .         name='MD.TypeOfRun', value='Force Constants')
        else if (leqi(dyntyp,'forces')) then
          idyn = 8
          write(6,'(a,a)')
     .     'redata: Dynamics option                  = ',
     .     '    Force evaluation'
          if (cml_p) call cmlAddParameter(xf=mainXML,
     .         name='MD.TypeOfRun', value='Force Evaluation')
        else
          call die('Invalid Option selected - value of MD.TypeOfRun'//
     .         'not recognised')
        endif
      endif
#ifdef MPI
      call MPI_Bcast(dyntyp,22,MPI_character,0,MPI_Comm_World,
     .     MPIerror)
      call MPI_Bcast(idyn,1,MPI_integer,0,MPI_Comm_World,
     .     MPIerror)
      call MPI_Bcast(varcel,1,MPI_logical,0,MPI_Comm_World,
     .     MPIerror)
      call MPI_Bcast(usesavecg,1,MPI_logical,0,MPI_Comm_World,
     .     MPIerror)
#endif
C ...


C     Maximum number of steps in CG coordinate optimization
      nmove = 0
      if (ionode) then
        nmove = fdf_integer('MD.NumCGsteps',nmove)
      endif
#ifdef MPI
      call MPI_Bcast(nmove,1,MPI_integer,0,MPI_Comm_World,
     .    MPIerror)
#endif
C ...

C     Maximum atomic displacement in one CG step
      dxmax = 0.2_dp
      if (ionode) then
        dxmax = fdf_physical('MD.MaxCGDispl',dxmax,'Bohr')
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(dxmax,1,MPI_double_precision,0,MPI_Comm_World,
     .    MPIerror)
#else
      call MPI_Bcast(dxmax,1,DAT_double,0,MPI_Comm_World,
     .    MPIerror)
#endif
#endif
C ...

C     Tolerance in the maximum atomic force (def 0.04 eV/Ang)
      ftol = 0.00155574_dp
      if (ionode) then
        ftol = fdf_physical('MD.MaxForceTol', ftol, 'Ry/Bohr')
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(ftol,1,MPI_double_precision,0,MPI_Comm_World,
     .    MPIerror)
#else
      call MPI_Bcast(ftol,1,DAT_double,0,MPI_Comm_World,
     .    MPIerror)
#endif
C ...
#endif

C     Tolerance in the maximum residual stress (var cell) def = 1 GPa 
      strtol = 6.79773e-5_dp
      if (ionode) then
        strtol = fdf_physical('MD.MaxStressTol', strtol, 'Ry/Bohr**3')
        strtol = abs(strtol)
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(strtol,1,MPI_double_precision,0,MPI_Comm_World,
     .    MPIerror)
#else
      call MPI_Bcast(strtol,1,DAT_double,0,MPI_Comm_World,
     .    MPIerror)
#endif
#endif
C ...

      if (idyn==0) then
        if (ionode) then
          write(6,'(a,i5)') 
     .         'redata: Maximum number of CG moves       = ', nmove
          write(6,'(a,f10.4,a)') 
     .         'redata: Max atomic displ per move        = ',
     .         dxmax,'  Bohr'
          write(6,'(a,f10.4,a)') 
     .         'redata: Force tolerance                  = ',
     .         ftol,'  Ry/Bohr'
          if (varcel) 
     .    write(6,'(a,f10.4,a)')
     .           'redata: Stress tolerance                 = ',
     .           strtol/6.79773e-5_dp, '  GPa'
        endif
        if (cml_p) then
          call cmlAddParameter(xf=mainXML,
     .         name='MD.NumCGSteps', value=nmove)
          call cmlAddParameter(xf=mainXML,
     .         name='MD.MaxCGDispl', value=dxmax, units='Bohr')
          call cmlAddParameter(xf=mainXML,
     .         name='MD.MaxForceTol', value=ftol, units='Ry/Bohr')
          if (varcel) call cmlAddParameter(xf=mainXML,
     .         name='MD.MaxStressTol', value=strtol, units='Ry/Bohr**3')
        endif ! cml_p
      endif ! idyn==0
C ...
  
C     Initial time step for MD
      istart = 1
      if (ionode) then
        istart = fdf_integer('MD.InitialTimeStep',istart)
      endif
#ifdef MPI
      call MPI_Bcast(istart,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
C ...

C     Final time step for MD
      ifinal = 1
      if (ionode) then
        ifinal = fdf_integer('MD.FinalTimeStep',ifinal)
      endif
#ifdef MPI
      call MPI_Bcast(ifinal,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
C ...

C     Length of time step for MD
      dt = 1._dp
      if (ionode) then
        dt = fdf_physical('MD.LengthTimeStep',dt,'fs')
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(dt,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(dt,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif
C ...

C     Quench Option
      qnch = .false.
      if (ionode) then
        qnch = fdf_boolean('MD.Quench',qnch)
        if (qnch .and. (idyn==2 .or. idyn==4)) 
     .    call die('redata: ERROR: You cannot quench and '//
     .         'use a Nose thermostat simultaneously')
        if (qnch) then
          iquench = 1
        else
          iquench = 0
        endif
      endif
#ifdef MPI
      call MPI_Bcast(iquench,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
C ...

C     Initial Temperature of MD simulation
C     (draws random velocities from the Maxwell-Boltzmann distribition
C      at the given temperature)
      tempinit = 0._dp
      if (ionode) then
        tempinit = fdf_physical('MD.InitialTemperature',tempinit,'K')
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(tempinit,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(tempinit,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif

      if (idyn .ge. 1 .and. idyn .le. 5) then
        if (ionode) then
          write(6,'(a,i5)') 
     .      'redata: Initial MD time step             = ',istart
          write(6,'(a,i5)') 
     .      'redata:   Final MD time step             = ',ifinal
          write(6,'(a,f10.4,a)') 
     .      'redata: Length of MD time step           = ',dt,'  fs'
          write(6,'(a,f10.4,a)') 
     .      'redata: Initial Temperature of MD run    = ',tempinit,'  K'
          if (idyn .ne. 5) write(6,'(a,4x,l1)') 
     .      'redata: Perform a MD quench              = ',qnch
        endif
        if (cml_p) then
          call cmlAddParameter(xf=mainXML, name='MD.InitialTimeStep',
     .         value=istart)
          call cmlAddParameter(xf=mainXML, name='MD.FinalTimeStep',
     .         value=ifinal)
          call cmlAddParameter(xf=mainXML, name='MD.LengthTimeStep',
     .         value=dt, units='fs')
          call cmlAddParameter(xf=mainXML, name='MD.InitialTemperature',
     .         value=tempinit, units='K')
          if (idyn/=5) call cmlAddParameter(xf=mainXML, 
     .         name='MD.Quench', value=qnch)
        endif
      endif

C     Target Temperature
      tt = 0._dp
      if (ionode) then
        tt = fdf_physical('MD.TargetTemperature',tt,'K')
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(tt,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(tt,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif
      
C     Target Pressure
      tp = 0._dp
      if (ionode) then
        tp = fdf_physical('MD.TargetPressure',tp,'Ry/Bohr**3')
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(tp,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(tp,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif

C     Mass of Nose variable
      mn = 1.e2_dp
      if (ionode) then
        mn = fdf_physical('MD.NoseMass',mn,'Ry*fs**2')
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(mn,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(mn,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif

C     Mass of Parrinello-Rahman variables
      mpr = 1.e2_dp
      if (ionode) then
        mpr = fdf_physical('MD.ParrinelloRahmanMass',
     .                    mpr,'Ry*fs**2')
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(mpr,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(mpr,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif

      if (idyn==2 .or. idyn==4) then
        if (ionode) write(6,'(a,f10.4,a)') 
     .  'redata: Nose mass                        = ',mn,'  Ry/fs**2'
        if (cml_p) call cmlAddParameter(xf=mainXML,
     .       name='MD.NoseMass', value=mn, units='Ry/fs**2')
      endif

      if (idyn==3 .or. idyn==4) then
        if (ionode) write(6,'(a,f10.4,a)') 
     .  'redata: Parrinello-Rahman mass           = ',mpr,'  Ry/fs**2'
        if (cml_p) call cmlAddParameter(xf=mainXML,
     .       name='MD.ParrinelloRahmanMass', value=mn, units='Ry/fs**2')
      endif

C     Annealing option
      annop = 'TemperatureAndPressure'
      if (ionode) then
        annop = fdf_string('MD.AnnealOption',annop)
      endif
#ifdef MPI
      call MPI_Bcast(annop,22,MPI_character,0,MPI_Comm_World,MPIerror)
#endif
      if (ionode) then
        if (idyn .eq. 5) then
          if (leqi(annop,'Temperature')) then
            ianneal = 1
            write(6,'(a,a)') 
     .       'redata: Annealing Option                 = ',
     .       'Temperature'
            if (cml_p) call cmlAddParameter(xf=mainXML,
     .           name='MD.AnnealOption', value='Temperature')
          else if (leqi(annop,'Pressure')) then
            ianneal = 2
            write(6,'(a,a)') 
     .       'redata: Annealing Option                 = ',
     .       'Pressure'
            if (cml_p) call cmlAddParameter(xf=mainXML,
     .           name='MD.AnnealOption', value='Pressure')
          else if (leqi(annop,'TemperatureAndPressure')) then
            ianneal = 3
            write(6,'(a,a)') 
     .       'redata: Annealing Option                 = ',
     .       'Temperature and Pressure'
            if (cml_p) call cmlAddParameter(xf=mainXML,
     .           name='MD.AnnealOption', value='TemperatureAndPressure')
          else
            call die('redata: ERROR: With annealing MD, you must '//
     .           'choose an appropriate value for MD.AnnealOption')
          endif !annop
        endif !idyn
      endif !ionode
      if (idyn==2 .or. idyn==4 .or. 
     .   (idyn==5 .and. (ianneal ==1 .or. ianneal==3))) then
        if (ionode) write(6,'(a,f10.4,a)') 
     .    'redata: Target Temperature               = ',tt,'  Kelvin'
        if (cml_p) call cmlAddParameter(xf=mainXML, 
     .       name='MD.TargetTemperature', value=tt, units='K')
      endif

      if (idyn==3 .or. idyn==4 .or. 
     .   (idyn==5 .and. (ianneal==2 .or. ianneal==3))) then
        if (ionode) write(6,'(a,f10.4,a)') 
     .       'redata: Target Pressure                  = ',
     .       tp, '  Ry/Bohr**3'
        if (cml_p) call cmlAddParameter(xf=mainXML, 
     .       name='MD.TargetPressure', value=tp, units='Ry/Bohr**3')
      endif

C     Relaxation Time for Annealing
      taurelax = 1.e2_dp
      if (ionode) then
        taurelax = fdf_physical('MD.TauRelax',taurelax,'fs')
        if (idyn==5) then
          write(6,'(a,f10.4,a)') 
     .         'redata: Annealing Relaxation Time        = ',
     .         taurelax,'  fs'
          if (cml_p) call cmlAddParameter(xf=mainXML, 
     .         name='MD.TauRelax', value=taurelax, units='fs')
        endif
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(taurelax,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(taurelax,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif
        
C     Estimated Bulk modulus (for Pressure annealing)
      bulkm = 1.e2_dp
      if (ionode) then
        bulkm = fdf_double('MD.BulkModulus',bulkm)
        if (idyn==5 .and. (ianneal==2 .or. ianneal==3)) then
          write(6,'(a,f10.4,a)') 
     .    'redata: Approx. Bulk Modulus             = ',
     .     bulkm,'  Ry/Bohr**3'
          if (cml_p) call cmlAddParameter(xf=mainXML, 
     .         name='MD.BulkModulus', value=bulkm, units='Ry/Bohr**3')
        endif
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(bulkm,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(bulkm,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif

C     Atomic displacement for force constant calculation
      dx = 0.04_dp
      if (ionode) then
        dx = fdf_physical('MD.FCDispl',dx,'Bohr')
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(dx,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(dx,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif

C     First and last atoms to displace for calculation of force constants
      ia1 = 1
      ia2 = na
      if (ionode) then
        ia1 = fdf_integer('MD.FCfirst',ia1)
        ia2 = fdf_integer('MD.FClast',ia2)
      endif
#ifdef MPI
      call MPI_Bcast(ia1,1,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
      call MPI_Bcast(ia2,1,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
#endif

      if (idyn==6) then
        if (ionode) then
          write(6,'(a,f10.4,a)')
     .  'redata: Atomic displ for force constants  = ',dx,'  Bohr'
          write(6,'(a,i8)')
     .  'redata: First atom to move               = ',ia1
          write(6,'(a,i8)')
     .  'redata: Last atom to move                = ',ia2
        endif
        if (cml_p) then
          call cmlAddParameter(xf=mainXML, name='MD.FCDispl',
     .         value=dx, units='Bohr')
          call cmlAddParameter(xf=mainXML, name='MD.FCFirst',
     .         value=ia1)
          call cmlAddParameter(xf=mainXML, name='MD.FCLast',
     .         value=ia2)
        endif
      endif
C ...

C Variable cell shape? Depending on input and type of dynamics
      varcel = varcel .or. (idyn==3) .or. (idyn==4) 
     .                .or. (idyn==5 .and. ianneal==1)
     .                .and. (idyn/=1) .and. (idyn/=2) 
     .                .and. (idyn/=6) .and. (idyn/=7)
     .                .and. (.not. (idyn==5 .and. ianneal/=1) )

      if (ionode) then
        write(6,'(2a)') 'redata: ', repeat('*', 71)
      endif

C Harris Forces?. Then DM.UseSaveDM should be ever false (use always
C Harris density in the first SCF step of each MD step), and
C MaxSCFIter should be ever 2, in the second one the SCF 
C Iteration are computed.

      if (ionode) then
         harrisfun = fdf_boolean('Harris_functional',.false.)
      endif
#ifdef MPI
      call MPI_Bcast(harrisfun,1,MPI_logical,0,MPI_Comm_World,
     $               MPIerror)
#endif

      if (harrisfun) then
        usesavedm = .false.
        nscf = 2
        mix = .false.
      endif

      if (cml_p) then
        call cmlEndParameterList(mainXML)
      endif

      end subroutine redata
