      module m_redata

      public :: redata
      private
      CONTAINS

      subroutine redata(na, ns, nspin, outlng, g2max,
     .  charnet, negl, nscf, dDtol, dEtol, mix, wmix,
     $  isolve, 
     .  temp, fixspin, ts, ncgmax, ftol, strtol, eta, 
     .  etol, rcoor, 
     .  ioptlwf, chebef, noeta, rcoorcp, beta, pmax,
     .  idyn, istart, ifinal, nmove, ianneal, iquench,
     .  dt, ia1, ia2, dx, dxmax, tt, tp, mn, mpr, 
     .  bulkm, taurelax, 
     .  usesavelwf, usesavedm, usesavecg,
     .  mullipop, inspn, maxsav, nkick, wmixkick, 
     .  pulfile, tempinit, dumpcharge, varcel, 
     .  harrisfun, occtol, broyden_maxit, require_energy_convergence,
     $  broyden_optim)
C *********************************************************************
C Subroutine to read the data for the SIESTA program
C
C     It uses the FDF (Flexible Data Format) package 
C     of J.M.Soler and A.Garcia
C
C Writen by P.Ordejon, December'96
C Modified for introduction of dynamic memory in SIESTA by JDG Sept 99
! Wrapping of most fdf and broadcast calls: A. Garcia, June 2005
!
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
C real*8 dDtol             : Maximum Density Matrix tolerance in SCF
C real*8 dEtol             : Maximum energy tolerance in SCF
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
C integer broyden_maxit    : Number of histories saved in Broyden mixing
C logical require_energy_convergence  : Impose E. conv. criterion?
C logical broyden_optim    : Broyden for forces instead of CG
C **********************************************************************

C
C  Modules
C
      use precision, only : dp
      use parallel,  only : IOnode, Nodes
      use fdf
      use sys
      use m_fdf_global, only: fdf_global_get
      use m_mpi_utils, only: broadcast
      use units,     only : eV
      use diagmemory, only: memoryfactor

      use siesta_cml

      implicit none

! Subroutine arguments:
      integer, intent(in)  :: na, nspin

      integer, intent(out) :: 
     .  ia1, ia2, ianneal, idyn, ifinal, ioptlwf,
     .  iquench, isolve, istart, maxsav, broyden_maxit,
     .  mullipop, ncgmax, nkick, nmove, ns, nscf, pmax

      real(dp), intent(out) ::
     .  beta, bulkm, charnet, dEtol,
     .  dDtol, dt, dx, dxmax, eta(2), etol, ftol, g2max,
     .  mn, mpr, occtol, rcoor, rcoorcp, strtol,
     .  taurelax, temp, tempinit, tp, ts, tt, wmix, wmixkick
      logical, intent(out) ::
     .  chebef, dumpcharge, fixspin, inspn, mix, negl, noeta, 
     .  outlng, pulfile, usesavecg, usesavelwf, usesavedm, varcel,
     .  harrisfun, require_energy_convergence, broyden_optim

C na_diag      : maximum number of atoms with diagon as default method
C g2max_default : Mesh cutoff default, in Ry
C temp_default  : Electronic temperature default, in Ry

      integer, parameter  :: na_diag = 100
      real(dp), parameter :: g2max_default = 100.e0_dp
      real(dp), parameter :: temp_default  = 1.900e-3_dp 
!
! Other defaults
!
      integer, parameter :: maxsav_default = 0
      integer, parameter :: nscf_default = 50
      integer, parameter :: ncgmax_default = 1000

      real(dp), parameter :: wmix_default = 0.25_dp
      real(dp), parameter :: wmixkick_default = 0.5_dp
      real(dp), parameter :: dDtol_default = 1.0e-4_dp
      real(dp), parameter :: dEtol_default = 1.0e-4_dp * eV
      real(dp), parameter :: occtol_default = 1.0e-12_dp
      real(dp), parameter :: etol_default = 1.0e-8_dp
      real(dp), parameter :: rcoor_default = 9.5_dp
      real(dp), parameter :: rcoorcp_default = 9.5_dp
      real(dp), parameter :: tcp_default = 0.05_dp
      integer, parameter  :: pmax_default = 100

      real(dp), parameter :: dxmax_default = 0.2_dp ! Bohr
      real(dp), parameter :: ftol_default =  0.00155574_dp ! Ry/Bohr 
                                                           ! 0.04 eV/Ang
      real(dp), parameter :: strtol_default = 6.79773e-5_dp ! 1 GPa
      real(dp), parameter :: dt_default = 1.0_dp ! 1 fs
      real(dp), parameter :: mn_default = 1.0e2_dp ! Nose mass in Ry*fs**2
      real(dp), parameter :: mpr_default = 1.0e2_dp ! PR mass in Ry*fs**2
      real(dp), parameter :: taurelax_default = 1.0e2_dp ! fs
      real(dp), parameter :: bulkm_default = 100*6.79773e-5_dp ! 100 GPa
      real(dp), parameter :: dx_default = 0.04_dp ! Bohr


C  Internal variables .................................................
      real(dp) :: tcp

      character
     .  sname*100, slabel*30,
     .  annop*22,  dyntyp*22, 
     .  method*6,  lwfopt*13, method_default*6

      logical  ::  DaC, leqi, qnch, usesaveddata
      external ::  leqi

! New template, using fdf and broadcast wrappers
!
!      call fdf_global_get(param,'ParamName', param_default)
!      if (ionode)  write(6,'(a,i)'),
!     .    'redata: ParamName           = ',param
!      if (cml_p) call cmlAddParameter(xf=mainXML, name='ParamName',
!     .                 value=param, dictref='siesta:param')

!
!      cml_p is only true in the master node
!
C Start of code
      if (cml_p) then
        call cmlStartParameterList(mainXML, title='Input Parameters')
      endif

C for cml output, find the system name & label

      if (cml_p) then
        sname = fdf_string('SystemName','')
        slabel = fdf_string('SystemLabel','')
        call cmlAddParameter(xf=mainXML, name='SystemName',
     .      value=trim(sname), dictref='siesta:sname')
        call cmlAddParameter(xf=mainXML, name='SystemLabel',
     .      value=trim(slabel), dictref='siesta:slabel')
      endif

C Type of output
      call fdf_global_get(outlng,"LongOutput", .false.)
      if (ionode)      write(6,'(a,4x,l1)')
     .   'redata: Long output                      = ',outlng
      if (cml_p) call cmlAddParameter(xf=mainXML, name='LongOutput',
     .       value=outlng, dictRef='siesta:verbosity')

C Define Number of species ...
      call fdf_global_get(ns,'NumberOfSpecies', 0)
      if (ionode)     write(6,'(a,3x,i2)') 
     .   'redata: Number of Atomic Species         = ',ns
        if (ns .le. 0) call die('redata: ERROR:'// 
     .       'Number of species must be larger than zero.')
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='NumberOfSpecies', value=ns, dictRef='siesta:ns')

C ...

C Dump information to plot charge contours
C by the external DENCHAR application program.

      call fdf_global_get(dumpcharge,'WriteDenchar',.false.)
      if (ionode)  write(6,'(a)') 
     .   'redata: Charge density info will appear in .RHO file'
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='WriteDenChar', value=dumpcharge)

C ...

C Perform Mulliken Population Analysis

      call fdf_global_get(mullipop,'WriteMullikenPop', 0)
      if (mullipop == 0 .and. outlng) mullipop = 1
      if (ionode) then
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
C ...

C Planewave cutoff of the real space mesh ...
      call fdf_global_get(g2max,'MeshCutoff',g2max_default,'Ry')
      if (ionode)   write(6,'(a,f10.4,a)') 
     .   'redata: Mesh Cutoff                      = ',g2max,'  Ry'
      if (cml_p) call cmlAddParameter(xf=mainXML, name='MeshCutOff',
     .     value=g2max, dictRef='siesta:g2max', units='Ry')
C ...

C Net charge in the cell ...
      call fdf_global_get(charnet,'NetCharge',0.0_dp)
      if (ionode)  write(6,'(a,f10.4,a)') 
     .   'redata: Net charge of the system         = ',charnet,' |e|'
      if (cml_p) call cmlAddParameter(xf=mainXML, name='NetCharge',
     .   value=charnet, dictRef='siesta:NetCharge', units='|e|')
C ...
       
C SCF Loop parameters ...
C     Maximum number of SCF iterations

      call fdf_global_get(nscf,'MaxSCFIterations',nscf_default)
      if (ionode)  write(6,'(a,i5)') 
     .   'redata: Max. number of SCF Iter          = ',nscf
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='MaxSCFIterations', value=nscf, dictRef='siesta:maxscf')

C Pulay mixing, number of iterations for one Pulay mixing (maxsav)
      call fdf_global_get(maxsav,'DM.NumberPulay', maxsav_default)

C Broyden mixing, number of iterations 
      call fdf_global_get(broyden_maxit,'DM.NumberBroyden',0)
      if (ionode) then
        if (broyden_maxit .gt. 0) then
           write(6,'(a,i5,a)') 
     .     'redata: Broyden mixing with ', broyden_maxit,
     $          ' saved histories.'
           if (maxsav > 1) then
              write(6,'(a)')
     $             'redata: Broyden supersedes Pulay!'
              maxsav = maxsav_default
           endif
        else if  (maxsav .gt. 1) then
          write(6,'(a,i5,a)') 
     .     'redata: One Pulay mixing every           = ',maxsav,
     .     ' iterations'
        else
           write(6,'(a)')'redata: Mixing is linear'
        endif
      endif
      if (cml_p) then
         call cmlAddParameter(xf=mainXML, name='DM.NumberPulay',
     .     value=maxsav, dictRef='siesta:maxsav')
         call
     $     cmlAddParameter(xf=mainXML, name='DM.NumberBroyden',
     .     value=broyden_maxit, dictRef='siesta:broyden_maxit')
      endif
C ...

C     Mix density matrix on first SCF step
C     (mix)

      call fdf_global_get(mix,'DM.MixSCF1',.false.)
      if (ionode)   write(6,'(a,4x,l1)')
     .    'redata: Mix DM in first SCF step ?       = ',mix
      if (cml_p) call cmlAddParameter(xf=mainXML, name='DM.MixSCF1',
     .     value=mix, dictRef='siesta:mix')


C Use disk or memory to store intermediate Pulay miximg vectors
C (pulfile)
      call fdf_global_get(pulfile,'DM.PulayOnFile',.false.)
      if (ionode) then
        if (pulfile.and.Nodes>1) then
          call die('redata: Cannot use DM.PulayOnFile=.true.'//
     .         'when running in parallel')
        endif
        write(6,'(a,4x,l1)')
     .    'redata: Write Pulay info on disk?        = ',pulfile
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, name='DM.PulayOnFile',
     .     value=pulfile, dictRef='siesta:pulfile')
C ...

C Density Matrix Mixing  (proportion of output DM in new input DM)
      call fdf_global_get(wmix,'DM.MixingWeight',wmix_default)
      if (ionode)  write(6,'(a,f10.4,a)') 
     .   'redata: New DM Mixing Weight             = ',wmix
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='DM.MixingWeight', value=wmix, dictRef='siesta:wmix')

C Density Matrix occupancy tolerance
      call fdf_global_get(occtol,
     $                    'DM.OccupancyTolerance',occtol_default)
      if (ionode)  write(6,'(a,f14.12)') 
     .   'redata: New DM Occupancy tolerance       = ',occtol
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='DM.OccupancyTolerance', value=occtol, 
     .     dictRef='siesta:occtol')


C Perform linear mixing each nkick SCF iterations (to kick system
C when it is pinned in a poorly convergent SCF loop)
      call fdf_global_get(nkick,'DM.NumberKick',0)
      if (ionode) then
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
 
C     Density Matrix Mixing each nkick SCF iterations
      call fdf_global_get( wmixkick,
     $                  'DM.KickMixingWeight',wmixkick_default)
      if (ionode) write(6,'(a,f10.4,a)')
     .   'redata: DM Mixing Weight for Kicks       = ',wmixkick
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='DM.KickMixingWeight', value=wmixkick, 
     .     dictRef='siesta:wmixkick')


C     Density Matrix Tolerance for achieving Self-Consistency
      call fdf_global_get(dDtol,'DM.Tolerance',dDtol_default)
      if (ionode)  write(6,'(a,f12.6,a)') 
     .   'redata: DM Tolerance for SCF             = ',dDtol
      if (cml_p) call cmlAddParameter(xf=mainXML, name='DM.Tolerance', 
     .     value=dDtol, dictRef='siesta:dDtol')

C  Require Energy convergence for achieving Self-Consistency?
      call fdf_global_get(require_energy_convergence,
     $     'DM.RequireEnergyConvergence',.false.)
      if (ionode)  write(6,'(a,4x,l1)') 
     .   'redata: Require Energy convergence for SCF = ',
     $            require_energy_convergence
      if (cml_p) call cmlAddParameter(xf=mainXML,
     $     name='DM.RequireEnergyConvergence', 
     .     value=require_energy_convergence,
     $     dictRef='siesta:ReqEnergyConv')

C Energy tolerance for achieving Self-Consistency
      call fdf_global_get(dEtol,'DM.EnergyTolerance',dEtol_default,"Ry")
      if (ionode)  write(6,'(a,f12.6,a)') 
     .   'redata: DM Energy tolerance for SCF      = ',dEtol/eV,' eV'
      if (cml_p) call cmlAddParameter(xf=mainXML,
     $     name='DM.EnergyTolerance', 
     .     value=dEtol, dictRef='siesta:dEtol')


C Initial spin density: Maximum polarization, Ferro (false), AF (true)
      if (nspin.eq.2) then
         call fdf_global_get(inspn,'DM.InitSpinAF',.false.)
        if (ionode)  write(6,'(a,4x,l1)')
     .     'redata: Antiferro initial spin density   = ',inspn
        if (cml_p) call cmlAddParameter(xf=mainXML, 
     .       name='DM.InitSpinAF', value=inspn, dictRef='siesta:inspn')
      endif
C ...

C Use Saved Data
      call fdf_global_get(usesaveddata,'UseSaveData',.false.)
      if (ionode)  write(6,'(a,4x,l1)') 
     .   'redata: Using Saved Data (generic)   = ', usesaveddata

C ...
C Use continuation DM files
      call fdf_global_get(usesavedm,'DM.UseSaveDM',usesaveddata)
      if (ionode) write(6,'(a,4x,l1)') 
     .   'redata: Use continuation files for DM    = ',  usesavedm
      if (cml_p) call cmlAddParameter(xf=mainXML, name='DM.UseSaveDM', 
     .     value=usesavedm, dictRef='siesta:usesavedm')
C ...

C Neglect Interactions between non-overlapping orbitals ...
      call fdf_global_get(negl,'NeglNonOverlapInt',.false.)
      if (ionode)  write(6,'(a,4x,l1)') 
     .   'redata: Neglect nonoverlap interactions  = ',negl
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='NeglNonOverlapInt', value=negl, dictRef='siesta:negl')
C ...

C Method to Solve LDA Hamiltonian ...
      if (na .le. na_diag) then
         method_default = 'diagon'
      else
         method_default = 'ordern'
      endif
      call fdf_global_get(method,'SolutionMethod',method_default)
      if (cml_p)
     $    call cmlAddParameter(xf=mainXML, name='SolutionMethod',
     $       value=method, dictRef='siesta:SCFmethod')

C
      if (leqi(method,'diagon')) then
         isolve = 0
C        DivideAndConquer is now the default
         call fdf_global_get(DaC,'Diag.DivideAndConquer',.true.)
         if (ionode)  then
            write(6,'(a,4x,a)') 
     .           'redata: Method of Calculation            = ',
     .           'Diagonalization'
            write(6,'(a,4x,l1)') 
     .           'redata: Divide and Conquer               = ',DaC
         endif

      else if (leqi(method,'ordern')) then
         isolve = 1
         DaC = .false.
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
      if (cml_p)
     $   call cmlAddParameter(xf=mainXML, name='Diag.DivideAndConquer',
     $     value=DaC, dictRef='siesta:DaC')

C ...

C Memory scaling factor for rdiag/cdiag - cannot be less than 1.0
      call fdf_global_get(MemoryFactor,'Diag.Memory', 1.0_dp )
      MemoryFactor = max(MemoryFactor,1.0_dp)
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='Diag.Memory', value=MemoryFactor, 
     .     dictRef='siesta:MemoryFactor')
C ...

C Electronic temperature for Fermi Smearing ...
      call fdf_global_get(temp,
     $                 'ElectronicTemperature',temp_default,'Ry')
      if (ionode .and. isolve == 0)  write(6,'(a,f10.4,a)') 
     .    'redata: Electronic Temperature           = ',temp,'  Ry'
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='ElectronicTemperature', value=temp, 
     .     dictRef='siesta:etemp', units='Ry')
C ...

C Fix the spin of the system to a given value ; and
C value of the Spin of the system (only used if fixspin = TRUE)

      call fdf_global_get(fixspin,'FixSpin',.false.)
      if (ionode)   write(6,'(a,4x,l1)') 
     .         'redata: Fix the spin of the system       = ',fixspin 
      if (fixspin) then
         if (nspin .ne. 2) 
     $      call die('redata: ERROR: You can only fix the spin of '//
     .         'the system for collinear spin polarized calculations.')

         call fdf_global_get(ts,'TotalSpin',0.0_dp)
         if (ionode) write(6,'(a,f10.4)') 
     .        'redata: Value of the Spin of the System  = ',ts
      else
         ts = 0.0_dp
      endif

      if (cml_p) then 
        call cmlAddParameter(xf=mainXML, name='FixSpin',
     .       value=fixspin, dictref='siesta:fixspin')
        call cmlAddParameter(xf=mainXML, name='TotalSpin',
     .       value=ts, dictref='siesta:ts')
      endif
C ...

C Order-N solution parameters ...
C     Maximum number of CG minimization iterations

      call fdf_global_get(ncgmax,'ON.MaxNumIter',ncgmax_default)
      if (ncgmax<1) then
      if (ionode)
     $   write(6,'(a)') "ON.MaxNumIter cannot be less than 1. "//
     .                  "Resetting to 1"
         ncgmax = 1 
      endif

C     Relative tolerance in total band structure energy

        call fdf_global_get(etol,'ON.etol',etol_default)
C ...

C     Fermi level parameter
      eta(1:2) = 0.0_dp
      if (ionode) then
        eta(1) = fdf_physical('ON.eta',eta(1),'Ry')
        eta(2) = eta(1)
        eta(1) = fdf_physical('ON.eta_alpha',eta(1),'Ry')
        eta(2) = fdf_physical('ON.eta_beta',eta(2),'Ry')
      endif
      call broadcast(eta(1:2))
C ...

C     Cutoff radius for Localized Wave Functions

      call fdf_global_get(rcoor,'On.RcLWF',rcoor_default,'Bohr')

C     Use continumation LWF files

      call fdf_global_get(usesavelwf,'ON.UseSaveLWF',usesaveddata)
C ...

C     Option on how to build LWF's (disk or functionals)
      call fdf_global_get(lwfopt,'ON.functional',"kim")
      if (leqi(lwfopt,'files')) then
         ioptlwf = 0
      else if (leqi(lwfopt,'kim')) then
         ioptlwf = 1
      else if (leqi(lwfopt,'ordejon-mauri')) then
         ioptlwf = 2
      else
         call die('redata: wrong ON.funcional option')
      endif

C     Option to calculate the Chemical potential in O(N)
C     Option to use the Chemical Potential calculated instead
C     of the eta variable of the input

      call fdf_global_get(noeta,'ON.ChemicalPotentialUse',.false.)
      if (noeta) then
        ! if so, we must (obviously) calculate the chemical potential
         chebef=.true.
      else
        ! otherwise, we may still want to calculate it but not use it.
         call fdf_global_get(chebef,'ON.ChemicalPotential',.false.)
      endif


C     Cutoff radius to calculate the Chemical Potential by projection

      call fdf_global_get(rcoorcp,'ON.ChemicalPotentialRc',
     .                        rcoorcp_default,'Bohr')
C ...

C     Temperature of the Fermi distribution to calculate the
C     Chemical potential by projection

      call fdf_global_get(tcp,'ON.ChemicalPotentialTemperature',
     .                   tcp_default,'Ry')
      beta = 1.0_dp/tcp

C     Order of the Chebishev expansion to calculate the Chemical
C     potential
      
      call fdf_global_get(pmax,
     $           'ON.ChemicalPotentialOrder',pmax_default)
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

      call fdf_global_get(varcel,'MD.VariableCell', .false. )
! NB reset below ...
C ...

C     Type of dynamics 
      call fdf_global_get(dyntyp,'MD.TypeOfRun',"verlet")
C ...

      if (leqi(dyntyp,'cg')) then
         idyn = 0
         call fdf_global_get(usesavecg,'MD.UseSaveCG', usesaveddata)
      else if (leqi(dyntyp,'verlet')) then
         idyn = 1
      else if (leqi(dyntyp,'nose')) then
         idyn = 2
      else if (leqi(dyntyp,'parrinellorahman')) then
         idyn = 3
      else if (leqi(dyntyp,'noseparrinellorahman')) then
         idyn = 4
      else if (leqi(dyntyp,'anneal')) then
         idyn = 5
      else if (leqi(dyntyp,'fc')) then
         idyn = 6
      else if (leqi(dyntyp,'phonon')) then
         idyn = 7
      else if (leqi(dyntyp,'forces')) then
         idyn = 8
      else
         call die('Invalid Option selected - value of MD.TypeOfRun'//
     .        'not recognised')
      endif

C     Broyden or Conjugate Gradients...

      call fdf_global_get(broyden_optim,'Optim.Broyden',.false.)

C     Maximum number of steps in CG coordinate optimization
      
      call fdf_global_get(nmove,'MD.NumCGsteps',0)

C     Maximum atomic displacement in one CG step

      call fdf_global_get(dxmax,'MD.MaxCGDispl',dxmax_default,'Bohr')

C     Tolerance in the maximum atomic force 

      call fdf_global_get(ftol,
     $            'MD.MaxForceTol', ftol_default, 'Ry/Bohr')

C     Tolerance in the maximum residual stress (var cell) def = 1 GPa 

      call fdf_global_get(strtol,
     $               'MD.MaxStressTol', strtol_default, 'Ry/Bohr**3')
      strtol = abs(strtol)

      if (ionode) then
        select case (idyn)
          case(0)
             write(6,'(a,a)') 
     .            'redata: Dynamics option                  = ',
     .            '    CG coord. optimization'
             write(6,'(a,4x,l1)')
     .            'redata: Variable cell                    = ', varcel
             write(6,'(a,4x,l1)')
     .            'redata: Use continuation files for CG    = ',
     .             usesavecg
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
             if (cml_p) then
                call cmlAddParameter(xf=mainXML, 
     .               name='MD.TypeOfRun', value='CG')
                call cmlAddParameter(xf=mainXML,
     .               name='MD.UseSaveCG', value=usesavecg)
                call cmlAddParameter(xf=mainXML,
     .               name='MD.NumCGSteps', value=nmove)
                call cmlAddParameter(xf=mainXML,
     .               name='MD.MaxCGDispl', value=dxmax, units='Bohr')
                call cmlAddParameter(xf=mainXML,
     .               name='MD.MaxForceTol', value=ftol, units='Ry/Bohr')
                if (varcel) call cmlAddParameter(xf=mainXML,
     .               name='MD.MaxStressTol',
     $               value=strtol, units='Ry/Bohr**3')
             endif
         case(1)
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    Verlet MD run'
          if (cml_p) call cmlAddParameter(xf=mainXML,
     .         name='MD.TypeOfRun', value='Verlet')
          case(2)
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    Nose thermostat MD run'
          if (cml_p) call cmlAddParameter(xf=mainXML,
     .         name='MD.TypeOfRun', value='Nose')

          case(3)
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    Parrinello-Rahman MD run'
          if (cml_p) call cmlAddParameter(xf=mainXML,
     .         name='MD.TypeOfRun', value='Parrinello-Rahman')
          case(4)
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    Nose-Parrinello-Rahman MD run'
          if (cml_p) call cmlAddParameter(xf=mainXML,
     .         name='MD.TypeOfRun', value='Nose-Parrinello-Rahman')
          case(5)
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    Annealing MD run'
          if (cml_p) call cmlAddParameter(xf=mainXML,
     .         name='MD.TypeOfRun', value='Annealing')

          case(6)
          write(6,'(a,a)')
     .     'redata: Dynamics option                  = ',
     .     '    Force Constants Matrix Calculation'
          if (cml_p) call cmlAddParameter(xf=mainXML,
     .         name='MD.TypeOfRun', value='Force Constants')

          case(7)
          write(6,'(a,a)') 
     .     'redata: Dynamics option                  = ',
     .     '    PHONON forces calculation'
          if (cml_p) call cmlAddParameter(xf=mainXML, 
     .         name='MD.TypeOfRun', value='Phonon')
          
          case(8)
          write(6,'(a,a)')
     .     'redata: Dynamics option                  = ',
     .     '    Force evaluation'
          if (cml_p) call cmlAddParameter(xf=mainXML,
     .         name='MD.TypeOfRun', value='Force Evaluation')

       end select
      endif

      if (ionode) write(6,'(a,l1)') 
     .     'redata: Broyden Optimization             = ', broyden_optim
      if (cml_p) call cmlAddParameter(xf=mainXML,
     .       name='Optim.Broyden', value=broyden_optim)
  
C     Initial and final time steps for MD
       call fdf_global_get(istart,'MD.InitialTimeStep',1)
       call fdf_global_get(ifinal,'MD.FinalTimeStep',1)

C     Length of time step for MD

       call fdf_global_get(dt,'MD.LengthTimeStep',dt_default,'fs')
C ...

C     Quench Option
       call fdf_global_get(qnch,'MD.Quench',.false.)
       if (qnch .and. (idyn==2 .or. idyn==4)) 
     .    call die('redata: ERROR: You cannot quench and '//
     .         'use a Nose thermostat simultaneously')
       if (qnch) then
          iquench = 1
       else
          iquench = 0
       endif
C ...

C     Initial Temperature of MD simulation
C     (draws random velocities from the Maxwell-Boltzmann distribition
C      at the given temperature)

       call fdf_global_get(tempinit,
     $             'MD.InitialTemperature',0.0_dp,'K')


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

C     Target Temperature and Pressure
      call fdf_global_get(tt,'MD.TargetTemperature',0.0_dp,'K')
      call fdf_global_get(tp,'MD.TargetPressure',0.0_dp,'Ry/Bohr**3')
      

C     Mass of Nose variable
      call fdf_global_get(mn,'MD.NoseMass',mn_default,'Ry*fs**2')

C     Mass of Parrinello-Rahman variables
      call fdf_global_get(mpr,
     $             'MD.ParrinelloRahmanMass',mpr_default,'Ry*fs**2')

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
      ianneal = 0
      call fdf_global_get(annop,
     $           'MD.AnnealOption','TemperatureAndPressure')

      if (idyn .eq. 5) then
         if (leqi(annop,'Temperature')) then
            ianneal = 1
         else if (leqi(annop,'Pressure')) then
            ianneal = 2
         else if (leqi(annop,'TemperatureAndPressure')) then
            ianneal = 3
         else
            call die('redata: ERROR: With annealing MD, you must '//
     .           'choose an appropriate value for MD.AnnealOption')
         endif                  !annop

         if (ionode) then
            select case (ianneal)
            case(1)
               write(6,'(a,a)') 
     .              'redata: Annealing Option                 = ',
     .              'Temperature'
               if (cml_p) call cmlAddParameter(xf=mainXML,
     .              name='MD.AnnealOption', value='Temperature')
            case(2)
               write(6,'(a,a)') 
     .              'redata: Annealing Option                 = ',
     .              'Pressure'
               if (cml_p) call cmlAddParameter(xf=mainXML,
     .              name='MD.AnnealOption', value='Pressure')
            case(3)
               write(6,'(a,a)') 
     .              'redata: Annealing Option                 = ',
     .              'Temperature and Pressure'
               if (cml_p) call cmlAddParameter(xf=mainXML,
     .          name='MD.AnnealOption', value='TemperatureAndPressure')
            end select
         endif                  !ionode
      endif                     !idyn


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
      call fdf_global_get(taurelax,
     $             'MD.TauRelax',taurelax_default,'fs')
      if (idyn==5) then
          if (ionode) write(6,'(a,f10.4,a)') 
     .         'redata: Annealing Relaxation Time        = ',
     .         taurelax,'  fs'
          if (cml_p) call cmlAddParameter(xf=mainXML, 
     .         name='MD.TauRelax', value=taurelax, units='fs')
       endif
        
C     Estimated Bulk modulus (for Pressure annealing)

      call fdf_global_get(bulkm,
     $              'MD.BulkModulus',bulkm_default,'Ry/Bohr**3')
      if (ionode) then
         if (idyn==5 .and. (ianneal==2 .or. ianneal==3))
     $            write(6,'(a,f10.4,a)') 
     .           'redata: Approx. Bulk Modulus             = ',
     .           bulkm,'  Ry/Bohr**3'
      endif
      if (cml_p) call cmlAddParameter(xf=mainXML, 
     .     name='MD.BulkModulus', value=bulkm, units='Ry/Bohr**3')

C     Atomic displacement for force constant calculation

      call fdf_global_get(dx,'MD.FCDispl',dx_default,'Bohr')

C     First and last atoms to displace for calculation of force constants

      call fdf_global_get(ia1,'MD.FCfirst',1)
      call fdf_global_get(ia2,'MD.FClast',na)

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


C Harris Forces?. Then DM.UseSaveDM should be false (use always
C Harris density in the first SCF step of each MD step), and
C MaxSCFIter should be  2, in the second one the SCF 
C Iteration are computed.

      call fdf_global_get(harrisfun,'Harris_functional',.false.)

      if (harrisfun) then
        usesavedm = .false.
        nscf = 2
        mix = .false.
      endif

      if (ionode)   write(6,'(2a)') 'redata: ', repeat('*', 71)
      if (cml_p)    call cmlEndParameterList(mainXML)

      end subroutine redata
      end module m_redata
