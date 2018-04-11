! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
subroutine read_options( na, ns, nspin )

  ! Subroutine to read the options for the SIESTA program
  !
  !     It uses the FDF (Flexible Data Format) package 
  !     of J.M.Soler and A.Garcia
  !
  ! Writen by P.Ordejon, December'96
  ! Modified for introduction of dynamic memory in SIESTA by JDG Sept 99
  ! Wrapping of most fdf and broadcast calls: A. Garcia, June 2005
  !

  use siesta_options
  use precision, only : dp
  use parallel,  only : IOnode, Nodes, ParallelOverK
  use fdf
  use files,     only : slabel
  use files,     only : filesOut_t   ! derived type for output file names
  use sys
  use units,     only : eV
  use diagmemory,   only: memoryfactor
  use siesta_cml
  use m_target_stress, only: set_target_stress

  implicit none
  !----------------------------------------------------------- Input Variables
  ! integer na               : Number of atoms
  ! integer ns               : Number of species
  ! integer nspin            : Spin polarization

  integer, intent(in)  :: na, ns, nspin

  ! This routine sets variables in the 'siesta_options' module

  ! g2max_default : Mesh cutoff default, in Ry
  ! temp_default  : Electronic temperature default, in Ry

  real(dp), parameter :: g2cut_default = 100.e0_dp
  real(dp), parameter :: temp_default  = 1.900e-3_dp 

  logical, parameter  :: mixH_def = .false.

  integer,  parameter :: maxsav_default = 0
  integer,  parameter :: nscf_default = 50
  integer,  parameter :: ncgmax_default = 1000

  real(dp), parameter :: wmix_default = 0.25_dp
  real(dp), parameter :: wmixkick_default = 0.5_dp
  real(dp), parameter :: dDtol_default = 1.0e-4_dp
  real(dp), parameter :: Energy_tolerance_default = 1.0e-5_dp * eV  ! Free energy...
  real(dp), parameter :: Harris_tolerance_default = 1.0e-5_dp * eV
  real(dp), parameter :: occtol_default = 1.0e-12_dp
  real(dp), parameter :: etol_default = 1.0e-8_dp
  real(dp), parameter :: rcoor_default = 9.5_dp
  real(dp), parameter :: rcoorcp_default = 9.5_dp
  real(dp), parameter :: tcp_default = 0.05_dp
  integer,  parameter :: pmax_default = 100

  real(dp), parameter :: dxmax_default = 0.2_dp             ! Bohr
  real(dp), parameter :: ftol_default =  0.00155574_dp      ! Ry/Bohr 
  ! 0.04 eV/Ang
  real(dp), parameter :: strtol_default = 6.79773e-5_dp     ! 1 GPa
  real(dp), parameter :: dt_default = 1.0_dp                ! 1 fs
  real(dp), parameter :: mn_default = 1.0e2_dp              ! Nose mass in Ry*fs**2
  real(dp), parameter :: mpr_default = 1.0e2_dp             ! PR mass in Ry*fs**2
  real(dp), parameter :: taurelax_default = 1.0e2_dp        ! fs
  real(dp), parameter :: bulkm_default = 100*6.79773e-5_dp  ! 100 GPa
  real(dp), parameter :: dx_default = 0.04_dp               ! Bohr


  ! The following are comment lines that should be merged into 'siesta_options'.

  ! real*8 charnet           : Net charge (in units of |e|)
  ! logical outlng           : Long (true) or short (false) output
  ! real*8 g2cut             : PW cutoff energy (Ry)
  ! logical negl             : True = Neglect interactions between
  !                            non-overlaping orbitals (coming from
  !                            KB projectors)
  ! integer nscf             : Maximum number of SCF cycles per time step
  ! real*8 dDtol             : Maximum Density Matrix tolerance in SCF
  ! real*8 Energy_tolerance  : Maximum Total energy tolerance in SCF
  ! real*8 Harris_tolerance  : Maximum Harris energy tolerance in SCF
  ! logical mix              : Perform mix in first SCF step
  ! real*8 wmix              : Amount of output DM for new DM
  ! integer isolve           : Method of solution.  0   = Diagonalization
  !                                                 1   = Order-N
  !                                                 2   = Transiesta
  !                                                 3   = OMM
  ! real*8 temp              : Temperature for Fermi smearing (Ry)
  ! logical fixspin          : Fix the spin of the system?
  ! real*8  ts               : Total spin of the system
  ! integer ncgmax           : Maximum number of CG steps for 
  !                            band structure energy minimization
  ! real*8 etol              : Relative tolerance in CG minimization
  !                            of band structure energy
  ! real*8 eta(2)            : Fermi level parameter of Kim functional
  ! real*8 rcoor             : Cutoff radius of LWF's (Bohr)
  ! integer ioptlwf          : Option to build LWF's according to:
  !                             0 = Read blindly from disk
  !                             1 = Functional of Kim et al.
  !                             2 = Functional of Ordejon-Mauri
  ! logical chebef          : Compute the chemical potential 
  ! logical noeta            : Use computed Chem.pot. instead of eta
  ! real*8 rcoorcp           : Cutoff (Bohr) to compute the chem.pot.
  ! real*8 beta              : Inverse temperature to compute chem.pot.
  ! integer pmax             : Order of Chebi expansion for chem.pot.
  ! integer idyn             : Atomic dynamics option:
  !                             0 = CG geometry optimization
  !                             1 = Standard MD run (Verlet)
  !                             2 = Nose thermostat MD
  !                             3 = Parrinello-Rahman MD
  !                             4 = Nose thermostat + Parrinello-Rahman MD
  !                             5 = Annealing MD
  !                             6 = Force constants
  !                             7 = Deprecated (Forces for PHONON program)
  !                             8 = Force evaluation
  ! integer istart           : Initial time step for MD
  ! integer ifinal           : Final time step for MD
  ! integer nmove            : Number of CG steps in CG optimization
  ! real*8 ftol              : Maximum force for CG structure optimization
  ! real*8 strtol            : Maximum stress for CG structure optimization
  ! integer ianneal          : Annealing option for idyn = 5
  !                             1 = Temperature 
  !                             2 = Pressure
  !                             3 = Temperature and Pressure
  ! integer iquench          : Quench option: 0 = No; 1 = Yes; 2 = Fire
  ! real*8 dt                : Length of time step (fs)
  ! real*8 dx                : Atomic displacement for Force Constants
  !                             calculation
  ! integer ia1              : First atom to displace for force constants
  ! integer ia2              : Last atom to displace for force constants
  ! real*8 dxmax             : Maximum atomic displacement in one CG move
  ! real*8 tt                : Target temperature (Kelvin)
  ! real*8 tp                : Target Pressure (Ry/Bohr**3)
  ! real*8 mn                : Mass of Nose variable (Ry/fs**2)
  ! real*8 mpr               : Mass of Parrinello-R. variable (Ry/fs**2)
  ! real*8 bulkm             : Estimate of bulk modulus (Ry/Bohr**3)
  ! real*8 taurelax          : Annealing time to reach targer T and P (fs)
  ! logical usesavelwf       : True = try to use continuation LWF files 
  !                              from disk
  ! logical usesavedm        : True = try to use continuation DM files 
  !                              from disk
  ! logical usesavecg        : True = try to use continuation CG files
  !                              from disk
  ! integer mullipop         : Option for Mulliken Pop. analysis
  ! logical inspn            : Spin initialization for spin-polarized
  !                              .true.  -> Antiferro
  !                              .false. -> Ferro
  ! integer maxsav           : Number of density-matrices stored for Pulay
  !                            mixing. .lt.2 => linear mixing only
  !                                    .ge.2 => pulay mixing
  ! integer nkick            : Perform a linear mixing eack nkick scf cycles
  ! real*8 wmixkick          : Mixing parameter for linear mixing each nkick scf
  !                            cycles
  ! logical pulfile          : Use file (.true.) or memory (.false.)
  !                            to store Pulay miximg intermediate vectors
  !                            Default: .false.
  ! real*8 tempinit          : Initial temperature (Kelvin) of the MD simulation
  ! logical dumpcharge       : True: Dump information to plot charge contours
  !                            by the external DENCHAR application program.
  !     (This is now obsolete: info will appear in .RHO file)
  ! logical varcel           : variable shape for CG optimization or dynamics
  ! logical harrisfun        : swith that indicates if harris functional will
  !                            be used or not
  ! real*8  occtol           : Occupancy threshold for DM build
  ! integer broyden_maxit    : Number of histories saved in Broyden SCF mixing
  ! logical require_energy_convergence  : Impose E. conv. criterion?
  ! logical broyden_optim    : Broyden for forces instead of CG
  ! logical want_domain_decomposition:  Use domain decomposition for orbitals in O(N)
  ! logical want_spatial_decomposition:  Use spatial decomposition for orbitals in O(N)


  !----------------------------------------------------------- Local Variables
  real(dp) :: tcp

  character annop*22,  dyntyp*22,  method*6,  lwfopt*13

  logical  ::  DaC, qnch, qnch2, usesaveddata

  !--------------------------------------------------------------------- BEGIN
  ! New template, using fdf
  !
  !      param = fdf_get('ParamName', param_default)
  !      if (ionode)  write(6,'(a,i)'),
  !     .    'redata: ParamName           = ',param
  !      if (cml_p) call cmlAddParameter(xf=mainXML, name='ParamName',
  !     .                 value=param, dictref='siesta:param')

  !
  !      cml_p is only true in the master node
  !
  ! Start of code

  if (cml_p) then
     call cmlStartParameterList(mainXML, title='Input Parameters')
  endif

  ! for cml output, find the system name & label
  if (cml_p) then
     call cmlAddParameter(xf=mainXML, name='SystemName',             &
          value=trim(sname), dictref='siesta:sname')
     call cmlAddParameter(xf=mainXML, name='SystemLabel',            &
          value=trim(slabel), dictref='siesta:slabel')
  endif

  ! H setup only
  h_setup_only = fdf_get('HSetupOnly', .false.)
  if (ionode .and. h_setup_only) then
     write(6,1) 'redata: H Setup Only                     = ', h_setup_only
  endif

  ! Type of output
  outlng = fdf_get('LongOutput', .false.)
  if (ionode) then
     write(6,1) 'redata: Long output                      = ', outlng
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='LongOutput',            &
          value=outlng, dictRef='siesta:verbosity' )
  endif

  ! Write about Number of species, as before
  if (ionode) then
     write(6,4) 'redata: Number of Atomic Species         = ', ns
  endif

  if (ns .le. 0) then
     call die( 'redata: ERROR: Number of species must be larger than zero.' )
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, title='NumberOfSpecies', &
          value=ns, dictRef='siesta:ns', units="cmlUnits:countable" )
  endif

  ! Dump information to plot charge contours
  ! by the external DENCHAR application program.
  dumpcharge = fdf_get('WriteDenchar',.false.)
  if (ionode) then
     write(6,2) 'redata: Charge density info will appear in .RHO file'
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='WriteDenChar', value=dumpcharge)
  endif

  ! Perform Mulliken Population Analysis
  mullipop = fdf_get('WriteMullikenPop', 0)
  if (mullipop == 0 .and. outlng) then
     mullipop = 1
  endif

  if (ionode) then
     select case (mullipop)
     case(0)
        write(6,2) 'redata: Write Mulliken Pop.              =     NO'
     case(1)
        write(6,2) 'redata: Write Mulliken Pop.              =     '//&
             'Atomic and Orbital charges'
     case(2)
        write(6,'(a,a/45x,a)') 'redata: Write Mulliken Pop.              =     ',&
             'Atomic and Orbital charges', 'plus Atomic Overlap Populations'

     case(3)
        write(6,'(a,a/45x,a/45x,a)') 'redata: Write Mulliken Pop.              =     ',&
             'Atomic and Orbital charges','plus Atomic Overlap Populations',        &
             'plus Orbital Overlap Populations'

     case default
        call die( 'redata: Invalid value for WriteMullikenPop' )
     end select
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='WriteMullikenPop', value=mullipop, &
          units="cmlUnits:dimensionless" )
  endif

  ! Perform Hirshfeld and/or Voronoi Population Analysis
  hirshpop= fdf_get('WriteHirshfeldPop',.false.)
  voropop=  fdf_get('WriteVoronoiPop',.false.)
  partial_charges_at_every_geometry =  &
       fdf_get('PartialChargesAtEveryGeometry',.false.)
  partial_charges_at_every_scf_step =  &
       fdf_get('PartialChargesAtEveryScfStep',.false.)

  
  if ( fdf_get('Compat.Matel.NRTAB', .false.) ) then
    matel_NRTAB = 128
  else
    matel_NRTAB = 1024
  end if
  if ( IONode ) then
    write(6,4) 'redata: Matel table size (NRTAB)         = ', matel_NRTAB
  end if
  if (cml_p) then
    call cmlAddParameter( xf=mainXML, name='MatelNRTAB',value=matel_NRTAB, &
        dictRef='siesta:matel_nrtab', units="cmlUnits:countable")
  end if

  ! Planewave cutoff of the real space mesh ...
  g2cut = fdf_get('MeshCutoff',g2cut_default,'Ry')
  if (ionode) then
     write(6,6)'redata: Mesh Cutoff                      = ', g2cut,' Ry'
   endif
   

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='MeshCutOff', value=g2cut,     &
          dictRef='siesta:g2max', units='siestaUnits:Ry' )
  endif

  ! Net charge in the cell ...
  charnet = fdf_get('NetCharge',0.0_dp)
  if (ionode) then
     write(6,6) 'redata: Net charge of the system         = ',charnet,' |e|'
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='NetCharge', value=charnet, &
          dictRef='siesta:NetCharge', units='siestaUnits:e__')
  endif

  ! SCF Loop parameters ...
  !     Minimum/Maximum number of SCF iterations
  min_nscf = fdf_get('MinSCFIterations',0)
  nscf     = fdf_get('MaxSCFIterations',nscf_default)
  SCFMustConverge = fdf_get('SCFMustConverge', .false.)
  if (ionode) then
     write(6,4) 'redata: Min. number of SCF Iter          = ',min_nscf
     write(6,4) 'redata: Max. number of SCF Iter          = ',nscf
     if (SCFMustConverge) then
        write(6,4) 'redata: SCF convergence failure will abort job'
     endif
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='MaxSCFIterations',  &
          value=nscf, dictRef='siesta:maxscf',  &
          units="cmlUnits:countable")
     call cmlAddParameter( xf=mainXML, name='MinSCFIterations',  &
          value=min_nscf, dictRef='siesta:minscf',  &
          units="cmlUnits:countable")
  endif

  mixH = fdf_get('MixHamiltonian',mixH_def)
  mixH = fdf_get('TS.MixH',mixH)   ! Catch old-style keyword
  mix_charge = fdf_get('MixCharge',.false.)

  if (mix_charge) then
     if (ionode) then
        write(6,1) 'redata: Mix charge density rho_g         = ', mix_charge
     endif
     if (mixH) then
        mixH = .false.
        if (ionode) then
           write(6,"(a)") 'redata: ***MixCharge takes precedence over MixH'
        endif
     endif
  endif
  if (mixH) then
     if (ionode) then
        write(6,1) 'redata: Mix Hamiltonian instead of DM    = ', mixH
     endif
  endif

  ! Options for pre-4.0 compatibility
  compat_pre_v4_DM_H  = fdf_get('Compat-pre-v4-DM-H',.false.)
  mix_after_convergence = fdf_get('SCF.MixAfterConvergence',compat_pre_v4_DM_H)
  recompute_H_after_scf = fdf_get('SCF.Recompute-H-After-Scf',compat_pre_v4_DM_H)

  if (ionode) then
     if (compat_pre_v4_DM_H) then
        write(6,"(a)") ':!:Next two options activated by pre-4.0 compat. switch'
     endif
     write(6,1) 'redata: Mix DM or H after convergence    = ',  &
          mix_after_convergence
     write(6,1) 'redata: Recompute H after scf cycle      = ',  &
          recompute_H_after_scf
  endif

  ! Pulay mixing, number of iterations for one Pulay mixing (maxsav)
  maxsav = fdf_get('DM.NumberPulay', maxsav_default)

  ! Broyden SCF mixing, number of iterations 
  broyden_maxit = fdf_get('DM.NumberBroyden',0)
  ! FIRE SCF mixing, no parameters
  fire_mix = fdf_get('DM.FIRE.Mixing',.false.)
  if (ionode) then
     if (fire_mix) then
        write(6,*) "Fire Mixing"
     else if (broyden_maxit .gt. 0) then
        write(6,5) 'redata: Broyden mixing with ', &
             broyden_maxit, &
             ' saved histories.'
        if (maxsav > 1) then
           write(6,2) 'redata: Broyden supersedes Pulay!'
           maxsav = maxsav_default
        endif
     else if (maxsav .gt. 1) then
        write(6,5) 'redata: Performing Pulay mixing using    = ',maxsav, &
             ' iterations'
     else
        write(6,2)'redata: Mixing is linear'
     endif
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='DM.NumberPulay',     &
          value=maxsav, dictRef='siesta:maxsav', &
          units="cmlUnits:countable" )
     call cmlAddParameter( xf=mainXML, name='DM.NumberBroyden',    &
          value=broyden_maxit, dictRef='siesta:broyden_maxit', &
          units="cmlUnits:countable" )
  endif

  ! Mix density matrix on first SCF step
  ! (mix)
  mix = fdf_get('DM.MixSCF1',.false.)
  !
  if (ionode) then
     write(6,1) 'redata: Mix DM in first SCF step ?       = ',mix
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='DM.MixSCF1',   &
          value=mix, dictRef='siesta:mix' )
  endif

  ! Use disk or memory to store intermediate Pulay mixing vectors
  ! (pulfile)
  pulfile = fdf_get('DM.PulayOnFile',.false.)
  if (ionode) then
     if (pulfile) then
        call die( 'redata: Cannot use DM.PulayOnFile=.true.'//&
             'in this version' )
     endif
     write(6,1) 'redata: Write Pulay info on disk?        = ',pulfile
  endif
  if (cml_p) then
     call cmlAddParameter(xf=mainXML, name='DM.PulayOnFile',      &
          value=pulfile, dictRef='siesta:pulfile')
  endif

  ! 
  avoid_first_after_kick = fdf_get (    &
       'DM.Pulay.Avoid.First.After.Kick',.false.)
  if (ionode) then
     write(6,1) 'redata: Discard 1st Pulay DM after  kick = ', &
          avoid_first_after_kick
  endif
  if (cml_p) then
     call cmlAddParameter(xf=mainXML,  &
          name='DM.Pulay.Avoid.First.After.Kick', &
          value=pulfile,   &
          dictRef='siesta:avoid_first_after_kick')
  endif

  ! Density Matrix Mixing  (proportion of output DM in new input DM)
  wmix = fdf_get('DM.MixingWeight',wmix_default)
  if (ionode) then
     write(6,6) 'redata: New DM Mixing Weight             = ',wmix
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML,name='DM.MixingWeight', &
          value=wmix, dictRef='siesta:wmix', &
          units="cmlUnits:dimensionless" )
  endif

  ! Density Matrix occupancy tolerance
  occtol = fdf_get('DM.OccupancyTolerance',occtol_default)
  if (ionode) then
     write(6,8) 'redata: New DM Occupancy tolerance       = ',occtol
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML,name='DM.OccupancyTolerance', &
          value=occtol, dictRef='siesta:occtol' ,  &
          units="cmlUnits:dimensionless" )
  endif

  ! Perform linear mixing each nkick SCF iterations (to kick system
  ! when it is pinned in a poorly convergent SCF loop)
  nkick = fdf_get('DM.NumberKick',0)
  if (ionode) then
     if (nkick .ge. 1) then
        write(6,5) 'redata: Kick with linear mixing every    = ',nkick,&
             ' iterations'
     else
        write(6,2)'redata: No kicks to SCF'
     endif
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML,name='DM.NumberKick',     &
          value=nkick, dictRef='siesta:nkick', &
          units="cmlUnits:countable" )
  endif

  ! Density Matrix Mixing each nkick SCF iterations
  wmixkick = fdf_get('DM.KickMixingWeight',wmixkick_default)
  if (ionode) then
     write(6,6) 'redata: DM Mixing Weight for Kicks       = ',wmixkick
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML,name='DM.KickMixingWeight',    &
          value=wmixkick, dictRef='siesta:wmixkick',&
          units="cmlUnits:dimensionless" )
  endif

  ! Density Matrix Tolerance for achieving Self-Consistency
  dDtol = fdf_get('DM.Tolerance',dDtol_default)
  if (ionode) then
     write(6,7) 'redata: DM Tolerance for SCF             = ',dDtol
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='DM.Tolerance',     &
          value=dDtol, dictRef='siesta:dDtol', &
          units='siestaUnits:eAng_3' )
  endif
  !--------------------------------------

  ! Require Energy convergence for achieving Self-Consistency?
  require_energy_convergence = fdf_get('DM.RequireEnergyConvergence', &
       .false.)
  if (ionode) then
     write(6,1) 'redata: Require (free) Energy convergence in SCF = ', &
          require_energy_convergence
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='DM.RequireEnergyConvergence', &
          value=require_energy_convergence,               &
          dictRef='siesta:ReqEnergyConv' )
  endif

  ! Energy tolerance for achieving Self-Consistency
  Energy_tolerance = fdf_get('DM.EnergyTolerance',    &
       Energy_tolerance_default, 'Ry' )
  if (ionode) then
     write(6,7) 'redata: DM (free)Energy tolerance for SCF = ', Energy_tolerance/eV, ' eV'
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='DM.EnergyTolerance', &
          value=Energy_tolerance/eV, dictRef='siesta:dEtol', &
          units="siestaUnits:eV" )
  endif

  !--------------------------------------
  ! Require Harris Energy convergence for achieving Self-Consistency?
  require_harris_convergence = fdf_get('DM.RequireHarrisConvergence', .false.)
  if (ionode) then
     write(6,1) 'redata: Require Harris convergence for SCF = ', &
          require_harris_convergence
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='DM.RequireHarrisConvergence', &
          value=require_harris_convergence,               &
          dictRef='siesta:ReqHarrisConv' )
  endif

  ! Harris energy tolerance for achieving Self-Consistency
  Harris_tolerance = fdf_get('DM.HarrisTolerance',    &
       Harris_tolerance_default, 'Ry' )
  if (ionode) then
     write(6,7) 'redata: DM Harris energy tolerance for SCF = ', Harris_tolerance/eV, ' eV'
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='DM.HarrisTolerance', units='siestaUnits:eV', &
          value=Harris_tolerance, dictRef='siesta:Harris_tolerance')
  endif

  ! Monitor forces and stresses during SCF loop
  monitor_forces_in_scf = fdf_get('MonitorForcesInSCF',.false.)

  !--------------------------------------
  ! Initial spin density: Maximum polarization, Ferro (false), AF (true)
  if (nspin.eq.2) then
     inspn = fdf_get('DM.InitSpinAF',.false.)
     if (ionode) then
        write(6,1) 'redata: Antiferro initial spin density   = ',inspn
     endif
     if (cml_p) then
        call cmlAddParameter( xf=mainXML, name='DM.InitSpinAF',   &
             value=inspn, dictRef='siesta:inspn')
     endif
  endif

  ! Use Saved Data
  usesaveddata = fdf_get('UseSaveData',.false.)
  if (ionode) then
     write(6,1) 'redata: Using Saved Data (generic)   = ', usesaveddata
  endif

  ! Use continuation DM files
  usesavedm = fdf_get('DM.UseSaveDM',usesaveddata)
  if (ionode) then
     write(6,1) 'redata: Use continuation files for DM    = ',  usesavedm
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='DM.UseSaveDM',            &
          value=usesavedm, dictRef='siesta:usesavedm')
  endif

  ! Neglect Interactions between non-overlapping orbitals ...
  negl = fdf_get('NeglNonOverlapInt',.false.)
  if (ionode) then
     write(6,1) 'redata: Neglect nonoverlap interactions  = ',negl
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='NeglNonOverlapInt', &
          value=negl, dictRef='siesta:negl' )
  endif

  ! Method to Solve LDA Hamiltonian ...
  method = fdf_get('SolutionMethod','diagon')
  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='SolutionMethod',        &
          value=method, dictRef='siesta:SCFmethod' )
  endif

  if (leqi(method,'diagon')) then
     isolve = SOLVE_DIAGON
     ! DivideAndConquer is now the default
     DaC = fdf_get('Diag.DivideAndConquer',.true.)
     if (ionode)  then
        write(6,'(a,4x,a)') 'redata: Method of Calculation            = ',&
             'Diagonalization'
        write(6,1) 'redata: Divide and Conquer               = ', DaC
     endif

  else if (leqi(method,'ordern')) then
     isolve = SOLVE_ORDERN
     DaC    = .false.
     if (ionode) then
        write(6,'(a,4x,a)') 'redata: Method of Calculation            = ', &
             'Order-N'
     endif
     if (nspin .gt. 2) then
        call die( 'redata: You chose the Order-N solution option '// &
             'together with nspin>2.  This is not allowed in '//&
             'this version of siesta' )
     endif
  else if (leqi(method,'omm')) then
     isolve = SOLVE_MINIM
     DaC    = .false.
     call_diagon_default=fdf_integer('OMM.Diagon',0)
     call_diagon_first_step=fdf_integer('OMM.DiagonFirstStep',call_diagon_default)
     minim_calc_eigenvalues=fdf_boolean('OMM.Eigenvalues',.false.)
     if (ionode) then
        write(6,'(a,4x,a)') 'redata: Method of Calculation            = ', &
             'Orbital Minimization Method'
     endif
#ifdef TRANSIESTA
     ! TSS Begin
  else if (leqi(method,'transi')) then
     isolve = SOLVE_TRANSI
     if (ionode) then
        write(*,'(a,4x,a)')                                &
             'redata: Method of Calculation            = ',  &
             '    Transiesta'
     endif
#endif /* TRANSIESTA */
  else
     call die( 'redata: The method of solution must be either '//&
#ifdef TRANSIESTA
          'Transiesta, '//&
#endif
          'OrderN, or Diagon' )
  endif

#ifdef DEBUG
  call write_debug( '    Solution Method: ' // method )
#endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='Diag.DivideAndConquer', &
          value=DaC, dictRef='siesta:DaC' )
  endif

  ! Memory scaling factor for rdiag/cdiag - cannot be less than 1.0
  MemoryFactor = fdf_get('Diag.Memory', 1.0_dp )
  MemoryFactor = max(MemoryFactor,1.0_dp)
  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='Diag.Memory', &
          value=MemoryFactor,             &
          dictRef='siesta:MemoryFactor',  &
          units="cmlUnits:dimensionless" )
  endif

  ! Electronic temperature for Fermi Smearing ...
  temp = fdf_get('ElectronicTemperature',temp_default,'Ry')
  if (ionode .and. isolve.eq.SOLVE_DIAGON) then
     write(6,6) 'redata: Electronic Temperature           = ',temp,'  Ry'
  endif

  if (cml_p) then
     call cmlAddParameter( xf=mainXML, name='ElectronicTemperature', &
          value=temp, dictRef='siesta:etemp',       &
          units='siestaUnits:Ry')
  endif

  ! Fix the spin of the system to a given value ; and
  ! value of the Spin of the system (only used if fixspin = TRUE)
  fixspin = fdf_get('FixSpin',.false.)
  if (ionode) then
     write(6,1) 'redata: Fix the spin of the system       = ',fixspin 
  endif

  if (fixspin) then
     if (nspin .ne. 2) then
        call die( 'redata: ERROR: You can only fix the spin of '//&
             'the system for collinear spin polarized calculations.' )
     endif
     ts = fdf_get('TotalSpin',0.0_dp)
     if (ionode) then
        write(6,9) 'redata: Value of the Spin of the System  = ',ts
     endif
  else
     ts = 0.0_dp
  endif

  if (cml_p) then 
     call cmlAddParameter( xf=mainXML, name='FixSpin', &
          value=fixspin, dictref='siesta:fixspin' )
     call cmlAddParameter( xf=mainXML, name='TotalSpin', &
          value=ts, dictref='siesta:ts',&
          units='siestaUnits:eSpin' )
  endif

  ! Order-N solution parameters ...
  !     Maximum number of CG minimization iterations
  ncgmax = fdf_get('ON.MaxNumIter',ncgmax_default)
  if (ncgmax<1) then
     if (ionode) then
        write(6,2) 'ON.MaxNumIter cannot be less than 1.  Resetting to 1'
     endif
     ncgmax = 1 
  endif

  ! Relative tolerance in total band structure energy
  etol = fdf_get('ON.etol',etol_default)

  ! Fermi level parameter
  eta(1:2) = 0.0_dp
  eta(1) = fdf_physical('ON.eta',eta(1),'Ry')
  eta(2) = eta(1)
  eta(1) = fdf_physical('ON.eta_alpha',eta(1),'Ry')
  eta(2) = fdf_physical('ON.eta_beta',eta(2),'Ry')

  ! Cutoff radius for Localized Wave Functions
  rcoor = fdf_get('On.RcLWF',rcoor_default,'Bohr')

  ! Use continumation LWF files
  usesavelwf = fdf_get('ON.UseSaveLWF',usesaveddata)

  ! Option on how to build LWF's (disk or functionals)
  lwfopt = fdf_get('ON.functional','kim')
  if (leqi(lwfopt,'files')) then
     ioptlwf = 0
  else if (leqi(lwfopt,'kim')) then
     ioptlwf = 1
  else if (leqi(lwfopt,'ordejon-mauri')) then
     ioptlwf = 2
  else
     call die('redata: wrong ON.funcional option')
  endif

  ! Option to calculate the Chemical potential in O(N)
  ! Option to use the Chemical Potential calculated instead
  ! of the eta variable of the input
  noeta = fdf_get('ON.ChemicalPotentialUse',.false.)
  ! NOTE: This does not yet work in parallel

  if (noeta) then
     ! if so, we must (obviously) calculate the chemical potential
     chebef=.true.
  else
     ! otherwise, we may still want to calculate it but not use it.
     chebef = fdf_get('ON.ChemicalPotential',.false.)
  endif

#ifdef MPI
  if (chebef) then
     call die("ON.ChemicalPotential(Use) options do not work with MPI")
  endif
#endif

  ! Cutoff radius to calculate the Chemical Potential by projection
  rcoorcp = fdf_get( 'ON.ChemicalPotentialRc', &
       rcoorcp_default, 'Bohr' )

  ! Temperature of the Fermi distribution to calculate the
  ! Chemical potential by projection
  tcp = fdf_get( 'ON.ChemicalPotentialTemperature', &
       tcp_default,'Ry' )
  beta = 1.0_dp/tcp

  ! Order of the Chebishev expansion to calculate the Chemical potential
  pmax = fdf_get('ON.ChemicalPotentialOrder',pmax_default)


  if (isolve==SOLVE_ORDERN) then
     if (ionode) then
        write(6,4) 'redata: Maximum number of iterations     = ',ncgmax
        write(6,'(a,d12.2)') 'redata: Relative tolerance               = ',etol
        if (nspin.eq.2) then
           write(6,6) 'redata: Eta (Fermi level) Alpha spin     = ',eta(1),'  Ry'
           write(6,6) 'redata: Eta (Fermi level) Beta spin      = ',eta(2),'  Ry'
        else
           write(6,6) 'redata: Eta (Fermi level parameter)      = ',eta(1),'  Ry'
        endif
        write(6,6) 'redata: Radius of LWFs                   = ',rcoor,'  Bohr'
        write(6,1) 'redata: Use continuation files for LWF   = ',usesavelwf
        write(6,'(a,a)') 'redata: Method to build LWFs             =     ',lwfopt
        if (chebef) then
           write(6,1) 'redata: Compute Chemical Potential       = ',chebef
           write(6,2) 'redata: Use the calculated Chemical ..'
           write(6,1) 'redata: ..Potential instead of eta       = ',noeta
           write(6,6) 'redata: Radius to compute the Chem. Pot. = ',rcoorcp,'  Bohr'
           write(6,2) 'redata: Temp. for Fermi distribution ..'
           write(6,6) 'redata: .. to compute the Chem. Pot.     = ',tcp,'    Ry'
           write(6,4) 'redata: Order of the Chebishev expansion = ',pmax
        endif
     endif
     if (cml_p) then
        call cmlAddParameter( xf      = mainXML,        &
             name    = 'ON.MaxNumIter',&
             value   = ncgmax,         &
             dictref = 'siesta:ncgmax', &
             units   = "cmlUnits:countable" )

        call cmlAddParameter( xf      = mainXML,       &
             name    = 'ON.etol',     &
             value   = etol,          &
             dictref = 'siesta:etol', &
             units   = "siestaUnits:eV" )
        if (nspin==2) then
           call cmlAddParameter( xf      = mainXML,          &
                name    = 'ON.eta_alpha',   &
                value   = eta(1),           &
                dictref = 'siesta:eta1',    &
                units   = 'siestaUnits:Ry' )

           call cmlAddParameter( xf      = mainXML,          &
                name    = 'ON.eta_beta',    &
                value   = eta(2),           &
                dictref = 'siesta:eta2',    &
                units   = 'siestaUnits:Ry' )
        else
           call cmlAddParameter( xf      = mainXML,         &
                name    = 'ON.eta',        &
                value   = eta(1),          &
                dictref = 'siesta:eta',    &
                units   = 'siestaUnits:Ry')
        endif
        call cmlAddParameter( xf      = mainXML,          &
             name    = 'On.RcLWF',       &
             value   = rcoor,            &
             dictref = 'siesta:rcoor',   &
             units   = 'siestaUnits:Bohr')

        call cmlAddParameter( xf=mainXML,                 &
             name='On.UseSaveLWF',       &
             value=usesavelwf,           &
             dictref='siesta:usesavelwf' )

        call cmlAddParameter( xf      = mainXML,        &
             name    = 'ON.functional',&
             value   = lwfopt,         &
             dictref = 'siesta:lwfopt' )
        if (chebef) then
           call cmlAddParameter( xf      = mainXML,                &
                name    = 'ON.ChemicalPotential', &
                value   = chebef,                 &
                dictref = 'siesta:chebef')

           call cmlAddParameter( xf      = mainXML,                   &
                name    = 'ON.ChemicalPotentialUse', &
                value   = noeta,                     &
                dictref = 'siesta:noeta')

           call cmlAddParameter( xf      = mainXML,                  &
                name    = 'ON.ChemicalPotentialRc', &
                value   = rcoorcp,                  &
                dictref = 'siesta:rcoorcp',         &
                units   = 'siestaUnits:Bohr')

           call cmlAddParameter( xf    = mainXML,                           &
                name  = 'ON.ChemicalPotentialTemperature', &
                value = tcp, dictref='siesta:tcp',         &
                units = 'siestaUnits:Ry' )

           call cmlAddParameter( xf      = mainXML,                     &
                name    = 'ON.ChemicalPotentialOrder', &
                value   = pmax,                        &
                dictref = 'siesta:pmax',               &
                units   = 'cmlUnits:dimensionless')
        endif
     endif
  endif

  ! Dynamics parameters ...
  varcel = fdf_get('MD.VariableCell', .false. )

  ! NB reset below ...
  ! Type of dynamics 

  compat_pre_v4_dynamics = fdf_get('compat-pre-v4-dynamics', .false. )
  if (compat_pre_v4_dynamics) then
     dyntyp = fdf_get('MD.TypeOfRun','verlet')
  else
     dyntyp = fdf_get('MD.TypeOfRun','cg')
  endif

  if (leqi(dyntyp,'cg')) then
     idyn = 0
     usesavecg = fdf_get('MD.UseSaveCG', usesaveddata)
     ! Support the old Broyden switch  for now
     broyden_optim = fdf_get('Optim.Broyden',.false.)

     if (broyden_optim) then
        write(6,2) '**Note: FDF symbol Optim.Broyden is '//&
             'deprecated. (Accepted for now.) See manual.'
     endif

  else if (leqi(dyntyp,'broyden')) then
     idyn = 0
     broyden_optim = .true.
  else if (leqi(dyntyp,'fire')) then
     idyn = 0
     fire_optim = .true.
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
     call die('Dynamics type "PHONON" is no longer supported')
  else if (leqi(dyntyp,'forces').or.leqi(dyntyp,'master')) then
     idyn = 8
  else
     call die('Invalid Option selected - value of MD.TypeOfRun not recognised')
  endif

  ! Maximum number of steps in CG/Broyden coordinate optimization
  nmove = fdf_get('MD.NumCGsteps',0)

  ! Maximum atomic displacement in one CG step
  dxmax = fdf_get('MD.MaxCGDispl',dxmax_default,'Bohr')

  ! Tolerance in the maximum atomic force 
  ftol = fdf_get('MD.MaxForceTol', ftol_default, 'Ry/Bohr')

  ! Tolerance in the maximum residual stress (var cell) def = 1 GPa 
  strtol = fdf_get('MD.MaxStressTol', strtol_default, 'Ry/Bohr**3')
  strtol = abs(strtol)
  
  GeometryMustConverge = fdf_get('GeometryMustConverge', .false.)

  if (ionode) then
     select case (idyn)
     case(0)
        if (nmove > 0) then
           if (broyden_optim) then
              write(6,2) 'redata: Dynamics option                  =     '//&
                   'Broyden coord. optimization'
           elseif (fire_optim) then
              write(6,2) 'redata: Dynamics option                  =     '//&
                   'FIRE coord. optimization'
           else
              write(6,2) 'redata: Dynamics option                  =     '//&
                   'CG coord. optimization'
           endif
           write(6,1) 'redata: Variable cell                    = ', varcel
           if (.not. broyden_optim) then
              write(6,1) 'redata: Use continuation files for CG    = ', usesavecg
              write(6,6) 'redata: Max atomic displ per move        = ', dxmax,&
                   '  Bohr'
           endif
           write(6,4) 'redata: Maximum number of CG moves       = ', nmove
           write(6,6) 'redata: Force tolerance                  = ', ftol,&
                '  Ry/Bohr'
           if (varcel) then
              write(6,6) 'redata: Stress tolerance                 = ', &
                   strtol/6.79773e-5_dp, '  GPa'
           endif
           if (cml_p) then
              if (broyden_optim) then
                 call cmlAddParameter( xf   = mainXML,        &
                      name = 'MD.TypeOfRun', &
                      value= 'Broyden' )
              else if (fire_optim) then
                 call cmlAddParameter( xf   = mainXML,        &
                      name = 'MD.TypeOfRun', &
                      value= 'FIRE' )
              else
                 call cmlAddParameter( xf    =mainXML,        &
                      name  ='MD.TypeOfRun', &
                      value ='CG' )

                 call cmlAddParameter( xf    = mainXML,       &
                      name  = 'MD.UseSaveCG',&
                      value = usesavecg )
              endif

              call cmlAddParameter( xf    = mainXML,         &
                   name  = 'MD.NumCGSteps', &
                   value = nmove,           &
                   units = "cmlUnits:countable" )

              call cmlAddParameter( xf    = mainXML,           &
                   name  = 'MD.MaxCGDispl',   &
                   value = dxmax,             &
                   units = 'siestaUnits:Bohr' )

              call cmlAddParameter( xf=mainXML,                 &
                   name='MD.MaxForceTol',      &
                   value=ftol,                 &
                   units='siestaUnits:Ry_Bohr')
              if (varcel) then
                 call cmlAddParameter( xf=mainXML,                     &
                      name='MD.MaxStressTol',         &
                      value=strtol,                   &
                      units='siestaUnits:Ry_Bohr__3' )
              endif
           endif
        else
           write(6,2) 'redata: Dynamics option                  =     '//&
                'Single-point calculation'
           if (cml_p) then
              call cmlAddParameter( xf   = mainXML,        &
                   name = 'MD.TypeOfRun', &
                   value= 'Single-Point' )
           endif
        endif
     case(1)
        write(6,2) 'redata: Dynamics option                  =     '//&
             'Verlet MD run'
        if (cml_p) then
           call cmlAddParameter( xf    = mainXML,        &
                name  = 'MD.TypeOfRun', &
                value = 'Verlet' )
        endif

     case(2)
        write(6,2) 'redata: Dynamics option                  =     '//&
             'Nose thermostat MD run'
        if (cml_p) then
           call cmlAddParameter( xf=mainXML, name='MD.TypeOfRun', value='Nose')
        endif

     case(3)
        write(6,2) 'redata: Dynamics option                  =     '//&
             'Parrinello-Rahman MD run'
        if (cml_p) then
           call cmlAddParameter( xf= mainXML,              &
                name= 'MD.TypeOfRun',     &
                value= 'Parrinello-Rahman' )
        endif

     case(4)
        write(6,2) 'redata: Dynamics option                  =     '//&
             'Nose-Parrinello-Rahman MD run'
        if (cml_p) then
           call cmlAddParameter( xf    = mainXML,                &
                name  = 'MD.TypeOfRun',         &
                value = 'Nose-Parrinello-Rahman' )
        endif

     case(5)
        write(6,2) 'redata: Dynamics option                  =     '//&
             'Annealing MD run'
        if (cml_p) then
           call cmlAddParameter( xf    = mainXML,        &
                name  = 'MD.TypeOfRun', &
                value = 'Annealing' )
        endif

     case(6)
        write(6,2) 'redata: Dynamics option                  =     '//&
             'Force Constants Matrix Calculation'
        if (cml_p) then
           call cmlAddParameter( xf    = mainXML,          &
                name  = 'MD.TypeOfRun',   &
                value = 'Force Constants' )
        endif

     case(7)
        ! deprecated

     case(8)
        write(6,2) 'redata: Dynamics option                  =     Force evaluation'
        if (cml_p) then
           call cmlAddParameter( xf    = mainXML,            &
                name  = 'MD.TypeOfRun',     &
                value = 'Force Evaluation' )
        endif
     end select
  endif

  ! Initial and final time steps for MD
  istart = fdf_get('MD.InitialTimeStep',1)
  ifinal = fdf_get('MD.FinalTimeStep',1)

  ! Length of time step for MD
  dt = fdf_get('MD.LengthTimeStep',dt_default,'fs')

  ! Quench Option
  qnch  = fdf_get('MD.Quench',.false.)
  qnch2 = fdf_get('MD.FireQuench',.false.)
  if ((qnch .or. qnch2) .and. (idyn==2 .or. idyn==4)) then 
     call die( 'redata: ERROR: You cannot quench and '//&
          'use a Nose thermostat simultaneously')
  endif

  iquench = 0
  if (qnch) then
     iquench = 1
  endif
  if (qnch2) then
     iquench = 2
  endif

  ! Initial Temperature of MD simulation
  ! (draws random velocities from the Maxwell-Boltzmann distribition
  !  at the given temperature)
  tempinit = fdf_get('MD.InitialTemperature',0.0_dp,'K')

  if (idyn .ge. 1 .and. idyn .le. 5) then
     if (ionode) then
        write(6,4) 'redata: Initial MD time step             = ',istart
        write(6,4) 'redata:   Final MD time step             = ',ifinal
        write(6,6) 'redata: Length of MD time step           = ',dt,'  fs'
        write(6,6) 'redata: Initial Temperature of MD run    = ',tempinit,'  K'
        if (idyn .ne. 5) then 
           if (qnch2) then
              write(6,1) 'redata: Perform a MD Fire quench         = ',qnch2
           else
              write(6,1) 'redata: Perform a MD quench              = ',qnch
           endif
        endif
     endif

     if (cml_p) then
        call cmlAddParameter( xf    = mainXML,             &
             name  = 'MD.InitialTimeStep',&
             value = istart,              &
             units = 'cmlUnits:countable' )

        call cmlAddParameter( xf    = mainXML,             &
             name  = 'MD.FinalTimeStep',  &
             value = ifinal,              &
             units = 'cmlUnits:countable' )

        call cmlAddParameter( xf=mainXML,              &
             name='MD.LengthTimeStep',&
             value=dt,                &
             units='siestaUnits:fs' )

        call cmlAddParameter( xf=mainXML,                  &
             name='MD.InitialTemperature',&
             value=tempinit,              &
             units='siestaUnits:K' )
        if (idyn/=5) then 
           if(qnch2) then
              call cmlAddParameter( xf    = mainXML,         &
                   name  = 'MD.FireQuench', &
                   value = qnch2 )
           else
              call cmlAddParameter( xf    = mainXML,     &
                   name  = 'MD.Quench', &
                   value = qnch )
           endif
        endif
     endif
  endif

  ! Target Temperature and Pressure
  tt = fdf_get('MD.TargetTemperature',0.0_dp,'K')
  tp = fdf_get('MD.TargetPressure',0.0_dp,'Ry/Bohr**3')
  !
  ! Used for now for the call of the PR md routine if quenching
  if (idyn == 3 .AND. iquench > 0) call set_target_stress()


  ! Mass of Nose variable
  mn = fdf_get('MD.NoseMass',mn_default,'Ry*fs**2')

  ! Mass of Parrinello-Rahman variables
  mpr = fdf_get('MD.ParrinelloRahmanMass',mpr_default,'Ry*fs**2')

  if (idyn==2 .or. idyn==4) then
     if (ionode) then
        write(6,6) 'redata: Nose mass                        = ',mn,'  Ry/fs**2'
     endif
     if (cml_p) then
        call cmlAddParameter( xf    = mainXML,              &
             name  = 'MD.NoseMass',        &
             value = mn,                   &
             units = 'siestaUnits:Ry_fs__2')
     endif
  endif

  if (idyn==3 .or. idyn==4) then
     if (ionode) then
        write(6,6) 'redata: Parrinello-Rahman mass           = ',mpr,'  Ry/fs**2'
     endif
     if (cml_p) then
        call cmlAddParameter( xf    = mainXML,                   &
             name  = 'MD.ParrinelloRahmanMass', &
             value = mn,                        &
             units = 'siestaUnits:Ry_fs__2' )
     endif
  endif

  ! Annealing option
  ianneal = 0
  annop = fdf_get( 'MD.AnnealOption','TemperatureAndPressure' )

  if (idyn .eq. 5) then
     if (leqi(annop,'Temperature')) then
        ianneal = 1
     else if (leqi(annop,'Pressure')) then
        ianneal = 2
     else if (leqi(annop,'TemperatureAndPressure')) then
        ianneal = 3
     else
        call die( 'redata: ERROR: With annealing MD, you must '//&
             'choose an appropriate value for MD.AnnealOption' )
     endif

     if (ionode) then
        select case (ianneal)
        case(1)
           write(6,2) 'redata: Annealing Option                 = Temperature'
           if (cml_p) then
              call cmlAddParameter( xf    = mainXML,           &
                   name  = 'MD.AnnealOption', &
                   value = 'Temperature' )
           endif

        case(2)
           write(6,2) 'redata: Annealing Option                 = Pressure'
           if (cml_p) then
              call cmlAddParameter( xf    = mainXML,           &
                   name  = 'MD.AnnealOption', &
                   value = 'Pressure')
           endif

        case(3)
           write(6,2) 'redata: Annealing Option                 = '//&
                'Temperature and Pressure'
           if (cml_p) then
              call cmlAddParameter( xf    = mainXML,                &
                   name  = 'MD.AnnealOption',      &
                   value = 'TemperatureAndPressure')
           endif
        end select
     endif
  endif


  if (idyn==2 .or. idyn==4 .or. (idyn==5 .and. (ianneal ==1 .or. ianneal==3))) then
     if (ionode) then
        write(6,6) 'redata: Target Temperature               = ',tt,'  Kelvin'
     endif

     if (cml_p) then
        call cmlAddParameter( xf    = mainXML,                &
             name  = 'MD.TargetTemperature', &
             value = tt,                     &
             units = 'siestaUnits:K' )
     endif
  endif

  if (idyn==3 .or. idyn==4 .or. (idyn==5 .and. (ianneal==2 .or. ianneal==3))) then
     if (ionode) then
        write(6,6) 'redata: Target Pressure                  = ', tp, '  Ry/Bohr**3'
     endif

     if (cml_p) then
        call cmlAddParameter( xf=mainXML,                     &
             name= 'MD.TargetPressure',      &
             value= tp,                      &
             units= 'siestaUnits:Ry_Bohr__3' )
     endif
  endif

  ! Relaxation Time for Annealing
  taurelax = fdf_get( 'MD.TauRelax',taurelax_default,'fs' )
  if (idyn==5) then
     if (ionode) then
        write(6,6) 'redata: Annealing Relaxation Time        = ', taurelax,'  fs'
     endif
     if (cml_p) then
        call cmlAddParameter( xf    = mainXML,          &
             name  = 'MD.TauRelax',    &
             value = taurelax,         &
             units = 'siestaUnits:fs' )
     endif
  endif

  ! Estimated Bulk modulus (for Pressure annealing)
  bulkm = fdf_get( 'MD.BulkModulus',bulkm_default,'Ry/Bohr**3' )
  if (ionode) then
     if (idyn==5 .and. (ianneal==2 .or. ianneal==3)) then
        write(6,6) 'redata: Approx. Bulk Modulus             = ', bulkm,&
             '  Ry/Bohr**3'
     endif
  endif

  if (cml_p) then
     call cmlAddParameter( xf    = mainXML,                &
          name  = 'MD.BulkModulus',       &
          value = bulkm,                  &
          units = 'siestaUnits:Ry_Bohr__3')
  endif

  ! Atomic displacement for force constant calculation
  dx = fdf_get('MD.FCDispl',dx_default,'Bohr')

  ! First and last atoms to displace for calculation of force constants
  ia1 = fdf_get('MD.FCfirst',1)
  ia2 = fdf_get('MD.FClast',na)

  ! Check that last atom doesn't exceed total number
  if (idyn.eq.6.and.ia2.gt.na) then
     call die( 'redata: ERROR:'//&
          'Last atom index for FC calculation is > number of atoms.')
  endif

  if (idyn==6) then
     if (ionode) then
        write(6,6) 'redata: Atomic displ for force constants  = ',dx,'  Bohr'
        write(6,4) 'redata: First atom to move               = ',ia1
        write(6,4) 'redata: Last atom to move                = ',ia2
     endif
     if (cml_p) then
        call cmlAddParameter( xf    = mainXML,           &
             name  = 'MD.FCDispl',      &
             value = dx,                &
             units = 'siestaUnits:Bohr' )

        call cmlAddParameter( xf= mainXML,        &
             name= 'MD.FCFirst', &
             value= ia1,         &
             units= 'cmlUnits:countable' )

        call cmlAddParameter( xf= mainXML,       &
             name= 'MD.FCLast', &
             value= ia2,        &
             units= 'cmlUnits:countable' )

     endif
  endif

  ! Variable cell shape? Depending on input and type of dynamics
  varcel = varcel .or. (idyn==3) .or. (idyn==4)             &
       .or. (idyn==5 .and. ianneal==1)           &
       .and. (idyn/=1) .and. (idyn/=2)           &
       .and. (idyn/=6) .and. (idyn/=7)           &
       .and. (.not. (idyn==5 .and. ianneal/=1) )

  want_spatial_decomposition = fdf_get('UseSpatialDecomposition', .false.)
  want_domain_decomposition = fdf_get('UseDomainDecomposition', .false.)
#ifndef ON_DOMAIN_DECOMP
#ifdef MPI
  if (want_domain_decomposition) then
     call die("Need to compile with ON_DOMAIN_DECOMP support")
  endif
#endif
#endif

  ! Harris Forces?. Then DM.UseSaveDM should be false (use always
  ! Harris density in the first SCF step of each MD step), and
  ! MaxSCFIter should be  2, in the second one the Harris 
  ! forces are computed. Also, should not exit if SCF did 
  ! not converge.

  harrisfun = fdf_get('Harris_functional',.false.)

  if (harrisfun) then
     usesavedm = .false.
     nscf      = 1  ! Note change from tradition, since siesta_forces        
     ! now explicitly separates the "compute_forces"        
     ! phase from the rest of the scf cycle.          
     mix       = .false.
     SCFMustConverge = .false.
  endif

  if (ionode) then
     write(6,'(2a)') 'redata: ', repeat('*', 71)
  endif

  if (cml_p) then
     call cmlEndParameterList(mainXML)
  endif

  ! Warn the user about deprecated symbols...
  call deprecated('Optim.Broyden.History.Steps')
  call deprecated('Optim.Broyden.Cycle.On.Maxit')
  call deprecated('Optim.Broyden.Variable.Weight')
  call deprecated('Optim.Broyden.Debug')
  call deprecated('Optim.Broyden.Initial.Inverse.Jacobian')


  ! Find some switches 
  writek                = fdf_get( 'WriteKpoints', outlng )
  writef                = fdf_get( 'WriteForces', outlng )
  writedm               = fdf_get( 'WriteDM', .true. )
  writedm_cdf           = fdf_get('WriteDM.NetCDF', .false. )
  writedm_cdf_history   = fdf_get('WriteDM.History.NetCDF', .false. )
  writedmhs_cdf         = fdf_get('WriteDMHS.NetCDF', .false. )
  writedmhs_cdf_history = fdf_get('WriteDMHS.History.NetCDF', .false.)
  read_charge_cdf       = fdf_get('SCF.Read.Charge.NetCDF' , .false. )
  read_deformation_charge_cdf = &
       fdf_get('SCF.Read.Deformation.Charge.NetCDF', .false. )

  if (read_charge_cdf .or. read_deformation_charge_cdf) then
     mix = .false.
  endif

  save_initial_charge_density = fdf_get(    &
       'SaveInitialChargeDensity' , .false.)

  analyze_charge_density_only = fdf_get(    &
       'AnalyzeChargeDensityOnly' , .false.)

  new_diagk              = fdf_get( 'UseNewDiagk', .false. )
  writb                  = fdf_get( 'WriteBands', outlng )
  writbk                 = fdf_get( 'WriteKbands', outlng )
  writeig                = fdf_get('WriteEigenvalues', outlng )
  writec                 = fdf_get( 'WriteCoorStep', outlng )
  writmd                 = fdf_get( 'WriteMDhistory', .false. )
  writpx                 = fdf_get( 'WriteMDXmol', .not. writec )
  default                = fdf_get( 'UseSaveData', .false. )
  savehs                 = fdf_get( 'SaveHS', .false. )
  fixauxcell             = fdf_get( 'FixAuxiliaryCell', .false. )
  naiveauxcell           = fdf_get( 'NaiveAuxiliaryCell', .false. )
  initdmaux              = fdf_get( 'ReInitialiseDM', .TRUE. )
  allow_dm_reuse         = fdf_get( 'DM.AllowReuse', .TRUE. )
  allow_dm_extrapolation = fdf_get( 'DM.AllowExtrapolation', .TRUE. )
  dm_normalization_tol   = fdf_get( 'DM.NormalizationTolerance',1.0d-5)
  normalize_dm_during_scf= fdf_get( 'DM.NormalizeDuringSCF',.true.)
  muldeb                 = fdf_get( 'MullikenInSCF'   , .false.)
  spndeb                 = fdf_get( 'SpinInSCF'   , (nspin>1) )
  rijmin                 = fdf_get( 'WarningMinimumAtomicDistance', &
       1.0_dp, 'Bohr' )
  bornz                  = fdf_get( 'BornCharge'   , .false. )
  if (idyn.ne.6) then
     bornz = .false.
  endif
  change_kgrid_in_md           = fdf_get('ChangeKgridInMD', .false.)
  ParallelOverK                = fdf_get('Diag.ParallelOverK', .false.)
  ! If non-collinear spin, it *MUST* be false.
  if ( nspin > 2 ) ParallelOverK = .false.
  RelaxCellOnly                = fdf_get('MD.RelaxCellOnly', .false.)
  RemoveIntraMolecularPressure = fdf_get( &
       'MD.RemoveIntraMolecularPressure', .false.)
  !
  !   COOP-related flags
  !
  write_coop = fdf_get('COOP.Write', .false.)
  !
  saverho  = fdf_get( 'SaveRho', dumpcharge)
  savedrho = fdf_get( 'SaveDeltaRho',       .false. )
  saverhoxc= fdf_get( 'SaveRhoXC', .false.)
  savevh   = fdf_get( 'SaveElectrostaticPotential', .false. )
  savevna  = fdf_get( 'SaveNeutralAtomPotential', .false. )
  savevt   = fdf_get( 'SaveTotalPotential', .false. )
  savepsch = fdf_get( 'SaveIonicCharge', .false. )
  savebader= fdf_get( 'SaveBaderCharge',  .false.)
  savetoch = fdf_get( 'SaveTotalCharge', savebader )

  !
  !   Siesta2Wannier90 -related flags
  !
  w90_write_mmn = fdf_get( 'Siesta2Wannier90.WriteMmn',   .false. )
  w90_write_unk = fdf_get( 'Siesta2Wannier90.WriteUnk',   .false. )
  w90_write_amn = fdf_get( 'Siesta2Wannier90.WriteAmn',   .false. )
  w90_write_eig = fdf_get( 'Siesta2Wannier90.WriteEig',   .false. )

  w90_processing = ( w90_write_mmn .or. w90_write_unk .or. &
       w90_write_amn .or. w90_write_eig )

  hasnobup   = fdf_defined( 'Siesta2Wannier90.NumberOfBandsUp'   )
  hasnobdown = fdf_defined( 'Siesta2Wannier90.NumberOfBandsDown' )
  hasnob     = fdf_defined( 'Siesta2Wannier90.NumberOfBands'     )

  nobup      = fdf_get( 'Siesta2Wannier90.NumberOfBandsUp',   0)
  nobdown    = fdf_get( 'Siesta2Wannier90.NumberOfBandsDown', 0)
  nob        = fdf_get( 'Siesta2Wannier90.NumberOfBands',     0)

  RETURN
  !-------------------------------------------------------------------- END
1 format(a,4x,l1)
2 format(a)
4 format(a,i8)
5 format(a,i5,a)
6 format(a,f10.4,a)
7 format(a,f12.6,a)
8 format(a,f14.12)
9 format(a,f10.4)

CONTAINS
  subroutine deprecated( str )

    implicit none

    character(len=*), intent(in) :: str

    if (ionode) then
       if (fdf_defined(trim(str))) then
          write(6,'(a)') '**Warning: FDF symbol '//trim(str)//&
               ' is deprecated. See the manual.'
       endif
    endif

  end subroutine deprecated

end subroutine read_options


