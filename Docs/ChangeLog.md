# ChangeLog for Siesta

## Master version updates

### Backward-compatibility issues:

* The labels in the Mulliken analysis CML blocks have been changed to use "population" instead of "charge".

### Changes

* Some fixes for library operation (avoid stopping when 'onlyS' is set; re-opening of unit 6)

* Cleanup of vibra/fcbuild code: use dynamic arrays, and LAPACK solver by default

* Extended Hirshfeld and Voronoi partition analysis to spin. Update of output, including new CML blocks.

*(Initial list of updates taken on April 8, 2020 from the output of `git see --first-parent rel-4.1..`,
filtering the direct merges from 4.1 series, but kept mostly in raw form,
with commit hash and date, for further editing, if needed.  Even
though the 4.1 and `master` progress was in parallel, one needs to
look at the 4.1 ChangeLog to get a complete picture of the changes in
`master` -- for example: GPU support)*

* 147db2b3 2020-03-05 Nick Papior   enh: performance increase for tall+skinny matrices in TS

* 586b553e 2019-10-21 Nick Papior   New 'random' experimental option to initialize the DM

* e4c44a01 2019-09-16 Jose M Soler   Implementation of a band unfolding utility

* 55b5f907 2019-03-06 Alberto Garcia   Benchmark support: dummy solver and scf-step granularity

* 2de79ac6 2018-12-07 Alberto Garcia   Compliance of code to F2003 standard plus further cleaning

* 8ed60c2e 2018-11-23 Nick Papior   Added class_*3D to allow for NC/SOC transiesta implementation

* 7f5ade48 2018-07-09 Nick Papior   Output gradient of S to nc file

* 17bee605 2018-07-09 Nick Papior   Enabled pjnl_j in ion.nc files

* 707d71f6 2018-06-18 Nick Papior   Enabled kgrid.MonkhorstPack as list input (only if block is not present)

* 2abadff2 2018-05-27 Nick Papior   Restructured k-point sampling, allowing user-supplied k-points

* f17ad461 2018-04-30 Ramon Cuadrado   Merge the 'offsite spin-orbit' (full) implementation by Ramon Cuadrado.

* 4b52c778 2018-01-03 Nick Papior   Initialize some components in the species_info derived type

* 6150be62 2017-11-20 Alberto Garcia   Remove old data structures and routines related to 'atom'

* c142809d 2017-10-10 Nick Papior   Update TD-DFT: k-points, H extrapolation, more parallelization, code readability

* 01356192 2017-08-22 Nick Papior   Interface to the CheSS linear-scaling library

* 8d69f1cf 2017-02-02 Nick Papior   Update of the siesta_forces logic wrt TDDFT

* e2f7686c 2016-09-28 Rafi Ullah   Implementation of Real-Time TDDFT

* cad6b290 2016-08-29 Nick Papior   Added OpenMP nesting information

## 4.1 (2020-   )   FUTURE Feature release

### Changes to defaults

* TranSiesta/TBtrans eta values are now default to 1meV

* TranSiesta equilibrium contours now default to using the "right" scheme

* Maximum l for KB projectors is now automatically set to the highest
  l in the PS file. This may result in slight changes.

### Changes

* Siesta is now developed on the GitLab platform: www.gitlab.com/siesta-project/siesta
  A number of Siesta-related packages are developed here: www.gitlab.com/siesta-project

* Document the setting of 'neigwanted' and print them if the diag solver allows it.

* Fix computation of NC/SOC occupations when the (optional) number of
  eigenstates handled ('neigwanted') is less than the number of orbitals.

* Fix reading of wave-functions in Util/COOP/fat.f90

* Added Obj/ARCH-EXPERIMENTAL for suggested more modular building scheme

* Added interface code to use GPUs with the ELPA library.

* Added Docs/ChangeLog.md file

* Removed a small memory leak (from siesta_init)

* Enabled parallel-over-k for NC/SOC calculations

* Fixed bug for writing eigenvalues when NumberOfOrbitals was used (together with NC/SOC)

* Fixed printing Voronoi and Hirshfeld charges if user requested LDOS calculations

* Fixed cell transpose when using socket calculations (only important for skewed cells)

* Added citation information to output (end of run)

* Allowed Mesh.Sizes as list input so users can specify their own Mesh size

* Fixed *.nc writes, fixes lp:1810279.

* Fixes polarization issue with Bessel orbitals

* Made the optional BSC_CELLXC code accessible at run time, instead of
  through pre-processing at compile time.

* Lots of updates for documentation

* Precision problem in the vibra utility is fixed, lp:1816719

* Added wavefunction tools for spin-orbit calculations.
  Now denchar, COOP/COHP, fatbands, and stm utilities can process spin-orbit output.

* Removed *ALL* OMP collapse statements, Intel 2019 is buggy.

### TranSiesta/TBtrans changes:

* Much more efficient dq implementation for fixing charge fluctuations

* Fixed TBtrans calculations with spin-flags (i.e. TB-only)

* Allow E-field for 1-electrode calculations, this allows 1-electrode
   capacitor setups

* Fixed bug in voltage potential which was severe for capacitor like setups.
   It also changes regular TS runs.

* Added TSFA[C] files which show forces on atoms in the device region

* Added TS-only energies, that is energies calculated from DM/H are
   now only using the updated elements of the respective regions.

* Fixed a bug in TBtrans when piping input

* Much faster Bloch expansions in transiesta/tbtrans calculations
   Using tiling will now result in extremely fast self-energy calculations

* Allowed complex contours in tbtrans calculations via negative eta values
   and user defined energy contours, the precision of outputted contours
   is also increased to the maximum field width.

* Changed Eta defaults in electrodes

* Enabled usage of external truly infinite electrodes (real-space self-energies)

* Fixed contour file output which is now usable for various post-processing
   utilities.

* Fixed tbtrans command line arguments, lp:1829974

* Added delta-Ef to the electrode block to specify off-set in electrodes.
   Mainly useful for semi-conducting electrodes.

* Added separate energies for NEGF calculation, now the TranSiesta energies
   are more divided and should be more comperable since they are calculated
   on the updated sub-set.


## 4.1-b4 (    )  Bugfix beta release

### Changes to default operational parameters

* MeshCutoff has been increased to 300 Ry (from 100)
* MaxSCFIterations has been increased to 1000 (from 50)
* SCFMustConverge is now default true (from false)

### Changes

* Added developer documentation found in Docs/developer
    Both ford (preferred) and Doxygen may be used

* Generally increased precision in many output files

* Lots of fixes and updates for the Lua/flook interaction

* Auxiliary supercell handling when reading DM matrices:
    Siesta can now read and convert nearly *any* DM matrix and make it
    match the used sparse pattern.

* Fixed minor inconsistencies when handling Bessel basis

* Updated all diagonalization routines
    - ELPA and MRRR for k-point sampling.
    - Less memory usage

* Fixed bug on reading *.ion* files (lp:1751723)

* Updated internal integration table sizes (slightly increased precision)

* PDOS files now also contain the fermi-level such that tools may easily
    align the energy to 0.

* Added more digits to dDmax which may be relevant when performing
    Spin-Orbit/Non-Collinear calculations.

* Fixed bug related to writing out non-collinear spin eigenvalues,
    and also for spin-orbit. (lp:1708634)

* Fixed parallel PDOS calculations of non-colinear and spin-orbit.
    (lp:1718162)

* Added calculated charges to the Lua interface (check the charges
    while running).

* Fixed lots of compilation issues related to the utilities
    (lp:1712317, lp:1712319, lp:1711850)

* Fix for reading a ghost basis (lp:1736455, lp:1738425)

* Fix when fdf-input lines are too long. Instead of discarding the
    remaining line, fdf now "dies" to inform users of possible erroneous
    input. (lp:1728281)

* Fixed Monkhorst-Pack displacements when the displacement was larger
    than 1 (lp:1721479)

* Fix for possible heap allocated arrays (Intel compilers) (lp:1704370)

* Ensured many files to be closed properly at the end of the runs.

* Added basic compiler information to the siesta/transiesta/tbtrans
    header (compiler output)

* Performing SOC calculations does not not require all species
    to have SOC contributions.

### TranSiesta / TBtrans changes:

* Disk-space reduction when mixing non-periodic and periodic electrodes

* Now tiling is also enabled for Bloch expansions. This is actually faster
   than repetitions, so users should prefer tiling the electrodes

* TranSiesta is now intrinsic to the Siesta executable. An
   electrode should now be calculated using 'siesta --electrode'
   The TranSiesta executable still exists but is nothing but 'siesta --electrode'

* Many bug-fixes related to pivoting tables; this should only
   change the effective BTD matrices, and has no relevance to the
   accuracy of the calculations

* Huge performance increase in TBtrans in many different parts of the code

* Bug-fix for out-of-core calculations for spin-polarized TBT.Spin 2 calculations

* Fixed the default 2-terminal Poisson ramp. The ramp is now
   defaulted to be positioned in the central region.
     TS.Poisson ramp-central

* Small memory reduction by de-allocating unused siesta memory when
   entering transiesta.
   
* Fixed the box Poisson for N-electrode calculations when using
   skewed electrodes. Thanks to Azar Ostovan and Thomas Frederiksen.

* Fixed tbtrans setup for bias-window-only calculations. Now the contours
   are correctly interpreted.

* Fixed tbtrans AVCEIG output.

* Change TBtrans DOS output such that there is no normalization

* Enabled tbtrans 1-orbital calculations in the BTD matrices.

* Fixed sign-convention changes in orbital-currents. Now they are
   checked and works together with sisl (>0.9).

* Allowed external GF files for the self-energies. This is mainly beneficial
   for TBtrans as we can add external electrodes *anywhere* in the device.
   Say Buttiker-probes.

* Bugfix when the left electrode was set to -|V|/2 (the default |V|/2 is
   unaffected).

* Added much more output to the TBT*.nc files; electrode information is now
   complete, and also the BTD matrices are written.

* Enabled tbtrans -fdf TBT.Analyze which runs all pivoting schemes, this
   may be very beneficial to run with tbtrans before performing calculations.
   Choosing the correct pivoting scheme can be really important!

* Enabled output file on tbtrans command line:
     tbtrans --out TBT.out RUN.fdf
   is (more or less) equivalent to:   
     tbtrans RUN.fdf > TBT.out

* Made Fermi charge correction more aggressive for faster convergence.

* TBtrans can now calculate DM, COOP and COHP curves. They are calculated
   in the supercell picture and can thus be analyzed cross-boundary as well.
   They are calculated both from the Green function and the spectral function.
   The coming >0.9.3 release of sisl will enable this analysis.

* Fixed TBtrans DOS (Green) calculations when performing k-point calculations. There
   can be small differences when comparing Green function DOS between this version
   and prior versions. The bug is only present when time-reversal-symmetry is applied.

## 4.1-b3  () Bugfix beta release


* Manual greatly overhauled and updated in various parts

* Fixed DOS and PDOS for non-colinear and spin-orbit

* Fixed bug when printing initial spin-configuration

* Enabled restarting calculations with different spin-configurations,
    i.e. one may go from an unpolarized calculation to a polarized, or
    from a polarized to an unpolarized (also non-colinear and spin-orbit).

* Lots of bug-fixes for transiesta and tbtrans

* Bug-fix for spin-orbit coupling normalization

* Fixed minor memory leaks

* Many improvements for Lua enabled runs

* Added installation scripts of
    netcdf/hdf5/zlib/flook

* Fixes to the <>.nc file for high spin configuration >= non-colinear

## 4.1-b2  () Bugfix beta release

* The configure script has been removed.  Its use was discouraged and
  would often yield erroneous arch.make files.  To circumvent any
  confusions it has been obsoleted until further notice.

* Instead of the configure script two default arch.make files now
  exist in the Obj directory (gfortran.make, intel.make) which should
  be guidelines for creating one's own arch.make file.

* Several fixes for bugs reported for the b1 release. See Docs/CHANGES

## 4.1-b1  () Beta release  

Please see the Manual for full details

### Backward-compatibility issues:

* The mixing routines have completely changed, hence the same
    convergence path cannot be expected. This, unfortunately, makes
    comparison difficult with prior versions. However, the final
    converged system should be transferable.

* SIESTA now defaults to mixing the Hamiltonian instead of the
    density matrix. To revert to density-matrix mixing, use
    "SCF.Mix DM". The option to mix after the initial scf step is now
    'on' by default.

* SIESTA now defaults to monitoring convergence for both the
    density matrix AND the Hamiltonian. To revert to only density
    matrix convergence, use: "SCF.Converge.H false"

* A major number of fdf-flags concerning mixing
    parameters have changed to a more consistent naming scheme.
    However, all previous flags are still in effect but the newer
    flags have precedence. The previous flags are the default values
    for the newer flag-names.

* Two additional files are created (H_DMGEN and H_MIXED), which
      contain the Hamiltonian at various stages through the SCF.
      Currently they are intended for developers and may be removed in
      the final 4.1 release.  You may delete these without problems.
      
### New features
    
* LDA+U (Javier Junquera)

    * Full incorporation of the LDA+U implementation in SIESTA
    * Two different LDA+U projectors available
    * Estimate the best U according to: Cococcioni and Gironcoli in PRB, 2005

* Spin-Orbit coupling (Ramon Cuadrado)

    * On-site approximation for spin-orbit-coupling matrix elements.

* MRRR method for diagonalization (Alberto Garcia)

    * This will typically be faster than divide-and-conquer algorithm
    and may be the future default. For Gamma-point calculations.

* ELPA method for diagonalization (Alberto Garcia)

    * This provides better scalability compared to ScaLAPACK for large
    # of processors. For Gamma-point calculations.

* Added interface to the PEXSI library for calculating the density
   matrix, DOS, and LDOS (Alberto Garcia)

    * This library provides massive parallelism and better
    scalability, but should only be used for very large systems.

* SIESTA is now hybrid-parallelised (Nick R. Papior)

    * One may compile Siesta/Transiesta in serial, OpenMP, MPI, or
    MPI+OpenMP mode.

* Re-write of non-equilibrium Green function code (Nick R. Papior)

    * N>=1 terminal support in transiesta
    * Improved convergence
    * Different ways of handling charge-reductions in SCF
    * All electrodes have settings _per electrode_ for full customization
    * Greatly reduced memory usage
    * Skewed axis are enabled for further reduction of complex systems
    * Implemented MUMPS method for inversion
    * Implemented LAPACK for inversion
    * Implemented BTD method for inversion (extremely fast)
    * Fully OpenMP threaded
    * Start directly from transiesta enabled
    * Temperature gradients as non-equilibrium a possibility
    
* Complete rewrite of tbtrans utility (Nick R. Papior)

    * Made tbtrans a stand-alone utility for user-defined tight-binding method
    * EXTREME SCALE version (BTD-only)
      - Memory consumption _only_ dependent on "device" region
    * N>=1 electrode support
    * Region projections for transmission through "eigenstates"
    * Custom change of the self-energies
    * k -> k' transmissions
    * Interpolation of bias Hamiltonians
    * Bond-currents
    * Fully OpenMP threaded and/or MPI parallelized
    * DOS and bulk transmission for electrodes
    * Gf-DOS/spectral-DOS for device region
      

* Mixing routines rewritten (Nick R. Papior)
    
    * New mixing schemes Pulay (Guarenteed Reduction)
    * Custom mixing restart options (full user customizability)

* Added more constraints (Nick R. Papior)

    * Constraints are verbose and many new ways of using constraints exists

* NetCDF4 file format for siesta -> parallel IO (Nick R. Papior)

    * Provides a standard intrinsically grouped file for retaining
    nearly all siesta related information in one file.

* Enabled convergence control of density-, energy density matrices,
   Hamiltonian and energy. 

* LUA scripting in siesta (Nick R. Papior)
    
    * This is an experimental feature

* Gate calculations (Nick R. Papior)

    * Charge and Hartree gate

* Utilities (Nick R. Papior)

    * All make-files are now prepared to enable parallel builds
      - this makes compilation *MUCH* faster. For example:
         make -j4
        will compile using 4 cores.
    * Grimme utility
      - easy creation of FDF block for Grimme parameters
    * SpPivot, pivot sparsity pattern using several different routines
    * TS/** different utilities related to transiesta


## 4.0.0 (2016-06-23) Feature release

This version includes the van der Waals functionals, the new
load-balancing code for real-space grid operations, a Wannier90
interface, a new Orbital-Minimization-Method solver, and other
improvements and bug fixes that have been part of the development
version for some time. 

### Backward compatibility issues

Please take into account the following changes in behavior (more
details in the TECHNICAL NOTES section in Docs/release_notes.4.0)

* The grid functions (charge densities, potentials, etc) were in
single precision by default in the 3.X versions, but are in double
precision by default for post-3.X versions. The 'phi' array that holds
the values of the basis orbitals on the real-space grid is kept in
single precision. Please take this into account if you compare the
results with those of siesta-3.X runs. See the manual in both versions
for more information.

* Changes in the geometry used for the analysis of the electronic
structure, as well as in the handling of the density matrix (DM) and
hamiltonian (H). This will _slightly_ change the output of most
calculations and the detailed results of any post-processing. Keep
this in mind if you need to maintain coherency within a project.

* The default 'dynamics' option has been changed from 'verlet' to 'CG'.
  There should really be a new 'single-point' default which completely
  avoids 'siesta_move'. The old behavior can be recovered by using the
  'compat-pre-v4-dynamics' switch.

* Single-point calculations do not write .STRUCT_NEXT_ITER files, and
  the coordinates in the XV file are the current ones, unmoved.
  Extra output in siesta_options is avoided for this case.

* Electric field and dipole correction for slab calculations

  - Older versions applied an incorrect dipole correction when also
    using an electric field (old behavior may be recovered by forcing
    SlabDipoleCorrection to .false.)

  - Older versions over-estimated the energy contribution from the
    dipole correction by a factor of 2 (old behavior cannot be
    recovered).

### New Features

Please see the relevant section of the manual for more information.

* New SiestaXC library for exchange and correlation, implementing several
  van der Waals functionals and some newer GGA ones.

* New code to improve the load-balancing of the operations in the
  real-space grid when running in parallel.

* A new interface to the Wannier90 code for the generation of
  maximally localized wannier functions.

* New electronic-structure solver implementing the the Orbital-Minimization-Method

* New mixing options, including hamiltonian and charge-density mixing

* Charge-confinement and "filteret" basis-set generation options.

* Improved MPI version of the Siesta-as-subroutine code for executing
  independent calculations.

* New code for 'server' operation via sockets and i-PI interface

* New JobList utility to organize and run multiple jobs.

* Timer with call-tree awareness.

* Performance and usability enhancements in TranSiesta.

* Enhancements to the restart capabilities in molecular-dynamics runs.

* New options for wave-function output. New WFSX format.

* 'Fatbands' analysis.

* Enhancements to the COOP/COHP analyzer

* Hirshfeld and Voronoi charges. Bader analysis output

* Force-convergence diagnostics.

* New HSX file format for Hamiltonian/Overlap files

* Calculation of the vacuum level for non-bulk systems.

* Updates to other analysis tools in Util/

* Replacement of license-encumbered routines by new ones.

* Enhancements to the manual and build system.

### Other changes

* A number of bugs have been fixed, and there have been
numerous cosmetic changes to the output and the code itself. See Docs/CHANGES
for a full list.

## 3.X 2.X 1.X 0.X  Older releases

See the files in `Docs/older_release_notes`

<!--
Local Variables:
mode: fundamental
fill-column: 70
indent-tabs-mode: nil
coding: utf-8
End:
-->