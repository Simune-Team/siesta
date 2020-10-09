# ChangeLog for Siesta

## 4.0.3 (2020- ) FUTURE Bug fix release

### Changes

* Removed some minor memory leaks in mesh-subs

* Update and add .md extension to main README

* Update Docs/REPORTING_BUGS

* Add ChangeLog.md since 4.0

* Fixed cell transpose for socket communication, fixes lp:1835196

* Make matel registry pool re-sizeable

* Fixed grid output when using cell-sampling, lp:1799991

* Fixed vibra utility precision issue, fixed lp:1816719

* Increase flexibility in the handling of pseudopotentials

* Fix some issues with polarization orbitals. Bessel orbs improvements.

* Fix C10 computation in molecularmechanics


## 4.0.2 (2018-07-19) Bug fix release

### Backward compatibility issues

* This release increases the size of the internal tables for
  two-center integrals used in some matrix element calculations. This
  means that calculations are slightly more heavy, but the accuracy is
  also superior. One can regain the *old* less accurate behaviour by
  setting Compat.Matel.NRTAB to true in the fdf input file (this is
  ONLY recommended for testing purposes).

### Changes

* Enabled ion.nc files for ghost atoms (#1738425)

* Enabled ghost atoms to read ion.nc files (#1736455)

* Forced libfdf to "die" when too long input strings are passed (#1728281)

* Monkhorst-Pack grids not properly shifted to [0;1[ when user specified large displacements (#1721479)

* Updated README content in Util directory (#1712319)

* Fixed building all utilities (#1712317)

* Fixed EIG file format for non-collinear spin (#1708634)

* Fixed possible segfault when using too large unit-cells (#1704370)

* Fixed ghost orbital energies (#1695130)

* Removed vpsa2bin and vpsb2asc codes (part of ATOM)

* Added new units in fdf (Hartree [Ha] and milli-Hartree [mHa])

* Added more tests

* Fixed a missing close when using Siesta in "master" mode

* Fixed several memory accounting errors and a couple of missing close
  statements. Also removed all deprecated PASTE calls.

* Added compiler version information to the version.F90 file

* Fixed writing the Bessel ghost atom to ion.xml/nc files

* Extended information in PDOS xml files. Now the atomic number
  as well as whether the orbital is a polarization shell is present

* Increased the default size of two-center integrals tabulation arrays from 128 to 1024.
  This is a change that results in more accurate values (of basically everything). Set
  'Compat.Matel.NRTAB true' to use previous value

* Added spin-monitoring in the SCF: for spin-polarized and non-collinear calculations
  the total spin-moment is written

* Fixed non-collinear Mulliken populations in parallel

* Updated gnubands to the new code that has been present since the 4.0 release.
  The new updated gnubands provides more functionalities

* Extended precision output in EIG, KP and PDOS files

* Fixed possible non-optimal DM initialization when unit-cell folding of orbital
  connections is performed.

* Bugfix for writing out too many wavefunction coefficients,
  see https://www.mail-archive.com/siesta-l@uam.es/msg10291.html

## 4.0.1 (2017-07-04) Bug fix release

### Changes

* Better standard compliance in code structure

* Fix bug related to SlabDipoleCorrection which couldn't be turned off (lp:1630827)

* Fix non-collinear bandstructure calculations (lp:1636100)

* Fix for 'nodes' basis generation option (lp:1625725)

* Fixes VCA mixing of pseudos (lp:1633039)

* Fixes integer energy specifications in ProjectedDensityOfStates block (lp:1657584)

* Added print-outs when GridCellSampling is used

* TBtrans_rep changed the written DOS units to 1/eV (they were in 1/Ry)

* Fix for Bader charge analysis (lp:1656273)

* Fix memory problem when memory usage close to limit (lp:1665294)

* Forced Diag.ParallelOverK to false if non-collinear spin configuration (lp:1666428)

* Updated Eig2DOS to be more like gnubands (options are the same)

* Added Geometry.Must.Converge flag to ensure the relaxation has converged

* Enabled internal walltime check to forcefully stop SIESTA after a
  certain limit.

* Updates and fixes for Util/STM/ol-STM
   - Fix wrong fftw call
   - Enabled direct reading of WFSX files

* Added new interpolation option to Util/Macroave


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