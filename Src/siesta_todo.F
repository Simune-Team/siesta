      Module SIESTA_TODO

      use precision, only: dp
      !use parallel, only: Node, Nodes, IOnode, ParallelOverK
      !use atomlist, only: indxuo, elem, indxua, iza, lastkb, lasto
      !use atomlist, only: qa, amass, qtot, no_u, rmaxo, rmaxkb
      !use atomlist, only: zvaltot
      !use atomlist, only: rmaxv, iphorb, iphkb, rco, rckb
      !use atomlist, only: datm, no_s, iaorb
      !use atomlist, only: superx, superc, initatomlists

      !use band, only: setup_bands, bands
      !use fdf, only : fdf_block, fdf_convfac, fdf_boolean
      !use fdf, only : fdf_physical
      !use m_fdf_global, only: fdf_global_get
      !use sys, only: die, bye
      !use parse
      !use periodic_table, only : symbol
      !use xcmod, only: setXC
      !use molecularmechanics, only : inittwobody, twobody
      !use metaforce, only: lMetaForce, initmeta, meta
#ifdef MPI
      use mpi_siesta, only: mpi_comm_world
      use m_mpi_utils, only : globalize_sum, globalize_max
#endif
      !use m_mpi_utils, only : broadcast
      !use alloc, only: re_alloc, alloc_report
      !use phonon, only: phonon_num_disps, phonon_setup
      !use phonon, only: phonon_write_forces, phonon_restore_coords
      !use phonon, only: phonon_set_coords

      !use densematrix

      !use m_ordern,     only: ordern
      !use m_hsparse,    only: hsparse

      !use parallelsubs, only: getnodeorbs
      !use writewave,    only: initwave, wwave
      !use iopipes,      only: forcestopipe, coordsfrompipe
      !use m_iostruct,   only: write_struct, read_struct
      !use m_nlefsm,     only: nlefsm
      !use m_overfsm,    only: overfsm
      !use m_check_supercell, only: check_sc_factors

      !use m_broyden_mixing

      !use siesta_cml

      !use files,       only : label_length

      !use m_timestamp, only : timestamp
      !use m_wallclock, only : wallclock
      !use units
      !use zmatrix,      only: lUseZmatrix, iofaZmat
      !use zmatrix,      only: CartesianForce_to_ZmatForce
      !use zmatrix,      only: write_Zmatrix
      !use m_broyden_optim,      only:broyden_optimizer

      !use m_dynamics, only: nose, verlet2, npr, anneal, pr
      !use md_out,     only: md_v_format
#ifdef CDF
      use md_out,     only: md_netcdf
#endif
      !use Kpoint_grid

      implicit none
      integer:: fincoor ! Final step in geommetry iteration
      integer:: i       ! General-purpose index
      integer:: ia      ! Atom index
      integer:: iadispl ! Displaced atom index
      integer:: ifa     ! Calculate forces? 0=>no, 1=>yes
      integer:: ihmat   ! Calculate hamiltonian? 0=>no, 1=>yes
      integer:: ihuge=1073741823 ! Huge integer
      integer:: iiscf
      integer:: ik      ! k-vector index
      integer:: ind     ! Index of sparse matrix elements
      integer:: inicoor ! Initial step in geommetry iteration
      integer:: io      ! Basis-orbital index
      integer:: iord    ! Order of density matrix extrapolation between geom. steps
      integer:: is      ! Species index
      integer:: isel    ! Auxiliary to call neighb: avoid double count of atom pairs?
      integer:: ispin   ! Spin index
      integer:: istp    ! Geommetry iteration step starting in istp=1
      integer:: istpsave! Temporary to save and restore istp when reading DM from file
      integer:: istr    ! Calculate stress? 0=>no, 1=>yes
      integer:: iu      ! I/O file unit
      integer:: iunit   ! Physical-units option: 1=>(eV,Ang), 2=>(Ry,Bohr)
      integer:: iv      ! Vector index. Not used
      integer:: ix      ! Cartesian coordinate index
      integer:: ixdispl ! Displaced cartesian coordinate
      integer:: j       ! General-purpose index
      integer:: ja      ! Atom index
      integer:: jo      ! Orbital index
      integer:: jx      ! Cartesian coordinate index
      integer:: maxwk   ! Number of k-points for wavefunction printout
      integer:: nbcell  ! Auxiliary to call shaper routine
      integer:: neigmin ! Min. number of eigenstates (per k point)
      integer:: nh      ! Max. nonzero elements in a row of the hamiltonian matrix
      integer:: nhist   ! Number of histogram intervals in projected density of states
      integer:: nnia    ! Number of neighbors of a given atom
      integer:: ntm(3)  ! Number of integration mesh divisions of each cell vector
      integer:: ntcon   ! Total number of position constraints imposed

      ! maxnh is (initially) the minimum size of all the sparse arrays.
      ! It must be 1 (not 0) since we frequently use the F77 idiom of
      ! passing the first element of the array instead of all of it.
      integer :: maxnh=1   ! Max. number of nonzero H matrix elements
      integer :: maxna=200 ! Max. number of neighbor atoms
      integer :: no_l=1    ! Number of basis orbitals stored in one node
      integer :: nnamax=1  ! Max. required number of neighbor atoms
      integer :: nauxpul=1 ! Max. required size of auxiliary sparse Pulay matrices
      integer :: ns=1      ! AG** do not initialize Number of chemical species
      integer :: nspin=1   ! Number of spin components
      integer :: nxij=1    ! Number of neighbors of a given basis orbital

      integer :: nscold(3) = 0 ! Previous diagonal elements of supercell vectors
      integer :: nsc(3) = 1    ! Diagonal elements of supercell vectors
      integer :: mscell(3,3) = 
     .                    reshape ((/1,0,0,0,1,0,0,0,1/), (/3,3/)) ! Matrix of auxiliary supercell

      integer, allocatable :: jna(:) ! Index of neighbor atoms

      real(dp):: bcell(3,3) ! Reciprocal unit vectors
      real(dp):: cfmax
      real(dp):: cftem
      real(dp):: const      ! Temporary factor (used in calculation of Harris forces)
      real(dp):: cstress(3,3) ! Constrained stress tensor
      real(dp):: dipol(3)   ! Electric dipole vector
      real(dp):: dDmax      ! Max. change in density matrix elements
      real(dp):: dEmax
      real(dp):: DEharr
      real(dp):: DEna       ! Neutral-atom energy term, calculated in dnaefs
      real(dp):: DUext      ! Interaction energy with external electric field,
                            ! calculated in dhscf
      real(dp):: DUscf      ! Electrostatic energy of (rhoscf-rhoatm), calc. in dhscf
      real(dp):: Dxc        ! Integral((epsxc-Vxc)*Rho), calculated in dhscf (not used)
      real(dp):: sigma      ! Energy width used to convolute partial density of states
      real(dp):: e1=1       ! Lower bound of calculated eigenvalues
      real(dp):: e2=-1      ! Upper bound of calculated eigenvalues
      real(dp):: Ecorrec    ! Energy term eta*DeltaQ, calculated in ordern
      real(dp):: ef         ! Fermi energy
      real(dp):: Eharrs     ! Harris-functional total energy
      real(dp):: Eharrs1    ! Same as Eharrs, but preserved in grid cell sampling
      real(dp):: Eions      ! Self-energy of isolated ions
      real(dp):: Ekin       ! Kinetic energy of electrons, calculated in kinefsm
      real(dp):: Ekinion    ! Kinetic energy of ions
      real(dp):: Elast      ! Total energy in the previous SCF iteration
      real(dp):: Emad       ! Madelung energy term, calculated in madelung
      real(dp):: Ena        ! Neutral-atom term in the total energy, calculated in naefs
      real(dp):: Enaatm     ! Integral of Vna * rhoatm, calculated in dhscf
      real(dp):: Enascf     ! Integral of Vna * rhoscf, calculated in dhscf
      real(dp):: Enl        ! Non-local pseudopot. energy term, calculated in nlefsm
      real(dp):: Emeta      ! Metadynamics energy contribution calculated in meta
      real(dp):: Entrop     ! Temporary to call diagon
      real(dp):: Entropy    ! Entropy due to electron state occupations, calc. in diagon
      real(dp):: Etot       ! Total electronic energy
      real(dp):: Etot_output ! Temporary to select among different total-energy opts
      real(dp):: Exc         ! Exchange-correlation energy, calculated in dhscf
      real(dp):: E0         ! Non-SCF part of total energy
      real(dp):: Emm        ! Classical two-body term, calculated in twobody
      real(dp):: factor     ! Temporary variable
      real(dp):: fmax       ! Maximum atomic force
      real(dp):: fmean      ! Mean atomic force (total force / number of atoms)
      real(dp):: FreeE      ! Free energy
      real(dp):: fres       ! Root mean squared atomic force
      real(dp):: ftem       ! Temporary force value
      real(dp):: ftot(3)    ! Total force
      real(dp):: g2max      ! Actual planewave cutoff of real-space integration mesh
      real(dp):: getot      ! Total energy including thermostat dynamic variables
      real(dp):: kn         ! Kinetic energy of the Nose' thermostat
      real(dp):: kpr        ! Kinetic energy of the Parrinello-Rahman variables
      real(dp):: Pint       ! Instantaneous pressure
      real(dp):: Pmol       ! "Molecular pressure" (discounting Virial term)
      real(dp):: Psol       ! Solid (periodic system) pressure
      real(dp):: rmax_bonds ! Cutoff length for bond definition
      real(dp):: tiny=1.e-15_dp ! Small tolerance used to compare k points
      real(dp):: tempion=0.0_dp ! Ionic temperature
      real(dp):: qaux       ! Auxiliary to call spnvec (not used later)
      real(dp):: qspin(4)   ! Total spin populations
      real(dp):: qsol       ! Total unnormalized electron charge
      real(dp):: rmax       ! Max. required interatomic distance
      real(dp):: rmaxh      ! Max. distance between non-zero hamiltonian matrix orbitals
      real(dp):: stot       ! Total spin polarization
      real(dp):: stress(3,3) ! Stress tensor = d_E/d_strain
      real(dp):: svec(3)    ! Total spin vector in non-collinear-spin calculations
      real(dp):: tstress(3,3) ! Total stress tensor, includin kinetic components
      real(dp):: Uatm       ! Harris hartree electron energy, calculated in dhscf
      real(dp):: Uscf       ! SCF hartree electron energy, calculated in dhscf
      real(dp):: vcell(3,3) ! Unit cell "velocity" (time derivative)
      real(dp):: virial     ! Virial = Sum(x*f)
      real(dp):: vn         ! Velocity (time derivative) of the Nose thermostat
      real(dp):: volume     ! Unit cell volume
      real(dp):: vpr        ! Velocity (time derivative) of the PR variables
      real(dp):: stressl(3,3) ! Local-mode part of the stress tensor
      real(dp):: veclen     ! Length of a unit-cell vector

      real(dp), external      :: volcel

      real(dp), allocatable   :: auxpul(:,:)! Auxilialy matrices for Pulay mixing
      real(dp), allocatable   :: cfa(:,:)   ! Atomic forces orthogonalized to geometry constraints
      real(dp), allocatable   :: ebk(:,:,:)! Band energies
      real(dp), allocatable   :: efs(:)     ! Fermi energy for each spin (for fixed spin calculations)
      real(dp), allocatable   :: eo(:,:,:) ! Hamiltonian eigenvalues
      real(dp), allocatable   :: fa(:,:)    ! Atomic forces
      real(dp), allocatable   :: fal(:,:)   ! Local-node part of atomic forces
      real(dp), allocatable   :: qo(:,:,:) ! Occupations of eigenstates
      real(dp), allocatable   :: qs(:)      ! Spin charges
      real(dp), allocatable   :: r2ij(:)    ! Interatomic distances squared
      real(dp), allocatable   :: polR(:,:)  ! Macroscopic polarization in unit-cell coordinates
      real(dp), allocatable   :: polxyz(:,:)! Macroscopic polarization in cartesian coordinates
      real(dp), allocatable   :: va(:,:)    ! Atomic velocities
      real(dp), allocatable   :: xij(:,:)   ! Ineratomic vectors
      real(dp), allocatable   :: xijo(:,:)  ! Vectors between basis-orbital centers


#     ifdef MPI
      integer :: MPIerror  ! Return error code in MPI routines
      integer :: ntmp      ! Temporary to call an MPI routine

      real(dp):: buffer1        ! Temporary for the result of a globalize operation
      real(dp):: qtmp(4)        ! Temporary for globalized spin charge
      real(dp):: stresstmp(3,3) ! Temporary for globalized stress

      real(dp), allocatable:: fatmp(:,:) ! Temporary for globalized atomic forces
#     endif

      logical:: auxchanged ! Has auxiliary supercell changed?
      logical:: dminit     ! Initialize density matrix?
      logical:: found      ! Found density matrix file?
      logical:: foundxv    ! Found position-velocity file?
      logical:: foundzm    ! Found Z-matrix file?
      logical:: gamma      ! Gamma point (k=0) in BZ only?
      logical:: genlogic   ! Generate local density of states?
      logical:: itest      ! Temporary used to check if DM have just been read
      logical:: lastst     ! Temporary to check if this is last geometry step
      logical:: mmix       ! Same as mix. Used in pulayx
      logical:: overflow   ! Too short array sizes?
      logical:: overflowed ! Have array sizes been increased?
    
      logical:: eggbox_block=.true. ! Read eggbox info from data file?
      logical:: final=.false.       ! Last geometry step?

      character(len=150) :: line    ! Data line read
      character(len=10)  :: message ! Not used
      character(len=150) :: sname   ! System name, used to initialise read
      character(len=150) :: shape*10


      !character(len=label_length+5), external :: paste

      !type(parsed_line), pointer  :: p


      CONTAINS
      
      subroutine siesta_write_forces()
      USE siesta_options
      use siesta_geom
      use siesta_cml
      use units
      implicit none
      ! Almost the same forces output whether during simulation
      ! or at the end. Unfortunately not quite, therefore slightly
      ! tortuous logic below. If we are content to change format
      ! of output file slightly, this can be simplified.
      if (.not.final) then
        ! print forces to xml every step. 
        ! output forces to stdout depending on writef
        if (cml_p) then
          call cmlStartPropertyList(mainXML, title='Forces')
          call cmlAddProperty(xf=mainXML, property=fa*Ang/eV,
     .         dictref='siesta:forces', units='siestaUnits:evpa')
          call cmlAddProperty(xf=mainXML, property=ftot,
     .         dictref='siesta:ftot')
          call cmlAddProperty(xf=mainXML, property=fmax, 
     .         dictref='siesta:fmax')
          call cmlAddProperty(xf=mainXML, property=fres,
     .         dictref='siesta:fres')
          call cmlAddProperty(xf=mainXML, property=cfmax, 
     .         dictref='siesta:cfmax')
          call cmlEndPropertyList(mainXML)
        endif
        write(6,'(/,a)') 'siesta: Atomic forces (eV/Ang):'
        if (writef) then
          write(6,'(i6,3f12.6)')(ia,(fa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
        else
          call iofa( na_u, fa )
        endif
        write(6,'(40("-"),/,a6,3f12.6)') 'Tot',(ftot(ix)*Ang/eV,ix=1,3)
        write(6,'(40("-"),/,a6, f12.6)') 'Max',fmax*Ang/eV
        write(6,'(a6,f12.6,a)')'Res',fres*Ang/eV,
     .       '    sqrt( Sum f_i^2 / 3N )'
        write(6,'(40("-"),/,a6, f12.6,a)') 'Max',cfmax*Ang/eV, 
     .       '    constrained'
      else !not final
! In finalization, only print forces if sufficiently large.
        fmax = maxval(abs(fa))
        ftot = sum(fa, dim=2)
        if (fmax .gt. ftol) then
          write(6,'(/,a)') 'siesta: Atomic forces (eV/Ang):'
          write(6,'(a,i6,3f12.6)')
     .         ('siesta: ', ia,(fa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
          write(6,'(a,40("-"),/,a,a6,3f12.6)')
     .         'siesta: ','siesta: ','Tot',(ftot(ix)*Ang/eV,ix=1,3)
          if (cml_p) then
            call cmlStartPropertyList(mainXML, title='Force Summary')
            call cmlAddProperty(xf=mainXML, property=fa*Ang/eV,
     .           dictref='siesta:forces', units='siestaUnits:evpa')
            call cmlAddProperty(xf=mainXML, property=ftot*Ang/eV, 
     .           dictref='siesta:ftot', units='siestaUnits:evpa')
            call cmlEndPropertyList(mainXML)
          endif !cml_p
        endif
        if (Any(cfa /= fa)) then
          fmax = maxval(abs(cfa))
          ftot = sum(cfa, dim=2)
          if (fmax .gt. ftol) then
            write(6,'(/,a)') 'siesta: Constrained forces (eV/Ang):'
            write(6,'(a,i6,3f12.6)')
     .           ('siesta: ',ia,(cfa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
            write(6,'(a,40("-"),/,a,a4,3f12.6)')
     .           'siesta: ','siesta: ','Tot',(ftot(ix)*Ang/eV,ix=1,3)
            if (cml_p) then
              call cmlStartPropertyList(mainXML,
     .             title='Constrained Force Summary')
              call cmlAddProperty(xf=mainXML, property=cfa*Ang/eV, 
     .             dictref='siesta:cforces', units='siestaUnits:evpa')
              call cmlAddProperty(xf=mainXML, property=ftot*Ang/eV, 
     .             dictref='siesta:cftot', units='siestaUnits:evpa')
              call cmlEndPropertyList(mainXML)
            endif !cml_p
          endif
        endif
      endif !final for forces

      end subroutine siesta_write_forces

      subroutine siesta_write_stress_pressure()
      use parallel, only: IOnode
      USE siesta_options
      use siesta_geom
      use atomlist, only: iza
      use m_iostruct,   only: write_struct
      use siesta_cml
      use units
      implicit none
! Stress tensor and pressure:
      

      if (.not.final) then
!
!           Write Voigt components of total stress tensor 
!
            write(6,'(/,a,6f12.2)')
     .           'Stress-tensor-Voigt (kbar):',
     .           (tstress(jx,jx)/kbar,jx=1,3),
     $            tstress(1,2)/kbar,
     $            tstress(2,3)/kbar,
     $            tstress(1,3)/kbar
!
!      Write "target enthalpy" (E + pV, where p is the *target* pressure)
            write(6,"(a,f14.4)") "Target enthalpy (eV/cell)",
     $           (FreeE + tp*volume)/eV
!
!      Write out structural information in "crystallography" format
!
            call write_struct( ucell, na_u, isa, iza, xa )

      ! Output depends on dynamics option
        select case (idyn)
        case(0:5,8)
          if (idyn==0 .and. (.not.varcel)) then
            continue
          else
            write(6,'(/,a,3(/,a,3f12.6))')
     .           'siesta: Stress tensor (static) (eV/Ang**3):',
     .           ('     ',(stress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
            Psol = - ((stress(1,1) + stress(2,2) + stress(3,3))/3.0_dp)
            write(6,'(/,a,f20.8,a)')
     .           'siesta: Pressure (static):', Psol/kBar, '  kBar'
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, property=stress*Ang**3, 
     .             dictref='siesta:stress')
              call cmlAddProperty(xf=mainXML, property=Psol, 
     .             dictref='siesta:psol', title='Pressure (Static)')
            endif !cml_p
!
            write(6,'(/,a,3(/,a,3f12.6))')
     .           'siesta: Stress tensor (total) (eV/Ang**3):',
     .           ('     ',(tstress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
            Psol = - ((tstress(1,1)+tstress(2,2) +tstress(3,3))/3.0_dp)
            write(6,'(/,a,f20.8,a)')
     .           'siesta: Pressure (total):', Psol/kBar, '  kBar'
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, property=tstress*Ang**3, 
     .             dictref='siesta:tstress')
              call cmlAddProperty(xf=mainXML, property=Psol,
     .             dictref='siesta:tpsol', title='Pressure (Total)')
            endif !cml_p
          endif !varcel
        ! Write Force Constant matrix if FC calculation ...
        case(6)
          call ofc(fa,dx,na_u)
        case(7)
!          call phonon_write_forces(fa,na_u,ns,ucell,istep)
           if (IOnode) write(*,*) 'phonon suport desactivated'
        end select !idyn

      else !final

! Print stress tensor unconditionally
        write(6,'(/,a,3(/,a,3f12.6))')
     .       'siesta: Stress tensor (static) (eV/Ang**3):',
     .       ('siesta: ',(stress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
        if (cml_p) then
          call cmlAddProperty(xf=mainXML, property=stress*Ang**3/eV, 
     .         dictref='siesta:stress', units='siestaUnits:eV_Ang__3')
        endif !cml_p

! Print constrained stress tensor if different from unconstrained
        if (Any(cstress /= stress )) then
             write(6,'(/,a,3(/,a,3f12.6))')
     .       'siesta: Constrained stress tensor (static) (eV/Ang**3):',
     .       ('siesta: ',(cstress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
          if (cml_p) then
            call cmlAddProperty(xf=mainXML, property=cstress*Ang**3/eV, 
     .           dictref='siesta:cstress', 
     .           units='siestaUnits:eV_Ang__3')
          endif !cml_p
        endif

! Find pressure
        virial = 0.0_dp
        do ix = 1,3
          fmean = 0.0_dp
          do ia = 1,na_u
            fmean = fmean + fa(ix,ia) / na_u
          enddo
          do ia = 1,na_u
            virial = virial + xalast(ix,ia) * (fa(ix,ia) - fmean)
          enddo
        enddo
        Psol = - (( stress(1,1) + stress(2,2) + stress(3,3) )/3.0_dp)
        Pmol = Psol - virial / volume / 3.0_dp
        write(6,'(/,a,f18.6,a)')
     .       'siesta: Cell volume =', volume/Ang**3, ' Ang**3'
        write(6,'(/,a,/,a,2a20,a,3(/,a,2f20.8,a))')
     .       'siesta: Pressure (static):',
     .       'siesta: ','Solid',        'Molecule',      '  Units',
     .       'siesta: ', Psol,           Pmol,           '  Ry/Bohr**3',
     .       'siesta: ', Psol*Ang**3/eV, Pmol*Ang**3/eV, '  eV/Ang**3',
     .       'siesta: ', Psol/kBar,      Pmol/kBar,      '  kBar'
        if (cml_p) then
          call cmlStartPropertyList(mainXML, title='Final Pressure')
          call cmlAddProperty(xf=mainXML, property=volume/Ang**3, 
     .         title='cell volume', dictref='siesta:cellvol', 
     .         units='siestaUnits:Ang__3')
          call cmlAddProperty(xf=mainXML, property=Psol/kBar, 
     .         title='Pressure of Solid', dictref='siesta:pressSol', 
     .         units='siestaUnits:kbar')
          call cmlAddProperty(xf=mainXML, property=Pmol/kBar,       
     .         title='Pressure of Molecule', dictref='siesta:pressMol', 
     .         units='siestaUnits:kbar')
          call cmlEndPropertyList(mainXML)
        endif !cml_p
        
      endif !final for stress & pressure

      end subroutine siesta_write_stress_pressure
      
      
      subroutine siesta_write_energies( first, iscf )
      USE siesta_options
      use siesta_cml
      use units
      implicit none
      logical :: first
      integer :: iscf
      ! Only print out full decomposition at very beginning and end.
      if ((istp==1.and.first).or.final) then
        write(6,'(/,a,/,(a,f17.6))')
     .     'siesta: Program''s energy decomposition (eV):',
     .     'siesta: Eions   =', Eions/eV,
     .     'siesta: Ena     =', Ena/eV,
     .     'siesta: Ekin    =', Ekin/eV,
     .     'siesta: Enl     =', Enl/eV,
     .     'siesta: DEna    =', DEna/eV,
     .     'siesta: DUscf   =', DUscf/eV,
     .     'siesta: DUext   =', DUext/eV,
     .     'siesta: Exc     =', Exc/eV,
     .     'siesta: eta*DQ  =', Ecorrec/eV,
     .     'siesta: Emadel  =', Emad/eV,
     .     'siesta: Emeta   =', Emeta/eV,
     .     'siesta: Emolmec =', Emm/eV,
     .     'siesta: Ekinion =', Ekinion/eV,
     .     'siesta: Eharris =', (Eharrs1+Ekinion)/eV,
     .     'siesta: Etot    =', (Etot+Ekinion)/eV,
     .     'siesta: FreeEng =', (FreeE+Ekinion)/eV
        if (cml_p) then
          call cmlStartPropertyList(mainXML,
     .         title='Energy Decomposition')
          call cmlAddProperty(xf=mainXML, property=Eions/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Eions', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Ena/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ena', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Ekin/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ekin', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Enl/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Enl', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=DEna/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:DEna', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=DUscf/eV, 
     .         units='siestaUnits:eV',
     .         dictref='siesta:DUscf', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=DUext/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:DUext', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Exc/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Exc', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML,property=Ecorrec/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ecorrec', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Emad/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Emad', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Emeta/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Emeta', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Emm/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Emm', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML,property=Ekinion/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ekinion', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=(Eharrs1+Ekinion)/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:EharrsK', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=(Etot+Ekinion)/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:EtotK', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=(FreeE+Ekinion)/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:FreeEK', fmt='(f17.6)')
          call cmlEndPropertyList(mainXML)
        endif
      endif
      ! On all SCF steps, print out the current energy (format depending on type of run)
      if (.not.final) then
        ! Print total energy and density matrix error .........................
        if (cml_p) then
          call cmlStartPropertyList(mainXML, title='SCF Cycle')
          ! Eharrs is always output
          call cmlAddProperty(xf=mainXML, property=Eharrs/eV, 
     .         units="siestaUnits:eV", 
     .         dictRef="siesta:Eharrs", fmt="(f14.7)")
        endif
        ! This chain of if statements determines which properties are output.
        if (harrisfun) then 
          write(6,"(/a,f14.6,/)") 'siesta: Eharris(eV) = ', Eharrs/eV
            ! No need for further cml output
        elseif (isolve==0) then
          if (cml_p) 
     .         call cmlAddProperty(xf=mainXML, property=FreeE/eV, 
     .         units="siestaUnits:eV", 
     .         dictRef="siesta:FreeE",  fmt="(f14.7)")
          if (fixspin) then
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, property=Etot/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:Etot",   fmt="(f14.7)")
              call cmlAddProperty(xf=mainXML, property=FreeE/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:FreeE", fmt="(f14.7)")
              call cmlAddProperty(xf=mainXML, property=dDmax/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:dDmax",  fmt="(f14.7)")
            endif
            if ((iscf .eq. 1).or.muldeb)
     .           write(6,'(/,a12,3a14,a8,a7,a11)')
     .           'siesta: iscf', '   Eharris(eV)', 
     .           '      E_KS(eV)', '   FreeEng(eV)', 
     .           '   dDmax', '  Ef_up', '  Ef_dn(eV)'
            write(6,'(a8,i4,3f14.4,f8.4,2f9.4)')
     .           'siesta: ',iscf, Eharrs/eV, Etot/eV, FreeE/eV, dDmax, 
     .           (Efs(i)/eV,i=1,2)
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, property=Efs(1)/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:Efs",    fmt="(f14.7)")
              call cmlAddProperty(xf=mainXML, property=Efs(2)/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:Efs",    fmt="(f14.7)")
            endif
          else !fixspin
            if ((iscf .eq. 1).or.muldeb)
     .           write(6,'(/,a12,3a14,2a8)')
     .           'siesta: iscf', '   Eharris(eV)', 
     .           '      E_KS(eV)', '   FreeEng(eV)', 
     .           '   dDmax', '  Ef(eV)'
            write(6,'(a8,i4,3f14.4,2f8.4)')
     .           'siesta: ',iscf, Eharrs/eV, Etot/eV, FreeE/eV, 
     .           dDmax, Ef/eV
            if (cml_p) then
              call cmlAddProperty(xf=mainXML,property=Ef/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:Ef", fmt="(f14.7)")
            endif !cml_p
          endif !fixspin
        elseif (isolve==1) then
          write(6,'(/,a15,i4)') 'siesta: iscf = ',iscf
          write(6,'(a14,f15.4,a13,f15.4,a10,f7.4/)') 
     .         'Eharris(eV) = ',Eharrs/eV,
     .         '  E_KS(eV) = ',Etot/eV,'  dDmax = ',dDmax
          if (cml_p) then
            call cmlAddProperty(xf=mainXML, property=Etot/eV, 
     .           units="siestaUnits:eV", 
     .           dictRef="siesta:Etot",   fmt="(f14.7)")
            call cmlAddProperty(xf=mainXML, property=dDmax/eV, 
     .           units="siestaUnits:eV", 
     .           dictRef="siesta:dDmax",  fmt="(f14.7)")
          endif
        endif !harrisfun/isolve
        
        if (cml_p) then
          call cmlEndPropertyList(mainXML)
        endif
        
      else !final
      ! Print out additional information in finalization.

        write(6,'(/,a)') 'siesta: Final energy (eV):'
        write(6,'(a,a15,f15.6)')
     .    'siesta: ',      'Kinetic =', Ekin/eV,
     .    'siesta: ',      'Hartree =', Uscf/eV,
     .    'siesta: ',   'Ext. field =', DUext/eV,
     .    'siesta: ',  'Exch.-corr. =', Exc/eV,
     .    'siesta: ', 'Ion-electron =', (Enascf+Enl+DUscf-Uscf-Uatm)/eV,
     .    'siesta: ',      'Ion-ion =', (Ena+Uatm-Enaatm-Eions)/eV,
     .    'siesta: ',      'Ekinion =', Ekinion/eV,
     .    'siesta: ',        'Total =', (Etot+Ekinion)/eV
        if (cml_p) then
          call cmlStartPropertyList(xf=mainXML, title='Final Energy')
          call cmlAddProperty(xf=mainXML, property=Ekin/eV,  
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ekin',  fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Uscf/eV,  
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Uscf',  fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=DUext/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:DUext', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Exc/eV,   
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Exc',   fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, 
     .         property=(Enascf+Enl+DUscf-Uscf-Uatm)/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:I-e', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, 
     .         property=(Ena+Uatm-Enaatm-Eions)/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:I-I', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Ekinion/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ekinion', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=(Etot+Ekinion)/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Etot', fmt='(f17.6)')
          call cmlEndPropertyList(mainXML)
        endif !cml_p
      endif !final

      end subroutine siesta_write_energies

      subroutine siesta_write_positions()
      use siesta_geom
      use atomlist, only: elem
      use siesta_cml
      use units
      implicit none
        if (cml_p) then
          call cmlAddMolecule(xf=mainXML, natoms=na_u, elements=elem,
     .           refs=cisa, coords=xa/Ang, style='x3', fmt='(f12.6)')
          call cmlAddLattice(xf=mainXML, cell=ucell/Ang, 
     .           units='siestaUnits:Ang', dictref='siesta:ucell')
        endif
      end subroutine siesta_write_positions

      END MODULE SIESTA_TODO
