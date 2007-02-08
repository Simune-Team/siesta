      Module SIESTA_TODO

      use precision, only: dp
      use parallel, only: Node, Nodes, IOnode, ParallelOverK
      use m_cell,   only: ucell, ucell_last
      use atmfuncs, only: rcut, uion, izofis
      use atomlist, only: xa, xalast, indxuo, na_u
      use atomlist, only: isa, cisa, elem, indxua, iza, lastkb, lasto
      use atomlist, only: qa, amass, qtot, no_u, rmaxo, rmaxkb
      use atomlist, only: zvaltot
      use atomlist, only: rmaxv, na_s, iphorb, iphkb, rco, rckb
      use atomlist, only: datm, no_s, iaorb
      use atomlist, only: superx, superc, initatomlists
!
      use band, only: setup_bands, bands
      use fdf, only : fdf_block, fdf_convfac, fdf_boolean
      use fdf, only : fdf_physical
      use m_fdf_global, only: fdf_global_get
      use sys, only: die, bye
      use parse
      use periodic_table, only : symbol
      use xcmod, only: setXC
      use molecularmechanics, only : inittwobody, twobody
      use metaforce, only: lMetaForce, initmeta, meta
#ifdef MPI
      use mpi_siesta, only: mpi_comm_world
      use m_mpi_utils, only : globalize_sum, globalize_max
#endif
      use m_mpi_utils, only : broadcast
      use alloc, only: re_alloc, alloc_report
      use phonon, only: phonon_num_disps, phonon_setup
      use phonon, only: phonon_write_forces, phonon_restore_coords
      use phonon, only: phonon_set_coords

      use densematrix

      use m_ordern,     only: ordern
      use m_hsparse,    only: hsparse

      use parallelsubs, only: getnodeorbs
      use writewave,    only: initwave, wwave
      use iopipes,      only: forcestopipe, coordsfrompipe
      use m_iostruct,   only: write_struct, read_struct
      use m_nlefsm,     only: nlefsm
      use m_overfsm,    only: overfsm
      use m_check_supercell, only: check_sc_factors
 
      use m_broyden_mixing

      use siesta_cml

      use files,       only : slabel, label_length

      use m_timestamp, only : timestamp
      use m_wallclock, only : wallclock
      use units
      use zmatrix,      only: lUseZmatrix, iofaZmat
      use zmatrix,      only: CartesianForce_to_ZmatForce
      use zmatrix,      only: write_Zmatrix
      use m_broyden_optim,      only:broyden_optimizer
      use m_redata,     only: redata
      use m_ioxv,       only: ioxv 

      use m_smearing, only: temp
      use m_dynamics, only: nose, verlet2, npr, anneal, pr
      use md_out,     only: md_v_format
#ifdef CDF
      use md_out,     only: md_netcdf
#endif
      use Kpoint_grid

      implicit none

      integer
     .  fincoor,
     .  i, ia, ia1, ia2, iadispl, ianneal,
     .  idyn, ifa, ifinal, ihmat, iiscf, ik, 
     .  ind, inicoor, io, ioptlwf, iord,
     .  iquench, is, isel, 
     .  isolve, ispin, istp, istpsave, istart, istr,
     .  iu, iunit,  iv, ix, ixdispl, 
     .  j, ja, jo, jx,
     .  level,  
     .  maxsav, broyden_maxit, maxwk, mullipop, 
     .  nkick, nbcell, ncgmax, nh,
     .  nmove, nnia, nscf, 
     .  ntm(3), ntcon, nwk, pmax, nkpol, nhist, 
     .  neigwanted, neigmin

      integer :: ihuge=1073741823
      ! maxnh is (initially) the minimum size of all the sparse arrays.
      ! It must be 1 (not 0) since we frequently use the F77 idiom of
      ! passing the first element of the array instead of all of it.
      integer :: maxnh=1
      integer :: maxna=200
      integer :: no_l=1
      integer :: nnamax=1
      integer :: nauxpul=1
      integer :: ns=1   ! AG** do not initialize
      integer :: nspin=1
      integer :: nxij=1

      integer, dimension(3)   :: nscold = 0
      integer, dimension(3)   :: nsc = 1
      integer, dimension(3,3) :: mscell = 
     .                    reshape ((/1,0,0,0,1,0,0,0,1/), (/3,3/))

#ifdef MPI
      integer
     .  MPIerror, ntmp
#endif

      integer, dimension(:), allocatable ::
     .  jna

      integer, pointer, save ::
     .  listh(:), listhold(:), listhptr(:), listhptrold(:),
     .  numh(:), numhold(:)

      real(dp)
     .  bcell(3,3), beta, bulkm, 
     .  charnet, cfmax, cftem, const, cstress(3,3),
     .  dipol(3), dDmax, dDtol, dEmax, dEtol, DEharr, DEna,
     .  dt, DUext, DUscf, Dxc, dx, dxmax, sigma,
     .  Ecorrec, ef, Eharrs, Eharrs1, Eions, Ekin, Ekinion, Emad, Ena, 
     .  Elast, Enaatm, Enascf, Enl, Entrop, Entropy, Emeta,
     .  eta(2), etol, Etot, Etot_output, Exc, E0, Emm,
     .  factor, fmax, fmean, FreeE, fres, ftem, ftol, ftot(3),
     .  g2cut, g2max, getot, kn, kpr, 
     .  mn, mpr, occtol, Pint, Pmol, Psol, rmax_bonds

      real(dp) :: e1=-1.0_dp, e2=-1.0_dp
      real(dp) :: tiny = 1.0E-15_dp
      real(dp) :: tempion=0.0_dp


      real(dp)
     .  qaux, qspin(4), qsol,
     .  rcoor, rcoorcp, rijmin, rmax, rmaxh,
     .  scell(3,3), scell_last(3,3), stot, stress(3,3), strtol, svec(3),
     .  taurelax, tempinit, tp, ts, tt,
     .  tstress(3,3), Uatm, Uscf,
     .  vcell(3,3), virial, vn,
     .  volume, vpr, wmix, wmixkick, 
     .  stressl(3,3), veclen

      real(dp), external:: volcel

      real(dp), dimension(:), allocatable ::
     .  efs, qs, r2ij

      real(dp), pointer, save :: 
     .  H0(:), S(:), wgthpol(:)

      real(dp), dimension(:,:), allocatable :: 
     .  auxpul, cfa, fa, fal,
     .  polR, polxyz, va, xij, xijo, wfk

      real(dp), pointer, save ::
     .  Dold(:,:), Dscf(:,:), Dscfsave(:,:), Eold(:,:), Escf(:,:),
     .  H(:,:), kpol(:,:)

#ifdef MPI
      real(dp)
     .  buffer1, stresstmp(3,3), qtmp(4)

      real(dp), dimension(:,:), allocatable :: 
     .  fatmp
#endif

      real(dp), dimension(:,:,:), allocatable :: 
     .  ebk, eo, qo

      logical
     .  auxchanged, chebef, dminit, default, dumpcharge, 
     .  fixauxcell, fixspin, found, foundxv, foundzm,
     .  gamma, naiveauxcell,
     .  initdmaux, inspn, itest, lastst, mix, mmix, negl, noeta, 
     .  outlng, overflow, overflowed, pulfile,
     .  savehs, savevh, savevt, savdrh, savrho,
     .  savepsch, savetoch,
     .  usesavecg, usesavelwf, usesavedm, usesavedmloc, usesavexv, 
     .  usesavezm, writeig,
     .  writbk, writmd, writpx, writb, writec, writef, 
     .  writek, writic, varcel, genlogic, do_pdos, writedm, atmonly,
     .  harrisfun, muldeb, use_struct_file,
     $  require_energy_convergence, broyden_optim, struct_only, bornz

    
      logical :: eggbox_block=.true.
      logical :: final=.false.

      character
     .  line*150, sname*150, shape*10, message*79

      character(len=label_length+5) :: fildrh
      character(len=label_length+5) :: filepsch
      character(len=label_length+5) :: filetoch
      character(len=label_length+5) :: filevh
      character(len=label_length+5) :: filevt
      character(len=label_length+5) :: filrho
      character(len=label_length+5), external :: paste

      type(parsed_line), pointer  :: p


      CONTAINS
      
      subroutine siesta_write_k_points()
      if ( writek ) then
        write(6,'(/,a)')
     .       'siesta: k-point coordinates (Bohr**-1) and weights:'
        write(6,'(a,i4,3f12.6,3x,f12.6)')
     .       ('siesta: ', ik, (kpoint(ix,ik),ix=1,3), kweight(ik),
     .       ik=1,nkpnt)
        if (cml_p) then
          call cmlAddProperty(xf=mainXML, property=kpoint, 
     .         dictref='siesta:kpoint')
          call cmlAddProperty(xf=mainXML, property=kweight,
     .         dictref='siesta:kweight')
        endif
      else
        call iokp( nkpnt, kpoint, kweight )
      endif
      write(6,'(/a,i6)')
     .     'siesta: k-grid: Number of k-points =', nkpnt
      write(6,'(a,f10.3,a)')
     .     'siesta: k-grid: Cutoff (effective) =', kcutof/Ang, ' Ang'
      write(6,'(a)')
     .     'siesta: k-grid: Supercell and displacements'
      write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .     (kscell(i,1),i=1,3), kdispl(1)
      write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .     (kscell(i,2),i=1,3), kdispl(2)
      write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .     (kscell(i,3),i=1,3), kdispl(3)
      if (cml_p) then
        call cmlStartPropertyList(mainXML, title='k-points')
        call cmlAddProperty(xf=mainXML, property=nkpnt, 
     .       dictref='siesta:nkpnt')
        call cmlAddProperty(xf=mainXML, property=kcutof/Ang, 
     .       dictref='siesta:kcutof', units='siestaUnits:angstrom')
        call cmlEndPropertyList(mainXML)
        call cmlAddProperty(xf=mainXML, property=kscell,
     .       dictref='siesta:kscell')
        call cmlAddProperty(xf=mainXML, property=kdispl,
     .       dictref='siesta:kdispl')
      endif
      end subroutine siesta_write_k_points


      subroutine siesta_write_forces()

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
! Stress tensor and pressure:
      

      if (.not.final) then
!
!           Write Voigt components of total stress tensor 
!
            write(6,'(/,a,6f12.2))')
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

      subroutine siesta_write_positions
        if (cml_p) then
          call cmlAddMolecule(xf=mainXML, natoms=na_u, elements=elem,
     .           refs=cisa, coords=xa/Ang, style='x3', fmt='(f12.6)')
          call cmlAddLattice(xf=mainXML, cell=ucell/Ang, 
     .           units='siestaUnits:Ang', dictref='siesta:ucell')
        endif
      end subroutine siesta_write_positions

      END MODULE SIESTA_TODO
