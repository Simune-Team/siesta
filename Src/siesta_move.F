      MODULE m_siesta_move
      private
      public :: siesta_move

      CONTAINS

      subroutine siesta_move( istep, relaxd )
      USE SIESTA_TODO
      use m_born_charge, only: born_charge

      implicit none

      integer  :: istep
      logical  :: relaxd

      external :: automatic_cell,
     .  bonds, cgvc, cgvc_zmatrix, fixed,
     .  dhscf, diagon, dnaefs, extrapol, initatom,
     .  iodm, iozm,
     .  kinefsm, mulliken, naefs, neighb, !kgrid, kgridinit, 
     .  pulayx, 
     .  reinit, shaper, spnvec, 
     .  timer, xijorb, memory,
     .  ioeig, iofa, iokp, iomd, prversion !, eggbox


!------------------------------------------------------------------------- BEGIN
! Write final Kohn-Sham Energy ........................................
      if (cml_p) call cmlStartPropertyList(mainXML,
     .                                   title='Final KS Energy')
      if (IOnode) then
        if ( .not. harrisfun) 
     .    write(6,"(/a,f14.4)")  'siesta: E_KS(eV) =        ', Etot/eV
        if (cml_p) call cmlAddProperty(xf=mainXML, property=Etot/eV,
     .       dictref='siesta:E_KS', units='siestaUnits:eV', 
     .       fmt='(f14.6)')
      endif

! Substract egg box effect form energy ................................
      if (eggbox_block) then
        call eggbox('energy',ucell,na_u,isa,ntm,xa,fa,Etot,
     .               eggbox_block)
        if (IOnode)
     .    write(6,"(/a,f14.4)") 'siesta: E_KS - E_eggbox = ',Etot/eV
        if (cml_p) call cmlAddProperty(xf=mainXML, property=Etot/eV,
     .         dictref='siesta:E_KS_egg', units='siestaUnits:eV', 
     .         fmt='(f14.6)')
      endif
      if (cml_p) call cmlEndPropertyList(mainXML)

#ifdef MPI
! Global reduction of forces and stresses
      allocate(fatmp(3,na_u))
      call memory('A','D',3*na_u,'siesta')
      call globalize_sum(stressl(1:3,1:3),stresstmp(1:3,1:3))
      call globalize_sum(fal(1:3,1:na_u),fatmp(1:3,1:na_u))
      stress(1:3,1:3) = stress(1:3,1:3) + stresstmp(1:3,1:3)
      fa(1:3,1:na_u) = fa(1:3,1:na_u) + fatmp(1:3,1:na_u)
      call memory('D','D',size(fatmp),'siesta')
      deallocate(fatmp)
#else
      stress(1:3,1:3) = stress(1:3,1:3) + stressl(1:3,1:3)
      fa(1:3,1:na_u) = fa(1:3,1:na_u) + fal(1:3,1:na_u)
#endif

! Substract egg box effect from the forces ............................
      if (eggbox_block) then
        call eggbox('forces',ucell,na_u,isa,ntm,xa,fa,Etot,eggbox_block)
      endif
! ...................

! Impose constraints to atomic movements by changing forces ...........
      call fixed(ucell,stress,na_u,isa, amass, xa, fa, cstress, cfa, 
     .           ntcon )
! ...................

! Write atomic forces .................................................
      fmax = 0.0_dp
      cfmax = 0.0_dp
      fres = 0.0_dp
      do ix = 1,3
        ftot(ix) = 0.0_dp
        do ia = 1,na_u
          ftem = fa(ix,ia)
          cftem = cfa(ix,ia)
          ftot(ix) = ftot(ix) + ftem
          fres = fres + ftem*ftem
          fmax = max( fmax, dabs(ftem) )
          cfmax = max( cfmax, dabs(cftem) )
        enddo
      enddo
      fres = dsqrt( fres / (3.0_dp*na_u) )

! Calculate and output Zmatrix forces
      if (lUseZmatrix) then
        call CartesianForce_to_ZmatForce(na_u,xa,fa)
        if (IOnode) call iofaZmat()
      endif

! Add kinetic term to stress tensor ..................................
      do ia = 1,na_u
        do jx = 1,3
          do ix = 1,3
            tstress(ix,jx) = stress(ix,jx) -
     .             amu * amass(ia) * va(ix,ia) * va(jx,ia) / volume
          enddo
        enddo
      enddo

! Force output .......................................................
      if (IOnode) then
        call siesta_write_forces()
        call siesta_write_stress_pressure()
        call wallclock('--- end of geometry step')
      endif

! Mulliken population analysis .......................................
      call mulliken( mullipop, nspin, na_u, no_u, maxnh,
     .               numh, listhptr, listh, S, Dscf, isa, 
     .               lasto, iaorb, iphorb )

!     Call the born effective charge routine only in those steps (even) 
!     in which the dx  is positive.
      if (bornz .and. (mod(istep,2) .eq. 0)) then
         call born_charge()
      endif

! Save the last coordinates for which the density matrix has been calculated
      if ( (istep.eq.fincoor)) then
        xalast(1:3,1:na_s)=xa(1:3,1:na_s)
      endif

      Ekinion  = 0.0_dp
      vn       = 0.0_dp
      vpr      = 0.0_dp
      kn       = 0.0_dp
      kpr      = 0.0_dp

      iunit = 2

! Move atoms ..........................................................
      select case(idyn)
      case(0)
        if (nmove .ne. 0) then
          if (lUseZmatrix) then
             if (broyden_optim)
     $          write(6,'(a)') 'siesta: Broyden method ' //
     $            'not implemented yet for Zmatrix relaxation'
             call cgvc_zmatrix( na_u, xa, cfa, ucell, cstress,
     $            volume, dxmax, tp, ftol, strtol, varcel,
     $            relaxd, usesavecg )
          else
             if (broyden_optim) then
                call broyden_optimizer( na_u, xa, cfa, ucell,
     $               cstress, volume, dxmax, tp, ftol, strtol,
     $               varcel, relaxd )
             else
                call cgvc( na_u, xa, cfa, ucell, cstress, volume,
     $               dxmax, tp, ftol, strtol, varcel,
     $               relaxd, usesavecg )
             endif
          endif
          ! Propagate the new structure to the virtual supercell
          call superx( ucell, nsc, na_u, na_s, xa, scell )
          if (relaxd) goto 60
            ! Exit coordinate relaxation loop
        endif

      case(1)
         call verlet2(istp, iunit, iquench, na_u, cfa, dt,
     .       amass, ntcon, va, xa, Ekinion, tempion)
        ! Propagate the new structure to the virtual supercell
         call superx( ucell, nsc, na_u, na_s, xa, scell )
        ! Check convergence for quenching runs (which are
        ! really relaxations)
         if ( iquench .ne. 0 ) then
            relaxd = .true.
            do ia = 1, na_u
              do i = 1, 3
                relaxd = relaxd .and. ( abs(cfa(i,ia)) .lt. ftol )
              enddo
            enddo
         endif

      case (2)
         call nose(istp, iunit, na_u, cfa, tt, dt, amass, mn,
     .       ntcon, va, xa, Ekinion, kn, vn, tempion)
        ! Propagate the new structure to the virtual supercell
        call superx( ucell, nsc, na_u, na_s, xa, scell )

      case (3)
         call pr(istp, iunit, iquench, na_u, cfa, cstress, tp, dt,
     .           amass, mpr, ntcon, va, xa, vcell, ucell, Ekinion, 
     .           kpr, vpr, tempion, Pint)
        ! Propagate the new structure to the virtual supercell
        call superx( ucell, nsc, na_u, na_s, xa, scell )
        if (IOnode) write(6,'(/,a,f12.3,a)')
     .       'siesta: E_kin PR =', kpr/Kelvin, ' K'

      case (4)
         call npr(istp, iunit, na_u, cfa, cstress, tp, tt, dt,
     .            amass, mn, mpr, ntcon, va, xa, vcell, ucell, 
     .            Ekinion, kn, kpr, vn, vpr, tempion, Pint)
        ! Propagate the new structure to the virtual supercell
        call superx( ucell, nsc, na_u, na_s, xa, scell )

      case (5)
         call anneal(istp, iunit, ianneal, taurelax, bulkm,
     .       na_u, cfa, cstress, tp, tt, dt, amass, ntcon,
     .       va, xa, ucell, Ekinion, tempion, Pint)
        ! Propagate the new structure to the virtual supercell
        call superx( ucell, nsc, na_u, na_s, xa, scell )

      case (6)

        ! Save last atomic positions and velocities 
        call ioxv( 'write', ucell, vcell, na_u, isa, iza,  
     &             xa, va, foundxv)
        if (lUseZmatrix) then
          call iozm('write',ucell,vcell,xa,foundzm)
        endif
        call siesta_write_positions()

        ! Save atomic positions and velocities accumulatively and
        ! accumulate coor in Xmol file for animation 
        call write_md_record( istep )

        ! Undo the last atom displacement
        if (istep .gt. inicoor) then
          iadispl = (istep-mod(istep-1,6))/6+ia1
          ix = mod(istep-1,6)+1
          ixdispl = (ix - mod(ix-1,2) +1)/2
          xa(ixdispl,iadispl) = xa(ixdispl,iadispl) - dx
        endif

        ! Displace atom by dx
        if (istep .lt. fincoor) then
          iadispl = ((istep+1)-mod(istep,6))/6+ia1
          ix = mod(istep,6)+1
          ixdispl = (ix - mod(ix-1,2) +1)/2
          dx=-dx
          xa(ixdispl,iadispl) = xa(ixdispl,iadispl) + dx
        endif

      case (7)
          ! Save last atomic positions and velocities 
          call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, 
     &                foundxv)
          if (lUseZmatrix) then
            call iozm('write',ucell,vcell,xa,foundzm)
          endif
          call siesta_write_positions()

          ! Save atomic positions and velocities accumulatively and
          ! accumulate coor in Xmol file for animation 
          call write_md_record( istep )

          ! Restore original coordinates after FC displacements
          if (istep.gt.inicoor) then
            call phonon_restore_coords(istep,xa,ucell)
          endif

          ! Modify coordinates
          if (istep .lt. fincoor) then
            call phonon_set_coords(istep+1,xa,ucell)
          endif

      case (8)
          call forcesToPipe( na_u, Etot, cfa, cstress )

          ! Save last atomic positions and velocities 
          call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, 
     &             foundxv)
          if (lUseZmatrix) then
            call iozm('write',ucell,vcell,xa,foundzm)
          endif
          call siesta_write_positions()

          ! Save atomic positions and velocities accumulatively and
          ! accumulate coor in Xmol file for animation 
          call write_md_record( istep )

          ! Get coordinates from driver program trough pipe
          call coordsFromPipe( na_u, xa, ucell )
          if (volcel(ucell) < 1.0e-8_dp) then
            call automatic_cell(ucell,scell,na_u,xa,isa,charnet)
          endif
      end select

      if (IOnode) then
        if (idyn .gt. 0 .and. idyn .lt. 6) then
          write(6,'(/,a,f12.3,a)')
     .      'siesta: Temp_ion =', tempion, ' K'
        endif
      endif

      ! Save last atomic positions and velocities 
      if (idyn.ne.6 .and. idyn.ne.7 .and. idyn.ne.8) then
        call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, 
     &             foundxv)
        if (lUseZmatrix)
     .    call iozm('write',ucell,vcell,xa,foundzm)
        call siesta_write_positions()


        ! Save atomic positions and velocities accumulatively and
        ! accumulate coor in Xmol file for animation 
        call write_md_record( istep )
      endif


   60 continue
! Output memory use at the end of this geometry step
      if (cml_p) call cmlEndStep(mainXML)
      call printmemory( 6, 0 )
      call timer( 'IterMD', 2 )

! End of one MD step - flush stdout
      if (ionode) call pxfflush(6)

!--------------------------------------------------------------------------- END
      END subroutine siesta_move

      END MODULE m_siesta_move
