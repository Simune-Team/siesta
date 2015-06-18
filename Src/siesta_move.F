      MODULE m_siesta_move

      private
      public :: siesta_move

      CONTAINS

      subroutine siesta_move( istep, relaxd )
!      use siesta_geom
      use zmatrix,         only: lUseZmatrix, iofaZmat, write_Zmatrix,
     .                           CartesianForce_to_ZmatForce
      use atomlist,        only: iaorb, iphorb, iza, amass, no_u, 
     .                           lasto, superx
      use m_ioxv,          only: ioxv
      USE siesta_options
      use parallel,        only: IOnode
      use siesta_cml,      only: cml_p, cmlEndStep, mainXML
      use sys,             only: die

      use m_broyden_optim, only: broyden_optimizer
      use m_fire_optim, only: fire_optimizer
      use m_zm_broyden_optim,      only: zm_broyden_optimizer
      use m_zm_fire_optim,      only: zm_fire_optimizer
      use m_cell_broyden_optim,      only: cell_broyden_optimizer
      use m_cell_fire_optim,      only: cell_fire_optimizer

      use m_dynamics,      only: nose, verlet2, npr, anneal, pr
      use phonon,          only: phonon_setup, phonon_restore_coords,
     .                           phonon_set_coords
!      use iopipes,         only: forcestopipe, coordsfrompipe
      use m_energies,      only: Ekinion, Etot
      use write_subs
      use m_steps
      use m_kinetic
      use siesta_geom,     only: na_u       ! Number of atoms in unit cell
      use siesta_geom,     only: na_s       ! Number of atoms in super-cell
      use siesta_geom,     only: nsc        ! Unit cells in each supercell dir
      use siesta_geom,     only: isa        ! Atomic species
      use siesta_geom,     only: xa         ! Atomic positions
      use siesta_geom,     only: xa_last    ! Previous atomic positions
      use siesta_geom,     only: va         ! Atomic velocities
      use siesta_geom,     only: ucell      ! Unit cell vectors
      use siesta_geom,     only: ucell_last ! Previous unit cell vectors
      use siesta_geom,     only: scell      ! Super-cell vectors
      use siesta_geom,     only: scell_last ! Previous super-cell vectors
      use siesta_geom,     only: vcell      ! Time derivative of unit cell vecs
      use m_forces,        only: cfa        ! Constrained forces
      use m_forces,        only: ntcon      ! Number of geometry constraints
      use m_stress,        only: cstress    ! Constrained stress tensor
      use siesta_master,   only: coordsFromMaster ! Get coords from master prog

#ifdef NCDF_4
      use m_exp_coord, only : exp_coord_next
#endif
#ifdef SIESTA__FLOOK
      use flook_siesta, only : slua_call, LUA_MOVE
#endif

      implicit none

      integer, intent(in)   :: istep
      logical, intent(out)  :: relaxd
      integer               :: i, ia, ix, iadispl, ixdispl
      integer               :: iunit   ! Physical-units option: 
                                       ! 1=>(eV,Ang), 2=>(Ry,Bohr)
      real(dp)              :: Pint    ! Instantaneous pressure 
                
      logical               :: foundxv  ! dummy for call to ioxv
      logical               :: foundzm  ! dummy for call to iozm

      real(dp), external :: volcel


!------------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( 'PRE siesta_move' )
#endif
      ! Save the last geometry for which the density matrix 
      ! (electronic structure) has been calculated
      call timer( 'SIMOVE', 1 )


!$OMP parallel workshare default(shared)
      xa_last(1:3,1:na_s) = xa(1:3,1:na_s)
      ucell_last(1:3,1:3) = ucell(1:3,1:3)
      scell_last(1:3,1:3) = scell(1:3,1:3)
!$OMP end parallel


      Ekinion = 0.0_dp   ! Fossil: already initialized in siesta_init

      iunit = 2

! Move atoms ..........................................................
      select case(idyn)
      case(0)
        if (nmove .ne. 0) then
          if (RelaxCellOnly) then
             if (broyden_optim) then
                call cell_broyden_optimizer( na_u, xa, ucell,
     $            cstress, tp, strtol,
     $            varcel, relaxd)
             elseif (fire_optim) then
                call cell_fire_optimizer( na_u, xa, ucell,
     $            cstress, tp, strtol,
     $            varcel, relaxd)
             else
                call die("Cell-only optim needs Broyden or FIRE")
             endif

          else   ! Coordinate relaxation (and maybe cell)

           if (lUseZmatrix) then
             if (broyden_optim) then
                call zm_broyden_optimizer( na_u, xa, ucell, cstress, tp,
     &                                     strtol, varcel, relaxd )
             elseif (fire_optim) then
                call zm_fire_optimizer( na_u, xa, cfa, ucell,
     $               cstress, dxmax, tp, ftol, strtol,
     $               varcel, relaxd)
             else
                call cgvc_zmatrix( na_u, xa, cfa, ucell, cstress,
     $               dxmax, tp, ftol, strtol, varcel,
     $               relaxd, usesavecg )
             endif
           else
             if (broyden_optim) then
                call broyden_optimizer( na_u, xa, cfa, ucell,
     $               cstress, tp, ftol, strtol, varcel, relaxd )
             elseif (fire_optim) then
                call fire_optimizer( na_u, xa, cfa, ucell,
     $               cstress, dxmax, tp, ftol, strtol,
     $               varcel, relaxd )
             else
                call cgvc( na_u, xa, cfa, ucell, cstress,
     $               dxmax, tp, ftol, strtol, varcel,
     $               relaxd, usesavecg )
             endif
           endif
          endif ! RelaxCellOnly

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
         call nose( istp, iunit, na_u, cfa, tt, dt, amass, mn,
     .              ntcon, va, xa, Ekinion, kn, vn, tempion )
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
        if ( istep > inicoor ) then ! inicoor == 0
          iadispl = (istep-mod(istep-1,6))/6+ia1
          ix = mod(istep-1,6)+1
          ixdispl = (ix - mod(ix-1,2) +1)/2
          xa(ixdispl,iadispl) = xa(ixdispl,iadispl) - dx
        end if

        ! Displace atom by dx
        ! NOTE: MOVE is before istep increment, hence
        ! fincoor == istep will be the last step and 
        ! the initial geometry should be retained.
        if ( istep < fincoor ) then ! fincoor = (ia2-ia1+1)*6
          iadispl = ((istep+1)-mod(istep,6))/6+ia1
          ix = mod(istep,6)+1
          ixdispl = (ix - mod(ix-1,2) +1)/2
          dx = -dx
          xa(ixdispl,iadispl) = xa(ixdispl,iadispl) + dx
        end if

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
!          Moved to siesta_forces. JMS Oct.2010
!          call forcesToPipe( na_u, Etot, cfa, cstress )

          ! Save last atomic positions and velocities 
          call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, 
     &             foundxv)
!!          if (lUseZmatrix) then
!!            call iozm('write',ucell,vcell,xa,foundzm)
!!          endif
          call siesta_write_positions()

          ! Save atomic positions and velocities accumulatively and
          ! accumulate coor in Xmol file for animation 
          call write_md_record( istep )

          ! Get coordinates from driver program trough pipe
!          Changed by JMS. Oct.2010
!          call coordsFromPipe( na_u, xa, ucell )
          call coordsFromMaster( na_u, xa, ucell )
          if (volcel(ucell) < 1.0e-8_dp) then
            call automatic_cell(ucell,scell,na_u,xa,isa,charnet)
          endif
#ifdef NCDF_4
       case (9)
          ! Save the md-record
          call write_md_record( istep )
          if ( istep < fincoor ) then
             call exp_coord_next( istep + 1 , na_u, xa)
          end if
          ! Force the velocity to zero to not move the 
          ! atoms
          va = 0._dp
          goto 60
#endif
#ifdef SIESTA__FLOOK
        case (10)
           ! Communicate with lua
           call slua_call(LUA, LUA_MOVE)
#endif
      end select

      if (IOnode) then
        if (idyn .gt. 0 .and. idyn .lt. 6) then
          write(6,'(/,a,f12.3,a)')
     .      'siesta: Temp_ion =', tempion, ' K'
        endif
      endif

      ! Save last atomic positions and velocities 
      ! NOTE that this means "the positions predicted for next step"
      ! and "the current velocities"
      if (idyn.ne.6 .and. idyn.ne.7 .and. idyn.ne.8) then
        call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, 
     &             foundxv)
        if (lUseZmatrix .and. (idyn.eq.0))
     .    call iozm('write',ucell,vcell,xa,foundzm)
        call siesta_write_positions()


        ! Save atomic positions and velocities accumulatively and
        ! accumulate coor in Xmol file for animation 
        call write_md_record( istep )
      endif

   60 continue
! Output memory use at the end of this geometry step
      if (cml_p) then
         call cmlEndStep(mainXML)
      endif
!!      call printmemory( 6, 0 )
      call timer( 'SIMOVE', 2 )
      call timer( 'IterMD', 2 )
! End of one MD step - flush stdout
      if (ionode) call pxfflush(6)

#ifdef DEBUG
      call write_debug( 'POS siesta_move' )
#endif
!--------------------------------------------------------------------------- END
      END subroutine siesta_move

      END MODULE m_siesta_move
