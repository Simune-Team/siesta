! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_siesta_tddft

  implicit none
  private

  public :: siesta_tddft

contains

  subroutine siesta_tddft( istep )
#ifdef MPI
    use mpi_siesta
#endif
    use sys, only : die, bye
    use precision, only: dp
    use parallel,     only : IOnode, SIESTA_worker

    use siesta_cml
    use siesta_options

    use m_state_init
    use m_setup_hamiltonian
    use m_setup_H0
    use m_steps
    use sparse_matrices, only:H, Dscf, Escf, maxnh, numh, listhptr
    use m_eo,          only: eo
    use m_energies,    only: Etot           ! Total energy
    use atomlist,      only: no_s, no_l, no_u, indxuo
    use m_spin,        only: nspin
    use m_gamma
    use Kpoint_grid,    only: nkpnt, kpoint, kweight

    use m_compute_energies, only: compute_energies
    use m_state_analysis, only: state_analysis

    use alloc
    use m_initwf,         only: initwf
    use wavefunctions,    only: wavef_ms, iowavef, compute_tddm
    use m_evolve,         only: evolve
    use m_iotddft,        only: write_tddft
    use m_overfsm,        only: overfsm
    use m_final_H_f_stress, only: final_H_f_stress

    use m_mpi_utils, only: broadcast
    use fdf

    integer, intent(inout)  :: istep
    
    integer :: mod_tded, mod_md
    integer :: istpp, ispin
    real(dp) :: dt_tded, G2max
#ifdef DEBUG
    call write_debug( '    PRE siesta_tddft' )
#endif

#ifdef SIESTA__PEXSI
    ! Broadcast relevant things for program logic
    ! These were set in read_options, called only by "SIESTA_workers".
    call broadcast(nscf,comm=true_MPI_Comm_World)
#endif

    if ( SIESTA_worker )  then
       ! Initialization tasks for a given geometry
       call state_init( istep )
    end if

#ifdef SIESTA__PEXSI
    if (ionode) call memory_snapshot("after state_init")
#endif

    if ( fdf_get("Sonly",.false.) ) then
       if ( SIESTA_worker ) then
          call timer( 'all', 2 )
          call timer( 'all', 3 )
       end if
       call bye("S only")
    end if

    ! The current structure of the loop tries to reproduce the
    ! historical Siesta usage. It should be made more clear.
    ! Two changes:
    !
    ! -- The number of scf iterations performed is exactly
    !    equal to the number specified (i.e., the "forces"
    !    phase is not counted as a final scf step)
    !
    ! -- At the change to a TranSiesta GF run the variable "first"
    !    is implicitly reset to "true".

    ! This call computes the non-scf part of H and initializes the
    ! real-space grid structures.  It might be better to split the two,
    ! putting the grid initialization into state_init and moving the
    ! calculation of H_0 to the body of the loop, done if first=.true.  This
    ! would suit "analysis" runs in which nscf = 0
    if ( SIESTA_worker ) call setup_H0( G2max )

#ifdef SIESTA__PEXSI
    if (ionode) call memory_snapshot("after setup_H0")
#endif

    ! The first call to change basis only calculates S^+1/2
    call chgbasis(nspin, gamma, nkpnt, kpoint, kweight, no_u,istep)
    
    ! Read the initial wavefunctions.
    ! In future reading wavefunctions can be moved before
    ! changebasis and then the first DM can be computed within
    ! changesbasis. Moving up iowavef would require a call to 
    ! ms_scalapack_setup. Which is now implicit in changebasis. 

    if ( istep == 1 ) then
       allocate(wavef_ms(nkpnt,nspin))
       call iowavef('read',wavef_ms,no_u,nkpnt,nspin, istpp, &
            rstart_time)
       IF (IONode) THEN
       write(6,'(a)') 'Computing DM from initial KS wavefunctions'
       END IF 
       DO ispin=1,nspin
         call compute_tddm(ispin, Dscf)
       END DO
       istep = istpp + 1
    end if
    ! To save the TD-wavefunctions for a restart.
    mod_md = mod(istep, tdednwrite)

    if( mod_md == 0 .or. istep == fincoor ) then
       call iowavef('write',wavef_ms,no_u,nkpnt,nspin, istep,totime)
    end if

    do itded = 1 , ntded ! TDED loop
       
       dt_tded = dt/ntded
       if(IONode) then
          write(*,'(/a)') &
               '                     ************************       '
          write(*,*)'                 TDED Step    =', itded
          write(*,'(a)') &
               '                     ************************       '
       end if
       
       ! To save the TD-wavefunctions to restart
       mod_tded = mod(itded, tdednwrite)
       
       if ( mod_tded == 0 ) then
          call iowavef('write',wavef_ms,no_u,nkpnt,nspin, istep, totime)
       end if
       
       call setup_hamiltonian( itded )
       
       call evolve(no_s,nspin,nspin,no_l,maxnh,maxnh,no_u, &
            gamma, indxuo,nkpnt,kpoint,kweight, &
            Escf, no_u, dt_tded,istep,itded)
       
       ! The total simulation time mainly for plotting
       totime = (istep*dt - dt) + (itded*dt_tded-dt_tded)
       
       call compute_energies (itded)
       call write_tddft(totime, istep, itded, ntded, rstart_time, &
            etot, eo, no_u,nspin,nkpnt)
       
    end do ! TDED loop

    call final_H_f_stress(istep, 1, .false.)
    call state_analysis( istep )
    
#ifdef DEBUG
    call write_debug( '    POS siesta_tddft' )
#endif

  end subroutine siesta_tddft
  
end module m_siesta_tddft
