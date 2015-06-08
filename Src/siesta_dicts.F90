
module siesta_dicts

  use precision, only : dp
#ifdef FLOOK
  use dictionary
#endif

  implicit none

#ifdef FLOOK

  ! A dictionary for all the options
  type(dict) :: options

#ifdef TRANSIESTA
  ! A dictionary for all transiesta options
  type(dict) :: ts_options
#endif

  ! A dictionary for all variables
  type(dict) :: variables

  private :: dict_variable_add_b_0d
  private :: dict_variable_add_i_0d
  private :: dict_variable_add_d_0d
  private :: dict_variable_add_d_1d
  private :: dict_variable_add_d_2d
  interface dict_variable_add
     module procedure dict_variable_add_i_0d, dict_variable_add_b_0d
     module procedure dict_variable_add_d_0d
     module procedure dict_variable_add_d_1d, dict_variable_add_d_2d
  end interface dict_variable_add

contains

  subroutine dict_clean()
    call delete(options,dealloc=.false.)
#ifdef TRANSIESTA
    call delete(ts_options,dealloc=.false.)
#endif
    call delete(variables,dealloc=.false.)
  end subroutine dict_clean

  subroutine dict_populate()
    call dict_populate_options()
    call dict_populate_variables()
  end subroutine dict_populate

  subroutine dict_populate_options()
    use siesta_options

    ! We simply re-create the options, (note the 
    ! de-allocation by "nullification")
    call delete(options,dealloc=.false.)

    ! unluckily the dictionary does not
    ! implement a stringent way of doing characters
    ! by pointers (as we do not know their initial length).

    options = &
         ('DM.HistoryDepth'.kvp.DM_history_depth)

    ! Output options
    options = options // &
         ('Write.DenChar'.kvp.dumpcharge)
    options = options // &
         ('Write.MullikenPopDenChar'.kvp.mullipop)
    options = options // &
         ('Write.HirshfeldPop'.kvp.hirshpop)
    options = options // &
         ('Write.VoronoiPop'.kvp.voropop)

    ! SCF options
    options = options // &
         ('SCF.MixFirst'.kvp.mix_first_scf_step)
    options = options // &
         ('SCF.MinIterations'.kvp.min_nscf)
    options = options // &
         ('SCF.MaxIterations'.kvp.nscf)
    options = options // &
         ('SCF.MixHamiltonian'.kvp.mixH)
    options = options // &
         ('SCF.MixCharge'.kvp.mix_charge)
    options = options // &
         ('SCF.NumberPulay'.kvp.maxsav)
    options = options // &
         ('SCF.NumberBroyden'.kvp.broyden_maxit)
    options = options // &
         ('SCF.MixingWeight'.kvp.wmix)
    options = options // &
         ('SCF.NumberKick'.kvp.nkick)
    options = options // &
         ('SCF.KickMixingWeight'.kvp.wmixkick)
    options = options // &
         ('SCF.DM.Tolerance'.kvp.dDtol)
    options = options // &
         ('SCF.H.Tolerance'.kvp.dHtol)
    options = options // &
         ('SCF.MonitorForces'.kvp.monitor_forces_in_scf)
    options = options // &
         ('ElectronicTemperature'.kvp.temp)

    options = options // &
         ('MD.NumSteps'.kvp.nmove)
    options = options // &
         ('MD.MaxDispl'.kvp.dxmax)
    options = options // &
         ('MD.MaxForceTol'.kvp.ftol)
    options = options // &
         ('MD.FinalTimeStep'.kvp.ifinal)
    options = options // &
         ('MD.FC.Displ'.kvp.dx)
    options = options // &
         ('MD.FC.First'.kvp.ia1)
    options = options // &
         ('MD.FC.Last'.kvp.ia2)
    options = options // &
         ('MD.Temperature.Target'.kvp.tt)
    
  end subroutine dict_populate_options

  subroutine dict_populate_variables()

    use siesta_geom
    use m_forces
    use m_energies
    use atomlist

    ! We simply re-create the options, (note the 
    ! de-allocation by "nullification")
    call delete(variables,dealloc=.false.)

    ! Add geometries (lets do with this for now)
    variables = &
         ('geom.na_u'.kvp.na_u)
    variables = variables // &
         ('geom.cell'.kvp.ucell)
    variables = variables // &
         ('geom.vcell'.kvp.vcell)
    variables = variables // &
         ('geom.nsc'.kvp.nsc)
    variables = variables // &
         ('geom.xa'.kvp.xa)
    variables = variables // &
         ('geom.xa_last'.kvp.xa_last)
    variables = variables // &
         ('geom.va'.kvp.va)

    ! This is an abstraction made
    ! easy for the user.
    ! The forces that are used internally
    ! are actually the constrained ones.
    ! Hence any constrainst imposed will
    ! be visible to the user via this handle.
    variables = variables // &
         ('geom.fa'.kvp.cfa)
    ! This will let the user interact more
    ! freely by retrieval of all instances of the forces.
    variables = variables // &
         ('geom.fa_pristine'.kvp.fa)
    variables = variables // &
         ('geom.fa_constrained'.kvp.cfa)
    variables = variables // &
         ('geom.species'.kvp.isa)
    variables = variables // &
         ('geom.z'.kvp.iza)
    variables = variables // &
         ('geom.last_orbital'.kvp.lasto)
    variables = variables // &
         ('geom.mass'.kvp.amass)
    variables = variables // &
         ('geom.neutral_charge'.kvp.qa)

    ! Add energies
    variables = variables // &
         ('E.neutral_atom'.kvp.DEna)
    variables = variables // &
         ('E.electrostatic'.kvp.DUscf)
    variables = variables // &
         ('E.fermi'.kvp.Ef)
    variables = variables // &
         ('E.harris'.kvp.Eharrs)
    variables = variables // &
         ('E.kinetic'.kvp.Ekin)
    variables = variables // &
         ('E.total'.kvp.Etot)
    variables = variables // &
         ('E.exchange_correlation'.kvp.Exc)
    variables = variables // &
         ('E.free'.kvp.FreeE)

    ! Add the number of charges to the system
    variables = variables // &
         ('charge'.kvp.qtot)

  end subroutine dict_populate_variables

  subroutine dict_variable_add_b_0d(name,val)
    character(len=*), intent(in) :: name
    logical, intent(inout) :: val
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_b_0d
  subroutine dict_variable_add_i_0d(name,val)
    character(len=*), intent(in) :: name
    integer, intent(inout) :: val
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_i_0d
  subroutine dict_variable_add_d_0d(name,val)
    character(len=*), intent(in) :: name
    real(dp), intent(inout) :: val
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_d_0d
  subroutine dict_variable_add_d_1d(name,val)
    character(len=*), intent(in) :: name
    real(dp), intent(inout) :: val(:)
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_d_1d
  subroutine dict_variable_add_d_2d(name,val)
    character(len=*), intent(in) :: name
    real(dp), intent(inout) :: val(:,:)
    if ( name.in.variables ) call delete(variables,name,dealloc=.false.)
    variables = variables // (name.kvp.val)
  end subroutine dict_variable_add_d_2d

#endif

end module siesta_dicts
  
