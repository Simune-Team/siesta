! Also the mixing container
module m_mixing_scf

  use class_Fstack_dData1D
  use m_mixing, only: tMixer
  use m_mixing, only: mixing_reset, mixing_history_clear

  implicit none

  private
  save

  type(tMixer), pointer :: scf_mixs(:) => null()
  type(tMixer), pointer :: scf_mix => null()

  public :: scf_mixs, scf_mix

  public :: reset_mixing_scf
  public :: mixing_scf_converged
  public :: mixing_scf_history_clear

contains

  subroutine mixing_scf_converged( SCFconverged )

    use parallel, only: IONode

    logical, intent(inout) :: SCFconverged
    integer :: i

    ! Return if no convergence
    if ( .not. SCFconverged ) return

    if ( scf_mix%n_itt < 0 ) then

       ! this means that we skip to the 
       ! following algorithm
       scf_mix => scf_mix%next
       SCFconverged = .false.

       if ( allocated(scf_mix%stack) ) then

          do i = 1 , size(scf_mix%stack)

             ! delete all but one history
             ! This should be fine
             call reset(scf_mix%stack(i), -1)

          end do

       end if

       if ( IONode ) then
         write(*,'(a,a)') ':!: SCF cycle continuation mixer: ', &
              trim(scf_mix%name)
       end if

    end if

  end subroutine mixing_scf_converged


  subroutine reset_mixing_scf()

    nullify(scf_mix)
    call mixing_reset( scf_mixs )

  end subroutine reset_mixing_scf


  subroutine mixing_scf_history_clear( )
    
    call mixing_history_clear( scf_mixs )
    scf_mix => scf_mixs(1)
    
  end subroutine mixing_scf_history_clear

end module m_mixing_scf
