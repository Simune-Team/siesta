
module m_ncps_froyen_ps_t

  implicit none

  private

  integer, parameter  :: dp = selected_real_kind(14)
      
  public :: froyen_ps_t
  public :: pseudo_init_constant

  type froyen_ps_t
        character(len=2)        :: name
        integer                 :: nr
        integer                 :: nrval
        real(dp)                :: zval
        real(dp)                :: gen_zval  ! Generation valence charge
        logical                 :: relativistic
        character(len=10)       :: correlation
        character(len=2)        :: icorr
        character(len=3)        :: irel
        character(len=4)        :: nicore
        real(dp)                :: a
        real(dp)                :: b
        character(len=10)       :: method(6)
        character(len=70)       :: text
        integer                 :: npotu
        integer                 :: npotd
        real(dp), pointer       :: r(:)        => null()
        real(dp), pointer       :: chcore(:)   => null()
        real(dp), pointer       :: chval(:)    => null()
        real(dp), pointer       :: vdown(:,:)  => null()
        real(dp), pointer       :: vup(:,:)    => null()
        integer, pointer        :: ldown(:)    => null()
        integer, pointer        :: lup(:)      => null()
     end type froyen_ps_t

      CONTAINS

      subroutine pseudo_init_constant(p)
      type(froyen_ps_t), intent(inout) :: p

      p%nr = 0
      p%nrval = 0
      p%zval = 0._dp
      p%gen_zval = 0._dp
      p%relativistic = .false.
      p%correlation = ' '
      p%icorr = ' '
      p%irel = ' '
      p%nicore = ' '
      p%a = 0._dp
      p%b = 0._dp
      p%method(:) = ' '
      p%text = ' '
      p%npotu = 0
      p%npotd = 0

      end subroutine pseudo_init_constant

   end module m_ncps_froyen_ps_t



