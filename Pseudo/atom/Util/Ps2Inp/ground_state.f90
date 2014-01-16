module m_ground_state

      integer, parameter, private :: dp = selected_real_kind(16,100)

      type, public ::  ground_state_t
          integer                   ::  lmax_valence
          integer                   ::  n(0:3)
          real(dp)                  ::  occupation(0:3)
          logical                   ::  occupied(0:4)   ! note 0..4
          real(dp)                  ::  z_valence
      end type ground_state_t

      public :: ground_state

      private

      character(len=1), dimension(0:3) :: sym = (/ "s", "p", "d", "f" /)

CONTAINS

      subroutine ground_state(z,gs)

        use periodic_table

      integer, intent(in)               ::  z
      type(ground_state_t), intent(out) :: gs
!
!     Determines ground state valence configuration from Z
!
      integer l, latm

      gs%z_valence = 0.d0
      do l=0,3
        gs%occupation(l)=0.0d0
      enddo

      call lmxofz(z,gs%lmax_valence,latm)
      call qvlofz(z,gs%occupation(:))
      do l=0,gs%lmax_valence
        gs%z_valence = gs%z_valence + gs%occupation(l)
      enddo
      call cnfig(z,gs%n(0:3))

      write(6,'(a,i2)',advance='no') 'Ground state valence configuration: '
      gs%occupied(4) = .false.         !! always
      do l=0,3
        gs%occupied(l) =  (gs%occupation(l).gt.0.1d0)
!        if (gs%occupied(l)) then
           write(6,'(2x,i1,a1,i2.2)',advance='no') &
                gs%n(l),sym(l),nint(gs%occupation(l))
!        endif
      enddo
      write(6,'(a)') ''

      end subroutine ground_state

    end module m_ground_state
