! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module m_cell
        use precision, only: dp
        use siesta_geom, only: ucell
        implicit none
        real(dp), public, save  :: celli(3,3)

        public :: cart2frac, frac2cart
        private

        CONTAINS

        subroutine cart2frac(cart,frac)
        real(dp), intent(in)  :: cart(3)
        real(dp), intent(out)  :: frac(3)

        frac =  matmul(transpose(celli),cart)
        end subroutine cart2frac

        subroutine frac2cart(frac,cart)
        real(dp), intent(in)  :: frac(3)
        real(dp), intent(out)  :: cart(3)

        cart =  matmul(ucell,frac)
        end subroutine frac2cart

      end module m_cell
!---------------------------------------------------
