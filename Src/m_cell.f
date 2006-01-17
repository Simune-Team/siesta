      module m_cell
        use precision, only: dp
        real(dp), public, save  :: ucell(3,3)
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
