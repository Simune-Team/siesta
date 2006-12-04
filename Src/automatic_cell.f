      subroutine automatic_cell(ucell,scell,na_u,xa,isa,charnet)

      use precision, only: dp
      use atmfuncs,  only: rcut
      use parallel,  only: IONode
      use units,     only: Ang

      implicit none

      real(dp), dimension(3,3), intent(inout) :: ucell
      real(dp), dimension(3,3), intent(inout) :: scell
      integer, intent(in)                     :: na_u
      real(dp), dimension(3,*), intent(in)    :: xa
      integer, dimension(*), intent(in)       :: isa
      real(dp), intent(in)                    :: charnet


      integer  :: ix, ia, is, iv
      real(dp) :: rc, xmin, xmax

        ucell(1:3,1:3) = 0.0_dp
        scell(1:3,1:3) = 0.0_dp
        do ix = 1,3
          xmin =  huge(1._dp)
          xmax = -xmin
          do ia = 1,na_u
            is = isa(ia)
            rc = rcut(is,0)
            xmin = min( xmin, xa(ix,ia)-rc )
            xmax = max( xmax, xa(ix,ia)+rc )
          enddo
C         Use a 10% margin for atomic movements
          ucell(ix,ix) = 1.10_dp * (xmax - xmin)
          scell(ix,ix) = ucell(ix,ix)
        enddo
C build cubic cell if system is charged
        if (charnet .ne. 0.0_dp) then
          xmax = -huge(1._dp)
          do ix = 1,3
            if (ucell(ix,ix) .gt. xmax) xmax = ucell(ix,ix)
          enddo
          do ix = 1,3
            ucell(ix,ix) = xmax
            scell(ix,ix) = xmax
          enddo
        endif
C
        if (IOnode) then
          write(6,'(/,a,3(/,a,3f12.6))')
     .      'siesta: Automatic unit cell vectors (Ang):',
     .      ('siesta:', (ucell(ix,iv)/Ang,ix=1,3), iv =1,3)
        endif

        end subroutine automatic_cell

