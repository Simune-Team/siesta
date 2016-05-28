!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been improved or fully created by:
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!
module m_ts_voltage

! Module for containing various ways of implementing the voltage in a
! TranSIESTA run.
! As a standard method TranSIESTA applies the voltage ramp across the
! entire unit cell.
! As a future progress the voltage ramp might be situated in the contact region
! only. This makes more physical sense as the voltage drop does not occur
! in the electrode regions.
!
! Created and copyrighted by: Nick Papior Andersen, 2012
! The use of this program is allowed for not-for-profit research only.
! Copy or disemination of all or part of this package is not
! permitted without prior and explicit authorization by the author.
!

  implicit none
  
  private

  ! The idea is to have sub routines in this module to do
  ! various voltage layouts
  public :: ts_voltage
  public :: print_ts_voltage

contains

  subroutine ts_voltage(cell,meshG,nsm,v)
    use precision,    only : dp, grid_p
! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: cell(3,3)
    integer,       intent(in) :: meshG(3), nsm
! ***********************
! * OUTPUT variables    *
! ***********************
    real(grid_p), intent(inout) :: v(*)
    call timer('ts_volt',1)

    ! Here one should implement different calls...
    call ts_ramp_cell(cell,meshG,nsm,v)

    call timer('ts_volt',2)
  end subroutine ts_voltage

  subroutine ts_ramp_cell(cell, meshG, nsm, v)
    use precision,    only : dp, grid_p
    use fdf,          only : fdf_physical
    use parallel,     only : IONode
    use units,        only : eV
    use mesh,         only : meshLim
    use m_ts_options, only : VoltFDF, VoltL, VoltR
! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: cell(3,3)
    integer,       intent(in) :: meshG(3), nsm
! ***********************
! * OUTPUT variables    *
! ***********************
    real(grid_p), intent(inout) :: v(*)

! ***********************
! * LOCAL variables     *
! ***********************
    integer                   :: ix, i1, i2, i3, i10, i20, i30, imesh
    integer                   :: meshl(3), Xoffset, Yoffset, Zoffset
    real(dp)                  :: Lvc, dLvc
    real(dp)                  :: f(3)
    integer,        parameter :: ivc=3              ! The voltage-direction
    real(dp)                  :: dot
    external                  :: dot

    Lvc = sqrt(dot(cell(1,ivc),cell(1,ivc),3))

    dLvc = Lvc/max( meshG(ivc), 1 ) !

! field in [0;Lvc]: v = e*x = f*index
    do ix = 1,3
       f(ix) = 0d0         !init
    enddo
    f(ivc) = -VoltFDF*dLvc/Lvc

! Find quantities in mesh coordinates
    meshl(1) = (meshLim(2,1) - meshLim(1,1)+1)*nsm
    meshl(2) = (meshLim(2,2) - meshLim(1,2)+1)*nsm
    meshl(3) = (meshLim(2,3) - meshLim(1,3)+1)*nsm

! Calculate starting point for grid
    Xoffset = (meshLim(1,1)-1)*nsm
    Yoffset = (meshLim(1,2)-1)*nsm
    Zoffset = (meshLim(1,3)-1)*nsm

! Add the electric field potential to the input potential
    imesh = 0
    i30 = Zoffset - 1
    do i3 = 0,meshl(3)-1
       i30 = i30 + 1
       i20 = Yoffset - 1
       do i2 = 0,meshl(2)-1
          i20 = i20 + 1
          i10 = Xoffset - 1
          do i1 = 0,meshl(1)-1
             i10   = i10 + 1
             imesh = imesh + 1
             v(imesh) = v(imesh) + VoltL + f(1)*i10 + f(2)*i20 + f(3)*i30
          enddo
       enddo
    enddo

  end subroutine ts_ramp_cell

! Print out the voltage direction dependent on the cell parameters.

  subroutine print_ts_voltage(cell)
    use precision,    only : dp
    use parallel,     only : IONode
    use units,        only : eV
    use m_ts_options, only : VoltFDF, isVolt
! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: cell(3,3)

! ***********************
! * LOCAL variables     *
! ***********************
    integer                   :: i
    real(dp)                  :: Lvc, vcdir(3)
    integer,        parameter :: ivc=3              ! The voltage-direction
    real(dp)                  :: dot
    external                  :: dot
    
    Lvc = sqrt(dot(cell(1,ivc),cell(1,ivc),3))
    do i=1,3
       vcdir(i) = cell(i,ivc)/Lvc
    end do

    if(IONode) then
       write(6,*)
       write(6,'(a,f6.3,1x,a)')'ts_voltage: Bias ', VoltFDF/eV,'V'
       write(6,'(a,3(f6.3,a))')'ts_voltage: In unit cell direction = {', &
            vcdir(1),',',vcdir(2),',',vcdir(3),'}'
    end if
  end subroutine print_ts_voltage

end module m_ts_voltage
