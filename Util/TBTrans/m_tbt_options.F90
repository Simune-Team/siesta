! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_tbt_options

use precision, only : dp
use fdf
use files, only : slabel


implicit none
public

real(dp) Emin               !Calc. transmission on [Emin;Emax]
real(dp) Emax
real(dp) Volt
real(dp) GFeta              ! state broadening

integer ncontour
integer neigch            ! No. eigenchannels to calculate
integer NBUFATL,NBUFATR   ! No. buffer atoms, L/R

logical USEBULK           ! If true use bulk params. in L/R regions
logical sppol

character*33 hsfile       !name of HS-input file
integer isoat1,isoat2

integer NA1L,NA2L,nqL 
integer NA1R,NA2R,nqR

character*33 Lhsfile,Rhsfile

integer Lnucuse, Rnucuse

logical CalcIeig


! PARAMETERS

integer RNode   ! Root Node

contains

subroutine tbt_read_options()

use parallel, only : IOnode

character*33 paste
external paste

real(dp) eV
parameter ( eV = 1_dp / 13.60580_dp )

real(dp), parameter :: Emin_def = -2.0_dp*eV
real(dp), parameter :: Emax_def = 2.0d0*eV
real(dp), parameter :: Volt_def = 0._dp
integer, parameter :: ncontour_def = 100
real(dp), parameter :: GFeta_def = 0.00001_dp
logical, parameter :: UseBulk_def = .true.
integer, parameter :: neigch_def = 0
logical, parameter :: sppol_def = .false.
integer, parameter :: NBUFATL_def=0
integer, parameter :: NBUFATR_def=0
logical, parameter :: CalcIeig_def = .false.
character*33 :: hsfile_def

hsfile_def=paste(slabel,'.TSHS')



! Read from fdf

 Emin = fdf_get('TS.TBT.Emin',Emin_def,'Ry')
 Emax = fdf_get('TS.TBT.Emax',Emax_def,'Ry')
 Volt = fdf_get( 'TS.Voltage', Volt_def)

 ncontour = fdf_get( 'TS.TBT.NPoints', ncontour_def)
 GFeta = fdf_get('TS.TBT.Eta',GFeta_def,'Ry')

 USEBULK = fdf_get( 'TS.UseBulkInElectrodes', UseBulk_def)
 neigch = fdf_get( 'TS.TBT.NEigen', neigch_def)

 sppol = fdf_get( 'SpinPolarized', sppol_def)

 NBUFATL = fdf_get( 'TS.BufferAtomsLeft', NBUFATL_def)
 NBUFATR = fdf_get( 'TS.BufferAtomsRight', NBUFATR_def) 

 hsfile = fdf_get( 'TS.TBT.HSFile', hsfile_def)

 Lhsfile = fdf_get( 'TS.HSFileLeft', 'LeftELEC.TSHS')
 Rhsfile = fdf_get( 'TS.HSFileRight', 'RightELEC.TSHS')

 Lnucuse = fdf_get( 'TS.NumUsedAtomsLeft', 0)
 Rnucuse = fdf_get( 'TS.NumUsedAtomsRight', 0)

 isoat1 = fdf_get( 'TS.TBT.PDOSFrom', 0)
 isoat2 = fdf_get( 'TS.TBT.PDOSTo', 0)

nqL=1
nqR=1
NA1L=1
NA2L=1
NA1R=1
NA2R=1

 CalcIeig = fdf_get( 'TS.TBT.CalcIeig', CalcIeig_def)




! Print Values in output
if (IOnode) then
  write(*,'(a,f10.4,a)') &
  'tbtrans: Voltage                         = ',Volt
  write(*,'(a,2F14.6)') 'Emin,Emax = ',Emin/eV,Emax/eV
  write(*,'(a,I10)') 'Energy points in tbtrans = ',ncontour
  write(*,'(a,7L)') 'tbtrans: Spin-polarized = ', sppol
end if ! IOnode


end subroutine tbt_read_options

end module m_tbt_options
