!============================================================
!= Sample program using the f90 FDF module : September 2007 =
!============================================================
!
!     Shows FDF capabilities..
!
PROGRAM SAMPLE
  USE fdf
  USE prec
  implicit none
!--------------------------------------------------------------- Local Variables
  integer, parameter         :: maxa = 100
  logical                    :: doit, debug
  character(20)              :: fname, axis, status
  character(2)               :: symbol(maxa)
  integer(sp)                :: i, ia, na, external_entry
  integer(sp)                :: isa(maxa)
  real(sp)                   :: wmix
  real(dp)                   :: cutoff, phonon_energy, factor
  real(dp)                   :: xa(3, maxa)
  type(block_fdf)            :: bfdf
  type(parsed_line), pointer :: pline


!------------------------------------------------------------------------- BEGIN

! Initialize
  call fdf_init('sample.fdf', 'sample.out')

! Handle/Use fdf structure
  if (fdf_defined('new-style')) write(6,*) 'New-style stuff'

  na = fdf_integer('NumberOfAtoms', 0)
  write(6,*) 'Examples: na =', na

  fname = fdf_string('NameOfFile', 'whatever')
  write(6,*) 'Name of file:', fname
 
  cutoff = fdf_physical('MeshCutoff', 8.d0, 'Ry')
  write(6,*) 'MeshCutOff:', cutoff

  phonon_energy = fdf_physical('phonon-energy', 0.01d0, 'eV')
  write(6,*) 'Phonom Energy:', phonon_energy

  i = fdf_integer('SomeInt', 34)
  write(6,*) '#elems:', i

  wmix = fdf_single('WmixValue', 0.55)
  write(6,*) 'WmixValue:', wmix

  factor = fdf_double('FactorValue', 1.d-10)
  write(6,*) 'Factor:', factor

  debug = fdf_boolean('Debug', .TRUE.)
  write(6,*) 'Debug:', debug

  doit = fdf_boolean('DoIt', .FALSE.)
  write(6,*) 'Doit:', doit

  if (fdf_block('AtomicCoordinatesAndAtomicSpecies', bfdf)) then
    ia = 1
    do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))
      do i= 1, 3
        xa(i,ia) = fdf_breals(pline, i)
      enddo
      isa(ia) = fdf_bintegers(pline, 1)
      ia = ia + 1
    enddo
  endif

  write(6,*) 'AtomicCoordinatesAndAtomicSpecies:'
  do ia= 1, na
    write(6,'(3F10.6,I5)') (xa(i,ia),i=1,3), isa(ia)
  enddo

  if (fdf_block('AtomicInfo', bfdf)) then
    ia = 1
    do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))
      do i= 1, 3
        xa(i,ia) = fdf_breals(pline, i)
      enddo
      ia = ia + 1
    enddo
  endif

  write(6,*) 'AtomicInfo:'
  do ia= 1, na
    write(6,'(3F10.6)') (xa(i,ia),i=1,3)
  enddo

  if (fdf_block('Other-Block', bfdf)) then

!   Forward reading
    ia = 1
    do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))
      symbol(ia) = fdf_bnames(pline, 1)
      do i= 1, na
        xa(i,ia) = fdf_breals(pline, i)
      enddo
      ia = ia + 1
    enddo

    write(6,*) 'Other-Block (Forward):'
    do ia= 1, na
      write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
    enddo

!   Backward reading
    ia = 1
    do while((fdf_bbackspace(bfdf, pline)) .and. (ia .le. na))
      symbol(ia) = fdf_bnames(pline, 1)
      do i= 1, na
        xa(i,ia) = fdf_breals(pline, i)
      enddo
      ia = ia + 1
    enddo

    write(6,*) 'Other-Block (Backward):'
    do ia= 1, na
      write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
    enddo

!   Forward reading
    ia = 1
    do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))
      symbol(ia) = fdf_bnames(pline, 1)
      do i= 1, na
        xa(i,ia) = fdf_breals(pline, i)
      enddo
      ia = ia + 1
    enddo

    write(6,*) 'Other-Block (Forward):'
    do ia= 1, na
      write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
    enddo

!   Forward reading with rewind
    call fdf_brewind(bfdf)
    ia = 1
    do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))
      symbol(ia) = fdf_bnames(pline, 1)
      do i= 1, na
        xa(i,ia) = fdf_breals(pline, i)
      enddo
      ia = ia + 1
    enddo

    write(6,*) 'Other-Block (Forward-with-rewind):'
    do ia= 1, na
      write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
    enddo
  endif

  external_entry = fdf_integer('externalentry', 60)
  write(6,*) 'ExternalEntry:', external_entry

  axis   = fdf_string('AxisXY', 'Cartesian')
  status = fdf_string('StatusXY', 'Enabled')
  write(6,*) 'Axis: ', TRIM(axis), ' | ', TRIM(status)

! Shutdown and deallocates fdf structure
  call fdf_shutdown()

  RETURN
!----------------------------------------------------------------------------END
END PROGRAM SAMPLE
