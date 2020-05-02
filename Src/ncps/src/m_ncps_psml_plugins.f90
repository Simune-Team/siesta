module m_ncps_psml_plugins

  use m_psml
  public :: check_atom_grid

  private
  integer, parameter :: dp = selected_real_kind(10,100)
  
CONTAINS

  subroutine check_atom_grid(ps,atom_log_grid_ok,nrval,a,b)
    type(ps_t), intent(in) :: ps
    logical, intent(out)     :: atom_log_grid_ok
    integer, intent(out)     :: nrval
    real(dp), intent(out)    :: a
    real(dp), intent(out)    :: b

    type(ps_annotation_t)    :: grid_annotation
    integer                  :: status
    character(len=40)        :: strvalue
    
    grid_annotation = ps_GridAnnotation(ps)
    call get_annotation_value(grid_annotation,  &
         "type",strvalue,status)

    atom_log_grid_ok = .false.
    a = 0.0_dp
    b = 0.0_dp
    nrval = 0
    
    if (status == 0) then
       atom_log_grid_ok = ( (trim(strvalue) == "log-atom") .or. &
                            (trim(strvalue) == "sampled-log-atom") )
    endif

    if (.not. atom_log_grid_ok) RETURN

    call get_annotation_value(grid_annotation,  &
         "nrval",strvalue,status)
    if (status == 0) then
       read(strvalue,*) nrval
    else
       nrval = 0   ! For sampled log grid
    endif
    
    ! Note that a and b are interchanged in Siesta!
    call get_annotation_value(grid_annotation,  &
         "scale",strvalue,status)
    if (status /= 0) RETURN
    read(strvalue,*) b
    
    call get_annotation_value(grid_annotation,  &
         "step",strvalue,status)
    if (status /= 0) RETURN
    read(strvalue,*) a

    atom_log_grid_ok = .true.

  end subroutine check_atom_grid
end module m_ncps_psml_plugins
