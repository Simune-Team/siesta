module m_verbosity

  ! Default verbosity
  integer, save :: verbosity = 5

contains

  subroutine init_verbosity(fdf_name,default,verb)
    use fdf, only : fdf_get
    character(len=*), intent(in) :: fdf_name
    integer, intent(in) :: default
    integer, intent(out), optional :: verb
    
    integer :: v

    v = fdf_get(fdf_name,default)
    if ( present(verb) ) then
       verb = v
    else
       verbosity = v
    end if

  end subroutine init_verbosity

end module m_verbosity

    
    
