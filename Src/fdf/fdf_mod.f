      module fdf

c     Declarations for fdf procedures

      interface
         function fdf_integer(label,default)
         integer fdf_integer
         character(len=*), intent(in) :: label
         integer, intent(in) ::  default
         end function fdf_integer

         function fdf_single(label,default)
         real fdf_single
         character(len=*), intent(in) :: label
         real, intent(in) ::  default
         end function fdf_single

         function fdf_double(label,default)
         real*8 fdf_double
         character(len=*), intent(in) :: label
         real*8, intent(in) ::  default
         end function fdf_double

         function fdf_physical(label,default,unit)
         real*8 fdf_physical
         character(len=*), intent(in) :: label, unit
         real*8, intent(in) ::  default
         end function fdf_physical

         function fdf_boolean(label,default)
         logical fdf_boolean
         character(len=*), intent(in) :: label
         logical, intent(in) ::  default
         end function fdf_boolean

         function fdf_string(label,default)
         character(len=80) fdf_string
         character(len=*), intent(in) :: label
         character(len=*), intent(in) ::  default
         end function fdf_string

         function fdf_defined(label)
         logical fdf_defined
         character(len=*), intent(in) :: label
         end function fdf_defined

         function fdf_enabled()
         logical fdf_enabled
         end function fdf_enabled

         function fdf_block(label,unit)
         logical fdf_block
         character(len=*), intent(in) :: label
         integer, intent(out)  :: unit
         end function fdf_block

         function fdf_convfac(unit1,unit2)
         real*8 fdf_convfac
         character(len=*), intent(in) :: unit1, unit2
         end function fdf_convfac

         subroutine fdf_init(filein,fileout)
         character(len=*), intent(in) :: filein, fileout
         end subroutine fdf_init

      end interface
      end module fdf



