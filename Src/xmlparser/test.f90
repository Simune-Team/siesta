!Test the xml output of siesta 
!A xml file (called modified!) is compared against a  xml reference file.

program test_xml
use flib_dom
use compare_m
use f2kcli

implicit none

integer :: nargs, iostat

type(fnode), pointer     :: reference,modified 
character(len=1000) reference_file,modified_file

      nargs = command_argument_count()
      if (nargs /= 2)  STOP "Usage: test-xml file1 file2"

      call get_command_argument(1,value=reference_file,status=iostat)
      if (iostat /= 0) then
          STOP "Cannot get first argument"
      endif
      call get_command_argument(2,value=modified_file,status=iostat)
      if (iostat /= 0) then
         STOP "Cannot get second argument"
      endif


!Parse the files
reference => parsefile(trim(reference_file),verbose=.false.)
modified  => parsefile(trim(modified_file))

call compare(reference,modified)

end program 


