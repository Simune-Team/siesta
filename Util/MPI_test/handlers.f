      subroutine timer_mpi( prog, iOpt )

      implicit none
      character(len=*),intent(in):: prog ! Name of program to time
      integer,         intent(in):: iOpt ! Action option

         call timer(prog,iOpt)

      end subroutine timer_mpi


      subroutine timer(str,i)

      character(len=*), intent(in)  :: str
      integer,  intent(in)  :: i
      end subroutine timer

