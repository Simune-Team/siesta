      subroutine sig(unit,str,x)
c
c Writes entries in a SIG file
c A SIG file is a series of lines of the form:
c
c                 [50-character-long string, 1x, g25.15]
c
      use precision, only : dp

      implicit none
c
      integer  unit
      real(dp) x
      character*(*) str

      character(len=50) string

      string = str

      write(unit,'(a50,1x,g25.15)') string, x

      end

      subroutine sig_setup(io_sig)
      integer io_sig
      open(unit=io_sig,file='SIG',
     $     form='formatted',status='unknown')
      end

