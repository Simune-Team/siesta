      subroutine sig(unit,str,x)
c
c Writes entries in a SIG file
c A SIG file is a series of lines of the form:
c
c                 [50-character-long string, 1x, g25.15]
c
      implicit none
c
      integer unit
      double precision x
      character*(*) str

      character*50 string
      logical first

      string = str

      write(unit,'(a50,1x,g25.15)') string, x

      end

      subroutine sig_setup(io_sig)
      integer io_sig
      open(unit=io_sig,file='SIG',
     $     form='formatted',status='unknown')
      end

