      module version_info

      implicit none
c
c     This file MUST be updated after every self-consistent commit,
c     and the PL ("patch level") number increased by one, unless the
c     modification involves raising a minor or major version number,
c     in which case the PL should be reset to zero.
c
c     A self-consistent commit is a group of changes that fix a bug
c     or implement a new feature, in such a way that the program can
c     be compiled (no loose ends left). An update to the ChangeLog file
c     should be an integral part of a commit (the PL number should be
c     included for reference.)
c
c     After it is done, this file should be commited.
c
      integer, dimension(3), save  :: num_version = (/1,0,27/)
      character(len=80), parameter ::  version_str =
     .  "SIESTA 1.0.27 -- [new sparse form] (Apr 12, 2000)"

      character(len=60), parameter :: arch =
     $  'SYS'
      character(len=60), parameter :: fflags =
     $  'FFLAGS'
      character(len=60), parameter :: mpilib =
     $  'MPILIB'


      end module version_info

      subroutine prversion
c
c     Simple routine to print the version string. Could be extended to
c     provide more information, if needed.
c
      use version_info
      implicit none
      write(6,'(a)') trim(version_str)
      write(6,'(2a)') 'Architecture  : ', trim(arch)
      write(6,'(2a)') 'Compiler flags: ', trim(fflags)
      if (len_trim(mpilib) .ne. 0) then
         write(6,'(a)') 'PARALLEL version'
      else
         write(6,'(a)') 'SERIAL version'
      endif

      end subroutine prversion

      subroutine get_version(v)
      use version_info
      implicit none
      integer, intent(out)  :: v(3)
      v = num_version
      end subroutine get_version


