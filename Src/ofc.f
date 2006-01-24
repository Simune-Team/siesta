      subroutine ofc(fa, dx, na)
C *******************************************************************
C Writes force constants matrix to file
C Input forces are in Ry/Bohr and input displacements are in Bohr.
C Force Constants written in file are in eV / Ang
C Written by P.Ordejon. August'98.
C Dynamic memory and save attribute for fres introduced by J.Gale
C Sept'99.
C ********* INPUT ***************************************************
C real*8 fa(3,na)             : atomic forces (in Ry / Bohr)
C real*8 dx                   : atomic displacements (in Bohr)
C integer na                  : number of atoms
C ********** BEHAVIOUR **********************************************
C On the first call (undisplaced coordinates), the forces should be 
C zero (relaxed structure).
C However, since the relaxation is usually not perfect, some
C residual forces are obtained. These residual forces are 
C substracted from the forces on other steps, to calculate the
C force constants matrix
C *******************************************************************

      use precision
      use fdf
      use files, only : slabel, label_length

      implicit          none

      integer           na
      real(dp)          dx, fa(3,na)
      external          io_assign, io_close, paste, timer,
     .                  memory

C Internal variables and arrays
      character(len=label_length+3), save :: fname
      character(len=label_length+3)       :: paste
      character(len=132)                  :: line
      logical,                       save :: frstme = .true.
      integer,                       save :: nwritten = 0
      integer                             :: i, ix, unit1, n
      real(dp)                            :: Ang, eV, rdummy
      real(dp), dimension(:,:), allocatable, save :: fres

C Allocate local array for storing residual forces
      if (.not.allocated(fres)) then
        allocate(fres(3,na))
        call memory('A','D',3*na,'ofc')
      endif

C Define conversion factors
      Ang = 1.0d0 / 0.529177d0
      eV  = 1.0d0 / 13.60580d0

C Find file name
      if (frstme) then
        fname = paste(slabel,'.FC')
      endif

      call io_assign(unit1)
      open( unit1, file=fname, status='unknown' )
      rewind(unit1)

      if (frstme) then
C Write header message if frstime
        write(unit1,'(a)') 'Force constants matrix'
C Set values of residual forces
        do i = 1,na
          do ix = 1,3
            fres(ix,i) = fa(ix,i)
          enddo
        enddo
        frstme = .false.
        call io_close(unit1)
        return
      endif

C Read file written so far to put pointer for write in the correct place
      read(unit1,'(a)') line
      do n = 1,nwritten
        do i=1,na
          read(unit1,'(3f15.7)') (rdummy, ix=1,3)
        enddo
      enddo

      do i=1,na
        write(unit1,'(3f15.7)') ((-fa(ix,i)+fres(ix,i))*
     .                              Ang**2/eV/dx, ix=1,3)
      enddo
      nwritten = nwritten + 1

      call io_close(unit1)

      return
      end
