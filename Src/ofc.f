C $Id: ofc.f,v 1.1 1999/02/28 19:08:35 ordejon Exp $

      subroutine ofc(fa, dx, na)
C *******************************************************************
C Writes force constants matrix to file
C Input forces are in Ry/Bohr and input displacements are in Bohr.
C Force Constants written in file are in eV / Ang
C Written by P.Ordejon. August'98.
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
      
      implicit          none
      integer           na
      double precision  dx, fa(3,na)
      external          chkdim, io_assign, io_close, paste, timer
      include          'fdf/fdfdefs.h'

c     Internal variables and arrays
      integer maxa
      parameter (maxa = 2000)
      character fname*33, sname*30, line*132, paste*33
      logical   frstme
      integer   i, ix, unit1
      double precision Ang, eV, fres(3,maxa)
      save      frstme, fname
      data      frstme /.true./

c     check dimensions
      call chkdim('ofs','maxa',maxa,na,1)

c     Define conversion factors
      Ang = 1.d0 / 0.529177d0
      eV  = 1.d0 / 13.60580d0

c     Find file name
      if (frstme) then
        sname = fdf_string('SystemLabel','siesta')
        fname = paste(sname,'.FC')
      endif

      call io_assign(unit1)
      open( unit1, file=fname, status='unknown' )
      rewind(unit1)

      if (frstme) then
c     Write header message if frstime
        write(unit1,'(a)') 'Force constants matrix'
c     Set values of residual forces
        do i=1,na
          do ix=1,3
            fres(ix,i) = fa(ix,i)
          enddo
        enddo
        frstme = .false.
        call io_close(unit1)
        return
      endif

c     Goto end of file if not frstime
10    read(unit1,end=100,err=100,fmt='(a)') line
      goto 10

100   do i=1,na
        write(unit1,'(3f15.7)') ((-fa(ix,i)+fres(ix,i))*
     .                              Ang**2/eV/dx, ix=1,3)
      enddo

      call io_close(unit1)

      return
      end




