c
      program cmpsig
c
c     Prints out the number of significant figures of agreement
c     between the data in two ASCII "SIG" files.
c
      implicit none
c
      double precision x1, x2
      character*50 s1, s2

      integer cmp
      external cmp

      integer iargc, nargs
      character*70 file1, file2, precision, form

      integer ndigits
c
c     Let's get the files from the command line:
c
      nargs = iargc()
      if (nargs .ne. 2) then
         write(0,*) 'Usage: cmpsig file1 file2'
         stop 
      endif
c
      call getarg(1,file1)
      call getarg(2,file2)
c
c      open files
c
      open(unit=1,file=file1,form='formatted',status='old')
      rewind(1)
      open(unit=2,file=file2,form='formatted',status='old')
      rewind(2)
c
c     Loop to process records
c
 1    continue
        read(1,9000,end=999) s1, x1
        read(2,9000,end=999) s2, x2
        write(6,'(i2,1x,a)') cmp(x1,x2), s1
        goto 1
c
 999   continue
 9000  format(a50,1x,g25.15)
       end

      integer function cmp(x1,x2)
c
c     Crude but good enough estimation of number of significant
c     figures which are common to x1 and x2
c
      real*8 x1, x2
      real*8 diff, relat
      diff = x2 - x1
      if (diff .eq. 0d0) then
         cmp = 15
         return
      endif
      relat = abs(diff)/x2
      cmp=int(-log10(abs(relat)))
      end






