      subroutine iolwf( task, maxc, maxo, nbasis, nspin,
     .                  numc, listc, c, found )

C *******************************************************************
C Reads/writes localized wave functions from/to file
C Written by P.Ordejon and J.M.Soler. May 1997.
C ********* INPUT ***************************************************
C character task*3 : 'read' or 'write'
C integer   maxc   : First dimension of listd and c
C integer   maxo   : Second dimension of c (must be maxo.ge.nbasis)
C integer   nbasis : Number of atomic orbitals
C integer   nspin  : Number of spins (1 or 2)
C ********* INPUT OR OUTPUT (depending on task) *********************
C integer numc(nbasis)       : Control vector of c matrix
C                              (number of nonzero elements of each row)
C integer listc(maxc,nbasis) : Control vector of c matrix
C                              (list of nonzero elements of each row)
C real*8  c(maxc,maxo,nspin) : Density matrix
C ********* OUTPUT *************************************************
C logical found : Have LWF's been found in disk? (Only when task='read')
C ******************************************************************

      implicit          none
      character         task*(*), paste*33
      logical           found
      integer           maxc, maxo, nbasis, nspin
      integer           listc(maxc,maxo), numc(nbasis)
      double precision  c(maxc,maxo,nspin)
      external          chkdim, io_assign, io_close, paste, timer
      include          'fdf/fdfdefs.h'

c     Internal variables and arrays
      character fname*33, sname*30
      logical   exist1, exist2, exist3, frstme
      integer   im, is, unit1, unit2, m, nb, ncmax, ns
      save      frstme, fname
      data      frstme /.true./

*     call timer( 'iolwf', 1 )
      
c     Find file name
      if (frstme) then
        sname = fdf_string('SystemLabel','siesta')
        fname = paste(sname,'.LWF')
        frstme = .false.
      endif

      if (task.eq.'read' .or. task.eq.'READ') then
        inquire (file='SAVE.LWF',     exist=exist1)
        inquire (file='SAVE.ctrlLWF', exist=exist2)
        inquire (file=fname,          exist=exist3)

c       First look for old-format files
        if (exist1 .and. exist2) then
          write(6,'(/,a)') 'iolwf: Reading LWFs from file'

          call io_assign(unit1)
          call io_assign(unit2)
          open( unit1, file='SAVE.LWF', status='unknown')
          open( unit2, file='SAVE.ctrlLWF', status='unknown')
          rewind(unit1)
          rewind(unit2)
          ncmax = 0
          do m = 1,nbasis
            read(unit2,*) numc(m)
            ncmax = max( ncmax, numc(m) )
          enddo
          call chkdim( 'iolwf', 'maxc', maxc, ncmax, 1 )
          do m = 1,nbasis
            do im = 1,numc(m)
              read(unit2,*) listc(im,m)
            enddo
          enddo
          do is = 1,nspin
            do m = 1,nbasis
              do im = 1,numc(m)
                read(unit1,*) c(im,m,is)
              enddo
            enddo
          enddo
          call io_close(unit1)
          call io_close(unit2)
          found = .true.

c       Look now for new-format files
        elseif (exist3) then
          write(6,'(/,a)') 'iolwf: Reading LWFs from files'
          call io_assign(unit1)
          open( unit1, file=fname,
     .          form='unformatted', status='unknown' )
          rewind(unit1)
          read(unit1) nb, ns
          call chkdim( 'iolwf', 'nbasis', nbasis, nb, 0 )
          call chkdim( 'iolwf', 'nspin',  nspin,  ns, 0 )
          read(unit1) (numc(m), m=1,nbasis)
          ncmax = 0
          do m = 1,nbasis
            ncmax = max( ncmax, numc(m) )
          enddo
          call chkdim( 'iolwf', 'maxc', maxc, ncmax, 1 )
          read(unit1) ( (listc(im,m), im=1,numc(m)), m=1,nbasis )
          read(unit1) (( (c(im,m,is), im=1,numc(m)), m=1,nbasis ),
     .                                               is=1,nspin )
          call io_close(unit1)
          found = .true.

        else
          found = .false.
        endif

      elseif (task.eq.'write' .or. task.eq.'WRITE') then

        call io_assign(unit1)
        open( unit1, file=fname,
     .        form='unformatted', status='unknown' )
        rewind(unit1)
        write(unit1) nbasis, nspin
        write(unit1) (numc(m), m=1,nbasis)
        write(unit1) ( (listc(im,m), im=1,numc(m)), m=1,nbasis )
        write(unit1) (( (c(im,m,is), im=1,numc(m)), m=1,nbasis ),
     .                                              is=1,nspin )
        call io_close(unit1)

      else
        stop 'iolwf: incorrect task'
      endif

*     call timer( 'iolwf', 2 )
      end

