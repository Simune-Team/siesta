C $Id: iodm.f,v 1.4 1999/01/31 11:14:56 emilio Exp $

      subroutine iodm( task, maxnd, maxo, nbasis, nspin,
     .                 numd, listd, dm, found )
C *******************************************************************
C Reads/writes density matrix from/to file
C Written by P.Ordejon and J.M.Soler. May 1997.
C ********* INPUT ***************************************************
C character task*3 : 'read' or 'write'
C integer   maxnd  : First dimension of listd and dm
C integer   maxo   : Second dimension of dm (must be maxo.ge.nbasis)
C integer   nbasis : Number of atomic orbitals
C integer   nspin  : Number of spins (1 or 2)
C ********* INPUT OR OUTPUT (depending on task) *********************
C integer numd(nbasis)        : Control vector of DM matrix
C                               (number of nonzero elements of each row)
C integer listd(maxnd,nbasis) : Control vector of DM matrix
C                               (list of nonzero elements of each row)
C real*8  dm(maxnd,maxo,nspin): Density matrix
C ********* OUTPUT *************************************************
C logical found : Has DM been found in disk? (Only when task='read')
C ******************************************************************
      
      implicit          none
      character         task*(*), paste*33
      logical           found
      integer           maxnd, maxo, nbasis, nspin
      integer           listd(maxnd,maxo), numd(nbasis)
      double precision  dm(maxnd,maxo,nspin)
      external          chkdim, io_assign, io_close, paste, timer
      include          'fdf/fdfdefs.h'

c     Internal variables and arrays
      character fname*33, sname*30
      logical   exist1, exist2, exist3, frstme
      integer   im, is, unit1, unit2, m, nb, ndmax, ns
      save      frstme, fname
      data      frstme /.true./

*     call timer( 'iodm', 1 )

c     Find file name
      if (frstme) then
        sname = fdf_string('SystemLabel','siesta')
        fname = paste(sname,'.DM')
        frstme = .false.
      endif

      if (task.eq.'read' .or. task.eq.'READ') then
        inquire (file='SAVE.DM',     exist=exist1)
        inquire (file='SAVE.ctrlDM', exist=exist2)
        inquire (file=fname,         exist=exist3)

c       First look for old-format files
        if (exist1 .and. exist2) then
          write(6,'(/,a)') 'iodm: Reading Density Matrix from files'
          call io_assign(unit1)
          call io_assign(unit2)
          open( unit1, file='SAVE.DM', status='unknown')
          open( unit2, file='SAVE.ctrlDM', status='unknown')
          rewind(unit1)
          rewind(unit2)
          ndmax = 0
          do m = 1,nbasis
            read(unit2,*) numd(m)
            ndmax = max( ndmax, numd(m) )
          enddo
          call chkdim( 'iodm', 'maxnd', maxnd, ndmax, 1 )
          do m = 1,nbasis
            do im = 1,numd(m)
              read(unit2,*) listd(im,m)
            enddo
          enddo
          do is = 1,nspin
            do m = 1,nbasis
              do im = 1,numd(m)
                read(unit1,*) dm(im,m,is)
              enddo
            enddo
          enddo
          call io_close(unit1)
          call io_close(unit2)
          found = .true.

c       Look now for new-format files
        elseif (exist3) then
          write(6,'(/,a)') 'iodm: Reading Density Matrix from files'
          call io_assign(unit1)
          open( unit1, file=fname,
     .          form='unformatted', status='unknown' )
          rewind(unit1)
          read(unit1) nb, ns
          call chkdim( 'iodm', 'nbasis', nbasis, nb, 0 )
          call chkdim( 'iodm', 'nspin',  nspin,  ns, 0 )
          read(unit1) (numd(m), m=1,nbasis)
          ndmax = 0
          do m = 1,nbasis
            ndmax = max( ndmax, numd(m) )
          enddo
          call chkdim( 'iodm', 'maxnd', maxnd, ndmax, 1 )
          read(unit1) (  (listd(im,m), im=1,numd(m)), m=1,nbasis )
          read(unit1) (( (dm(im,m,is), im=1,numd(m)), m=1,nbasis ),
     .                                                is=1,nspin )
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
        write(unit1) (numd(m), m=1,nbasis)
        write(unit1) (  (listd(im,m), im=1,numd(m)), m=1,nbasis )
        write(unit1) (( (dm(im,m,is), im=1,numd(m)), m=1,nbasis ),
     .                                               is=1,nspin )
        call io_close(unit1)

      else
        stop 'iodm: incorrect task'
      endif

*     call timer( 'iodm', 2 )
      end




