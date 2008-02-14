      module debugMpi
      logical, save     :: debugON = .false.
      integer, save     :: mpiUnit  = 23
      end module

      subroutine debugMpiOn( )
      use parallel
      use debugMpi
      implicit none
      character*128         :: fname, cnode
!------------------------------------------------------------------------- BEGIN
      write(*,*) 'Que esta pasando RGT'
      if (.not.debugON) then
        debugON = .true.
        if (Node.lt.10) then
          write(cnode,'(I1)') Node
        else if (Node.lt.100) then
          write(cnode,'(I2)') Node
        else if (Node.lt.1000) then
          write(cnode,'(I3)') Node
        else if (Node.lt.10000) then
          write(cnode,'(I4)') Node
        else
          write(cnode,'(I5)') Node
        endif
        fname = TRIM('debug.'//cnode)
        open( mpiUnit, file=fname, status='unknown' )
      endif
!--------------------------------------------------------------------------- END
      end subroutine debugMpiOn

      subroutine write_debug( mesg )
      use debugMpi
      implicit none
      character*(*)        :: mesg
!------------------------------------------------------------------------- BEGIN
      if (debugON) then
        write(mpiUnit,*) mesg
        call pxfflush(mpiUnit)
      endif
!--------------------------------------------------------------------------- END
      end subroutine write_debug
