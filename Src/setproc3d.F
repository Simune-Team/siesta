      subroutine setProc3d( Nodes, nm, ProcX, ProcY, ProcZ )
      use precision, only: dp
      use parallel, only: Node
      implicit none
C     Passed arguments
      integer, intent(in)  :: Nodes, nm(3)
      integer, intent(out) :: ProcX, ProcY, ProcZ
C     Local variables
      integer              :: NN, DIV, NDI, LProcX, LProcY, LProcZ
      real(dp)             :: nmpp(3), GMIN, GMAX, LMAX, LMIN
      logical              :: notfound, notend
      integer, allocatable :: VDI(:), STATUS(:)
#ifdef MPI
      integer              :: istat
#endif
!------------------------------------------------------------------------- BEGIN
      ProcX = 1
      ProcY = 1
      ProcZ = 1

      if (Nodes.EQ.1) return

C     Count the number of divisors
      NN  = Nodes
      DIV = 2
      NDI = 0
      do while (NN.gt.1)
        if (MOD(NN,DIV).NE.0) then
          notfound = .true.
          do while (notfound)
            DIV = DIV + 1
            if (MOD(NN,DIV).EQ.0) notfound = .false.
          enddo
        endif
        NDI = NDI + 1
        NN = NN/DIV
      enddo

C     Allocate local memory
      allocate( VDI(NDI), STATUS(NDI) )
      call memory( 'A', 'I', 2*NDI, 'setProc3d' )

C     Make a list of divisors: VDI(1:NDI)
      NN  = Nodes
      DIV = 2
      NDI = 0
      do while (NN.gt.1)
        if (MOD(NN,DIV).NE.0) then
          notfound = .true.
          do while (notfound)
            DIV = DIV + 1
            if (MOD(NN,DIV).EQ.0) notfound = .false.
          enddo
        endif
        NDI = NDI + 1
        VDI(NDI) = DIV
        NN = NN/DIV
      enddo

      notend = .true.
      STATUS = 1
      GMIN   = max(max(nm(1),nm(2)),nm(3)) + 1.0
      GMAX   = 1
      do while( notend )
        LProcX = 1
        LProcY = 1
        LProcZ = 1
        call distriDIV( NDI, VDI, STATUS, LProcX, LProcY, LProcZ )

        nmpp = nm
        nmpp(1) = nmpp(1)/LProcX
        nmpp(2) = nmpp(2)/LProcY
        nmpp(3) = nmpp(3)/LProcZ
        LMAX     = max(max(nmpp(1),nmpp(2)),nmpp(3))
        if (LMAX.LT.GMIN) then
          GMIN  = LMAX
          ProcX = LProcX
          ProcY = LProcY
          ProcZ = LProcZ
          LMIN  = min(min(nmpp(1),nmpp(2)),nmpp(3))
          GMAX  = LMIN
        else if (LMAX.EQ.GMIN) then
          LMIN  = min(min(nmpp(1),nmpp(2)),nmpp(3))
          if (LMIN.GT.GMAX) then
            GMAX = LMIN
            ProcX = LProcX
            ProcY = LProcY
            ProcZ = LProcZ
          endif
        endif
        call updSTATUS( NDI, STATUS, notend )
      enddo

C     Free local memory
      call memory( 'D', 'I', 2*NDI, 'setProc3d' )
      deallocate( VDI, STATUS )

      nmpp    = nm
      nmpp(1) = nmpp(1)/ProcX
      nmpp(2) = nmpp(2)/ProcY
      nmpp(3) = nmpp(3)/ProcZ
      if (nmpp(1).lt.1.0 .or. nmpp(2).lt.1.0 .or. nmpp(3).lt.1.0 ) then
        if (Node.eq.0) then
          write(*,*) 'ERROR: Cannot distribute the mesh'
          write(*,*) 'ERROR: Bad number of processors'
        endif
#ifdef MPI
        call MPI_Finalize( istat )
#endif
        STOP        
      endif
!--------------------------------------------------------------------------- END
      end subroutine setProc3d

      recursive subroutine updSTATUS( NDI, STATUS, notend )
      implicit none
C     Input variables
      integer, intent(in)  :: NDI
C     Output variables
      integer, intent(out) :: STATUS(NDI)
      logical, intent(out) ::notend
      
      STATUS(1) = STATUS(1) + 1
      if (STATUS(1).EQ.4) then
        if (NDI.GT.1) then
          STATUS(1) = 1
          call updSTATUS( NDI-1, STATUS(2), notend )
        else
          notend = .false.
        endif
      endif
      end subroutine updSTATUS

      recursive subroutine distriDIV( NDI, VDI, STATUS,
     &                                ProcX, ProcY, ProcZ )
      implicit none
C     Input variables
      integer, intent(in)  :: NDI, VDI(NDI), STATUS(NDI)
C     Output variables
      integer, intent(out) :: ProcX, ProcY, ProcZ

      if (STATUS(1).EQ.1) then
        ProcX = ProcX*VDI(1)
      else if (STATUS(1).EQ.2) then
        ProcY = ProcY*VDI(1)
      else
        ProcZ = ProcZ*VDI(1)
      endif
      if (NDI.GT.1) then
        call distriDIV( NDI-1, VDI(2), STATUS(2), ProcX, ProcY, ProcZ )
      endif
      end subroutine distriDIV
