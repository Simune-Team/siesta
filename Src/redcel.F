      subroutine redcel( alat, ucell, scell, mscell )

c *********************************************************************
c Reads the unit cell (and supercell, if present) vectors from the
c data file, or generates them automatically if absent there.
c Written by J.M.Soler
c ********* Output ****************************************************
c real*8 alat       : Lattice constant
c real*8 ucell(3,3) : Unit cell vectors ucell(ixyz,ivec)
c real*8 scell(3,3) : Supercell vectors (equal to ucell if no supercell)
c integer mscell(3,3): Supercell vectors in units of ucell, i.e.
c                      scell(ix,i) = Sum_j( ucell(ix,j)*mscell(j,i) )
c *********************************************************************
C
C  Modules
C
      use precision
      use fdf
#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

      integer           mscell(3,3)
      double precision  alat, scell(3,3), ucell(3,3)

C Internal variables
      integer           i, iunit, ix, j
      double precision  alp, alplp, betlp, blp, cell(3,3), clp,
     .                  gamlp, pi, xxx
#ifdef MPI
      integer           MPIerror
#endif

C Lattice constant
      alat = fdf_physical('LatticeConstant',0.d0,'Bohr')

C Lattice vectors

      if ( fdf_block('LatticeParameters',iunit) .and.
     .     fdf_block('LatticeVectors',iunit) ) then
         write(6,'(2a)')'redcel: ERROR: Lattice vectors doubly ',
     .     'specified: by LatticeVectors and by LatticeParameters.' 
#ifdef MPI
         write(6,'(a)')
     .     'redcel: ERROR: Double input for lattice vectors'
         call MPI_Abort(MPI_Comm_World,1,MPIerror)
#else
         stop 'redcel: ERROR: Double input for lattice vectors'
#endif
      endif

      if ( fdf_block('LatticeParameters',iunit) ) then
         read(iunit,*) alp, blp, clp, alplp, betlp, gamlp
*        write(6,'(a)')
*    .    'redcel: Lattice Parameters (units of Lattice Constant) ='
*        write(6,'(a,3f10.5,3f9.3)')
*    .    'redcel: ',alp,blp,clp,alplp,betlp,gamlp
         pi = acos(-1.d0)
         alplp = alplp * pi/180.d0
         betlp = betlp * pi/180.d0
         gamlp = gamlp * pi/180.d0
         cell(1,1) = alp
         cell(2,1) = 0.d0
         cell(3,1) = 0.d0
         cell(1,2) = blp * cos(gamlp)
         cell(2,2) = blp * sin(gamlp)
         cell(3,2) = 0.d0
         cell(1,3) = clp * cos(betlp)
         xxx = (cos(alplp) - cos(betlp)*cos(gamlp))/sin(gamlp)
         cell(2,3) = clp * xxx
         cell(3,3) = clp * sqrt(sin(betlp)*sin(betlp) - xxx*xxx)
      elseif ( fdf_block('LatticeVectors',iunit) ) then
        do i = 1,3
          read(iunit,*) (cell(j,i), j=1,3)
        enddo
      else
        do i = 1,3
          do j  = 1,3
            cell(i,j) = 0.d0
          enddo
          cell(i,i) = 1.d0
        enddo
      endif
C ...

C Multiply cell vectors by by lattice constant ........................
      do i = 1,3
        do ix = 1,3
          cell(ix,i) = alat * cell(ix,i)
        enddo
      enddo
C ..................

C Super cell ..........................................
      if ( fdf_block('SuperCell',iunit) ) then
        if ( alat .eq. 0.d0 ) then
          write(6,*) 'redcel: ERROR: LatticeConstant required',
     .                      ' to define SuperCell'
#ifdef MPI
          write(6,'(a)')
     .      'redcel: ERROR: LatticeConstant required'
          call MPI_Abort(MPI_Comm_World,1,MPIerror)
#else
          stop 'redcel: ERROR: LatticeConstant required'
#endif
        endif
        do i = 1,3
          read(iunit,*) (mscell(j,i),j=1,3)
        enddo
      else
        do i = 1,3
          do j = 1,3
            mscell(j,i) = 0
          enddo
          mscell(i,i) = 1
        enddo
      endif
      do i = 1,3
        do ix = 1,3
          ucell(ix,i) = cell(ix,i)
        enddo
      enddo
      do i = 1,3
        do ix = 1,3
          scell(ix,i) = ucell(ix,1) * mscell(1,i) +
     .                  ucell(ix,2) * mscell(2,i) +
     .                  ucell(ix,3) * mscell(3,i)
        enddo
      enddo

      end

