      subroutine resize(na, nspin, cell, ncells)
C *********************************************************************
C Subroutine to read the data that influences the size of the main
C arrays for the SIESTA program. In particular the values that control
C the size of arrays that are filled in redata must be obtained here.
C
C     It uses the FDF (Flexible Data Fromat) package 
C     of J.M.Soler and A.Garcia
C
C Writen by J.Gale, September'99 based on redata.F
C **************************** OUTPUT *********************************
C integer na               : Number of atoms
C integer nspin            : Spin polarization
C real*8 cell(3,3)         : (Super) lattice vectors CELL(ixyz,ivector)
C                            (in Bohr)
C integer ncells           : Number of unit cells in supercell
C **********************************************************************

C
C  Modules
C
      use precision
      use fdf
#ifdef MPI
      use mpi
#endif

      implicit none

      integer
     .  na, nspin, ncells
 
      double precision
     .  cell(3,3)

C  Internal variables .................................................
      integer 
     .  na_default, i, ix, Node, mscell(3,3)

#ifdef MPI
      integer
     .  MPIerror
#endif

      double precision
     .  alat, ucell(3,3), volcel

      logical
     .  noncol, sppol

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

      na = 1
      nspin = 1

C Read Number of Atoms ...
      na_default = 0
      if (Node.eq.0) then
        na = fdf_integer('NumberOfAtoms',na_default)
        write(6,'(a,i5)') 
     .   'redata: Number of Atoms                  = ',na
      endif
C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(na,1,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
#endif
      if (na .le. 0) then
        if (Node.eq.0) then
          write(6,100) 
          write(6,101) 
          write(6,'(a)') 
     .   'redata: ERROR: Number of atoms must be larger than zero.'
          write(6,102)
        endif
        stop
      endif
C ..................

C Lattice constant and lattice vectors .........
      if (Node.eq.0) then
        call redcel( alat, ucell, cell, mscell )
        if (alat .ne. 0.d0) then
          write(6,'(a,f10.4,a)')
     .   'redata: Lattice Constant = ',alat,'  Bohr'
          write(6,'(a,/,(8x,3f12.6))')
     .     'redata: Lattice vectors (in units of Lattice Constant) =',
     .     ((cell(ix,i)/alat,ix=1,3),i=1,3)
        endif
      endif
#ifdef MPI
      call MPI_Bcast(alat,1,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
      call MPI_Bcast(ucell(1,1),9,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
      call MPI_Bcast(cell(1,1),9,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
      call MPI_Bcast(mscell(1,1),9,MPI_integer,0,MPI_Comm_World,
     .  MPIerror)
#endif
      if (volcel(ucell) .lt. 1.d-8) then
        ncells = 1
      else
        ncells = nint( volcel(cell) / volcel(ucell) )
        if (Node.eq.0) then
          if (ncells.gt.1) then
            write(6,'(a,/,(8x,3f12.6))')
     .        'redata: Total-cell (supercell) vectors (Bohr)', cell
            write(6,'(a,i6)') 'redata: Number of unit cells  =', ncells
            write(6,'(a,i6)') 'redata: Total number of atoms =', 
     .        na*ncells
          else
            write(6,'(a,/,(a,1x,3f12.6))')
     .       'redata: Lattice vectors (in Bohr)', 
     .       ('redata:',(cell(ix,i),ix=1,3),i=1,3)
          endif
        endif
      endif
C ..................

C Spin Polarization ...
      if (Node.eq.0) then
        sppol  = fdf_boolean('SpinPolarized',.false.)
        noncol = fdf_boolean('NonCollinearSpin',.false.)
      endif
#ifdef MPI
      call MPI_Bcast(sppol,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(noncol,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
      if (noncol) then
        nspin = 4
      elseif (sppol) then
        nspin = 2
      else 
        nspin = 1
      endif
      if (Node.eq.0) then
        write(6,'(a,4x,i1)') 
     . 'redata: Number of spin components        = ',nspin
      endif
C ...
 
100   format(/,'redata: ',71(1h*))
101   format('redata:                  INPUT ERROR')
102   format('redata: ',71(1h*))

      return
      end
