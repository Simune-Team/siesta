      subroutine kgridinit( cell, kscell, displ,
     .                  cutoff, nk )

c **********************************************************************
c Finds Monkhost-Pack k-point coordinates and weights.
c This version assumes no symmetry except time reversal, i.e.
c inversion in reciprocal space.
c Refs: H.J.Monkhorst and J.D.Pack, Phys Rev B 13, 5188 (1976)
c       J.Moreno and J.M.Soler, Phys Rev B 45, 13891 (1992)
c Written by J.M.Soler. July 1997.
c ***************** INPUT **********************************************
c real*8  cell(3,3)  : Unit cell vectors in real space cell(ixyz,ivec)
c integer kscell(3,3): Supercell reciprocal of k-grid unit cell
c                      scell(ix,i) = sum_j cell(ix,j)*kscell(j,i)
c real*8  displ(3)   : Grid origin in k-grid-vector coordinates:
c                      origin(ix) = sum_j gridk(ix,j)*displ(j)
c real*8  cutoff     : Minimum k-grid cutoff required
c                      Not used unless det(kscell)=0
c integer nk         : Dimension of arrays points and weight
c ***************** OUTPUT *********************************************
c integer kscell(3,3)  : Supercell reciprocal of k-grid unit cell
c                        Only if det(kscell)=0 on input
c real*8  displ(3)     : Grid origin in k-grid-vector coordinates:
c                        Only if det(kscell)=0 on input
c real*8  cutoff       : Actual k-grid cutoff
c integer nk           : Actual number of irreducible k-points
c ***************** UNITS **********************************************
c cutoff must be in the same units as cell
c ***************** BEHAVIOUR ******************************************
c - If det(kscell).ne.0, input cutoff is not used
c - If det(kscell).eq.0, kscell and displ are generated according to
c   input cutoff
c - If det(kscell).eq.0 .AND. cutoff.le.0.d0, they are readed
c   from the input fdf data file.
c   The relevant fdf labels are kgrid_cutoff and kgrid_Monkhorst_Pack.
c   If both are present, kgrid_Monkhorst_Pack has priority. If none is
c   present, the cutoff defect is cero, producing only the gamma point.
c   Examples of fdf data specifications:
c     kgrid_cutoff  50. Bohr
c     %block kgrid_Monkhorst_Pack  # Defines kscell and displ
c     4  0  0   0.50               # (kscell(i,1),i=1,3), displ(1)
c     0  4  0   0.50               # (kscell(i,2),i=1,3), displ(2)
c     0  0  4   0.50               # (kscell(i,3),i=1,3), displ(3)
c     %endblock kgrid_Monkhorst_Pack
c **********************************************************************
C
C  Modules
C
      use precision
      use fdf
#ifdef MPI
      use mpi
#endif

      implicit          none

      integer           kscell(3,3), nk
      double precision  cell(3,3), cutoff, displ(3)
      external          idiag, reclat
c ----------------------------------------------------------------------

c Internal variables
      integer           i, ir, iu, ix, j,
     .                  kdsc(3,3), maux(3,3,2), ml(3,3), mr(3,3),
     .                  ng(3), ni, nkr(3), nktot, Node
#ifdef MPI
      integer           MPIerror
#endif
      double precision  defcut, dkg(3), dkx(3), dscell(3,3),
     .                  gridk(3,3), gscell(3,3), huge, pi,  
     .                  scell(3,3), scmin(3,3), tiny, vmod
      logical           genlogic, spiral
      parameter (defcut = 0.d0)
      parameter (huge   = 1.d30)
      parameter (tiny   = 1.d-12)

c Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

c Find out if there is spiral arrangement of spins 
      if (Node.eq.0) then
        spiral = fdf_block('SpinSpiral',iu)
      endif
#ifdef MPI
      call MPI_Bcast(spiral,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

c Find total number of points (determinant of kscell)
      nktot = abs( kscell(1,1) * kscell(2,2) * kscell(3,3) +
     .             kscell(2,1) * kscell(3,2) * kscell(1,3) +
     .             kscell(3,1) * kscell(1,2) * kscell(2,3) -
     .             kscell(1,1) * kscell(3,2) * kscell(2,3) -
     .             kscell(2,1) * kscell(1,2) * kscell(3,3) -
     .             kscell(3,1) * kscell(2,2) * kscell(1,3) )

c Look for kscell or cutoff in input fdf file
      if ( nktot.eq.0 .and. cutoff.lt.tiny ) then
        if (Node.eq.0) then
          genlogic = fdf_block('kgrid_Monkhorst_Pack',iu)
        endif
#ifdef MPI
        call MPI_Bcast(genlogic,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
        if ( genlogic ) then
          if (Node.eq.0) then
            do i = 1,3
              read(iu,*) (kscell(j,i),j=1,3), displ(i)
            enddo
          endif
#ifdef MPI
          call MPI_Bcast(kscell(1,1),9,MPI_integer,0,MPI_Comm_World,
     .      MPIerror)
#ifdef NODAT
          call MPI_Bcast(displ,3,MPI_double_precision,0,MPI_Comm_World,
     .      MPIerror)
#else
          call MPI_Bcast(displ,3,DAT_double,0,MPI_Comm_World,
     .      MPIerror)
#endif
#endif
          nktot = abs( kscell(1,1) * kscell(2,2) * kscell(3,3) +
     .                 kscell(2,1) * kscell(3,2) * kscell(1,3) +
     .                 kscell(3,1) * kscell(1,2) * kscell(2,3) -
     .                 kscell(1,1) * kscell(3,2) * kscell(2,3) -
     .                 kscell(2,1) * kscell(1,2) * kscell(3,3) -
     .                 kscell(3,1) * kscell(2,2) * kscell(1,3) )
        else
c         The second argument is the default value
          if (Node.eq.0) then
            cutoff = fdf_physical('kgrid_cutoff',defcut,'Bohr')
          endif
#ifdef MPI
#ifdef NODAT
          call MPI_Bcast(cutoff,1,MPI_double_precision,0,
     .      MPI_Comm_World,MPIerror)
#else
          call MPI_Bcast(cutoff,1,DAT_double,0,
     .      MPI_Comm_World,MPIerror)
#endif
#endif
        endif
      endif

c Find kscell from required cutoff
      if (nktot .eq. 0) then
        nktot = 1
        do j = 1,3
          do i = 1,3
            kscell(j,i) = 0
          enddo
          vmod = sqrt( cell(1,j)**2 + cell(2,j)**2 + cell(3,j)**2 )
          kscell(j,j) = int(cutoff/(vmod/2.d0)) + 1
          if (mod(kscell(j,j),2) .eq. 0) displ(j) = 0.5d0
          nktot = nktot * kscell(j,j)
        enddo
      endif

c Find k-grid supercell
      do i = 1,3
        do ix = 1,3
          scell(ix,i) = cell(ix,1) * kscell(1,i) +
     .                  cell(ix,2) * kscell(2,i) +
     .                  cell(ix,3) * kscell(3,i)
        enddo
      enddo

c Find actual cutoff
      call minvec( scell, scmin )
      cutoff = huge
      do i = 1,3
        vmod = sqrt( scmin(1,i)**2 + scmin(2,i)**2 + scmin(3,i)**2 )
        cutoff = min( cutoff, vmod/2.d0 )
      enddo

c Find equivalent diagonal supercell
      call idiag( 3, kscell, kdsc, ml, mr, maux )
      do i = 1,3
        do ix = 1,3
          dscell(ix,i) = scell(ix,1) * mr(1,i) +
     .                   scell(ix,2) * mr(2,i) +
     .                   scell(ix,3) * mr(3,i)
        enddo
      enddo

c Find k-grid unit vectors
      call reclat( dscell, gridk, 1 )

c Find grid origin in cartesian coordinates
      call reclat( scell, gscell, 1 )
      do ix = 1,3
        dkx(ix) = gscell(ix,1) * displ(1) +
     .            gscell(ix,2) * displ(2) +
     .            gscell(ix,3) * displ(3)
      enddo

c Find grid origin in gridk coordinates
      pi = 4.d0 * atan(1.d0)
      do i = 1,3
        dkg(i) = ( dkx(1) * dscell(1,i) +
     .             dkx(2) * dscell(2,i) +
     .             dkx(3) * dscell(3,i) ) / (2*pi)
      enddo

c Some printout for debugging
*     write(6,'(/,a,/,(3f12.6,i6,f12.6))') 'kgrid: gridK,ng,dg =',
*    .  ((gridk(ix,i),ix=1,3),kdsc(i,i),dkg(i),i=1,3)

c Find total range of grid indexes
      do j = 1,3
        ng(j) = kdsc(j,j)
      enddo

c Find number of points with time-reversal (inversion) symmetry,
c (if not spiral) after reflection on each alternative plane
      if (spiral) then
c       If spiral, use all k points
        nk = nktot
      else
        do j = 1,3
          ni = ng(j)
          if (abs(dkg(j)) .lt. tiny) then
            ni = ng(j)/2 + 1
          elseif (abs(dkg(j)-0.5d0) .lt. tiny) then
            ni = (ng(j)-1)/2 + 1
          endif
          nkr(j) = ni * nktot / kdsc(j,j)
        enddo

c       Select reflection plane
        ir = 3
        if (nkr(2) .lt. nkr(ir)) ir = 2
        if (nkr(1) .lt. nkr(ir)) ir = 1

c       Find k points and weights
        nk = nkr(ir)
      endif

      return
      end

