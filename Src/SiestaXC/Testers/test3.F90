PROGRAM siestaXCtest3

  ! Compares the energy and potential calculated by atomXC and cellXC.
  ! J.M.Soler. Sept.2009

  ! Used module procedures
  USE siestaXC, only: atomXC
  USE siestaXC, only: cellXC
  USE siestaXC, only: setXC

  ! Used module parameters
  USE siestaXC, only: dp => siestaXC_std_p
  USE siestaXC, only: gp => siestaXC_grid_p

! Used MPI types
#ifdef MPI
  USE mpi_siesta, only: MPI_Double_Precision
  USE mpi_siesta, only: MPI_Max
  USE mpi_siesta, only: MPI_Sum
  USE mpi_siesta, only: MPI_Comm_World
#endif

  implicit none

  ! Tester parameters
  integer, parameter:: irel  =  0 ! Relativistic? 0=>no, 1=>yes
  integer, parameter:: nSpin =  2 ! Number of spin components
  integer, parameter:: nfTot = 10 ! Number of functionals
  integer, parameter:: jf1 =  10  ! First functional tested
  integer, parameter:: jf2 =  10  ! Last functional tested
  integer, parameter:: nr = 501   ! Number of radial points
  integer, parameter:: nx = 60    ! Number of grid points per lattice vector
  integer, parameter:: n1cut = 8  ! Cutoff parameter
  integer, parameter:: n2cut = 2  ! Cutoff parameter:
                                  !    fCut(r)=(1-(r/rMax)**n1cut)**n2cut
  real(dp),parameter:: dWidth = 2._dp ! Width of density distribution, in Bohr
  real(dp),parameter:: Qtot = 10._dp  ! Integral of density distribution
  real(dp),parameter:: spinPol= 2._dp ! Integral of densUp - densDown
  real(dp),parameter:: rMax = 12._dp  ! Cutoff radius, in Bohr
  real(dp),parameter:: rBuff = 3._dp  ! Radial buffer of zero density, in Bohr

  ! Functionals to be tested
  character(len=3):: func(nfTot) = (/'LDA',   'LDA',   'GGA',   'GGA',    &
                                     'GGA',   'GGA',   'GGA',   'GGA',    &
                                     'GGA',   'VDW'/)
  character(len=6):: auth(nfTot) = (/'PZ    ','PW92  ','PW91  ','PBE   ', &
                                     'RPBE  ','revPBE','LYP   ','WC    ', &
                                     'PBESOL','DRSLL '/)

  ! Tester variables and arrays
  integer :: cellMesh(3) = (/nx,nx,nx/)
  integer :: i1, i1max, i2, i2max, i3, i3max, ir, irmax, &
             lb1, lb2, lb3, myNode, nf, nNodes, ub1, ub2, ub3
  real(dp):: atomDens(nr,nSpin), atomEc, atomEx, atomDc, atomDx, &
             atomVxc(nr,nSpin), avgDeltaVxc, &
             cell(3,3), cellEc, cellEx, cellDc, cellDx, &
             d0, d0s(nSpin), deltaVxc(nSpin), dr, dx, Ecut, kCut, &
             latConst, maxDeltaVxc, pi, r, r2, recCell(3,3), rMesh(nr), &
             stress(3,3), sumDeltaVxc, Vxc(nSpin), &
             wc(nfTot), wr, wx(nfTot), x(3), x0(3)
  real(gp),allocatable:: cellDens(:,:,:,:), cellVxc(:,:,:,:)

#ifdef MPI
  ! MPI-related variables
  integer :: MPIerror
  integer :: nLarger, nxNode
#endif

  ! Initialize hybrid XC functional
  nf = jf2-jf1+1
  wx(jf1:jf2) = 1._dp / nf
  wc(jf1:jf2) = 1._dp / nf
  call setXC( nf, func(jf1:jf2), auth(jf1:jf2), wx(jf1:jf2), wc(jf1:jf2) )

  ! Find radial mesh points and gaussian density
  pi = acos(-1._dp)
  d0 = Qtot / (2*pi*dWidth**2)**1.5_dp    ! Total density at origin
  d0s(1) = d0 * (Qtot + spinPol/2) / Qtot ! Spin up density at origin
  d0s(2) = d0 * (Qtot - spinPol/2) / Qtot ! Spin down density at origin
  dr = rmax / (nr-1)                      ! Interval between radial points
  do ir = 1,nr
    rMesh(ir) = dr * (ir-1)               ! Radial point values
    atomDens(ir,:) = d0s(:) * exp(-rMesh(ir)**2/2/dWidth**2)
  end do

  ! Impose a smooth radial cutoff
  do ir = 1,nr
    atomDens(ir,:) = atomDens(ir,:) * ( 1 - (rMesh(ir)/rMax)**n1cut )**n2cut
  end do


  ! Find exchange and correlation energy and potential from radial density
  call atomXC( irel, nr, nr, rMesh, nSpin, atomDens, &
               atomEx, atomEc, atomDx, atomDc, atomVxc )

  ! Define fcc unit cell, such that a sphere of radius rMax+rBuff fits in it
  latConst = (rMax+rBuff) * 2*sqrt(2._dp)
  cell(:,1) = (/ 0.0_dp, 0.5_dp, 0.5_dp /)
  cell(:,2) = (/ 0.5_dp, 0.0_dp, 0.5_dp /)
  cell(:,3) = (/ 0.5_dp, 0.5_dp, 0.0_dp /)
  cell(:,:) = cell(:,:) * latConst

  ! Define reciprocal unit cell
  recCell(:,1) = (/-1.0_dp, 1.0_dp, 1.0_dp /)
  recCell(:,2) = (/ 1.0_dp,-1.0_dp, 1.0_dp /)
  recCell(:,3) = (/ 1.0_dp, 1.0_dp,-1.0_dp /)
  recCell(:,:) = recCell(:,:) * 2*pi/latConst
  kCut = cellMesh(1) * sqrt(sum(recCell(:,1)**2)) / 2  ! Max. wave vector
  Ecut = kCut**2                                       ! Mesh cutoff, in Ry
  dx = pi / kCut                                 ! Dist. between mesh planes

#ifdef MPI
  ! Initialize MPI and get myNode and nNodes
  call MPI_Init( MPIerror )
  call MPI_Comm_Rank( MPI_Comm_World, myNode, MPIerror )
  call MPI_Comm_Size( MPI_Comm_World, nNodes, MPIerror )
#else
  myNode = 0
  nNodes = 1
#endif

  ! Find the box of mesh points own by my processor
#ifdef MPI
  ! Do simplest thing: divide only along first axis
  nxNode = Nx / nNodes          ! Points per node along first vector
  nLarger = nx - nxNode*nNodes  ! Number of nodes with one more point
  if (myNode<nLarger) then      ! My node has nx+1 points
    lb1 = (nxNode+1)*(myNode-1)
    ub1 = (nxNode+1)*(myNode-1) - 1
  else                          ! My node has nx points
    lb1 = (nxNode+1)*nLarger + nxNode*(myNode-nLarger)
    ub1 = (nxNode+1)*nLarger + nxNode*(myNode-nLarger+1) - 1
  end if
#else
  ! All points belong to the only processor
  lb1 = 0
  ub1 = nx-1
#endif
  lb2 = 0
  lb3 = 0
  ub2 = nx-1
  ub3 = nx-1

  ! Allocate arrays for density and potential
  allocate( cellDens(lb1:ub1,lb2:ub2,lb3:ub3,nSpin), &
             cellVxc(lb1:ub1,lb2:ub2,lb3:ub3,nSpin) )

  ! Find density at mesh points
  x0(:) = sum(cell,2) / 2     ! Center of cell
  do i3 = lb3,ub3
  do i2 = lb2,ub2
  do i1 = lb1,ub1
    x(:) = i1*cell(:,1)/cellMesh(1) &   ! Mesh point position
         + i2*cell(:,2)/cellMesh(2) &
         + i3*cell(:,3)/cellMesh(3)
    r2 = sum((x-x0)**2)                 ! Distance to center of cell squared
    cellDens(i1,i2,i3,:) = d0s(:) * exp(-r2/2/dWidth**2)
  end do ! i1
  end do ! i2
  end do ! i3

  ! Find exchange and correlation energy and potential from density in cell
  call cellXC( irel, cell, cellMesh, lb1, ub1, lb2, ub2, lb3, ub3, nSpin, &
               cellDens, cellEx, cellEc, cellDx, cellDc, stress, cellVxc )

  ! Print parameters
  print'(a,2i6)', 'jf1, jf2 = ', jf1, jf2
  print'(a,3f12.6)', 'dr, dx, Ecut = ', dr, dx, Ecut
  print'(a,2f12.6)', 'rMax, rBuff = ', rMax, rBuff

  ! Compare energies
  if (myNode==0) then
    print'(a,3f12.6)', 'atomEx, cellEx, diff =', atomEx, cellEx, atomEx-cellEx
    print'(a,3f12.6)', 'atomEc, cellEc, diff =', atomEc, cellEc, atomEc-cellEc
    print'(a,3f12.6)', 'atomDx, cellDx, diff =', atomDx, cellDx, atomDx-cellDx
    print'(a,3f12.6)', 'atomDc, cellDc, diff =', atomDc, cellDc, atomDc-cellDc
  end if

  ! Compare potentials
  sumDeltaVxc = 0
  maxDeltaVxc = 0
  do i3 = lb3,ub3
  do i2 = lb2,ub2
  do i1 = lb1,ub1
    x(:) = cell(:,1)*i1/cellMesh(1) &
         + cell(:,2)*i2/cellMesh(2) &
         + cell(:,3)*i3/cellMesh(3)
    r = sqrt(sum((x-x0)**2))
    if (r>=rMax) cycle
    ! Simplest thing: a linear interpolation of atomVxc (requires large nr)
    ir = r/dr + 1
    wr = (r - ir*dr) / dr
    Vxc(:) = atomVxc(ir,:)*(1-wr) + atomVxc(ir+1,:)*wr
    deltaVxc(:) = abs(cellVxc(i1,i2,i3,:)-Vxc(:))
    if (maxval(DeltaVxc(:)) > maxDeltaVxc) then
      i1max = i1
      i2max = i2
      i3max = i3
      irmax = ir
    end if
    sumDeltaVxc = sumDeltaVxc + sum(deltaVxc(:)**2)
    maxDeltaVxc = max( maxDeltaVxc, maxval(deltaVxc(:)) )
  end do ! i1
  end do ! i2
  end do ! i3
#ifdef MPI
  ! Find sumDeltaVxc and maxDeltaVxc accross all processor nodes
  call MPI_AllReduce( sumDeltaVxc, sumDeltaVxc, 1, MPI_double_precision, &
                      MPI_Sum, MPI_Comm_World, MPIerror )
  call MPI_AllReduce( maxDeltaVxc, maxDeltaVxc, 1, MPI_double_precision, &
                      MPI_Max, MPI_Comm_World, MPIerror )
#endif
  avgDeltaVxc = sumDeltaVxc / nSpin / nx**3
  print'(a,2f15.9)', 'avgDeltaVxc, maxDeltaVxc = ', avgDeltaVxc, maxDeltaVxc
  print'(a,4i6)', 'i1max,i2max,i3max, irmax = ', i1max, i2max, i3max, irmax

END PROGRAM siestaXCtest3

