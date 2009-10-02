! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine doping(cell,np,task,samestruc,rho,rhoatm,nspin)

c **********************************************************************
c Adds a background charge density to simulate doping.
c This routine calculates the net charge of the system, and
c adds a compensating background charge that makes the system 
c neutral. This background charge is constant at points of
c the mesh near the atoms, and zero at points far from the atoms.
c This simulates situations like doped slabs, where the extra
c electrons (holes) are compensated by oposite charges at the
c material (the ionized dopant impurities), but not at the vacuum.
c This serves to simulate properly doped systems in which there
c are large portions of vacuum, such as doped slabs.
c
c Written by P. Ordejon, July 2009
c ********* Input ******************************************************
c real*8  cell(3,3) : Unit cell vectors
c integer np        : Number of mesh points locally
c integer task      : 0 = add background;   1 = remove background added
c                     in a previous call
c logical samestruc : same mesh and atomic positions as in last call?
c real rhoatm(*)    : Sum of atomic charge densities
c integer nspin     : Number of spins (1 or 2)
c ********* Input and output *******************************************
c real    rho(*)    : Charge density, to which the background charge is
c                     added.
c ********* Units ******************************************************
c Charges in electros/borh**3
c Energies in Rydbergs
c ********* Behaviour **************************************************
c The routine must be first called with 'task=0', and later with 'task=1'
c - 'task=0':  the routine calculates the net charge of the system 
c   determines in which mesh points the compensating charge will be 
c   placed.  Then, it adds the background charge to these points.
c - 'task=1':  the background charge added in a previous call with task=0
c   is now removed.
c The points in which the background charge is removed are those for
c which the absolute value of the atomic charge rhoatm is larger than a 
c thresshold value, defined in this subroutine as thress=0.005  (empirically 
c found to work well). The user can exprement with this at will...
c **********************************************************************

C
C  Modules
C
      use parallel,     only : ionode
      use precision, only: dp, grid_p
      use mesh,         only : nsp
      use meshphi,      only : endpht
      use sys,          only: die
      use fdf


#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

      integer           np, task, nspin

      logical           samestruc

      real(grid_p)  ::  rho(nsp, np), rhoatm(nsp,np), charge
      real(dp)          cell(3,3)

C Internal parameters
C tol : tolerance for bulk components of the electric field
C thress:  thresshold charge for points where background is added
      real(dp), save :: tol = 1.0d-12, thress = 0.5d-2

C Internal variables
      logical           it
      logical,  save :: frstime = .true., found = .false.
      integer           i,ip, isp, nc
      integer, save, dimension(:), allocatable   ::  indip, indisp
      integer, save ::  ntloc, nt, nptot

#ifdef MPI
      integer           MPIerror, npl
#endif

      real(dp)          qtot,volcel,fact
      save              qtot
      external          volcel


      if (ionode) then
        if (frstime) then
          frstime = .false.
          found=fdf_boolean('SimulateDoping',.false.)
          if (found) then
            write(6,'(/,(a))')
     .        'doping: SimulateDoping = .true. in input file',
     .        'doping: Neutralizing background will be added to points',
     .        'doping: other than vacumm, to simulate doping'
            qtot = fdf_double('NetCharge',0.0_dp)
          endif
        endif
      endif

#ifdef MPI
      call MPI_Bcast(found,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(qtot,1,MPI_double_precision,0,
     .      MPI_Comm_World,MPIerror)

#endif


      if (.not. found) return

      if (.not. samestruc) then

        if(allocated(indip)) deallocate(indip)
        if(allocated(indisp)) deallocate(indisp)
        allocate(indip(np*nsp))
        allocate(indisp(np*nsp))


        ntloc = 0
C  Loop over grid points
        do ip = 1,np
C  Find number of nonzero orbitals at this point
          nc = endpht(ip) - endpht(ip-1)
C  Check if there is a nonzero orbital at that point
          if (nc .ne. 0) then
            do isp = 1,nsp
              it = .false.
C  Check if charge at that point is larger than thresshold
              charge = rhoatm(isp,ip)*nspin
              if(abs(charge).ge.thress) it = .true.
              if (it) then
                ntloc = ntloc+1
                indip(ntloc) = ip
                indisp(ntloc) = isp
              endif
            enddo
          endif
        enddo


#ifdef MPI
        call MPI_AllReduce(ntloc,nt,1,MPI_integer,MPI_sum,
     .                     MPI_Comm_World,MPIerror)
        call MPI_AllReduce(np,nptot,1,MPI_integer,MPI_sum,
     .                     MPI_Comm_World,MPIerror)
#else
        nt = ntloc
        nptot = np
#endif

      endif

      if (task .eq. 0) then
        fact = qtot*(nsp*nptot)/(nt*volcel(cell))
      else if (task .eq. 1) then
        fact = -qtot*(nsp*nptot)/(nt*volcel(cell))
      else 
        call die('ERROR: wrong task in doping.F')
      endif

C  Loop over grid points where charge has to be added/removed
      do i=1,ntloc
        ip = indip(i)
        isp = indisp(i)
        rho(isp,ip) = rho(isp,ip) + fact
      enddo


      end
