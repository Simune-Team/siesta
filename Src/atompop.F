      module m_atompop
      public :: atompop
      private
      CONTAINS
!
      subroutine atompop(hirshpop, voropop,
     .                  nspin, iaorb, iphorb, 
     .                  no_u, na_u, na_s, isa, ntm, ucell,
     .                  maxnd, numd, listdptr, listd, 
     .                  Dscf )

C integer na_u           : Number of atoms in unit cell
C
C Calculates the Hirshfeld and Voronoi charges 
C Coded by P. Ordejon, August 2004
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C logical hirshpop      : Perform Hirshfeld population analysis
C logical voropop       : Perform Voronoi population analysis
C integer nspin         : Number of different spin polarisations
C                         nspin=1 => Unpolarized, nspin=2 => polarized
C                         nspin=4 => Noncollinear spin -NOT IMPLEMENTED-
C integer iaorb(no_s)   : Atom to which each orbital belongs
C integer iphorb(no_s)  : Orbital index (within atom) of each orbital
C integer indxuo        : Index of equivalent orbital in unit cell
C integer na_u           : Number of atoms in unit cell
C integer na_s            : Number of atoms in supercell
C integer isa(na_s)       : Species index of all atoms in supercell
C integer ntm(3)        : Number of mesh divisions of each cell
C                         vector, including subgrid.
C real*8 ucell(3,3)     : Unit cell vectors
C integer maxnd             : First dimension of listd and Dscf
C integer numd(no_l)         : Number of nonzero density-matrix
C                             elements for each matrix row
C integer listdptr(no_l)     : Pointer to start of rows of density-matrix
C integer listd(maxnd)      : Nonzero-density-matrix-element column
C                             indexes for each matrix row
C real*8  Dscf(maxnd,nspin): SCF density-matrix elements
C ----------------------------------------------------------------------
C Output : None
C ----------------------------------------------------------------------

C
C  Modules
C
      use precision
      use parallel
      use parallelsubs
      use atomlist, only: indxua, indxuo, no_s, no_l, Datm
      use atmfuncs, only: rcut, rcore, symfio, cnfigfio, labelfis, 
     .                    nofis, zvalfis
      use fdf
      use sys
      use mesh
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      logical, intent(in)     ::   hirshpop, voropop
      integer, intent(in)     ::
     .  maxnd, na_u, na_s, nspin, ntm(3), no_u,
     .  iaorb(no_s), iphorb(no_s), isa(na_s), 
     .  listd(*), listdptr(no_l), numd(no_l)

      real(dp), intent(In)   ::
     .  ucell(3,3), Dscf(maxnd,nspin)

C ----------------------------------------------------------------------
C Routines called internally:
C ----------------------------------------------------------------------
C        rhooda(...)    : Finds Harris electron density in the mesh
C        rhoofd(...)    : Finds SCF electron density in the mesh
C        volcel(...)    : Returns volume of cell
C ----------------------------------------------------------------------

C
C Internal variables
C
C  rhoscf, rhoatm  : defined as in dhscf
C  real*8 q(na) : atomic charges 
C
C
      integer
     .  ia, ispin, ntpl, ntp, 
     .  meshnsm(3), meshl(3), npl, ns, ispec

#ifdef MPI
      integer
     .  MPIerror
#endif

      real, dimension(:), allocatable, save ::
     .  rhoatm

      real, dimension(:,:), allocatable, save ::
     .  DRho

      real*8, dimension(:), allocatable, save ::
     .  q, qtmp

      real*8
     .  dvol, volume, volcel

      character
     .  atm_label*20

      external
     .  rhooda, rhoofd, rhoofdsp,
     .  memory, volcel


C ----------------------------------------------------------------------
C General initialisation
C ----------------------------------------------------------------------

      if (nspin .eq. 4) then
        write(6,*) 'atompop: Hirsfeld and Voronoi populations for'
        write(6,*) 'atompop: non-collinear spin not implemented...'
        return
      endif

      ns=0
      do ia = 1,na_u
        ns=max(ns,isa(ia))
      enddo


C ----------------------------------------------------------------------
C Allocate memory that depends on the mesh
C ----------------------------------------------------------------------

C Find system volume and volume per mesh point

      meshnsm(1) = ntm(1)/nsm
      meshnsm(2) = ntm(2)/nsm
      meshnsm(3) = ntm(3)/nsm

      call HowManyMeshPerNode(meshnsm,Node,Nodes,npl,meshl)

      volume = volcel( ucell )
      ntp = ntm(1) * ntm(2) * ntm(3)
      dvol = volume / ntp

C Number of extended mesh points locally
      ntpl = npl*nsp

C Allocate memory related to ntpl

      allocate(rhoatm(ntpl))
      call memory('A','S',ntpl,'hirsh')
      allocate(DRho(ntpl,nspin))
      call memory('A','S',ntpl*nspin,'hirsh')

C ----------------------------------------------------------------------

C Find Harris (sum of atomic) electron density
        call rhooda( no_s, npl, datm, rhoatm, 
     .               iaorb, iphorb, isa )

C ----------------------------------------------------------------------
C Find SCF electron density at mesh points. Store it in array DRho
C ----------------------------------------------------------------------

      call rhoofd(no_s, npl, maxnd, numd, listdptr, listd,
     .              nspin, Dscf, DRho, no_l, no_u, iaorb, iphorb, isa)


      allocate(q(na_u))
      call memory('A','D',na_u,'hirsh')
      allocate(qtmp(na_u))
      call memory('A','D',na_u,'hirsh')

C ----------------------------------------------------------------------
C Find Hirshfeld charges
C ----------------------------------------------------------------------

      if (hirshpop) then
        call hirshfeld( no_s, na_s, na_u, nspin, indxuo, indxua, 
     .                  npl, datm, rhoatm, DRho, iaorb, iphorb, isa, 
     .                  q )

#ifdef MPI
#ifdef NODAT
        call MPI_AllReduce(q(1),qtmp(1),na_u,
     .    MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
#else
        call MPI_AllReduce(q(1),qtmp(1),na_u,
     .    DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
#endif
        q(1:na_u) = qtmp(1:na_u)
#endif

        if (Node.eq.0) then
          write(6,*)
          write(6,"(a)")
     .          'atompop: Hirshfeld Net Atomic Populations:'
        endif

        do ispec =1, ns

          atm_label=labelfis(ispec)
          if (Node.eq.0) then
            write(6,'(/2a)')'Species: ', atm_label
            write(6,'(a4,a7)') 'Atom', 'Qatom'
          endif

          do ia = 1,na_u
            if (isa(ia).eq.ispec) then
              if (Node.eq.0) then
                write(6,'(i4,f7.3)')
     .          ia, zvalfis(ispec)-q(ia)*dvol
              endif
            endif
          enddo
        enddo
      endif

C ----------------------------------------------------------------------
C Find Voronoi charges
C ----------------------------------------------------------------------

      if (voropop) then
        call voronoi( no_s, na_s, na_u, nspin, indxua, 
     .                npl, rhoatm, DRho, iaorb,
     .                q )


#ifdef MPI
#ifdef NODAT
        call MPI_AllReduce(q(1),qtmp(1),na_u,
     .    MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
#else
        call MPI_AllReduce(q(1),qtmp(1),na_u,
     .    DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
#endif

        q(1:na_u) = qtmp(1:na_u)
#endif

        if (Node.eq.0) then
          write(6,*)
          write(6,"(a)")
     .          'atompop: Voronoi Net Atomic Populations:'
        endif

        do ispec =1, ns

          atm_label=labelfis(ispec)
          if (Node.eq.0) then
            write(6,'(/2a)')'Species: ', atm_label
            write(6,'(a4,a7)') 'Atom', 'Qatom'
          endif

          do ia = 1,na_u
            if (isa(ia).eq.ispec) then
              if (Node.eq.0) then
                write(6,'(i4,f7.3)')
     .            ia, -q(ia)*dvol
              endif
            endif
          enddo
        enddo
      endif

      call memory('D','D',size(qtmp),'hirsh')
      deallocate(qtmp)
      call memory('D','D',size(q),'hirsh')
      deallocate(q)

C ----------------------------------------------------------------------
C Free locally allocated memory
C ----------------------------------------------------------------------
      call memory('D','S',size(DRho),'hirsh')
      deallocate(DRho)
      call memory('D','S',size(rhoatm),'hirsh')
      deallocate(rhoatm)

      end subroutine atompop



      subroutine hirshfeld( no_s, na_s, na_u, nspin, indxuo, indxua, np,
     .                   Datm, rhoatm, rhoscf, iaorb, iphorb, isa, qh )
C ********************************************************************
C Finds the Hirshfeld atomic charges  
C Hirshfeld, Theo Chem Acta 44, 129 (1977)
C See Fonseca et al, J. Comp. Chem. 25, 189 (2003)
C 
C Written by P.Ordejon. August'04.
C *********************** InpUT **************************************
C integer no_s                : Number of basis orbitals
C integer na_s                : Number of atoms
C integer na_u               : Number of atoms in unit cell
C integer nspin             : Number of spins
C integer indxuo(no_s)        : Index of equivalent orbital in unit cell
C integer indxua(na_s)        : Index of equivalent atom in unit cell
C integer np                : Number of mesh points
C real*8 Datm(no_s)           : Occupations of basis orbitals in free atom
C integer iaorb(*)          : Pointer to atom to which orbital belongs
C integer iphorb(*)         : Orbital index within each atom
C integer isa(*)            : Species index of all atoms
C real rhoatm(nsp,np)       : Harris (sum of atoms) density at mesh points
C real rhoscf(nsp,np,nspin) : Selfconsistent charge density at mesh points
C *********************** OUTPUT **************************************
C real*8 qh(na_u)            : Hirshfeld charge on each orbital
C *********************************************************************
C
C  Modules
C
      use atmfuncs, only: rcut, phiatm
      use mesh,     only: nsp, dxa, xdop, xdsp
      use meshphi
C
      implicit none

      integer          no_s, np, nspin, na_s, na_u
      integer          indxuo(no_s), indxua(na_s),
     .                 iaorb(*), iphorb(*), isa(*)

      real             rhoatm(nsp,np),rhoscf(nsp,np,nspin)

      real*8           Datm(no_s), phip, qh(na_u)

      integer          i, ip, isp, iu, kn, io, iop, is, iphi, ia, ix, 
     .                 ispin, iua
      real*8           Ci, gradCi(3), r2o, r2sp, dxsp(3), Qi



C  Initialise Hirshfeld charges

      do iua = 1, na_u
        qh(iua) = 0.0
      enddo

C  Loop on mesh points
      do ip = 1,np

C  Loop on orbitals of mesh point
        do kn = 1+endpht(ip-1), endpht(ip)
          i = lstpht(kn)
          iu = indxuo(i)

C  Generate phi value and loop on subpoints
          iphi = iphorb(i)
          if (iphi .gt. no_s) stop 'error!!!'
          ia = iaorb(i)
          iua = indxua(ia)
          is = isa(ia)
          r2o = rcut(is,iphi)**2
          iop = listp2(kn)
          do isp = 1,nsp
            do ix = 1,3
              dxsp(ix) = xdop(ix,iop) + xdsp(ix,isp) - dxa(ix,ia)
            enddo
            r2sp = dxsp(1)**2 + dxsp(2)**2 + dxsp(3)**2
            if (r2sp.lt.r2o) then
              call phiatm(is,iphi,dxsp,phip,gradCi)
              Ci = phip
              Qi = Datm(iu) * Ci * Ci
              do ispin = 1, nspin
                qh(iua) = qh(iua) 
     .                + Qi * rhoscf(isp,ip,ispin)/rhoatm(isp,ip)
C			write(6,*) rhoscf(isp,ip,ispin),rhoatm(isp,ip)
              enddo
            endif
          enddo

        enddo

      enddo

      end subroutine hirshfeld


      subroutine voronoi( no_s, na_s, na_u, nspin, indxua, np,
     .                   rhoatm, rhoscf, iaorb, qv )
C ********************************************************************
C Finds the Voronoi atomic charges
C Bickelhaupt et al, Organometallics 15, 2923 (1996)
C See Fonseca et al, J. Comp. Chem. 25, 189 (2003)
C 
C Written by P.Ordejon. August'04.
C *********************** InpUT **************************************
C integer no_s                : Number of basis orbitals
C integer na_s                : Number of atoms
C integer na_u               : Number of atoms in unit cell
C integer nspin             : Number of spins
C integer indxua(na_s)        : Index of equivalent atom in unit cell
C integer np                : Number of mesh points
C integer iaorb(*)          : Pointer to atom to which orbital belongs
C real rhoatm(nsp,np)       : Harris (sum of atoms) density at mesh points
C real rhoscf(nsp,np,nspin) : Selfconsistent charge density at mesh points
C *********************** OUTPUT **************************************
C real*8 qv(na_u)              : Voronoi charge on each orbital
C *********************************************************************
C
C  Modules
C
      use atmfuncs, only: rcut, phiatm
      use mesh,     only: nsp, dxa, xdop, xdsp
      use meshphi
C
      implicit none

      integer          no_s, np, nspin, na_s, na_u
      integer          indxua(na_s),
     .                 iaorb(no_s)

      real             rhoatm(nsp,np),rhoscf(nsp,np,nspin)

      real*8           phip, qv(na_u)

      integer          i, ip, isp, iu, kn, io, iop, is, ia, ix, 
     .                 ispin, iua, ntimes(na_u), maxeq, neq, ieq

      real*8           Ci, gradCi(3), r2o, r2sp, dxsp(3), Qi, huge, 
     .                 rmin, qtot, tol

      logical          eq

      parameter (huge = 1.d30)
      parameter (maxeq = 20)
      parameter (tol = 1.0d-2)

      integer          iatom(maxeq)


C  Initialise Hirshfeld charges

      do iua = 1, na_u
        ntimes(iua) = 0
        qv(iua) = 0.0
      enddo


C  Loop on mesh points
      do ip = 1,np

C  Loop on mesh subpoints
        do isp = 1,nsp
          rmin = huge
          eq = .false.
          neq = 1
          do ieq = 1, maxeq
            iatom(ieq) = 0
          enddo

C  Loop on orbitals of mesh point, to check which atom is closest to subpoint
          do kn = 1+endpht(ip-1), endpht(ip)
            i = lstpht(kn)
            ia = iaorb(i)

C  iua is the index of the atom in unit cell
            iua = indxua(ia)

            iop = listp2(kn)
            do ix = 1,3
              dxsp(ix) = xdop(ix,iop) + xdsp(ix,isp) - dxa(ix,ia)
            enddo
            r2sp = dxsp(1)**2 + dxsp(2)**2 + dxsp(3)**2
            r2sp = dsqrt(r2sp)

C  If distance is equal to the previous minimum, determine if it is
C  another atom, or the same. If it is a different one, then
C  handle multiplicity of nearest atoms to grip point
C  Consider distances equal if within tolerance

            if (dabs(rmin-r2sp) .lt. tol) then
              eq = .false.
              do ieq = 1,neq
                if (iua .eq. iatom(ieq)) eq = .true.
              enddo
              if (.not. eq) then
                neq = neq+1
                if (neq .gt. maxeq) stop 'voronoi: increase maxeq'
                iatom(neq) = iua
              endif
              goto 100
            endif

            if (r2sp .lt. rmin) then
              neq = 1
              iatom(neq) = iua
              rmin = r2sp
            endif

100         continue

          enddo

C  Assign charge to atom iatom; if no atom was found, then the charge
C  is zero, so move to next subpoint

          if (iatom(1) .ne. 0) then
            do ieq = 1,neq
              ntimes(iatom(ieq)) = ntimes(iatom(ieq))+1
            enddo
            do ispin = 1, nspin
              do ieq = 1,neq
                qv(iatom(ieq)) = qv(iatom(ieq)) 
     .            + (rhoscf(isp,ip,ispin) - (rhoatm(isp,ip)/nspin))/neq
              enddo
            enddo
          else
C  Check that charge is actually zero
            qtot=0.0d0
            do ispin=1,nspin
              qtot = qtot + rhoscf(isp,ip,ispin) + rhoatm(isp,ip)
            enddo
            if (qtot .gt. 0.0d0) stop 'voronoi: Error in grid charge'
          endif

        enddo

      enddo

      end subroutine voronoi
      end module m_atompop
