! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_partial_charges

  implicit none
  
  public :: compute_partial_charges
  
  logical, public :: want_partial_charges = .false.

  private

contains

  subroutine compute_partial_charges(DRho,rhoatm, &
      nspin, iaorb, iphorb, &
      isa, nmpl,dvol)

! Calculates the Hirshfeld and Voronoi charges 
! Coded by P. Ordejon, August 2004
! Integrated into dhscf workflow by A. Garcia, October 2012
! Added spin-resolved output+NC+SOC, Nick Papior, May 2020
!
! ----------------------------------------------------------------------
! Input :
! ----------------------------------------------------------------------
! integer nspin         : Number of different spin polarisations
!                         nspin=1 => Unpolarized, nspin=2 => polarized
!                         nspin=4 => Noncollinear spin -NOT IMPLEMENTED-
! integer iaorb(no_s)   : Atom to which each orbital belongs
! integer iphorb(no_s)  : Orbital index (within atom) of each orbital
! integer indxuo        : Index of equivalent orbital in unit cell
! integer isa(na_s)     : Species index of all atoms in supercell
! real(grid_p) Drho(:,:): Electron charge density (in cluster form)
! real(grid_p) rhoatm(:): Superposition of atomic Electron charge densities
! ----------------------------------------------------------------------
! Output : None
! ----------------------------------------------------------------------
    
    use precision, only:  dp, grid_p
    use parallel
    use atomlist, only: indxua, indxuo, no_s, no_u, Datm
    use siesta_geom, only: na_u, na_s
    use atmfuncs, only: zvalfis
    use sys
    use siesta_options, only: hirshpop, voropop
#ifdef MPI
    use mpi_siesta
#endif

    integer, intent(in) :: nspin, iaorb(no_s), iphorb(no_s), isa(na_s), nmpl

    real(dp), intent(in) :: dvol
    real(grid_p), intent(in) :: rhoatm(:)
    real(grid_p), intent(in) :: DRho(:,:)

! Internal variables
    integer :: ia, is, nsd

    ! Atomic charges
    real(dp), dimension(:,:), allocatable :: q

#ifdef MPI
    integer :: MPIerror
    real(dp), dimension(:,:), allocatable :: qtmp
#endif

    character(len=20) ::  atm_label

    external :: memory

! ----------------------------------------------------------------------
! General initialisation
! ----------------------------------------------------------------------
    allocate(q(nspin,na_u))
    call memory('A','D',na_u,'hirsh')
#ifdef MPI
    allocate(qtmp(nspin,na_u))
    call memory('A','D',na_u,'hirsh')
#endif

! ----------------------------------------------------------------------
! Find Hirshfeld charges
! ----------------------------------------------------------------------
    if (hirshpop) then
      call hirshfeld(no_s, na_s, na_u, nspin, indxuo, indxua, &
          nmpl, datm, rhoatm, DRho, iaorb, iphorb, isa, dvol, q)

#ifdef MPI
      call MPI_AllReduce(q(1,1),qtmp(1,1),nspin*na_u, &
          MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
      q(:,:) = qtmp(:,:)
#endif

      call write_partial_charges(nspin, na_u, q, method=1)

    end if

! ----------------------------------------------------------------------
! Find Voronoi charges
! ----------------------------------------------------------------------
    if (voropop) then
      call voronoi(no_s, na_s, na_u, nspin, indxua, &
          nmpl, rhoatm, DRho, iaorb, dvol, q)

#ifdef MPI
      call MPI_AllReduce(q(1,1),qtmp(1,1),nspin*na_u, &
          MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
      q(:,:) = qtmp(:,:)
#endif

      ! Note that the 'voronoi' routine is not really computing the
      ! partition of the scf density, but that of the 'deformation
      ! density' ("rho_scf - rho_atm", where rho_atm is the
      ! superposition of atomic charge densities)

      ! The routine to write the partial charges expects the total
      ! number of electrons as input, so we add the number of valence
      ! electrons in the 'charge' components of the generalized
      ! density array.
      
      select case ( nspin )
      case ( 1 )
        do ia = 1, na_u
          is = isa(ia)
          q(1,ia) = q(1,ia) + zvalfis(is)
        end do
      case ( 2, 4 )
        do ia = 1, na_u
          is = isa(ia)
          q(1:2,ia) = q(1:2,ia) + zvalfis(is) / 2
        end do
      end select

      call write_partial_charges(nspin, na_u, q, method=2)

    end if

    ! Add new-line
    if ( Node == 0 ) write(6, *) ''

#ifdef MPI
    call memory('D','D',size(qtmp),'hirsh')
    deallocate(qtmp)
#endif
    call memory('D','D',size(q),'hirsh')
    deallocate(q)

  contains

    subroutine write_partial_charges(nspin, na_u, q, method)
      use atmfuncs, only: labelfis
      use siesta_cml

      integer, intent(in) :: nspin, na_u
      ! On entry number of electrons for each spin-component.
      ! Note that this array sums up to the total number of electrons
      real(dp), intent(in) :: q(nspin,na_u)
      integer, intent(in) :: method ! 1 == Hirshfeld, 2 == Voronoi

      character(len=20) :: atm_label
      character(len=64) :: fmt
      integer :: ia, is
      real(dp) :: S, Stot, Svec(3)
      real(dp), allocatable :: qout(:,:)

      if ( Node /= 0 ) return

      ! Now write
      if ( method == 1 ) then
        write(6,"(/,a)") 'Hirshfeld Atomic Populations:'
        if ( cml_p ) then
          call cmlStartPropertyList(xf=mainXML, dictRef='siesta:hirshfeld', &
              title='Hirshfeld populations')
        end if
      else if ( method == 2 ) then
        write(6,"(/,a)") 'Voronoi Atomic Populations:'
        if ( cml_p ) then
          call cmlStartPropertyList(xf=mainXML, dictRef='siesta:voronoi', &
              title='Voronoi populations')
        end if
      end if

      ! Note that the first column is the net charge on the atom.
      ! Hence we do:
      !    zvalfis - q
      ! which amounts to a positive number for electron deficiency.
      ! The 2nd column is the total number of electrons on the atom.
      ! With both columns present it is clear how they relate to the
      ! number of electrons.

      select case ( nspin )
      case ( 1 )
        write(6,'(a6,1x,2(tr1,a9),2x,a)') 'Atom #', 'dQatom', 'Atom pop', 'Species'
        fmt = "(i6,1x,2(tr1,f9.5),2x,a)"

        allocate(qout(2,na_u))
        do ia = 1, na_u
          is = isa(ia)
          qout(1,ia) = zvalfis(is) - q(1,ia)
          qout(2,ia) = q(1,ia)
        end do

      case ( 2 )
        write(6,'(a6,1x,3(tr1,a9),2x,a)') 'Atom #', 'dQatom', 'Atom pop', 'Sz', 'Species'
        fmt = "(i6,1x,3(tr1,f9.5),2x,a)"

        allocate(qout(3,na_u))

        Stot = 0._dp
        do ia = 1, na_u
          is = isa(ia)
          ! this needs to be reversed (compared to mulliken)
          ! since the charges are with respect to atomic charge
          qout(1,ia) = zvalfis(is) - (q(1,ia) + q(2,ia))
          qout(2,ia) = q(1,ia) + q(2,ia)
          qout(3,ia) = q(1,ia) - q(2,ia)
          Stot = Stot + qout(3,ia)
        end do

        if ( cml_p ) then
          call cmlAddProperty(xf=mainXML, title="Atomic spin z", &
              dictRef='siesta:spin_z', value=qout(3,:), &
              units="siestaUnits:spin")
          call cmlAddProperty(xf=mainXML, title="Total atomic spin z", &
              dictRef='siesta:spin_z_tot', value=Stot, &
              units="siestaUnits:spin")
        end if

      case ( 4 )
        write(6,'(a6,1x,6(tr1,a9),2x,a)') 'Atom #', 'dQatom', 'Atom pop', &
            'S', 'Sx', 'Sy', 'Sz', 'Species'
        fmt = "(i6,1x,6(tr1,f9.5),2x,a)"

        ! Convert to spin-out
        allocate(qout(6,na_u))

        Svec(:) = 0._dp
        do ia = 1, na_u
          call spnvec(nspin, q(:,ia), qout(2,ia), qout(3,ia), qout(4:6,ia))
          is = isa(ia)
          qout(1,ia) = zvalfis(is) - qout(2,ia)
          Svec(:) = Svec(:) + qout(4:6,ia)
        end do
        Stot = sum(Svec(:) * Svec(:)) ** 0.5_dp

        if ( cml_p ) then
          call cmlAddProperty(xf=mainXML, title="Atomic spin", &
              dictRef='siesta:spin', value=qout(3,:), &
              units="siestaUnits:spin")
          call cmlAddProperty(xf=mainXML, title="Atomic spin x", &
              dictRef='siesta:spin_x', value=qout(4,:), &
              units="siestaUnits:spin")
          call cmlAddProperty(xf=mainXML, title="Atomic spin y", &
              dictRef='siesta:spin_y', value=qout(5,:), &
              units="siestaUnits:spin")
          call cmlAddProperty(xf=mainXML, title="Atomic spin z", &
              dictRef='siesta:spin_z', value=qout(6,:), &
              units="siestaUnits:spin")
          call cmlAddProperty(xf=mainXML, title="Total atomic spin", &
              dictRef='siesta:spin_tot', value=Stot, &
              units="siestaUnits:spin")
          call cmlAddProperty(xf=mainXML, title="Total atomic spin_xyz", &
              dictRef='siesta:spin_xyz_tot', value=Svec, &
              units="siestaUnits:spin")
        end if

      case default
        call die("partial_charges: could not figure out the spin configuration")
      end select

      if ( cml_p ) then
        call cmlAddProperty(xf=mainXML, title="Atomic net charge", &
            dictRef='siesta:atom_net_charge', value=qout(1,:), &
            units="siestaUnits:|e|")
        call cmlAddProperty(xf=mainXML, title="Atomic population", &
            dictRef='siesta:atom_charge', value=qout(2,:), &
            units="siestaUnits:e")
      end if

      ! Write out to std-out
      do ia = 1,na_u
        is = isa(ia)
        atm_label = labelfis(is)
        write(6,fmt) ia, qout(:,ia), trim(atm_label)
      end do

      deallocate(qout)

      if ( cml_p ) then
        ! Close property
        call cmlEndPropertyList(xf=mainXML)
      end if

      ! finaly print out total spin
      select case ( nspin )
      case ( 2 )
        write(6,'(a)') repeat('-',6+1+3*10)
        write(6,'(a6,22x,f9.5)') 'Total', Stot
      case ( 4 )
        write(6,'(a)') repeat('-',6+1+6*10)
        write(6,'(a6,21x,4(tr1,f9.5))') 'Total', Stot, Svec
      end select

    end subroutine write_partial_charges

  end subroutine compute_partial_charges


  subroutine hirshfeld( no_s, na_s, na_u, nspin, indxuo, indxua, np, &
      Datm, rhoatm, rhoscf, iaorb, iphorb, isa, dvol, q)
! ********************************************************************
! Finds the Hirshfeld atomic charges  
! Hirshfeld, Theo Chem Acta 44, 129 (1977)
! See Fonseca et al, J. Comp. Chem. 25, 189 (2003)
! 
! Written by P.Ordejon. August'04.
! *********************** InpUT **************************************
! integer no_s              : Number of basis orbitals
! integer na_s              : Number of atoms
! integer na_u              : Number of atoms in unit cell
! integer nspin             : Number of spins
! integer indxuo(no_s)      : Index of equivalent orbital in unit cell
! integer indxua(na_s)      : Index of equivalent atom in unit cell
! integer np                : Number of mesh points
! real*8 Datm(no_s)         : Occupations of basis orbitals in free atom
! integer iaorb(no_s)       : Pointer to atom to which orbital belongs
! integer iphorb(no_s)      : Orbital index within each atom
! integer isa(na_s)         : Species index of all atoms
! real rhoatm(nsp,np)       : Harris (sum of atoms) density at mesh points
! real rhoscf(nsp,np,nspin) : Selfconsistent charge density at mesh points
! *********************** OUTPUT **************************************
! real*8 q(nspin,na_u)      : Hirshfeld electron charge on each atom
! *********************************************************************
!
!  Modules

    use precision, only: dp, grid_p
    use atmfuncs, only: rcut, phiatm
    use mesh, only: nsp, dxa, xdop, xdsp
    use meshphi

    integer :: no_s, np, nspin, na_s, na_u
    integer :: indxuo(no_s), indxua(na_s), iaorb(no_s), iphorb(no_s), isa(na_s)
    real(grid_p) :: rhoatm(nsp,np),rhoscf(nsp,np,nspin)
    real(dp) :: Datm(no_s), phip, dvol, q(nspin,na_u)

    integer :: i, ip, isp, iu, kn, io, iop, is, iphi, ia, ix
    integer :: ispin, iua
    real(dp) :: gradCi(3), r2o, r2sp, dxsp(3), Qi

    ! A very small number to avoid division by zero
    real(grid_p), parameter :: rhoatm_tolerance = 1.0e-12_grid_p

! Initialise Hirshfeld charges
    q(:,:) = 0._dp

!  Loop on mesh points
    do ip = 1,np

!  Loop on orbitals of mesh point
      do kn = 1+endpht(ip-1), endpht(ip)
        i = lstpht(kn)
        iu = indxuo(i)

!  Generate phi value and loop on subpoints
        iphi = iphorb(i)
!!          if (iphi .gt. no_s) call die("hirshfeld: iphi error")
        ia = iaorb(i)
        iua = indxua(ia)
        is = isa(ia)
        r2o = rcut(is,iphi)**2
        iop = listp2(kn)
        do isp = 1,nsp
          if (abs(rhoatm(isp,ip)) <= rhoatm_tolerance) then
                ! Atomic charge is really too small here... so
                ! we do not count this point
            CYCLE
!$$$                 write(6,"(a,g14.6,/,10x,a,2g14.6)")
!$$$     $               "rhoatm: ", rhoatm(isp,ip), "rho: ",
!$$$     $                rhoscf(isp,ip,:)
          end if

          do ix = 1,3
            dxsp(ix) = xdop(ix,iop) + xdsp(ix,isp) - dxa(ix,ia)
          end do
          r2sp = dxsp(1)**2 + dxsp(2)**2 + dxsp(3)**2
          if ( r2sp < r2o ) then
            call phiatm(is,iphi,dxsp,phip,gradCi)

            ! This translates to the fraction of atomic charge
            ! associated with orbital `iu` on atom `iua`
            ! i.e. rhoatm(isp,ip) = sum_atoms Datm(atom) * phip**2
            Qi = Datm(iu) * phip * phip / rhoatm(isp,ip)
            do ispin = 1, nspin
              q(ispin,iua) = q(ispin,iua) + Qi * rhoscf(isp,ip,ispin)
            end do

          end if
        end do

      end do

    end do

    ! Correct charges so they account for voxel volume
    ! This sums up to be:
    !  q(ia) = \int_rr rho(rr) * rho_atom(rr, ia) / sum_atoms rho_atomic(rr, atom)
    do ia = 1 , na_u
      q(:,ia) = q(:,ia) * dvol
    end do

  end subroutine hirshfeld


  subroutine voronoi( no_s, na_s, na_u, nspin, indxua, np, &
      rhoatm, rhoscf, iaorb, dvol, q)
! ********************************************************************
! Finds the Voronoi deformation density
! Bickelhaupt et al, Organometallics 15, 2923 (1996)
! See Fonseca et al, J. Comp. Chem. 25, 189 (2003)
! 
! Written by P.Ordejon. August'04.
! *********************** InpUT **************************************
! integer no_s              : Number of basis orbitals
! integer na_s              : Number of atoms
! integer na_u              : Number of atoms in unit cell
! integer nspin             : Number of spins
! integer indxua(na_s)      : Index of equivalent atom in unit cell
! integer np                : Number of mesh points
! integer iaorb(no_s)       : Pointer to atom to which orbital belongs
! real rhoatm(nsp,np)       : Harris (sum of atoms) density at mesh points
! real rhoscf(nsp,np,nspin) : Selfconsistent charge density at mesh points
! *********************** OUTPUT **************************************
! real*8 q(nspin,na_u)      : Voronoi deformation density in units of electrons on each atom
! *********************************************************************

    use precision, only: dp, grid_p
    use atmfuncs, only: rcut, phiatm, zvalfis
    use mesh,     only: nsp, dxa, xdop, xdsp
    use meshphi

    integer :: no_s, np, nspin, na_s, na_u
    integer :: indxua(na_s), iaorb(no_s)

    real(grid_p) :: rhoatm(nsp,np),rhoscf(nsp,np,nspin)

    real(dp) :: dvol, q(nspin,na_u)
    
    integer :: i, ip, isp, iu, kn, io, iop, is, ia, ix
    integer :: ispin, iua, neq, ieq

    real(dp) :: phip, gradCi(3), r2o, r2sp, dxsp(3), Qi, rmin, qtot

    integer, parameter :: maxeq = 20
    real(dp), parameter :: huge = 1.e30_dp
    real(dp), parameter :: tol = 1.e-4_dp

    integer :: nsd
    logical :: eq
    integer :: iatom(maxeq)

!  Initialise Voronoi charges
    q(:,:) = 0._dp

    ! Number of diagonal spin-components
    nsd = min(2, nspin)

!  Loop on mesh points
    do ip = 1,np

!  Loop on mesh subpoints
      do isp = 1,nsp
        rmin = huge
        eq = .false.
        neq = 1
        ! We only need first element to check
        iatom(1) = 0

!  Loop on orbitals of mesh point, to check which atom is closest to subpoint
        do kn = 1+endpht(ip-1), endpht(ip)
          i = lstpht(kn)
          ia = iaorb(i)

!  iua is the index of the atom in unit cell
          iua = indxua(ia)

          iop = listp2(kn)
          do ix = 1,3
            dxsp(ix) = xdop(ix,iop) + xdsp(ix,isp) - dxa(ix,ia)
          end do
          r2sp = dxsp(1)**2 + dxsp(2)**2 + dxsp(3)**2

! If distance is equal to the previous minimum, determine if it is
! another atom, or the same. If it is a different one, then
! handle multiplicity of nearest atoms to grip point
! Consider distances equal if within tolerance

          if (dabs(rmin-r2sp) .lt. tol) then
            eq = .false.
            do ieq = 1,neq
              if (iua .eq. iatom(ieq)) eq = .true.
            end do
            if (.not. eq) then
              neq = neq+1
              if (neq .gt. maxeq) stop 'voronoi: increase maxeq'
              iatom(neq) = iua
            end if
            goto 100
          end if

          if (r2sp .lt. rmin) then
            neq = 1
            iatom(neq) = iua
            rmin = r2sp
          end if

100       continue

        end do

!  Assign charge to atom iatom; if no atom was found, then the charge
!  is zero, so move to next subpoint

        if (iatom(1) .ne. 0) then
          do ispin = 1, nsd
            Qi = (rhoscf(isp,ip,ispin) - rhoatm(isp,ip)/nsd)/neq
            do ieq = 1,neq
              q(ispin,iatom(ieq)) = q(ispin,iatom(ieq)) + Qi
            end do
          end do
          do ispin = 3, nspin
            Qi = rhoscf(isp,ip,ispin) / neq
            do ieq = 1,neq
              q(ispin,iatom(ieq)) = q(ispin,iatom(ieq)) + Qi
            end do
          end do
        else
!  Check that charge is actually zero (both should be 0)
          qtot = 0.0_dp
          do ispin = 1 , nsd
            qtot = qtot + rhoscf(isp,ip,ispin) + rhoatm(isp,ip)
          end do
          if ( qtot > 0._dp ) stop 'voronoi: Error in grid charge'
        end if

      end do

    end do

    ! Correct charges so they account for voxel volume
    ! This will sum up to 0
    do ia = 1 , na_u
      q(:,ia) = q(:,ia) * dvol
    end do

  end subroutine voronoi

end module m_partial_charges
