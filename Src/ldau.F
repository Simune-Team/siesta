! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module m_ldau
      use precision,     only : dp


!#ifdef MPI
      use m_mpi_utils, only: globalize_sum, globalize_max
!#endif

      implicit none

      private

      public :: hubbard_term
      logical, save, public :: conv_ldaupop  ! Are the populations of the LDA+U 
                                             !   projectors converged?
    
      integer, save :: iter_ldaupop = 0

      CONTAINS

      subroutine hubbard_term( scell, nua, na, isa, xa, indxua,
     .                         maxnh, maxnd, lasto, iphorb, 
     .                         numd, listdptr, listd, numh, 
     .                         listhptr, listh, nspin, Dscf, 
     .                         Eldau, DEldau,
     .                         fa, stress, H, iscf,
     .                         matrix_elements_only )      
C *********************************************************************
C Calculates Hubbard-like U term contribution to total 
C energy, atomic forces, stress and hamiltonian matrix elements.
C Energies in Ry. Lengths in Bohr.
C Merged into the trunk by J. Junquera, March 2016.
C Written by D. Sanchez-Portal, October 2008, 
C after subroutine nlefsm by J.Soler and P.Ordejon (June 1997).
C Based on a previous version by S. Riikonen and D. Sanchez-Portal (2005)
C **************************** INPUT **********************************
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C integer isa(na)          : Species index of each atom
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
C integer indxua(na)       : Index of equivalent atom in unit cell
C integer maxnh            : First dimension of H and listh
C integer maxnd            : Maximum number of elements of the
C                            density matrix
C integer lasto(0:na)      : Position of last orbital of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom,
C                            where no=lasto(na)
C integer numd(nuo)        : Number of nonzero elements of each row of the
C                            density matrix
C integer listdptr(nuo)    : Pointer to the start of each row (-1) of the
C                            density matrix
C integer listd(maxnd)     : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer numh(nuo)        : Number of nonzero elements of each row of the
C                            hamiltonian matrix
C integer listhptr(nuo)    : Pointer to the start of each row (-1) of the
C                            hamiltonian matrix
C integer listh(maxnh)     : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer nspin            : Number of spin components
C integer iscf             : Counter of the cycles of SCF iterations
C real*8  Dscf(maxnd,nspin): Density matrix
C logical matrix_elements_only: Determine whether the only the matrix elements
C                            of the Hamiltonian are computed, or also the
C                            forces and stresses
C ******************* INPUT and OUTPUT *********************************
C real*8 fa(3,na)          : NL forces (added to input fa)
C real*8 stress(3,3)       : NL stress (added to input stress)
C real*8 H(maxnh,nspin)    : NL Hamiltonian (added to input H)
C **************************** OUTPUT *********************************
C real*8 Eldau             : U-Hubbard energy 1 
C real*8 DEldau            : U-hubbard energy 2 
C *********************************************************************
C
C  Modules
C
      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : GetNodeOrbs, LocalToGlobalOrb
      use parallelsubs,  only : GlobalToLocalOrb
      use atmfuncs,      only : rcut, orb_gindex, ldau_gindex
      use neighbour    , only : iana=>jan, r2ki=>r2ij, xki=>xij
      use neighbour    , only : mneighb
      use alloc,         only : re_alloc, de_alloc
      use m_steps, only: istp
      use m_new_matel,   only : new_matel
      use atm_types,     only : species_info    ! Derived type with all the info
                                                !   about the radial functions
                                                !   (PAOs, KB projectors,
                                                !   LDA+U proj,
                                                !   VNA potentials, etc)
                                                !   for a given atomic specie
      use atm_types,   only : species           ! Actual array where the
                                                !   previous information is
                                                !   stored
      use basis_types, only : basis_def_t       ! Derived type where all the
                                                !   information relative to the
                                                !   definition of the basis set
                                                !   is defined
      use basis_types, only : ldaushell_t       ! Derived type where all the
                                                !   information relative to the
                                                !   atomic orbitals where the U
                                                !   correction will be applied
                                                !   is defined
      use basis_types, only : basis_parameters  ! Derived type where all the
                                                !   information about the
                                                !   - basis set
                                                !   - Kleinman-Bylander proj.
                                                !   - LDA+U proj.
                                                !   ...
                                                !   for all the species
                                                !   are defined
      use ldau_specs,  only : dDmax_threshold   ! Parameter that defines the
                                                !   criterium required to start
                                                !   or update the calculation of
                                                !   the populations of
                                                !   the LDA+U projections
      use ldau_specs,  only : dtol_ldaupop      ! Parameter that defines the
                                                !   convergence criterium 
                                                !   for the LDA+U local 
                                                !   population
      use ldau_specs,  only : ldau_init         ! Flag that determines whether 
                                                !   the local populations are 
                                                !   calculated on the 
                                                !   first iteration
      use ldau_specs,  only : ldau_shift        ! Flag that determines whether 
                                                !   the U parameter
                                                !   is interpreted
                                                !   as a local potential shift
      use m_compute_max_diff, only: dDmax_current

      integer, intent(in) ::
     .   maxnh, na, maxnd, nspin, nua, iscf

      integer, intent(in)  ::
     .  indxua(na), iphorb(*), isa(na),  
     .  lasto(0:na), listd(maxnd), listh(maxnh),
     .  numd(*), numh(*), listdptr(*), listhptr(*)

      real(dp), intent(in) :: scell(3,3), Dscf(maxnd,nspin),
     .                        xa(3,na)
      real(dp), intent(inout) :: fa(3,nua), stress(3,3)
      real(dp), intent(inout) :: H(maxnh,nspin)
      real(dp), intent(out)   :: Eldau, DEldau
      logical, intent(in)     :: matrix_elements_only

      real(dp) ::   volcel
      external ::   timer, volcel

C Internal variables ................................................
C maxno  = maximum number of basis orbitals overlapping a KB projector

      integer, save ::  maxno = 2000
      integer, save ::  maxldau       ! Maximum number of LDA+U projectors 
                                      ! on a given atom 
                                      ! (used to allocate some matrices)

      integer
     .  ia, ishldau, ina, ind, ino,
     .  io, iio, ioa, is, ispin, ix, ikb, 
     .  j, jno, jo, jx, ka, ko, koa, ks, kua,
     .  nldauproj, nna, nno, no, nuo, nuotot, ja,
     .  lko, lkb, nprin_ko, nprin_kb, kg, ig 

      integer, save :: istp_last      ! Last MD step when the Hldau 
                                      !   was reallocated
      integer, dimension(:), pointer :: iano, iono

      real(dp)
     .  Cijk, epsk, fik, rki, rmax, rmaxldau, rmaxo, 
     .  Sik, Sjk, volume,  oc(2), Ueff, dn, dnmax,
     .  Dij, stress_ldau(3,3)

      real(dp), dimension(:,:), pointer :: Vi, Di
      real(dp), dimension(:,:), pointer :: Ski, xno
      real(dp), dimension(:,:), pointer, save :: Hldau(:,:)    ! Hubbard term
                                                               ! to be included
                                                               ! in the 
                                                               ! Hamiltonian
      real(dp), dimension(:,:,:), pointer :: grSki
      real(dp), dimension(:,:,:,:), pointer, save :: occu      ! Array used to 
                                                               ! store the 
                                                               ! occupations of
                                                               ! the LDA+U proj.
      real(dp), dimension(:,:,:,:), pointer, save :: occu_old  ! Same as occu
                                                               ! but in the 
                                                               ! previous step

#ifdef MPI
      real(dp), dimension(:,:,:,:), pointer, save :: buffer1   ! Local copies 
                                                               ! of occu 
                                                               ! (only for 
                                                               ! parallel
                                                               ! execution)
#endif


      logical :: first                            ! For a given set of 
                                                  !   atomic positions
                                                  !   determine whether this 
                                                  !   is the first step in the
                                                  !   SCF iterations
      logical ::   within
      logical, dimension(:), pointer ::  listed, listedall
      logical, save :: firstime = .true.          !  First time that this 
                                                  !   subroutine is called?

      type(species_info),  pointer :: spp
      type(basis_def_t),   pointer :: basp
      type(ldaushell_t),   pointer :: ldaushell
C ......................

!     Start time counter
      call timer( 'hubbard_term', 1 )

      ! Nullify pointers
      nullify( grSki, Ski, xno, iono, iano )
      nullify( listedall, listed, Vi, Di )
      
      Eldau = 0.0_dp

!     Determine whether this is the first SCF step for a 
!     given atomic configuration
      first = (iscf == 1)

!!     For debugging
!      write(6,'(a,i5)')'hubbard_term: iscf  = ', iscf
!      write(6,'(a,l5)')'hubbard_term: first = ', first
!!     End debugging

!     Initialization and allocation of matrices
      if( firstime ) then 
!       Find maximum number of LDA+U projectors on a given atom 
        maxldau = 0
        do ka = 1, na
          is  = isa(ka)
          spp => species(is)
          maxldau = max(maxldau,spp%nprojsldau)
        enddo
!!       For debugging
!        write(6,'(a,i5)')'hubbard_term: maxldau = ', maxldau
!        write(6,'(a,i5)')'hubbard_term: nua     = ', nua
!        write(6,'(a,i5)')'hubbard_term: nspin   = ', nspin
!!       End debugging

!       Allocate local array to store the occupations of the 
!       LDA+U projectors
        nullify( occu, occu_old )
        allocate( occu(maxldau,maxldau,nua,nspin) )
        call memory( 'A', 'D', size(occu), 'hubbard_term' )
        allocate( occu_old(maxldau,maxldau,nua,nspin) )
        call memory( 'A', 'D', size(occu_old), 'hubbard_term' )
        occu_old = 0.0_dp
#ifdef  MPI
        nullify( buffer1 )
	allocate(buffer1(maxldau,maxldau,nua,nspin))
        call memory('A','D',size(buffer1),'hubbard_term')
#endif

        nullify( Hldau )
        call re_alloc( Hldau, 1, maxnh, 1, nspin, name='Hldau',
     &                 routine = 'hubbard_term', shrink=.false., 
     &                 copy=.false. )
        Hldau=0.0_dp 

        if( ldau_shift ) ldau_init = .true.
        istp_last = istp
        firstime  = .false.
      endif
!!     For debugging
!      write(6,'(a,l1)')'hubbard_term: ldau_shift = ', ldau_shift
!      write(6,'(a,l1)')'hubbard_term: ldau_init  = ', ldau_init
!      write(6,'(a,i5)')'hubbard_term: istp_last  = ', istp_last
!      write(6,'(a,l1)')'hubbard_term: firstime   = ', firstime
!!     End debugging
!     End initialization

!     If the system has moved, reallocate Hldau and set ldau_init to true
      if( istp_last .ne. istp ) then 
        call re_alloc( Hldau, 1, maxnh, 1, nspin, name='Hldau',
     &                 routine='hubbard_term', shrink=.false., 
     &                 copy=.false. )
        ldau_init = .true.
        istp_last = istp
      endif

!     Here we determine whether the occupations are computed on the first
!     SCF step or not.
      if( first .and. (.not. ldau_init) ) then
         if ( Node == 0 ) then
            write(6,'(a)')
     &           'hubbard_term: not computing the occupations '
            write(6,'(a)')
     &           '              in the first SCF step'
         end if
         call timer( 'hubbard_term', 2 )
         return
      endif
      
!!     For debugging
!      write(6,'(a,l1)')'hubbard_term: first = ', first
!      write(6,'(a,f12.5)')
!     &  'hubbard_term: dDmax_current   = ', dDmax_current
!      write(6,'(a,f12.5)')
!     &  'hubbard_term: dDmax_threshold = ', dDmax_threshold
!      write(6,'(a,l1)')'hubbard_term: compar = ', 
!     &  (dDmax_current .lt. dDmax_threshold)
!      write(6,'(a,l1)')'hubbard_term: ldau_shift = ', 
!     &  ldau_shift
!!     End debugging
      occupations: if( first .or.
     &                 (dDmax_current .lt. dDmax_threshold) .or.
     &                 ldau_shift ) then 


!       Find unit cell volume
        volume = volcel( scell ) * nua / na

!       Find maximum range of the atomic orbitals (rmaxo) 
!       and of the LDA+U projectors (rmaxldau)
        rmaxo    = 0.0_dp
        rmaxldau = 0.0_dp
        do ia = 1, nua
          is = isa(ia)
          basp => basis_parameters(is)
!!         For debugging
!          write(6,'(a,2i5)')' is, basp%nldaushells_tmp = ', 
!     &                        is, basp%nldaushells_tmp
!!         End debugging
!         Loop on all the LDA+U projectors for a given species
          do ishldau = 1, basp%nldaushells_tmp
            ldaushell => basp%tmp_ldaushell(ishldau)
            rmaxldau = max( rmaxldau, ldaushell%rc )
!!           For debugging
!            write(6,'(a,2i5,f12.5)')' is, ishldau, ldaushell%rc = ', 
!     &                                is, ishldau, ldaushell%rc
!!           End debugging
          enddo
          do io = lasto(ia-1)+1,lasto(ia)
            ioa = iphorb(io)
            rmaxo = max( rmaxo, rcut(is,ioa) )
          enddo
        enddo
        rmax = rmaxo + rmaxldau
!!       For debugging
!        write(6,'(a,3f12.5)')'hubbard_term: rmax, rmaxo, rmaxldau = ',
!     &                               rmax, rmaxo, rmaxldau
!!       End debugging

!       Initialize arrays Di and Vi only once
        no = lasto(na)
        nuotot = lasto(nua)
        call GetNodeOrbs(nuotot,Node,Nodes,nuo)

!       Allocate local memory
        call re_alloc( Di, 1, no, 1, nspin, 
     &                 name='Di', routine='hubbard_term' )
        Di = 0.0_dp

        call re_alloc( Vi, 1, no, 1, nspin, 
     &                 name='Vi', routine='hubbard_term' )
        Vi = 0.0_dp

        call re_alloc( listed, 1, no, 
     &                 name='listed', routine='hubbard_term')
        listed(1:no) = .false.

        call re_alloc( listedall, 1, no, 
     &                 name='listedall', routine='hubbard_term' )
        listedall(1:no) = .false.

!       Make list of all orbitals needed for this node
        do io = 1, nuo
          call LocalToGlobalOrb(io,Node,Nodes,iio)
          listedall(iio) = .true.
          do j = 1, numh(io)
            jo = listh(listhptr(io)+j)
            listedall(jo) = .true.
          enddo
        enddo

!!       For debugging
!        do io = 1, nuo
!          write(6,'(a,2i5)')'io, numh(io) = ', io, numh(io)
!          do j = 1,numh(io)
!            jo = listh(listhptr(io)+j)
!            write(6,'(a,2i5)')'  j, neigh = ', j, jo
!          enddo
!        enddo
!
!        do io = 1, no
!          write(6,'(a,i5,l3)')'io, listedall(io) = ', io, listedall(io)
!        enddo
!!       End debugging

!       Allocate local arrays that depend on saved parameters
        call re_alloc( iano, 1, maxno, 
     &                 name='iano', routine='hubbard_term' )
        call re_alloc( iono, 1, maxno, 
     &                 name='iono', routine='hubbard_term' )
        call re_alloc( xno, 1, 3, 1, maxno, 
     &                 name='xno', routine='hubbard_term' )
        call re_alloc( Ski, 1, maxldau, 1, maxno, 
     &                 name='Ski', routine='hubbard_term' )
        call re_alloc( grSki, 1, 3, 1, maxldau, 1, maxno, 
     &                 name='grSki', routine='hubbard_term' )
      
!       Counter for the SCF loops to converge the population of the 
!       LDA+U projectors
        iter_ldaupop = iter_ldaupop + 1
        if( Node == 0 ) write(6,'(a,i4)')
     &   'hubbard_term: recalculating local occupations ', iter_ldaupop

!       Initialize occupations
        occu = 0.0_dp
#ifdef  MPI
        buffer1 = 0.0_dp
#endif

!       Initialize neighb subroutine
        call mneighb( scell, rmax, na, xa, 0, 0, nna )

!       Loop on atoms with LDA+U projectors      
        do ka = 1, na
          kua = indxua(ka)
          ks = isa(ka)
          spp => species(ks)
          basp => basis_parameters(ks)
          nldauproj = spp%nprojsldau
          if( nldauproj == 0 ) cycle

!         Find neighbour atoms 
          call mneighb( scell, rmax, na, xa, ka, 0, nna )
!!         For debugging
!          do io = 1, 3
!            write(6,'(a,3f12.5)')' scell = ', scell(io,:)
!          enddo
!          write(6,'(a,f12.5)')' rmax = ', rmax
!          write(6,'(a,i5)')' na   = ', na
!          do io = 1, na
!            write(6,'(a,3f12.5)')' xa = ', xa(:,io)
!          enddo
!          write(6,'(a,i5)')' nna   = ', nna
!!         End debugging

!         Find neighbour orbitals
          nno = 0
          do ina = 1, nna
            ia = iana(ina)
            is = isa(ia)
            rki = sqrt(r2ki(ina))
            do io = lasto(ia-1)+1, lasto(ia)

!             Only calculate if needed locally
              if (listedall(io)) then
                ioa = iphorb(io)

!               Find if orbital is within range
                within = .false.
                do ko = 1, nldauproj
                  ldaushell => basp%tmp_ldaushell(spp%pjldau_index(ko))
!!                 For debugging
!                  write(6,'(a,3i5,f12.5,2i5,f12.5)')
!     .                      'ka, ks, ko, ldaushell%rc = ', 
!     .                       ka, ks, ko, ldaushell%rc, ina, io, 
!     .                       rcut(is,ioa)
!!                 End debugging
                  if (rki .lt. rcut(is,ioa)+ldaushell%rc) 
     .               within = .true.
                enddo

!               Find overlap between neighbour orbitals and LDA+U projectors
                if (within) then
!! For debugging
!                  write(6,'(a,3i5,f12.5)')'hubbard_term: ka, ina, io, rki = ',
!     .              ka, ina, io, rki
!! End debugging
!                 Check maxno - if too small then increase array sizes
                  if (nno.eq.maxno) then
                    maxno = maxno + 100
                    call re_alloc( iano, 1, maxno, name='iano',
     &                       copy=.true., routine='hubbard_term' )
                    call re_alloc( iono, 1, maxno, name='iono',
     &                       copy=.true., routine='hubbard_term' )
                    call re_alloc( xno, 1, 3, 1, maxno, name='xno', 
     &                       copy=.true., routine='hubbard_term' )
                    call re_alloc( Ski,1, maxldau, 1, maxno, name='Ski',
     &                       copy=.true., routine='hubbard_term' )
                    call re_alloc( grSki, 1, 3, 1, maxldau, 1, maxno,
     &                       name='grSki', routine='hubbard_term', 
     &                       copy=.true. )
                  endif
                  nno = nno + 1
!                 The nno-eme neighbour orbital of the LDA+U projector is io,
!                 where io runs between 1 and the total number of orbitals
!                 in the supercell
                  iono(nno) = io

!                 The nno-eme neighbour orbital of the LDA+U projector belongs
!                 to atom ia,
!                 where ia runs between 1 and the total number of atoms
!                 in the supercell
                  iano(nno) = ia

!                 The relative position between the center of the LDA+U proj.
!                 and the center of the nno-eme neighbour orbital 
!                 is xno 
                  do ix = 1,3
                    xno(ix,nno) = xki(ix,ina)
                  enddo

!                 Here we compute the overlap between a LDA+U projector
!                 and a neighbour atomic orbital
                  do ko = 1, nldauproj
                    ioa = iphorb(io)
                    kg  = ldau_gindex(ks,ko)
                    ig  = orb_gindex(is,ioa)
                    call new_matel( 'S', kg, ig, xki(1:3,ina),
     &                  Ski(ko,nno), grSki(1:3,ko,nno) )
!!                   For debugging
!                    write(6,'(4i5,9f10.5)')
!     &                ka, ko, io, ioa, xki(1:3,ina), 
!     &                Ski(ko,nno), grSki(1:3,ko,nno)  
!!                   End debugging
                  enddo

                endif   ! If on orbitals within range

              endif     ! If on orbitals within the local node

            enddo       ! End loop on neighbour orbitals

          enddo         ! End loop on neighbour atoms

!         Loop on neighbour orbitals
          do ino = 1,nno
            io = iono(ino)
            ia = iano(ino)

            call GlobalToLocalOrb(io,Node,Nodes,iio)
            if (iio.gt.0) then
!             Valid orbital
              if(ia .le. nua) then
!               Scatter density matrix row of orbital io
                do j = 1, numd(iio)
                  ind = listdptr(iio) + j
                  jo = listd(ind)
                  do ispin = 1, nspin
                    Di(jo,ispin) = Di(jo,ispin) + 
     &                             Dscf(ind,ispin)
                  enddo
                enddo

!               Scatter filter of desired matrix elements
                do j = 1,numh(iio)
                  jo = listh(listhptr(iio)+j)
                  listed(jo) = .true.
                enddo

!               Find matrix elements with other neighbour orbitals
                do jno = 1,nno
                  jo = iono(jno)
                  ja = iano(jno)
                  if (listed(jo)) then

!                   Loop on LDA+U projectors
                    do ko = 1, nldauproj
                      do ikb = 1, nldauproj
                        Sik = Ski(ko,ino)
                        Sjk = Ski(ikb,jno)
                        do ispin = 1, nspin
                          Dij = Di(jo,ispin) 
                          occu(ko,ikb,kua,ispin) = 
     &                    occu(ko,ikb,kua,ispin) + 
     &                    Dij * Sik * Sjk/(3.0_dp-dble(nspin))
!!                         For debugging
!                          write(6,'(a,4i5,4f12.5)')
!     &                     'ko, ikb, kua, ispin, occu = ',
!     &                      ko, ikb, kua, ispin, occu(ko,ikb,kua,ispin),
!     &                      Dij, Sik, Sjk
!!                         End debugging
                        enddo 
                      enddo
                    enddo
                  endif
                enddo

!               Restore  Di and listed
                do j = 1, numh(iio)
                  ind = listhptr(iio) + j
                  jo = listh(ind)
                  listed(jo) = .false.
                enddo
                do j = 1,numd(iio)
                  jo = listd(listdptr(iio)+j)
                  do ispin=1,nspin
                     Di(jo,ispin) = 0.0_dp
                  enddo
                enddo

              endif        ! Endif the neighbout atom is within the unit cell

            endif          ! Endif (iio .gt. 0)

          enddo            ! End loop on neighbour orbitals (ino)

        enddo              ! End of the loop on atoms with LDA+U projectors
!! For debugging        
!        do kua = 1, nua
!          ks = isa(kua)
!          spp => species(ks)
!          nldauproj = spp%nprojsldau
!          do ko = 1, nldauproj
!            do ikb = 1, nldauproj
!              do ispin = 1, nspin
!                write(6,'(a,4i5,f12.5)')
!     &            'ko, ikb, kua, ispin, occu = ',
!     &             ko, ikb, kua, ispin, occu(ko,ikb,kua,ispin)
!              enddo
!            enddo
!          enddo
!        enddo
!        stop
!! End debugging        
     
        dnmax=0.0_dp 
        oc=0.0_dp 
#ifdef MPI
!       Global reduction of occupation
        do ka=1,nua
          is=isa(ka)
          spp  => species(is)
          nldauproj = spp%nprojsldau
          if( nldauproj .ne. 0) then
            do ispin=1,nspin             
              call globalize_sum(occu(1:nldauproj,1:nldauproj,ka,ispin),
     &             buffer1(1:nldauproj,1:nldauproj,ka,ispin))
            enddo 
          endif
        enddo 
        occu=buffer1
#endif
        do ka=1,nua
          is=isa(ka)
          spp => species(is)
          nldauproj = spp%nprojsldau
          if( nldauproj .ne. 0) then
            oc=0.0_dp

            if( ldau_shift .and. Node == 0 ) then
              write(6,*) 'hubbard_term: projector occupations'
              write(6,*) 'hubbard_term: atom, species: ',ka, is
            endif
            do ko = 1, nldauproj
              do ikb = 1, nldauproj
!                 write(6,'(2i4,2f12.5)') ko, ikb,
!    &          (occu_old(ko,ikb,ka,ispin),ispin=1,nspin)
!
                if( ldau_shift .and. Node == 0 ) then
                   write(6,'(2i4,2f12.5)')  ko, ikb,
     &             (occu(ko,ikb,ka,ispin),ispin=1,nspin)
                endif

                do ispin = 1, nspin
                  dn=occu(ko,ikb,ka,ispin)-occu_old(ko,ikb,ka,ispin)
                  dnmax = max(dnmax,dabs(dn))
                enddo 
              enddo 
              if(ldau_shift) then
                oc(1) = oc(1) + occu(ko,ko,ka,1)
                oc(2) = oc(2) + occu(ko,ko,ka,2)
              endif
            enddo 
            if(ldau_shift.and.Node == 0) then
              write(6,'(a,/,a,3f12.6)') 
     &           'hubbard_term: Total projector shell',
     &           'Occupations: ', oc(1), oc(2),oc(1)+oc(2)
            endif
          endif
        enddo 
!!       For debugging
!        write(6,'(a,2f12.6)')'hubbard_term: occup ', oc(1), oc(2)
!!       End debugging
#ifdef MPI
        call globalize_max(dnmax,buffer1(1,1,1,1))
        dnmax = buffer1(1,1,1,1)
#endif
        if ( Node == 0 ) then
           write(6,'(a,f12.6)')
     &    'hubbard_term: maximum change in local occup.',dnmax
        endif

        conv_ldaupop = .false.
        if( dnmax .lt. dtol_ldaupop ) conv_ldaupop = .true. 
        if( .not. conv_ldaupop ) occu_old = occu

        recalc_hamilt: if( .not. conv_ldaupop .or. first .or. 
     &                     .not. matrix_elements_only ) then
          if(Node.eq.0) then 
            if(matrix_elements_only) then
              write(6,'(a)')'hubbard_term: recalculating Hamiltonian'
            else
              write(6,'(a)')
     &          'hubbard_term: recalculating Hamiltonian and forces'      
            endif
          endif

          Hldau = 0.0_dp

!         Initialize neighb subroutine
          call mneighb( scell, rmax, na, xa, 0, 0, nna )

          do ka = 1, na
            kua = indxua(ka)
            ks = isa(ka)
            spp => species(ks)
            nldauproj = spp%nprojsldau
            basp => basis_parameters(ks)
            if( nldauproj == 0 ) cycle

!           Find neighbour atoms
            call mneighb( scell, rmax, na, xa, ka, 0, nna )

!           Find neighbour orbitals
            nno = 0
            do ina = 1, nna
              ia = iana(ina)
              is = isa(ia)
              rki = sqrt(r2ki(ina))
              do io = lasto(ia-1)+1, lasto(ia)

!               Only calculate if needed locally
                if (listedall(io)) then
                  ioa = iphorb(io)

!                 Find if orbital is within range
                  within = .false.
                  do ko = 1, nldauproj
                    ldaushell=>basp%tmp_ldaushell(spp%pjldau_index(ko))      
                    if (rki .lt. rcut(is,ioa)+ldaushell%rc) 
     .                within = .true.
                  enddo

!                 Find overlap between neighbour orbitals and LDA+U projectors
                  if (within) then
!                   Check maxno - if too small then increase array sizes
                    if (nno.eq.maxno) then
                      maxno = maxno + 100
                      call re_alloc( iano, 1, maxno, name='iano',
     &                         copy=.true., routine='hubbard_term' )
                      call re_alloc( iono, 1, maxno, name='iono',
     &                         copy=.true., routine='hubbard_term' )
                      call re_alloc( xno, 1, 3, 1, maxno, name='xno', 
     &                         copy=.true., routine='hubbard_term' )
                      call re_alloc( Ski,1, maxldau, 1,maxno,name='Ski',       
     &                         copy=.true., routine='hubbard_term' )
                      call re_alloc( grSki, 1, 3, 1, maxldau, 1, maxno,
     &                         name='grSki', routine='hubbard_term', 
     &                         copy=.true. )
                    endif
                    nno = nno + 1
                    iono(nno) = io
                    iano(nno) = ia
                    do ix = 1,3
                      xno(ix,nno) = xki(ix,ina)
                    enddo
                    do ko = 1, nldauproj
                      ioa = iphorb(io)
                      kg  = ldau_gindex(ks,ko)
                      ig  = orb_gindex(is,ioa)
                      call new_matel( 'S', kg, ig, xki(1:3,ina),
     &                    Ski(ko,nno), grSki(1:3,ko,nno) )
!!                     For debugging
!                      write(6,'(4i5,7f10.5)')
!     &                  ka, ko, io, ioa, xki(1:3,ina), 
!     &                  Ski(ko,nno), grSki(1:3,ko,nno)  
!!                     End debugging
                    enddo

                  endif   ! If on orbitals within range

                endif     ! If on orbitals within the local node

              enddo       ! end loop on neighbour orbitals

            enddo         ! end loop on neighbour atoms


!           Loop on neighbour orbitals
            do ino = 1,nno
              io = iono(ino)
              ia = iano(ino)

              call GlobalToLocalOrb(io,Node,Nodes,iio)
              if (iio.gt.0) then
!             Valid orbital

                if (ia .le. nua) then

!                 Scatter density matrix row of orbital io
                  do j = 1,numd(iio)
                    ind = listdptr(iio)+j
                    jo = listd(ind)
                    do ispin = 1,nspin
                      Di(jo,ispin) = Di(jo,ispin) + 
     &                 Dscf(ind,ispin)
                    enddo
                  enddo

!                 Scatter filter of desired matrix elements
                  do j = 1,numh(iio)
                    jo = listh(listhptr(iio)+j)
                    listed(jo) = .true.
                  enddo

!                 Find matrix elements with other neighbour orbitals
                  do jno = 1,nno
                    jo = iono(jno)
                    ja = iano(jno)

                    if ( listed(jo) ) then
!                     Loop on LDA+U projectors
                      do ko = 1, nldauproj
                        lko      = spp%pjldau_l(ko)
                        nprin_ko = spp%pjldau_n(ko)
                        ldaushell => 
     &                     basp%tmp_ldaushell(spp%pjldau_index(ko))      
                        Ueff     = ldaushell%U - ldaushell%J
                        if(ldau_shift) Ueff = 2.0_dp * Ueff
                        Sik = Ski(ko,ino)
                        Sjk = Ski(ko,jno)
!!                       For debugging
!                        write(6,'(a,7i5,3f12.5)')
!     &                  'ka, ko, lko, nprin_ko, ino, jno, jo, U_eff = ',     
!     &                   ka, ko, lko, nprin_ko, ino, jno, jo, Ueff, 
!     &                   Sik, Sjk
!!                       End debugging
                        do ispin=1, nspin
                          Vi(jo,ispin) = Vi(jo,ispin) +
     &                                   0.5_dp * Sik * Sjk * Ueff
                          Dij  = Di(jo,ispin)
                          Cijk = Ueff * Dij * Sjk 
                          if(.not. matrix_elements_only) then 
                            do ix=1,3
                              fik = Cijk * grSki(ix,ko,ino)
                              fa(ix,ia)  = fa(ix,ia)  - fik
                              fa(ix,kua) = fa(ix,kua) + fik
                              do jx = 1, 3
                                stress(jx,ix) = stress(jx,ix) +
     &                                          xno(jx,ino)*fik/volume       
                              enddo 
                            enddo 
                          endif

                          if( .not. ldau_shift ) then 
                            do ikb = 1, nldauproj
                              lkb      = spp%pjldau_l(ikb)
                              nprin_kb = spp%pjldau_n(ikb)
                              if( lko      .eq. lkb        .and.
     &                            nprin_ko .eq. nprin_kb ) then 
!                                For the time being we will use the formulation
!                                of Dudarev and collaborators
                                 Cijk = occu(ko,ikb,kua,ispin) *
     &                                  Ueff * Ski(ikb,jno)
                                 Vi(jo,ispin) = Vi(jo,ispin) -
     &                                          Sik * Cijk
                                if(.not. matrix_elements_only) then 
                                  do ix=1,3
                                    fik = -2.0_dp * Dij * Cijk * 
     &                                    grSki(ix,ko,ino) 
                                    fa(ix,ia)  = fa(ix,ia)  - fik
                                    fa(ix,kua) = fa(ix,kua) + fik
                                    do jx = 1, 3
                                      stress(jx,ix)= stress(jx,ix) +
     &                                        xno(jx,ino) * fik / volume
                                    enddo
                                  enddo 
                                endif
                              endif 
                            enddo 
                          endif

                        enddo   ! Enddo ispin
                      enddo     ! Enddo ko, on LDA+U projectors
                    endif       ! Endif if the orbital jo is listed
                  enddo         ! Enddo on neighbour orbitals jno

!                 Restore  Di and listed
                  do j = 1,numh(iio)
                    ind = listhptr(iio)+j
                    jo = listh(ind)
                    do ispin = 1,nspin
                      Hldau(ind,ispin) = Hldau(ind,ispin)+
     &                                   Vi(jo,ispin)
                      Vi(jo,ispin) = 0.0_dp
                    enddo 
                    listed(jo) = .false.
                  enddo

                  do j = 1,numd(iio)
                    jo = listd(listdptr(iio)+j)
                    do ispin=1,nspin
                      Di(jo,ispin) = 0.0_dp
                    enddo
                  enddo

                endif    ! endif (ia .le. nua) 
              endif      ! endif (iio.gt.0) 
            enddo        ! enddo ino = 1,nno
          enddo          ! enddo loop on orbitals with LDA+U projectors(ka=1,na)

        endif recalc_hamilt

      endif occupations

      Eldau  = 0.0_dp
      DEldau = 0.0_dp
      add_hamilt: if(iter_ldaupop.gt.0) then 
        if(Node.eq.0) 
     &    write(6,'(a)') 'hubbard_term: Adding Hamiltonian contribution'

        do ia = 1,nua
          is = isa(ia)
          spp => species(is)
          nldauproj = spp%nprojsldau
          basp => basis_parameters(is)
          do io = lasto(ia-1)+1,lasto(ia)
            call GlobalToLocalOrb(io,Node,Nodes,iio)
            if ( iio .gt. 0 ) then
              do j = 1, numh(iio)
                ind = listhptr(iio) + j
                do ispin = 1, nspin
                   Eldau = Eldau +
     &                     0.5_dp * Hldau(ind,ispin) * Dscf(ind,ispin)       
                   H(ind,ispin) = H(ind,ispin) + Hldau(ind,ispin)
                enddo
              enddo
            endif
          enddo
          if( nldauproj .gt. 0 ) then
            do ispin = 1, nspin
              do ko = 1, nldauproj
                ldaushell => 
     &             basp%tmp_ldaushell(spp%pjldau_index(ko))      
                Ueff     = ldaushell%U - ldaushell%J
                if( ldau_shift ) Ueff = 2.0_dp * Ueff
                  DEldau = DEldau +
     &                     0.25_dp * (3.0_dp-dble(nspin)) * Ueff *
     &                     occu(ko,ko,ia,ispin)
              enddo
            enddo   
          endif
        enddo

!!       For debugging
!        do io = 1, nuo
!          do j = 1, numh(io)
!              ind = listhptr(io) + j
!              write(6,'(a,4i7,2f12.5)')'io, j, ind, listh, Hldau = ',
!     &              io, j, ind, listh(ind), Hldau(ind,:)
!          enddo 
!        enddo 
!        write(6,'(a,f12.6)') 
!     &    'hubbard_term: Eldau(eV)=',Eldau*13.6057_dp
!        write(6,'(a,f12.6)') 
!     &    'hubbard_term: DEldau(eV)=',DEldau*13.6057_dp      
!        if( .not. matrix_elements_only) then 
!          do io = 1, nua
!            write(6,'(a,i7,3f12.5)')'ia, fa = ',
!     &        io, fa(:,io)
!          enddo
!          do io = 1, 3
!            write(6,'(a,3f12.5)')'stress = ',
!     &        stress(io,:)
!          enddo
!        endif
!!       End debugging
      endif add_hamilt

!      Deallocate local memory

      call de_alloc( grSki, name='grSki' )
      call de_alloc( Ski, name='Ski' )
      call de_alloc( xno, name='xno' )
      call de_alloc( iono, name='iono' )
      call de_alloc( iano, name='iano' )

      call de_alloc( listedall, name='listedall' )
      call de_alloc( listed, name='listed' )
      call de_alloc( Vi, name='Vi' )
      call de_alloc( Di, name='Di' )

!      call de_alloc( occu,     name = 'occu'     )
!      call memory( 'D', 'D', size(occu), 'hubbard_term' )
!      call de_alloc( occu_old, name = 'occu_old' )
!      call memory( 'D', 'D', size(occu_old), 'hubbard_term' )

!     Stop time counter
      call timer( 'hubbard_term', 2 )

      end subroutine hubbard_term
!
      end module m_ldau
