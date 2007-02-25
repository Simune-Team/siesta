      Module SIESTA_TODO

      use precision, only: dp
#ifdef MPI
      use mpi_siesta, only: mpi_comm_world
      use m_mpi_utils, only : globalize_sum, globalize_max
#endif
#ifdef CDF
      use md_out,     only: md_netcdf
#endif

      implicit none

!      integer:: j       ! General-purpose index
!      integer:: ia      ! Atom index
!      integer:: ja      ! Atom index  NOT USED
!      integer:: jo      ! Orbital index
!      integer:: jx      ! Cartesian coordinate index
!      integer:: ix      ! Cartesian coordinate index

!      integer:: ind     ! Index of sparse matrix elements
!      integer:: ik      ! k-vector index
!      integer:: iv      ! Vector index. Not used
!      integer:: io      ! Basis-orbital index
      !integer:: iord    ! Order of density matrix extrapolation between geom. steps
!      integer:: is      ! Species index
!      integer:: ispin   ! Spin index
!      integer:: iadispl ! Displaced atom index
!      integer:: ixdispl ! Displaced cartesian coordinate

!      integer:: iu      ! I/O file unit
      !integer:: iunit   ! Physical-units option: 1=>(eV,Ang), 2=>(Ry,Bohr)

!      integer:: ifa     ! Calculate forces? 0=>no, 1=>yes
!      integer:: istr    ! Calculate stress? 0=>no, 1=>yes
!      integer:: ihmat   ! Calculate hamiltonian? 0=>no, 1=>yes


!      integer:: isel    ! Auxiliary to call neighb: avoid double count of atom pairs?
      !integer:: nbcell  ! Auxiliary to call shaper routine
!      integer:: istpsave! Temporary to save and restore istp when reading DM from file

!      integer:: neigmin ! Min. number of eigenstates (per k point)
      !integer:: ntcon   ! Total number of position constraints imposed

!      integer:: maxwk   ! Number of k-points for wavefunction printout
!      integer:: nhist   ! Number of histogram intervals in projected density of states
      !integer:: ntm(3)  ! Number of integration mesh divisions of each cell vector

      ! maxnh is (initially) the minimum size of all the sparse arrays.
      ! It must be 1 (not 0) since we frequently use the F77 idiom of
      ! passing the first element of the array instead of all of it.

      !integer              :: nnia      ! Number of neighbors of a given atom
      !integer              :: maxna=200 ! Max. number of neighbor atoms
      !integer, allocatable :: jna(:)    ! Index of neighbor atoms

      !integer :: nnamax  = 1 ! Max. required number of neighbor atoms
      !integer :: nauxpul = 1 ! Max. required size of auxiliary sparse Pulay matrices

      !integer :: no_l    = 1 ! Number of basis orbitals stored in one node
      !integer :: ns      = 1 ! Number of chemical species. AG** do not initialize
      !integer :: nspin   = 1 ! Number of spin components
      !integer :: nxij    = 1 ! Number of neighbors of a given basis orbital


!      real(dp):: bcell(3,3) ! Reciprocal unit vectors
      !real(dp):: cfmax
      !real(dp):: cftem
!      real(dp):: const      ! Temporary factor (used in calculation of Harris forces)
!      real(dp):: factor     ! Temporary variable
      !real(dp):: dipol(3)   ! Electric dipole vector

      !real(dp):: fmax       ! Maximum atomic force
      !real(dp):: fmean      ! Mean atomic force (total force / number of atoms)
      !real(dp):: fres       ! Root mean squared atomic force

      !real(dp):: ftem       ! Temporary force value
      !real(dp):: ftot(3)    ! Total force

      !real(dp):: g2max      ! Actual planewave cutoff of real-space integration mesh

!      real(dp):: getot      ! Total energy including thermostat dynamic variables

      !real(dp):: kn         ! Kinetic energy of the Nose' thermostat
      !real(dp):: kpr        ! Kinetic energy of the Parrinello-Rahman variables

      !real(dp):: Pint       ! Instantaneous pressure
      !real(dp):: Pmol       ! "Molecular pressure" (discounting Virial term)
      !real(dp):: Psol       ! Solid (periodic system) pressure

!      real(dp):: rmax_bonds ! Cutoff length for bond definition
!      real(dp):: tiny=1.e-15_dp ! Small tolerance used to compare k points
      !real(dp):: tempion=0.0_dp ! Ionic temperature
!      real(dp):: qaux       ! Auxiliary to call spnvec (not used later)
      !real(dp):: qspin(4)   ! Total spin populations
!      real(dp):: qsol       ! Total unnormalized electron charge
      !real(dp):: rmax       ! Max. required interatomic distance
      !real(dp):: rmaxh      ! Max. distance between non-zero hamiltonian matrix orbitals

!      real(dp):: sigma      ! Energy width used to convolute partial density of states

!      real(dp):: stot       ! Total spin polarization
!      real(dp):: svec(3)    ! Total spin vector in non-collinear-spin calculations

!      real(dp):: virial     ! Virial = Sum(x*f)
      !real(dp):: volume     ! Unit cell volume

      !real(dp):: vn         ! Velocity (time derivative) of the Nose thermostat
      !real(dp):: vpr        ! Velocity (time derivative) of the PR variables


      !real(dp), external                  :: volcel

      !real(dp), allocatable               :: auxpul(:,:)! Auxilialy matrices for Pulay mixing
      !real(dp),             pointer, save :: cfa(:,:)   ! Atomic forces orthogonalized to geometry constraints
      !real(dp),             pointer, save :: fa(:,:)    ! Atomic forces
      !real(dp),             pointer, save :: fal(:,:)   ! Local-node part of atomic forces
!      real(dp), allocatable               :: ebk(:,:,:) ! Band energies
      !real(dp), allocatable               :: efs(:)     ! Fermi energy for each spin (for fixed spin calculations)
      !real(dp),             pointer, save :: eo(:,:,:)  ! Hamiltonian eigenvalues
      !real(dp),             pointer, save :: qo(:,:,:)  ! Occupations of eigenstates
      !real(dp), allocatable               :: qs(:)      ! Spin charges
      !real(dp), allocatable               :: r2ij(:)    ! Interatomic distances squared
      !real(dp), allocatable               :: va(:,:)    ! Atomic velocities
      !real(dp), allocatable               :: xij(:,:)   ! Interatomic vectors
      !real(dp), allocatable               :: xijo(:,:)  ! Vectors between basis-orbital centers


#     ifdef MPI
      integer :: MPIerror  ! Return error code in MPI routines
      integer :: ntmp      ! Temporary to call an MPI routine

      real(dp):: buffer1   ! Temporary for the result of a globalize operation
      real(dp):: qtmp(4)   ! Temporary for globalized spin charge
#     endif

      logical:: auxchanged ! Has auxiliary supercell changed?
      logical:: dminit     ! Initialize density matrix?
      logical:: found      ! Found density matrix file?
      logical:: foundxv    ! Found position-velocity file?
      logical:: foundzm    ! Found Z-matrix file?
      logical:: gamma      ! Gamma point (k=0) in BZ only?
      logical:: genlogic   ! Generate local density of states?
      logical:: itest      ! Temporary used to check if DM have just been read
      logical:: mmix       ! Same as mix. Used in pulayx
      logical:: overflow   ! Too short array sizes?
      logical:: overflowed ! Have array sizes been increased?
    
      logical:: eggbox_block=.true. ! Read eggbox info from data file?

      character(len=150) :: line    ! Data line read
      character(len=10)  :: message ! Not used
      character(len=150) :: sname   ! System name, used to initialise read
      character(len=150) :: shape*10


      !character(len=label_length+5), external :: paste

      !type(parsed_line), pointer  :: p      

      END MODULE SIESTA_TODO
