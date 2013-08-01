
subroutine writemmn( ispin )
!
! In this subroutine, the overlap matrices between the periodic parts of the
! wavefunctions at neighbour k-points are written,
! $M_{m n} ^{( \vec{k} + \vec{b} )}$, see Eq. (27) of the paper by
! N. Marzari et al., Review of Modern Physics 84, 1419 (2012).
! 
! The structure of the outfile (seedname.mmn) is described below,
! according to the rules explained in page 45 of the Wannier90 User Guide,
! version 1.2
!
! Written by J. Junquera in July 2013, 
! based on a previous subroutine by R. Korytar
! 

  use m_siesta2wannier90, only: seedname ! Seed for the name of the file 
                                         !   where the Wannier90
                                         !   code, when used as a postprocessing
                                         !   tool, dumps the information.
  use m_siesta2wannier90, only: numbands ! Number of bands for wannierization
  use m_siesta2wannier90, only: numkpoints ! Total number of k-points
                                         !   for which the overlap of the
                                         !   periodic part of the wavefunct
                                         !   with a neighbour k-point will
                                         !   be computed
  use m_siesta2wannier90, only: nncount  ! Number of neighbour k-points
  use m_siesta2wannier90, only: nnlist   ! nnlist(ikp,inn) is the index of the 
                                         !   inn-neighbour of ikp-point
                                         !   in the Monkhorst-Pack grid 
                                         !   folded to the 
                                         !   first Brillouin zone
  use m_siesta2wannier90, only: nnfolding ! nnfolding(i,ikp,inn) is the 
                                         !   i-component 
                                         !   of the reciprocal lattice vector 
                                         !   (in reduced units) that brings
                                         !   the inn-neighbour specified in 
                                         !   nnlist (which is in the first BZ)  
                                         !   to the actual \vec{k} + \vec{b} 
                                         !   that we need.
                                         !   In reciprocal lattice units.
  use m_siesta2wannier90, only: Mmnkb    ! Matrix of the overlaps of 
                                         !   periodic parts of Bloch waves.
                                         !   <u_{m k}|u_{n k+b}>
                                         !   The first two indices refer to 
                                         !   the number of occupied bands
                                         !   (indices m and n in standard
                                         !   notation, see for instance, 
                                         !   Eq. (27) of the paper by 
                                         !   Marzari et al., RMP 84, 1419 (2012)
                                         !   The third index refer to the kpoint
                                         !   The fourth index refer to the neig

  implicit none

  integer, intent(in)     :: ispin      ! Spin component

!
! Local variables
!

  character(len=len_trim(seedname)+4) :: mmnfilename ! Name of the file where
                                                     !   the overlap Mmn
                                                     !   matrices will be
                                                     !   written

  integer      :: ik        ! Counter for the loop on k-points
  integer      :: inn       ! Counter for the loop on beighbours
  integer      :: nband     ! Counter for the loop on bands
  integer      :: mband     ! Counter for the loop on bands
  integer      :: mmnunit   ! Logical unit assigned to the file where the 
                            !   Mmn matrices will be written
  integer      :: g(3)      ! Auxiliary variable to write the folding vector
  integer      :: eof

  external     :: io_assign ! Assign a logical unit 
  external     :: io_close  ! Close a logical unit


! Asign a logical unit to the file where the Mmn matrices will be written
  call io_assign( mmnunit )

! Asign a name to the file where the Mmn matrices will be written
  mmnfilename = trim( seedname )// ".mmn"

! Open the output file where the Mmn matrices will be written
  open( unit=mmnunit, err=1984, file=mmnfilename, status="replace", &
 &      iostat=eof)

! The first line of the mmn file is a user comment
  write( unit=mmnunit, fmt="(a38)", err=1984 )                      &
 &  "siesta2wannier90: Siesta for Wannier90"

! The second line includes three integers:
!   numbands:   the number of bands for wannierization
!   numkpoints: the number of k-points for which the overlap of the
!      periodic part of the wavefunction with a neighbour k-point will be 
!      computed
!   nncount:    the number of neighbour k-points 
  write( unit=mmnunit, fmt="(i5,x,i5,x,i2)", err=1984 )             &
 &  numbands(ispin), numkpoints, nncount
! &  numIncBands(ispin),numKPoints,nnCount

! Then, there are numkpoints x nncount blocks of data.
  do ik = 1,numkpoints
    do inn = 1,nncount
   ! override type for the sake of formatting
      g(:) = nnfolding(:,ik,inn)

!     The first line of each block contains five integers:
!     The first specifies the k-points (i.e. gives the ordinal correspondence
!        to its position in the list of k-points in seedname.win).
!     The second to the fifth integers specify \vec{k} + \vec{b}
!     The second integer, in particular, points to the k-point on the list
!        that is a periodic image of \vec{k} + \vec{b}, and in particular
!        is the image that is actually mentioned in the list.
!     The last three integers specify the \vec{G} vector, 
!        in reciprocal lattice units, that brings the k-point specified by
!        the second integer, and that this lives inside the first BZ zone,
!        to the actual \vec{k} + \vec{b} that we need.

      write( unit=mmnunit, fmt="(i5,i5,3x,3i4)", err=1984 )         &
 &      ik, nnlist(ik,inn), g

!     Subsequent numbands x numbands lines of each block: 
!     two real numbers per line. 
!     These are the real and imaginary parts, respectively, of the actual
!     scalar product $M_{m n} ^{( \vec{k} + \vec{b} )} for
!     $m, n \in [1, numbands].
!     The order of these elements is such that the first index $m$ is fastest

!      do nband = 1,numIncBands(ispin)
!        do mband = 1,numIncBands(ispin)
      do nband = 1, numbands(ispin)
        do mband = 1,numbands(ispin)
        write( unit=mmnunit, fmt="(f12.5,2x,f12.5)", err=1984 )     &
 &         real(  Mmnkb(mband,nband,ik,inn) ),                      &
 &         aimag( Mmnkb(mband,nband,ik,inn) )
        enddo  ! End loop on bands      (mband)
      enddo    ! End loop on bands      (nband)
    enddo      ! End loop on neighbours (inn)
  enddo        ! End loop on k-points   (ik) 

! Close the output file where the Mmn matrices will be written
  call io_close(mmnunit)

  write(6,'(/,a)')  &
 &  'mmn: Overlap matrices between periodic part of wavefunctions'
  write(6,'(a)')  &
 &  'mmn: written in ' // trim(seedname) // '.mmn file'

  return

1984 call die("writemmn: Error writing to .mmn file")

end subroutine writemmn

!
! ------------------------------------------------------------------------------
!

subroutine writeamn( ispin )
!
! In this subroutine, the overlap matrices between the trial projection
! functions and the Hamiltonian eigenstates are written,
! $A_{m n}^{(\vec{k})} = 
!    \langle \psi_{m \vec{k}} \vbar g_{n} \rangle 
! (see paragraph after Eq. (16) of the review by
! N. Marzari et al., Review of Modern Physics 84, 1419 (2012),
! or Eq. (1.8) of the Wannier Users Guide, Version 1.2
! 
! The structure of the outfile (seedname.amn) is described below,
! according to the rules explained in page 45 of the Wannier90 User Guide,
! version 1.2
!
! Written by J. Junquera in July 2013, 
! based on a previous subroutine by R. Korytar
! 

  use m_siesta2wannier90, only: seedname ! Seed for the name of the file 
                                         !   where the Wannier90
                                         !   code, when used as a postprocessing
                                         !   tool, dumps the information.
  use m_siesta2wannier90, only: numbands ! Number of bands for wannierization
  use m_siesta2wannier90, only: numkpoints ! Total number of k-points
                                         !   for which the overlap of the
                                         !   periodic part of the wavefunct
                                         !   with a neighbour k-point will
                                         !   be computed
  use m_siesta2wannier90, only: numproj  ! Total number of projection centers,
                                         !   equal to the number of MLWF
  use m_siesta2wannier90, only: Amnmat   ! Projections of a trial function
                                         !   with a Bloch orbital
                                         !   <\psi_{m k}|g_n>

!
! Local variables
!

  character(len=len_trim(seedname)+4) :: amnfilename ! Name of the file where
                                                     !   the overlap Amn
                                                     !   matrices will be
                                                     !   written
  integer      :: ik        ! Counter for the loop on k-points
  integer      :: iproj     ! Counter for the loop on projectors
  integer      :: mband     ! Counter for the loop on bands
  integer      :: amnunit   ! Logical unit assigned to the file where the 
                            !   Amn matrices will be written
  integer      :: eof

  external     :: io_assign ! Assign a logical unit 
  external     :: io_close  ! Close a logical unit


! Asign a logical unit to the file where the Amn matrices will be written
  call io_assign( amnunit )

! Asign a name to the file where the Amn matrices will be written
  amnfilename = trim( seedname )// ".amn"

! Open the output file where the Amn matrices will be written
  open( unit=amnunit, err=1992, file=amnfilename, status="replace", &
 &      iostat=eof)

! The first line of the mmn file is a user comment
  write( unit=amnunit, fmt="(a38)", err=1992 )                      &
 &  "siesta2wannier90: Siesta for Wannier90"

! The second line includes three integers:
!   numbands:   the number of bands for wannierization
!   numkpoints: the number of k-points for which the overlap of the
!      periodic part of the wavefunction with a neighbour k-point will be 
!      computed
!   numproj:    the number of projections
  write( unit=amnunit, fmt="(i5,x,i5,x,i2)", err=1992 )             &
 &  numbands(ispin), numkpoints, numproj
! &  numIncBands(ispin),numKPoints,numproj

! Subsequent numbands x numproj x numkpoint lines:  
! three integers and two real numbers on each line
! The first two integers are the band indices m and n.
! The third integer specifies the \vec{k} by giving the ordinal 
! corresponding to its position in the list of k-points in seedname.win
! The real numbers are the real and the imaginary parts,
! respectively, of the actual $A_{m n} (\vec{k})$.

  do ik = 1, numkpoints
    do iproj = 1, numproj
!      do mband = 1,numIncBands(ispin)
      do mband = 1,numbands(ispin)
        write(unit=amnunit,fmt="(3i5,x,f12.5,2x,f12.5)",err=1992)      &
 &         mband, iproj, ik,                                           &
 &         real(Amnmat(mband,iproj,ik)),aimag(Amnmat(mband,iproj,ik))
      enddo
    enddo
  enddo

! Close the output file where the Amn matrices will be written
  call io_close(amnunit)

  write(6,'(/,a)')  &
 &  'amn: Overlap matrices between trial projection functions and wavefunctions'
  write(6,'(a)')  &
 &  'amn: written in ' // trim(seedname) // '.amn file'

  return

1992 call die("writeamn: Error writing to .amn file")

end subroutine writeamn

!
! ------------------------------------------------------------------------------
!

subroutine writeeig( ispin )
!
! In this subroutine, the Kohn-Sham eigenvalues $\epsilon_{n \vec{k}}$ 
! (in eV) at each point in the Monkhorst-Pack mesh are written
!
! This is required if any of disentanglement, plot_bands, plot_fermi_surface
! or hr_plot options are set to true in Wannier90.
! 
! The structure of the outfile (seedname.eig) is described below,
! according to the rules explained in page 45 of the Wannier90 User Guide,
! version 1.2
!
! Written by J. Junquera in July 2013, 
! based on a previous subroutine by R. Korytar
! 

  use m_siesta2wannier90, only: seedname  ! Seed for the name of the file 
                                          !   where the Wannier90
                                          !   code, when used as postprocessing
                                          !   tool, dumps the information.
  use m_siesta2wannier90, only: numkpoints! Total number of k-points
                                          !   for which the overlap of the
                                          !   periodic part of the wavefunct
                                          !   with a neighbour k-point will
                                          !   be computed
  use m_siesta2wannier90, only: numbands  ! Number of bands for wannierization
  use m_siesta2wannier90, only: eo        ! Eigenvalues of the Hamiltonian 
                                          !   at the numkpoints introduced in
                                          !   kpointsfrac 
                                          !   First  index: band index
                                          !   Second index: k-point index

  implicit none

  integer, intent(in)     :: ispin       ! Spin component

!
! Local variables
!

  character(len=len_trim(seedname)+4) :: eigfilename ! Name of the file where
                                                     !   the eigenvalues
                                                     !   written

  integer      :: ik        ! Counter for the loop on k-points
  integer      :: iband     ! Counter for the loop on bands
  integer      :: eigunit   ! Logical unit assigned to the file where the 
                            !   eigenvalues will be written
  integer      :: eof       ! Status of the output file

  external     :: io_assign ! Assign a logical unit 
  external     :: io_close  ! Close a logical unit

! Asign a logical unit to the file where the eignvalues will be written
  call io_assign(eigunit)

! Asign a name to the file where the eigenvalues will be written
  eigfilename = trim(seedname)// ".eig"

! Open the output file where the eigenvalues will be written
  open( unit=eigunit, err=1983, file=eigfilename, status="replace", &
 &      iostat=eof )

! Each line consists of two integers and a real number.
! The first integer is the band index, 
! The second integer gives the ordinal corresponding to the k-point in
!   the list of k-points in seedname.win.
! The real number is the eigenvalue (in eV).
  do ik = 1, numkpoints
!    do iband = 1,numIncBands(ispin)
    do iband = 1, numbands(ispin)
      write( unit=eigunit, fmt="(i5,i5,3x,f12.5)", err=1983 )       &
 &      iband, ik, eo(iband,ik)
    enddo
  enddo

  write(6,'(/,a)')  &
 &  'eig: Eigenvalues of the Hamiltonian '
   write(6,'(a)')  &
 &  'eig: written in ' // trim(seedname) // '.eig file'


! Close the output file where the eigenvalues will be written
  call io_close(eigunit)

  return

  1983 call die("writeeig: Error writing to .eig file")
end subroutine writeeig

!
! ------------------------------------------------------------------------------
!

subroutine writeunk( ispin )
!
! Produces UNKXXXXX.X files which contain the periodic
! part of a Bloch function in the unit cell on a grid given by
! global nx,ny,nz variables. Hamiltonian must be diagonalized
! for every k before the routine is invoked. The coefficients
! have to be stored in coeffs.
!
!@ See wannier_plot.tex, Eq. (4)
use neighbour,        only: maxnna
!use m_denchar_neighb, only: neighb ! Yields list of nonzero orbitals
                                   ! in a given point

! This is copied from constants.f90 in Wannier90-1.1
! We use it due to i/o compatibility of binary UNKXXX files.
  integer, parameter :: wannier90dp = selected_real_kind(15,300)


  integer,intent(in):: ispin
  character(len=11) :: unkfilename

! k point index
  integer           :: kpt

end subroutine writeunk
