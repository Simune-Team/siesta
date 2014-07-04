module m_ts_io_version

  use precision, only : dp
  implicit none

contains

  subroutine write_TSHS_0(filename, &
       onlyS, Gamma, TSGamma, &
       ucell, na_u, no_l, no_u, no_s, maxnh, nspin,  &
       kscell, kdispl, &
       xa, iza, lasto, &
       numh, listhptr, listh, xij, indxuo, &
       H, S, Ef, &
       Qtot, Temp, &
       istep, ia1)

! *********************************************************************
! Saves the hamiltonian and overlap matrices, and other data required
! to obtain the bands and density of states
! Writen by J.Soler July 1997.
! Note because of the new more compact method of storing H and S
! this routine is NOT backwards compatible
! Modified by M.Paulsson 2009 to:
! 1: To include information of which FC step for phonon calculations
! 2: To only save the overlap matrix if onlyS flag is set
!    (Used for e-ph coupling calculations)
! 3: File format changed to unify Copenhagen/Barcelona Transiesta vers.
! 4: Smaller files by writing arrays directly instead of element wise
! *************************** INPUT **********************************
! logical       Gamma         : Is only gamma point used?
! logical       TSGamma       : Is only TS gamma point used?
! logical       onlyS         : Should only overlap matrix be saved?
! ******************** INPUT or OUTPUT (depending on task) ***********
! integer no_u                : Number of basis orbitals per unit cell
! integer no_s                : Number of basis orbitals per supercell
! integer Enspin              : Spin polarization (1 or 2)
! integer indxuo(no_s)        : Index of orbitals in supercell
! integer maxnh               : First dimension of listh, H, S and
!                               second of xij
! integer numh(nuo)           : Number of nonzero elements of each row
!                               of hamiltonian matrix
! integer listhptr(nuo)       : Pointer to the start of each row (-1)
!                               of hamiltonian matrix
! integer listh(maxnh)        : Nonzero hamiltonian-matrix element column
!                               indexes for each matrix row
! real*8  H(maxnh,Enspin)     : Hamiltonian in sparse form
! real*8  S(maxnh)            : Overlap in sparse form
! real*8  qtot                : Total number of electrons
! real*8  temp                : Electronic temperature for Fermi smearing
! real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
!                               (not read/written if only gamma point)
! TSS Begin
! ********************* ADDED ARGUMENTS FOR TRANSIESTA ****************
! integer fnlength            : file name length
! character(fnlength)  fname  : file nema for input or output
! integer na_u                : Number of atoms per unit cell
! integer istep, ia1          : Force constant step and atom number
! logical check_kcell        : Do a check of the kscell and kdispl
! TSS End
! *************************** UNITS ***********************************
! Units should be consistent between task='read' and 'write'
! *********************************************************************

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: filename
    logical, intent(in) :: onlyS
    logical, intent(in) :: Gamma, TSGamma
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in) :: na_u, no_l, no_u, no_s, maxnh, nspin
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in) :: iza(na_u)
    integer, intent(in) :: numh(no_l), listhptr(no_l)
    integer, intent(in) :: listh(maxnh)
    real(dp), intent(in) :: xij(3,maxnh)
    integer, intent(in) :: indxuo(no_s)
    integer, intent(in) :: lasto(0:na_u)
    real(dp), intent(in) :: H(maxnh,nspin), S(maxnh)
    real(dp), intent(in) :: Ef
    integer, intent(in) :: kscell(3,3)
    real(dp), intent(in) :: kdispl(3)
    real(dp), intent(in) :: Qtot,Temp
    integer, intent(in) :: istep, ia1
    
! ************************
! * LOCAL variables      *
! ************************
    integer :: iu
    integer :: ispin, i,j

    external :: io_assign, io_close

    ! Open file
    call io_assign( iu )
    open( iu, file=filename, form='unformatted', status='unknown' )

    ! Write Dimensions Information
    write(iu) na_u, no_u, no_s, nspin, maxnh
       
    ! Write Geometry information
    write(iu) xa
    write(iu) iza
    write(iu) ucell
    
    ! Write k-point samplung information
    write(iu) Gamma
    ! Is this only an S containing object?
    write(iu) onlyS
    write(iu) TSGamma
    write(iu) kscell
    write(iu) kdispl
    write(iu) istep, ia1
    
    write(iu) lasto
    
    if ( .not. Gamma ) then
       write(iu) indxuo
    endif

    write(iu) numh

    ! Write Electronic Structure Information
    write(iu) Qtot,Temp
    write(iu) Ef
    
    ! Write listh
    do i = 1 , no_u
       write(iu) listh(listhptr(i)+1:listhptr(i)+numh(i))
    end do

    ! Write Overlap matrix
    do i = 1 , no_u
       write(iu) S(listhptr(i)+1:listhptr(i)+numh(i))
    end do
    
    if ( .not. onlyS ) then
       ! Write Hamiltonian 
       do ispin = 1 , nspin 
          do i = 1 , no_u
             write(iu) H(listhptr(i)+1:listhptr(i)+numh(i),ispin)
          end do
       end do
    end if

    if ( .not. Gamma ) then
       do i = 1 , no_u
          write(iu) (xij(j,listhptr(i)+1:listhptr(i)+numh(i)),j=1,3)
       end do
    end if

    ! Close file
    call io_close( iu )

  end subroutine write_TSHS_0

end module m_ts_io_version
    

