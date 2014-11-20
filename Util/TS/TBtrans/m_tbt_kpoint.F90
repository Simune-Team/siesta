module m_tbt_kpoint
  
  use precision, only : dp

  implicit none

  private
  save

  logical, public   :: Gamma
  integer, public   :: nkpnt             ! Total number of k-points

  real(dp), pointer, public :: kweight(:) 
  real(dp), pointer, public :: kpoint(:,:)

  integer,  public, dimension(3,3)  :: kscell = 0
  real(dp), public, dimension(3)    :: kdispl = 0.0_dp

  public :: setup_kscell, setup_kpoint_grid
  public :: write_k_points

contains

  subroutine setup_kscell( cell, N_Elec, Elecs )

! ***************** INPUT **********************************************
! real*8  cell(3,3)  : Unit cell vectors in real space cell(ixyz,ivec)
! ***************** OUTPUT *********************************************
! logical firm_displ   : User-specified displacements (firm)?

!   The relevant fdf labels are kgrid_cutoff and kgrid_Monkhorst_Pack.
!   If both are present, kgrid_Monkhorst_Pack has priority. If none is
!   present, the cutoff default is zero, producing only the gamma point.
!   Examples of fdf data specifications:
!     kgrid_cutoff  50. Bohr
!     %block kgrid_Monkhorst_Pack  # Defines kscell and kdispl
!     4  0  0   0.50               # (kscell(i,1),i=1,3), kdispl(1)
!     0  4  0   0.50               # (kscell(i,2),i=1,3), kdispl(2)
!     0  0  4   0.50               # (kscell(i,3),i=1,3), kdispl(3)
!     %endblock kgrid_Monkhorst_Pack
! **********************************************************************

    use fdf
    use sys,        only : die

    use m_ts_electype
    use intrinsic_missing, only : SPC_PROJ, VEC_PROJ, VNORM
    
    implicit          none

    ! Passed variables
    real(dp), intent(in) :: cell(3,3)
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)

    ! Local variables
    integer :: i, j

    ! For finding projected directions..
    real(dp) :: p(3), contrib

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
    logical :: block_found

    block_found = fdf_block('TBT.kgrid_Monkhorst_Pack',bfdf)
    if ( .not. block_found ) then
       block_found = fdf_block('kgrid_Monkhorst_Pack',bfdf)
    end if

    if ( block_found ) then

       do i = 1 , 3

          if ( fdf_bline(bfdf,pline) ) then
             kscell(1,i) = fdf_bintegers(pline,1)
             kscell(2,i) = fdf_bintegers(pline,2)
             kscell(3,i) = fdf_bintegers(pline,3)
             kdispl(i)   = fdf_breals(pline,1)
          else
             call die( 'setup_tbt_kscell: ERROR no data in' // &
                  'kgrid_Monkhorst_Pack block' )
          end if

       end do

    else

       ! We only calculate the gamma point transport
       kscell(:,:) = 0
       do i = 1 , 3
          kscell(i,i) = 1
       end do
       kdispl(:) = 0._dp

    end if

    do j = 1 , size(Elecs)
       ! project the electrode transport direction onto
       ! the corresponding unit-cell direction
       p = SPC_PROJ(cell,Elecs(j)%ucell(:,Elecs(j)%t_dir))
       ! See which unitcell direction has the highest contribution
       do i = 1 , 3 
          ! project the unit-cell vector onto each cell component
          contrib = VNORM(VEC_PROJ(cell(:,i),p))
          if ( contrib > 1.e-7_dp ) then ! TODO electrode k-points
             ! the contribution along this vector is too much
             ! to disregard the elongation along this
             ! direction.
             ! We *MUST* kill all k-points in this direction
             kscell(:,i) = 0
             kscell(i,:) = 0
             kscell(i,i) = 1
             kdispl(i)   = 0._dp
          end if
       end do
    end do
    
  end subroutine setup_kscell
  
  subroutine setup_kpoint_grid( cell , N_Elec, Elecs )
    
    use fdf, only       : fdf_get, leqi
    use m_find_kgrid, only : find_kgrid
    use parallel, only  : IONode

    use m_ts_electype

    ! Local Variables
    real(dp), intent(in)   :: cell(3,3)
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)

    real(dp) :: tmp, bkpt(3)
    integer :: i
    ! Whether we should apply time-reversal symmetry
    logical :: TRS
    character(len=250) :: user_kfile

    nullify(kweight,kpoint)

    ! In case the user requests to utilize a provided
    ! k-point file
    user_kfile = fdf_get('TBT.k.File','NONE')
    if ( .not. leqi(user_kfile,'NONE') ) then

       if ( IONode ) then
          write(*,'(a)') 'tbtrans: Reading user specified k-points.'
          write(*,'(2a)')'tbtrans: k-points found in file: ',trim(user_kfile)
          write(*,'(2a)')'tbtrans: *** Responsibility is on your side! ***'
       end if

       call tbt_iokp_read(user_kfile,nkpnt,kpoint,kweight)

       ! Convert the k-points to local units
       do i = 1 , nkpnt
          bkpt(:) = kpoint(:,i)
          call kpoint_convert(cell,bkpt(:),kpoint(:,i),-1)
       end do

    else    

       call setup_kscell(cell, N_Elec, Elecs)
       
       TRS = .not. fdf_get('SpinSpiral',.false.)
       TRS = fdf_get('TBT.Symmetry.TimeReversal',TRS)
       
       call find_kgrid(cell, kscell, kdispl, .true., &
            TRS, &
            nkpnt, kpoint, kweight, tmp)

    end if
    
    Gamma = (nkpnt == 1) .and. &
         dot_product(kpoint(:,1),kpoint(:,1)) < 1.0e-20_dp
    
    if (IONode) call write_k_points()
    
  end subroutine setup_kpoint_grid
  
  subroutine write_k_points()
    
    integer  :: ik, ix, i
    
    write(*,'(/,a)') 'tbtrans: k-point coordinates (Bohr**-1) and weights:'
    write(*,'(a,i4,3f12.6,3x,f12.6)') &
         ('tbtrans: ', ik, (kpoint(ix,ik),ix=1,3), kweight(ik), &
         ik=1,nkpnt)

    ! Always write the TranSIESTA k-points
    call tbt_iokp( nkpnt, kpoint, kweight )

    write(*,'(/,a,i0)')  'tbtrans: k-grid: Number of transport k-points = ', nkpnt
    write(*,'(a)') 'tbtrans: k-grid: Supercell and displacements'
    do ix = 1 , 3
       write(*,'(a,3i4,3x,f8.3)') 'tbtrans: k-grid: ',        &
            (kscell(i,ix),i=1,3), kdispl(ix)
    end do

  end subroutine write_k_points
  
  subroutine tbt_iokp( nk, points, weight )
! *******************************************************************
! Saves tbtrans k-points (only writing) Bohr^-1
! Emilio Artacho, Feb. 1999
! Modified by Nick Papior Andersen to not overwrite the SIESTA/TranSIESTA k-points
! ********** INPUT **************************************************
! integer nk           : Number of k-points
! real*8  points(3,nk) : k-point coordinates
! real*8  weight(3,nk) : k-point weight
! *******************************************************************
    use files, only : slabel

    integer  :: nk
    real(dp) :: points(3,nk), weight(nk)
    external :: io_assign, io_close

! Internal 
    integer :: iu, ik, ix

    call io_assign( iu )
    open( iu, file=trim(slabel)//'.TBTKP', form='formatted', status='unknown' ) 

    write(iu,'(i6)') nk
    write(iu,'(i6,3f12.6,3x,f12.6)') &
         (ik, (points(ix,ik),ix=1,3), weight(ik), ik=1,nk)

    call io_close( iu )
    
  end subroutine tbt_iokp

  ! The user can specify their own k-points
  subroutine tbt_iokp_read(fname,nkpt,kpt,wkpt)
    use parallel, only : Node
#ifdef MPI
    use mpi_siesta
#endif
    character(len=*), intent(in) :: fname
    integer, intent(inout) :: nkpt
    real(dp), pointer :: kpt(:,:), wkpt(:)

    real(dp) :: wsum

#ifdef MPI
    integer :: MPIerror
#endif

    ! The user has requested to read in 
    ! k-points from a specific file...
    ! We will do that for him/her.

    external :: io_assign, io_close

    integer :: iu, ik, ix, stat

    if ( Node == 0 ) then
       call io_assign( iu )
       open( iu, file=trim(fname), form='formatted', status='old' ) 

       ! Read number of k-points
       read(iu,*,iostat=stat) nkpt
       call kill_iokp(stat,0)

    end if

#ifdef MPI
    call MPI_Bcast(nkpt,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
#endif

    nullify(kpt,wkpt)
    allocate(kpt(3,nkpt),wkpt(nkpt))

    if ( Node == 0 ) then
       ! Read in the k-points
       wsum = 0._dp
       do ik = 1 , nkpt
          ! read current k-point
          read(iu,*,iostat=stat) ix, kpt(:,ik), wkpt(ik) ! (i6,3f12.6,3x,f12.6)
          call kill_iokp(stat,ik)
          wsum = wsum + wkpt(ik)
       end do
      
       if ( abs(wsum - 1._dp) > 1.e-7_dp ) then
          call die('TBT user specified k-grid does not sum weights to 1.')
       end if

       call io_close( iu )

    end if

#ifdef MPI
    call MPI_Bcast(nkpt,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
#endif

  contains
    
    subroutine kill_iokp(stat,line)
      integer, intent(in) :: stat, line
      if ( stat == 0 ) return
      write(*,*) 'TBtrans iokp could not read your input file'
      write(*,*) 'The k-points MUST be in units of reciprocal vectors!'
      write(*,*) 'TBtrans will convert the unit to correct units.'
      write(*,*) 'Also the sum of weights MUST equal 1.'
      write(*,*) !
      if ( line == 0 ) then
         write(*,*) 'Error occured on reading number of k-points (first line)'
      else
         write(*,'(a,i0,a)') 'Error occured on reading the ',line,' kpoint.'
      end if
      write(*,*) 'Please format your file like this:'
      write(*,*) ' $> cat '//trim(fname)
      write(*,*) ' <nkpt>'
      write(*,*) '     1  <kpt-A1> <kpt-A2> <kpt-A3> <w-kpt>'
      write(*,*) '     2  <kpt-A1> <kpt-A2> <kpt-A3> <w-kpt>'
      write(*,*) ' ....'
      write(*,*) ' <nkpt> <kpt-A1> <kpt-A2> <kpt-A3> <w-kpt>'
      
      call die('TBT reading user specified k-points')
      
    end subroutine kill_iokp
    
  end subroutine tbt_iokp_read
  
end module m_tbt_kpoint
