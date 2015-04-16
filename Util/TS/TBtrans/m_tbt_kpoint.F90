module m_tbt_kpoint
  
  use precision, only : dp

  implicit none

  private
  save

  logical, public :: Gamma
  integer, public :: nkpnt             ! Total number of k-points

  real(dp), pointer, public :: kweight(:) 
  real(dp), pointer, public :: kpoint(:,:)

  integer,  public :: kscell(3,3) = 0
  real(dp), public :: kdispl(3) = 0.0_dp

  public :: setup_kpoint_grid
  public :: write_k_points
  public :: read_kgrid
  public :: tbt_iokp

contains

  subroutine read_kgrid(bName,N_Elec,Elecs,TRS,cell,kpt,wkpt, &
       is_b, &
       kcell,kdispl)

    use parallel, only : IONode
    use fdf
    use fdf_extra, only : fdf_bnext
    use intrinsic_missing, only : EYE, VNORM, SPC_PROJ, VEC_PROJ
    use units, only : Pi
    use m_find_kgrid, only : find_kgrid

    use m_ts_electype

    ! INPUT
    character(len=*), intent(in) :: bName
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! Whether time-reversal symmetry applies
    logical, intent(in) :: TRS
    real(dp), intent(in) :: cell(3,3)
    ! OUTPUT
    real(dp), pointer :: kpt(:,:), wkpt(:)
    logical, intent(in), optional :: is_b
    ! Optional OUTPUT
    integer, intent(out), optional :: kcell(3,3)
    real(dp), intent(out), optional :: kdispl(3)

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
    integer :: i, ik, nkpt, iEl
    real(dp) :: rcell(3,3), displ(3), ksize(3), rtmp, p(3)
    real(dp) :: contrib
    integer :: kscell(3,3)
    character(len=50) :: ctmp

    logical :: is_list

    ! Initialize values
    ksize(:) = 1._dp
    displ(:) = 0._dp
    kscell(:,:) = 0
    kscell(1,1) = 1
    kscell(2,2) = 1
    kscell(3,3) = 1

    ! If the block does not exist, simply 
    ! create the Gamma-point
    nullify(kpt,wkpt)

    if ( .not. fdf_block(bName,bfdf) ) then
       ! the block does not exist, hence the user
       ! requests a Gamma-point.
       allocate(kpt(3,1),wkpt(1))
       kpt(:,:) = 0._dp
       wkpt(:) = 1._dp
       if ( present(kcell) ) kcell = kscell
       if ( present(kdispl) ) kdispl = displ
       return
    end if

    ! Read in the blocks
    ik = 0
    is_list = .false.
    do while ( fdf_bnext(bfdf,pline) )
       
       if ( fdf_bnnames(pline) > 0 ) then

          ctmp = fdf_bnames(pline,1)

          ! We have some kind of designation
          if ( leqi(ctmp(1:5),'diag-') ) then

             if ( fdf_bnintegers(pline) /= 1 ) then
                call die('Please correct your input, you have not supplied a number in diag-.')
             end if
             
             ctmp = ctmp(6:)
             if ( leqi(ctmp,'A1') .or. leqi(ctmp,'a') ) then
                ik = 1
             else if ( leqi(ctmp,'A2') .or. leqi(ctmp,'b') ) then
                ik = 2
             else if ( leqi(ctmp,'A3') .or. leqi(ctmp,'c') ) then
                ik = 3
             end if
             
             ! Get the diagonal
             kscell(ik,ik) = max(1,fdf_bnintegers(pline,1))
             
          else if ( leqi(ctmp,'diagonal') .or. &
               leqi(ctmp,'diag') ) then
             
             ik = fdf_bnintegers(pline)
             
             if ( ik < 3 .and. IONode ) then
                write(*,'(a)') 'tbtrans: POSSIBLE WARNING'
                write(*,'(a,i0,a)') 'tbtrans: You have only supplied ', &
                     ik,' of the 3 diagonal k-cell elements.'
                write(*,'(a)') 'tbtrans: Will assume this order A1-A2-A3'
             end if
             
             ! Set the diagonal
             do i = 1 , ik
                kscell(i,i) = max(1,fdf_bintegers(pline,i))
             end do

          else if ( leqi(ctmp,'displacement') .or. &
               leqi(ctmp,'displ') ) then

             displ(1) = fdf_bvalues(pline,1)
             displ(2) = fdf_bvalues(pline,2)
             displ(3) = fdf_bvalues(pline,3)

          else if ( leqi(ctmp,'size') ) then

             ksize(1) = fdf_bvalues(pline,1)
             ksize(2) = fdf_bvalues(pline,2)
             ksize(3) = fdf_bvalues(pline,3)
             if ( any(ksize > 1._dp) ) then
                call die('The size of the Brillouin zone MUST be &
                     &less than or equal to 1.')
             end if

          else if ( leqi(ctmp,'list') ) then

             ! Get number of k-points
             nkpt = fdf_bintegers(pline,1)
             ! allocate for the k-points
             allocate(kpt(3,nkpt),wkpt(nkpt))
             ! reset all weights to be equal
             wkpt(:) = 1._dp / real(nkpt,dp)
             is_list = .false.
             do ik = 1 , nkpt
                if ( .not. fdf_bnext(bfdf,pline) ) &
                     call die('Could not read correct number of k-points in list')

                if ( ik == 1 .and. fdf_bnvalues(pline) > 3 ) is_list = .true.
                kpt(1,ik) = fdf_bvalues(pline,1)
                kpt(2,ik) = fdf_bvalues(pline,2)
                kpt(3,ik) = fdf_bvalues(pline,3)

                if ( is_list ) then
                   if ( .not. fdf_bnvalues(pline) > 3 ) then
                      call die('Could not read weight for k-point, &
                           &either supply all, or none (which means equal weight)')
                   end if
                   wkpt(ik) = fdf_bvalues(pline,4)
                end if

             end do

             ! This tells the remaining algorithm
             ! to not create a Monkhorst-Pack grid.
             is_list = .true.

          end if

       else if ( fdf_bnintegers(pline) == 3 ) then

          ! There exists two variants
          ! 1. Either the user only supplies the diagonal,
          ! 2. or the full kscell is supplied

          do ik = 1 , 3
             kscell(1,ik) = fdf_bintegers(pline,1)
             kscell(2,ik) = fdf_bintegers(pline,2)
             kscell(3,ik) = fdf_bintegers(pline,3)
             if ( fdf_bnvalues(pline) > 3 ) then
                displ(ik) = fdf_bvalues(pline,4)
             end if
             ! To not error out of only 3 lines grids
             if ( ik == 3 ) cycle
             if ( .not. fdf_bnext(bfdf,pline) ) &
                  call die('Could not read kgrid from block: '//trim(bName))
          end do

       end if

    end do

    if ( .not. is_list ) then

       do iEl = 1 , N_Elec
          ! project the electrode transport direction onto
          ! the corresponding unit-cell direction
          p = SPC_PROJ(cell,Elecs(iEl)%ucell(:,Elecs(iEl)%t_dir))
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
                displ(i)   = 0._dp
             end if
          end do
       end do

       if ( present(kcell) ) kcell = kscell
       if ( present(kdispl) ) kdispl = displ

       call EYE(3, rcell, 2._dp * Pi)
       call find_kgrid(rcell, kscell, displ, .true., &
            TRS , &
            nkpt, kpt, wkpt, rtmp)

       ! Re-scale the k-points to the correct size
       do ik = 1 , nkpt
          kpt(:,ik) = kpt(:,ik) * ksize(:)
       end do

       ! Rescale the weights if the size of the
       ! kgrid zone has been narrowed
       ! This is simple, the basic size is 1.
       ! hence if the size of the k-cell is 1/2
       ! then the weights must also be 1/2
       ! I.e. the size is the weight scale.
       do ik = 1 , 3
          ! we do not need to rescale if the number
          ! of k-points in the i'th direction is
          ! 1. (regardless of what the user says ;))
          if ( kscell(ik,ik) == 1 ) cycle
          wkpt(:) = wkpt(:) * ksize(ik)
       end do

    else

       ! Print out warning if the sum of the weights
       ! does not equal 1
       contrib = 0._dp
       do ik = 1 , nkpt
          contrib = contrib + wkpt(ik)
       end do

       if ( IONode .and. abs(contrib - 1._dp) > 1.e-7_dp ) then
          write(*,'(a)')'WARNING: Weights for k-points in &
               & %block '//trim(bName)//' does not sum to 1.'
       end if

    end if

    if ( present(is_b) ) then
       if ( is_b ) return
    end if

    ! Transform the reciprocal units into length
    call reclat(cell,rcell,1)
    do ik = 1 , nkpt
       ksize(:) = kpt(:,ik)
       call kpoint_convert(rcell,ksize(:),kpt(:,ik),-2)
    end do

  end subroutine read_kgrid
  
  subroutine setup_kpoint_grid( cell , N_Elec, Elecs )
    
    use fdf, only       : fdf_get, leqi, block_fdf, fdf_block
    use parallel, only  : IONode

    use m_ts_electype

    ! Local Variables
    real(dp), intent(in)   :: cell(3,3)
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)

    real(dp) :: bkpt(3)
    integer :: i
    type(block_fdf) :: bfdf
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

       TRS = .not. fdf_get('SpinSpiral',.false.)
       TRS = fdf_get('TBT.Symmetry.TimeReversal',TRS)

       if ( fdf_block('TBT.kgrid_Monkhorst_Pack',bfdf) ) then
          call read_kgrid('TBT.kgrid_Monkhorst_Pack', &
               N_Elec,Elecs,TRS,cell,kpoint,kweight, &
               kcell=kscell,kdispl=kdispl)
       else
          call read_kgrid('kgrid_Monkhorst_Pack', &
               N_Elec,Elecs,TRS,cell,kpoint,kweight, &
               kcell=kscell,kdispl=kdispl)
       end if
       nkpnt = size(kweight)
       
    end if
    
    Gamma = (nkpnt == 1) .and. &
         dot_product(kpoint(:,1),kpoint(:,1)) < 1.0e-20_dp
    
    call write_k_points()
    
  end subroutine setup_kpoint_grid
  
  subroutine write_k_points()
    use parallel, only : IONode
    
    integer  :: ik, ix, i

    if ( .not. IONode ) return 
    
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
  
  subroutine tbt_iokp( nk, points, weight , fend)
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

    integer, intent(in) :: nk
    real(dp), intent(in) :: points(3,nk), weight(nk)
    character(len=*), intent(in), optional :: fend
    external :: io_assign, io_close

! Internal 
    integer :: iu, ik, ix

    call io_assign( iu )
    if ( present(fend) ) then
       open( iu, file=trim(slabel)//'.'//trim(fend), &
            form='formatted', status='unknown' ) 
    else
       open( iu, file=trim(slabel)//'.TBTKP', &
            form='formatted', status='unknown' ) 
    end if

    write(iu,'(i6)') nk
    write(iu,'(i6,3f12.6,3x,f12.6)') &
         (ik, (points(ix,ik),ix=1,3), weight(ik), ik=1,nk)

    call io_close( iu )
    
  end subroutine tbt_iokp

  ! The user can specify their own k-points
  subroutine tbt_iokp_read(fname,nkpt,kpt,wkpt)
    use parallel, only : Node
    use m_io_s, only : file_exist
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
       
       if ( .not. file_exist(trim(fname)) ) then
          call die('Could not locate file '//trim(fname)// &
               ' please ensure that the file exists.')
       end if

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
          write(*,'(a)')'WARNING: Weights for user specified k-points does &
               &not sum to 1.'
       end if

       call io_close( iu )

    end if

#ifdef MPI
    call MPI_Bcast(kpt(1,1),3*nkpt,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(wkpt(1),nkpt,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
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
