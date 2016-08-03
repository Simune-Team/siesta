MODULE  wavefunctions
  USE  precision
  USE  MatrixSwitch
  !
  implicit none 
  !    
  type(matrix), allocatable, save :: wavef_ms(:,:)
  ! 
CONTAINS
  !
  subroutine iowavef(task,wavef_rw,nuotot,nk,nspin,istpp,totime)
 
    ! To read and write the TDKS orbitals
    ! 
    ! Written by Daniel Sanchez-Portal
    ! Re-written by Rafi Ullah, CIC nanoGUNE, October 16, 2015 to work
    ! with parallel TDDFT flavor of siesta using MatrixSwtich.
    !
    ! Although it is prepared to work with parallel TDDFT-siesta the
    ! reading/writing process itself is not parallel.

    use parallel
    use fdf
    use MatrixSwitch
#ifdef MPI
    use mpi_siesta
#endif
    !
    implicit none
    !
    integer, intent(in)  :: nuotot, nk, nspin
    integer, intent(inout) :: istpp
    double precision, intent(inout) :: totime
    character*(*), intent(in):: task
    type(matrix), intent(inout)  :: wavef_rw(nk,nspin)
    ! Internal variables and arrays
    character :: fname*33, sname*30, paste*33, m_storage*5
    logical   ::  exist1, frstme
    integer   :: unit1, ie,nuototread,nkread,nspinread,dim2read
    integer   :: nwf, ik, ispin, mxnwf, io, ix,i,j
    external  :: chkdim, io_assign, io_close, paste, timer,memory
    complex(kind=dp) :: varaux
    save      :: frstme, fname
    data         frstme /.true./
    !
#ifdef MPI
    INTEGER :: MPIerror
#endif
    !
#ifdef MPI
    call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
    call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
    Node = 0
    Nodes = 1
#endif
    !
    ! Find file name
    if (frstme) then
      if (Node.eq.0) then
        sname = fdf_string('SystemLabel','siesta')
      endif
#ifdef MPI
      call MPI_Bcast(sname,30,MPI_character,0,MPI_Comm_World,MPIerror)
#endif
#ifdef MPI
      m_storage='pzdbc'
#else
      m_storage='szden'
#endif
      fname = paste(sname,'.TDWF')
      frstme = .false.
    endif
    !
    if (task.eq.'read' .or. task.eq.'READ') then      ! task read
      if (Node.eq.0) then
        inquire (file=fname, exist=exist1)
        if(.not.exist1) then 
          write(6,'(/,a)')                                            &
          'iowavef: file containing the KS orbitals not found'
          write(6,'(a)')                                              &
          'iowavef: simulation can not be started or restarted'
          stop 
        endif
        call io_assign(unit1)
        open( unit1, file=fname, form='unformatted',                 &
        status='unknown' )
        read(unit1) istpp, totime
      endif
#ifdef MPI 
      call MPI_Bcast(istpp,1,MPI_double_complex,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(totime,1,MPI_double_complex,0,MPI_Comm_World,MPIerror)
#endif
      if (Node.eq.0) then
        read(unit1) nuototread, nkread, nspinread
        if(nkread.ne.nk) stop 'iowavef: Nunber of K-points inconsistent '
        if(nuototread.ne.nuotot) stop 'iowavef: Number of KS orbitals inconsistent ' 
        if(nspinread.ne.nspin) stop 'iowavef: Spin inconsistent '
      endif
      !   
      do ispin=1,nspin                                 
        do ik=1,nk
          if (Node.eq.0) read(unit1) dim2read
#ifdef MPI
          call MPI_Bcast(dim2read,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
          call m_allocate(wavef_rw(ik,ispin),nuotot,dim2read,m_storage)
        enddo
      enddo
      !
      do ispin=1,nspin
        do ik=1,nk
          do i=1,nuotot
            do j=1,wavef_rw(ik,ispin)%dim2
              if (Node==0) read(unit1) varaux
#ifdef MPI 
              call MPI_Bcast(varaux,1,MPI_double_complex,0,MPI_Comm_World,MPIerror)
#endif
              call m_set_element(wavef_rw(ik,ispin),i,j,varaux,'lap')
            end do
          end do
        enddo
      enddo
      if (Node.eq.0)   call io_close(unit1)
      !......................................................................   
    elseif(task.eq.'write'.or.task.eq.'WRITE') then     ! task write
      if (Node.eq.0) then
        call io_assign(unit1)
        open( unit1, file=fname, form='unformatted',status='unknown' )
        rewind(unit1)
        write(unit1) istpp, totime
        write(unit1) nuotot,nk,nspin
        !
        do ispin=1,nspin
          do ik=1,nk
            write(unit1) (wavef_rw(ik,ispin)%dim2)
          enddo
        enddo
        !
      end if
      !
      do ispin=1,nspin
        do ik=1,nk
          do i=1,nuotot
            do j=1,wavef_rw(ik,ispin)%dim2
              call m_get_element(wavef_rw(ik,ispin),i,j,varaux,'lap')
              if (Node==0) write(unit1) varaux
            enddo
          enddo
        enddo 
      enddo 
      !
      if (Node .eq. 0) then
        call io_close(unit1)
      endif
    endif            ! end task 
    !
  end subroutine iowavef
  !
end module wavefunctions       
