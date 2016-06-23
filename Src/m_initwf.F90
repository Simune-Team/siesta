module m_initwf
  use precision
  use wavefunctions, only:iowavef, wavef_ms
  use MatrixSwitch
  use sparse_matrices, only: numh, listhptr, listh, H, S, xijo
  use m_eo,            only: eo, qo
  !
  implicit none
  !
  private
  !
  public :: initwf
  !
  !type(matrix), allocatable, save :: wavef_ms(:,:)
CONTAINS
  !
  subroutine initwf(no,nspin,maxspn,maxuo,maxnh,         &
                    maxo, qtot,gamma,indxuo,             &
                    nk,kpoint,wk,nuotot,ef,        &
                    istpp,totime)
! *********************************************************************
! Subroutine to calculate the eigenvalues and eigenvectors,
! for given Hamiltonian and Overlap matrices (including
! spin polarization), providing the initial wavefunctions 
! for a time dependent electonic simulations.
! Written by D. Sanchez-Portal, November 2002-March 2003 after
! subroutine diagon by J.Soler, P.Ordejon, and J. D. Gale (1997-1999)
! Modified by Rafi Ullah, CIC nanoGUNE, October 2015 
! to make it  work with parallel TDDFT-siesta using MatrixSwitch
! **************************** INPUT **********************************
! integer no                  : Number of basis orbitals
! integer nspin               : Spin polarization (1 or 2)
! integer maxspn              : Second dimension of eo and qo
! integer maxnh               : Maximum number of orbitals interacting  
! integer maxnd               : Maximum number of nonzero elements of 
!                               each row of density matrix
! integer maxo                : First dimension of eo and qo
! integer numh(nuo)           : Number of nonzero elements of each row 
!                               of hamiltonian matrix
! integer listhptr(nuo)       : Pointer to each row (-1) of the
!                               hamiltonian matrix
! integer listh(maxlh)        : Nonzero hamiltonian-matrix element  
!                               column indexes for each matrix row
! integer numd(nuo)           : Number of nonzero elements of each row 
!                               of density matrix
! integer listdptr(nuo)       : Pointer to each row (-1) of the
!                               density matrix
! integer listd(maxnh)        : Nonzero density-matrix element column 
!                               indexes for each matrix row
! real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
! real*8  S(maxnh)            : Overlap in sparse form
! logical gamma               : Only gamma point?
! real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
!                               (not used if only gamma point)
! integer indxuo(no)          : Index of equivalent orbital in unit cell
!                               Unit cell orbitals must be the first in
!                               orbital lists, i.e. indxuo.le.nuo, with
!                               nuo the number of orbitals in unit cell
! integer nk                  : Number of k points
! real*8  kpoint(3,nk)        : k point vectors
! real*8  wk(nk)              : k point weights (must sum one)
! integer nuotot              : total number of orbitals in unit cell 
!                               over all processors
! *************************** OUTPUT **********************************
! real*8 eo(maxo,maxspn,nk)   : Eigenvalues
! real*8 qo(maxo,maxspn,nk)   : Occupations of eigenstates
! *************************** UNITS ***********************************
! xij and kpoint must be in reciprocal coordinates of each other.
! temp and H must be in the same energy units.
! eo, Enew and ef returned in the units of H.
! *************************** Parallel ********************************
! When running in parallel some of the dimensions are now the 
! maximum per node and the corresponding number passed in as
! an argument is the number of locally stored values. The
! variables for which this is the case are:
!
! maxuo/no
!
! *********************************************************************
!
!  Modules
!
      use parallel,     only : Node, Nodes, BlockSize
      use parallelsubs, only : GlobalToLocalOrb, GetNodeOrbs
      use fdf
      use densematrix,  only : Haux, Saux, psi
      use alloc
      use m_memory
      use m_fermid,      only : fermid
      use sys, only: die
#ifdef MPI
      use mpi_siesta,   only : mpi_bcast, mpi_comm_world, &
                               mpi_logical, mpi_double_precision
#endif
      !
      implicit none
      !
      !
      integer, intent (in) :: no, nspin, maxspn, maxuo, maxnh, maxo, nk,       &
                              nuotot, indxuo(no)
      integer, intent (inout) :: istpp
      !
      real(dp), intent (in)    :: kpoint(3,nk), wk(nk), qtot
      real(dp), intent (inout) :: totime
      logical, intent (in)     :: gamma
      ! 
      external           :: io_assign, io_close, readsp, paste
      !
      character          :: sname*30, fname*33, paste*33, m_storage*5
      !
      logical            :: fixspin, ioifound, degen
      logical, save      :: spiral
      logical, save      :: frstme = .true.
      !
      !
      integer            :: io, iuo, iu, nhs, npsi, nuo, nocc(2), ispin,ik,     &
                            i, j, sumqo,ikmax,iomax,ioi,iof,ispinmax

      real(dp)           :: qspiral(3), ef, temp,nelect, entrp,qomax,qtol
      complex(dp)        :: varaux,varaux2,varaux3,varaux4
#ifdef MPI
      integer            :: MPIerror
      logical, save      :: ParallelOverK
#endif
!     Dynamic arrays
      integer, dimension(:),     allocatable, save :: muo
      integer, dimension(:,:),   allocatable, save :: nocck
      logical, dimension(:,:,:), allocatable, save :: occup
!     First call initialisation
      if (frstme) then
#ifdef MPI
        if (Node.eq.0) then
          ParallelOverK = fdf_boolean( 'Diag.ParallelOverK', .false. )
        end if
        if(ParallelOverk) then
          stop "initwf: is not parallelized over k-points."
        end if
        call MPI_Bcast(ParallelOverK,1,MPI_logical,0,MPI_Comm_World, MPIerror)
#endif
!       Read spin-spiral wavevector (if defined)
        call readsp( qspiral, spiral )
        if (spiral.and.Node.eq.0) then
          if (gamma) write(6,*) &
            'diagon: WARNING: spin-spiral requires k sampling'
          if (nspin.ne.4) write(6,*) &
            'diagon: WARNING: spin-spiral requires nspin=4'
        end if
        frstme = .false.
      end if
!     Get Node number and calculate local orbital range
#ifdef MPI
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)
#else
      nuo = nuotot
#endif
!     Start time counter ................................................
      call timer( 'initwf', 1 )
!     Check internal dimensions ..........................................
      if (nspin.le.2 .and. gamma) then
        nhs  = nuotot * nuo
        npsi = nuotot * maxuo * nspin
      else if (nspin.le.2 .and. .not.gamma) then
        nhs  = 2 * nuotot * nuo
        npsi = 2 * nuotot * nuo
#ifdef MPI
        if (ParallelOverK) then
          nhs  = 2 * nuotot * nuotot
          npsi = 2 * nuotot * nuotot
        end if
#endif
      else if (nspin.eq.4) then 
        if(Node.eq.0 ) write(6,'(a,/,a)') &
          'initwf: Electron-Ion dynamics not yet', &
          'initwf: implemented for non-collinear spin'  
        stop
        nhs  = 2 * (2*nuotot) * (2*nuo)
        npsi = 2 * (2*nuotot) * (2*maxuo)
      else
        call die('diagon: ERROR: incorrect value of nspin')
      end if
!     Allocate local arrays
      call re_alloc(Haux,1,nhs,name='Haux',routine='initwf')
      call re_alloc(Saux,1,nhs,name='Saux',routine='initwf')
      call re_alloc(psi,1,npsi,name='psi',routine='initwf')
      allocate(muo(nuo),stat=mem_stat)
      call memory('A','I',nuo,'initwf',stat=mem_stat)
      allocate(nocck(nk,nspin),stat=mem_stat)
      call memory('A','I',nk*nspin,'initwf',stat=mem_stat)
      allocate(occup(nuotot,nspin,nk),stat=mem_stat)
      call memory('A','L',nuo*nk*nspin,'initwf',stat=mem_stat)
!     Check indxuo .......................................................
      do iuo = 1,nuo
        muo(iuo) = 0
      end do
      do io = 1,no
        iuo = indxuo(io)
        if (iuo.le.0 .or. iuo.gt.nuotot) then
          if (Node.eq.0) then
            write(6,*) 'initwf: ERROR: invalid index: io, indxuo =',io, indxuo(io)
            stop 'initwf: ERROR: invalid indxuo'
          else
            stop
          end if
        end if
        call GlobalToLocalOrb(indxuo(io),Node,Nodes,iuo)
        if (iuo.gt.0) then
          muo(iuo) = muo(iuo) + 1
        end if
      end do
      do iuo = 1,nuo
        if (muo(iuo) .ne. muo(1)) then
          if (Node.eq.0) then
            write(6,'(/,2a,3i6)') 'initwf: ERROR: inconsistent indxuo', &
             '. iuo, muo(iuo), muo(1) =', iuo, muo(iuo), muo(1)
            stop 'initwf: ERROR: inconsistent indxuo.'
          else
            stop
          end if
        end if
      end do
! ............................................................................!
!     Determine the number of occupied states this is not necessarily trivial !
!     if we have a metal. Only the occupied KS are saved and subsequently     !
!     evolved by integrating TDKS equations.                                  !
! ............................................................................!
      temp=1.0d-6
      call fermid( nspin, maxspn, nk, wk, maxo, nuotot, eo, &
                   temp, qtot, qo, ef, entrp )
      nocc(1) = 0
      nocc(2) = 0
      nelect=0.0d0
      degen= .false.
      !
      if (Node .eq. 0) then
      write(6,fmt="(a,tr3,a,tr3,a,tr3,a)")   ,"initwf:","ik", "occupancy","maximum occupancy"
      end if
      !
      do ik=1,nk
        do ispin=1,nspin
          nocck(ik,ispin)=0
          do io=1,nuotot
            occup(io,ispin,ik)=.false.
            if(dabs(qo(io,ispin,ik)-2.0d0*wk(ik)/nspin).le.    &
               1.0d-2*dabs(2.0d0*wk(ik)/nspin))  then
              nocc(ispin)=nocc(ispin)+1
              nocck(ik,ispin)=nocck(ik,ispin)+1
!             Accounting the number of electrons corresponding the states being marked
!             as occupied.
              nelect=nelect+dabs(2.0d0*wk(ik)/nspin)
              occup(io,ispin,ik)=.true.
            else
              if ( dabs( qo(io,ispin,ik)) .gt.1.0d-2*dabs(2.0d0*wk(ik)/nspin)) then
                write(6,"(tr2,I10,tr3,f8.6,tr4,f8.6)") ik, qo(io,ispin,ik), &
                  2.0d0*wk(ik)/nspin
                degen = .true.  
              end if
            end if
          end do
        end do
      end do
!-------------------------------------------------------------------------------!
! Systems with odd number of electrons in spin-unpolarized calctions            !
! may encounter situations with partial occupations. In such case we            !
! try to print the occupation map around Fermi level to have a clear picture.   !
! However the program will stop after printing this information in case of      !
! partial occupations.                                                          !
!-------------------------------------------------------------------------------!
      if(Node .eq. 0) then
        write(6,*) "initwf: No. of electrons corresponding occupied states =  ", nelect
        write(6,*) "initwf: (Total charge - charge in selected states) =  ",qtot-nelect
      end if
!     Tolrance for charge discrepacy 
      qtol = 1.0d-2
      if(nelect .gt. 0.d0 .and. (qtot-nelect) .gt. qtol) then
!       Find the partially filled orbitals
        ioifound=  .false.
        do ispin=1,nspin
          do io=1,nuo
            sumqo   =  0 ! integer
            do ik=1,nk
              if (occup(io,ispin,ik) ) sumqo=sumqo+1
           end do
           if (sumqo .lt. nk ) then ! Unoccupied orbitals at some k-points.
    	        if (Node .eq. 0)   write(6,*)                                          &
                   "initwf: k-points with occupied orbital ,total k-points = ",sumqo,nk
			        if(.not. ioifound .and. sumqo .gt. 0) then
			          ioi=io
			          ioifound=.true.
		            if(Node .eq. 0) write(6,*) "initwf: First partially filled orbital = ",ioi
			        end if
			        if (sumqo .eq. 0) then
			          iof=io
		            if(Node .eq. 0)  write(6,*) "initwf: First completely unoccupied orbital = ",iof
			          exit ! exit the loop on io becaseu occupations have been found
			        end if
		        end if ! sumqo .lt. nk 	
		      end do  ! io
          if(sumqo .eq. 0) exit 
		    end do  !ispin
        !..............
		    if(Node.eq.0) then  
		      if(nspin.eq.2) then
            write(6,*) 'initwf: number of occupied wave functions'
            write(6,*) 'initwf: spin up    ',nocc(1)
            write(6,*) 'initwf: spin down  ',nocc(2)
		      else
            write(6,*) 'initwf: number of occupied wave functions ',               &
               nocc(1)
		      end if
		    end if
      end if  ! qtol
!--------------------------------------------------------------------------------------!
!     Stop if the system has degenracy.
      if (degen) then
        STOP "initwf: System has degenracy. Change spin polarization or Fermi &
        temperature to avoid it"
      end if
      IF(Node.eq.0) write(6,*) 'Debug01: Occupations complete!!!!!!!!!!!!!!'
      !..............
#ifdef MPI
      call ms_scalapack_setup(Nodes,1,'c',BlockSize)
      m_storage='pddbc'
#else
      m_storage='sdden'
#endif
      allocate(wavef_ms(1:nk,1:nspin)) ! allocate (nk*npsin) matrices inside wavef_ms
      do i=1,nk !for every value of nk and nspin, allocate a matrix of size (nuotot x nocck(i,j))
        do j=1,nspin
          call m_allocate(wavef_ms(i,j),nuotot,nocck(i,j),m_storage)
        end do
      end do

      IF(Node.eq.0) write(6,*) 'Debug02: MatrixSwitch allocations!!!!!!!!!!!!!!'
!     Call apropriate routine .............................................
      if (nspin.le.2 .and. gamma) then
        call diaggiwf( nspin, nuo, maxuo, maxnh, maxo,                     &
                    Haux, Saux, psi, nuotot, occup)
      else if (nspin.le.2 .and. .not.gamma) then
          call diagkiwf( nspin, nuo, no, maxspn, maxuo, maxnh,                 &
                         maxo, indxuo, nk, kpoint, Haux, Saux,                 &
                         psi, nuotot, occup)
      else 
         stop                                                                  &
         'initwf: ERROR: non-collinear spin options for TDDFT not yet implemented'
      end if
      IF(Node.eq.0) write(6,*) 'Debug03: Gamma-diagonalization !!!!!!!!!!!!!!'
!     Write/save wavefunction in .TDWF file to use for TDDFT calculation.
      call  iowavef('write',wavef_ms,nuotot,nk,nspin,istpp,totime)
!     Free local arrays
      call memory('D','I',size(muo),'initwf',stat=mem_stat)
      deallocate(muo,stat=mem_stat)
      call memory('D','I',size(nocck),'initwf',stat=mem_stat)
      deallocate(nocck,stat=mem_stat)
      call memory('D','L',size(occup),'initwf',stat=mem_stat)
      deallocate(occup,stat=mem_stat)

      call de_alloc( Haux, 'Haux', 'initwf')
      call de_alloc( Saux, 'Saux', 'initwf')
      call de_alloc( psi,  'phi',  'initwf')

!     Stop time counter ...................................................
      call timer( 'initwf', 2 )
      !
      IF(Node.eq.0) write(6,*) 'Debug04: Free-local arrays !!!!!!!!!!!!!!'
  end subroutine initwf
  ! Gamma point: solve KS by diagonalisation and store the occupied wavefunctions in wavef
  subroutine diaggiwf(nspin,nuo,maxuo,maxnh, maxo,Haux,Saux,psi,           &
                      nuotot,occup)
#ifdef MPI
      use parallel, only : BlockSize,Node
      use m_diagon, only : ictxt
#endif
      !
      implicit none
      !
      integer, intent(in)         :: maxnh, maxuo, maxo, nuo, nspin, nuotot
      real(dp), intent(inout)  :: Haux(nuotot,nuo), Saux(nuotot,nuo), psi(nuotot,maxuo,nspin)
      logical, intent(inout)      :: occup(nuotot,nspin,1)
      ! Internal variables
      integer                     :: ie, io, ispin, j, jo, ind, ierror, ioc, indwf
      real(dp)                    :: element
#ifdef MPI
      integer                     :: desch(9)
#else
      integer                     :: Node=0
#endif
      !
#ifdef MPI
        call descinit(desch,nuotot,nuotot,BlockSize,BlockSize,0,0,ictxt,nuotot,ierror)
#endif
      !
      indwf=0
      do ispin=1,nspin
        do ie=1,10
          Saux(1:nuotot,1:nuo) = 0.0d0
          Haux(1:nuotot,1:nuo) = 0.0d0
          do io=1,nuo
            do j=1,numh(io)
              ind=listhptr(io)+j
              jo=listh(ind)
              Saux(jo,io)=Saux(jo,io)+S(ind)
              Haux(jo,io)=Haux(jo,io)+H(ind,ispin)
            end do
          end do
          IF(Node.eq.0) write(6,*) 'Gamam-diag: sparse2dense !!!!!!!!!!!!!!'
          IF(Node.eq.0) print*, shape(Haux)
          IF(Node.eq.0) print*, shape(Saux)
          IF(Node.eq.0) print*, nuotot
          IF(Node.eq.0) print*, nuo
          IF(Node.eq.0) print*, shape(eo)
          IF(Node.eq.0) print*, maxo, ispin
          IF(Node.eq.0) print*, shape(psi)
          IF(Node.eq.0) print*, ispin
          call rdiag(Haux,Saux,nuotot,nuo,nuotot,eo,psi(1,1,ispin),nuotot,1,ierror)
          IF(Node.eq.0) write(6,*) 'Gamam-diag: after-diag !!!!!!!!!!!!!!'
          if (ierror .eq. 0) then
            exit
          else if ((ierror .ne. -1) .or. (ie .eq. 10)) then
          call die('Terminating due to failed diagonalisation')
          end if
        end do     ! ie
!.............................       
        ioc=0
        do ie=1,nuotot
          if (occup(ie,ispin,1)) then 
            ioc=ioc+1
            indwf=indwf+1
            do j=1,nuotot
#ifdef MPI
              call pdelget('a',' ',element,psi(:,:,ispin),j,ie,desch)
#else
              element=psi(j,ie,ispin)
#endif
              call m_set_element(wavef_ms(1,ispin),j,ioc,element,'lap')
            end do          ! j=1,nuotot
          end if            ! occup
        end do              ! ie=1,nuotot
      end do                ! ispin
      !
          IF(Node.eq.0) write(6,*) 'Gamam-diag: wavef2M_S !!!!!!!!!!!!!!'
  end subroutine diaggiwf
  ! k points: solve KS by diagonalisation and store the occupied wavefunctions in wavef
  subroutine diagkiwf(nspin,nuo,no,maxspn,maxuo,maxnh, maxo, indxuo,nk,        &
                      kpoint, Haux,Saux,psi,nuotot,occup)
#ifdef MPI
      use parallel, only : BlockSize
      use m_diagon, only : ictxt
#endif
      !
      implicit none
      !
      integer, intent(in)      :: maxnh, maxuo, maxo, no, nspin, nuo,          &
                                  nuotot, nk, maxspn, indxuo(no)
      real(dp), intent(in)     :: kpoint(3,nk)
      real(dp), intent(inout)  :: Haux(2,nuotot,nuo), Saux(2,nuotot,nuo), psi(2,nuotot,nuo)
      logical, intent(inout)   :: occup(nuotot,nspin,nk)
      ! Internal variables
      integer                  :: ispin, ie, ierror, ik, ind, iuo, j, jo, juo, indwf, ioc
      real(dp)                 :: ckxij, kxij, skxij, varaux(2)
      complex(dp)              :: element
#ifdef MPI
      integer                  :: desch(9)
#endif
      !
#ifdef MPI
      call descinit(desch,nuotot,nuotot,BlockSize,BlockSize,0,0,ictxt,nuotot,ierror)
#endif
      !
      indwf=0
      do ik=1,nk
        do ispin=1,nspin
          Saux(1:2,1:nuotot,1:nuo)=0.0d0
          Haux(1:2,1:nuotot,1:nuo)=0.0d0
          do iuo=1,nuo
            do j=1,numh(iuo)
              ind=listhptr(iuo)+j
              jo=listh(ind)
              juo=indxuo(jo)
              kxij=kpoint(1,ik)*xijo(1,ind)+                      &
              kpoint(2,ik)*xijo(2,ind)+                           &
              kpoint(3,ik)*xijo(3,ind)
              ckxij=cos(kxij)
              skxij=sin(kxij)
!              Note: sign of complex part changed to match change in order of iuo/juo
              Saux(1,juo,iuo)=Saux(1,juo,iuo)+S(ind)*ckxij
              Saux(2,juo,iuo)=Saux(2,juo,iuo)-S(ind)*skxij
              Haux(1,juo,iuo)=Haux(1,juo,iuo)+H(ind,ispin)*ckxij
              Haux(2,juo,iuo)=Haux(2,juo,iuo)-H(ind,ispin)*skxij
            end do
          end do
          !
          call cdiag(Haux,Saux,nuotot,nuo,nuotot,eo(1,ispin,ik),psi,nuotot,1,ierror)
          if (ierror .ne. 0) then
          call die('Terminating due to failed diagonalisation')
          end if
          !.....................
          ioc=0
          do ie=1,nuotot
            if (occup(ie,ispin,ik)) then 
              ioc=ioc+1
              indwf=indwf+1
              do j=1,nuotot
#ifdef MPI
                do iuo=1,2
                call pdelget('a',' ',varaux(iuo),psi(iuo,:,:),j,ie,desch)
                enddo
                element=cmplx(varaux(1),varaux(2),dp)
#else
                element=cmplx(psi(1,j,ie),psi(2,j,ie),dp)
#endif
                call m_set_element(wavef_ms(ik,ispin),j,ioc,element,'lap')
              end do ! j
            end if
          end do ! ie
        end do  ! do ispin
      end do    ! do ikmax
  end subroutine diagkiwf
end module m_initwf
