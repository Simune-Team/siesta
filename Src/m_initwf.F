      module m_initwf 
      implicit none
      private

      public :: initwf
      

      CONTAINS
      
      
      
      
      subroutine initwf(no, nspin, maxspn, maxuo, maxnh, 
     .                  maxo, numh, listhptr, listh, 
     .                  H, S, qtot,
     .                  gamma, xij, indxuo, nk, kpoint, wk,
     .                  eo, qo, nuotot,ef, istpp,totime)
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors,
C for given Hamiltonian and Overlap matrices (including
C spin polarization), providing the initial wavefunctions 
C for a time dependent electonic simulations.
C Written by D. Sanchez-Portal, November 2002-March 2003 after
C subroutine diagon by J.Soler, P.Ordejon, and J. D. Gale (1997-1999)
C **************************** INPUT **********************************
C integer no                  : Number of basis orbitals
C integer nspin               : Spin polarization (1 or 2)
C integer maxspn              : Second dimension of eo and qo
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : Maximum number of nonzero elements of 
C                               each row of density matrix
C integer maxo                : First dimension of eo and qo
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C integer listh(maxlh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer numd(nuo)           : Number of nonzero elements of each row 
C                               of density matrix
C integer listdptr(nuo)       : Pointer to each row (-1) of the
C                               density matrix
C integer listd(maxnh)        : Nonzero density-matrix element column 
C                               indexes for each matrix row
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C logical gamma               : Only gamma point?
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C integer nk                  : Number of k points
C real*8  kpoint(3,nk)        : k point vectors
C real*8  wk(nk)              : k point weights (must sum one)
C integer nuotot              : total number of orbitals in unit cell 
C                               over all processors
C *************************** OUTPUT **********************************
C real*8 eo(maxo,maxspn,nk)   : Eigenvalues
C real*8 qo(maxo,maxspn,nk)   : Occupations of eigenstates
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C temp and H must be in the same energy units.
C eo, Enew and ef returned in the units of H.
C *************************** Parallel ********************************
C When running in parallel some of the dimensions are now the 
C maximum per node and the corresponding number passed in as
C an argument is the number of locally stored values. The
C variables for which this is the case are:
C
C maxuo/no
C
C *********************************************************************
C
C  Modules
C
      use precision
      use parallel,     only : Node, Nodes
      use parallelsubs, only : GlobalToLocalOrb, GetNodeOrbs
      use fdf
      use densematrix,  only : Haux, Saux, psi
      use alloc
      use m_memory
      use m_fermid,      only : fermid
      use sys, only: die
#ifdef MPI
      use mpi_siesta,   only : mpi_bcast, mpi_comm_world,
     $                         mpi_logical, mpi_double_precision
#endif
      use wavefunctions
c maz. added ClProj.
      use elec_dyn_options, only: usesaveTDWF,ClProj

      implicit none

      integer
     .  maxnd, maxnh, maxspn, maxuo, maxo, nk, no, nspin, nuotot,
     . sumqo,ikmax,iomax,ioi,iof,ispinmax, istpp

      integer 
     .  indxuo(no), listh(maxnh), numh(*), 
     .  listhptr(*)

      double precision
     .  eo(maxo,maxspn,nk), H(maxnh,nspin), kpoint(3,nk), 
     .  qo(maxo,maxspn,nk), S(maxnh), wk(nk), 
     .  xij(3,maxnh), qspiral(3), qtot,  ef, temp, nelect, entrp,
     . qomax,qtol, totime
       
      character 
     .   sname*30, fname*33, paste*33

      logical
     .  fixspin, gamma, ioifound

      logical, save :: spiral
      logical, save :: frstme = .true.

      external
     .  io_assign, io_close,
     .  readsp, paste
     
#ifdef MPI
      integer
     .  MPIerror
      external
     .  diagkp
#endif


      integer
     .  io, iuo, iu, naux, nhs, npsi, nuo, nocc(2), ispin,
     .  ik

C     Dynamic arrays
      double precision, dimension(:), allocatable, save :: aux
      integer         , dimension(:), allocatable, save :: muo
      integer         , dimension(:,:), allocatable, save :: nocck
      logical         , dimension(:,:,:), allocatable, save :: occup
#ifdef MPI
      logical, save ::
     .  ParallelOverK
#endif

C  ....................

C First call initialisation
      if (frstme) then
#ifdef MPI
        if (Node.eq.0) then
          ParallelOverK = fdf_boolean( 'Diag.ParallelOverK', .false. )
        endif
        call MPI_Bcast(ParallelOverK,1,MPI_logical,0,MPI_Comm_World,
     .    MPIerror)
#endif

C Read spin-spiral wavevector (if defined)
        call readsp( qspiral, spiral )
        if (spiral.and.Node.eq.0) then
          if (gamma) write(6,*)
     .      'diagon: WARNING: spin-spiral requires k sampling'
          if (nspin.ne.4) write(6,*)
     .      'diagon: WARNING: spin-spiral requires nspin=4'
        endif

        frstme = .false.
      endif

C Get Node number and calculate local orbital range
#ifdef MPI
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)
#else
      nuo = nuotot
#endif

C Start time counter ................................................
      call timer( 'initwf', 1 )
         
C Read wavefunctions from a previous electron-dynamics simulation 
      
!Debug - Rafi      
	 	
       !  call  iowavef('read',wavef,nuotot,nk,nspin)
!End - Rafi

C If we do not read information from a previous simulation
      
C Check internal dimensions ..........................................
      if (nspin.le.2 .and. gamma) then
        nhs  = nuotot * nuo
        npsi = nuotot * maxuo * nspin
        naux = nuotot * 5
      elseif (nspin.le.2 .and. .not.gamma) then
        nhs  = 2 * nuotot * nuo
        npsi = 2 * nuotot * nuo
        naux = 2 * nuotot * 5            !!!! AG: Put factor 2 !!
#ifdef MPI
        if (ParallelOverK) then
          nhs  = 2 * nuotot * nuotot
          npsi = 2 * nuotot * nuotot
        endif
#endif
      elseif (nspin.eq.4) then 
        if(Node.eq.0 ) write(6,'(a,/,a)')
     .    'initwf: Electron-Ion dynamics not yet', 
     .    'initwf: implemented for non-collinear spin'  
        stop
        nhs  = 2 * (2*nuotot) * (2*nuo)
        npsi = 2 * (2*nuotot) * (2*maxuo)
        naux = (2*nuotot) * 5
      else
        call die('diagon: ERROR: incorrect value of nspin')
      endif

C Allocate local arrays
      call re_alloc(Haux,1,nhs,name='Haux',routine='initwf')
      call re_alloc(Saux,1,nhs,name='Saux',routine='initwf')
      call re_alloc(psi,1,npsi,name='psi',routine='initwf')
      allocate(aux(naux),stat=mem_stat)
      call memory('A','D',naux,'initwf',stat=mem_stat)
      allocate(muo(nuo),stat=mem_stat)
      call memory('A','I',nuo,'initwf',stat=mem_stat)
      allocate(nocck(nk,nspin),stat=mem_stat)
      call memory('A','I',nk*nspin,'initwf',stat=mem_stat)
      allocate(occup(nuo,nspin,nk),stat=mem_stat)
      call memory('A','L',nuo*nk*nspin,'initwf',stat=mem_stat)
C Check indxuo .......................................................
      do iuo = 1,nuo
        muo(iuo) = 0
      enddo
      do io = 1,no
        iuo = indxuo(io)
        if (iuo.le.0 .or. iuo.gt.nuotot) then
          if (Node.eq.0) then
            write(6,*) 'initwf: ERROR: invalid index: io, indxuo =',
     .        io, indxuo(io)
            stop 'initwf: ERROR: invalid indxuo'
          else
            stop
          endif
        endif
        call GlobalToLocalOrb(indxuo(io),Node,Nodes,iuo)
        if (iuo.gt.0) then
          muo(iuo) = muo(iuo) + 1
        endif
      enddo
      do iuo = 1,nuo
        if (muo(iuo) .ne. muo(1)) then
          if (Node.eq.0) then
            write(6,'(/,2a,3i6)') 'initwf: ERROR: inconsistent indxuo',
     .       '. iuo, muo(iuo), muo(1) =', iuo, muo(iuo), muo(1)
            stop 'initwf: ERROR: inconsistent indxuo.'
          else
            stop
          endif
        endif
      enddo
C ......................

C Determine the number of occupied states 
C this is not necessarily trivial if we have a metal
c (if ClProj, only Ef is needed).
      temp=1.0d-6
      call fermid( nspin, maxspn, nk, wk, maxo, nuotot, eo,
     .             temp, qtot, qo, ef, entrp )
      nocc(1) = 0
      nocc(2) = 0
      nelect=0.0d0

C	********
C	MAZ . if ClProj, then make nocc=nuotot for every k point. 
	if (ClProj) then

	do ik=1,nk
	do ispin=1,nspin
         nocck(ik,ispin)=0
         do io=1,nuo
                   nocc(ispin)=nocc(ispin)+1
                   nocck(ik,ispin)=nocck(ik,ispin)+1

		if (dabs(qo(io,ispin,ik)
     .              - 2.0d0*wk(ik)/nspin).le.
     .        1.0d-2*dabs(2.0d0*wk(ik)/nspin))  then
                   nelect=nelect+dabs(2.0d0*wk(ik)/nspin)
		end if
                   occup(io,ispin,ik)=.true.
	  enddo
	 enddo
	enddo

	if(Node.eq.0) then  
	if(nspin.eq.2) then
         write(6,*) 'initwf: number of total wave functions'
         write(6,*) 'initwf: spin up    ',nocc(1)
         write(6,*) 'initwf: spin down  ',nocc(2)
	else
         write(6,*) 'initwf: number of total wave functions ',
     .       nocc(1)
	endif
	write(6,*) 'initwf: number of electrons  ',nelect
	endif


	else 
C	********
      do ik=1,nk
       do ispin=1,nspin
         nocck(ik,ispin)=0
         do io=1,nuo
              occup(io,ispin,ik)=.false.
              if (dabs(qo(io,ispin,ik)
     .              - 2.0d0*wk(ik)/nspin).le.
     .        1.0d-2*dabs(2.0d0*wk(ik)/nspin))  then
                   nocc(ispin)=nocc(ispin)+1
                   nocck(ik,ispin)=nocck(ik,ispin)+1
                   nelect=nelect+dabs(2.0d0*wk(ik)/nspin)
                   occup(io,ispin,ik)=.true.
              else
               if ( dabs( qo(io,ispin,ik)) .gt.
     .           1.0d-2*dabs(2.0d0*wk(ik)/nspin)) then
           write(6,*) ik, qo(io,ispin,ik),
     .         qo(io,ispin,ik)/(2.0d0*wk(ik)/nspin)
c                  stop 'initwf: ERROR system with degeneracy'
               endif
              endif
	  enddo
	 enddo
	enddo
 


c--------------c--------------c--------------c--------------c--------------
! maz: Nov 9, 2011. 
! to select only lowest lying energy states that can be filled with electrons.
! the same as suing a sharp step distribution function, or FD function at T_e=0K. 
! only for spin unpolarised calcualtions!

c qtol = tolerance for q!	
	qtol=1.0d-2
	write(6,*) "initwf: nelect =  ", nelect
	if(nelect .gt. 0.d0 .and. (qtot-nelect) .gt. qtol) then

	write(6,*) "initwf: (qtot-nelect) =  ",qtot-nelect

! find the partial filled range of io's.

	ioifound=.false.

	do ispin=1,nspin
	do io=1,nuo


		sumqo=0 ! integer

		do ik=1,nk
		if (occup(io,ispin,ik) ) sumqo=sumqo+1
		end do



		if (sumqo .lt. nk ) then ! there are empty states at some k points.

	write(6,*) "initwf: sumqo,nk = ",sumqo,nk

		if(.not. ioifound) then
		ioi=io
		ioifound=.true.
	write(6,*) "initwf: maz: ioifound! ioi= ",ioi
		end if

		if (sumqo .eq. 0) then
		iof=io
	write(6,*) "initwf: maz: ioffound! iof= ",iof
		go to 21 ! exit the loop on io becaseu we have found ioi and iof.
		end if

		end if ! sumqo .lt. nk 	


	end do  ! io
	end do  !ispin

 21 	continue

! 	now find the lowest energy empty state, i.e., with largest qo.
	qomax=0.0d0

	do ispin=1,nspin
	do io=max(1,ioi-20),min(iof+20,nuo)
	do ik=1,nk
	
		if(.not. occup(io,ispin,ik)) then

		if (qo(io,ispin,ik) .gt. qomax) then 
		qomax=qo(io,ispin,ik)
		ispinmax=ispin
		iomax=io
		ikmax=ik

		end if

		end if



	end do ! ik
	end do
	end do

	write(6,*) "initwf: qomax,iomax,ikmax= ",qomax,iomax,ikmax


!	set this state as occupied!
!	add fraction of electron in it to total no. of electrons filling the states.

                   nocc(ispinmax)=nocc(ispinmax)+1
                   nocck(ikmax,ispinmax)=nocck(ikmax,ispinmax)+1
                   nelect=nelect+dabs(2.0d0*wk(ikmax)/nspin)
                   occup(iomax,ispinmax,ikmax)=.true.
	write(6,*)"initwf: occup=false -->> true ; io,ik=",iomax,ikmax
	write(6,*)"initwf: occup(iomax,ispinmax,ikmax)=",
     .   occup(iomax,ispinmax,ikmax)
	write(6,*)"initwf: nocck(ikmax,ispin)=",nocck(ikmax,ispinmax)
	write(6,*)"initwf: nocc(1)=",nocc(1)
	write(6,*)"initwf: nelec,qtot=",nelect,qtot


!	check whether the no. of electrons filling the states is equal to total electrons?


c	if(nspin.eq.2) then
c
c	  if ( int(qtot) .gt. (nocc(1)+nocc(2)) ) then ! still empty states.
c
c		go to 21 ! repeat above cycle to see if 
!c		 there is any other state that can be set to occupied
c		else ! no. of electrons = no. of electrons in occupied states. ?
c		go to 22
c	  end if
c
c 
c
c	else ! nspin=1
c	   if (qtot .gt. 2*nocc(1)/nk ) then
c		go to 21 ! repeat above cycle to see if 
!c		 there is any other state that can be set to occupied
c		else ! no. of electrons = no. of electrons in occupied states. ?
c		go to 22
c	   end if
c
c	endif




		if ((qtot-nelect) .gt. qtol) then ! still empty states.
		go to 21 ! repeat above cycle to see if 
!c		 there is any other state that can be set to occupied
		else ! no. of electrons = no. of electrons in occupied states. ?
		go to 22
		end if

	end if ! if nelect .lt. qtot

 22	write(6,*)"electrons in selected orb (tddft),qtot=",nelect,qtot

c--------------c--------------c--------------c--------------c--------------   




	if(Node.eq.0) then  
	if(nspin.eq.2) then
         write(6,*) 'initwf: number of occupied wave functions'
         write(6,*) 'initwf: spin up    ',nocc(1)
         write(6,*) 'initwf: spin down  ',nocc(2)
	else
         write(6,*) 'initwf: number of occupied wave functions ',
     .       nocc(1)
	endif
	write(6,*) 'initwf: number of electrons  ',nelect
	endif 

	end if
C	***********
C	MAZ 16 March 2011.
C	remember:if (ClProj),  nocc=nuotot; its total number of KS-orbital then.


c maz. 26 april, 2011.
c to get the file "qo-file" that would contain percentage occupancy of states:
c in search of missing norm problem.

	call io_assign(iu)
	open(iu,file="qo-file",status="unknown")

	write(iu,*) "((100*qo(io,ispin,ik)/wk,"//
     .             "occup(io,ispin,ik)),ik=1,nk)"
	do ispin=1,nspin
	do io=1,nuo
	write(iu,*) (50.0d0*qo(io,ispin,ik)*nspin/wk(ik),
     .       ik=1,nk)
	write(iu,*) ( occup(io,ispin,ik)  ,ik=1,nk)

	end do
	end do 

	call io_close(iu)



C Allocate space for the wavefunctions storage.........................
      call alloc_wavef(wavef,nuotot,nocck,nk,nspin)
      wavef%nocck(1:nk,1:nspin)=nocck(1:nk,1:nspin)
      wavef%iwf(:,:,:)=0
C......................................................................

C Call apropriate routine .............................................
      if (nspin.le.2 .and. gamma) then
        call diaggiwf( nspin, nuo, maxuo, maxnh, maxo,
     .              numh, listhptr, listh, H, S, eo, 
     .              Haux, Saux, psi, aux, 
     .              nuotot, occup)
      elseif (nspin.le.2 .and. .not.gamma) then
          call diagkiwf( nspin, nuo, no, maxspn, maxuo, maxnh, 
     .              maxo, numh, listhptr, listh,
     .              H, S, xij, indxuo, nk, kpoint,
     .              eo, Haux, Saux, psi, aux, 
     .              nuotot, occup)
      else 
         stop 
     .   'initwf: ERROR: non-collinear spin options not yet implemented'
      endif
C ....................

C Storing a copy of the initial wavefunctions for debugging purpose...
c     call alloc_wavef(init,nuotot,nocck,nk,nspin)
c     init%nocck(:,:)=wavef%nocck(:,:)
c     init%iwf(:,:,:)=wavef%iwf(:,:,:)
c     init%phi(:,:,:)=wavef%phi(:,:,:)
C ....................................................................


C Free local arrays
      call memory('D','I',size(muo),'initwf',stat=mem_stat)
      deallocate(muo,stat=mem_stat)
      call memory('D','D',naux,'initwf',stat=mem_stat)
      deallocate(aux,stat=mem_stat)
      call memory('D','I',size(nocck),'initwf',stat=mem_stat)
      deallocate(nocck,stat=mem_stat)
      call memory('D','L',size(occup),'initwf',stat=mem_stat)
      deallocate(occup,stat=mem_stat)

C If block on reading versus calculate wavefunctions   
      

      call  iowavef('write',wavef,nuotot,nk,nspin,istpp,totime)
C Stop time counter ...................................................
      call timer( 'initwf', 2 )
C .......................

      end subroutine initwf



      subroutine diaggiwf( nspin, nuo, maxuo, maxnh,
     .                  maxo, numh, listhptr, listh, 
     .                  H, S, eo, Haux, Saux, psi, aux,
     .                  nuotot, occup)
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization). Gamma-point version.
C Writen by J.Soler, August 1998.
C **************************** INPUT **********************************
C integer nspin               : Number of spin components (1 or 2)
C integer nuo                 : Number of basis orbitals local to node
C integer maxuo               : Last dimension of xij
C                               Must be at least max(indxuo)
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxo                : First dimension of eo 
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C *************************** OUTPUT **********************************
C real*8 eo(maxo,nspn)        : Eigenvalues
C real*8 psi(nuotot,maxuo,nspin) : Eigenvectors
C *************************** AUXILIARY *********************************
C real*8 Haux(nuotot,nuo)     : Auxiliary space for the hamiltonian matrix
C real*8 Saux(nuotot,nuo)     : Auxiliary space for the overlap matrix
C real*8 aux(nuotot*5)        : Extra auxiliary space
C *************************** UNITS ***********************************
C eo returned in the units of H.
C *************************** PARALLEL ********************************
C The auxiliary arrays are now no longer symmetry and so the order
C of referencing has been changed in several places to reflect this.
C *********************************************************************
C
C  Modules
C
      use precision
      use sys
      use parallel
      use wavefunctions
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

#ifdef MPI
      integer 
     .  MPIerror
#endif

      integer
     .  maxnh, maxuo, maxo, nuo, nspin, nuotot, 
     .  indwf

      integer 
     .  listh(maxnh), numh(nuo), listhptr(nuo),
     .  numd(nuo), nocc(nspin)

      double precision
     .  eo(maxo,nspin), H(maxnh,nspin),
     .  S(maxnh)
     
      double precision
     .  Haux(nuotot,nuo), Saux(nuotot,nuo), psi(nuotot,maxuo,nspin), 
     .  aux(5*nuotot)
      
      logical 
     .  occup(nuo,nspin)
      external
     .  rdiag

C  Internal variables .............................................
      integer           ie, io, iio, ispin, j, jo, BNode, iie, ind,
     .                  BTest, ierror, nd, ioc
      double precision  ee, pipj, qe, t
C  ....................

C Solve eigenvalue problem .........................................
      indwf=0
      do ispin = 1,nspin
   10   do io = 1,nuo
          do jo = 1,nuotot
            Saux(jo,io) = 0.0d0
            Haux(jo,io) = 0.0d0
          enddo
        enddo
        do io = 1,nuo
          do j = 1,numh(io)
            ind = listhptr(io) + j
            jo = listh(ind)
            Saux(jo,io) = Saux(jo,io) + S(ind)
            Haux(jo,io) = Haux(jo,io) + H(ind,ispin)
          enddo
        enddo
        call rdiag( Haux, Saux, nuotot, nuo, nuotot, 
     .              eo(1,ispin), psi(1,1,ispin), 
     .              nuotot, 1, ierror )
        if (ierror.gt.0) then
          call die('Terminating due to failed diagonalisation')
        elseif (ierror.eq.-1) then
          goto 10
        endif

        BNode = 0
        iie = 0
        ioc=0
        do ie = 1,nuotot

          if (Node.eq.BNode) iie = iie + 1

            if(occup(ie,ispin)) then 
             if(Node.eq.BNode) then 
              do j = 1,nuotot
                aux(j) = psi(j,iie,ispin)
              enddo
             endif
#ifdef MPI
              call MPI_Bcast(aux,nuotot,mpi_double_precision,BNode,
     .           MPI_Comm_World,MPIerror)
#endif   
          
          indwf=indwf+1
          ioc=ioc+1
          wavef%iwf(ioc,1,ispin)=indwf
          do j=1,nuotot
             wavef%phi(1,j,indwf)=aux(j)
             wavef%phi(2,j,indwf)=0.0d0
          enddo
          BTest = ie/BlockSize
          if (BTest*BlockSize.eq.ie) then
            BNode = BNode + 1
            if (BNode .gt. Nodes-1) BNode = 0
          endif
         endif
       enddo
       if(ioc.ne.wavef%nocck(1,ispin)) stop 'occup. inconsistency'
      enddo

      end subroutine diaggiwf



      subroutine diagkiwf( nspin, nuo, no, maxspn, maxuo, maxnh, 
     .                  maxo, numh, listhptr, listh, H, S,
     .                  xij, indxuo, nk, kpoint, eo,
     .                  Haux, Saux, psi, aux, nuotot, occup)
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization). K-sampling version.
C Writen by J.Soler, August 1998.
C **************************** INPUT **********************************
C integer nspin               : Number of spin components (1 or 2)
C integer nuo                 : Number of basis orbitals in unit cell
C                               local to this processor
C integer no                  : Number of basis orbitals in supercell
C integer maxuo               : Last dimension of xij
C                               Must be at least max(indxuo)
C integer maxo                : First dimension of eo and qo
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : First dimension of listd and DM
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C real*8  kpoint(3,nk)           : k point vectors
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C *************************** OUTPUT **********************************
C real*8 eo(maxo)             : Eigenvalues
C real*8 psi(2,nuotot,nuo)    : Eigenvectors
C *************************** AUXILIARY *******************************
C real*8 Haux(2,nuotot,nuo) : Auxiliary space for the hamiltonian matrix
C real*8 Saux(2,nuotot,nuo) : Auxiliary space for the overlap matrix
C real*8 aux(2,nuotot*5)    : Extra auxiliary space
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C eo returned in the units of H.
C *************************** PARALLEL ********************************
C The auxiliary arrays are now no longer symmetry and so the order
C of referencing has been changed in several places to reflect this.
C *********************************************************************
C
C  Modules
C
      use precision
      use sys
      use parallel
      use wavefunctions
#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

#ifdef MPI
      integer 
     .  MPIerror
#endif

      integer           maxnh, maxuo, maxo, no, ispin,
     .                  nspin, nuo, nuotot, nk, maxspn
      integer           indxuo(no), listh(maxnh), numh(nuo),
     .                  listhptr(nuo)
      double precision  
     .                  eo(maxo,maxspn,nk), H(maxnh,nspin),
     .                  kpoint(3,nk), S(maxnh), xij(3,maxnh)
      double precision  
     .                  Haux(2,nuotot,nuo), Saux(2,nuotot,nuo),
     .                  psi(2,nuotot,nuo), aux(2,nuotot*5)
      logical           occup(nuotot,nspin,nk)
      external          cdiag

C  Internal variables .............................................
      integer
     .  BNode, BTest, ie, ierror, iie, ik, ind, io, iio,
     .  iuo, j, jo, juo, nd, indwf, ioc
      double precision
     .  ckxij, ee, kxij, pipj1, pipj2, qe, skxij, t
C  ....................

C Find eigenvalues ................................................
     
      indwf=0
      do ik = 1,nk
        do ispin = 1,nspin
          Saux(1:2,1:nuotot,1:nuo) = 0.0d0
          Haux(1:2,1:nuotot,1:nuo) = 0.0d0
          do iuo = 1,nuo
            do j = 1,numh(iuo)
              ind = listhptr(iuo) + j
              jo = listh(ind)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xij(1,ind) +
     .               kpoint(2,ik) * xij(2,ind) +
     .               kpoint(3,ik) * xij(3,ind)
              ckxij = cos(kxij)
              skxij = sin(kxij)
C Note : sign of complex part changed to match change in order of iuo/juo
              Saux(1,juo,iuo)=Saux(1,juo,iuo)+S(ind)*ckxij
              Saux(2,juo,iuo)=Saux(2,juo,iuo)-S(ind)*skxij
              Haux(1,juo,iuo)=Haux(1,juo,iuo)+H(ind,ispin)*ckxij
              Haux(2,juo,iuo)=Haux(2,juo,iuo)-H(ind,ispin)*skxij
            enddo
          enddo
          call cdiag( Haux, Saux, nuotot, nuo, nuotot, 
     .        eo(1,ispin,ik), psi, 
     .        nuotot, 1, ierror )
          if (ierror.ne.0) then
            call die('Terminating due to failed diagonalisation')
          endif
C ....................
C Global operation to store wavefunctions.........................
C this is actually not very useful at the moment since TDDFT
C option only work in serial.
          BNode = 0
          iie = 0
          ioc=0
c         do ie = 1,wavef%nocck(ik,ispin) 
          do ie = 1, nuotot
             if (Node.eq.BNode) iie=iie+1
             if(occup(ie,ispin,ik)) then 
              if(Node.eq.BNode) then 
               do j = 1,nuotot
                aux(1,j) = psi(1,j,iie)
                aux(2,j) = psi(2,j,iie)
               enddo
              endif
#ifdef MPI
            call MPI_Bcast(aux(1,1),2*nuotot,mpi_double_precision,
     .        BNode,MPI_Comm_World,MPIerror)
#endif
           ioc=ioc+1
           indwf=indwf+1
           wavef%iwf(ioc,ik,ispin)=indwf
           do j=1,nuotot
             wavef%phi(1,j,indwf)=aux(1,j)
             wavef%phi(2,j,indwf)=aux(2,j)
           enddo
           BTest = ie/BlockSize
             if (BTest*BlockSize.eq.ie) then
              BNode = BNode + 1
              if (BNode .gt. Nodes-1) BNode = 0
             endif
            endif
          enddo
          if(ioc.ne.wavef%nocck(ik,ispin)) stop 'occup. inconsistency'
         enddo
        enddo

      end  subroutine diagkiwf

      END module m_initwf


