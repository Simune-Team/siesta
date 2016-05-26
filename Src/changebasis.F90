      subroutine chgbasis(no, nspin, maxspn, maxuo, maxnh, maxnd,&
                          maxo, numh, listhptr, listh, numd,&
                          listdptr, listd, S,&
                          gamma, xij, indxuo, nk, kpoint, wk,&
                          Dnew, nuotot,istpmove,psi,Saux)

!******************************************************************************
!Modified  by D. Sanchez-Portal, Feb 2009
!Modified by Adiran Garaizar, June 2015
!Modified by Rafi Ullah, October 2015
!**************************** INPUT ********************************************
!integer no                  : Number of basis orbitals the supercell
!integer nspin               : Spin polarization (1 or 2)
!integer maxspn              : Maximum number of spin orentations 
!integer maxuo               : Maximum number of orbitals stored in a 
!                              given Node
!integer maxnh               : Maximum number of orbitals interacting
!integer maxnd               : Maximum number of nonzero elements of
!                              each row of density matrix
!integer maxo                : Maximum number of orbitals in the unit cell
!integer numh(nuo)           : Number of nonzero elements of each row
!                              of hamiltonian matrix
!integer listhptr(nuo)       : Pointer to each row (-1) of the
!                              hamiltonian matrix
!integer listh(maxlh)        : Nonzero hamiltonian-matrix element
!                              column indexes for each matrix row
!integer numd(nuo)           : Number of nonzero elements of each row
!                              of density matrix
!integer listdptr(nuo)       : Pointer to each row (-1) of the
!                              density matrix
!integer listd(maxnh)        : Nonzero density-matrix element column
!                              indexes for each matrix row
!real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
!real*8  S(maxnh)            : Overlap in sparse form
!logical gamma               : Only gamma point?
!real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
!                              (not used if only gamma point)
!integer indxuo(no)          : Index of equivalent orbital in unit cell
!                              Unit cell orbitals must be the first in
!                              orbital lists, i.e. indxuo.le.nuo, with
!                              nuo the number of orbitals in unit cell
!integer nk                  : Number of k points
!real*8  kpoint(3,nk)        : k point vectors
!real*8  wk(nk)              : k point weights (must sum one)
!integer nuotot              : total number of orbitals in unit cell
!                              over all processors
!integer nuo		      : The number of orbitals stored on this Node - if zero
!                              on input then calculated otherwise left unchanged
!
!integer                     : istpmove
!*************************** OUTPUT **********************************
!real*8 Dnew(maxnd,nspin)    : Output Density Matrix in the new basis
!*************************** UNITS ***********************************
!xij and kpoint must be in reciprocal coordinates of each other.
!Enew returned in the units of H.
!delt in femtoseconds
!*************************** Parallel ********************************
! Very important!!!!!!!!!!!!
! DSP: This subroutine is not yet prepared to run in parallel!!!!
! Sorry!
!When running in parallel some of the dimensions are now the
!maximum per node and the corresponding number passed in as
!an argument is the number of locally stored values. The
!variables for which this is the case are:
!
!maxuo/no
!
!*********************************************************************



  use precision
  use parallel,     only : Node, Nodes,BlockSize
  use parallelsubs, only : GlobalToLocalOrb, GetNodeOrbs
  use fdf
  use sys, only: die
#ifdef MPI
  use mpi_siesta,   only : mpi_bcast, mpi_comm_world, mpi_logical
#endif
  use wavefunctions
  use m_diagon, only : ictxt
  use MatrixSwitch
  use matdiagon,    only: geteigen 
  !
  implicit none
  !
  integer                 :: maxnd, maxnh, maxspn, maxuo, maxo, nk, no, nspin, nuotot,istpmove,ierror
  integer                 :: indxuo(no), listh(maxnh), numh(*), listd(maxnd), numd(*)
  integer                 :: listhptr(*), listdptr(*)
  double precision        :: Dnew(maxnd,nspin),kpoint(3,nk),S(maxnh), wk(nk),xij(3,maxnh)
  logical                 :: gamma
  external                :: io_assign, io_close

#ifdef MPI
  integer                 :: MPIerror,desch(9)
  logical, save           :: ParallelOverK
  external                :: diagkp
#endif

  logical, save           :: frstme = .true.
  integer                 :: io, iuo, iu, naux, nhs,  nuo, juo, jo, ind, ispin, nocc, nwf, ik, j, nd
  real(dp)                :: skxij,ckxij, kxij, qe,Ddense(nuotot,nuotot)
  complex(dp)             :: pipj
  !
  type(matrix)                     :: Maux,invsqS,phi, Sauxms
  type(matrix),allocatable,save    :: sqrtS(:)
  character(3)                     :: m_operation
  character(5)                          :: m_storage
  complex(dp)                      :: varaux,varaux2,varaux3, Saux(nuotot,maxuo), psi (nuotot,maxuo)
  !
#ifdef MPI
  call GetNodeOrbs(nuotot,Node,Nodes,nuo)
  if (frstme) then
    if (Node.eq.0) then
      ParallelOverK = fdf_boolean( 'Diag.ParallelOverK', .false. )
    endif
    if(ParallelOverK) then
      stop "chgbasis: tddft-siesta not parallelized over k-points."
    else
      call MPI_Bcast(ParallelOverK,1,MPI_logical,0,MPI_Comm_World,MPIerror)
    endif
      call ms_scalapack_setup(Nodes,1,'c',BlockSize)
  endif
#else
  Node = 0
  Nodes = 1
  nuo = nuotot
#endif
  nhs=nuotot*nuotot*2
  call timer( 'chgbasis', 1 )
#ifdef MPI
  if (ParallelOverK) then
    nhs  = 2 * nuotot * nuotot
  endif
#endif
  !
#ifdef MPI
  m_storage='pzdbc'
  m_operation='lap'
#else
  m_storage='szden'
  m_operation='lap'
#endif
  ! 
  call m_allocate(Sauxms,nuotot,nuotot,m_storage)
  call m_allocate(Maux,nuotot,nuotot,m_storage)
  call m_allocate(invsqS,nuotot,nuotot,m_storage)
  if(frstme) then
    allocate(sqrtS(nk))
    do ik=1,nk
      call m_allocate(sqrtS(ik),nuotot,nuotot,m_storage)
    end do
    frstme=.false.
  endif
  if(istpmove.gt.0) then
    nd = listdptr(nuo) + numd(nuo)
    Dnew(1:nd,1:nspin) = 0.d0
  endif
#ifdef MPI
  call descinit(desch,nuotot,nuotot,BlockSize,BlockSize,0,0,ictxt,nuotot,ierror)
#endif
  ! 
  do ik = 1,nk
    Saux=0.0_dp
    call m_set(Sauxms,'a',cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)
    call timer( 'SsparseTodense', 1 )
    do iuo = 1,nuo
      do j = 1,numh(iuo)
        ind = listhptr(iuo) + j
        jo = listh(ind)
        juo = indxuo(jo)
        if(.not.gamma) then 
          kxij = kpoint(1,ik) * xij(1,ind) +&
          kpoint(2,ik) * xij(2,ind) +&
          kpoint(3,ik) * xij(3,ind)
          ckxij = cos(kxij)
          skxij = -sin(kxij)
        else 
          ckxij=1.0_dp
          skxij=0.0_dp
        endif
        ! Saux=S*e^-ikx, and passing sparse to dense
        Saux(juo,iuo)=Saux(juo,iuo) + cmplx(S(ind)*ckxij,S(ind)*skxij,dp)
      enddo
    enddo
    !
    call timer( 'SsparseTodense', 2 )
    !
    call timer( 'SdenseToMS', 1 )
    !
    do io=1,nuotot
      do j=1,nuotot
#ifdef MPI
        call pzelget('a',' ',varaux2,Saux,j,io,desch)
#else
        varaux2=Saux(j,io)
#endif
        call m_set_element( Sauxms,j,io,varaux2,m_operation)
      enddo
    enddo
    !
    call timer( 'SdenseToMS', 2 )
    !
    if(istpmove.eq.1) then   ! istpmove 
      ! If first step calculate S0^1/2 and save for next step. 
      call calculatesqrtS(Sauxms,invsqS,sqrtS(ik),m_storage,m_operation)
    elseif(istpmove.gt.1) then 
      ! Calculate both Sn^1/2 and Sn^-1/2 where Sn^1/2 is used in n+1 step. 
      call timer( 'S2halfs', 1 )
      call calculatesqrtS(Sauxms,invsqS,Maux,m_storage,m_operation)
      !Saux= Sn-1^1/2*Sn^-1/2
      call mm_multiply(invsqS,'n',sqrtS(ik),'n',&
      Sauxms,cmplx(1.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)
      ! Passing Sn^1/2 from Maux to sqrtS for next step usage.
      call m_add(Maux,'n',sqrtS(ik),cmplx(1.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)
      call timer( 'S2halfs', 2 )
      ! C1=S0^1/2*S1^1/2*C0
      qe=2.0d0*wk(ik)/dble(nspin)
      do ispin=1,nspin
        ! Cn = Saux*Cn-1 where Saux= Sn-1^1/2*Sn^-1/2
        call m_allocate (phi,wavef_ms(ik,ispin)%dim1,wavef_ms(ik,ispin)%dim2,m_storage)
        call m_add(wavef_ms(ik,ispin),'n',phi,cmplx(1.0,0.0,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)
        call timer( 'SauxCn', 1 )
        call mm_multiply(Sauxms,'n',phi,'n',               &
             wavef_ms(ik,ispin),cmplx(1.0_dp,0.0_dp,dp),     &
             cmplx(0.0_dp,0.0_dp,dp),m_operation)
        call timer( 'SauxCn', 2 )
        call m_deallocate(phi)
        ! Constructing  the new DM
        call timer('DMinMS-CB',1)
        call mm_multiply(wavef_ms(1,ispin),'n',wavef_ms(1,ispin),'c',Maux,cmplx(1.0_dp,0.0_dp,dp),&
             cmplx(0.0_dp,0.0_dp,dp),m_operation)
        call timer('DMinMS-CB', 2)
        call timer( 'dmMStodense', 1)
        !  
        do io=1,nuotot
          do jo = 1,nuotot
            call m_get_element(Maux,io,jo,varaux,m_operation)
#ifdef MPI
            call  pzelset(psi,io,jo,desch,varaux)
#else
            psi(io,jo)=varaux
#endif
          end do
        end do
        !
        call timer( 'dmMStodense', 2)
        call timer( 'dmDensetoSparse',1)
        do iuo = 1,nuo
          do j = 1,numd(iuo)
            ind = listdptr(iuo) + j
            jo = listd(ind)
            juo = indxuo(jo)
            if(.not.gamma) then 
              kxij = kpoint(1,ik) * xij(1,ind) +&
              kpoint(2,ik) * xij(2,ind) +&
              kpoint(3,ik) * xij(3,ind)
              ckxij = cos(kxij)
              skxij = -sin(kxij)
            else
              ckxij=1.0d0
              skxij=0.0d0
            endif
            varaux=real(psi(juo,iuo))*ckxij+ aimag(psi(juo,iuo))*skxij
            Dnew(ind,ispin)=Dnew(ind,ispin)+varaux
          enddo
        enddo
        !        
        call timer( 'dmDensetoSparse',2)
      enddo  
    endif   !istpmove 
  enddo          ! ik 
  !
  call timer('chgbasis',2)
  end subroutine chgbasis

 subroutine calculatesqrtS(S,invsqS,sqrtS,m_storage,m_operation)
 
 use precision 
 use matdiagon
 use MatrixSwitch
 ! 
 implicit none
 ! 
 character(5), intent(in)                      :: m_storage
 character(3), intent(in)                      :: m_operation
 type(matrix), intent(inout)                   :: S,invsqS,sqrtS
 complex(kind=dp),dimension(S%dim1,S%dim1)     :: A,B
 type(matrix)                                  :: StimesD
 complex(dp)                                   :: eig, varaux
 real(dp), allocatable                         :: eigen(:)
 integer                                       :: no, info, i, j
 real(dp)  tiny
 data tiny  /1.0d-10/
 ! 
 no=S%dim1 
 allocate(eigen(no))
 call m_allocate(StimesD,no,no,m_storage)
 call geteigen(S,eigen,m_operation)
 call m_set(StimesD,'a',cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)
 ! 
 do j=1,no
   eig=dsqrt(dabs(eigen(j)))
   do i=1,no
     call m_get_element(S,i,j,varaux,m_operation)
     call m_set_element(StimesD,i,j,eig*varaux,m_operation)
   enddo
 enddo 
 !
 call mm_multiply(StimesD,'n',S,'c',sqrtS,cmplx(1.0_dp,0.0_dp,dp),&
      cmplx(0.0_dp,0.0_dp,dp),m_operation)
 do j=1,no
   eig=dsqrt(dabs(eigen(j)))
   eig=1.0d0/(eig+tiny)
   do i=1,no
     call m_get_element(S,i,j,varaux,m_operation)
     call m_set_element(StimesD,i,j,eig*varaux,m_operation)
   enddo
 enddo
 !
 call mm_multiply(StimesD,'n',S,'c',invsqS,cmplx(1.0_dp,0.0_dp,dp),&
      cmplx(0.0_dp,0.0_dp,dp),m_operation)
 call m_deallocate(StimesD)
 deallocate(eigen)
 !
 end subroutine calculatesqrtS
