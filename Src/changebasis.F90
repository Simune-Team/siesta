      subroutine chgbasis(no, nspin, maxspn, maxuo, maxnh, maxnd,            &
                          maxo, gamma, indxuo, nk, kpoint, wk, Dnew,         &
                          nuotot,istpmove)

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
  use parallel,            only : Node, Nodes,BlockSize
  use parallelsubs,        only : GlobalToLocalOrb, GetNodeOrbs,               &
                                  LocalToGlobalOrb
  use fdf
  use alloc
  use sys, only: die
#ifdef MPI
  use mpi_siesta,          only : mpi_bcast, mpi_comm_world, mpi_logical
#endif
  use wavefunctions
  use sparse_matrices,     only : numh, listhptr, listh, S, xijo
!  use densematrix,         only : Saux, psi
  use m_diagon,            only : ictxt
  use MatrixSwitch
  use matdiagon,           only: geteigen 
  !
  implicit none
  !
  integer, intent(in)     :: no, nspin, maxspn, maxuo, maxnh, maxnd, maxo
  integer, intent(in)     :: indxuo(no), nk, nuotot, istpmove
  logical, intent(in)     :: gamma
  !
  real(dp), intent(in)       :: kpoint(3,nk), wk(nk)
  real(dp), intent(out)      :: Dnew(maxnd,nspin) 
!  real(dp), intent(inout)    :: Saux(nuotot,maxuo), psi(nuotot, maxuo)
  !
#ifdef MPI
  integer                 :: MPIerror,desch(9)
  logical, save           :: ParallelOverK
  external                :: diagkp
#endif
  !
  external                :: io_assign, io_close
  !
  logical, save           :: frstme = .true.
  integer                 :: io, iuo, iu, naux, nhs,  nuo, juo, jo, ind, &
                             ispin, nocc, nwf, ik, j,jio, nd, ierror, npsi
  real(dp)                :: skxij,ckxij, kxij, qe
  complex(dp)             :: pipj, varaux,varaux2,varaux3
  !complex(dp), allocatable :: Sx(:,:), psix(:,:)
  !
  type(matrix)             :: Maux,invsqS,phi
  type(matrix)             :: Sauxms
  type(matrix),allocatable,save    :: sqrtS(:)
  character(3)                     :: m_operation
  character(5)                     :: m_storage
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
  if(nspin .le. 2 .and. gamma) then
    nhs  = nuotot * nuo 
    npsi = nuotot * maxuo * nspin
  else if (nspin .le. 2 .and. .not. gamma) then
    nhs  = 2 * nuotot * nuotot
    npsi = 2 * nuotot * nuotot
  else if (nspin .eq. 4) then
    call die ('chgbasis: ERROR: EID not yet prepared for non-collinear spin')
  else 
    call die ('chgbasis: ERROR: incorrect value of nspin')
  end if 
! Allocate local arrays
!  allocate (Sx(nuotot, nuo))
!  allocate (psix(nuotot, nuo))

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
    nd = listhptr(nuo) + numh(nuo)
    Dnew(1:nd,1:nspin) = 0.d0
  endif
#ifdef MPI
  call descinit(desch,nuotot,nuotot,BlockSize,BlockSize,0,0,ictxt,nuotot,ierror)
#endif
  ! 
  do ik = 1,nk
!    Sx(1:nuotot,1:nuo)=0.0_dp
!    call m_set(Sauxms,'a',cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)
    call timer( 'S2MSdense', 1 )
    do iuo = 1,nuo
      call LocalToGlobalOrb(iuo, Node, Nodes, io)
      do j = 1,numh(iuo)
        ind = listhptr(iuo) + j
        jo = listh(ind)
        juo = indxuo(jo)
        if(.not.gamma) then 
          kxij = kpoint(1,ik) * xijo(1,ind) +&
          kpoint(2,ik) * xijo(2,ind) +&
          kpoint(3,ik) * xijo(3,ind)
          ckxij = cos(kxij)
          skxij = -sin(kxij)
        else 
          ckxij=1.0_dp
          skxij=0.0_dp
        endif
        ! Saux=S*e^-ikx, and passing sparse to dense
        !Sx(juo,iuo)=Sx(juo,iuo) + cmplx(S(ind)*ckxij,S(ind)*skxij,dp)
        varaux2 = cmplx(S(ind)*ckxij,S(ind)*skxij)
        call m_set_element(Sauxms, jo, io, varaux2, m_operation)
      enddo
    enddo
    !
    call timer( 'S2MSdense', 2 )
    !
!    call timer( 'SdenseToMS', 1 )
    !
  !  do io=1,nuo
  !    do j=1,nuotot
#ifdef MPI
   !     call LocalToGlobalOrb(io,Node, Nodes, jo)
   !     varaux2 = Sx(j,io)
        !call pzelget('a',' ',varaux2,Sx,j,io,desch)
#else
   !     jo = io
   !     varaux2=Sx(j,io)
#endif
    !    call m_set_element( Sauxms,j,jo,varaux2,m_operation)
    !  enddo
    !enddo
    !
 !   call timer( 'SdenseToMS', 2 )
    !
    if(istpmove.eq.1) then   ! istpmove 
      ! If first step calculate S0^1/2 and save for next step. 
      call calculatesqrtS(Sauxms,invsqS,sqrtS(ik),nuo,m_storage,m_operation)
    elseif(istpmove.gt.1) then 
      ! Calculate both Sn^1/2 and Sn^-1/2 where Sn^1/2 is used in n+1 step. 
      call timer( 'S2halfs', 1 )
      call calculatesqrtS(Sauxms,invsqS,Maux,nuo,m_storage,m_operation)
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
        call mm_multiply(wavef_ms(1,ispin),'n',wavef_ms(1,ispin),'c',          &
                         Maux,cmplx(1.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp), &
                         m_operation)
        call timer('DMinMS-CB', 2)
!        call timer( 'dmMStodense', 1)
        !  
!        do io=1,nuotot
!#ifdef MPI
!            call LocalToGlobalOrb(io, Node, Nodes, j)
!#endif
!          do jo = 1,nuotot
#ifdef MPI
!            call m_get_element(Maux,jo,io,varaux,m_operation)
!            call  pzelset(psix,jo,io,desch,varaux)
#else
!            call m_get_element(Maux,jo,io,varaux,m_operation)
!            psix(jo,io)=varaux
#endif
!         end do
!        end do
        !
!        call timer( 'dmMStodense', 2)
        call timer( 'dmDensetoSparse',1)
        do iuo = 1,nuo
          do j = 1,numh(iuo)
            ind = listhptr(iuo) + j
            jo = listh(ind)
            juo = indxuo(jo)
            if(.not.gamma) then 
              kxij = kpoint(1,ik) * xijo(1,ind) +&
              kpoint(2,ik) * xijo(2,ind) +&
              kpoint(3,ik) * xijo(3,ind)
              ckxij = cos(kxij)
              skxij = -sin(kxij)
            else
              ckxij=1.0d0
              skxij=0.0d0
            endif
            varaux2=real(Maux%zval(jo,iuo))*ckxij+ aimag(Maux%zval(jo,iuo))*skxij
            Dnew(ind,ispin)=Dnew(ind,ispin)+varaux2
          enddo
        enddo
        !        
        call timer( 'dmDensetoSparse',2)
      enddo  
    endif   !istpmove 
  enddo          ! ik 
  !
 ! deallocate(Sx)
 ! deallocate(psix)
  call m_deallocate(Sauxms)
  call m_deallocate(Maux)
  call m_deallocate(invsqS)

  call timer('chgbasis',2)
  end subroutine chgbasis

 subroutine calculatesqrtS(S,invsqS,sqrtS,nu,m_storage,m_operation)
 
 use precision 
 use matdiagon
 use MatrixSwitch
 use parallelsubs,          only: LocalToGlobalOrb
 use parallel,              only: Node, Nodes
 ! 
 implicit none
 ! 
 character(5), intent(in)                      :: m_storage
 character(3), intent(in)                      :: m_operation
 type(matrix), intent(inout)                   :: S,invsqS,sqrtS
 type(matrix)                                  :: SD01, SD02
 complex(dp)                                   :: varaux
 real(dp)                                      :: eig01, eig02
 real(dp), allocatable                         :: eigen(:)
 integer                                       :: no,nu, info, i, j,jo
 real(dp)  tiny
 data tiny  /1.0d-10/
 ! 
 no=S%dim1 
 allocate(eigen(no))
 call m_allocate(SD01,no,no,m_storage)
 call m_allocate(SD02,no,no,m_storage)
 !call m_set(SD01,'a',cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)
 !call m_set(SD02,'a',cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)
 ! 
 ! Takes overlap matrix S in dense form and returns its eigenvalues
 ! in eigen(*) and eigenvectors in S.
 call geteigen(S,eigen,m_operation)
 !
 do j=1,nu
   call LocalToGlobalOrb(j,Node,Nodes,jo)
   eig01=dsqrt(dabs(eigen(jo)))
   eig02=1.0d0/(eig01+tiny)
   do i=1,no
     varaux = S%zval(i,j)
     !call m_get_element(S,i,j,varaux,m_operation)
     call m_set_element(SD01,i,jo,eig01*varaux,m_operation)
     call m_set_element(SD02,i,jo,eig02*varaux,m_operation)
   enddo
 enddo 
 !
 deallocate(eigen)
 
 call mm_multiply(SD01,'n',S,'c',sqrtS,cmplx(1.0_dp,0.0_dp,dp),&
                  cmplx(0.0_dp,0.0_dp,dp),m_operation)
 call mm_multiply(SD02,'n',S,'c',invsqS,cmplx(1.0_dp,0.0_dp,dp),&
                  cmplx(0.0_dp,0.0_dp,dp),m_operation)
 call m_deallocate(SD01)
 call m_deallocate(SD02)
 !
 end subroutine calculatesqrtS
