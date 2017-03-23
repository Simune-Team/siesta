      subroutine evolg( nspin, nuo, no, maxo, maxnh, maxnd,           &
                        Enew, nuotot, delt, ist,itd)
! ********************************************************************
! Subroutine to calculate the eigenvalues and eigenvectors, density
! and energy-density matrices, and occupation weights of each 
! eigenvector, for given Hamiltonian and Overlap matrices (including
! spin polarization). Gamma-point version.
! Written by A. Tsolakidis, May 2000 after a subroutine
! by J. M. Soler.
! Rewritten by D. Sanchez-Portal, November 2002-March 2003
! Rewritten by Rafi Ullah and Adiran Garaizar June-October 2015
! Making it parallel using Matrix Switch.
! **************************** INPUT **********************************
! integer nspin               : Number of spin components (1 or 2)
! integer nuo                 : Number of basis orbitals local to node
! integer no                  : Number of basis orbitals
! integer maxo                : Maximum number of orbitals in the unit cell
! integer maxnh               : Maximum number of orbitals interacting  
! integer maxnd               : Maximum number of nonzero elements of 
!                               each row of density matrix
! integer numh(nuo)           : Number of nonzero elements of each row 
!                               of hamiltonian matrix
! integer listhptr(nuo)       : Pointer to each row (-1) of the
!                               hamiltonian matrix
! integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
!                               column indexes for each matrix row
! integer numd(nuo)           : Number of nonzero elements of each row 
!                               ofdensity matrix
! integer listdptr(nuo)       : Pointer to each row (-1) of the
!                               density matrix
! integer listd(maxnd)        : Nonzero density-matrix element column 
!                               indexes for each matrix row
! real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
! real*8  S(maxnh)            : Overlap in sparse form
! integer nuotot              : total number of orbitals per unit cell
!                               over all processors
! real*8  delt                : length of the time step
! ******************** OUTPUT **************************************
! real*8 Dnew(maxnd,nspin)    : Output New Density Matrix
! real*8 Enew(maxnd,nspin)    : Output New Energy-Density Matrix
! real*8 eo(maxo,nspin,1)       : Output instantaneous eigenvalues
!                              (only calculated if explicitly required
!                               by user and in the last MD step)
! New wavefunctions are calculated and stored for the next step.
! *************************** AUXILIARY *******************************
! real*8 Haux(nuotot,nuo)     : Auxiliary space for the hamiltonian matrix
! real*8 Saux(nuotot,nuo)     : Auxiliary space for the overlap matrix
! real*8 aux(2*nuotot)        : Extra auxiliary space
! *************************** UNITS ***********************************
! Enew returned in the units of H.
! *************************** PARALLEL ********************************
! The auxiliary arrays are now no longer symmetry and so the order
! of referencing has been changed in several places to reflect this.
! *********************************************************************
!
!  Modules
!
      use precision
      use sys
      use parallel    
      use parallelsubs,          only: LocalToGlobalOrb
      use alloc
      use wavefunctions
      use MatrixSwitch
      use siesta_options,        only: eigen_time, ntded, extrapol_H_tdks, ntded_sub
      use sparse_matrices,       only: H, S, numh, listh, listhptr, Dscf 
      use m_eo,                  only: eo
      use m_steps,               only: fincoor, final
#ifdef MPI
      use mpi_siesta
#endif
      !
      implicit none
#ifdef MPI
      INTEGER :: MPIerror
#endif
      !
      integer              :: itd, ist, maxnd, maxnh, nuo, no, nspin, nuotot
      integer              :: ncounter, maxo, asn,desch(9)
      real(dp)             :: Enew(maxnd,nspin), delt
      !
      type(matrix)         :: Hauxms,Sauxms
      character(3)         :: m_operation
      character(5)         :: m_storage
      complex(dp)          :: varaux, varaux2,varaux3, pipj, pipj2, varaux4  
       
      integer              :: ie, io, iio,iee, ispin, j, jo, BNode, iie, ind, BTest
      integer              :: mm, maxnuo, ierror, nd, nocc, nstp,i,npsi
      real(dp)             :: qe, t, eigv, dnrm,el1,el2,el3,el4
      logical              :: calculateEnew ! Not sure if it is really needed?
      !
#ifdef MPI
      m_storage='pzdbc'
      m_operation='lap'
#else
      m_storage='szden'
      m_operation='lap'
#endif
      npsi=2*nuotot*nuotot
      calculateEnew = .true.        ! Why? 
      !
      call timer( 'Evolg', 1 )
      !
      call m_allocate( Hauxms,nuotot,nuotot,m_storage)
      call m_allocate( Sauxms,nuotot,nuotot,m_storage)
      !
      nstp=1 ! I guess this determines the Hamiltonia extrapolation.
             ! Needs to be checked and should be made user defined.
      !
      !if(calculateEnew) Enew(1:nd,1:nspin) = 0.d0
      ! Evolve wavefunctions.............................................
      do ispin = 1,nspin
        ncounter=0
        if(ispin.eq.2) ncounter=wavef_ms(1,1)%dim2
        ! One can use the dense overlap matrix constructed in changebais? 
        call timer( 'HSSparseToDense', 1 )
        !
        do io = 1,nuo
          call LocalToGlobalOrb (io,Node, Nodes, i)
          do j = 1,numh(io)
            ind = listhptr(io) + j
            jo = listh(ind)
             call m_set_element(Hauxms, jo, i, H(ind,ispin), m_operation)
             call m_set_element(Sauxms, jo, i, S(ind), m_operation)
          enddo
        enddo
        !
        call timer( 'HSSparseToDense', 2 )
        !
        call timer( 'CntoCn1', 1 )
        !
        call evol1new(Hauxms, Sauxms, nuotot, nuo, nspin,ispin, ncounter,              &
             delt,extrapol_H_tdks,ntded_sub)
        nocc=wavef_ms(1,ispin)%dim2
        !
        call timer( 'CntoCn1', 2 )
        !
        ! Calculating denisty matrix.
        IF (IONode .and. ispin .eq. 1) THEN
          WRITE(*,'(a)') 'evolg: Computing DM after evolving TDKS wavefunctions'
        END IF
        call compute_tddm(ispin, Dscf)
        ! This needs to be fixes as it is not being properly done.
        call timer( 'Eigenvalue', 1 )
        if (eigen_time) then 
          nocc=wavef_ms(1,ispin)%dim2
          do ie = 1,nocc
            eigv=0.0d0
            do io = 1,nuotot
              do jo = 1,nuotot
                call m_get_element(wavef_ms(1,ispin),io,ie,varaux,m_operation)
                call m_get_element(wavef_ms(1,ispin),jo,ie,varaux2,m_operation)
                pipj = real(varaux)*real(varaux2) + aimag(varaux)*aimag(varaux2)
                call m_get_element(Hauxms,jo,io,varaux3,m_operation)
                eigv=eigv+real(varaux3)*pipj
              enddo
            enddo
            eo(ie,ispin,1)=eigv
          enddo
        endif
        !
        call timer( 'Eigenvalue', 2 )
      enddo ! ispin
      call m_deallocate(Hauxms)
      call m_deallocate(Sauxms)
      !
      call timer( 'Evolg', 2 )
    END SUBROUTINE evolg
!---------------------------------------------------------------------------------!
    subroutine Uphi(H, S, phi, no, nocc, deltat)
!*************************************************************************
!Subroutine that calculates the new wavefunction, given the old 
!wavefunction by using the formula for the time evolution. Gamma-point 
!version. Written by A. Tsolakidis, May 2000 
!Modified by D. Sanchez-Portal, July 2002
!Modified by D. Sanchez-Portal,  2008. 
!This version is limited to first order expansion
!and avoids the inversion of the overlap.
!Modified by Rafi Ullah, October, 2015.
!Parallelized by using Matrix Swtich.
!*************************** INPUT ***************************************
!integer no                  : Number of basis orbitals
!integer nol                 : Local number of basis orbitals
!integer noccok              : Number of occupied wavefunctions per spin
!real*8 H(no,nol)             : Hamiltonian matrix
!real*8 S(no,nol)             : Overlap matrix
!complex*16  Phi(no,nocc) : Old wavefunctions
!*************************** OUTPUT *************************************
!complex*16  Phi(no,nocc) : New wavefunctions
!*************************** AUXILIARY ********************************** 
!complex*16 Q_1(no_max,no_max)           : Auxiliary Matrix
!complex*16 Q_2(no_max,no_max)           : Auxiliary Matrix
!complex*16 Q_3(no_max,no_max)           : Auxiliary Matrix
!real*8 deltat                           : Duration of the time step 
!*************************************************************************
! Modules
!
      use precision
      use MatrixSwitch 
      use matswinversion,   only: getinverse
#ifdef MPI
      use mpi_siesta
#endif
      use parallel
!*************************************************************************   
  
 implicit none 
 !     
 integer               :: no, ispin, nocc,nuo
 complex(kind=dp)      :: pi, pj
 real(kind=dp)         :: deltat
 type(matrix)          :: H,S,phi
 integer   ,  dimension(:), allocatable, save   :: ipiv
 ! Internal variables 
 integer               :: i, j , k, info, no2, l
 type(matrix)          :: aux1,aux2,aux3,aux4
 complex(kind=dp)      :: alpha, ss, hh,varaux,varaux2
 character             :: m_storage*5, m_operation*3
 logical, save         :: onlyelectrons 
 logical, save         :: frsttime = .true. 
 !
#ifdef MPI
 m_storage='pzdbc'
 m_operation='lap'
#else
 m_storage='szden'
 m_operation='lap'
#endif
 no2=no*no
 !  
 call m_allocate(aux1,no,no,m_storage)
 call m_allocate(aux2,no,no,m_storage)
 call m_allocate(aux3,no,no,m_storage)
 call m_allocate(aux4,phi%dim1,phi%dim2,m_storage)
 ! First order expansion for the evolution operator
 alpha=-0.5_dp*cmplx(0.0_dp,1.0_dp,dp)*deltat
 ! Copying S to aux1 and aux2
 call m_add(S,'n',aux1,cmplx(1.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)
 call m_add(S,'n',aux2,cmplx(1.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)
 ! copying phi_0 (wavefunction) to aux4
 call m_add(phi,'n',aux4,cmplx(1.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)
 ! Calculating S - alpha * H
 call m_add(H,'n',aux1,alpha,cmplx(1.0_dp,0.0_dp,dp),m_operation)
 ! Calculating S + alpha * H
 call m_add(H,'n',aux2,-1.0_dp*alpha,cmplx(1.0_dp,0.0_dp,dp),m_operation)
 ! Calculating inverse of (S + alpha * H)
 call getinverse(aux2)
 ! aux3 = (S + alpha * H)^-1 * (S - alpha * H)
 call mm_multiply(aux2,'n',aux1,'n',aux3,cmplx(1.0,0.0,dp),cmplx(0.0,0.0,dp),m_operation)
 ! phi_1 = aux3 * phi_0     
 call mm_multiply(aux3,'n',aux4,'n',phi,cmplx(1.0,0.0,dp),cmplx(0.0,0.0,dp),m_operation)
 !
 call m_deallocate(aux1)
 call m_deallocate(aux2)
 call m_deallocate(aux3)
 call m_deallocate(aux4)
 !
 END SUBROUTINE Uphi
 !------------------------------------------------------------------------------------!
 SUBROUTINE evol1new(Hauxms, Sauxms, no, nol, nspin, ispin,         & 
            ncounter, delt, extrapol, nstp)
!*************************************************************************
!Subroutine that calculates the new wavefunction, given the old 
!wavefunction by using the formula for the time evolution. Gamma-point 
!version. Written by A. Tsolakidis, May 2000 
!Modified by D. Sanchez-Portal, July 2002
!This version is limited to first order expansion
!and avoids the inversion of the overlap
!Re-written by Rafi Ullah and Adiran Garaizar June-October 2015.
!Making it parallel using Matrix Swtich.
!*************************************************************************
! Modules
!
      use precision
      use fdf
      use wavefunctions, only: wavef_ms
      use parallel
      use MatrixSwitch
#ifdef MPI
      use mpi_siesta
#endif
!*************************************************************************   
  implicit none 
  !
  integer                 :: no,  ispin, ncounter, nol, nstp, nspin,nuo
  complex(kind=dp)        :: pi, pj
  real(dp)                :: deltat, delt, varaux
  !
  type(matrix),intent(in)         :: Hauxms, Sauxms 
  type(matrix),allocatable,save   :: Hsve(:)
  character(5)                    :: m_storage
  character(3)                    :: m_operation
  logical                         :: extrapol
  ! Internal variables ...
  integer                :: i, j , k, info, no2, nocc, l
  complex(dp)            ::  alpha,hh, sum
  logical, save          :: fsttim(2) = (/.true. , .true./)
  logical, save          :: frsttime = .true.
  logical, save          :: onlyelectrons  = .false.
  save                   ::  deltat
  !
#ifdef MPI
  m_storage='pzdbc'
  m_operation='lap'
#else
  m_storage='szden'
  m_operation='lap'
#endif
  no2=no*no
  if (frsttime) then
    !nstp is the number of "substeps" in the electronic evolution
    !the evolution operator is applied in each substep although
    !an extrapolated Hamiltonian is used "rather" than 
    !a SCF Hamiltonian
    deltat=delt/0.04837d0/dble(nstp)
    if (Node.eq.0) then
      write(6,*) 'evol1: time step (Ry**-1) ',deltat
    end if
    allocate(Hsve(nspin))
    do i=1, nspin
      call m_allocate(Hsve(i),no,no,m_storage)
    end do 
    frsttime=.false.
  endif     ! frsttime
  nocc=wavef_ms(1,ispin)%dim2
  call timer('evol1.xtpl',1)
  !
  do l=1,nstp
    if(fsttim(ispin).or..not.extrapol) then
      call Uphi(Hauxms, Sauxms, wavef_ms(1,ispin), no, nocc, deltat)
    else
      varaux=(l-0.5_dp)/dble(nstp)
      call m_add(Hauxms,'n',Hsve(ispin),cmplx(1.0,0.0,dp),cmplx(-1.0,0.0,dp),m_operation)
      call m_add(Hauxms,'n',Hsve(ispin),cmplx(1.0,0.0,dp),cmplx(varaux,0.0,dp),m_operation) 
      call Uphi(Hsve(ispin), Sauxms, wavef_ms(1,ispin),no, nocc, deltat)
    endif
  enddo 
  fsttim(ispin)=.false.
  !
  !Storing Hamitonian for extrapolation and later correction    
  !
  call m_add(Hauxms,'n',Hsve(ispin),cmplx(1.0,0.0_dp,dp),cmplx(0.0,0.0,dp),m_operation)
  call timer('evol1.xtpl',2)
  END SUBROUTINE evol1new
!---------------------------------------------------------------------------------------!
      subroutine applyinverSH(S,H,no,nol,psi,psi2)
!**********************************************************************************!
! Re-written by Rafi Ullah and Adiran Garaizar June-October 2015                   !
! to make it parallel using Matrix Switch.                                         !
!**********************************************************************************!
      use parallel
      use precision
      use MatrixSwitch
      use matswinversion,   only: getinverse
      !
      implicit none
      !
      integer                :: i, j , k, info, no2,no, nol
      logical, save          :: frstime = .true.
      character              :: m_storage*5, m_operation*3
      !
      type(matrix)           :: S, H,psi,psi2
      type(matrix)           ::  aux1, S_1
      !
#ifdef MPI
      m_storage='pzdbc'
      m_operation='lap'
#else
      m_storage='szden'
      m_operation='lap'
#endif
      no2=no*no
        call m_allocate(aux1,no,no,m_storage)
        call m_allocate(S_1,no,no,m_storage)
      ! Invert the overlap matrix.
      call m_add (S,'n',S_1,cmplx(1.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)
      call getinverse(S_1,m_operation)
      call mm_multiply(H,'n',S_1,'n',aux1,cmplx(1.0_dp,0.0_dp,dp),                         &
           cmplx(0.0_dp,0.0_dp,dp),m_operation)
      call mm_multiply(aux1,'n',psi,'t',psi2,cmplx(1.0_dp,0.0_dp,dp),                      &
           cmplx(0.0_dp,0.0_dp,dp),m_operation)
      call m_deallocate(S_1)
      call m_deallocate(aux1)
      end subroutine applyinverSH
