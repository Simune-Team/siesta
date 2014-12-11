
      subroutine chgbasis(no, nspin, maxspn, maxuo, maxnh, maxnd,
     .                  maxo, numh, listhptr, listh, numd,
     .                  listdptr, listd, S, 
     .                  gamma, xij, indxuo, nk, kpoint, wk,
     .                  Dnew, nuotot,istpmove)

C *********************************************************************
C Modified by D. Sanchez-Portal, Feb 2009
C **************************** INPUT **********************************
C integer no                  : Number of basis orbitals the supercell
C integer nspin               : Spin polarization (1 or 2)
C integer maxspn              : Maximum number of spin orentations 
C integer maxuo               : Maximum number of orbitals stored in a 
C                               given Node
C integer maxnh               : Maximum number of orbitals interacting
C integer maxnd               : Maximum number of nonzero elements of
C                               each row of density matrix
C integer maxo                : Maximum number of orbitals in the unit cell
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
C integer                     : istpmove
C *************************** OUTPUT **********************************
C real*8 Dnew(maxnd,nspin)    : Output Density Matrix in the new basis
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C Enew returned in the units of H.
C delt in femtoseconds
C *************************** Parallel ********************************
C  Very important!!!!!!!!!!!!
C  DSP: This subroutine is not yet prepared to run in parallel!!!!
C  Sorry!
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
      use densematrix,  only : Haux, Saux
      use alloc
      use m_memory
      use sys, only: die
#ifdef MPI
      use mpi_siesta,   only : mpi_bcast, mpi_comm_world,
     $                         mpi_logical
#endif
      use wavefunctions

      implicit none

      integer
     .  maxnd, maxnh, maxspn, maxuo, maxo, nk, no, nspin, nuotot,
     .   istpmove

      integer
     .  indxuo(no), listh(maxnh), numh(*), listd(maxnd), numd(*),
     .  listhptr(*), listdptr(*)

      double precision
     .  Dnew(maxnd,nspin),
     .  kpoint(3,nk),
     .  S(maxnh), wk(nk),
     .  xij(3,maxnh)

      logical
     .  gamma

      external
     . io_assign, io_close

#ifdef MPI
      integer
     .  MPIerror
      external
     .  diagkp
#endif
      logical
     .  frstme

C  Internal variables .............................................
      integer
     .  io, iuo, iu, naux, nhs,  nuo, juo, jo, ind, ispin, nocc, nwf,
     .  ik, j, nd

      real(dp) skxij,ckxij, kxij, qe

C Dynamic arraysdd
      complex*16, dimension(:,:,:), allocatable, save :: sqrtS
      complex*16, dimension(:,:), allocatable, save :: invsqS
      complex*16 pipj
#ifdef MPI
      logical, save ::
     .  ParallelOverK
#endif
C  ....................

      save frstme
      data frstme /.true./


C Get Node number and calculate local orbital range
#ifdef MPI
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)
      if (frstme) then
        if (Node.eq.0) then
          ParallelOverK = fdf_boolean( 'Diag.ParallelOverK', .false. )
        endif
        call MPI_Bcast(ParallelOverK,1,MPI_logical,0,MPI_Comm_World,
     .    MPIerror)
      endif
#else
      Node = 0
      Nodes = 1
      nuo = nuotot
#endif
C  ....................

C Start time counter ................................................
      call timer( 'chgbasis', 1 )
C ......................

C Check internal dimensions ..........................................
        nhs  =  2 * nuotot * nuo
#ifdef MPI
        if (ParallelOverK) then
          nhs  = 2 * nuotot * nuotot
        endif
#endif

C Allocate local arrays
      call re_alloc(Haux,1,nhs,name='Haux',routine='chgbasis')
      call re_alloc(Saux,1,nhs,name='Saux',routine='chgbasis')
      
      if(frstme) then
C  Allocate memory to store the inverse matrix ................. 
       allocate(sqrtS(nuotot,nuotot,nk))
       call memory('A','Z',nuotot*nuotot*nk,'chgbasis')
       allocate(invsqS(nuotot,nuotot))
       call memory('A','Z',nuotot*nuotot,'chgbasis')
       frstme=.false.
      endif
 
C.....................................................................
      if(istpmove.gt.0) then
         nd = listdptr(nuo) + numd(nuo)
         Dnew(1:nd,1:nspin) = 0.d0
      endif
          
      do ik = 1,nk
        do iuo=1,nhs
            Saux(iuo) = 0.0d0
        enddo
        do iuo = 1,nuo
          do j = 1,numh(iuo)
            ind = listhptr(iuo) + j
            jo = listh(ind)
            juo = indxuo(jo)
            io=2*(juo-1)+2*nuotot*(iuo-1)
            if(.not.gamma) then 
              kxij = kpoint(1,ik) * xij(1,ind) +
     .               kpoint(2,ik) * xij(2,ind) +
     .               kpoint(3,ik) * xij(3,ind)
              ckxij = cos(kxij)
              skxij = -sin(kxij)
            else 
              ckxij=1.0d0
              skxij=0.0d0
            endif
              Saux(io+1)=Saux(io+1)+S(ind)*ckxij
              Saux(io+2)=Saux(io+2)+S(ind)*skxij
          enddo
        enddo
        if(istpmove.eq.0) then 

          call calculatesqrtS(Saux,nuotot,nuo,invsqS,sqrtS(1,1,ik))
        
        elseif(istpmove.gt.0) then 

          
          call calculatesqrtS(Saux,nuotot,nuo,invsqS,Haux)

          call zhemm('l','u',nuotot,nuotot,(1.0d0,0.0d0),
     .       invsqS,nuotot,sqrtS(1,1,ik),nuotot,
     .       (0.0d0,0.0d0),Saux,nuotot)
     
          do iuo=1,nuo
              do juo=1,nuotot
               io=2*(juo-1)+2*nuotot*(iuo-1)
c              write(6,*) iuo,juo, Saux(io+1),Saux(io+2)
               sqrtS(juo,iuo,ik)=cmplx(Haux(io+1),Haux(io+2))
              enddo 
          enddo 

          qe=2.0d0*wk(ik)/dble(nspin)
          do ispin=1,nspin
             nocc=wavef%nocck(ik,ispin)
             do iu=1,nocc
               nwf=wavef%iwf(iu,ik,ispin)
               do iuo=1,nuotot
                 io=2*(iuo-1)+2*nuotot*(iu-1) 
                 Haux(1+io)= wavef%phi(1,iuo,nwf)
                 Haux(2+io)= wavef%phi(2,iuo,nwf)
               enddo 
             enddo 

             call zgemm('n','n',nuotot,nocc,nuotot,
     .           (1.0d0,0.0d0), Saux,nuotot,
     .           Haux,nuotot,(0.0d0,0.0d0),invsqS,nuotot)

          do iu=1,nocc
             nwf=wavef%iwf(iu,ik,ispin)
             do iuo=1,nuotot
               wavef%phi(1,iuo,nwf)=dreal(invsqS(iuo,iu))
               wavef%phi(2,iuo,nwf)=dimag(invsqS(iuo,iu))
             enddo
          enddo

           do iuo =1, nhs
              Haux(iuo)=0.0d0
           enddo

        do iu = 1,nocc
          do iuo = 1,nuo
           do juo = 1,nuotot
              io=2*(juo-1)+2*nuotot*(iuo-1)
              pipj = qe*conjg(invsqS(iuo,iu))*invsqS(juo,iu)
              Haux(io+1) = Haux(io+1) + dreal(pipj)
              Haux(io+2) = Haux(io+2) + dimag(pipj)
           enddo
         enddo
        enddo

          do iuo = 1,nuo
            do j = 1,numd(iuo)
              ind = listdptr(iuo) + j
              jo = listd(ind)
              juo = indxuo(jo)
              io=2*(juo-1)+2*nuotot*(iuo-1)
              if(.not.gamma) then 
                kxij = kpoint(1,ik) * xij(1,ind) +
     .               kpoint(2,ik) * xij(2,ind) +
     .               kpoint(3,ik) * xij(3,ind)
                ckxij = cos(kxij)
                skxij = -sin(kxij)
              else
                ckxij=1.0d0
                skxij=0.0d0
              endif
              Dnew(ind,ispin)=Dnew(ind,ispin)+
     .              Haux(io+1)*ckxij +
     .              Haux(io+2)*skxij
            enddo
          enddo

         enddo  

        endif  

        enddo

        call timer('chgbasis',2)

       end

      subroutine calculatesqrtS(S,no,nol,invsqS,sqrtS)

      use precision

      complex*16 S(no,nol), eig, sqrtS(no,nol),invsqS(no,nol)
      integer lwork, info,i,j
      complex*16,dimension(:),allocatable,save :: work
      real(dp),dimension(:),allocatable,save :: rwork, eigen
      complex*16, dimension(:,:), allocatable, save :: StimesD, aux
      real(dp)  tiny
      logical  frstime
      data frstime /.true./
      save frstime 
      data tiny  /1.0d-10/

      lwork=2*no
      if(frstime) then 
          allocate(eigen(no))
          call memory('A','D',no,'calculatesqrtS')
          allocate(work(lwork))
          call memory('A','Z',lwork,'calculatesqrtS')
          allocate(rwork(3*no))
          call memory('A','D',3*no,'calculatesqrtS')
          allocate(StimesD(no,no))
          call memory('A','Z',no*no,'calculatesqrtS')
c         allocate(aux(no,no))
c         call memory('A','Z',no*no,'calculatesqrtS')
          frstime=.false.
      endif

c     do j=1,no
c        do i=1,no
c           aux(i,j)=S(i,j)
c        enddo 
c     enddo 
c     call zheev('V','U',no,aux,no,eigen,work,lwork,rwork,info)
      call zheev('V','U',no,S,no,eigen,work,lwork,rwork,info)
       
      if(info.ne.0) stop 'Error calculatesqrtS'
      StimesD(1:no,1:no)=(0.0d0,0.0d0)
      do j=1,no
         eig=dsqrt(dabs(eigen(j)))
         do i=1,no
c           StimesD(i,j)=eig*aux(i,j)
            StimesD(i,j)=eig*S(i,j)
         enddo 
      enddo 
       
c     call zgemm('n','t',no,no,no,(1.0d0,0.0d0),StimesD,no,
c    .       aux,no,(0.0d0,0.0d0),sqrtS,no)
      call zgemm('n','c',no,no,no,(1.0d0,0.0d0),StimesD,no,
     .       S,no,(0.0d0,0.0d0),sqrtS,no)

      do j=1,no
         eig=dsqrt(dabs(eigen(j)))
         eig=1.0d0/(eig+tiny)
         do i=1,no
c           StimesD(i,j)=eig*aux(i,j)
            StimesD(i,j)=eig*S(i,j)
         enddo
      enddo

c     call zgemm('n','t',no,no,no,(1.0d0,0.0d0),StimesD,no,
c    .       aux,no,(0.0d0,0.0d0),sqrtS,no)
      call zgemm('n','c',no,no,no,(1.0d0,0.0d0),StimesD,no,
     .       S,no,(0.0d0,0.0d0),invsqS,no)
      
   
      end

