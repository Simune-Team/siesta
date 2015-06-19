
      subroutine evolg( nspin, nuo, no, maxo, maxnh, maxnd, 
     .                  numh, listhptr, listh, numd, 
     .                  listdptr, listd, H, S, eo,
     .                  Dnew, Enew, 
     .                  Haux, Saux, nuotot,
     .                  delt)
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization). Gamma-point version.
C Written by A. Tsolakidis, May 2000 after a subroutine
C by J. M. Soler.
C rewritten by D. Sanchez-Portal, November 2002-March 2003
C Modified by M. Ahsan Zeb, 11 March 2011.
C **************************** INPUT **********************************
C integer nspin               : Number of spin components (1 or 2)
C integer nuo                 : Number of basis orbitals local to node
C integer no                  : Number of basis orbitals
C integer maxo                : Maximum number of orbitals in the unit cell
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : Maximum number of nonzero elements of 
C                               each row of density matrix
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer numd(nuo)           : Number of nonzero elements of each row 
C                               ofdensity matrix
C integer listdptr(nuo)       : Pointer to each row (-1) of the
C                               density matrix
C integer listd(maxnd)        : Nonzero density-matrix element column 
C                               indexes for each matrix row
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C real*8  delt                : length of the time step
C ******************** OUTPUT **************************************
C real*8 Dnew(maxnd,nspin)    : Output New Density Matrix
C real*8 Enew(maxnd,nspin)    : Output New Energy-Density Matrix
C real*8 eo(maxo,nspin)       : Output instantaneous eigenvalues
C                              (only calculated if explicitly required
C                               by user and in the last MD step)
C New wavefunctions are calculated and stored for the next step.
C *************************** AUXILIARY *******************************
C real*8 Haux(nuotot,nuo)     : Auxiliary space for the hamiltonian matrix
C real*8 Saux(nuotot,nuo)     : Auxiliary space for the overlap matrix
C real*8 aux(2*nuotot)        : Extra auxiliary space
C *************************** UNITS ***********************************
C Enew returned in the units of H.
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
      use elec_dyn_options     

      implicit none

#ifdef MPI

#endif
      
      integer
     .  maxnd, maxnh, nuo, no, nspin, nuotot, 
     .  ncounter, maxo, asn

      integer 
     .  listh(maxnh), numh(nuo), listhptr(nuo),
     .  listd(maxnd), numd(nuo), listdptr(nuo)

      double precision
     .  Dnew(maxnd,nspin), Enew(maxnd,nspin),
     .  H(maxnh,nspin),
     .  S(maxnh), delt, eo(maxo,nspin)
     
      complex*16
     .  Haux(nuotot,nuo), Saux(nuotot,nuo)
     
      double precision
     .  aux(2,nuotot), aux2(2,nuotot)

	double precision, allocatable ::

     .  bix(:,:)
	
	complex*16, allocatable :: 

     .  bix2(:,:)



C  Internal variables .............................................
      integer           ie, io, iio, ispin, j, jo, BNode, iie, ind,
     .                  BTest, ierror, nd, nocc, nstp
      double precision  ee, pipj, qe, t, eigv, dnrm
C  ....................

C Number of "substeps" for electron dynamics (in between 
C Hamiltonian reevaluations)
       nstp=Nelecsubsteps

C New density and energy-density matrices of unit-cell orbitals .......
      nd = listdptr(nuo) + numd(nuo)
      Dnew(1:nd,1:nspin) = 0.d0
      if(calculateEnew) Enew(1:nd,1:nspin) = 0.d0

C Evolve wavefunctions.............................................
C Attention!! At the moment this is only prepared for serial work 
      do ispin = 1,nspin
        ncounter=0
        if(ispin.eq.2) ncounter=wavef%nocck(1,1)
   10   do io = 1,nuo
          do jo = 1,nuotot
            Saux(jo,io) = (0.0d0,0.0d0)
            Haux(jo,io) = (0.0d0,0.0d0)
          enddo
        enddo
        do io = 1,nuo
          do j = 1,numh(io)
            ind = listhptr(io) + j
            jo = listh(ind)
            Saux(jo,io) = Saux(jo,io) + cmplx(S(ind),0.0d0)
            Haux(jo,io) = Haux(jo,io) + cmplx(H(ind,ispin),0.0d0)
          enddo
        enddo
        call evol1new(Haux, Saux, nuotot, nuo, nspin,
     .              ispin, ncounter, delt,nstp,ElecDynPC)
C ....................

C New Density Matrix this again is not prepared for parallel operation
        nocc=wavef%nocck(1,ispin)

C	allocate bix & bix2. it's inside spin-loop, so deallocate before next iteration.
	allocate(bix(2*nocc,nuotot))
	allocate(bix2(nuotot,nocc))

         do ie = 1,nocc
          iie=wavef%iwf(ie,1,ispin)
          do j = 1,nuotot
             bix(2*ie-1,j) = wavef%phi(1,j,iie)
             bix(2*ie,j) = wavef%phi(2,j,iie)
          enddo
         enddo

	  if(calculateEnew) then
          call applyinverSH(Saux,Haux,nuotot,nuo,nocc,bix,bix2)
          endif

         do io = 1,nuo
            do j = 1,numd(io)
              ind = listdptr(io) + j
              jo = listd(ind)

	do asn=1,nocc

              pipj = bix(2*asn-1,io)*bix(2*asn-1,jo)+
     .      bix(2*asn,io)*bix(2*asn,jo)
              pipj = pipj * dble(3 - nspin)
              Dnew(ind,ispin) = Dnew(ind,ispin) +  pipj

              if(calculateEnew) then 
                pipj = bix(2*asn-1,io)*real(bix2(jo,asn))+
     .      bix(2*asn,io)*aimag(bix2(jo,asn))
                pipj = pipj * dble(3 - nspin)
                Enew(ind,ispin) = Enew(ind,ispin) + pipj
              endif 

            enddo

          enddo
        enddo
C ....................
C Plot instantaneous "eigenvalues" if required
       if (EigenvsTime.or.(finalstep.eq.itimestep)) then 
        nocc=wavef%nocck(1,ispin)
        do ie = 1,nocc
          iie=wavef%iwf(ie,1,ispin)
          do j = 1,nuotot
             aux(1,j) = wavef%phi(1,j,iie)
             aux(2,j) = wavef%phi(2,j,iie)
          enddo
          eigv=0.0d0
c         dnrm=0.0d0
          do io = 1,nuotot
            do jo = 1,nuotot
              pipj = aux(1,io) * aux(1,jo) + aux(2,io) * aux(2,jo)
              eigv=eigv+dreal(Haux(jo,io))*pipj
c             dnrm=dnrm+dreal(Saux(io,jo))*pipj
            enddo
          enddo
          eo(ie,ispin)=eigv
        enddo
      endif



C	deallocate bix & bix2 before the next iteration on spin index.
	deallocate(bix,bix2)
C	below is end do loop on spins.
      enddo 



CSave wavefunctions when required
CZEB. write when wrttdfiles=true
C	if(wrttdfiles) then
C        call  iowavef('write',wavef,nuotot,1,nspin)
C	end if





  
      end



       subroutine Uphi(H, S, phi, no, nocc, deltat)


C *************************************************************************
C Subroutine that calculates the new wavefunction, given the old 
C wavefunction by using the formula for the time evolution. Gamma-point 
C version. Written by A. Tsolakidis, May 2000 
C Modified by D. Sanchez-Portal, July 2002
C DSP 2008. This version is limited to first order expansion
C and avoids the inversion of the overlap
C *************************** INPUT ***************************************
C integer no                  : Number of basis orbitals
C integer nol                 : Local number of basis orbitals
C integer noccok              : Number of occupied wavefunctions per spin
C real*8 H(no,nol)             : Hamiltonian matrix
C real*8 S(no,nol)             : Overlap matrix
C complex*16  Phi(no,nocc) : Old wavefunctions
C *************************** OUTPUT *************************************
C complex*16  Phi(no,nocc) : New wavefunctions
C *************************** AUXILIARY ********************************** 
C complex*16 Q_1(no_max,no_max)           : Auxiliary Matrix
C complex*16 Q_2(no_max,no_max)           : Auxiliary Matrix
C complex*16 Q_3(no_max,no_max)           : Auxiliary Matrix
C real*8 deltat                           : Duration of the time step 
C*************************************************************************
C  Modules
C
      use precision
C*************************************************************************   
      implicit none 
     
      integer     no, ispin, nocc


      complex*16 H(no,no), S(no,no),  pi, pj, phi(no,nocc)
      
   
      external zgetrf, zgetri, zgemm

      double precision deltat

C  Internal variables .........................................


      integer i, j , k, info, no2, l
      complex*16,  dimension(:,:), allocatable, save :: aux1,aux2,aux3
      integer   ,  dimension(:), allocatable, save   :: ipiv
      complex*16  alpha, ss, hh

      logical  frsttime, onlyelectrons 
      

      data  frsttime / .true. /
      save  frsttime



      no2=no*no
      if (frsttime) then
C  Allocate memory for auxiliar storage
       allocate(aux1(no,no))
       call memory('A','Z',no2,'Uphi')
       allocate(aux2(no,no))
       call memory('A','Z',no2,'Uphi')
       allocate(aux3(no,no))
       call memory('A','Z',no2,'Uphi')
       allocate(ipiv(no))
       call memory('A','I',no,'Uphi')
       frsttime=.false.
      endif

C     first order expansion for the evolution operator
         alpha=-0.5d0*(0.d0,1.d0)*deltat
         do i=1,no
           do j=1,no
             hh=H(i,j)
             ss=S(i,j)
             aux1(i,j) = ss+alpha*hh
             aux2(i,j) = ss-alpha*hh
           enddo
         enddo 


C Invert Q_2 .........................................................

      call zgetrf(no,no,aux2,no,ipiv,info)
      call zgetri(no,aux2,no,ipiv,aux3,no2,info)
      if(info.ne.0) stop 'problem inverting aux2'

C Multiply Q_2 inverse  and Q_1 ...............................................
       call zgemm('n','n',no,no,no,(1.d0,0.d0),aux2,no,
     .           aux1,no,(0.d0,0.d0),aux3,no)

C Multiply Q_3 and psi ............................................
        call zgemm('n','n',no,nocc,no,(1.d0,0.d0),aux3,no,
     .                 phi,no,(0.d0,0.d0),aux2,no)
          do i=1,nocc
            do j=1,no
              phi(j,i)=aux2(j,i)
             enddo 
          enddo 
C ....................

      end


       subroutine evol1new(H, S, no, nol, nspin, ispin, 
     .       ncounter, delt, nstp, ElecDynPC)


C *************************************************************************
C Subroutine that calculates the new wavefunction, given the old 
C wavefunction by using the formula for the time evolution. Gamma-point 
C version. Written by A. Tsolakidis, May 2000 
C Modified by D. Sanchez-Portal, July 2002
C DSP 2008. This version is limited to first order expansion
C and avoids the inversion of the overlap
C*************************************************************************
C  Modules
C
      use precision
      use fdf
      use wavefunctions
C*************************************************************************   
      implicit none 
     
      integer     no,  ispin, ncounter, nol, nstp, nspin


      complex*16 H(no,nol), S(no,nol),  pi, pj
      

      double precision deltat, delt

      logical ElecDynPC

C  Internal variables .........................................


      integer i, j , k, info, no2, nocc, l
      complex*16,dimension(:,:,:),allocatable,save :: Hsave,Hextrpol
      complex*16,dimension(:,:),allocatable,save :: Hcorr
      complex*16  alpha,hh, sum
      
      logical fsttim(2)
      logical  frsttime, onlyelectrons, extrapol
C use linear extrapolation of the Hamiltonian 
      parameter(extrapol=.true.)

      data  fsttim(1), fsttim(2) /.true., .true. /
      data  frsttime / .true. /
      data onlyelectrons / .false. /
      save  frsttime, deltat, onlyelectrons, fsttim



      no2=no*no
      if (frsttime) then
C  Transform dt  to Ry**-1
C nstp is the number of "substeps" in the electronic evolution
C the evolution operator is applied in each substep although
C an extrapolated Hamiltonian is used "rather" than 
C a SCF Hamiltonian
       deltat=delt/0.04837d0/dble(nstp)
       write(6,*) 'evol1: time step (Ry**-1) ',deltat

       onlyelectrons=fdf_boolean('MD.OnlyElectrons',.false.)

C  Allocate memory for auxiliar storage
       allocate(Hsave(no,no,nspin))
       call memory('A','Z',no2*nspin,'evol1')
       if(ElecDynPC) then
         allocate(Hextrpol(no,no,nspin))
         call memory('A','Z',no2*nspin,'evol1')
       endif
       allocate(Hcorr(no,no))
       call memory('A', 'Z', no2,'evol1')
       frsttime=.false.
      endif
      nocc=wavef%nocck(1,ispin)

C Correction step
      if(.not.fsttim(ispin).and.ElecDynPC) then
      call timer('evol1.corr',1)
C Propagation backwards with the Hamiltonian extrapolated in the previous step
      do l=1,nstp
       do i=1,no
          do j=1,no
             hh=Hextrpol(i,j,ispin) 
             hh=hh+(l-0.5d0)*(Hsave(i,j,ispin)-hh)/dble(nstp)
             Hcorr(i,j)=-1.0d0*hh
           enddo
       enddo 
       call Uphi(Hcorr,S, wavef%phi(1,1,ncounter+1),
     .   no, nocc, deltat)
      enddo 
C Forward propagation with newly calculated Hamiltonian
      do l=1,nstp
       do i=1,no
          do j=1,no
             hh=H(i,j)
             hh=hh+(l-0.5d0)*(hh-Hsave(i,j,ispin))/dble(nstp)
             Hcorr(i,j)=hh
           enddo
       enddo
       call Uphi(Hcorr,S, wavef%phi(1,1,ncounter+1),
     .   no, nocc, deltat)
      enddo
      call timer('evol1.corr',2)
      endif

C Propagatio step into future times (using Hamilt. extrapolation)
      call timer('evol1.xtpl',1)
      do l=1,nstp
      if(fsttim(ispin).or..not.extrapol) then 
        call Uphi(H, S, wavef%phi(1,1,ncounter+1), 
     .   no, nocc, deltat)
      else
        do i=1,no
           do j=1,no
             hh=H(i,j)
             Hcorr(i,j)=hh+(l-0.5d0)*(hh-Hsave(i,j,ispin))/dble(nstp)
           enddo
        enddo 
         call Uphi(Hcorr, S, wavef%phi(1,1,ncounter+1),
     .   no, nocc, deltat)
      endif
      enddo 
      fsttim(ispin)=.false.

C Checking overlap with initial wavefunctions
c       do i=1,nocc
c           sum=(0.0d0,0.0d0)
c           do j=1,no
c             do k=1,no
c       pi=cmplx(wavef%phi(1,j,ncounter+i), wavef%phi(2,j,ncounter+i))
c       pj=cmplx( init%phi(1,k,ncounter+i), -init%phi(2,k,ncounter+i))  
c              sum=sum+pi*S(k,j)*pj
c            enddo 
c          enddo 
c          write(44+ncounter+i,*) dreal(sum), dimag(sum)
c       enddo 
C....................................................

C Storing Hamitonian for extrapolation and later correction    
        do i=1,no
           do j=1,no
              hh=H(i,j)
              if(ElecDynPC) 
     .          Hextrpol(i,j,ispin)=2.0d0*hh-Hsave(i,j,ispin)
              Hsave(i,j,ispin)=hh
           enddo 
        enddo 
        call timer('evol1.xtpl',2)
C ....................

      end

C	maz modified applyinverSH. 13 March 2011.
      subroutine applyinverSH(S,H,no,nol,nocc,psi,psi2)

      use precision

      complex*16 S(no,nol), H(no,nol)
      complex*16 psi(nocc,no), psi2(no,nocc)
      integer i, j , k, info, no2
      complex*16,  dimension(:,:), allocatable, save :: S_1
      complex*16,  dimension(:,:), allocatable, save :: aux1
      integer   ,  dimension(:), allocatable, save   :: ipiv
      logical  frstime
      data frstime /.true./
      save frstime


      no2=no*no
      if (frstime) then
C  Allocate memory for auxiliar storage
       allocate(aux1(no,no))
       call memory('A','Z',no2,'applyinverSH')
       allocate(ipiv(no))
       call memory('A','I',no,'applyinverSH')
C  Allocate memory to store the inverse matrix ................. 
       allocate(S_1(no,no))
       call memory('A','Z',no2,'applyinverSH')
       frstime=.false.
      endif

C  Invert the overlap matrix ...................................
      do i = 1, no
       do j = 1, nol
        S_1(i,j) = S(i,j)
       enddo
      enddo

      call zgetrf(no,no,S_1,no,ipiv,info)
      if(info.gt.0) write(6,*) 'info1 =',info 

      call zgetri(no,S_1,no,ipiv,aux1,no2,info)

      if(info.gt.0) then
          write(6,*) 'applyinverSH: info= ',info
          stop
      endif

      call zhemm('r','u',no,no,(1.0d0,0.0d0),H,no,S_1,no,
     .            (0.d0,0.d0),aux1,no)

      call zgemm('n','t',no,nocc,no,(1.0d0,0.0d0),aux1,no,
     .             psi,nocc,(0.d0,0.d0),psi2,no)


       end
