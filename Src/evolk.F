

      subroutine evolk( nspin, maxspn, nuo,no,maxo,maxnh,maxnd, 

     .                  numh, listhptr, listh, numd, 

     .                  listdptr, listd, H, S, eo,

     .                  xij, indxuo, nk, kpoint, wk,

     .                  Dnew, Enew, 

     .                  Haux, Saux, Dk, Ek, nuotot,

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

C integer maxspn              : Maximum number of spin components

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

C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)

C                               (not used if only gamma point)

C integer indxuo(no)          : Index of equivalent orbital in unit cell

C                               Unit cell orbitals must be the first in

C                               orbital lists, i.e. indxuo.le.nuo, with

C                               nuo the number of orbitals in unit cell

C integer nk                  : Number of k points

C real*8  kpoint(3,nk)        : k point vectors

C real*8  wk(nk)              : k point weights (must sum one)

C integer nuotot              : total number of orbitals per unit cell

C                               over all processors

C real*8  delt                : length of the time step

C ******************** OUTPUT **************************************

C real*8 Dnew(maxnd,nspin)    : Output New Density Matrix

C real*8 Enew(maxnd,nspin)    : Output New Energy-Density Matrix

C real*8 eo(maxo,maxspn,nk)   : Output instantaneous eigenvalues

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

      use mpi

#endif

      use elec_dyn_options





      implicit none



#ifdef MPI

      integer 

     .  MPIerror

#endif

      

      integer

     .  maxnd, maxnh, nuo, no, nspin, nuotot, 

     .  ncounter, nk, maxspn, maxo



      integer 

     .  listh(maxnh), numh(nuo), listhptr(nuo),

     .  listd(maxnd), numd(nuo), listdptr(nuo),

     .  indxuo(no) 



      double precision

     .  Dnew(maxnd,nspin), Enew(maxnd,nspin),

     .  H(maxnh,nspin),

     .  S(maxnh), delt, wk(nk),

     .  xij(3,maxnh), kpoint(3,nk), eo(maxo,maxspn,nk)

     

	complex*16

     .  Haux(nuotot,nuo), Saux(nuotot,nuo), 

     .  Dk(nuotot,nuo), Ek(nuotot,nuo)

     

	double precision 

     .  aux(2,nuotot), aux2(2,nuotot)


	double precision, allocatable ::

     .  bix(:,:)
	
	complex*16, allocatable :: 

     .  bix2(:,:)





C  Internal variables .............................................

      integer           ie, io,  ispin, j, jo, BNode, iie, ind,

     .                  BTest, ierror, nd, nocc, ik, iuo, juo, nstp

      double precision

     .  ckxij, ee, kxij, qe, skxij, t, eigv, pipj1, pipj2

C MAZ		
	integer asn

C  ....................





C Number of "substeps" for electron dynamics (in between 

C Hamiltonian reevaluations)
       
       nstp=1
       !nstp=Nelecsubsteps



C New density and energy-density matrices of unit-cell orbitals .......

      nd = listdptr(nuo) + numd(nuo)

      Dnew(1:nd,1:nspin) = 0.d0

      if(calculateEnew) Enew(1:nd,1:nspin) = 0.d0



C Evolve wavefunctions.............................................

C Attention!! At the moment this is only prepared for serial work 

      ncounter=0

      do ik = 1,nk

       do ispin = 1,nspin

        nocc=wavef%nocck(ik,ispin)

C	allocate bix & bix2.
	allocate(bix(2*nocc,nuotot))
	allocate(bix2(nuotot,nocc))


   10   do io = 1,nuo

          do jo = 1,nuotot

            Saux(jo,io) = (0.0d0,0.0d0)

            Haux(jo,io) = (0.0d0,0.0d0)

          enddo

        enddo

        do iuo = 1,nuo

          do j = 1,numh(iuo)

            ind = listhptr(iuo) + j

            jo = listh(ind)

            juo = indxuo(jo)

              kxij = kpoint(1,ik) * xij(1,ind) +

     .               kpoint(2,ik) * xij(2,ind) +

     .               kpoint(3,ik) * xij(3,ind)

              ckxij = cos(kxij)

              skxij = -sin(kxij)

            Saux(juo,iuo) = Saux(juo,iuo) + 

     .                   cmplx(S(ind)*ckxij,S(ind)*skxij)

            Haux(juo,iuo) = Haux(juo,iuo) + 

     .       cmplx(H(ind,ispin)*ckxij,H(ind,ispin)*skxij)

          enddo

        enddo
   

	call evol2new( Haux, Saux, nuotot, nuo, nspin, nk,

     .              ispin, ik, ncounter, delt,

     .              nstp, ElecDynPC )


C Calculate "eigenvalue" if required

         if (EigenvsTime.or.(finalstep.eq.itimestep)) then

          nocc=wavef%nocck(ik,ispin)


          do ie=1,nocc

           iie=wavef%iwf(ie,1,ispin)

           do j = 1,nuotot

             aux(1,j) = wavef%phi(1,j,iie)

             aux(2,j) = wavef%phi(2,j,iie)

           enddo


           eigv=0.0d0

           do io = 1,nuotot

             do jo = 1,nuotot

               pipj1 = aux(1,io) * aux(1,jo) + aux(2,io) * aux(2,jo)

               eigv=eigv+dreal(Haux(jo,io))*pipj1

               pipj2 = aux(1,io) * aux(2,jo) - aux(2,io) * aux(1,jo)

               eigv=eigv+dimag(Haux(jo,io))*pipj2

             enddo

           enddo

           eo(ie,ispin,ik)=eigv

         enddo



        endif


        ncounter=ncounter+nocc

C.....................



C New Density Matrix this again is not prepared for parallel operation

C    Add contribution to density matrices of unit-cell orbitals


        qe=2.0d0*wk(ik)/dble(nspin)

        do iuo = 1,nuo

           do juo = 1,nuotot

             Dk(juo,iuo) = (0.d0,0.0d0)

             Ek(juo,iuo) = (0.d0,0.0d0)

          enddo

        enddo





C loop starts here.


        do ie = 1,nocc

          iie=wavef%iwf(ie,ik,ispin)

          do j = 1,nuotot

             bix(2*ie-1,j) = wavef%phi(1,j,iie)

             bix(2*ie,j) = wavef%phi(2,j,iie)

          enddo

         enddo



	  if(calculateEnew) then

          call applyinverSH(Saux,Haux,nuotot,nuo,nocc,bix,bix2)

          endif

          do iuo = 1,nuo
           do juo = 1,nuotot
		do asn=1,nocc

		pipj1 = bix(2*asn-1,iuo)*bix(2*asn-1,juo)+
     .      bix(2*asn,iuo)*bix(2*asn,juo)

		pipj2 = bix(2*asn-1,iuo)*bix(2*asn,juo)-
     .      bix(2*asn,iuo)*bix(2*asn-1,juo)

		Dk(juo,iuo) = Dk(juo,iuo) + qe*cmplx(pipj1,pipj2)

		if(calculateEnew) then
		
		pipj1 = bix(2*asn-1,iuo)*real(bix2(juo,asn))+
     .      bix(2*asn,iuo)*aimag(bix2(juo,asn))

		pipj2 = bix(2*asn-1,iuo)*aimag(bix2(juo,asn))-
     .      bix(2*asn,iuo)*real(bix2(juo,asn))

		Ek(juo,iuo) = Ek(juo,iuo) + qe*cmplx(pipj1,pipj2)
              endif

                enddo
           enddo
         enddo


          do iuo = 1,nuo

            do j = 1,numd(iuo)

              ind = listdptr(iuo) + j

              jo = listd(ind)

              juo = indxuo(jo)

              kxij = kpoint(1,ik) * xij(1,ind) +

     .               kpoint(2,ik) * xij(2,ind) +

     .               kpoint(3,ik) * xij(3,ind)

              ckxij = cos(kxij)

              skxij = -sin(kxij)

              Dnew(ind,ispin)=Dnew(ind,ispin)+ 

     .              dreal(Dk(juo,iuo))*ckxij +

     .              dimag(Dk(juo,iuo))*skxij

              if(calculateEnew)then

                  Enew(ind,ispin)=Enew(ind,ispin)+

     .              dreal(Ek(juo,iuo))*ckxij +

     .              dimag(Ek(juo,iuo))*skxij

              endif 

            enddo

          enddo


C MAZ	deallocate bix & bix2 before the next iteration of k/spins loops (end loops below).
	deallocate(bix,bix2)

       enddo

      enddo

C ....................


CSave wavefunctions when required
CZEB. write when wrttdfiles=true
	if(wrttdfiles) then
         call  iowavef('write',wavef,nuotot,nk,nspin)
	end if


      end





       subroutine evol2new(H, S, no, nol, nspin, nk, ispin, ik, 

     .       ncounter, delt,nstp,ElecDynPC)





C *************************************************************************

C Subroutine that calculates the new wavefunction, given the old 

C wavefunction by using the formula for the time evolution. Gamma-point 

C version. Written by A. Tsolakidis, May 2000 

C Modified by D. Sanchez-Portal, July 2002

C DSP 2008. This version is limited to first order expansion and

C avoids the inversion of the overlap

C*************************************************************************

C  Modules

C

      use precision

      use fdf

      use wavefunctions

C*************************************************************************   

      implicit none 

     

      integer     no, nstp, ispin, ncounter, nol, ik, nk, nspin



 



      complex*16 H(no,nol), S(no,nol), sum, pi, pj

      

      double precision deltat, delt



      logical ElecDynPC



C  Internal variables .........................................





      integer i, j , k, info, no2, nocc, l, no2k, jk

      complex*16,  dimension(:,:,:), allocatable, save :: Hsave,Hextrpol

      complex*16,  dimension(:,:), allocatable, save :: Hcorr

      logical   ,  dimension(:,:), allocatable, save   :: fsttimk

      complex*16  alpha, hh



      logical  frsttime, onlyelectrons, extrapol

      parameter(extrapol=.true.) 

      



      data  frsttime / .true. /

      data onlyelectrons /.false./

      save  frsttime, deltat, onlyelectrons







      no2=no*no

      no2k=no2*nk



      if (frsttime) then


C  Transform dt  to Ry**-1

C nstp is the number of "substeps" in the electronic evolution

C the evolution operator is applied in each substep although

C an extrapolated Hamiltonian is used "rather" than 

C a SCF Hamiltonian

       deltat=delt/0.04837d0/dble(nstp)
       write(6,*) 'evol2: n time step = ',nstp
       write(6,*) 'evol2: time step (Ry**-1) ',deltat



       onlyelectrons=fdf_boolean('MD.OnlyElectrons',.false.)



C  Allocate memory for auxiliar storage

C We store copies of the Hamiltonian in a clumsy way,

C this is OK if we have a very small number of k points

C and/or not too large unit cell



       allocate(Hsave(no,no,nk*nspin))

       call memory('A','Z',no2k*nspin,'evol2')

       if(ElecDynPC) then

         allocate(Hextrpol(no,no,nk*nspin))

         call memory('A','Z',no2k*nspin,'evol2')

       endif

       allocate(Hcorr(no,no))

       call memory('A', 'Z', no2,'evol2')

       allocate(fsttimk(nk,nspin))

       call memory('A','L',nk,'evol2')

       fsttimk(1:nk,1:nspin)=.true.

       frsttime=.false.


      endif

      nocc=wavef%nocck(ik,ispin)



C Correction step

      if(.not.fsttimk(ik,ispin).and.ElecDynPC) then

      call timer('evol2.corr',1)

C Propagation backwards with the Hamiltonian extrapolated in the previous step

      do l=1,nstp

       do i=1,no

          do j=1,no

             jk=ispin+nspin*(ik-1)

             hh=Hextrpol(i,j,jk)

             hh=hh+(l-0.5d0)*(Hsave(i,j,jk)-hh)/dble(nstp)

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

             jk=ispin+nspin*(ik-1)

             hh=H(i,j)

             hh=hh+(l-0.5d0)*(hh-Hsave(i,j,jk))/dble(nstp)

             Hcorr(i,j)=hh

           enddo

       enddo

       call Uphi(Hcorr,S, wavef%phi(1,1,ncounter+1),

     .   no, nocc, deltat)

      enddo

      call timer('evol2.corr',2)

      endif



C Propagation step into future times

      call timer('evol2.xtpl',1)

      do l=1,nstp

      if(fsttimk(ik,ispin).or..not.extrapol) then

        call Uphi(H, S, wavef%phi(1,1,ncounter+1),

     .   no, nocc, deltat)

      else

        do i=1,no

           do j=1,no

             jk=ispin+nspin*(ik-1)

             hh=H(i,j)

             Hcorr(i,j)=hh+(l-0.5d0)*(hh-Hsave(i,j,jk))/dble(nstp)

           enddo

        enddo

         call Uphi(Hcorr, S, wavef%phi(1,1,ncounter+1),

     .   no, nocc, deltat)

      endif

      enddo

      fsttimk(ik,ispin)=.false.



C Checking overlap with initial state

c       do i=1,nocc

c          sum=(0.0d0,0.0d0)

c          do j=1,no

c            do k=1,no

c       pi=cmplx(wavef%phi(1,j,ncounter+i), wavef%phi(2,j,ncounter+i))

c       pj=cmplx( init%phi(1,k,ncounter+i), -init%phi(2,k,ncounter+i))  

c              sum=sum+pi*S(k,j)*pj

c            enddo 

c          enddo 

c          write(44+ncounter+i,*) dreal(sum), dimag(sum)

c       enddo    

C Checking norm of the wavefunctions

c       do i=1,nocc

c          sum=(0.0d0,0.0d0)

c          do j=1,no

c            do k=1,no

c       pi=cmplx(wavef%phi(1,j,ncounter+i), wavef%phi(2,j,ncounter+i))

c       pj=cmplx(wavef%phi(1,k,ncounter+i), -wavef%phi(2,k,ncounter+i))

c              sum=sum+pi*S(k,j)*pj

c            enddo

c          enddo

c          write(60+i,*) dreal(sum), dimag(sum)

c       enddo

C.......................................................



C Storing Hamitonian for extrapolation and later correction    

        do i=1,no

           do j=1,no

              jk=ispin+nspin*(ik-1)

              hh=H(i,j)

              if(ElecDynPC) 

     .     Hextrpol(i,j,jk)=2.0d0*hh-Hsave(i,j,jk)

              Hsave(i,j,jk)=hh

           enddo

        enddo

        call timer('evol2.xtpl',2)

C ....................



      end



