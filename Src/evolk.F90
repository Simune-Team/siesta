      subroutine evolk( nspin, maxspn, nuo,no,maxo,maxnh,maxnd, numh,     &
                      listhptr,listh, numd, listdptr, listd, H, S, eo,    &
                      xij, indxuo, nk, kpoint, wk, Dnew, Enew,      &
                      Dk, Ek, nuotot,delt,Haux,Saux,psi)



! *********************************************************************
! Subroutine to calculate the eigenvalues and eigenvectors, density
! and energy-density matrices, and occupation weights of each 
! eigenvector, for given Hamiltonian and Overlap matrices (including
! spin polarization). Gamma-point version.
! Written by A. Tsolakidis, May 2000 after a subroutine
! by J. M. Soler.
! rewritten by D. Sanchez-Portal, November 2002-March 2003
! Modified by M. Ahsan Zeb, 11 March 2011.
! **************************** INPUT **********************************
! integer nspin               : Number of spin components (1 or 2)
! integer maxspn              : Maximum number of spin components
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
! real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
!                               (not used if only gamma point)
! integer indxuo(no)          : Index of equivalent orbital in unit cell
!                               Unit cell orbitals must be the first in
!                               orbital lists, i.e. indxuo.le.nuo, with
!                               nuo the number of orbitals in unit cell
! integer nk                  : Number of k points
! real*8  kpoint(3,nk)        : k point vectors
! real*8  wk(nk)              : k point weights (must sum one)
! integer nuotot              : total number of orbitals per unit cell
!                               over all processors
! real*8  delt                : length of the time step
! ******************** OUTPUT **************************************
! real*8 Dnew(maxnd,nspin)    : Output New Density Matrix
! real*8 Enew(maxnd,nspin)    : Output New Energy-Density Matrix
! real*8 eo(maxo,maxspn,nk)   : Output instantaneous eigenvalues
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

      use precision
      use sys
      use parallel
      use wavefunctions
      use MatrixSwitch
      use siesta_options,    only: eigen_time
      use m_diagon_opt, only : ictxt
#ifdef MPI
      use mpi_siesta
#endif
      implicit none

      integer   :: maxnd, maxnh, nuo, no, nspin, nuotot,ncounter, nk, maxspn, maxo
      integer   :: listh(maxnh), numh(nuo), listhptr(nuo), listd(maxnd), numd(nuo), listdptr(nuo), indxuo(no) 
      REAL (kind=dp)     :: Dnew(maxnd,nspin), Enew(maxnd,nspin),            &
                            H(maxnh,nspin), S(maxnh), delt, wk(nk),          &
                            xij(3,maxnh), kpoint(3,nk), eo(maxo,maxspn,nk)
      complex(kind=dp)  :: Dk(nuotot,nuo), Ek(nuotot,nuo), varaux, varaux2, varaux3, varaux4, varaux5
      type(matrix)      ::  Hauxms,Sauxms, aux,aux2,bix,bix2
      character         :: m_operation*3, m_storage*5
! Internal variables .............................................
      complex(dp) :: Haux(nuotot,nuo),Saux(nuotot,nuo),psi(nuotot,nuo)
      integer   :: ie, io,  ispin, j, jo, BNode, iie, ind, BTest,            &
                   ierror, nd, nocc, ik, iuo, juo, nstp,desch(9)
      real(kind=dp)     :: ckxij, ee, kxij, qe, skxij, t, eigv, pipj1, pipj2
      logical           :: calculateEnew
!MAZ		
	  integer           :: asn
! ....................

!Number of "substeps" for electron dynamics (in between 
#ifdef MPI
 m_storage='pzdbc'
 m_operation='lap'
#else
 m_storage='szden'
 m_operation='lap'
#endif
 call m_allocate( Hauxms,nuotot,nuotot,m_storage)
 call m_allocate( Sauxms,nuotot,nuotot,m_storage)


#ifdef MPI
    call descinit(desch,nuotot,nuotot,BlockSize,BlockSize,0,0,ictxt,nuotot,ierror)
#endif
 

!Hamiltonian reevaluations)
nstp=1       
       !nstp=Nelecsubsteps

!New density and energy-density matrices of unit-cell orbitals .......

      nd = listdptr(nuo) + numd(nuo)
      Dnew(1:nd,1:nspin) = 0.d0

      if(calculateEnew) Enew(1:nd,1:nspin) = 0.d0

!Evolve wavefunctions.............................................

!Attention!! At the moment this is only prepared for serial work 

       call m_allocate(aux,1,nuotot,m_storage)

      ncounter=0
      do ik = 1,nk
       do ispin = 1,nspin
        !nocc=wavef%nocck(ik,ispin)
         nocc=wavef_ms(ik,ispin)%dim2
!allocate bix & bix2.
        call m_allocate(bix,nocc,nuotot,m_storage)
        call m_allocate(bix2,nuotot,nocc,m_storage)

        call m_set(Hauxms,'a',cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)
       call m_set(Sauxms,'a',cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)


      do io = 1,nuo
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
 				kxij = kpoint(1,ik) * xij(1,ind) +&
 				kpoint(2,ik) * xij(2,ind) +&
 				kpoint(3,ik) * xij(3,ind)
 				ckxij = cos(kxij)
 				skxij = -sin(kxij)

            Saux(juo,iuo) = Saux(juo,iuo) + cmplx(S(ind)*ckxij,S(ind)*skxij,dp)

            Haux(juo,iuo) = Haux(juo,iuo) + cmplx(H(ind,ispin)*ckxij,H(ind,ispin)*skxij,dp)

             !Haux(juo,iuo)=Haux(juo,iuo)+cmplx(H(ind,ispin)*ckxij,H(ind,ispin)*skxij,dp)
             !Saux(juo,iuo)=Saux(juo,iuo)+cmplx(S(ind)*ckxij,S(ind)*skxij,dp)
    enddo
 enddo





do io=1,nuotot
  do j=1,nuotot

#ifdef MPI

             call pzelget('a',' ',varaux,Haux,j,io,desch)
           
  !           call pzelget('a',' ',varaux2,Saux,j,io,desch)
#else
varaux=Haux(j,io)
varaux2=Saux(j,io)
#endif
  call m_set_element( Hauxms,j,io,varaux,m_operation)
 call m_set_element( Sauxms,j,io,varaux2,m_operation)
enddo
enddo

call m_get_element(wavef_ms(ik,ispin),1,1,varaux,m_operation)

     call evol2new( Hauxms, Sauxms, nuotot, nuo, nspin, nk, ispin, ik, ncounter, delt, nstp)
 
!
        ncounter=ncounter+nocc

!New Density Matrix this again is not prepared for parallel operation

!Add contribution to density matrices of unit-cell orbitals


        qe=2.0d0*wk(ik)/dble(nspin)

        do iuo = 1,nuo

           do juo = 1,nuotot

             Dk(juo,iuo) = (0.d0,0.0d0)

             Ek(juo,iuo) = (0.d0,0.0d0)

          enddo

        enddo

!loop starts here.

             call m_add(wavef_ms(ik,ispin),'c',bix,cmplx(1.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)
    
     


	  if(calculateEnew) then

          call applyinverSH(Sauxms,Hauxms,nuotot,nuo,nocc,bix,bix2)

          endif

do io=1,nuotot
  do j=1,nuotot 

#ifdef MPI
   call pzelset(psi,io,j,desch,cmplx(0.0_dp,0.0_dp,dp))
#else
              psi(io,j)=cmplx(0.0_dp,0.0_dp,dp)
#endif

end do
end do


pipj1= 0.0_dp
pipj2=0.0_dp

          do iuo = 1,nuotot
           do juo = 1,nuotot
               varaux3 = cmplx(0.0_dp,0.0_dp,dp)
               varaux5 = cmplx(0.0_dp,0.0_dp,dp)
		do asn=1,nocc

               call m_get_element(bix,asn,iuo,varaux,m_operation)
               call m_get_element(bix,asn,juo,varaux2,m_operation)
               
		pipj1 =real(varaux)*real(varaux2)+aimag(varaux)*aimag(varaux2) 
      
		pipj2 = real(varaux)*aimag(varaux2)-aimag(varaux)*real(varaux2)
              
                varaux3 = varaux3 + qe*cmplx(pipj1,pipj2)
               
               IF (calculateEnew) THEN
               call m_get_element(bix2,juo,asn,varaux4,m_operation)
               pipj1 = real(varaux)*real(varaux4) + aimag(varaux)*aimag(varaux4)
               pipj2 = real(varaux)*aimag(varaux4)+aimag(varaux)*real(varaux4)
               varaux5 = varaux5 + qe*cmplx(pipj1,pipj2)
               END IF
              
                enddo
#ifdef MPI
                call pzelset(Dk,juo,iuo,desch,varaux3)
                IF (calculateEnew) call pzelset (Ek,juo,iuo,desch,varaux5)
#else
		Dk(juo,iuo) = varaux3
                IF (calculateEnew) Ek(juo,iuo) = varaux5
#endif
           enddo
         enddo

! 		
do iuo=1,nuo
  do j = 1,numh(iuo)
 			ind = listhptr(iuo) + j
 			jo = listh(ind)
 			juo = indxuo(jo)
 				kxij = kpoint(1,ik) * xij(1,ind) +&
 				kpoint(2,ik) * xij(2,ind) +&
 				kpoint(3,ik) * xij(3,ind)
 				ckxij = cos(kxij)
 				skxij = -sin(kxij)
              Dnew(ind,ispin)=Dnew(ind,ispin)+ dreal(Dk(juo,iuo))*ckxij +&
                dimag(Dk(juo,iuo))*skxij

              if(calculateEnew)then

                 Enew(ind,ispin)=Enew(ind,ispin)+ dreal(Ek(juo,iuo))*ckxij     &
                                + dimag(Ek(juo,iuo))*skxij

             endif 

            enddo

          enddo


!MAZ	deallocate bix & bix2 before the next iteration of k/spins loops (end loops below).
 call m_deallocate(bix)
 call m_deallocate(bix2)
       enddo

      enddo
call m_deallocate(aux)
call m_deallocate(Hauxms)
!call m_deallocate(Sauxms)

!ave wavefunctions when required
!EB. write when wrttdfiles=true
!if(wrttdfiles) then
!        call  iowavef('write',wavef,nuotot,nk,nspin)
!end if


      END SUBROUTINE evolk





      SUBROUTINE  evol2new(Hauxms, Sauxms, no, nol, nspin, nk, ispin, ik, ncounter, delt,nstp)





!*************************************************************************

!Subroutine that calculates the new wavefunction, given the old 

!wavefunction by using the formula for the time evolution. Gamma-point 

!version. Written by A. Tsolakidis, May 2000 

!Modified by D. Sanchez-Portal, July 2002

!DSP 2008. This version is limited to first order expansion and

!avoids the inversion of the overlap

!************************************************************************

! Modules



      use  matswinversion
      use precision
      use fdf
      use wavefunctions
      use MatrixSwitch
!************************************************************************   

      implicit none 

      integer                                           :: no, nstp, ispin, ncounter, nol, ik, nk, nspin
      complex(kind=dp)                                  :: sum, pi, pj,varaux3,varaux2
      type(matrix)                                      :: Hauxms, Sauxms
      real(kind=dp)                                     :: delt, varaux

! Internal variables .........................................

      integer                                           :: i, j , k, info, no2, nocc, l, no2k, jk
      type(matrix),allocatable,save                     :: Hsave(:)
      logical,dimension(:,:), allocatable, save         :: fsttimk
      complex(kind=dp)                                  :: al, hh
      logical, parameter                                :: extrapol=.true.
      character(5)                                      :: m_storage
      character(3)                                      :: m_operation
      logical, save                                     :: frsttime = .true. 
      logical, save                                     :: onlyelectrons =.false.
      real(kind=dp), save                               :: deltat
      logical                                           :: calculateEnew
#ifdef MPI
 m_storage='pzdbc'
 m_operation='lap'
#else
 m_storage='szden'
 m_operation='lap'
#endif

      no2=no*no
      no2k=no2*nk

      if (frsttime) then


! Transform dt  to Ry**-1

!nstp is the number of "substeps" in the electroni!evolution

!the evolution operator is applied in each substep although

!an extrapolated Hamiltonian is used "rather" than 

!a S! Hamiltonian

       deltat=delt/0.04837d0/dble(nstp)
       
       write(6,*) 'evol2: n time step = ',nstp
       write(6,*) 'evol2: time step in (fm) ',delt
       write(6,*) 'evol2: time step (Ry**-1) ',deltat



       onlyelectrons=fdf_boolean('MD.OnlyElectrons',.false.)



! Allocate memory for auxiliar storage

!We store copies of the Hamiltonian in a clumsy way,

!this is OK if we have a very small number of k points

!and/or not too large unit cell

       allocate(Hsave(nk*nspin))
       do i=1,nk*nspin
           call m_allocate(Hsave(i),no,no,m_storage)
       end do
           call memory('A','Z',no2k*nspin,'evol2')



       call memory('A', 'Z', no2,'evol2')

       allocate(fsttimk(nk,nspin))

       call memory('A','L',nk,'evol2')

       fsttimk(1:nk,1:nspin)=.true.

       frsttime=.false.


      endif
      nocc=wavef_ms(ik,ispin)%dim2

      !nocc=wavef%nocck(ik,ispin)


!correction step
!Propagation step into future times

      call timer('evol2.xtpl',1)

      do l=1,nstp

      if(fsttimk(ik,ispin).or..not.extrapol) then
        !call Uphi(H, S, wavef%phi(1,1,ncounter+1), no, nocc, deltat)
        ! call Uphi(Hauxms, Sauxms,wavef_ms(ik,ispin), no, nocc, deltat)
      else
             jk=ispin+nspin*(ik-1)
        varaux=(l-0.5d0)/dble(nstp)
        call m_add(Hauxms,'n',Hsave(jk),cmplx(1.0,0.0,dp),cmplx(-1.0,0.0,dp),m_operation)
        call m_add(Hauxms,'n',Hsave(jk),cmplx(1.0,0.0,dp),cmplx(varaux,0.0,dp),m_operation) 
       ! call Uphi(Hsave(jk), Sauxms,wavef_ms(ik,ispin), no, nocc, deltat)
      endif
      enddo

      
      fsttimk(ik,ispin)=.false.





!Storing Hamitonian for extrapolation and later correction    



              jk=ispin+nspin*(ik-1)

       call m_add(Hauxms,'n',Hsave(jk),cmplx(1.0,0.0_dp,dp),cmplx(0.0,0.0,dp),m_operation)

        call timer('evol2.xtpl',2)

!....................


 END SUBROUTINE evol2new


