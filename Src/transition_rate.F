      subroutine transition_rate(psi, ek, temp, efermi, Sr, 
     .                           Aux, Aux2, numh, listhptr, listh, 
     .                           indxuo, no, nuo, nuotot, xij, 
     .                           maxnh, nbands, kpoint, matrix )
C *********************************************************************
C Finds the matrix element for the dipolar transition between to states
C Written by DSP. August 1999
C Restyled for f90 version by JDG. June 2004
C **************************** INPUT ********************************** 
C real*8  psi(2,nuotot,nuo)   : Wavefunctions in current k point
C real*8  ek(nuotot)          : Eigenvalues
C real*8  temp                : Electronic temperature
C real*8  efermi              : Fermi level 
C real*8  Sr(maxnh)           : Position operator matrix elements (sparse)
C real*8  Sr2(maxnh)          : Position operator matrix elements (sparse)
C real*8  Sr3(maxnh)          : Position operator matrix elements (sparse)
C real*8  Aux(2,nuotot,nuo)   : Auxiliary space
C real*8  Aux2(2,nuotot,nuo)  : Auxiliary space
C integer numh(nuo)           : Number of nonzero elements of each row
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to start of row in listh
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element
C                               column indexes for each matrix row 
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C integer no                  : Number of basis orbitals in supercell
C integer nuo                 : Number of orbitals in the cell (locally)
C integer nuotot              : Number of orbitals in the cell (globally)
C integer maxnh               : Maximum dimension of listh
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)  
C integer nocc                : number of occupied states
C real*8  kpoint(3)           : Current kpoint
C real*8  dk(3)               : Vector joining the previous and current 
C                               kpoint
C character matrix*1          : 'R' or 'P' for position or momentum operator
C *************************** OUTPUT ****************************
C real*8  Aux(2,nuotot,nuo)   : matrix elements of the dipolar transition
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C *********************************************************************

      use precision
      use parallel,     only : BlockSize, Node, Nodes
      use parallelsubs, only : LocalToGlobalOrb
      use m_fermid,     only : stepf
#ifdef MPI
      use mpi_siesta
#endif

      implicit none
 
C Passed variables
      integer nuo, nuotot, maxnh, no,
     .  listh(maxnh), numh(nuo), listhptr(nuo),
     .  indxuo(no), nbands

      real(dp)
     .  xij(3,maxnh), Sr(maxnh), 
     .  psi(2,nuotot,nuo), kpoint(3), Aux(2,nuotot,nuo), 
     .  Aux2(2,nuotot,nuo), ek(nuo), temp, efermi 

      character
     .  matrix*1

C Internal variables
      integer 
     .  ind, iuo, juo, j, ie, iie, iio, jje, je, jo, 
     .  BNodei, Bnodej, BTest
      real(dp)
     .  kxij, skxij, ckxij, pipj1, pipj2,
     .  f1,  f2

#ifdef MPI
      integer                     :: MPIerror
      real(dp), allocatable, save :: AuxLocal2(:,:)
#endif
      real(dp), allocatable, save :: AuxLocal(:,:)
      real(dp), allocatable, save :: psibandi(:,:)
      real(dp), allocatable, save :: psibandj(:,:)
      real(dp),              save :: ediffmin = 1.0d-3
      real(dp),              save :: tiny = 1.0d-9


C Start timer
      call timer('transrate',1)

C Initialise matrix elements to zero
      do iuo = 1,nuo
        do juo = 1,nuotot
          Aux(1,juo,iuo) = 0.0d0
          Aux(2,juo,iuo) = 0.0d0
          Aux2(1,juo,iuo) = 0.0d0
          Aux2(2,juo,iuo) = 0.0d0
        enddo 
      enddo 
          
C Compute matrix elements
      do iuo = 1,nuo
        do j = 1,numh(iuo)
          ind = listhptr(iuo) + j
          jo = listh(ind)
          juo = indxuo(jo)
          kxij = kpoint(1) * xij(1,ind) +
     .           kpoint(2) * xij(2,ind) +
     .           kpoint(3) * xij(3,ind) 
          ckxij = dcos(kxij)
          skxij = dsin(kxij)
          Aux2(1,juo,iuo) = Aux2(1,juo,iuo) + Sr(ind)*ckxij
          Aux2(2,juo,iuo) = Aux2(2,juo,iuo) + Sr(ind)*skxij
        enddo 
      enddo 

C Allocate workspace array
#ifdef MPI
      allocate(AuxLocal2(2,nuotot))
#endif
      allocate(AuxLocal(2,nuotot))
      allocate(psibandi(2,nuotot))
      allocate(psibandj(2,nuotot))
 
      BNodei = 0
      iie = 0
      do ie = 1,nbands
        f1 = 2.0d0*stepf((ek(ie)-efermi)/temp)
        if (Node.eq.BNodei) then
          iie = iie + 1 
          do j = 1,nuotot 
            psibandi(1,j) = psi(1,j,iie)
            psibandi(2,j) = psi(2,j,iie)
          enddo
        endif
#ifdef MPI
#ifdef NODAT
        call MPI_Bcast(psibandi(1,1),2*nuotot,MPI_double_precision,
     .    BNodei,MPI_Comm_World,MPIerror)
#else       
        call MPI_Bcast(psibandi(1,1),2*nuotot,DAT_double,
     .    BNodei,MPI_Comm_World,MPIerror)
#endif          
#endif

        AuxLocal(1:2,1:nuotot) = 0.0d0

        BNodej = 0
        jje = 0
        do je = 1,nbands
          if (Node.eq.BNodej) then
            jje = jje + 1 
          endif
          if (dabs(ek(ie)-ek(je)).gt.ediffmin) then   
            f2 = 2.0d0*stepf((ek(je)-efermi)/temp)
                  
            if (f1*(2.0d0-f2).gt.tiny) then 
              if (Node.eq.BNodej) then
                do j = 1,nuotot 
                  psibandj(1,j) = psi(1,j,jje)
                  psibandj(2,j) = psi(2,j,jje)
                enddo
              endif
#ifdef MPI
#ifdef NODAT
              call MPI_Bcast(psibandj(1,1),2*nuotot,
     $        MPI_double_precision,BNodej,MPI_Comm_World,MPIerror)
#else       
              call MPI_Bcast(psibandj(1,1),2*nuotot,DAT_double,
     .          BNodej,MPI_Comm_World,MPIerror)
#endif          
#endif

              do iuo = 1,nuo
                call LocalToGlobalOrb(iuo,Node,Nodes,iio)
                do juo = 1,nuotot
                
                  pipj1 = psibandi(1,iio)*psibandj(1,juo) +
     .                    psibandi(2,iio)*psibandj(2,juo)
                  pipj2 = psibandi(1,iio)*psibandj(2,juo) -
     .                    psibandi(2,iio)*psibandj(1,juo) 

                  AuxLocal(1,je) = AuxLocal(1,je)
     .               + pipj1*Aux2(1,juo,iuo)
     .               - pipj2*Aux2(2,juo,iuo)

                  AuxLocal(2,je) = AuxLocal(2,je)
     .               + pipj1*Aux2(2,juo,iuo)
     .               + pipj2*Aux2(1,juo,iuo)

         
                enddo 
              enddo  
              if (matrix.eq.'P') then 
                AuxLocal(1,je) = AuxLocal(1,je)/(ek(je)-ek(ie))
                AuxLocal(2,je) = AuxLocal(2,je)/(ek(je)-ek(ie)) 
              endif 
            endif 
          endif 
          BTest = je/BlockSize
          if (BTest*BlockSize.eq.je) then
            BNodej = BNodej + 1
            if (BNodej .gt. Nodes-1) BNodej = 0
          endif
        enddo 
#ifdef MPI
#ifdef NODAT
        call MPI_Reduce(AuxLocal(1,1),AuxLocal2(1,1),2*nuotot,
     .    MPI_double_precision,MPI_sum,BNodei,MPI_Comm_World,MPIerror)
#else
        call MPI_Reduce(AuxLocal(1,1),AuxLocal2(1,1),2*nuotot,
     .    DAT_double,MPI_sum,BNodei,MPI_Comm_World,MPIerror)
#endif
        if (Node.eq.BNodei) then
          Aux(1:2,1:nuotot,iie) = AuxLocal2(1:2,1:nuotot)
        endif
#else
        if (Node.eq.BNodei) then
          Aux(1:2,1:nuotot,iie) = AuxLocal(1:2,1:nuotot)
        endif
#endif
        BTest = ie/BlockSize
        if (BTest*BlockSize.eq.ie) then
          BNodei = BNodei + 1
          if (BNodei .gt. Nodes-1) BNodei = 0
        endif
      enddo 

C Free workspace array
      deallocate(psibandj)
      deallocate(psibandi)
      deallocate(AuxLocal)
#ifdef MPI
      deallocate(AuxLocal2)
#endif
         
C Stop timer
      call timer('transrate',2)

      end
