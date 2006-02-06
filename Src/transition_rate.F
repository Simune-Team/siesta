      subroutine transition_rate( ng, psi, ek, temp, efermi, Sr, 
     .                           Aux, Aux2, numh, listhptr, listh, 
     .                           indxuo, no, nuo, nuotot, xij, 
     .                           maxnh, nbands, kpoint, matrix )
C *********************************************************************
C Finds the matrix element for the dipolar transition between to states
C Written by DSP. August 1999
C Restyled for f90 version by JDG. June 2004
C **************************** INPUT ********************************** 
C integer ng                  : first dimension of psi, Aux and Aux2
C real*8  psi(ng,nuotot,nuo)  : Wavefunctions in current k point
C real*8  ek(nuotot)          : Eigenvalues
C real*8  temp                : Electronic temperature
C real*8  efermi              : Fermi level 
C real*8  Sr(maxnh)           : Position operator matrix elements (sparse)
C real*8  Aux(ng,nuotot,nuo)  : Auxiliary space
C real*8  Aux2(ng,nuotot,nuo)  : Auxiliary space
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
      use sys,          only : die
#ifdef MPI
      use mpi_siesta
#endif

      implicit none
 
C Passed variables
      integer nuo, nuotot, maxnh, no,
     .  listh(maxnh), numh(nuo), listhptr(nuo),
     .  indxuo(no), nbands, ng

      real(dp)
     .  xij(3,maxnh), Sr(maxnh), 
     .  psi(ng,nuotot,nuo), kpoint(3), Aux(ng,nuotot,nuo), 
     .  Aux2(ng,nuotot,nuo), ek(nuo), temp, efermi 

      character
     .  matrix*1

C Internal variables
      integer 
     .  ind, iuo, juo, j, ie, iie, iio, jje, je, jo, 
     .  BNodei, Bnodej, BTest, ig
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

C Check input matrix
      if(matrix.ne.'P'.and.matrix.ne.'R')
     $  call die('transrate: matrix only can take values R or P')

C Initialise matrix elements to zero
      do iuo = 1,nuo
        do juo = 1,nuotot
         do ig=1,ng
          Aux(ig,juo,iuo) = 0.0d0
          Aux2(ig,juo,iuo) = 0.0d0
         enddo 
        enddo 
      enddo 
          
C Compute matrix elements
      do iuo = 1,nuo
        do j = 1,numh(iuo)
          ind = listhptr(iuo) + j
          jo = listh(ind)
          juo = indxuo(jo)
          if(ng.gt.1) then 
           kxij = kpoint(1) * xij(1,ind) +
     .            kpoint(2) * xij(2,ind) +
     .            kpoint(3) * xij(3,ind) 
           ckxij = dcos(kxij)
           skxij = dsin(kxij)
           Aux2(1,juo,iuo) = Aux2(1,juo,iuo) + Sr(ind)*ckxij
           Aux2(2,juo,iuo) = Aux2(2,juo,iuo) + Sr(ind)*skxij
          else
           Aux2(1,juo,iuo) = Aux2(1,juo,iuo) + Sr(ind)
          endif
        enddo 
      enddo 

C Allocate workspace array
#ifdef MPI
      allocate(AuxLocal2(ng,nuotot))
#endif
      allocate(AuxLocal(ng,nuotot))
      allocate(psibandi(ng,nuotot))
      allocate(psibandj(ng,nuotot))
 
      BNodei = 0
      iie = 0
      do ie = 1,nbands
        f1 = 2.0d0*stepf((ek(ie)-efermi)/temp)
        if (Node.eq.BNodei) then
          iie = iie + 1 
          do j = 1,nuotot 
            do ig=1,ng
              psibandi(ig,j) = psi(ig,j,iie)
            enddo 
          enddo
        endif
#ifdef MPI
        call MPI_Bcast(psibandi(1,1),ng*nuotot,MPI_double_precision,
     .    BNodei,MPI_Comm_World,MPIerror)
#endif

        AuxLocal(1:ng,1:nuotot) = 0.0d0

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
                  do ig=1,ng
                    psibandj(ig,j) = psi(ig,j,jje)
                  enddo 
                enddo
              endif
#ifdef MPI
              call MPI_Bcast(psibandj(1,1),ng*nuotot,
     $        MPI_double_precision,BNodej,MPI_Comm_World,MPIerror)
#endif

              do iuo = 1,nuo
                call LocalToGlobalOrb(iuo,Node,Nodes,iio)
                do juo = 1,nuotot
                 if(ng.eq.2) then  
                   pipj1 = psibandi(1,iio)*psibandj(1,juo) +
     .                     psibandi(2,iio)*psibandj(2,juo)
                   pipj2 = psibandi(1,iio)*psibandj(2,juo) -
     .                     psibandi(2,iio)*psibandj(1,juo) 

                   AuxLocal(1,je) = AuxLocal(1,je)
     .                + pipj1*Aux2(1,juo,iuo)
     .                - pipj2*Aux2(2,juo,iuo)

                   AuxLocal(2,je) = AuxLocal(2,je)
     .                + pipj1*Aux2(2,juo,iuo)
     .                + pipj2*Aux2(1,juo,iuo)
                  else 
                   pipj1 = psibandi(1,iio)*psibandj(1,juo) 
                      
                   AuxLocal(1,je) = AuxLocal(1,je)
     .                + pipj1*Aux2(1,juo,iuo)
                  endif

                enddo 
              enddo  
              if (matrix.eq.'P') then 
                do ig=1,ng
                 AuxLocal(ig,je) = AuxLocal(ig,je)/(ek(je)-ek(ie))
                enddo 
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
        call MPI_Reduce(AuxLocal(1,1),AuxLocal2(1,1),ng*nuotot,
     .    MPI_double_precision,MPI_sum,BNodei,MPI_Comm_World,MPIerror)
        if (Node.eq.BNodei) then
          Aux(1:ng,1:nuotot,iie) = AuxLocal2(1:ng,1:nuotot)
        endif
#else
        if (Node.eq.BNodei) then
          Aux(1:ng,1:nuotot,iie) = AuxLocal(1:ng,1:nuotot)
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
