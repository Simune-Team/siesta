      subroutine diagpol( ispin, nspin, nuo, no, nuotot,
     .                    maxnh, numh, listhptr, listh, H, S,
     .                    xij, indxuo, kpoint, eo, psi,
     .                    Haux, Saux )
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, 
C for given Hamiltonian and Overlap matrices (including
C spin polarization), and for a given k_point
C Written by DSP, March 1999. From the routine diagk of J.Soler.
C Modified for parallel execution by J.D.Gale, March 2000.
C **************************** INPUT **********************************
C integer ispin               : Spin component which will be calculated
C integer nspin               : Number of spin components (1 or 2)
C integer nuo                 : Number of basis orbitals in unit cell
C integer no                  : Number of basis orbitals in supercell
C integer nuotot              : First dimension of eo, qo, last of xij
C integer maxnh               : Maximum number of orbitals interacting  
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to the start of each row
C                               of hamiltonian matrix
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
C real*8  kpoint(3)           : k point vectors
C real*8  psi                 : Workspace array
C real*8  Haux(2,nuotot,nuo)  : Workspace for dense H
C real*8  Saux(2,nuotot,nuo)  : Workspace for dense S
C *************************** OUTPUT **********************************
C real*8 eo(nuotot)           : Eigenvalues
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C eo H.
C *********************************************************************

      use precision
      use sys

      implicit          none

      integer
     .  maxnh, nuotot, no, nspin, nuo, indxuo(no), listh(maxnh), 
     .  listhptr(nuo), numh(nuo)

      real(dp)
     .  eo(nuotot), H(maxnh,nspin), kpoint(3), S(maxnh), 
     .  xij(3,maxnh), psi(2,nuotot,nuo), Haux(2,nuotot,nuo),
     .  Saux(2,nuotot,nuo)

      external          cdiag, memory

C  Internal variables .............................................
      integer
     .  ierror, ind, ispin, iuo, j, jo, juo
      real(dp)
     .  ckxij, kxij, skxij
      real(dp), dimension(:), allocatable, save :: aux
C  ....................

C Allocate local memory
      allocate(aux(2*nuotot*5))
      call memory('A','D',2*nuotot*5,'diagpol')

C Solve eigenvalue problem .........................................
      Saux = 0.0d0
      Haux = 0.0d0
      do iuo = 1,nuo
        do j = 1,numh(iuo)
          ind = listhptr(iuo) + j
          jo = listh(ind)
          juo = indxuo(jo)
          kxij = kpoint(1) * xij(1,ind) +
     .           kpoint(2) * xij(2,ind) +
     .           kpoint(3) * xij(3,ind)
          ckxij = cos(kxij)
          skxij = sin(kxij)
          Saux(1,juo,iuo) = Saux(1,juo,iuo) + S(ind)*ckxij
          Saux(2,juo,iuo) = Saux(2,juo,iuo) - S(ind)*skxij
          Haux(1,juo,iuo) = Haux(1,juo,iuo) + H(ind,ispin)*ckxij
          Haux(2,juo,iuo) = Haux(2,juo,iuo) - H(ind,ispin)*skxij
        enddo
      enddo
      call cdiag( Haux, nuotot, Saux, nuotot, nuo,
     .            eo, psi, nuotot, nuo, ierror )

C Deallocate local memory
      call memory('D','D',size(aux),'diagpol')
      deallocate(aux)

C Trap error flag from cdiag
      if (ierror.gt.0) then
        call die('Terminating due to failed diagonalisation')
      endif

      end
