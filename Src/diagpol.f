      subroutine diagpol( ispin, nspin, nuo, no, 
     .                  maxuo, maxno, numh, listh, H, S,
     .                  xij, indxuo, kpoint, eo, psi)
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, 
C for given Hamiltonian and Overlap matrices (including
C spin polarization), and for a given k_point
C Written by DSP, March 1999. From the routine diagk de J.Soler.
C **************************** INPUT **********************************
C integer ispin               : Spin component which will be calculated
C integer nspin               : Number of spin components (1 or 2)
C integer nuo                 : Number of basis orbitals in unit cell
C integer no                  : Number of basis orbitals in supercell
C integer maxuo               : First dimension of eo, qo, last of xij
C                               Must be at least max(indxuo)
C integer maxno               : Maximum number of orbitals interacting  
C                               with any orbital
C integer numh(no)            : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listh(maxno,no)     : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C real*8  H(maxno,maxuo,nspin): Hamiltonian in sparse form
C real*8  S(maxno,maxuo)      : Overlap in sparse form
C real*8  xij(3,maxno,maxuo)  : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C real*8  kpoint(3)           : k point vectors
C real*8  psi                 : Workspace array
C *************************** OUTPUT **********************************
C real*8 eo(maxuo)     : Eigenvalues
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C eo H.
C *********************************************************************
      implicit          none
      integer           maxno, maxuo, no,
     .                  nspin, nuo
      integer           indxuo(no), listh(maxno,no), numh(no)

      double precision  
     .                  eo(maxuo), H(maxno,maxuo,nspin),
     .                  kpoint(3), 
     .                  S(maxno,maxuo), 
     .                  xij(3,maxno,maxuo),
     .                  psi(2,maxuo,nuo)
      external          cdiag, memory

C  Internal variables .............................................
      integer
     .  io, ispin, iuo, j, jo, juo, ierror
      double precision
     .  ckxij, kxij, skxij
      double precision, dimension(:,:,:), allocatable ::
     .  Haux, Saux
      double precision, dimension(:), allocatable ::
     .  aux
C  ....................

C Allocate local memory
      allocate(Haux(2,nuo,nuo))
      call memory('A','D',2*nuo*nuo,'diagpol')
      allocate(Saux(2,nuo,nuo))
      call memory('A','D',2*nuo*nuo,'diagpol')
      allocate(aux(2*nuo*5))
      call memory('A','D',2*nuo*5,'diagpol')

C Solve eigenvalue problem .........................................
      do iuo = 1,nuo
        do juo = 1,nuo
          Saux(1,juo,iuo) = 0.d0
          Saux(2,juo,iuo) = 0.d0
          Haux(1,juo,iuo) = 0.d0
          Haux(2,juo,iuo) = 0.d0
        enddo
      enddo
      do io = 1,nuo
        do j = 1,numh(io)
          jo = listh(j,io)
          iuo = indxuo(io)
          juo = indxuo(jo)
          kxij = kpoint(1) * xij(1,j,io) +
     .           kpoint(2) * xij(2,j,io) +
     .           kpoint(3) * xij(3,j,io)
          ckxij = cos(kxij)
          skxij = sin(kxij)
          Saux(1,iuo,juo)=Saux(1,iuo,juo)+S(j,io)*ckxij
          Saux(2,iuo,juo)=Saux(2,iuo,juo)+S(j,io)*skxij
          Haux(1,iuo,juo)=Haux(1,iuo,juo)+H(j,io,ispin)*ckxij
          Haux(2,iuo,juo)=Haux(2,iuo,juo)+H(j,io,ispin)*skxij
        enddo
      enddo
      do iuo = 1,nuo
        do juo = 1,iuo-1
          Saux(1,juo,iuo) = 0.5d0 * ( Saux(1,juo,iuo) +
     .                                Saux(1,iuo,juo) )
          Saux(1,iuo,juo) = Saux(1,juo,iuo)
          Saux(2,juo,iuo) = 0.5d0 * ( Saux(2,juo,iuo) -
     .                                Saux(2,iuo,juo) )
          Saux(2,iuo,juo) = -Saux(2,juo,iuo)
          Haux(1,juo,iuo) = 0.5d0 * ( Haux(1,juo,iuo) +
     .                                Haux(1,iuo,juo) )
          Haux(1,iuo,juo) = Haux(1,juo,iuo)
          Haux(2,juo,iuo) = 0.5d0 * ( Haux(2,juo,iuo) -
     .                                Haux(2,iuo,juo) )
          Haux(2,iuo,juo) = -Haux(2,juo,iuo)
        enddo
        Saux(2,iuo,iuo) = 0.d0
        Haux(2,iuo,iuo) = 0.d0
      enddo
      call cdiag( Haux, nuo, Saux, nuo, nuo,
     .            eo, psi, nuo, aux, ierror )

C Deallocate local memory
      call memory('D','D',size(Haux),'diagpol')
      deallocate(Haux)
      call memory('D','D',size(Saux),'diagpol')
      deallocate(Saux)
      call memory('D','D',size(aux),'diagpol')
      deallocate(aux)

C Trap error flag from cdiag
      if (ierror.gt.0) then
        call die('Terminating due to failed diagonalisation')
      endif

      end
