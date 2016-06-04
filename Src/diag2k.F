      
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine diag2k( nuo, no, maxnh, maxnd, maxo,
     .                   numh, listhptr, listh, numd, listdptr, 
     .                   listd, H, S, getD, qtot, temp, e1, e2,
     .                   xij, indxuo, nk, kpoint, wk,
     .                   eo, qo, Dnew, Enew, ef, Entropy,
     .                   psi, nuotot, occtol, iscf)

!
!     Contributed by Volodymyr Maslyuk
!      
      use precision
      use sys
      use parallel,      only : Node, Nodes, BlockSize
      use parallelsubs,  only : LocalToGlobalOrb
      use m_fermid,      only : fermid, stepf
#ifdef MPI
      use mpi_siesta
#endif
      implicit           none
C *********************************************************************
C Calculates the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices.
C This version is for non-collinear spin with k-sampling and time 
C reversal symmetry.
C Written by J.Soler, August 1998.
C Modified by V.M.Garcia, June 2002.
C Density matrix build and computation of eigenvectors restricted
C by J.D. Gale November 2004.
C Corrected by V.V. Maslyuk for parallel case, Mai 2007
C **************************** INPUT **********************************
C integer nuo                 : Number of basis orbitals in unit cell
C integer no                  : Number of basis orbitals in supercell
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : First dimension of listd / DM
C integer maxo                : First dimension of eo and qo
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer numd(nuo)           : Number of nonzero elements of each row 
C                               of density matrix
C integer listdptr(nuo)       : Pointer to each row (-1) of the
C                               density matrix
C integer listd(maxnd)        : Nonzero density-matrix element column 
C                               indexes for each matrix row
C real*8  H(maxnh,4)          : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C logical getD                : Find occupations and density matrices?
C real*8  qtot                : Number of electrons in unit cell
C real*8  temp                : Electronic temperature 
C real*8  e1, e2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
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
C real*8  occtol              : Occupancy threshold for DM build
C integer iscf                : SCF cycle number
C *************************** OUTPUT **********************************
C real*8 eo(maxo*4,nk)        : Eigenvalues
C real*8 qo(maxo*4,nk)        : Occupations of eigenstates
C real*8 Dnew(maxnd,4)        : Output Density Matrix
C real*8 Enew(maxnd,4)        : Output Energy-Density Matrix
C real*8 ef                   : Fermi energy
C real*8 Entropy              : Electronic entropy
C *************************** AUXILIARY *******************************
C real*8 Haux(2,nuotot,2,nuo) : Aux. space for the hamiltonian matrix
C real*8 Saux(2,nuotot,2,nuo) : Aux. space for the overlap matrix
C real*8 psi(2,2,nuotot,2*nuo)  : Aux. space for the eigenvectors
C real*8 Dk(2,nuotot,2,nuo)   : Aux. space that may be the same as Haux
C real*8 Ek(2,nuotot,2,nuo)   : Aux. space that may be the same as Saux
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C temp and H must be in the same energy units.
C eo, Enew and ef returned in the units of H.
C *************************** PARALLEL ********************************
C The auxiliary arrays are now no longer symmetry and so the order
C of referencing has been changed in several places to reflect this.
C *********************************************************************
!
!     INPUT / OUTPUT
      
      integer maxuo, maxnd, maxnh, maxo, nk, no, nuo, nuotot, iscf
      
      integer indxuo(no), numh(nuo), numd(nuo)
      integer listh(maxnh), listd(maxnd)
      integer listhptr(*), listdptr(*)
      
      real(dp) Dnew(maxnd,4), Enew(maxnd,4), H(maxnh,4), S(maxnh)
      real(dp) kpoint(3,nk), wk(nk), xij(3,maxnh)
      real(dp) e1, e2, ef, eo(maxo*2,nk), qo(maxo*2,nk)
      real(dp) occtol, qtot, temp, Entropy
      
      logical               getD
      
!     TEMPOS, INTERNAL VARIABLES etc.
      
      real(dp) psi(2,2,nuotot,nuo*2)
      
      integer  BNode, BTest, ie, ierror, iie, ik, ind, io, iio
      integer  iuo, j, jo, juo, neigneeded
      real(dp) ee, qe, t
      
!     Haux(js,juo,is,iuo) = <js,juo|H|is,iuo>
!     Indices is and js are for spin components
!     Indices iuo and juo are for orbital components
      complex(dp),dimension(:,:,:,:), allocatable :: Haux, Saux
      
      complex(dp),dimension(:,:,:,:), allocatable :: Dkc, Ekc
      complex(dp)                                 :: cicj
      complex(dp)                                 :: D11, D22, D12, D21
      complex(dp),dimension(:,:),     allocatable :: caux
      complex(dp)                                 :: kphs
      real(dp)                                    :: kxij
#ifdef MPI
      integer            :: MPIerror
#endif

      external              cdiag, memory

#ifdef DEBUG
      call write_debug( '    PRE diag2k' )
#endif

!***********************************************************************
!     B E G I N
!***********************************************************************
!     define the arrays

      allocate(Haux(2,nuotot,2,nuo))
      call memory('A','Z',2*nuotot*2*nuo,'diag2k')
      allocate(Saux(2,nuotot,2,nuo))
      call memory('A','Z',2*nuotot*2*nuo,'diag2k')
      allocate(Dkc(2,nuotot,2,nuo))
      call memory('A','Z',2*nuotot*2*nuo,'diag2k')
      allocate(Ekc(2,nuotot,2,nuo))
      call memory('A','Z',2*nuotot*2*nuo,'diag2k')
      allocate(caux(2,nuotot))  ! CC RC Chnaged for osSO
      call memory('A','Z',2*nuotot,'diag2k')
      
      eo = 0.0_dp

!****************************************************************************
! Build Hamiltonian      
!****************************************************************************
! (Notes by Jaime Ferrer)
! Indices iuo and juo are for orbital components in the unit cell
! while jo is for orbital components in the supercell:
! d_{jo} = R + d_{juo}, x_{ij} = d_{jo} - d_{iuo}
! Since <jo|H|iuo> is computed, the Bloch phase factor is e^(-i k x_{ij}
! )
!
! The different subroutines build
! H_{jo,iuo)^{js,is} = H_{juo,iuo}^{js,is}(R) = <js,jo|H|is,iuo>
!
! The spin notation is as follows:
!
!               | H_{jo,iuo}^{u,u}  H_{jo,iuo}^{u,d} |
!  H_{jo,iuo} = |                                    |
!               | H_{jo,iuo}^{d,u}  H_{jo,iuo}^{d,d} |
!
!            | H(ind,1) + i H(ind,5)   H(ind,3) - i H(ind,4) |
!          = |                                               |
!            | H(ind,7) + i H(ind,8)   H(ind,2) + i H(ind,6) |
!
! 1. Hermiticity imposes H_{i,j}^{is,js}=H_{j,i}^{js,is}^*
! 2. Since wave functions are real, if there are no single P or L
! operators:
!      (a) H_{i,j}^{is,js}=H_{j,i}^{is,js}
!      (b) These imply spin-box hermiticity:
!      H_{i,j}^{is,js}=H_{i,j}^{js,is}^*
!
!      (c) Hence
!
!                               | H(ind,1)                H(ind,3) - i H(ind,4) |
!           H_{j,i} = H_{i,j} = |                                               |
!                               | H(ind,3) + i H(ind,4)   H(ind,2)              |
!
!

      
!     Find eigenvalues at every k point
      do ik = 1,nk
      
        Saux = dcmplx(0.0_dp,0.0_dp)
        Haux = dcmplx(0.0_dp,0.0_dp)
      
!
!       Transfer S,H matrices from sparse format in supercell to
!                                       full format in unit cell
!
!       Convention: ispin=1 => H11,       ispin=2 => H22,
!                   ispin=3 => Real(H12), ispin=4 => Imag(H12)
      
        do iuo = 1,nuo
          do j = 1,numh(iuo)
            ind = listhptr(iuo) + j
            jo = listh(ind)
            juo = indxuo(jo)
            kxij = kpoint(1,ik) * xij(1,ind) +
     .             kpoint(2,ik) * xij(2,ind) +
     .             kpoint(3,ik) * xij(3,ind)
! JFR: sign change in phase because of orbital ordering j,i & convention for xij = d_j - d_i
            kphs = cdexp(dcmplx(0.0_dp, -1.0_dp)*kxij)
      
            Saux(1,juo,1,iuo) = Saux(1,juo,1,iuo) + S(ind)   * kphs
            Saux(2,juo,2,iuo) = Saux(2,juo,2,iuo) + S(ind)   * kphs
            Haux(1,juo,1,iuo) = Haux(1,juo,1,iuo) + H(ind,1) * kphs
            Haux(2,juo,2,iuo) = Haux(2,juo,2,iuo) + H(ind,2) * kphs
            Haux(1,juo,2,iuo) = Haux(1,juo,2,iuo)
     .                        + dcmplx(H(ind,3), - H(ind,4)) * kphs
            Haux(2,juo,1,iuo) = Haux(2,juo,1,iuo)
     .                        + dcmplx(H(ind,3), + H(ind,4)) * kphs
          enddo
        enddo
      
!       Find eigenvalues
!       Note duplication since the spins are mixed
        call cdiag(Haux,Saux,2*nuotot,2*nuo,2*nuotot,eo(1,ik),psi,
     .             0,iscf,ierror, 2*BlockSize)

        if (ierror.ne.0) then
          call die('Terminating due to failed diagonalisation')
        endif
      enddo
      
!-----------------------------------------------------------------------
!     Check if we are done
!-----------------------------------------------------------------------
      if (.not.getD) goto 1001
!-----------------------------------------------------------------------
      
!     Find new Fermi energy and occupation weights
      call fermid(2,2,nk,wk,maxo,nuotot,eo,temp,qtot,qo,ef,Entropy)
      
!     Find weights for local density of states
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.d-6 )
        do ik = 1,nk
          do io = 1,nuotot*2
            qo(io,ik) = wk(ik) * ( stepf( (eo(io,ik)-e2)/t ) -
     .                             stepf( (eo(io,ik)-e1)/t ) )
          enddo
        enddo
      endif

!     New density and energy-density matrices of unit-cell orbitals
      Dnew(:,:) = 0.0_dp
      Enew(:,:) = 0.0_dp
      
!     Find maximum eigenvector that is required for this k point 
      do ik = 1,nk
        neigneeded = 0
        ie = 2*nuotot
        do while (ie.gt.0.and.neigneeded.eq.0)
          qe = qo(ie,ik)
          if (abs(qe).gt.occtol) neigneeded = ie
          ie = ie - 1
        enddo
      
!       Find eigenvectors 
        Saux = dcmplx(0.0_dp,0.0_dp)
        Haux = dcmplx(0.0_dp,0.0_dp)
        do iuo = 1,nuo
          do j = 1,numh(iuo)

            ind = listhptr(iuo) + j
            jo = listh(ind)
            juo = indxuo(jo)
            kxij = kpoint(1,ik) * xij(1,ind) +
     .             kpoint(2,ik) * xij(2,ind) +
     .             kpoint(3,ik) * xij(3,ind)
            kphs = cdexp(dcmplx(0.0_dp, -1.0_dp)*kxij)
      
            Saux(1,juo,1,iuo) = Saux(1,juo,1,iuo) + S(ind)   * kphs
            Saux(2,juo,2,iuo) = Saux(2,juo,2,iuo) + S(ind)   * kphs
            Haux(1,juo,1,iuo) = Haux(1,juo,1,iuo) + H(ind,1) * kphs
            Haux(2,juo,2,iuo) = Haux(2,juo,2,iuo) + H(ind,2) * kphs
            Haux(1,juo,2,iuo) = Haux(1,juo,2,iuo)
     .                        + dcmplx(H(ind,3), - H(ind,4)) * kphs
            Haux(2,juo,1,iuo) = Haux(2,juo,1,iuo)
     .                        + dcmplx(H(ind,3), + H(ind,4)) * kphs

          enddo
        enddo

        call cdiag(Haux,Saux,2*nuotot,2*nuo,2*nuotot,caux,psi,
     .             neigneeded,iscf,ierror, 2*BlockSize)
      
!       Check error flag and take appropriate action
        if (ierror.gt.0) then
          call die('Terminating due to failed diagonalisation')
        elseif (ierror.lt.0) then
!         Repeat diagonalisation with increased memory to handle clustering
          Saux = dcmplx(0.0_dp,0.0_dp)
          Haux = dcmplx(0.0_dp,0.0_dp)
          do iuo = 1,nuo
            do j = 1,numh(iuo)
              ind = listhptr(iuo) + j
              jo = listh(ind)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xij(1,ind) +
     .               kpoint(2,ik) * xij(2,ind) +
     .               kpoint(3,ik) * xij(3,ind)
              kphs = cdexp(dcmplx(0.0_dp, -1.0_dp)*kxij)
      
              Saux(1,juo,1,iuo) = Saux(1,juo,1,iuo) + S(ind)   * kphs
              Saux(2,juo,2,iuo) = Saux(2,juo,2,iuo) + S(ind)   * kphs
              Haux(1,juo,1,iuo) = Haux(1,juo,1,iuo) + H(ind,1) * kphs
              Haux(2,juo,2,iuo) = Haux(2,juo,2,iuo) + H(ind,2) * kphs
              Haux(1,juo,2,iuo) = Haux(1,juo,2,iuo)
     .                          + dcmplx(H(ind,3), - H(ind,4)) * kphs
              Haux(2,juo,1,iuo) = Haux(2,juo,1,iuo)
     .                          + dcmplx(H(ind,3), + H(ind,4)) * kphs
            enddo
          enddo
          call cdiag(Haux,Saux,2*nuotot,2*nuo,2*nuotot,caux,psi,
     .               neigneeded,iscf,ierror,2*BlockSize)
        endif
      
!-----------------------------------------------------------------------	     
!     Store the products of eigenvectors in matrices Dk and Ek
!     WARNING: Dk and Ek may be EQUIVALENCE'd to Haux and Saux
!-----------------------------------------------------------------------	     
!***********************************************************************
! BUILD NEW DENSITY MATRIX
!***********************************************************************
!
!                 | ------- 1,1 -------     ------- 1,2 ------- |
!                 | c_{i,up} c_{j,up}^*     c_{i,up} c_{j,down)^* |
!     D_{i,j} =   |                                             |
!                 | ------- 2,1 -------     ------- 2,2 ------- |
!                 | c_{i,down} c_{j,up}^*     c_{i,dn} c_{j,dn)^* |
!
!             =   | D_{i,j}(1)              D_{i,j}(3)-i D_{i,j}(4) |
!                 | D_{i,j}(7)+i D_{i,j}(8) D_{i,j}(2)              |

! The Energy is computed as E = Tr [ D_{i,j} H_{j,i} ]
!
! The Density Matrix is not "spin box hermitian" even if H does not
! contain P or
! L operators.
!
! Spin-box symmetrization of D inside this present subroutine:
! D_{i,j}(1,2) = 0.5 ( D_{i,j}(1,2) + D_{i,j}(2,1)^* )
! does not affect the results in any case, since H is spin-box
! hermitian.
!
! The wave function coefficients c_{i,sigma}(ie,k) only satisfy
! "TR"-related 
! symmetry relations if the state is collinear (nspin = 2). In this case
! c_{i,sigma}(ie,k) = c_{i,sigma}(ie, -k)^*
!
! For a non-collinear state, the above relation does not hold
! generically. Hence the
! full BZ must be used to compute the DM.

      
        Dkc  = dcmplx(0.0_dp,0.0_dp)
        Ekc  = dcmplx(0.0_dp,0.0_dp)
        caux = dcmplx(0.0_dp,0.0_dp)
      
        BNode = 0
        iie = 0
        do ie = 1,2*nuotot
          qe = qo(ie,ik)
          if (Node.eq.BNode) then
            iie = iie + 1
          endif
      
          caux(:,:) = dcmplx(0.0_dp,0.0_dp)
          if (qe.gt.occtol) then
            if (Node.eq.BNode) then
              do j = 1,nuotot
                caux(1,j)=dcmplx(psi(1,1,j,iie),psi(2,1,j,iie)) ! c_{i,up}
                caux(2,j)=dcmplx(psi(1,2,j,iie),psi(2,2,j,iie)) ! c_{i,dn}
              enddo
            endif
#ifdef MPI
            call MPI_Bcast(caux(1,1),2*nuotot,MPI_double_complex,BNode,
     .           MPI_Comm_World,MPIerror)
#endif
            ee = qo(ie,ik) * eo(ie,ik)
            do iuo = 1,nuo
              call LocalToGlobalOrb(iuo,Node,Nodes,iio)
              do juo = 1,nuotot
CC RC
C    In the JF 3.1.1 version the following matrix elements are
C    calculated by means of dconj(caux(1,juo)) instead of caux(1,juo)
C    I will leave it as it is. maybe it needs further revision...
      
!------- 1,1 -----------------------------------------------------------
               cicj = caux(1,iio) * dconjg(caux(1,juo))
               Dkc(1,juo,1,iuo) = Dkc(1,juo,1,iuo) + qe * cicj
               Ekc(1,juo,1,iuo) = Ekc(1,juo,1,iuo) + ee * cicj
!------- 2,2 -----------------------------------------------------------
               cicj = caux(2,iio) * dconjg(caux(2,juo))
               Dkc(2,juo,2,iuo) = Dkc(2,juo,2,iuo) + qe * cicj
               Ekc(2,juo,2,iuo) = Ekc(2,juo,2,iuo) + ee * cicj
!------- 2,1 -----------------------------------------------------------
               cicj = caux(1,iio) * dconjg(caux(2,juo))
               Dkc(1,juo,2,iuo) = Dkc(1,juo,2,iuo) + qe * cicj
               Ekc(1,juo,2,iuo) = Ekc(1,juo,2,iuo) + ee * cicj
!------- 1,2 -----------------------------------------------------------
               cicj = caux(2,iio) * dconjg(caux(1,juo))
               Dkc(2,juo,1,iuo) = Dkc(2,juo,1,iuo) + qe * cicj
               Ekc(2,juo,1,iuo) = Ekc(2,juo,1,iuo) + ee * cicj
              enddo
            enddo
          endif
          BTest = ie/(2*BlockSize)
          if (BTest*2*BlockSize.eq.ie) then
            BNode = BNode + 1
            if (BNode .gt. Nodes-1) BNode = 0
          endif
        enddo
      
!       Add contribution to density matrices of unit-cell orbitals
!       ------------------------------------------------------
!       | D_{j,i}(1)                 D_{j,i}(3)-i D_{j,i}(4) |
!       | D_{i,j}(3)+i D_{i,j}(4)    D_{j,i}(2)              |
!       ------------------------------------------------------
      
        do iuo = 1,nuo
        do j = 1,numd(iuo)
          ind = listdptr(iuo) + j
          jo = listd(ind)
          juo = indxuo(jo)
          kxij = kpoint(1,ik) * xij(1,ind) +
     .           kpoint(2,ik) * xij(2,ind) +
     .           kpoint(3,ik) * xij(3,ind)
! JFR: the sign in the phase of H_ji is reflected in the DM
          kphs = cdexp(dcmplx(0.0_dp,-1.0_dp)*kxij)
     
CC RC  
C      In the JF's version D11, D22, etc are not time-reversal averaged
C      I will leave as it is here but maybe it will need major revision
 
!         Average k and -k solutions because time-reversal symetry
!          D11 = 0.5_dp * (Dkc(1,juo,1,iuo) * kphs
!     .        +    dconjg(Dkc(1,juo,1,iuo) * kphs))
!      
!          D22 = 0.5_dp * (Dkc(2,juo,2,iuo) * kphs
!     .        +    dconjg(Dkc(2,juo,2,iuo) * kphs))
!          
!          D12 = 0.5_dp * (Dkc(1,juo,2,iuo) * kphs
!     .        +    dconjg(Dkc(2,juo,1,iuo) * kphs))
!      
!          D21 = 0.5_dp * (Dkc(2,juo,1,iuo) * kphs
!     .        +    dconjg(Dkc(1,juo,2,iuo) * kphs))

          D11 = Dkc(1,juo,1,iuo) * kphs
          D22 = Dkc(2,juo,2,iuo) * kphs
          D12 = Dkc(1,juo,2,iuo) * kphs
          D21 = Dkc(2,juo,1,iuo) * kphs

      
! Make DM Spin-box hermitian 
          D12 = 0.5_dp * (D12 + dconjg(D21))
          D21 = dconjg(D12)

          Dnew(ind,1) = Dnew(ind,1) + dreal(D11)
          Dnew(ind,2) = Dnew(ind,2) + dreal(D22)
          Dnew(ind,3) = Dnew(ind,3) + dreal(D12)
          Dnew(ind,4) = Dnew(ind,4) + dimag(D21)

CC RC  
C      In the JF's version D11, D22, etc are not time-reversal averaged
C      I will leave as it is here but maybe it will need major revision

!         Average k and -k solutions because time-reversal symetry
!          D11 = 0.5_dp * (Ekc(1,juo,1,iuo) * kphs
!     .        +    dconjg(Ekc(1,juo,1,iuo) * kphs))
!      
!          D22 = 0.5_dp * (Ekc(2,juo,2,iuo) * kphs
!     .        +    dconjg(Ekc(2,juo,2,iuo) * kphs))
!          
!          D12 = 0.5_dp * (Ekc(1,juo,2,iuo) * kphs
!     .        +    dconjg(Ekc(2,juo,1,iuo) * kphs))
!      
!          D21 = 0.5_dp * (Ekc(2,juo,1,iuo) * kphs
!     .        +    dconjg(Ekc(1,juo,2,iuo) * kphs))

          D11 = Ekc(1,juo,1,iuo) * kphs
          D22 = Ekc(2,juo,2,iuo) * kphs
          D12 = Ekc(1,juo,2,iuo) * kphs
          D21 = Ekc(2,juo,1,iuo) * kphs

     
! Make DM Spin-box hermitian  
          D12 = 0.5_dp * (D12 + dconjg(D21))
          D21 = dconjg(D12)

          Enew(ind,1) = Enew(ind,1) + dreal(D11)
          Enew(ind,2) = Enew(ind,2) + dreal(D22)
          Enew(ind,3) = Enew(ind,3) + dreal(D12) 
          Enew(ind,4) = Enew(ind,4) + dimag(D21) 

        enddo
        enddo
      enddo

 1001 continue
      
!     deallocate the arrays
      call memory('D','Z',2*nuotot*2*nuo,'diag2k')
      deallocate(Haux)
      call memory('D','Z',2*nuotot*2*nuo,'diag2k')
      deallocate(Saux)
      call memory('D','Z',2*nuotot*2*nuo,'diag2k')
      deallocate(Dkc)
      call memory('D','Z',2*nuotot*2*nuo,'diag2k')
      deallocate(Ekc)
      call memory('D','Z',2*2*nuotot,'diag2k')
      deallocate(caux)
      
!        write(160+node,*)'ITERATION'
!        write(170+node,*)'ITERATION'
!        do iuo = 1,nuo
!          do j = 1,numd(iuo)
!            ind = listdptr(iuo) + j
!            jo = listd(ind)
!            juo = indxuo(jo)
!             write(160+node,'(3i5,4f15.9)')iuo,juo,jo,Dnew(ind,1:4)
!             write(170+node,'(3i5,4f15.9)')iuo,juo,jo,
!     .    (Dnew(ind,1)+Dnew(ind,2))*0.5,
!     .    (Dnew(ind,1)-Dnew(ind,2))*0.5,Dnew(ind,3:4)
!          enddo
!        enddo
!      
!      call die('maslyuk')
!***********************************************************************
#ifdef DEBUG
      call write_debug( '    POS diag2k' )
#endif

      return
      end subroutine diag2k
!***********************************************************************      
      
