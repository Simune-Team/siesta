!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine diag3k( nuo, no, maxnh, maxnd, maxo,
     .                   numh, listhptr, listh, numd, listdptr, 
     .                   listd, H, S, getD, getPSI, qtot, temp, e1, e2,
     .                   xij, indxuo, nk, kpoint, wk,
     .                   eo, qo, Dnew, Enew, ef, Entropy,
     .                   Haux, Saux, psi, Dk, Ek, caux, 
     .                   nuotot, occtol, iscf, neigwanted)

!
      use precision
      use sys
      use parallel,      only : Node, Nodes, BlockSize
      use parallelsubs,  only : LocalToGlobalOrb
      use writewave,     only : writew
      use m_fermid,      only : fermid, stepf
      use iso_c_binding, only : c_f_pointer, c_loc
#ifdef MPI
      use mpi_siesta
#endif
      implicit           none
C *********************************************************************
C Calculates the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices.
C This version is for spin-orbit coupling with k-sampling and time 
C reversal symmetry.
C Written by J.Soler, August 1998.
C Modified by V.M.Garcia, June 2002.
C Density matrix build and computation of eigenvectors restricted
C by J.D. Gale November 2004.
C Reduced memory requirements, Nick, Aug. 2017
C **************************** INPUT **********************************
C integer nuo                 : Number of basis orbitals in unit cell
C integer no                  : Number of basis orbitals in supercell
C integer maxuo               : Maximum number of basis orbitals
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
C logical getPSI              : Find & print wavefunctions?
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
C integer neigwanted          : Number of eigenvalues wanted
C *************************** OUTPUT **********************************
C real*8 eo(maxo*2,nk)        : Eigenvalues
C real*8 qo(maxo*2,nk)        : Occupations of eigenstates
C real*8 Dnew(maxnd,4)        : Output Density Matrix
C real*8 Enew(maxnd,4)        : Output Energy-Density Matrix
C real*8 ef                   : Fermi energy
C real*8 Entropy              : Electronic entropy
C *************************** AUXILIARY *******************************
C complex*16 Haux(2,nuotot,2,nuo): Auxiliary space for the hamiltonian matrix
C complex*16 Saux(2,nuotot,2,nuo): Auxiliary space for the overlap matrix
C complex*16 psi(2,nuotot,2*nuo) : Auxiliary space for the eigenvectors
C complex*16 Dk(2,nuotot,nuo)    : Aux. space that may be the same as Haux
C complex*16 Ek(2,nuotot,nuo)    : Aux. space that may be the same as Saux
C complex*16 caux(2,nuotot)      : Extra auxiliary space
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
      
      real(dp) Dnew(maxnd,8), Enew(maxnd,4), H(maxnh,8), S(maxnh)
      real(dp) kpoint(3,nk), wk(nk), xij(3,maxnh)
      real(dp) e1, e2, ef, eo(maxo*2,nk), qo(maxo*2,nk)
      real(dp) occtol, qtot, temp, Entropy
      integer :: neigwanted
      
      complex(dp), dimension(2,nuotot,2,nuo) :: Haux, Saux, Dk, Ek
      complex(dp), dimension(2,nuotot,2*nuo), target :: psi
      complex(dp), dimension(2,nuotot), target :: caux
      logical, intent(in) ::               getD, getPSI
      

      real(dp),    pointer :: aux(:)
      real(dp), pointer  :: psi_real(:)

      integer  BNode, BTest, ie, ierror, iie, ik, ind, io, iio
      integer  iuo, j, jo, juo, neigneeded
      real(dp) ee, qe, t
      
!     Haux(js,juo,is,iuo) = <js,juo|H|is,iuo>
!     Indices is and js are for spin components
!     Indices iuo and juo are for orbital components
      complex(dp)                                 :: cicj
      complex(dp)                                 :: D11, D22, D12, D21
      complex(dp)                                 :: kphs
      real(dp)                                    :: kxij
#ifdef MPI
      integer            :: MPIerror
#endif

      external              cdiag
!***********************************************************************
!     B E G I N
!***********************************************************************

!****************************************************************************
! Build Hamiltonian      
!****************************************************************************
! (Notes by Jaime Ferrer)
! Indices iuo and juo are for orbital components in the unit cell
! while jo is for orbital components in the supercell:
! d_{jo} = R + d_{juo}, x_{ij} = d_{jo} - d_{iuo}
! Since <jo|H|iuo> is computed, the Bloch phase factor is e^(-i k x_{ij} )
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
! 2. Since wave functions are real, if there are no single P or L operators:
!      (a) H_{i,j}^{is,js}=H_{j,i}^{is,js}
!      (b) These imply spin-box hermiticity: H_{i,j}^{is,js}=H_{i,j}^{js,is}^*
!
!      (c) Hence
!
!                               | H(ind,1)                H(ind,3) - i H(ind,4) |
!           H_{j,i} = H_{i,j} = |                                               |
!                               | H(ind,3) + i H(ind,4)   H(ind,2)              |
!
!

      ! Jump point for bands and wfs output only
      if (getPSI .and. (.not. GetD)) goto 10
      
!     Find eigenvalues at every k point

      do ik = 1,nk
      
        Saux = cmplx(0.0_dp,0.0_dp,dp)
        Haux = cmplx(0.0_dp,0.0_dp,dp)
      
        do iuo = 1,nuo
          do j = 1,numh(iuo)
            ind = listhptr(iuo) + j
            jo = listh(ind)
            juo = indxuo(jo)
            kxij = kpoint(1,ik) * xij(1,ind) +
     .             kpoint(2,ik) * xij(2,ind) +
     .             kpoint(3,ik) * xij(3,ind)
            kphs = exp(cmplx(0.0_dp, -1.0_dp, dp) * kxij)
      
            Saux(1,juo,1,iuo) = Saux(1,juo,1,iuo) + S(ind)   * kphs
            Saux(2,juo,2,iuo) = Saux(2,juo,2,iuo) + S(ind)   * kphs
            Haux(1,juo,1,iuo) = Haux(1,juo,1,iuo) + 
     .                            cmplx(H(ind,1), H(ind,5),dp) * kphs
            Haux(2,juo,2,iuo) = Haux(2,juo,2,iuo) + 
     .                            cmplx(H(ind,2), H(ind,6),dp) * kphs
            Haux(1,juo,2,iuo) = Haux(1,juo,2,iuo)
     .                        + cmplx(H(ind,3), - H(ind,4),dp) * kphs
            Haux(2,juo,1,iuo) = Haux(2,juo,1,iuo)
     .                        + cmplx(H(ind,7), + H(ind,8),dp) * kphs
          enddo
        enddo
      
!       Find eigenvalues
!       Note duplication since the spins are mixed
        call cdiag(Haux,Saux,2*nuotot,2*nuo,2*nuotot,eo(1,ik),psi,
     .             -2*neigwanted,iscf,ierror, 2*BlockSize)
      
        if (ierror.ne.0) then
          call die('Terminating due to failed diagonalisation')
        endif
      enddo

      ! Check if we are done (i.e.: pure bands only)
      if (.not.getD .and. .not. getPSI) then
         RETURN
      end if
!-----------------------------------------------------------------------
      
!     Find new Fermi energy and occupation weights
      call fermid(2,1,nk,wk,maxo*2,neigwanted*2,eo,
     &    temp,qtot,qo,ef,Entropy)
      
!     Find weights for local density of states
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.d-6 )
        do ik = 1,nk
          do io = 1,neigwanted*2
            qo(io,ik) = wk(ik) * ( stepf( (eo(io,ik)-e2)/t ) -
     .                             stepf( (eo(io,ik)-e1)/t ) )
          enddo
        enddo
      endif
      
!     New density and energy-density matrices of unit-cell orbitals
      Dnew(:,:) = 0.0_dp
      Enew(:,:) = 0.0_dp

 10   continue   ! ---- Start here for (getPsi and (not getDM))
      
      do ik = 1,nk
         !  Find maximum eigenvector index that is required for this k point 
         if (getPSI) then
             ! We do not have info about occupation
             neigneeded = neigwanted*2
         else
             neigneeded = 0
             ie = 2*neigwanted ! for mixed spins
             do while (ie.gt.0.and.neigneeded.eq.0)
                if (abs(qo(ie,ik)).gt.occtol) neigneeded = ie
                ie = ie - 1
             enddo
         endif
      
!       Find eigenvectors 
        Saux = cmplx(0.0_dp,0.0_dp,dp)
        Haux = cmplx(0.0_dp,0.0_dp,dp)
        do iuo = 1,nuo
          do j = 1,numh(iuo)
            ind = listhptr(iuo) + j
            jo = listh(ind)
            juo = indxuo(jo)
            kxij = kpoint(1,ik) * xij(1,ind) +
     .             kpoint(2,ik) * xij(2,ind) +
     .             kpoint(3,ik) * xij(3,ind)
            kphs = exp(cmplx(0.0_dp, -1.0_dp, dp) * kxij)
      
            Saux(1,juo,1,iuo) = Saux(1,juo,1,iuo) + S(ind)   * kphs
            Saux(2,juo,2,iuo) = Saux(2,juo,2,iuo) + S(ind)   * kphs
            Haux(1,juo,1,iuo) = Haux(1,juo,1,iuo) + 
     .                            cmplx(H(ind,1), H(ind,5),dp) * kphs
            Haux(2,juo,2,iuo) = Haux(2,juo,2,iuo) + 
     .                            cmplx(H(ind,2), H(ind,6),dp) * kphs
            Haux(1,juo,2,iuo) = Haux(1,juo,2,iuo)
     .                        + cmplx(H(ind,3), - H(ind,4),dp) * kphs
            Haux(2,juo,1,iuo) = Haux(2,juo,1,iuo)
     .                        + cmplx(H(ind,7), + H(ind,8),dp) * kphs
          enddo
        enddo

        ! Note that neigneeded is explicitly calculated
        call cdiag(Haux,Saux,2*nuotot,2*nuo,2*nuotot,caux,psi,
     .             neigneeded,iscf,ierror, 2*BlockSize)
      
!       Check error flag and take appropriate action
        if (ierror.gt.0) then
          call die('Terminating due to failed diagonalisation')
        elseif (ierror.lt.0) then
!         Repeat diagonalisation with increased memory to handle clustering
          Saux = cmplx(0.0_dp,0.0_dp,dp)
          Haux = cmplx(0.0_dp,0.0_dp,dp)
          do iuo = 1,nuo
            do j = 1,numh(iuo)
              ind = listhptr(iuo) + j
              jo = listh(ind)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xij(1,ind) +
     .               kpoint(2,ik) * xij(2,ind) +
     .               kpoint(3,ik) * xij(3,ind)
              kphs = exp(cmplx(0.0_dp, -1.0_dp, dp) * kxij)
      
              Saux(1,juo,1,iuo) = Saux(1,juo,1,iuo) + S(ind)   * kphs
              Saux(2,juo,2,iuo) = Saux(2,juo,2,iuo) + S(ind)   * kphs
              Haux(1,juo,1,iuo) = Haux(1,juo,1,iuo) + 
     .                              cmplx(H(ind,1), H(ind,5),dp) * kphs
              Haux(2,juo,2,iuo) = Haux(2,juo,2,iuo) + 
     .                              cmplx(H(ind,2), H(ind,6),dp) * kphs
              Haux(1,juo,2,iuo) = Haux(1,juo,2,iuo)
     .                          + cmplx(H(ind,3), - H(ind,4),dp) * kphs
              Haux(2,juo,1,iuo) = Haux(2,juo,1,iuo)
     .                          + cmplx(H(ind,7), + H(ind,8),dp) * kphs
            enddo
          enddo
          call cdiag(Haux,Saux,2*nuotot,2*nuo,2*nuotot,caux,psi,
     .               neigneeded,iscf,ierror, 2*BlockSize)
        endif
      
          if (getPSI) then
             if (.not. getD) then
                ! Fill in the missing eigenvalue information
                ! This is useful if we are requesting WFS in a "bands" setting
                call c_f_pointer(c_loc(caux),aux,[2*nuotot])
                call dcopy(2*nuotot, aux(1), 1, eo(1,ik), 1)
             endif
             call c_f_pointer(c_loc(psi),psi_real,[2*size(psi)])
             call writew(nuotot,nuo,ik,kpoint(1,ik),1,     ! last '1' is 'spin-block number'
     .                   eo(1,ik),                         ! energies
     .                   psi_real,gamma=.false.,non_coll=.true.,
     $                   blocksize=2*BlockSize)
          endif

          if (.not. getD) CYCLE   ! process new k-point

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
! The Density Matrix is not "spin box hermitian" even if H does not contain P or
! L operators.
!
! Spin-box symmetrization of D:
! D_{i,j}(1,2) = 0.5 ( D_{i,j}(1,2) + D_{i,j}(2,1)^* )
! can not be enforced here.

        Dk  = cmplx(0.0_dp,0.0_dp,dp)
        Ek  = cmplx(0.0_dp,0.0_dp,dp)
      
        BNode = 0
        iie = 0
        do ie = 1,2*nuotot
          qe = qo(ie,ik)
          if (Node.eq.BNode) then
            iie = iie + 1
          endif
      
          caux(:,:) = cmplx(0.0_dp,0.0_dp,dp)
          if (abs(qe).gt.occtol) then
            if (Node.eq.BNode) then
              do j = 1,nuotot
                caux(1,j)= psi(1,j,iie) ! c_{i,up}
                caux(2,j)= psi(2,j,iie) ! c_{i,dn}
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
      
!------- 1,1 -----------------------------------------------------------
               cicj = caux(1,iio) * conjg(caux(1,juo))
               Dk(1,juo,1,iuo) = Dk(1,juo,1,iuo) + qe * cicj
               Ek(1,juo,1,iuo) = Ek(1,juo,1,iuo) + ee * cicj
!------- 2,2 -----------------------------------------------------------
               cicj = caux(2,iio) * conjg(caux(2,juo))
               Dk(2,juo,2,iuo) = Dk(2,juo,2,iuo) + qe * cicj
               Ek(2,juo,2,iuo) = Ek(2,juo,2,iuo) + ee * cicj
!------- 1,2 -----------------------------------------------------------
               cicj = caux(1,iio) * conjg(caux(2,juo))
               Dk(1,juo,2,iuo) = Dk(1,juo,2,iuo) + qe * cicj
               Ek(1,juo,2,iuo) = Ek(1,juo,2,iuo) + ee * cicj
!------- 2,1 -----------------------------------------------------------
               cicj = caux(2,iio) * conjg(caux(1,juo))
               Dk(2,juo,1,iuo) = Dk(2,juo,1,iuo) + qe * cicj
               Ek(2,juo,1,iuo) = Ek(2,juo,1,iuo) + ee * cicj
              enddo
            enddo
          endif
          BTest = ie/(2*BlockSize)
          if (BTest*2*BlockSize.eq.ie) then
            BNode = BNode + 1
            if (BNode .gt. Nodes-1) BNode = 0
          endif
        enddo
      
        do iuo = 1,nuo
        do j = 1,numd(iuo)
          ind = listdptr(iuo) + j
          jo = listd(ind)
          juo = indxuo(jo)
          kxij = kpoint(1,ik) * xij(1,ind) +
     .           kpoint(2,ik) * xij(2,ind) +
     .           kpoint(3,ik) * xij(3,ind)
          kphs = exp(cmplx(0.0_dp,-1.0_dp,dp) * kxij)
      
          D11 = Dk(1,juo,1,iuo) * kphs
          D22 = Dk(2,juo,2,iuo) * kphs
          D12 = Dk(1,juo,2,iuo) * kphs
          D21 = Dk(2,juo,1,iuo) * kphs
      
          Dnew(ind,1) = Dnew(ind,1) + real(D11,dp)
          Dnew(ind,2) = Dnew(ind,2) + real(D22,dp)
          Dnew(ind,3) = Dnew(ind,3) + real(D12,dp)
          Dnew(ind,4) = Dnew(ind,4) - aimag(D12)
          Dnew(ind,5) = Dnew(ind,5) + aimag(D11) 
          Dnew(ind,6) = Dnew(ind,6) + aimag(D22)
          Dnew(ind,7) = Dnew(ind,7) + real(D21,dp) 
          Dnew(ind,8) = Dnew(ind,8) + aimag(D21)

      
          D11 = Ek(1,juo,1,iuo) * kphs
          D22 = Ek(2,juo,2,iuo) * kphs
          D12 = Ek(1,juo,2,iuo) * kphs
          D21 = Ek(2,juo,1,iuo) * kphs
      
          Enew(ind,1) = Enew(ind,1) + real(D11,dp)
          Enew(ind,2) = Enew(ind,2) + real(D22,dp)
          Enew(ind,3) = Enew(ind,3) + real(D12,dp) 
          Enew(ind,4) = Enew(ind,4) - aimag(D12) 
      
        enddo
        enddo
      enddo
      
      end subroutine diag3k
!***********************************************************************      

      
