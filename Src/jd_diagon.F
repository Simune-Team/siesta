      module m_jd_diagon
      implicit none
      public :: jd_diagon

      CONTAINS
      subroutine jd_diagon( nspin, nuotot, nuo, neig, numh, listh, H, S,
     $                      Dscf, Escf, ef, efs, Entrop)

!     Work in progress ...

!     Diagonalizes
!     Builds density matrix

      integer,  intent(in) :: nspin, nuotot, nuo, neig
      integer,  intent(in) :: numh(*), listh(*)
      real(dp), intent(in) :: H(*,nspin), S(*)
      real(dp), intent(out)          :: Dscf(:,nspin), Escf(:,nspin)
      real(dp), intent(out)          :: ef, efs, Entrop

      real(dp),    pointer :: eigenvalues(:,:), psi(:,:,:)

      nullify( eigenvalues, psi )
      call re_alloc( eigenvalues, 1, neig, 1, nspin, 'w', 'jd_diag' )
      call re_alloc( psi, 1, nuo, 1, neig, 1, nspin, 'psi', 'jd_diagon' )

      do ispin = 1, nspin
         call jd_diag_wZ( nuotot, nuo, neig, numh, listh, H(:,ispin), S,
     $                    eigenvalues(:,ispin), psi(:,:,ispin))

      enddo

!     Use eigenvalues and psi to build Dscf and Escf
         
C Find new Fermi energy and occupation weights ........................
      if (fixspin) then
        call fermispin( nspin, nspin, nk, wk, maxo, neig, 
     .                  eigenvalues, temp, qs, qo, efs, Entropy )
      else
        call fermid( nspin, nspin, nk, wk, maxo, neig, eigenvalues, 
     .               temp, qtot, qo, ef, Entropy )
      endif


C New density and energy-density matrices of unit-cell orbitals .......
      nd = listdptr(nuo) + numd(nuo)
      Dnew(1:nd,1:nspin) = 0.0d0
      Enew(1:nd,1:nspin) = 0.0d0

C Global operation to form new density matrix
      nullify( nuo_LOC )
      call re_alloc( nuo_LOC, 0, Nodes-1, 'nuo_LOC', 'diagg' )
#ifdef MPI
      call MPI_Allgather( nuo, 1, MPI_INTEGER, nuo_LOC, 1,
     &                    MPI_INTEGER, MPI_Comm_World, MPIerror )
#else
      nuo_LOC(0) = nuo
#endif
      maxnuo = 0
      do BNode=0, Nodes-1
        maxnuo = max(nuo_LOC(BNode),maxnuo)
      enddo
      nullify( PSI_TMP )
      call re_alloc( PSI_TMP, 1, nuotot*maxnuo*nspin, 'PSI_TMP',
     &               'diagg' )

      do Bnode=0, Nodes-1
        if (BNode.eq.Node) then
          mm = 0
          do ispin= 1, nspin
            do io= 1, nuo
              PSI_TMP(mm+1:mm+nuotot) = PSI(1:nuotot,io, ispin)
              mm = mm + nuotot
            enddo
          enddo
        endif
#ifdef MPI
        call MPI_Bcast( PSI_TMP, nuotot*nuo_LOC(Bnode)*nspin,
     &                  MPI_double_precision, BNode, MPI_Comm_World,
     &                  MPIerror )
#endif
        do ispin = 1,nspin
          do io = 1,nuo
            call LocalToGlobalOrb(io,Node,Nodes,iio)
            mm   = (ispin-1)*nuo_LOC(Bnode)*nuotot
            do iie = 1,nuo_LOC(Bnode)
              call  LocalToGlobalOrb(iie,BNode,Nodes,ie)
              qe = qo(ie,ispin)
              if (qe.gt.occtol) then
                paux => PSI_TMP(mm+1:mm+nuotot)
                ee  = qe*eo(ie,ispin)
                qei = qe*paux(iio)
                eei = ee*paux(iio)
                do j = 1,numd(io)
                  ind = listdptr(io) + j
                  jo  = listd(ind)
                  Dnew(ind,ispin) = Dnew(ind,ispin) + qei*paux(jo)
                  Enew(ind,ispin) = Enew(ind,ispin) + eei*paux(jo)
                enddo
              endif
              mm = mm + nuotot
            enddo
          enddo
        enddo
      enddo
      call de_alloc( nuo_LOC, 'nuo_LOC', 'diagg' )
      call de_alloc( PSI_TMP, 'PSI_TMP', 'diagg' )


      call de_alloc( psi, 'psi', 'jd_diagon' )
      call de_alloc( eigenvalues, 'eigenvalues', 'jd_diagon' )

      end subroutine jd_diagon

      end module m_jd_diagon
