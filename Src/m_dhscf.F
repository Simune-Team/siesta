! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module m_dhscf
C     Modules
      use precision
      implicit none

C     MODULE variables
C     Arrays that hold data which do not change during the SCF loop
      real(grid_p),   pointer :: rhopcc(:), rhoatm(:), Vna(:)
      real(dp)                :: Uharrs     ! Harris energy
      logical                 :: IsDiag, IsEfld, spiral
      character(len=10)       :: shape
      integer                 :: nm(3), nml(3), ntml(3), npcc,
     &                           nmpl, ntpl, ntopl
      real(dp)                :: qspiral(3), bcell(3,3), cell(3,3),
     &                           dvol, field(3), rmax, scell(3,3)
      real(dp)                :: G2mesh =  0.0_dp

      interface
      subroutine dhscf_init(nspin, norb, iaorb, iphorb,
     &                      nuo, nuotot, nua, na,
     &                      isa, xa, indxua, ucell,
     &                      mscell, G2max, ntm,
     &                      maxnd, numd, listdptr, listd, datm,
     &                      Fal, stressl)
      use precision
      implicit none
      integer, intent(in)     :: nspin, norb, iaorb(norb), iphorb(norb),
     &                           nuo, nuotot, nua, na, isa(na),
     &                           indxua(na), mscell(3,3), maxnd,
     &                           numd(nuo), listdptr(nuo), listd(maxnd)
      real(dp), intent(in)    :: xa(3,na), ucell(3,3), G2max, datm(norb)
      integer, intent(inout)  :: ntm(3)
      real(dp), intent(inout) :: Fal(3,nua), stressl(3,3)
      end subroutine dhscf_init
      END interface

      interface
      subroutine dhscf( nspin, norb, iaorb, iphorb, nuo,
     &                  nuotot, nua, na, isa, xa, indxua,
     &                  ntm, ifa, istr, iHmat, filrho, fildrh,
     &                  filevh, filevt, filepsch, filetoch, maxnd, numd,
     &                  listdptr, listd, Dscf, datm, maxnh, Hmat,
     &                  Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     &                  Exc, Dxc, dipol, stress, Fal, stressl )
      use precision
      implicit none
      integer
     &  maxnd, maxnh, nua, na, norb, nspin, nuo, nuotot,
     &  iaorb(norb), ifa, iHmat, indxua(na), iphorb(norb), isa(na), 
     &  istr, listd(*), listh(maxnh), listdptr(nuo), listhptr(nuo),
     &  mscell(3,3), ntm(3), numd(nuo), numh(nuo)

      real(dp)
     &  datm(norb), dipol(3), Dscf(maxnd,nspin),
     &  DUext, DUscf, Dxc, Enaatm, Enascf, Exc,
     &  Hmat(maxnh,nspin), Uatm, ucell(3,3), Uscf, xa(3,na),
     &  stress(3,3), Fal(3,nua), stressl(3,3)

      character
     &  filevh*(*), filevt*(*), fildrh*(*), filrho*(*), 
     &  filepsch*(*), filetoch*(*)
      end subroutine dhscf
      END interface

      end module m_dhscf
