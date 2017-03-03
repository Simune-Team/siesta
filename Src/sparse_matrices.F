! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE sparse_matrices
      use precision
      use class_dSpData1D
      use class_dSpData2D
      use class_Sparsity
      use class_OrbitalDistribution
      use class_Fstack_Pair_Geometry_dSpData2D

      implicit none

      private

!     Max. number of nonzero H matrix elements    
      integer, public :: maxnh = 10 

      integer, public, pointer :: listh(:)=>null(), listhptr(:)=>null(),
     $                            numh(:)=>null()

      real(dp), public, pointer :: Dold(:,:)=>null(), Dscf(:,:)=>null()
      real(dp), public, pointer :: Eold(:,:)=>null(), Escf(:,:)=>null()
      real(dp), public, pointer :: Hold(:,:)=>null(), H(:,:)=>null()

      real(dp), public, pointer :: S(:)=>null()

      real(dp), public, pointer :: xijo(:,:)=>null()

! CC RC  Added for the offSpOrb
      complex(dp), public, pointer :: H0_SO(:,:) => null()
      integer, public, pointer :: listht(:)=>null()
! CC RC  Added for the offSpOrb

      ! Pieces of H that do not depend on the SCF density matrix
      ! Formerly there was a single array H0 for this
      type(dSpData1D), public, save :: H_vkb_1D, H_kin_1D
      ! LDA+U and spin-orbit coupling Hamiltonian
      type(dSpData2D), public, save :: H_ldau_2D, H_so_2D

      !  New interface data
      type(Sparsity), public, save :: sparse_pattern
      type(dSpData1D), public, save :: S_1D
      type(dSpData2D), public, save :: DM_2D, EDM_2D, xij_2D, H_2D
      type(Fstack_Pair_Geometry_dSpData2D), public, save :: DM_history

      ! Orbital distribution descriptor
      type(OrbitalDistribution), public, save :: block_dist
      ! Always store a single distributor
      type(OrbitalDistribution), public, save :: single_dist

      public :: resetSparseMatrices

! CC RC  Added for the offSpOrb
      public :: MM_HDM, herm_HDM, euler
! CC RC  Added for the offSpOrb

      CONTAINS

      subroutine resetSparseMatrices( )
      use alloc, only : de_alloc

      implicit none

      call delete( block_dist )
      call delete( sparse_pattern )
      nullify(numh,listhptr,listh)

      call delete( H_kin_1D )
      call delete( H_vkb_1D )
      call delete( H_ldau_2D )
      call delete( H_so_2D )

      call delete( DM_2D ) ; nullify(Dscf)
      call delete( EDM_2D ) ; nullify(Escf)
      call delete( S_1D ) ; nullify(S)
      call delete( H_2D ) ; nullify(H)
      call delete( xij_2D ) ; nullify(xijo)

      call delete( DM_history )

! CC RC  Added for the offSpOrb
      deallocate( H0_SO )
      call de_alloc( listht,'listht','sparseMat' ) ! CC RC 
! CC RC  Added for the offSpOrb

      
      ! Using MixH is bad as several utilities
      ! are not relying on FoX, better to leave out
      ! Note that it is nullified from start, 
      ! and hence de_alloc returns immediately
      
      call de_alloc( Hold, 'Hold', 'sparseMat' )
      call de_alloc( Dold, 'Dold', 'sparseMat' )
      call de_alloc( Eold, 'Eold', 'sparseMat' )

      end subroutine resetSparseMatrices

! CC RC  Added for the offSpOrb
      subroutine MM_HDM( cHDM )

      use m_spin,           only : spin
      use atomlist,         only : lasto
      use siesta_geom,      only : na_u
      use precision

      implicit none

      character(len=*), intent(in) :: cHDM

      integer    :: isp, ia, i, nh, ioff, it
      real(dp)   :: Qs(2), Q_tot, MM(0:3), theta, phi
      real(dp)   :: r_xy, r_mod, r_xy2

      real(dp), pointer       :: HDM(:,:)
      real(dp), allocatable   :: MM_at(:,:)

      character(len=1) :: cindent = ' '

      real(dp), parameter :: pi = 4.0d0*atan(1.0_dp)
      real(dp), parameter :: fac= 1.8d2/pi

      real(dp), external :: ddot

c-----------------------------------------------------------------------

      if ( spin%Grid.eq.1 ) return

c---- set pointers
      select case( cHDM )
      case( 'DM' )
       HDM   => Dscf
      case( 'H' )
       HDM   => H
      case default
       write(6,'(a)')
     &  'MM_HDM: internal ERROR, wrong type='//trim(cHDM)
       stop
      end select

      allocate( MM_at(0:3,na_u) )

!c--- atom resolved
c      write(spin%iout_SO,'(a)') trim(cindent)//
c     & '    :  atom  MM_total   MM_x     MM_y     MM_z   THETA   PHI'
      write(6,'(a)') trim(cindent)//
     & '    :  atom  MM_total   MM_x     MM_y     MM_z   THETA   PHI'
      MM_at(:,:)= 0.0_dp; ioff= 0
      theta= 0.0_dp; phi= 0.0_dp
      do ia = 1, na_u
       Qs(:)= 0.0_dp
       nh = sum(numh(lasto(ia-1)+1:lasto(ia)))

       do isp= 1, min(2,spin%Grid)
        Qs(isp) = ddot( nh, HDM(ioff+1:,isp), 1, S(ioff+1:), 1 )
       enddo
       MM_at(3,ia) = Qs(1)-Qs(2)

       if ( spin%Grid.le.2 ) then
        MM_at(0,ia)= abs(MM_at(3,ia))
       else
        MM_at(1,ia)= 2.0d0*ddot(nh,HDM(ioff+1:,3),1,S(ioff+1:),1)
        MM_at(2,ia)=-2.0d0*ddot(nh,HDM(ioff+1:,4),1,S(ioff+1:),1)
        MM_at(0,ia) = sqrt( sum(MM_at(1:3,ia)**2) )
        call euler( MM_at(1:3,ia), theta, phi, MM_at(0,ia) )
       endif

       ioff = ioff + nh
c       write(spin%iout_SO,'(a,5x,i5,4f9.3,2f7.1)') trim(cindent),
c     &            ia,MM_at(0:3,ia), theta*fac, phi*fac
       write(6,'(a,5x,i5,4f9.3,2f7.1)') trim(cindent),
     &            ia,MM_at(0:3,ia), theta*fac, phi*fac
      enddo ! ia

c---- total MMs
      MM(3)= sum(MM_at(3,:))
      if ( spin%Grid.eq.2 ) then
       MM(0) = abs(MM(3))
      else
       MM(1)= sum(MM_at(1,:))
       MM(2)= sum(MM_at(2,:))
       call euler( MM(1:3), theta, phi, MM(0) )
      endif
C      write(spin%iout_SO,'(a,4f9.3,a,2f7.1,a)') trim(cindent)//
C     &  '    : '//trim(cHDM)//' Total Magnetic Moment=',MM(0:3),
C     &         ' [theta/phi=',theta*fac,phi*fac,']'
      write(6,'(a,4f9.3,a,2f7.1,a)') trim(cindent)//
     &  '    : '//trim(cHDM)//' Total Magnetic Moment=',MM(0:3),
     &         ' [theta/phi=',theta*fac,phi*fac,']'

      return
      END subroutine MM_HDM
c-----------------------------------------------------------------------
c
!> Finds Euler angles
c
c-----------------------------------------------------------------------
      subroutine euler( vec, theta, phi, r_mod )

      use precision

      implicit none

      real(dp), intent(in)  :: vec(3)
      real(dp), intent(out) :: r_mod, theta, phi

      real(dp) :: r_xy2, r_xy

c-----------------------------------------------------------------------

      theta= 0.0_dp; phi= 0.0_dp
      r_xy2= vec(1)**2 + vec(2)**2
      r_mod= sqrt( r_xy2 + vec(3)**2 )
      if ( r_mod .lt. 1.0d-6 ) return

      r_xy = sqrt( r_xy2 )
      theta= atan2( r_xy,vec(3) )
      if ( r_xy.gt.1.0d-6 ) phi = atan2( vec(2), vec(1) )

      return
      end subroutine euler
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      subroutine herm_HDM( cHDM )

      use units,     only: eV
      use m_spin,    only: spin
      USE atomlist,  only: no_u

      implicit none

      character(len=*), intent(in) :: cHDM

      real(dp)  :: err_HDM(4)

      complex(dp):: H_tmp(maxnh,4)

      real(dp), pointer :: HDM(:,:)

      integer :: ij, ji, i, j, js

      select case( trim(cHDM) )
      case( 'H' )
       HDM => H
      case( 'DM' )
       HDM => Dscf
      end select

      H_tmp(:,1) = cmplx( HDM(:,1), HDM(:,5) )
      H_tmp(:,2) = cmplx( HDM(:,2), HDM(:,6) )
      H_tmp(:,3) = cmplx( HDM(:,3), HDM(:,4) )
      H_tmp(:,4) = cmplx( HDM(:,7), HDM(:,8) )

      err_HDM = 0.0_dp

!! CC RC
!      do i = 1, no_u
!       do js = 1, numh(i)
!        ij = listhptr(i)+js
!        j = listh(ij)
!        ji = listhptr(j)+i
!
!       if ( ji.lt.ij ) then
!        err_HDM(1) = max( err_HDM(1),
!     &                  cdabs(H_tmp(ij,1)-conjg(H_tmp(ji,1))) )
!        err_HDM(2) = max( err_HDM(2),
!     &                  cdabs(H_tmp(ij,2)-conjg(H_tmp(ji,2))) )
!       endif
!       err_HDM(3) = max( err_HDM(3),
!     &                cdabs(H_tmp(ij,3)-conjg(H_tmp(ji,4))) )
!       err_HDM(4) = max( err_HDM(4),
!     &                cdabs(H_tmp(ij,4)-conjg(H_tmp(ji,3))) )
!
!       enddo
!      enddo

      do ij = 1, maxnh
       ji = listht(ij)
       if ( ji.lt.ij ) then
        err_HDM(1) = max( err_HDM(1),
     &                  cdabs(H_tmp(ij,1)-conjg(H_tmp(ji,1))) )
        err_HDM(2) = max( err_HDM(2),
     &                  cdabs(H_tmp(ij,2)-conjg(H_tmp(ji,2))) )
       endif
       err_HDM(3) = max( err_HDM(3),
     &                cdabs(H_tmp(ij,3)-conjg(H_tmp(ji,4))) )
       err_HDM(4) = max( err_HDM(4),
     &                cdabs(H_tmp(ij,4)-conjg(H_tmp(ji,3))) )
      enddo

c---- Print errors
!      write(spin%iout_SO,'(4x,a,7f8.5)')
      write(6,'(4x,a,7f8.5)')
     &   trim(cHDM)//'  max/sym error  =',
     &   maxval(cdabs(H_tmp))/eV,err_HDM
      if ( maxval(err_HDM).gt.1.0d-4 ) then
       write(6,*) 'ERROR, '//cHDM//' not hermitian'
       stop
      endif

      end subroutine herm_HDM
! CC RC  Added for the offSpOrb

      END MODULE sparse_matrices
