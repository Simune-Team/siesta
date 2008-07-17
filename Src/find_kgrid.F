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
      module m_find_kgrid
      private
      public :: find_kgrid
      CONTAINS
      subroutine find_kgrid ( cell, kscell, displ, firm_displ,
     $                        time_reversal_symmetry,
     $                        nk, points, weight, eff_kgrid_cutoff )

c ***************** INPUT **********************************************
c real*8  cell(3,3)  : Unit cell vectors in real space cell(ixyz,ivec)
c integer kscell(3,3): Supercell reciprocal of k-grid unit cell
c                      scell(ix,i) = sum_j cell(ix,j)*kscell(j,i)
c logical time_reversal_symmetry : Flag to fold k to -k if possible. 
c ***************** INPUT/OUTPUT  **************************************
c real*8  displ(3)   : Grid origin in k-grid-vector coordinates:
c                      origin(ix) = sum_j gridk(ix,j)*displ(j)
c ***************** OUTPUT *********************************************
c integer nk           : Actual number of irreducible k-points
c real(dp) points(3:nk): Kpoints
c real(dp) weight(nk)  : weights
c real(dp) eff_kgrid_cutoff : actual equivalent kgrid cutoff 

C  Modules

      use precision,  only : dp
      use units,      only : pi
      use m_minvec,   only : minvec
      use alloc,      only : re_alloc
      use parallel,   only : Node

      implicit          none

C Passed variables
      integer, intent(in)      :: kscell(3,3)
      real(dp), intent(in)     :: cell(3,3)
      logical, intent(in)      :: firm_displ
      logical, intent(in)      :: time_reversal_symmetry

      real(dp), intent(inout)  :: displ(3)

      integer, intent(out)     :: nk
      real(dp), intent(out)    :: eff_kgrid_cutoff
      real(dp), pointer        :: points(:,:)
      real(dp), pointer        :: weight(:)

      external                 ::  idiag, reclat

C Internal variables
      integer           i, ir, ix, j, i1, i2, i3,
     .                  kdsc(3,3), maux(3,3,2), ml(3,3), mr(3,3),
     .                  ng(3), ni, nkr(3), nktot,
     $                  proj(3,3), igmin(3), igmax(3)

      real(dp)          d(3), dkg(3), dkx(3), dscell(3,3),
     .                  gridk(3,3), gscell(3,3), 
     .                  scell(3,3), scmin(3,3), tiny, vmod, w1, wtot
      real(dp)          ctransf(3,3)

      parameter (tiny   = 1.d-12)

C Find total number of points (determinant of kscell)
      nktot = abs( kscell(1,1) * kscell(2,2) * kscell(3,3) +
     .             kscell(2,1) * kscell(3,2) * kscell(1,3) +
     .             kscell(3,1) * kscell(1,2) * kscell(2,3) -
     .             kscell(1,1) * kscell(3,2) * kscell(2,3) -
     .             kscell(2,1) * kscell(1,2) * kscell(3,3) -
     .             kscell(3,1) * kscell(2,2) * kscell(1,3) )

C 
C Find k-grid supercell
C
      do i = 1,3
        do ix = 1,3
          scell(ix,i) = cell(ix,1) * kscell(1,i) +
     .                  cell(ix,2) * kscell(2,i) +
     .                  cell(ix,3) * kscell(3,i)
        enddo
        vmod = sqrt( scell(1,i)**2 + scell(2,i)**2 + scell(3,i)**2 )
      enddo

C Find actual cutoff
      call minvec( scell, scmin, ctransf )
      eff_kgrid_cutoff = huge(1.0_dp)
      do i = 1,3
        vmod = sqrt( scmin(1,i)**2 + scmin(2,i)**2 + scmin(3,i)**2 )
        eff_kgrid_cutoff = min( eff_kgrid_cutoff, vmod/2.d0 )
      enddo

C     Equivalent supercell DA with the property that there exists
C     a primitive cell (pa') such that DA_i = N_i*pa'_i
C     (See Moreno and Soler)
!!    Direct route
!!    call DIGCEL( ucell, kSCELL, new_ucell, dscell, NSC, ISDIAG )

      call idiag( 3, kscell, kdsc, ml, mr, maux )
      proj(:,:) = 0  ! Possible sign changes
      do i = 1, 3
         proj(i,i) = 1
         if (kdsc(i,i) < 0) proj(i,i) = -1
      enddo
      kdsc = matmul(kdsc,proj)
      mr = matmul(mr,proj)
C
C     Set the displacements if not firm (i.e., specified by the
C     user). Even if firm, warn if a better choice is possible.
C
      do j = 1, 3
         if (mod(kdsc(j,j),2) .eq. 0) then
            if (firm_displ .and. displ(j) /= 0.5d0) then
               if (Node .eq. 0)
     $           write(6,"(a,i4,a,2f8.2)")
     $            "k-point displ. along", j, " input, could be: ",
     $           displ(j), 0.5d0
            else
               displ(j) = 0.5d0
            endif
         else
            if (firm_displ .and. displ(j) /= 0.0d0) then
               if (Node .eq. 0)
     $           write(6,"(a,i4,a,2f8.2)")
     $            "k-point displ. along", j, " input, could be: ",
     $           displ(j), 0.0d0
            else
               displ(j) = 0.0d0
            endif
         endif
      enddo
C
      dscell = matmul(scell,mr)

C Find k-grid unit vectors
      call reclat( dscell, gridk, 1 )

C Find grid origin in cartesian coordinates
      call reclat( scell, gscell, 1 )
      do ix = 1,3
        dkx(ix) = gscell(ix,1) * displ(1) +
     .            gscell(ix,2) * displ(2) +
     .            gscell(ix,3) * displ(3)
      enddo

C Find grid origin in gridk coordinates

      do i = 1,3
        dkg(i) = ( dkx(1) * dscell(1,i) +
     .             dkx(2) * dscell(2,i) +
     .             dkx(3) * dscell(3,i) ) / (2*pi)
      enddo

C Find total range of grid indexes
      do j = 1,3
        ng(j) = kdsc(j,j)
        igmin(j) = -( (ng(j)-1) / 2)
        igmax(j) = ng(j) / 2
      enddo

C Find number of points with time-reversal (inversion) symmetry,
C (if possible) after reflection on each alternative plane
      if (.not. time_reversal_symmetry) then
        ! Use all k points
        nk = nktot
      else
        do j = 1,3
          ni = ng(j)
          if (abs(dkg(j)) .lt. tiny) then
            ni = ng(j)/2 + 1
          elseif (abs(dkg(j)-0.5d0) .lt. tiny) then
            ni = (ng(j)-1)/2 + 1
          endif
          nkr(j) = ni * nktot / kdsc(j,j)
        enddo

C Select reflection plane
        ir = 3
        if (nkr(2) .lt. nkr(ir)) ir = 2
        if (nkr(1) .lt. nkr(ir)) ir = 1
        igmin(ir) = 0
        if (abs(dkg(ir)-0.5d0) .lt. tiny)
     .    igmax(ir) = (ng(ir)-1)/2
        nk = nkr(ir)
      endif

      call re_alloc(points,1,3,1,nk,name='points',
     $              routine='find_kgrid',copy=.false.)
      call re_alloc(weight,1,nk,name='weight',
     $              routine='find_kgrid',copy=.false.)

C Find k points and weights
      w1 = 1.0d0 / nktot
      nk = 0
      do i3 = igmin(3),igmax(3)
      do i2 = igmin(2),igmax(2)
      do i1 = igmin(1),igmax(1)
        nk = nk + 1
        d(1) = i1 + dkg(1)
        d(2) = i2 + dkg(2)
        d(3) = i3 + dkg(3)
        if (d(1) .gt. 0.5d0*ng(1)+tiny) d(1) = d(1) - ng(1)
        if (d(2) .gt. 0.5d0*ng(2)+tiny) d(2) = d(2) - ng(2)
        if (d(3) .gt. 0.5d0*ng(3)+tiny) d(3) = d(3) - ng(3)
        do ix = 1,3
          points(ix,nk) = gridk(ix,1)*d(1) + 
     .                    gridk(ix,2)*d(2) +
     .                    gridk(ix,3)*d(3)
        enddo
        if (.not. time_reversal_symmetry) then
          weight(nk) = w1
        else
          if ( abs(d(ir))              .lt. tiny .or.
     .         abs(d(ir)-0.5d0*ng(ir)) .lt. tiny) then
            weight(nk) = w1
          else
            weight(nk) = 2.0d0 * w1
          endif
        endif
      enddo
      enddo
      enddo

C Check that weight is normalised
      wtot = sum(weight(1:nk))
      if (abs(wtot-1.0d0) .gt. nk*tiny) then
        w1 = dble(nk)/wtot
        weight(1:nk) =  w1*weight(1:nk)
      endif

      end subroutine find_kgrid
      end module m_find_kgrid
