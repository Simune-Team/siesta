!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2014, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

module m_ts_mumps_scat

  use precision, only : dp

#ifdef SIESTA__MUMPS

  private

  public :: GF_Gamma_GF
#ifdef USE_GEMM3M
# define GEMM zgemm3m
#else
# define GEMM zgemm
#endif

contains

  subroutine GF_Gamma_GF(El, mum, no_u_TS, no, GF)
    use iso_c_binding, only: c_loc, c_f_pointer
    use m_ts_electype

    implicit none

    include 'zmumps_struc.h'

! *********************
! * INPUT variables   *
! *********************
    ! electrode self-energy
    type(Elec), intent(in) :: El
    type(zMUMPS_STRUC), intent(inout) :: mum
    integer, intent(in) :: no_u_TS ! no. states in contact region
    integer, intent(in) :: no      ! no. states for all electrodes
    ! The Green function (it has to be the column that corresponds to the electrode)
    complex(dp), intent(inout) :: GF(no_u_TS,no)

! *********************
! * LOCAL variables   *
! *********************
    complex(dp), parameter :: z0 = cmplx(0._dp, 0._dp,dp)
    complex(dp), parameter :: z1 = cmplx(1._dp, 0._dp,dp)
    complex(dp), parameter :: zi = cmplx(0._dp, 1._dp,dp)

    complex(dp), pointer :: rows(:,:), ztmp(:)
    integer :: io, jo, ind, indG, SB, CB

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE GFGammaGF' )
#endif

    call timer("GFGGF",1)

    ! re-use the work-array in mum%S
    if ( .not. associated(mum%S) ) &
        call die('Error in MUMPS, S not allocated')

    ! Calculate maximum blocking size
    SB = size(mum%S) / (no_u_TS + no)
    ! The block must not be larger than no_u_TS
    ! MUMPS working array could be larger than the entire matrix :(
    SB = min(SB,no_u_TS)
    if ( SB <= 0 ) call die('We cannot calculate the Gf.G.Gf^\dagger &
        &product for your system.')

    ! Setup the different segments
    call c_f_pointer(c_loc(mum%S),rows,[SB, no_u_TS])
    ztmp => mum%S(no_u_TS*SB+1:no_u_TS*SB+no*SB) ! Gf.Gamma

    ! Loop over rows
    ind = 1
    CB = 1
    do while ( CB <= no_u_TS )
      call c_f_pointer(c_loc(mum%S),rows,[SB, no_u_TS])

      ! We will now the the MM of Gf[CB:CB-1+SB,<electrode>].Gamma.Gf^\dagger = A
      ! with bounds A[CB:CB-1+SB,:]

      ! Copy over rows of Gf
      do jo = 1 , no
        rows(:,jo) = Gf(CB:CB-1+SB,jo)
      end do
      ! now rows(:,1:no) contains the Gf rows at CB:CB-1+SB

      ! This routine will save the transposed of: Gf.Gamma.Gf^\dagger
      ! as that fits with the siesta sparsity pattern

      ! Do Gf.Gamma for these rows
      call GEMM ('N','T',SB,no,no,z1, &
          rows(1,1), SB, &
          El%Gamma, no, &
          z0, ztmp(1), SB)

      ! Calculate the Gf.Gamma.Gf^\dagger product for the entire rows
      call GEMM ('N','C',SB,no_u_TS,no,z1, &
          ztmp(1), SB, &
          Gf(1,1), no_u_TS, &
          z0, rows(1,1),  SB)

      ! Backtrack until we have the correct row
      do while ( CB <= mum%IRN(ind) )
        ind = ind - 1
        if ( ind == 0 ) exit
      end do

      ! We now have the full Gf.G.Gf^\dagger rows at CB:CB-1+SB
      ind = ind + 1
      do while ( mum%IRN(ind) <= CB + SB - 1)
        io = mum%IRN(ind) - CB + 1
        jo = mum%JCN(ind)

        ! save the row
        mum%A(ind) = rows(io,jo)

        ind = ind + 1 ! update index
        if ( ind > mum%NZ ) exit

      end do

      ! Update what we have calculated
      CB = CB + SB ! this is the current reached row

      ! Calculate the next block size
      ! Only change the blocking size if
      ! we need fewer than SB rows
      if ( CB + SB - 1 > no_u_TS ) then
        SB = no_u_TS - CB + 1
        call c_f_pointer(c_loc(mum%S),rows,[SB, no_u_TS])
      end if

    end do
    
    ! Check that we actually calculated all entries
    if ( ind < mum%NZ ) then
      call die('We did not complete calculation')
    end if

    call timer("GFGGF",2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS GFGammaGF' )
#endif

  end subroutine GF_Gamma_GF

#undef GEMM

#endif

end module m_ts_mumps_scat
