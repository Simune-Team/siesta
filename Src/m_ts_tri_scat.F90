!
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
!
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

module m_ts_tri_scat

  use precision, only : dp

  implicit none

  private

  public :: GF_Gamma_GF
  public :: GFGGF_needed_worksize
  public :: ts_needed_mem

  ! Used for BLAS calls (local variables)
  complex(dp), parameter :: z0  = dcmplx( 0._dp, 0._dp)
  complex(dp), parameter :: z1  = dcmplx( 1._dp, 0._dp)
  complex(dp), parameter :: zm1 = dcmplx(-1._dp, 0._dp)

contains

  ! The problem of this routine is that we wish not to
  ! overwrite the old half-inverted matrix, that would mean
  ! that we need to do the full calculation for each electrode
  ! (which can be quite time-consuming!)
  subroutine GF_Gamma_GF(Gf_tri, El, nwork, work)

    use alloc, only : re_alloc, de_alloc

    use class_zTriMat
    use m_ts_trimat_invert, only : TriMat_Bias_idxs
    use m_ts_electype

    implicit none

! *********************
! * INPUT variables   *
! *********************
    ! The Green's function column
    type(zTriMat), intent(inout) :: Gf_tri
    type(Elec), intent(in) :: El ! contains: i (Sigma - Sigma^dagger)/2 ^T

! *********************
! * OUTPUT variables  *
! *********************
    integer, intent(in) :: nwork
    complex(dp), intent(inout), target :: work(nwork)

    ! local variables
    complex(dp), pointer :: fGf(:), Gf(:), GGG(:), oW(:)
    integer :: nr, np, no, tn
    integer :: sIdx, eIdx
    integer :: last_eIdx
    logical :: need_alloc
    integer :: ip, cp, n
    integer :: sN, sNc

    integer :: lsPart, lePart
    integer :: BsPart, BePart

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE GFGammaGF' )
#endif

    call timer("GFGGF",1)

    ! tri-diagonal parts information
    nr = nrows_g(Gf_tri)
    no = TotUsedOrbs(El)
    np = parts(Gf_tri)

    ! Which parts are needed
    lsPart = 1
    lePart = np

    ! Capture the full elements
    fGf => val(Gf_tri)

    ip = 0
    do n = lsPart , lePart
       ip = max(ip,no * nrows_g(Gf_tri,n))
    end do

    if ( nwork < ip ) &
         call die('Work size not big enough')
    oW => work(:)
    
    ! to track the last overwrite place
    tn = elements(Gf_tri) + 1

    do n = lePart , lsPart , -1

       ! Calculate the \Gamma Gf^\dagger sPart,1
       sN = nrows_g(Gf_tri,n)

       if ( size(oW) < no * sN ) then
          call die('Something went wrong with calculating &
               &the maximum work size')
       end if

       ! find the index of the thing that we don't want
       ! to overwrite...
       call TriMat_Bias_idxs(Gf_tri,no,min(n+1,lePart),sIdx,last_eIdx)
       
       call TriMat_Bias_idxs(Gf_tri,no,n,sIdx,eIdx)
       ! obtain the Gf in the respective column
       Gf => fGf(sIdx:eIdx)
       call zgemm('T','C',no,sN,no, z1, El%Gamma, no, &
            Gf, sN, z0, oW, no)
       
       ! Now we are ready to perform the multiplication
       ! for the requested region

       ! correct to the quantities that is available
       BsPart = max(n-1,lsPart)
       BePart = min(n+1,lePart)

#ifdef TRANSIESTA_DEBUG
       write(*,'(a,2(tr1,i0),a,2(tr1,i0))')'GfGGf at:',BsPart,ip,' --',BePart,ip
#endif
       
       ! this will populate in decreasing column major order
       do cp = BePart , BsPart , -1

          sNc = nrows_g(Gf_tri,cp)
          
          ! Update the index of which we will update last
          tn = tn - sNc * sN

          if ( tn <= last_eIdx ) then
             ! transfer to the work-array

             ! Retrieve Gf block
             call TriMat_Bias_idxs(Gf_tri,no,max(cp,n),sIdx,eIdx)
             
             ! copy over the elements in the end
             work(nwork-eIdx+1:nwork) = fGf(1:eIdx)
             ! point to the new place of the Gf-column
             fGf => work(nwork-eIdx+1:nwork)

             ! restrict work-array to be the remaining size
             ! lets us check that what we do is correct
             oW => work(1:nwork-eIdx)
             if ( size(oW) < no * sN ) then
                call die('Something went wrong with calculating &
                     &the maximum work size')
             end if

          end if

          ! Retrieve Gf block
          call TriMat_Bias_idxs(Gf_tri,no,cp,sIdx,eIdx)
          Gf => fGf(sIdx:eIdx)

          ! retrieve the GGG block
          GGG => val(Gf_tri,cp,n)

          ! We need only do the product in the closest
          ! regions (we don't have information anywhere else)
          call zgemm('N','N', sNc, sN, no, z1, &
               Gf, sNc, oW, no, z0, GGG, sNc)
          
       end do
       
    end do
    
    call timer('GFGGF',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS GFGammaGF' )
#endif

  end subroutine GF_Gamma_GF

  subroutine GFGGF_needed_worksize(N_tri_part, tri_parts, &
       N_Elec, Elecs, worksize)
    use m_ts_electype
    integer, intent(in) :: N_tri_part
    integer, intent(in) :: tri_parts(N_tri_part)
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(out) :: worksize

    integer :: n, idx, sIdx, eIdx, tn, np

    integer :: sN, sNc
    integer :: cp, scp, ecp

    logical :: hasEl

    integer :: no_max, max_n

    ! calculate the maximum electrode size
    no_max = maxval(TotUsedOrbs(Elecs))
    
    ! find at which point they will cross...
    tn = tri_parts(N_tri_part)**2
    do n = 1 , N_tri_part - 1
       tn = tn + tri_parts(n)*( tri_parts(n) + 2 * tri_parts(n+1) )
    end do

    tn = tn + 1

    find: do n = N_tri_part , 1 , - 1

       sN = tri_parts(n)

       scp = max(1 ,n-1)
       ecp = min(N_tri_part,n+1)
       
       do cp = ecp , scp , -1

          ! Update the index of which we will update last
          sNc = tri_parts(cp)
          tn = tn - sNc * sN

          call TriMat_Bias_idxs(N_tri_part, tri_parts, &
               no_max,cp,sIdx,eIdx)
          if ( tn <= eIdx ) then
             max_n = max(cp,n)
             exit find
          end if
          
       end do

    end do find
    
    ! first calculate the naive worksize
    sN = 0
    do n = 1 , max_n
       sN = max(sN,tri_parts(n))
    end do
    sNc = sN
    do n = max_n + 1 , N_tri_part
       sN = max(sN,tri_parts(n))
    end do

    ! Update the index
    call TriMat_Bias_idxs(N_tri_part, tri_parts, &
         no_max,max_n,sIdx,eIdx)

    ! the worksize must be the
    ! maximum matrix created by the zgemm-routines
    ! plus the size of the remaining elements
    tn   = no_max * sN
    sIdx = no_max * sNc + eIdx
    worksize = max(tn, sIdx)

    ! Check whether we can utilize a diagonal
    ! block from the tri-diagonal matrix which is
    ! not entering any of the electrodes.
    ! If able, we can reduce the memory needed!
    sNc = 0
    idx = 0
    do n = 1 , N_tri_part
       hasEl = .false.
       sN = tri_parts(n)
       do cp = sNc + 1 , sNc + sN
          if ( any(OrbInElec(Elecs,cp)) ) then
             hasEl = .true.
             exit
          end if
       end do

       if ( .not. hasEl ) then
          if ( worksize <= sN ** 2 ) then
             idx = n
          end if
       end if
       ! update the next columns
       sNc = sNc + sN
    end do

    ! tell transiesta to point to this diagonal
    ! block in the tri-diagonal matrix instead of allocating
    if ( idx > 0 ) worksize = -idx

  contains
    
    ! We will partition the system by:
    ! 1. nrows_g(tri,1) x no
    ! 2. nrows_g(tri,2) x no
    ! ...
    subroutine TriMat_Bias_idxs(N_tri_part,tri_parts,no,p,sIdx,eIdx)
      integer, intent(in) :: N_tri_part
      integer, intent(in) :: tri_parts(N_tri_part)
      ! no is the number of orbitals we wish to take out
      ! p is the part that we wish to point to
      integer, intent(in) :: no, p
      integer, intent(out) :: sIdx, eIdx
      integer :: cum
      
      ! the size of the partition
      cum = 0
      ! we are requesting the first column,
       ! hence we order the matrix in from the
       ! beginning...
      do eIdx = 1 , p - 1
         cum = cum + tri_parts(eIdx)
      end do
      ! This is the number of elements already occupied
      sIdx = no * cum + 1
      eIdx = sIdx + no * tri_parts(p) - 1
      
    end subroutine TriMat_Bias_idxs

  end subroutine GFGGF_needed_worksize

  subroutine ts_needed_mem(N_tri_part, tri_parts, worksize)
    use m_ts_electype
    use m_ts_options, only : N_Elec, Elecs, IsVolt
    integer, intent(in) :: N_tri_part
    integer, intent(in) :: tri_parts(N_tri_part)
    integer, intent(out) :: worksize

    integer :: n

    ! find at which point they will cross...
    worksize = tri_parts(N_tri_part)**2
    do n = 1 , N_tri_part - 1
       worksize = worksize + &
            tri_parts(n)*( tri_parts(n) + 2 * tri_parts(n+1) )
    end do
    ! for the two arrays
    worksize = worksize * 2

    if ( IsVolt ) then
       call GFGGF_needed_worksize(N_tri_part, tri_parts, &
            N_Elec, Elecs, n)
       worksize = max(0,n) + worksize
    end if
    
  end subroutine ts_needed_mem

end module m_ts_tri_scat
