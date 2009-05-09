      MODULE m_compute_dm
      private
      public :: compute_dm

      CONTAINS

      subroutine compute_dm( iscf )
      use precision
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use atomlist, only: qa, lasto, iphorb, iaorb, no_u, no_s, indxuo, 
     .                    qtot, no_l
      use sys, only: die
      use Kpoint_grid
      use m_energies
      use m_rmaxh
      use m_eo
      use m_gamma
      use parallel, only: IONode
#ifdef MPI
      use m_mpi_utils,     only:  globalize_sum
#endif
#ifdef CDF
      use iodmhs_netcdf, only: write_dmh_netcdf
#endif
      use m_ordern,        only : ordern
      use m_steps,         only : istp
      use m_spin,          only : nspin, qs, efs
      use domain_decom,    only : gatherGlobalMatrix
      use jacobi_davidson, only : jd_diag
      use modi_lanczos_m
      implicit none

      integer, intent(in) :: iscf

      integer :: io, ispin
      real(dp):: const  ! Temporary factor used in calculation of
      ! Harris forces 
      ! e1>e2 To signal that we do not want DOS weights
      real(dp), parameter  :: e1 = 1.0_dp, e2 = -1.0_dp
#ifdef MPI
      real(dp):: buffer1
#endif
!------------------------------------------------------------------------- BEGIN
#if defined(MPI_TRACE) && defined(MPI)
      call mpitrace_eventandcounters( 1000, 6 )
#endif
! Solve eigenvalue problem 
      if (isolve .eq. 0) then
        call write_debug( '    PRE diagon' )
        call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     .              numh, listhptr, listh, numh, listhptr, listh, 
     .              H, S, qtot, fixspin, qs, temp, e1, e2,
     .              gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     .              eo, qo, Dscf, Escf, ef, efs, Entrop, no_u,
     .              occtol, iscf, neigwanted)
        Ecorrec = 0.0_dp
!
      elseif (isolve .eq. 1) then
        if (.not. gamma) call die("Cannot do O(N) with k-points.")
        call write_debug( '    PRE ordern' )
            call ordern(usesavelwf,ioptlwf,na_u,no_u,no_l,lasto,iaorb,
     .                  isa,qa,rcoor,rmaxh,ucell,xa,iscf,istp,ncgmax,
     .                  etol,eta,qtot,maxnh,numh,listhptr,listh,H,S,
     .              chebef,noeta,rcoorcp,beta,pmax,Dscf,Escf,
     .              Ecorrec,nspin,qs)
        Entrop = 0.0_dp
      elseif (isolve .eq. 2) then
        call write_debug( '    PRE Jacobi' )
        call gatherGlobalMatrix( no_s, no_l, maxnh, S, 'S.matrix' )
        call gatherGlobalMatrix( no_s, no_l, maxnh, H(1:maxnh,1),
     &                           'H.matrix' )

        call jd_diag( no_s, no_l, nspin, neigwanted, fixspin, temp,
     &                qtot, qs, occtol, numh, listh, H, S,
     &                ef, efs, Entrop, Dscf, Escf )

      elseif (isolve .eq. 3) then
        call write_debug( '    PRE ITERATIVE' )
        call modi_lanczos( no_s, no_l, nspin, fixspin, temp, qtot,
     &                     qs, occtol, numh, listh, H, S, neigwanted,
     &                     ef, efs, Entrop, Dscf, Escf )
        call write_debug( '    POS ITERATIVE' )
        call MPI_Finalize( io )
        STOP
      else
        call die('siesta: ERROR: wrong solution method')
      endif

! Harris-functional energy 
      DEharr = 0.0_dp
      do ispin = 1,nspin
! const factor takes into account that there are two nondiagonal
! elements in non-collinear spin density matrix, stored as
! ispin=1 => D11; ispin=2 => D22, ispin=3 => Real(D12);
! ispin=4 => Imag(D12)
        const = 1._dp
        if (ispin .gt. 2) const = 2._dp
        do io = 1,nh
          DEharr = DEharr + H(io,ispin) * const * 
     .                 ( Dscf(io,ispin) - Dold(io,ispin) )
        enddo
      enddo
#ifdef MPI
! Global reduction of DEharr
      call globalize_sum(DEharr,buffer1)
      DEharr = buffer1
#endif
! ..................
#ifdef CDF
          if ( writedmhs_cdf_history) then
             call write_dmh_netcdf( no_l, nh, nspin,
     $                              Dold, H, Dscf)
          else if ( writedmhs_cdf) then
             call write_dmh_netcdf( no_l, nh, nspin,
     $                              Dold, H, Dscf, overwrite=.true.)
          endif
#endif

! Print populations at each SCF step if requested before mixing ......

      if (muldeb) then
         write (6,"(/a)")
     .         'siesta: Mulliken populations before mixing'
         call mulliken( mullipop, nspin, na_u, no_u, maxnh,
     .                  numh, listhptr, listh, S, Dscf, isa,
     .                  lasto, iaorb, iphorb )
      endif
#if defined(MPI_TRACE) && defined(MPI)
      call mpitrace_eventandcounters( 1000, 0 )
#endif
!--------------------------------------------------------------------------- END
      END subroutine compute_dm
      END MODULE m_compute_dm
