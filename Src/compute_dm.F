      MODULE m_compute_dm
      private
      public :: compute_dm

      CONTAINS

      subroutine compute_dm( iscf )
      USE SIESTA_TODO
      USE siesta_options
      use sparse_matrices
      implicit none

      integer :: iscf

      external :: automatic_cell,
     .  bonds, cgvc, cgvc_zmatrix, fixed,
     .  dhscf, diagon, dnaefs, extrapol, initatom,
     .  iodm, iozm,
     .  kgrid, kgridinit, kinefsm, mulliken, naefs, neighb,
     .  pulayx, 
     .  reinit, shaper, spnvec, 
     .  timer, xijorb, memory,
     .  ioeig, iofa, iokp, iomd, prversion, eggbox


!------------------------------------------------------------------------- BEGIN

! Solve eigenvalue problem .........................................
          if (isolve .eq. 0) then
            call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     .                  numh, listhptr, listh, numh, listhptr, listh, 
     .                  H, S, qtot, fixspin, qs, temp, e1, e2,
     .                  gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     .                  eo, qo, Dscf, Escf, ef, efs, Entrop, no_u,
     .                  occtol, iscf, neigwanted)
            Ecorrec = 0.0_dp
!
          elseif (isolve .eq. 1) then
            if (.not. gamma) call die("Cannot do O(N) with k-points.")
            call ordern(usesavelwf,ioptlwf,na_s,no_s,no_l,lasto,isa,qa,
     .                  rcoor,rmaxh,ucell,xa,iscf,istp,ncgmax,etol,eta,
     .                  qtot,maxnh,numh,listhptr,listh,H,S,
     .                  chebef,noeta,rcoorcp,beta,pmax,Dscf,Escf,
     .                  Ecorrec,nspin,qs)
            Entrop = 0.0_dp
          else
            call die('siesta: ERROR: wrong solution method')
          endif

! Harris-functional energy ............................................
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
     .                     ( Dscf(io,ispin) - Dold(io,ispin) )
            enddo
          enddo
#ifdef MPI
! Global reduction of DEharr
          call globalize_sum(DEharr,buffer1)
          DEharr = buffer1
#endif
! ..................

! Print populations at each SCF step if requested before mixing ......

          if (muldeb) then
             write (6,"(/a)")
     .             'siesta: Mulliken populations before mixing'
             call mulliken( mullipop, nspin, na_u, no_u, maxnh,
     .                      numh, listhptr, listh, S, Dscf, isa,
     .                      lasto, iaorb, iphorb )
          endif
!--------------------------------------------------------------------------- END
      END subroutine compute_dm
      END MODULE m_compute_dm
