      MODULE m_mixer
      private
      public :: mixer

      CONTAINS

      subroutine mixer( iscf )
      USE SIESTA_TODO
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use atomlist, only: iaorb, iphorb, lasto, no_s, no_u
      use m_broyden_mixing
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
! Mix input and output energy-density and density matrices ............
! Following line for using and saving the density matrix without mix ..
          if (wmix.ne.0._dp) then
! Pulay or Broyden mixing
            mmix  = mix
            iiscf = iscf
            if (maxsav .le. 0) then
              iiscf = 1
              if (iscf .ne. 1) mmix = .true.
            endif
            if (broyden_maxit == 0) then
               call pulayx( pulfile, iiscf, mmix, no_l, no_s, maxnh,
     .                   numh, listhptr, nspin, maxsav, wmix, nkick, 
     .                   wmixkick, auxpul(1,1), auxpul(1,2), nauxpul,
     .                   Dscf, Dold, dDmax)
            else
               call broyden_mixing(iscf, mix, no_l, maxnh,
     .                 numh(1:no_l), listhptr(1:no_l), nspin,
     $                 wmix, nkick, wmixkick, Dscf, Dold, dDmax)
            endif
          endif

! Ensure that dDmax is the same on all nodes for convergence test/output
#ifdef MPI
          call globalize_max(dDmax,buffer1)
          dDmax = buffer1
#endif
! ...................

! Print populations at each SCF step, if requested, after mixing ......

          if (muldeb) then 
             write (6,"(/a)")
     .             'siesta: Mulliken populations after mixing'
             call mulliken( mullipop, nspin, na_u, no_u, maxnh,
     .                      numh, listhptr, listh, S, Dscf, isa,
     .                      lasto, iaorb, iphorb )
          endif
! ..................

! Save density matrix on disk, after mixing ...........................
          if (writedm) then
            if ((idyn .eq. 6) .or. (idyn .eq. 7)) then
              if (istp.eq.1)
     .        call iodm( 'write', maxnh, no_l, nspin,
     .                   numh, listhptr, listh, Dscf, found )
            else
              call iodm( 'write', maxnh, no_l, nspin,
     .                   numh, listhptr, listh, Dscf, found )
            endif
          endif
!--------------------------------------------------------------------------- END
      END subroutine mixer

      End MODULE m_mixer


