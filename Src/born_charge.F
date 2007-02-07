      module m_born_charge
      use siesta_todo

      private
      public :: born_charge

      CONTAINS
      
      subroutine born_charge()
!
!     Computes Born-effective charges by finite differences,
!     using the "force constant calculation" mode of operation.
!     It is called only for those steps (even) in which the dx
!     is positive.
!
      if (nkpol.lt.1) then
         if (IOnode) write(6,'(/,a,f12.6)')
     .        'siesta: specify polarization grid for BC calculation'
         if (IOnode) write(6,'(a,f12.6)')
     .        'siesta: The Born charge matrix will not be calculated'
         RETURN
      endif
      if (IOnode) write(6,'(/,a,f12.6)')
     .     'siesta: Calculating polarization. '

!     Find total population of spin up and down
      if (nspin .ge. 2) then
         do ispin = 1,nspin
            qspin(ispin) = 0.0_dp
            do io = 1,no_l
               do j = 1,numh(io)
                  ind = listhptr(io) + j
                  qspin(ispin) = qspin(ispin)
     .                 + Dscf(ind,ispin)*S(ind)
               enddo
            enddo
         enddo
#ifdef MPI
!     Global reduction of spin components
         call globalize_sum(qspin(1:nspin),qtmp(1:nspin))
         qspin(1:nspin) = qtmp(1:nspin)
#endif
      endif
      if (nkpol.gt.0) then
         ! Note use of xa instead of xalast now
         call KSV_pol(na_u, na_s, xa, rmaxo, scell, ucell,
     .        no_u, no_l, no_s, nspin, qspin, maxna,
     .        maxnh, nkpol, numh, listhptr, listh,
     .        H, S, H0, xijo, indxuo, isa, iphorb,
     .        iaorb, lasto, jna, xij, r2ij,shape,
     .        nkpol,kpol,wgthpol, polR, polxyz)
      endif
      if (nkpol.gt.0.and.IOnode) then
         call obc( polxyz, polR, ucell, dx, nspin, node )
      endif

      end subroutine born_charge

      end module m_born_charge
