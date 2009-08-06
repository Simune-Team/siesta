      module m_siesta_analysis
      use write_subs

      private
      public :: siesta_analysis

      CONTAINS

      subroutine siesta_analysis( relaxd )
      USE band, only: nbk, bk, maxbk, bands
      USE writewave, only: nwk, wfk, wwave
      USE m_ksvinit, only: nkpol, kpol, wgthpol
      use m_ksv
      USE m_projected_DOS, only: projected_DOS
      USE m_local_DOS, only: local_DOS
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use m_dhscf, only: dhscf
      use atomlist, only: indxuo, iaorb, lastkb, lasto, datm, no_l,
     &                    iphkb, no_u, no_s, iza, iphorb, rmaxo, indxua
      use fdf, only : fdf_boolean
      use writewave,    only: wwave
      use siesta_cml
      use files,        only : slabel
      use zmatrix,      only: lUseZmatrix, write_zmatrix
      use Kpoint_grid
      use parallel, only: IOnode
      use files,       only : label_length
      use m_energies
      use m_steps
      use m_ntm
      use m_spin,         only: nspin
      use m_dipol
      use m_eo
      use m_forces,          only: fa
      use m_gamma 
      use alloc,             only: re_alloc, de_alloc
      use basis_enthalpy,    only: write_basis_enthalpy
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif

      implicit none

      character(len=label_length+5), external :: paste

      logical :: relaxd

      character(len=label_length+5) :: fildrh
      character(len=label_length+5) :: filepsch
      character(len=label_length+5) :: filetoch
      character(len=label_length+5) :: filevh
      character(len=label_length+5) :: filevt
      character(len=label_length+5) :: filrho

      real(dp)          :: dummy_str(3,3), dummy_strl(3,3)
      real(dp)          :: qspin(4)         ! Local

      real(dp)          :: polxyz(3,nspin)  ! Autom., small
      real(dp)          :: polR(3,nspin)    ! Autom., small
      real(dp)          :: qaux
      real(dp)          :: stot             ! Total spin polar.
      real(dp)          :: svec(3)          ! Total spin vector
      real(dp), pointer :: ebk(:,:,:)       ! Band energies

      integer           :: j, ix, ind, ik, io, ispin
      integer           :: maxwk   ! Number of k-points/wavefunction
      real(dp)          :: g2max
#ifdef MPI
      real(dp)          :: qtmp(4)
#endif

!------------------------------------------------------------------------- BEGIN
      final = .true.
!     We want xalast to equal xa for coordinate relaxation only.
!     Possible ERROR !!!!!!!!!!!!!!!!!!!!!!!!
!     xlast stores the previos atom possition and it is needed further
!     if (idyn==0) xalast(1:3,1:na_s)=xa(1:3,1:na_s)
      if (cml_p) then
        call cmlStartModule(xf=mainXML, title='Finalization')
      endif

      if (IOnode) then
!       Print atomic coordinates (and also unit cell for ParrRah.)
        if (nmove .ne. 0) then
          if (relaxd) 
     &      call outcoor(ucell, xa, na_u, 'Relaxed', .true. )
          if (.not.relaxd) 
     &      call outcoor(ucell, xa, na_u,
     &                 'Final (unrelaxed)', .true. )
        endif
        call siesta_write_positions()
        if (lUseZmatrix) call write_Zmatrix
        if ( varcel .or. (idyn.eq.8)) call outcell(ucell)

!       Print coordinates in xmol format in a separate file
        if (fdf_boolean('WriteCoorXmol',.false.)) 
     &     call coxmol(iza, xa, na_u )

!       Print coordinates in cerius format in a separate file
        if (fdf_boolean('WriteCoorCerius',.false.))
     &     call coceri(iza, xa, ucell, na_u, sname )

!       Find interatomic distances (output in file BONDS_FINAL)
        call bonds( ucell, na_u, isa, xa,
     &       rmax_bonds, trim(slabel) // ".BONDS_FINAL" )

       endif ! IONode

!     Find and print wavefunctions at selected k-points
      if (nwk.gt.0) then
        call wwave( no_s, nspin, nspin, no_u, no_l, maxnh, maxwk,
     &              numh, listhptr, listh, H, S, Ef, xijo, indxuo,
     &              nwk, wfk, no_u, gamma, occtol )
      endif

!     Compute bands
      nullify( ebk )
      call re_alloc( ebk, 1, no_u, 1, nspin, 1, maxbk,
     &               'ebk', 'siesta_analysis' )

      if (nbk.gt.0) then
#ifdef DEBUG
        call write_debug( '  PRE compute bands' )
#endif
        call bands( no_s, nspin, nspin, no_u, maxnh, maxbk,
     &              numh, listhptr, listh, H, S, Ef, xijo, indxuo,
     &              .true., nbk, bk, ebk, no_u, occtol )
        if (IOnode) then
          if ( writbk ) then
            write(6,'(/,a,/,a4,a12)')
     &       'siesta: Band k vectors (Bohr**-1):', 'ik', 'k'
            do ik = 1,nbk
              write(6,'(i4,3f12.6)') ik, (bk(ix,ik),ix=1,3)
            enddo
          endif
        
          if ( writb ) then
            write(6,'(/,a,/,a4,a3,a7)')
     &       'siesta: Band energies (eV):', 'ik', 'is', 'eps'
            do ispin = 1,min(nspin,2)
              do ik = 1,nbk
                write(6,'(i4,i3,10f7.2)')
     &            ik, ispin, (ebk(io,ispin,ik)/eV,io=1,min(10,no_u))
                if (no_u.gt.10) write(6,'(7x,10f7.2)')
     &              (ebk(io,ispin,ik)/eV,io=11,no_u)
              enddo
            enddo
          endif
        endif
      endif

!     Print eigenvalues
      if (IOnode .and. writeig) then
        if (isolve.eq.0 .and. no_l.lt.1000) then
          if (nspin .le. 2) then
            write(6,'(/,a,/,a4,a3,a7)')
     &       'siesta: Eigenvalues (eV):', 'ik', 'is', 'eps'
            do ik = 1,nkpnt
              do ispin = 1,nspin
                write(6,'(i4,i3,10f7.2)')
     &            ik,ispin,(eo(io,ispin,ik)/eV,io=1,min(10,neigwanted))
                if (no_u.gt.10) write(6,'(7x,10f7.2)')
     &            (eo(io,ispin,ik)/eV,io=11,neigwanted)
              enddo
            enddo
          else
            write(6,'(/,a)') 'siesta: Eigenvalues (eV):'
            do ik = 1,nkpnt
              write(6,'(a,i6)') 'ik =', ik
              write(6,'(10f7.2)')
     &          ((eo(io,ispin,ik)/eV,io=1,neigwanted),ispin=1,2)
            enddo
          endif
          write(6,'(a,f15.6,a)') 'siesta: Fermi energy =', ef/eV, ' eV'
        endif
      endif

      if ((isolve.eq.0.or.isolve.eq.3).and.IOnode)
     &     call ioeig(eo,ef,neigwanted,nspin,nkpnt,no_u,nspin,maxk,
     &                kpoint, kweight)

      call projected_DOS()

!     Print program's energy decomposition and final forces
      if (IOnode) then
        call siesta_write_energies( first=.false., iscf=0, dDmax=0.0_dp)
        call siesta_write_forces()
        call siesta_write_stress_pressure()
        call write_basis_enthalpy(FreeE,FreeEHarris)
      endif

!     Print spin polarization
      if (nspin .ge. 2) then
        do ispin = 1,nspin
          qspin(ispin) = 0.0_dp
          do io = 1,no_l
            do j = 1,numh(io)
              ind = listhptr(io)+j
              qspin(ispin) = qspin(ispin) + Dscf(ind,ispin) * S(ind)
            enddo
          enddo
        enddo

#ifdef MPI
!       Global reduction of spin components
        call globalize_sum(qspin(1:nspin),qtmp(1:nspin))
        qspin(1:nspin) = qtmp(1:nspin)
#endif
        if (nspin .eq. 2) then
          if (IOnode) then
            write(6,'(/,a,f12.6)')
     &       'siesta: Total spin polarization (Qup-Qdown) =', 
     &       qspin(1) - qspin(2)
          endif
          if (cml_p) call cmlAddProperty(xf=mainXML,
     &         property=qspin(1)-qspin(2), dictref='siesta:qspin')
        elseif (nspin .eq. 4) then
          call spnvec( nspin, qspin, qaux, stot, svec )
          if (IOnode) then
            write(6,'(/,a,f12.6)')
     &       'siesta: Total spin polarization (Qup-Qdown) =', stot
            write(6,'(a,3f12.6)') 'siesta: Spin vector =', svec
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, property=stot,
     &             dictref='siesta:stot')
              call cmlAddProperty(xf=mainXML, property=svec,
     &             dictref='siesta:svec')
            endif !cml_p
          endif
        endif
      endif

!     Print electric dipole
      if (shape .ne. 'bulk') then
        if (IOnode) then
          write(6,'(/,a,3f12.6)')
     &      'siesta: Electric dipole (a.u.)  =', dipol
          write(6,'(a,3f12.6)')
     &      'siesta: Electric dipole (Debye) =', 
     &      (dipol(ix)/Debye,ix=1,3)
        endif
        if (cml_p) then
          call cmlAddProperty(xf=mainXML, property=dipol,
     &         title='Electric dipole', dictref='siesta:dipol',
     &         units='siestaUnits:atomic')
        endif !cml_p
      endif

!     Calculation of the bulk polarization using the Berry phase
!     formulas by King-Smith and Vanderbilt
!     Attention H0 is used as an auxiliary array
      if (nkpol.gt.0 .and. .not.bornz) then
        call KSV_pol(na_u, na_s, xalast, rmaxo, scell_last, ucell_last,
     &               no_u, no_l, no_s, nspin, qspin, 
     &               maxnh, nkpol, numh, listhptr, listh, 
     &               H, S, H0, xijo, indxuo, isa, iphorb, 
     &               iaorb, lasto, shape,
     &               nkpol,kpol,wgthpol, polR, polxyz ) 
      endif

!     Calculation of the optical conductivity
!     Attention H0, Eold, Dold are used as auxiliary arrays
#ifdef DEBUG
      call write_debug( '  PRE optical' )
#endif
      call optical(na_u, na_s, xalast, scell_last, ucell_last,
     &             no_u, no_l, no_s, nspin, qspin,
     &             maxnh, numh, listhptr, listh, H, S, H0,
     &             Eold(1,1), Dold(1,1), xijo, indxuo, ebk, ef, temp,
     &             isa, iphorb, iphKB, lasto, lastkb, shape )

      call de_alloc( ebk, 'ebk', 'siesta_analysis' )
!...................................

!     Save electron density and potential
      if (savrho .or. savdrh .or. savevh .or. savevt .or.
     &    savepsch .or. savetoch ) then
        filrho = ' '
        fildrh = ' '
        filevh = ' '
        filevt = ' '
        filepsch = ' '
        filetoch = ' '
        if (savrho) filrho = paste( slabel, '.RHO' )
        if (savdrh) fildrh = paste( slabel, '.DRHO' )
        if (savevh) filevh = paste( slabel, '.VH'  )
        if (savevt) filevt = paste( slabel, '.VT'  )
        if (savepsch) filepsch = paste( slabel, '.IOCH'  )
        if (savetoch) filetoch = paste( slabel, '.TOCH'  )
        g2max = g2cut
#ifdef DEBUG
        call write_debug( '  PRE dhscf' )
#endif
        call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     &              no_u, na_u, na_s, isa, xalast, indxua, 
     &              ntm, 0, 0, 0, filrho, fildrh, filevh, filevt,
     &              filepsch, filetoch, 
     &              maxnh, numh, listhptr, listh, Dscf, Datm,
     &              maxnh, H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     &              Exc, Dxc, dipol, dummy_str, fa, dummy_strl )
                    ! next to last argument is dummy here,
                    ! as no forces are calculated
                    ! todo: make all these optional
#ifdef DEBUG
        call write_debug( '  POS dhscf' )
#endif
      endif

C     Find local density of states
#ifdef DEBUG
      call write_debug( '  PRE local_DOS' )
#endif
      call local_DOS()

!--------------------------------------------------------------------------- END
      END subroutine siesta_analysis
      END module m_siesta_analysis
