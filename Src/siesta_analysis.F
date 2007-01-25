      subroutine siesta_analysis( first, relaxd, iscf )
      USE SIESTA_TODO
      implicit none
      logical :: first, relaxd
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
      final = .true.

      ! We want xalast to equal xa for coordinate relaxation only.
      if (idyn==0) xalast(1:3,1:na_s)=xa(1:3,1:na_s)

      if (cml_p) then
        call cmlStartModule(xf=mainXML, title='Finalization')
      endif

      if (IOnode) then
! Print atomic coordinates (and also unit cell for ParrRah.)
        if (nmove .ne. 0) then
          if (relaxd) 
     .      call outcoor(ucell, xa, na_u, 'Relaxed', .true. )
          if (.not.relaxd) 
     .      call outcoor(ucell, xa, na_u,
     .                 'Final (unrelaxed)', .true. )
        endif
        call siesta_write_positions()
        if (lUseZmatrix) call write_Zmatrix
        if ( varcel .or. (idyn.eq.8)) call outcell(ucell)

! Print coordinates in xmol format in a separate file

        if (fdf_boolean('WriteCoorXmol',.false.)) 
     .     call coxmol(iza, xa, na_u )

! Print coordinates in cerius format in a separate file

        if (fdf_boolean('WriteCoorCerius',.false.))
     .     call coceri(iza, xa, ucell, na_u, sname )

! Find interatomic distances (output in file BONDS_FINAL)

        call bonds( ucell, na_u, isa, xa,
     $       rmax_bonds, trim(slabel) // ".BONDS_FINAL" )

       endif ! IONode

! Find and print wavefunctions at selected k-points
      if (nwk.gt.0) then
        call wwave( no_s, nspin, nspin, no_u, no_l, maxnh, maxwk,
     .              numh, listhptr, listh, H, S, Ef, xijo, indxuo,
     .              nwk, wfk, no_u, gamma, occtol )
      endif

! Find and print bands
      if (nbk.gt.0) then
        call bands( no_s, nspin, nspin, no_u, no_l, maxnh, maxbk,
     .              numh, listhptr, listh, H, S, Ef, xijo, indxuo,
     .              .true., nbk, bk, ebk, no_u, occtol )
        if (IOnode) then
          if ( writbk ) then
            write(6,'(/,a,/,a4,a12)')
     .       'siesta: Band k vectors (Bohr**-1):', 'ik', 'k'
            do ik = 1,nbk
              write(6,'(i4,3f12.6)') ik, (bk(ix,ik),ix=1,3)
            enddo
          endif
        
          if ( writb ) then
            write(6,'(/,a,/,a4,a3,a7)')
     .       'siesta: Band energies (eV):', 'ik', 'is', 'eps'
            do ispin = 1,min(nspin,2)
              do ik = 1,nbk
                write(6,'(i4,i3,10f7.2)')
     .            ik, ispin, (ebk(io,ispin,ik)/eV,io=1,min(10,no_u))
                if (no_u.gt.10) write(6,'(7x,10f7.2)')
     .              (ebk(io,ispin,ik)/eV,io=11,no_u)
              enddo
            enddo
          endif
        endif
      endif

! Print eigenvalues
      if (IOnode .and. writeig) then
        if (isolve.eq.0 .and. no_l.lt.1000) then
          if (nspin .le. 2) then
            write(6,'(/,a,/,a4,a3,a7)')
     .       'siesta: Eigenvalues (eV):', 'ik', 'is', 'eps'
            do ik = 1,nkpnt
              do ispin = 1,nspin
                write(6,'(i4,i3,10f7.2)')
     .            ik,ispin,(eo(io,ispin,ik)/eV,io=1,min(10,neigwanted))
                if (no_u.gt.10) write(6,'(7x,10f7.2)')
     .            (eo(io,ispin,ik)/eV,io=11,neigwanted)
              enddo
            enddo
          else
            write(6,'(/,a)') 'siesta: Eigenvalues (eV):'
            do ik = 1,nkpnt
              write(6,'(a,i6)') 'ik =', ik
              write(6,'(10f7.2)')
     .          ((eo(io,ispin,ik)/eV,io=1,neigwanted),ispin=1,2)
            enddo
          endif
          write(6,'(a,f15.6,a)') 'siesta: Fermi energy =', ef/eV, ' eV'
        endif
      endif

      if (isolve.eq.0.and.IOnode)
     .     call ioeig(eo,ef,neigwanted,nspin,nkpnt,no_u,nspin,maxk,
     .                kpoint, kweight)

! Compute the projected density of states
      if (IOnode) then
        do_pdos = fdf_block('ProjectedDensityOfStates',iu)
        if (isolve.ne.0.and.do_pdos) then
          write(6,*)
     .         'siesta: ERROR: PDOS implemented only with diagon'
          do_pdos = .false.
        endif
      endif
      call broadcast(do_pdos)

      if (do_pdos) then
! Find the desired energy range
        if (IOnode) then
          read(iu,'(a)') line
          p=>digest(line)
          if (nvalues(p).lt.3 .or. nnames(p).ne.1)
     $          call die("Wrong format in PDOS block")
          factor = fdf_convfac( names(p,1), 'Ry' )
          e1 = values(p,1) * factor
          e2 = values(p,2) * factor
          sigma = values(p,3) * factor
          nhist = integers(p,1)
          write(6,'(a)') 'siesta: PDOS info: '
          write(6,'(a,3(f8.2,a),2x,i5)')
     $           'siesta: e1, e2, sigma, nhist: ',
     $           e1/eV,' eV',e2/eV,' eV',sigma/eV,' eV', nhist
        endif

        call broadcast(e1)
        call broadcast(e2)
        call broadcast(sigma)
        call broadcast(nhist)

        call pdos( no_s, nspin, nspin, no_l, maxnh,
     .             no_u, numh, listhptr, listh, H, S,
     .             e1, e2, sigma, nhist,
     .             gamma, xijo, indxuo, nkpnt, kpoint, kweight, eo,
     .             no_u)

      endif                     ! PDOS calc (do_pdos)

! Print program's energy decomposition and final forces
      if (IOnode) then
        call siesta_write_energies( first, iscf )
        call siesta_write_forces()
        call siesta_write_stress_pressure()
      endif

! Print spin polarization
      if (nspin .ge. 2) then
        do ispin = 1,nspin
          qspin(ispin) = 0.0_dp
          do io = 1,no_l
            do j = 1,numh(io)
              ind = listhptr(io)+j
              jo = listh(ind)
              qspin(ispin) = qspin(ispin) + Dscf(ind,ispin) * S(ind)
            enddo
          enddo
        enddo

#ifdef MPI
! Global reduction of spin components
      call globalize_sum(qspin(1:nspin),qtmp(1:nspin))
      qspin(1:nspin) = qtmp(1:nspin)
#endif
        if (nspin .eq. 2) then
          if (IOnode) then
            write(6,'(/,a,f12.6)')
     .       'siesta: Total spin polarization (Qup-Qdown) =', 
     .       qspin(1) - qspin(2)
          endif
          if (cml_p) call cmlAddProperty(xf=mainXML,
     .         property=qspin(1)-qspin(2), dictref='siesta:qspin')
        elseif (nspin .eq. 4) then
          call spnvec( nspin, qspin, qaux, stot, svec )
          if (IOnode) then
            write(6,'(/,a,f12.6)')
     .       'siesta: Total spin polarization (Qup-Qdown) =', stot
            write(6,'(a,3f12.6)') 'siesta: Spin vector =', svec
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, property=stot,
     .             dictref='siesta:stot')
              call cmlAddProperty(xf=mainXML, property=svec,
     .             dictref='siesta:svec')
            endif !cml_p
          endif
        endif
      endif

! Print electric dipole
      if (shape .ne. 'bulk') then
        if (IOnode) then
          write(6,'(/,a,3f12.6)')
     .      'siesta: Electric dipole (a.u.)  =', dipol
          write(6,'(a,3f12.6)')
     .      'siesta: Electric dipole (Debye) =', 
     .      (dipol(ix)/Debye,ix=1,3)
        endif
        if (cml_p) then
          call cmlAddProperty(xf=mainXML, property=dipol,
     .         title='Electric dipole', dictref='siesta:dipol',
     .         units='siestaUnits:atomic')
        endif !cml_p
      endif

! Calculation of the bulk polarization using the Berry phase
! formulas by King-Smith and Vanderbilt
! Attention H0 is used as an auxiliary array
      if (nkpol.gt.0 .and. .not.bornz) then
        call KSV_pol(na_u, na_s, xalast, rmaxo, scell, ucell,
     .               no_u, no_l, no_s, nspin, qspin, maxna, 
     .               maxnh, nkpol, numh, listhptr, listh, 
     .               H, S, H0, xijo, indxuo, isa, iphorb, 
     .               iaorb, lasto, jna, xij, r2ij,shape,
     .               nkpol,kpol,wgthpol, polR, polxyz ) 
      endif

! Calculation of the optical conductivity
! Attention H0, Eold, Dold are used as auxiliary arrays
      call optical(na_u, na_s, xa, scell, ucell,
     .             no_u, no_l, no_s, nspin, qspin,
     .             maxna, maxnh, numh, listhptr, listh, H, S, H0,
     .             Eold(1,1), Dold(1,1),
     .             xijo, indxuo, indxua, ebk, ef, temp,
     .             isa, iphorb, iphKB, iaorb, lasto, lastkb,
     .             jna, xij, r2ij, shape )

!...................................

! Save electron density and potential
      call fdf_global_get(savrho,'SaveRho',
     $                      dumpcharge .or. .false.)
      call fdf_global_get(savdrh,'SaveDeltaRho',       .false.)
      call fdf_global_get(savevh,'SaveElectrostaticPotential',
     $                      .false.)
      call fdf_global_get(savevt,'SaveTotalPotential', .false.)
      call fdf_global_get(savepsch,'SaveIonicCharge',  .false.)
      call fdf_global_get(savetoch,'SaveTotalCharge',  .false.)

      if (savrho .or. savdrh .or. savevh .or. savevt .or.
     .    savepsch .or. savetoch ) then
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
        call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .              no_u, na_u, na_s, isa, xa, indxua, 
     .              ucell, mscell, g2max, ntm,
     .              0, 0, 0, filrho, fildrh, filevh, filevt,
     .              filepsch, filetoch, 
     .              maxnh, numh, listhptr, listh, Dscf, Datm,
     .              maxnh, numh, listhptr, listh, H,
     .              Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, Exc, Dxc,
     .              dipol, fa, stress, fal, stressl )
      endif

! Find local density of states
      if (IOnode) then
        genlogic = fdf_block('LocalDensityOfStates',iu)
      endif
      call broadcast(genlogic)

      if ( genlogic ) then

! Find the desired energy range
        if (IOnode) then
          read(iu,'(a)') line
          p=>digest(line)
          if (.not. match(p,"vvn"))
     .       call die("Wrong format in LocalDensityofStates")
          factor = fdf_convfac( names(p,1), 'Ry' )
          e1 = values(p,1)*factor
          e2 = values(p,2)*factor
          call destroy(p)
        endif
        call broadcast(e1)
        call broadcast(e2)

!       Find the density matrix for states between e1 and e2
        if (isolve .eq. 0) then
          call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     .                numh, listhptr, listh, numh, listhptr, listh, 
     .                H, S, qtot, fixspin, qs, temp, e1, e2,
     .                gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     .                eo, qo, Dscf, Escf, ef, efs, Entrop, no_u,
     .                occtol, iscf, neigwanted)

!       Find the LDOS in the real space mesh
          filrho = paste( slabel, '.LDOS' )
          g2max = g2cut
          call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .              no_u, na_u, na_s, isa, xa, indxua, 
     .              ucell, mscell, g2max, ntm,
     .              0, 0, 0, filrho, ' ', ' ', ' ', ' ', ' ',
     .              maxnh, numh, listhptr, listh, Dscf, Datm,
     .              maxnh, numh, listhptr, listh, H,
     .              Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, Exc, Dxc,
     .              dipol, fa, stress, fal, stressl )
        else
          if (IOnode)  write(6,*)
     .       'siesta: ERROR: LDOS implemented only with diagon'
        endif

      endif ! genlogic
!--------------------------------------------------------------------------- END
      END subroutine siesta_analysis
