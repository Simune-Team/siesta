      Module write_subs
      use siesta_todo
      use m_stress
      use m_steps

      private
      public :: siesta_write_forces, siesta_write_stress_pressure,
     &          siesta_write_energies, siesta_write_positions
      CONTAINS
      
      subroutine siesta_write_forces()
      USE siesta_options
      use siesta_geom
      use siesta_cml
      use units
      implicit none
      ! Almost the same forces output whether during simulation
      ! or at the end. Unfortunately not quite, therefore slightly
      ! tortuous logic below. If we are content to change format
      ! of output file slightly, this can be simplified.
      if (.not.final) then
        ! print forces to xml every step. 
        ! output forces to stdout depending on writef
        if (cml_p) then
          call cmlStartPropertyList(mainXML, title='Forces')
          call cmlAddProperty(xf=mainXML, property=fa*Ang/eV,
     .         dictref='siesta:forces', units='siestaUnits:evpa')
          call cmlAddProperty(xf=mainXML, property=ftot,
     .         dictref='siesta:ftot')
          call cmlAddProperty(xf=mainXML, property=fmax, 
     .         dictref='siesta:fmax')
          call cmlAddProperty(xf=mainXML, property=fres,
     .         dictref='siesta:fres')
          call cmlAddProperty(xf=mainXML, property=cfmax, 
     .         dictref='siesta:cfmax')
          call cmlEndPropertyList(mainXML)
        endif
        write(6,'(/,a)') 'siesta: Atomic forces (eV/Ang):'
        if (writef) then
          write(6,'(i6,3f12.6)')(ia,(fa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
        else
          call iofa( na_u, fa )
        endif
        write(6,'(40("-"),/,a6,3f12.6)') 'Tot',(ftot(ix)*Ang/eV,ix=1,3)
        write(6,'(40("-"),/,a6, f12.6)') 'Max',fmax*Ang/eV
        write(6,'(a6,f12.6,a)')'Res',fres*Ang/eV,
     .       '    sqrt( Sum f_i^2 / 3N )'
        write(6,'(40("-"),/,a6, f12.6,a)') 'Max',cfmax*Ang/eV, 
     .       '    constrained'
      else !not final
! In finalization, only print forces if sufficiently large.
        fmax = maxval(abs(fa))
        ftot = sum(fa, dim=2)
        if (fmax .gt. ftol) then
          write(6,'(/,a)') 'siesta: Atomic forces (eV/Ang):'
          write(6,'(a,i6,3f12.6)')
     .         ('siesta: ', ia,(fa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
          write(6,'(a,40("-"),/,a,a6,3f12.6)')
     .         'siesta: ','siesta: ','Tot',(ftot(ix)*Ang/eV,ix=1,3)
          if (cml_p) then
            call cmlStartPropertyList(mainXML, title='Force Summary')
            call cmlAddProperty(xf=mainXML, property=fa*Ang/eV,
     .           dictref='siesta:forces', units='siestaUnits:evpa')
            call cmlAddProperty(xf=mainXML, property=ftot*Ang/eV, 
     .           dictref='siesta:ftot', units='siestaUnits:evpa')
            call cmlEndPropertyList(mainXML)
          endif !cml_p
        endif
        if (Any(cfa /= fa)) then
          fmax = maxval(abs(cfa))
          ftot = sum(cfa, dim=2)
          if (fmax .gt. ftol) then
            write(6,'(/,a)') 'siesta: Constrained forces (eV/Ang):'
            write(6,'(a,i6,3f12.6)')
     .           ('siesta: ',ia,(cfa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
            write(6,'(a,40("-"),/,a,a4,3f12.6)')
     .           'siesta: ','siesta: ','Tot',(ftot(ix)*Ang/eV,ix=1,3)
            if (cml_p) then
              call cmlStartPropertyList(mainXML,
     .             title='Constrained Force Summary')
              call cmlAddProperty(xf=mainXML, property=cfa*Ang/eV, 
     .             dictref='siesta:cforces', units='siestaUnits:evpa')
              call cmlAddProperty(xf=mainXML, property=ftot*Ang/eV, 
     .             dictref='siesta:cftot', units='siestaUnits:evpa')
              call cmlEndPropertyList(mainXML)
            endif !cml_p
          endif
        endif
      endif !final for forces

      end subroutine siesta_write_forces

      subroutine siesta_write_stress_pressure()

      use parallel, only: IOnode
      USE siesta_options
      use siesta_geom
      use atomlist, only: iza
      use m_iostruct,   only: write_struct
      use siesta_cml
      use units
      use m_energies, only: FreeE

      implicit none

! Stress tensor and pressure:
      

      if (.not.final) then
!
!           Write Voigt components of total stress tensor 
!
            write(6,'(/,a,6f12.2)')
     .           'Stress-tensor-Voigt (kbar):',
     .           (tstress(jx,jx)/kbar,jx=1,3),
     $            tstress(1,2)/kbar,
     $            tstress(2,3)/kbar,
     $            tstress(1,3)/kbar
!
!      Write "target enthalpy" (E + pV, where p is the *target* pressure)
            write(6,"(a,f14.4)") "Target enthalpy (eV/cell)",
     $           (FreeE + tp*volume)/eV
!
!      Write out structural information in "crystallography" format
!
            call write_struct( ucell, na_u, isa, iza, xa )

      ! Output depends on dynamics option
        select case (idyn)
        case(0:5,8)
          if (idyn==0 .and. (.not.varcel)) then
            continue
          else
            write(6,'(/,a,3(/,a,3f12.6))')
     .           'siesta: Stress tensor (static) (eV/Ang**3):',
     .           ('     ',(stress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
            Psol = - ((stress(1,1) + stress(2,2) + stress(3,3))/3.0_dp)
            write(6,'(/,a,f20.8,a)')
     .           'siesta: Pressure (static):', Psol/kBar, '  kBar'
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, property=stress*Ang**3, 
     .             dictref='siesta:stress')
              call cmlAddProperty(xf=mainXML, property=Psol, 
     .             dictref='siesta:psol', title='Pressure (Static)')
            endif !cml_p
!
            write(6,'(/,a,3(/,a,3f12.6))')
     .           'siesta: Stress tensor (total) (eV/Ang**3):',
     .           ('     ',(tstress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
            Psol = - ((tstress(1,1)+tstress(2,2) +tstress(3,3))/3.0_dp)
            write(6,'(/,a,f20.8,a)')
     .           'siesta: Pressure (total):', Psol/kBar, '  kBar'
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, property=tstress*Ang**3, 
     .             dictref='siesta:tstress')
              call cmlAddProperty(xf=mainXML, property=Psol,
     .             dictref='siesta:tpsol', title='Pressure (Total)')
            endif !cml_p
          endif !varcel
        ! Write Force Constant matrix if FC calculation ...
        case(6)
          call ofc(fa,dx,na_u)
        case(7)
!          call phonon_write_forces(fa,na_u,ns,ucell,istep)
           if (IOnode) write(*,*) 'phonon suport deactivated'
        end select !idyn

      else !final

! Print stress tensor unconditionally
        write(6,'(/,a,3(/,a,3f12.6))')
     .       'siesta: Stress tensor (static) (eV/Ang**3):',
     .       ('siesta: ',(stress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
        if (cml_p) then
          call cmlAddProperty(xf=mainXML, property=stress*Ang**3/eV, 
     .         dictref='siesta:stress', units='siestaUnits:eV_Ang__3')
        endif !cml_p

! Print constrained stress tensor if different from unconstrained
        if (Any(cstress /= stress )) then
             write(6,'(/,a,3(/,a,3f12.6))')
     .       'siesta: Constrained stress tensor (static) (eV/Ang**3):',
     .       ('siesta: ',(cstress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
          if (cml_p) then
            call cmlAddProperty(xf=mainXML, property=cstress*Ang**3/eV, 
     .           dictref='siesta:cstress', 
     .           units='siestaUnits:eV_Ang__3')
          endif !cml_p
        endif

! Find pressure
        virial = 0.0_dp
        do ix = 1,3
          fmean = 0.0_dp
          do ia = 1,na_u
            fmean = fmean + fa(ix,ia) / na_u
          enddo
          do ia = 1,na_u
            virial = virial + xalast(ix,ia) * (fa(ix,ia) - fmean)
          enddo
        enddo
        Psol = - (( stress(1,1) + stress(2,2) + stress(3,3) )/3.0_dp)
        Pmol = Psol - virial / volume / 3.0_dp
        write(6,'(/,a,f18.6,a)')
     .       'siesta: Cell volume =', volume/Ang**3, ' Ang**3'
        write(6,'(/,a,/,a,2a20,a,3(/,a,2f20.8,a))')
     .       'siesta: Pressure (static):',
     .       'siesta: ','Solid',        'Molecule',      '  Units',
     .       'siesta: ', Psol,           Pmol,           '  Ry/Bohr**3',
     .       'siesta: ', Psol*Ang**3/eV, Pmol*Ang**3/eV, '  eV/Ang**3',
     .       'siesta: ', Psol/kBar,      Pmol/kBar,      '  kBar'
        if (cml_p) then
          call cmlStartPropertyList(mainXML, title='Final Pressure')
          call cmlAddProperty(xf=mainXML, property=volume/Ang**3, 
     .         title='cell volume', dictref='siesta:cellvol', 
     .         units='siestaUnits:Ang__3')
          call cmlAddProperty(xf=mainXML, property=Psol/kBar, 
     .         title='Pressure of Solid', dictref='siesta:pressSol', 
     .         units='siestaUnits:kbar')
          call cmlAddProperty(xf=mainXML, property=Pmol/kBar,       
     .         title='Pressure of Molecule', dictref='siesta:pressMol', 
     .         units='siestaUnits:kbar')
          call cmlEndPropertyList(mainXML)
        endif !cml_p
        
      endif !final for stress & pressure

      end subroutine siesta_write_stress_pressure
      
      
      subroutine siesta_write_energies( first, iscf, dDmax )
      USE siesta_options
      use siesta_cml
      use units
      use m_energies

      implicit none

      logical :: first
      integer :: iscf
      real(dp), intent(in) :: dDmax     ! Max. change in density matrix elements

      ! Only print out full decomposition at very beginning and end.
      if ((istp==1.and.first).or.final) then
        write(6,'(/,a,/,(a,f17.6))')
     .     'siesta: Program''s energy decomposition (eV):',
     .     'siesta: Eions   =', Eions/eV,
     .     'siesta: Ena     =', Ena/eV,
     .     'siesta: Ekin    =', Ekin/eV,
     .     'siesta: Enl     =', Enl/eV,
     .     'siesta: DEna    =', DEna/eV,
     .     'siesta: DUscf   =', DUscf/eV,
     .     'siesta: DUext   =', DUext/eV,
     .     'siesta: Exc     =', Exc/eV,
     .     'siesta: eta*DQ  =', Ecorrec/eV,
     .     'siesta: Emadel  =', Emad/eV,
     .     'siesta: Emeta   =', Emeta/eV,
     .     'siesta: Emolmec =', Emm/eV,
     .     'siesta: Ekinion =', Ekinion/eV,
     .     'siesta: Eharris =', (Eharrs1+Ekinion)/eV,
     .     'siesta: Etot    =', (Etot+Ekinion)/eV,
     .     'siesta: FreeEng =', (FreeE+Ekinion)/eV
        if (cml_p) then
          call cmlStartPropertyList(mainXML,
     .         title='Energy Decomposition')
          call cmlAddProperty(xf=mainXML, property=Eions/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Eions', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Ena/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ena', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Ekin/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ekin', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Enl/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Enl', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=DEna/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:DEna', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=DUscf/eV, 
     .         units='siestaUnits:eV',
     .         dictref='siesta:DUscf', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=DUext/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:DUext', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Exc/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Exc', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML,property=Ecorrec/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ecorrec', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Emad/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Emad', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Emeta/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Emeta', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Emm/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Emm', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML,property=Ekinion/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ekinion', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=(Eharrs1+Ekinion)/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:EharrsK', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=(Etot+Ekinion)/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:EtotK', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=(FreeE+Ekinion)/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:FreeEK', fmt='(f17.6)')
          call cmlEndPropertyList(mainXML)
        endif
      endif
      ! On all SCF steps, print out the current energy (format depending on type of run)
      if (.not.final) then
        ! Print total energy and density matrix error .........................
        if (cml_p) then
          call cmlStartPropertyList(mainXML, title='SCF Cycle')
          ! Eharrs is always output
          call cmlAddProperty(xf=mainXML, property=Eharrs/eV, 
     .         units="siestaUnits:eV", 
     .         dictRef="siesta:Eharrs", fmt="(f14.7)")
        endif
        ! This chain of if statements determines which properties are output.
        if (harrisfun) then 
          write(6,"(/a,f14.6,/)") 'siesta: Eharris(eV) = ', Eharrs/eV
            ! No need for further cml output
        elseif (isolve==0) then
          if (cml_p) 
     .         call cmlAddProperty(xf=mainXML, property=FreeE/eV, 
     .         units="siestaUnits:eV", 
     .         dictRef="siesta:FreeE",  fmt="(f14.7)")
          if (fixspin) then
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, property=Etot/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:Etot",   fmt="(f14.7)")
              call cmlAddProperty(xf=mainXML, property=FreeE/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:FreeE", fmt="(f14.7)")
              call cmlAddProperty(xf=mainXML, property=dDmax/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:dDmax",  fmt="(f14.7)")
            endif
            if ((iscf .eq. 1).or.muldeb)
     .           write(6,'(/,a12,3a14,a8,a7,a11)')
     .           'siesta: iscf', '   Eharris(eV)', 
     .           '      E_KS(eV)', '   FreeEng(eV)', 
     .           '   dDmax', '  Ef_up', '  Ef_dn(eV)'
            write(6,'(a8,i4,3f14.4,f8.4,2f9.4)')
     .           'siesta: ',iscf, Eharrs/eV, Etot/eV, FreeE/eV, dDmax, 
     .           (Efs(i)/eV,i=1,2)
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, property=Efs(1)/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:Efs",    fmt="(f14.7)")
              call cmlAddProperty(xf=mainXML, property=Efs(2)/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:Efs",    fmt="(f14.7)")
            endif
          else !fixspin
            if ((iscf .eq. 1).or.muldeb)
     .           write(6,'(/,a12,3a14,2a8)')
     .           'siesta: iscf', '   Eharris(eV)', 
     .           '      E_KS(eV)', '   FreeEng(eV)', 
     .           '   dDmax', '  Ef(eV)'
            write(6,'(a8,i4,3f14.4,2f8.4)')
     .           'siesta: ',iscf, Eharrs/eV, Etot/eV, FreeE/eV, 
     .           dDmax, Ef/eV
            if (cml_p) then
              call cmlAddProperty(xf=mainXML,property=Ef/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:Ef", fmt="(f14.7)")
            endif !cml_p
          endif !fixspin
        elseif (isolve==1) then
          write(6,'(/,a15,i4)') 'siesta: iscf = ',iscf
          write(6,'(a14,f15.4,a13,f15.4,a10,f7.4/)') 
     .         'Eharris(eV) = ',Eharrs/eV,
     .         '  E_KS(eV) = ',Etot/eV,'  dDmax = ',dDmax
          if (cml_p) then
            call cmlAddProperty(xf=mainXML, property=Etot/eV, 
     .           units="siestaUnits:eV", 
     .           dictRef="siesta:Etot",   fmt="(f14.7)")
            call cmlAddProperty(xf=mainXML, property=dDmax/eV, 
     .           units="siestaUnits:eV", 
     .           dictRef="siesta:dDmax",  fmt="(f14.7)")
          endif
        endif !harrisfun/isolve
        
        if (cml_p) then
          call cmlEndPropertyList(mainXML)
        endif
        
      else !final
      ! Print out additional information in finalization.

        write(6,'(/,a)') 'siesta: Final energy (eV):'
        write(6,'(a,a15,f15.6)')
     .    'siesta: ',      'Kinetic =', Ekin/eV,
     .    'siesta: ',      'Hartree =', Uscf/eV,
     .    'siesta: ',   'Ext. field =', DUext/eV,
     .    'siesta: ',  'Exch.-corr. =', Exc/eV,
     .    'siesta: ', 'Ion-electron =', (Enascf+Enl+DUscf-Uscf-Uatm)/eV,
     .    'siesta: ',      'Ion-ion =', (Ena+Uatm-Enaatm-Eions)/eV,
     .    'siesta: ',      'Ekinion =', Ekinion/eV,
     .    'siesta: ',        'Total =', (Etot+Ekinion)/eV
        if (cml_p) then
          call cmlStartPropertyList(xf=mainXML, title='Final Energy')
          call cmlAddProperty(xf=mainXML, property=Ekin/eV,  
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ekin',  fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Uscf/eV,  
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Uscf',  fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=DUext/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:DUext', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Exc/eV,   
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Exc',   fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, 
     .         property=(Enascf+Enl+DUscf-Uscf-Uatm)/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:I-e', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, 
     .         property=(Ena+Uatm-Enaatm-Eions)/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:I-I', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Ekinion/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ekinion', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=(Etot+Ekinion)/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Etot', fmt='(f17.6)')
          call cmlEndPropertyList(mainXML)
        endif !cml_p
      endif !final

      end subroutine siesta_write_energies

      subroutine siesta_write_positions()
      use siesta_geom
      use atomlist, only: elem
      use siesta_cml
      use units
      implicit none
        if (cml_p) then
          call cmlAddMolecule(xf=mainXML, natoms=na_u, elements=elem,
     .           refs=cisa, coords=xa/Ang, style='x3', fmt='(f12.6)')
          call cmlAddLattice(xf=mainXML, cell=ucell/Ang, 
     .           units='siestaUnits:Ang', dictref='siesta:ucell')
        endif
      end subroutine siesta_write_positions

      END MODULE write_subs
