      Module write_subs_pressure

      use m_steps

      private
      public :: siesta_write_stress_pressure

      CONTAINS
      
      subroutine siesta_write_stress_pressure()

      use parallel, only: IOnode
      use precision
      USE siesta_options
      use siesta_geom
      use atomlist, only: iza
      use m_iostruct,   only: write_struct
      use siesta_cml
      use units
      use m_energies, only: FreeE
      use m_stress, only: stress, kin_stress, mstress, cstress, tstress
      use m_forces,   only: fa   ! This does not really belong here  

      implicit none

      integer :: ia, jx, ix

      real(dp):: Pmol       ! "Molecular pressure" (discounting Virial term)
      real(dp):: Psol       ! Pressure of "solid"
      real(dp):: Press      ! Pressure
      real(dp):: ps(3,3)    ! Auxiliary array

! Stress tensor and pressure:
      
      if (.not.final) then
!
!           Write Voigt components of total stress tensor 
!
            ps = stress + kin_stress
            write(6,'(/,a,6f12.2)')
     .           'Stress-tensor-Voigt (kbar):',
     .           (ps(jx,jx)/kbar,jx=1,3),
     $            ps(1,2)/kbar,
     $            ps(2,3)/kbar,
     $            ps(1,3)/kbar
            Press = - ((ps(1,1) + ps(2,2) + ps(3,3))/3.0_dp)
            write(6,"(a,f14.4)") "(Free)E + p*V (eV/cell)",
     $           (FreeE + Press*volume_of_some_cell)/eV

            if (RemoveIntraMolecularPressure) then
               ps = mstress + kin_stress
               write(6,'(/,a,6f12.2)')
     .           'Inter-Molecular-Stress-Voigt (kbar):',
     .           (ps(jx,jx)/kbar,jx=1,3),
     $            ps(1,2)/kbar,
     $            ps(2,3)/kbar,
     $            ps(1,3)/kbar
               Press = - ((ps(1,1) + ps(2,2) + ps(3,3))/3.0_dp)
               write(6,"(a,f14.4)")
     $              "(Free)E + p_inter_molec * V  (eV/cell)",
     $              (FreeE + Press*volume_of_some_cell)/eV
            endif
!
!      This use of the volume is OK, as it is called from state_analysis,
!      before possibly changing the cell.
!      Write "target enthalpy" (E + pV, where p is the *target* pressure)
            write(6,"(a,f14.4)") "Target enthalpy (eV/cell)",
     $           (FreeE + tp*volume_of_some_cell)/eV

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
              call cmlAddProperty(xf=mainXML, value=stress*Ang**3, 
     .             dictref='siesta:stress', title='Stress',
     .             units='siestaUnits:evpa3')
              call cmlAddProperty(xf=mainXML, value=Psol, 
     .             dictref='siesta:psol', title='Pressure (Static)',
     .             units='siestaUnits:kBar')
            endif !cml_p
!
            write(6,'(/,a,3(/,a,3f12.6))')
     .           'siesta: Stress tensor (total) (eV/Ang**3):',
     .           ('     ',(tstress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
            Psol = - ((tstress(1,1)+tstress(2,2) +tstress(3,3))/3.0_dp)
            write(6,'(/,a,f20.8,a)')
     .           'siesta: Pressure (total):', Psol/kBar, '  kBar'
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, value=tstress*Ang**3, 
     .             dictref='siesta:tstress', title='Total Stress',
     .             units='siestaUnits:evpa3')
              call cmlAddProperty(xf=mainXML, value=Psol,
     .             dictref='siesta:tpsol', title='Pressure (Total)',
     .             units='siestaUnits:kBar')
            endif !cml_p

            if (RemoveIntraMolecularPressure) then
             ps = mstress
             write(6,'(/,a,3(/,a,3f12.6))')
     .           'siesta: Stress tensor (nonmol) (eV/Ang**3):',
     .           ('     ',(ps(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
             Psol = - ((ps(1,1) + ps(2,2) + ps(3,3))/3.0_dp)
             write(6,'(/,a,f20.8,a)')
     .           'siesta: Pressure (nonmol):', Psol/kBar, '  kBar'
             if (cml_p) then
              call cmlAddProperty(xf=mainXML, value=ps*Ang**3, 
     .          dictref='siesta:mstress', 
     .          title='Stress tensor (normal)',
     .          units='siestaUnits:evpa3')
              call cmlAddProperty(xf=mainXML, value=Psol, 
     .          dictref='siesta:pmol', title='Pressure (Nonmol)',
     .          units='siestaUnits:kBar')
             endif                !cml_p
!
             ps = mstress + kin_stress
             write(6,'(/,a,3(/,a,3f12.6))')
     .           'siesta: Stress tensor (nonmol+kin) (eV/Ang**3):',
     .           ('     ',(ps(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
             Psol = - ((ps(1,1)+ps(2,2) +ps(3,3))/3.0_dp)
             write(6,'(/,a,f20.8,a)')
     .           'siesta: Pressure (nonmol+kin):', Psol/kBar, '  kBar'
             if (cml_p) then
              call cmlAddProperty(xf=mainXML, value=ps*Ang**3, 
     .             dictref='siesta:tmstress', 
     .             title='Stress tensor (nonmol+kin)',
     .             units='siestaUnits:evpa3')
              call cmlAddProperty(xf=mainXML, value=Psol,
     .         dictref='siesta:tpmol', title='Pressure (Nonmol+Kin)',
     .         units='siestaUnits:kBar')
             endif                !cml_p
            endif                 ! Remove intramolecular pressure

          endif                  !varcel

        ! Write Force Constant matrix if FC calculation ...
        case(6)
          call ofc(fa,dx,na_u)
        case(7)
!          call phonon_write_forces(fa,na_u,ucell,istep)
           if (IOnode) write(*,*) 'phonon suport deactivated'
        end select !idyn

      else !final

           ! AG: Possible BUG 
           ! The volume here refers to the "old" cell, and
           ! might be out of date if the cell has changed.

! Print stress tensor unconditionally
        write(6,'(/,a,3(/,a,3f12.6))')
     .       'siesta: Stress tensor (static) (eV/Ang**3):',
     .       ('siesta: ',(stress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
        if (cml_p) then
          call cmlAddProperty(xf=mainXML, value=stress*Ang**3/eV, 
     .         dictref='siesta:stress', units='siestaUnits:eV_Ang__3')
        endif !cml_p

! Print constrained stress tensor if different from unconstrained
        if (Any(cstress /= stress )) then
             write(6,'(/,a,3(/,a,3f12.6))')
     .       'siesta: Constrained stress tensor (static) (eV/Ang**3):',
     .       ('siesta: ',(cstress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
          if (cml_p) then
            call cmlAddProperty(xf=mainXML, value=cstress*Ang**3/eV, 
     .           dictref='siesta:cstress', 
     .           units='siestaUnits:eV_Ang__3')
          endif !cml_p
        endif

! Find pressure

        Psol = - (( stress(1,1) + stress(2,2) + stress(3,3) )/3.0_dp)
        Pmol = - (( mstress(1,1) + mstress(2,2) + mstress(3,3) )/3.0_dp)

        write(6,'(/,a,f18.6,a)')
     .       'siesta: Cell volume =',
     $       volume_of_some_cell/Ang**3, ' Ang**3'
        write(6,'(/,a,/,a,2a20,a,3(/,a,2f20.8,a))')
     .       'siesta: Pressure (static):',
     .       'siesta: ','Solid',        'Molecule',      '  Units',
     .       'siesta: ', Psol,           Pmol,           '  Ry/Bohr**3',
     .       'siesta: ', Psol*Ang**3/eV, Pmol*Ang**3/eV, '  eV/Ang**3',
     .       'siesta: ', Psol/kBar,      Pmol/kBar,      '  kBar'
        if (cml_p) then
          call cmlStartPropertyList(mainXML, title='Final Pressure')
          call cmlAddProperty(xf=mainXML,
     $         value=volume_of_some_cell/Ang**3, 
     .         title='cell volume', dictref='siesta:cellvol', 
     .         units='siestaUnits:Ang__3')
          call cmlAddProperty(xf=mainXML, value=Psol/kBar, 
     .         title='Pressure of Solid', dictref='siesta:pressSol', 
     .         units='siestaUnits:kbar')
          call cmlAddProperty(xf=mainXML, value=Pmol/kBar,       
     .         title='Pressure of Molecule', dictref='siesta:pressMol', 
     .         units='siestaUnits:kbar')
          call cmlEndPropertyList(mainXML)
        endif !cml_p
        
      endif !final for stress & pressure

      end subroutine siesta_write_stress_pressure
      
      
      End Module write_subs_pressure

