      Module write_subs_energies

      use m_steps

      private
      public :: siesta_write_energies
      CONTAINS
      
      subroutine siesta_write_energies( first, iscf, dDmax )
      USE siesta_options
      use siesta_cml
      use units
      use m_energies 
      use m_spin
#ifdef TRANSIESTA
! TSS Begin
! Transiesta Modules
      use m_ts_options, only: TSinit, TSrun, TSiscf
! TSS End
#endif /* TRANSIESTA */
      implicit none

      logical :: first
      integer :: iscf
      real(dp), intent(in) :: dDmax     ! Max. change in density matrix elements
#ifdef TRANSIESTA
! TSS Begin
      character(11) :: tmpname
! TSS End
#endif /* TRANSIESTA */

      integer :: i 


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
          call cmlAddProperty(xf=mainXML, value=Eions/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Eions', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=Ena/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ena', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=Ekin/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ekin', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=Enl/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Enl', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=DEna/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:DEna', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=DUscf/eV, 
     .         units='siestaUnits:eV',
     .         dictref='siesta:DUscf', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=DUext/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:DUext', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=Exc/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Exc', fmt='r6')
          call cmlAddProperty(xf=mainXML,value=Ecorrec/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ecorrec', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=Emad/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Emad', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=Emeta/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Emeta', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=Emm/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Emm', fmt='r6')
          call cmlAddProperty(xf=mainXML,value=Ekinion/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ekinion', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=(Eharrs1+Ekinion)/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:EharrsK', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=(Etot+Ekinion)/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:EtotK', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=(FreeE+Ekinion)/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:FreeEK', fmt='r6')
          call cmlEndPropertyList(mainXML)
        endif
      endif
      ! On all SCF steps, print out the current energy (format depending on type of run)
      if (.not.final) then
        ! Print total energy and density matrix error .........................
        if (cml_p) then
          call cmlStartPropertyList(mainXML, title='SCF Cycle')
          ! Eharrs is always output
          call cmlAddProperty(xf=mainXML, value=Eharrs/eV, 
     .         units="siestaUnits:eV", 
     .         dictRef="siesta:Eharrs", fmt="r7")
        endif
        ! This chain of if statements determines which properties are output.
        if (harrisfun) then 
          write(6,"(/a,f14.6,/)") 'siesta: Eharris(eV) = ', Eharrs/eV
            ! No need for further cml output
        elseif (isolve==0) then
          if (cml_p) 
     .         call cmlAddProperty(xf=mainXML, value=FreeE/eV, 
     .         units="siestaUnits:eV", 
     .         dictRef="siesta:FreeE", fmt="r7")
          if (fixspin) then
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, value=Etot/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:Etot", fmt="r7")
              call cmlAddProperty(xf=mainXML, value=FreeE/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:FreeE", fmt="r7")
              call cmlAddProperty(xf=mainXML, value=dDmax/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:dDmax",  fmt="r7")
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
              call cmlAddProperty(xf=mainXML, value=Efs(1)/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:Efs", fmt="r7")
              call cmlAddProperty(xf=mainXML, value=Efs(2)/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:Efs", fmt="r7")
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
              call cmlAddProperty(xf=mainXML,value=Ef/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:Ef", fmt="r7")
            endif !cml_p
          endif !fixspin
        elseif (isolve==1) then
          write(6,'(/,a15,i4)') 'siesta: iscf = ',iscf
          write(6,'(a14,f15.4,a13,f15.4,a10,f7.4/)') 
     .         'Eharris(eV) = ',Eharrs/eV,
     .         '  E_KS(eV) = ',Etot/eV,'  dDmax = ',dDmax
          if (cml_p) then
            call cmlAddProperty(xf=mainXML, value=Etot/eV, 
     .           units="siestaUnits:eV", 
     .           dictRef="siesta:Etot", fmt="r7")
            call cmlAddProperty(xf=mainXML, value=dDmax/eV, 
     .           units="siestaUnits:eV", 
     .           dictRef="siesta:dDmax", fmt="r7")
          endif
        endif !harrisfun/isolve

#ifdef TRANSIESTA
! TSS Begin
        if (TSinit) then     !TSS
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
         else
          if (iscf .eq. 1) write(6,'(/,a12,3a14,2a8)')
     .             'siesta: iscf', '   Eharris(eV)',
     .             '      E_KS(eV)', '   FreeEng(eV)',
     .             '   dDmax', '  Ef(eV)'
          write(6,'(a8,i4,3f14.4,2f8.4)')
     .             'siesta: ',iscf, Eharrs/eV, Etot/eV, FreeE/eV,
     .             dDmax, Ef/eV
         end if !fixspin
        endif                !TSinit
        if (TSrun) then      !TSS
              tmpname = 'transiesta:'
              if (fixspin) then
                 if (TSiscf .eq. 1)
     .                write(6,'(/,a11,a7,3a14,a8,a7,a11)')
     .                tmpname,' TSiscf', '   Eharris(eV)',
     .                '      E_KS(eV)', '   FreeEng(eV)',
     .                '   dDmax', '  Ef_up', '  Ef_dn(eV)'
                 write(6,'(a11,i4,3f14.4,f8.4,2f9.4)')
     .                tmpname,TSiscf, Eharrs/eV, Etot/eV, FreeE/eV,
     .                dDmax, (Efs(i)/eV,i=1,2)
              else
                 if (TSiscf .eq. 1) write(6,'(/,a11,a7,3a14,2a8)')
     .                tmpname,' TSiscf', '   Eharris(eV)',
     .                '      E_KS(eV)', '   FreeEng(eV)',
     .                '   dDmax', '  Ef(eV)'
                 write(6,'(a11,i4,3f14.4,2f8.4)')
     .                tmpname,TSiscf, Eharrs/eV, Etot/eV,
     .                FreeE/eV, dDmax, Ef/eV
              endif             !fixspin
        endif                !TSrun
! TSS End
#endif /* TRANSIESTA */
        
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
          call cmlAddProperty(xf=mainXML, value=Ekin/eV,  
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ekin', fmt='r7')
          call cmlAddProperty(xf=mainXML, value=Uscf/eV,  
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Uscf', fmt='r7')
          call cmlAddProperty(xf=mainXML, value=DUext/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:DUext', fmt='r7')
          call cmlAddProperty(xf=mainXML, value=Exc/eV,   
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Exc', fmt='r7')
          call cmlAddProperty(xf=mainXML, 
     .         value=(Enascf+Enl+DUscf-Uscf-Uatm)/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:I-e', fmt='r7')
          call cmlAddProperty(xf=mainXML, 
     .         value=(Ena+Uatm-Enaatm-Eions)/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:I-I', fmt='r7')
          call cmlAddProperty(xf=mainXML, value=Ekinion/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ekinion', fmt='r7')
          call cmlAddProperty(xf=mainXML, value=(Etot+Ekinion)/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Etot', fmt='r7')
          call cmlEndPropertyList(mainXML)
        endif !cml_p
      endif !final

      end subroutine siesta_write_energies

      END MODULE write_subs_energies
