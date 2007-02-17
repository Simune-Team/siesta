      MODULE m_state_analysis
      use write_subs

      private
      public :: state_analysis

      CONTAINS

      subroutine state_analysis( istep  )
      USE SIESTA_TODO
      use siesta_cml
      use m_born_charge, only: born_charge
      use parallel, only: IOnode
      use m_wallclock, only : wallclock
      use zmatrix,      only: lUseZmatrix, iofaZmat, write_Zmatrix,
     .                        CartesianForce_to_ZmatForce
      use sparse_matrices
      use atomlist, only: iaorb, iphorb, iza, amass, no_u, lasto, superx
      use siesta_geom
      USE siesta_options
      use m_stress
      implicit none

      integer  :: istep

!       external :: automatic_cell,
!      .  bonds, cgvc, cgvc_zmatrix, fixed,
!      .  dhscf, diagon, dnaefs, extrapol, initatom,
!      .  iodm, iozm,
!      .  kinefsm, mulliken, naefs, neighb, !kgrid, kgridinit, 
!      .  pulayx, 
!      .  reinit, shaper, spnvec, 
!      .  timer, xijorb, memory,
!      .  ioeig, iofa, iokp, iomd, prversion !, eggbox


!------------------------------------------------------------------------- BEGIN
! Write final Kohn-Sham Energy 
      if (cml_p) call cmlStartPropertyList(mainXML,
     .                                   title='Final KS Energy')
      if (IOnode) then
        if ( .not. harrisfun) 
     .    write(6,"(/a,f14.4)")  'siesta: E_KS(eV) =        ', Etot/eV
        if (cml_p) call cmlAddProperty(xf=mainXML, property=Etot/eV,
     .       dictref='siesta:E_KS', units='siestaUnits:eV', 
     .       fmt='(f14.6)')
      endif

! Substract egg box effect form energy 
      if (eggbox_block) then
        call eggbox('energy',ucell,na_u,isa,ntm,xa,fa,Etot,
     .               eggbox_block)
        if (IOnode)
     .    write(6,"(/a,f14.4)") 'siesta: E_KS - E_eggbox = ',Etot/eV
        if (cml_p) call cmlAddProperty(xf=mainXML, property=Etot/eV,
     .         dictref='siesta:E_KS_egg', units='siestaUnits:eV', 
     .         fmt='(f14.6)')
      endif
      if (cml_p) call cmlEndPropertyList(mainXML)

#ifdef MPI
! Global reduction of forces and stresses
      allocate(fatmp(3,na_u))
      call memory('A','D',3*na_u,'siesta')
      call globalize_sum(stressl(1:3,1:3),stresstmp(1:3,1:3))
      call globalize_sum(fal(1:3,1:na_u),fatmp(1:3,1:na_u))
      stress(1:3,1:3) = stress(1:3,1:3) + stresstmp(1:3,1:3)
      fa(1:3,1:na_u) = fa(1:3,1:na_u) + fatmp(1:3,1:na_u)
      call memory('D','D',size(fatmp),'siesta')
      deallocate(fatmp)
#else
      stress(1:3,1:3) = stress(1:3,1:3) + stressl(1:3,1:3)
      fa(1:3,1:na_u) = fa(1:3,1:na_u) + fal(1:3,1:na_u)
#endif

! Substract egg box effect from the forces 
      if (eggbox_block) then
        call eggbox('forces',ucell,na_u,isa,ntm,xa,fa,Etot,eggbox_block)
      endif

! Impose constraints to atomic movements by changing forces
      call fixed(ucell,stress,na_u,isa, amass, xa, fa, cstress, cfa, 
     .           ntcon )

! Write atomic forces 
      fmax = 0.0_dp
      cfmax = 0.0_dp
      fres = 0.0_dp
      do ix = 1,3
        ftot(ix) = 0.0_dp
        do ia = 1,na_u
          ftem = fa(ix,ia)
          cftem = cfa(ix,ia)
          ftot(ix) = ftot(ix) + ftem
          fres = fres + ftem*ftem
          fmax = max( fmax, dabs(ftem) )
          cfmax = max( cfmax, dabs(cftem) )
        enddo
      enddo
      fres = dsqrt( fres / (3.0_dp*na_u) )

! Calculate and output Zmatrix forces
      if (lUseZmatrix) then
        call CartesianForce_to_ZmatForce(na_u,xa,fa)
        if (IOnode) call iofaZmat()
      endif

! Add kinetic term to stress tensor
      do ia = 1,na_u
        do jx = 1,3
          do ix = 1,3
            tstress(ix,jx) = stress(ix,jx) -
     .             amu * amass(ia) * va(ix,ia) * va(jx,ia) / volume
          enddo
        enddo
      enddo

! Force output 
      if (IOnode) then
        call siesta_write_forces()
        call siesta_write_stress_pressure()
        call wallclock('--- end of geometry step')
      endif

! Mulliken population analysis 
      call mulliken( mullipop, nspin, na_u, no_u, maxnh,
     .               numh, listhptr, listh, S, Dscf, isa, 
     .               lasto, iaorb, iphorb )

!     Call the born effective charge routine only in those steps (even) 
!     in which the dx  is positive.
      if (bornz .and. (mod(istep,2) .eq. 0)) then
         call born_charge()
      endif
!--------------------------------------------------------------------------- END
      END subroutine state_analysis

      END MODULE m_state_analysis
