      module m_siesta_tddft
      implicit none
      private

      public :: siesta_tddft



      CONTAINS

      subroutine siesta_tddft(istep, rstart)
      use sys,             only: die
      use files,           only: slabel
      use siesta_options
      use wavefunctions
      use m_spin,          only: nspin
      use Kpoint_grid,     only: nkpnt, kweight, kpoint
      use atomlist,        only: no_u,no_l,no_s, indxuo
      use sparse_matrices, only: H, Hold, Dold, Dscf, maxnh, numh,
     &                           listh, listhptr, S, xijo, Escf,
     &                           H_kin, H_vkb
      use m_iodm,          only: read_dynamic_dm
      use m_state_init,    only: state_init
      use m_setup_H0
      use m_evolve,        only: evolve
      use m_eo
      use m_gamma
      use m_energies
      use units,           only: eV
      use m_compute_energies
      use m_final_H_f_stress,    only: final_H_f_stress
      use elec_dyn_options,       only: ioetot, etotfile
      implicit none
      integer, intent(in)  :: istep
      integer              :: ispn, nio
      integer              :: itded = 1 ! a TDDFT counterpart of iscf
                                      ! should be read from input file
      logical              :: rstart, found      
      double precision     :: totime
      !character*7          :: etotfile
      character*70, external  :: paste





      call state_init(istep)
      
     
     
     
      call setup_H0()
      
      
      if(rstart) then
     
      call iowavef('read',wavef,no_u,nkpnt,nspin)
      etotfile  =  paste( slabel, '.etot_vs_time'  )
      !write(6,'(2A15)') 'Debug.....!!!!!     ',etotfile
      ! call read_dynamic_dm(maxnh,no_l,nspin,numh,listhptr,listh,
      !&                      Dscf,found)
      
      rstart=.false. !To read the density matrix and TDWF only for the
      
      end if !rstart if
     


      call evolve(no_s,nspin,nspin,no_l,maxnh,maxnh,no_u,
     &            numh,listhptr,listh,numh,listhptr,listh,
     &            H,S,eo,
     &            gamma, xijo, indxuo,nkpnt,kpoint,kweight,
     &            Dscf, Escf, no_u, dt)

      totime = istep*dt
            
      Ekin = 0.0_dp
      Enl  = 0.0_dp
      do ispn = 1,min(nspin,2)
        do nio = 1,maxnh
          Ekin = Ekin + H_kin(nio) * Dscf(nio,ispn)
          Enl  = Enl  + H_vkb(nio)* Dscf(nio,ispn)
        enddo
      enddo
      !E0 = Ekin + Enl - Eions + Ena
      
      call final_H_f_stress (itded)
      DEna = Enascf - Enaatm
      Etot = E0 + DEna + DUscf + DUext + Exc + Ecorrec + Emad +
     &        Emm + Emeta


      write(6,*) 'E0     = ',E0/eV
      write(6,*)'DEna    = ',DEna/eV
      write(6,*)'DUscf   = ',DUscf/eV
      write(6,*)'DUext   = ',DUext/eV
      write(6,*)'Exc     = ',Exc/eV
      write(6,*)'Ecorrec = ',Ecorrec/eV
      write(6,*)'Emad    = ',Emad/eV
      

      call ioetot(totime, Etot)

      END subroutine siesta_tddft
      
      END module m_siesta_tddft
