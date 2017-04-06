! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine atm_transfer()

      use atm_types, only: maxnorbs, nspecies
      use atm_types, only: species, species_info

      use radial
      use atmparams, only: NTBMAX
      use m_spin,    only: spin ! offS-SO
      use parallel,  only : IONode
!----------------------------------------------------------------
      use old_atmfuncs, only: nsmax
!
!     old_atmfuncs arrays
!
      use old_atmfuncs, only: tabpol, table, tab2
      use old_atmfuncs, only: table_offsiteSO, tab2_offsiteSO ! offS-SO
      use old_atmfuncs, only: coretab, tab2pol
      use old_atmfuncs, only: qtb, slfe
      use old_atmfuncs, only: chloctab, vlocaltab
      use old_atmfuncs, only: lmxosave, npolorbsave
      use old_atmfuncs, only: nzetasave, nsemicsave, nkblsave
!
!     old_atmfuncs procedures
!
      use old_atmfuncs, only: labelfis, izofis, zvalfis
      use old_atmfuncs, only: massfis, lomaxfis, nofis
      use old_atmfuncs, only: cnfigfio, lofio, mofio, jofio
      use old_atmfuncs, only: atmpopfio, epskb, rcut, rcut_offsiteSO
      use old_atmfuncs, only: epskb_offsiteSO 
      use old_atmfuncs, only: lmxkbfis, nkbfis


!----------------------------------------------------------------
      use ldau_specs,     only: populate_species_info_ldau

      use periodic_table, only: symbol
      use sys,            only: die

      implicit none

      type(species_info), pointer        :: spp
      type(rad_func), pointer            :: op
      type(rad_func), pointer            :: pp

      integer is, io, i , n, ntot, l, m, max_npjnl
      integer max_norbnl, nsm, izeta, j, nj_offsiteSO, ij
      integer norb, indx, ipol, num_normal, num_pol
!      integer, allocatable :: j_offsiteSO(:)  ! OffS-SO
      integer j_offsiteSO  ! OffS-SO

      real :: aj, amj 

      integer, dimension(maxnorbs) :: index_normal, z_normal,
     $     nsm_normal, index_pol, z_pol, nsm_pol
      integer, dimension(maxnorbs) :: mark  ! overkill

!     Allocate main structures for new atmfuncs -----------

      nspecies = nsmax           ! From old_atmfuncs

      allocate(species(nspecies))
!-----------------------------------------------------------

      max_npjnl = 0
      max_norbnl = 0
      do is=1,nspecies
         spp => species(is)

         spp%read_from_file = .false.

         spp%label = labelfis(is)
         spp%z     = izofis(is)
         if (spp%z.eq.-100) then
            spp%symbol = 'BS'
         else
            ! The function 'symbol' knows how to deal
            ! with (ghost) synthetics
            spp%symbol = symbol(spp%z)
         endif
         spp%zval  = zvalfis(is)
         spp%mass  = massfis(is)
         spp%self_energy  = slfe(is)


         spp%lmax_basis = lomaxfis(is)
         spp%norbs = nofis(is)

!        Check that number of orbitals is below maximum 
         if (spp%norbs.gt.maxnorbs) 
     .     call die("atm_transfer: Increase maxnorbs in atm_types.f")

         do io = 1,  spp%norbs
            spp%orb_n(io) = cnfigfio(is,io)  !! Not sure about this
            spp%orb_l(io) = lofio(is,io)
            spp%orb_m(io) = mofio(is,io)
            spp%orb_pop(io) = atmpopfio(is,io)
         enddo
!
!           Double check
!

         index_normal = 0
         index_pol = 0
         do io = 1,  spp%norbs
            norb=0 
            indx=0
            do  l=0,lmxosave(is)
               do nsm=1,nsemicsave(l,is)+1
                  do izeta=1,nzetasave(l,nsm,is)
                     norb=norb+(2*l+1)
                     indx=indx+1
                     if(norb.ge.io) then
                        index_normal(io) = indx
                        z_normal(io) = izeta
                        nsm_normal(io) = nsm
                        goto 20
                     endif
                  enddo 
               enddo 
            enddo 

            indx=0
            do  l=0,lmxosave(is)
               do nsm=1,nsemicsave(l,is)+1
                  do ipol=1, npolorbsave(l,nsm,is)
                     norb=norb+(2*(l+1)+1)
                     indx=indx+1
                     if(norb.ge.io) then
                        index_pol(io) = indx
                        z_pol(io) = ipol
                        nsm_pol(io) = nsm
                        goto 20
                     endif
                  enddo 
               enddo 
            enddo  

 20         continue

         enddo

         num_normal = maxval(index_normal)
         num_pol = maxval(index_pol)

         spp%n_orbnl = num_normal + num_pol


         allocate(spp%orbnl(spp%n_orbnl))
!
         mark = 0
         do io = 1, spp%norbs
            i = index_normal(io)
            if (i .eq. 0) cycle     ! not a normal orbital

            spp%orb_index(io) = i

            if (mark(i) .ne. 0) cycle   ! nl orb already set up
            mark(i) = 1

            l = spp%orb_l(io)
            spp%orbnl_l(i) = l
            spp%orbnl_n(i) = spp%orb_n(io)
!!            spp%orbnl_n(i) = nsm_pol(io)
            spp%orbnl_z(i) = z_normal(io)
            spp%orbnl_pop(i) = qtb(io,is) * (2*l+1)
            spp%orbnl_ispol(i) = .false.

            op => spp%orbnl(i)

            call rad_alloc(op,NTBMAX)
            op%delta  =   table(1,i,is)
            op%cutoff =   op%delta *  (NTBMAX - 1)
            op%f(1:)    = table(3:,i,is)
            op%d2(1:)   = tab2(1:,i,is)

         enddo
!
!        Polarization orbitals
!
         do io = 1, spp%norbs
            j = index_pol(io)
            if (j .eq. 0) cycle     ! not a polarization orbital
            i = num_normal + j      ! pol orbs come after normal ones
            spp%orb_index(io) = i

            if (mark(i) .ne. 0) cycle   ! nl orb already set up
            mark(i) = 1

            l = spp%orb_l(io)
            spp%orbnl_l(i) = l
            spp%orbnl_n(i) = spp%orb_n(io)
!!            spp%orbnl_n(i) = nsm_pol(io)
            spp%orbnl_z(i) = z_pol(io)
            spp%orbnl_pop(i) = qtb(io,is) * (2*l+1)
            spp%orbnl_ispol(i) = .true.

            op => spp%orbnl(i)
            call rad_alloc(op,NTBMAX)
            op%delta  =   tabpol(1,j,is)
            op%cutoff =   op%delta *  (NTBMAX - 1)
            op%f(1:)    = tabpol(3:,j,is)
            op%d2(1:)   = tab2pol(1:,j,is)
         enddo

!
!        KB projectors (relatively easy...)
!
         spp%nprojs = nkbfis(is)  ! Total number of projs per specie
         spp%lmax_projs = lmxkbfis(is)
        
CC         write(6,*) ' nprojs=',spp%nprojs
!         allocate(j_offsiteSO(spp%nprojs))

         do i = 1, spp%nprojs
          io = - i                        !! Old convention
!!!!      spp%pj_n(i)  ???????????? useful??
          spp%pj_l(i) = lofio(is,io)

           spp%pj_m(i) = mofio(is,io)
CC           write(6,'(i4,2(a,i5))') i, ' l=', spp%pj_l(i),
CC     .                                ' m=', spp%pj_m(i)  
         enddo

!
!        This piece of code assumes that the projectors are ordered
!        in their usual manner
!
         if ( .not.spin%SO_offsite ) then
!          write(spin%iout_offsiteSO,'(a)') ' passing by not offSPOrb..'

          n = 0
          ntot = 0
          do l = 0, spp%lmax_projs
             do i = 1, nkblsave(l,is)
                n = n + 1
                spp%pjnl_n(n) = i
                spp%pjnl_l(n) = l
                do m = 1, 2*l+1
                   ntot = ntot + 1
                   spp%pj_index(ntot) = n
                enddo
             enddo
          enddo
          spp%n_pjnl = n
          if (ntot .ne. spp%nprojs) call die('KB indexing...')
!          write(6,'(a)') '            '
!          write(6,'(a,i4)') ' spp%n_pjnl=', spp%n_pjnl
!          write(6,'(a,i4)') ' spp%nprojs=', spp%nprojs
!          write(6,'(a)') '            '
 
 
          allocate(spp%pjnl(spp%n_pjnl))
 
          do i = 1, spp%n_pjnl
             pp => spp%pjnl(i)
             call rad_alloc(pp,NTBMAX)
             pp%delta  =   table(1,-i,is)
             pp%f(1:)    = table(3:,-i,is)
             pp%d2(1:)   = tab2(1:,-i,is)
             if ( spin%SO_offsite .and. IONode ) then
              write(spin%iout_offsiteSO,'(a,i5,a,f12.6)') 'NTBMAX=',
     .            NTBMAX, ' delta=', pp%delta 
             endif
          enddo
!
!         Fill in the KB energy array and the cutoffs
!         A bit redundant

          do i = 1, spp%nprojs
             io = -i
             indx = spp%pj_index(i)
             spp%pjnl_ekb(indx) = epskb(is,io)
             pp => spp%pjnl(indx)
             pp%cutoff = rcut(is,io)
             if ( spin%SO_offsite .and. IONode ) then
              write(spin%iout_offsiteSO,'(2(a,i5),2(3x,a,f14.8))') 
     .          ' indx=',spp%pj_index(i),
     .          '   io=',io, ' ekb=',spp%pjnl_ekb(indx),
     .         ' rcut=',pp%cutoff
             endif
          enddo
         else
          n = 1
          ntot = 0
          nj_offsiteSO = 1
          do l = 0, spp%lmax_projs
!           if ( offSpOrb .and. l.ne.0) nj_offsiteSO = 2
!           do j_offsiteSO = 1, nj_offsiteSO
           do i = 1, nkblsave(l,is)
            if ( spin%SO_offsite .and. l.ne.0 ) nj_offsiteSO = 2
            do j_offsiteSO = 1, nj_offsiteSO
             spp%pjnl_n(n) = i  ! n of pjnl
            spp%pjnl_l(n) = l  ! l of pjnl
             if ( .not.spin%SO_offsite .or. l.eq.0 ) then 
              aj=l
             else
              aj=dble(l)+(2*j_offsiteSO-3)*0.5d0    
              spp%pj_j(n) = aj
             endif
!             do m = 1, 2*aj+1
             do m = 1, 2*l+1
              ntot = ntot + 1
              spp%jso(ntot)=j_offsiteSO
              amj = -aj + dfloat(m-1)
 
              if (.not.spin%SO_offsite) then
               spp%pj_m(ntot) = int(amj)
               spp%pj_index(ntot) = n 
 
!               write(6,'(6(a,i3))') 
!     .         ' ntot=',ntot,
!     .         ' spp%pj_index(ntot)=', spp%pj_index(ntot), 
!     .         ' l=',l, 
!     .         ' jso(ntot)=',spp%jso(ntot),
!     .         ' ml=',spp%pj_m(ntot),
!     .         ' m=',m
              else
!               spp%pj_mj(ntot) = amj
!               spp%pj_j(ntot) = jofio(is,ntot) 
               spp%pj_index(ntot) = n 
 
!               write(6,'(4(3x,a,i3))') 
!     .         ' ntot=',ntot,
!     .         ' spp%pj_index(ntot)=', spp%pj_index(ntot), 
!     .         ' l=',l, 
!     .         ' jso(ntot)=',spp%jso(ntot)
              endif

             enddo
             n = n + 1
            enddo
           enddo
          enddo
          spp%n_pjnl = n-1
! CC RC
          write(6,'(a)') '            '
          write(6,'(a,i4)') ' spp%n_pjnl=', spp%n_pjnl
          write(6,'(a,i4)') ' spp%nprojs=', spp%nprojs
          write(6,'(a)') '            '

          if (ntot .ne. spp%nprojs) call die('KB indexing...')
! CC RC  
          allocate(spp%pjnl(spp%n_pjnl))
          n = 1
          nj_offsiteSO = 1
          do l = 0, spp%lmax_projs
           do i = 1, nkblsave(l,is)
            if ( spin%SO_offsite .and. l.ne.0 ) nj_offsiteSO = 2
            do j_offsiteSO = 1, nj_offsiteSO
             if ( .not.spin%SO_offsite .or. l.eq.0 ) then 
              pp => spp%pjnl(n)
              call rad_alloc_offsiteSO(pp,NTBMAX)
              pp%delta_offsiteSO(j_offsiteSO)=
     .               table_offsiteSO(1,-spp%pjnl_n(n),l,j_offsiteSO,is)
              pp%f_offsiteSO(1:,j_offsiteSO)=
     .               table_offsiteSO(3:,-spp%pjnl_n(n),l,j_offsiteSO,is)
              pp%d2_offsiteSO(1:,j_offsiteSO)=
     .               tab2_offsiteSO(1:,-spp%pjnl_n(n),l,j_offsiteSO,is)
!              write(spin%iout_offsiteSO,'(4(a,i5),a,f12.6)') 'NTBMAX=',NTBMAX, 
!     .          ' n=', n, ' spp%pjnl_n(n)=', spp%pjnl_n(n),
!     .          ' j_offsiteSO=', j_offsiteSO, ' delta=', pp%delta_offsiteSO(j_offsiteSO) 
             else
              pp => spp%pjnl(n)
              call rad_alloc_offsiteSO(pp,NTBMAX)
              pp%delta_offsiteSO(j_offsiteSO)= 
     .               table_offsiteSO(1,-spp%pjnl_n(n),l,j_offsiteSO,is)
              pp%f_offsiteSO(1:,j_offsiteSO) = 
     .               table_offsiteSO(3:,-spp%pjnl_n(n),l,j_offsiteSO,is)
              pp%d2_offsiteSO(1:,j_offsiteSO)= 
     .               tab2_offsiteSO(1:,-spp%pjnl_n(n),l,j_offsiteSO,is)
!              write(spin%iout_offsiteSO,'(4(a,i5),a,f12.6)') 'NTBMAX=',NTBMAX, 
!     .          ' n=', n, ' spp%pjnl_n(n)=', spp%pjnl_n(n),
!     .          ' j_offsiteSO=', j_offsiteSO, ' delta=', pp%delta_offsiteSO(j_offsiteSO) 
             endif
             n = n + 1
            enddo
           enddo
          enddo
!
!        Fill in the KB energy array and the cutoffs
!        A bit redundant

! CC RC
! 
!        In the following loop, the j_offsiteSO assigment is as it will be used
!        in nlefsm, i.e., for each l it will be l +/- 1/2: 
!        l=0 --> one proj
!        l=1 --> First J_-, second J_+ and if we have semicore states 
!                for the same l, this will be repeated the number of 
!                times needed...
!
          do i = 1, spp%nprojs ! All the projs, including m_l/m_J
             io = -i
             indx = spp%pj_index(i)
!   This will be indexed as before: 1, 2 2, 3 3 3 3, etc, ...   
             spp%pjnl_ekb(indx) = epskb_offsiteSO(is,io)
             pp => spp%pjnl(indx)
             pp%cutoff_offsiteSO = rcut_offsiteSO(is,io)
!             write(spin%iout_offsiteSO,'(2(a,i5),2(3x,a,f14.8))') 
!     .         ' indx=',spp%pj_index(i),
!     .         '   io=',io, ' ekb=',spp%pjnl_ekb(indx),
!     .         ' rcut=',pp%cutoff_offsiteSO
          enddo
         endif

         call rad_alloc(spp%vna,NTBMAX)
         spp%vna%f(1:)       = table(3:,0,is)
         spp%vna%cutoff      = table(2,0,is)
         spp%vna%delta       = table(1,0,is)
         spp%vna%d2(1:)      = tab2(1:,0,is)

         call rad_alloc(spp%chlocal,NTBMAX)
         spp%chlocal%delta      = chloctab(1,1,is)
         spp%chlocal%cutoff     = chloctab(1,1,is)*(NTBMAX-1)
         spp%chlocal%f(1:)      = chloctab(2:,1,is)
         spp%chlocal%d2(1:)     = chloctab(2:,2,is)

         call rad_alloc(spp%reduced_vlocal,NTBMAX)
         spp%reduced_vlocal%delta      = vlocaltab(1,1,is)
         spp%reduced_vlocal%cutoff     = vlocaltab(1,1,is)*(NTBMAX-1)
         spp%reduced_vlocal%f(1:)      = vlocaltab(2:,1,is)
         spp%reduced_vlocal%d2(1:)     = vlocaltab(2:,2,is)

         spp%there_is_core      = (coretab(1,2,is) .eq. 1)

         call rad_alloc(spp%core,NTBMAX)
         spp%core%delta         = coretab(1,1,is)
         spp%core%cutoff        = coretab(1,1,is)*(NTBMAX-1)
         spp%core%f(1:)         = coretab(2:,1,is)
         spp%core%d2(1:)        = coretab(2:,2,is)

      enddo

      call populate_species_info_ldau

      end subroutine atm_transfer







