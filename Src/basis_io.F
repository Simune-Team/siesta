      module basis_io
!
!     Support for dumping and reading PAO and KB information from
!     ASCII or NetCDF files.
!
!     Alberto Garcia, 2000, 2001
!
      use chemical
      use sys
      use precision
      use atmparams, only:NTBMAX
      use atm_types
      use basis_types
      use pseudopotential
      use radial
      use fdf

      implicit none

      private
      public :: dump_basis_ascii, read_basis_ascii
      public :: dump_basis_xml
      public :: dump_basis_netcdf, read_basis_netcdf

      type(species_info), pointer        :: spp
      type(rad_func), pointer            :: op
      type(rad_func), pointer            :: pp
      type(rad_func), pointer            :: func

      CONTAINS

      subroutine read_basis_netcdf

#ifndef CDF
      call die(
     $ '*** You need netCDF to read the new user-defined basis files...'
     $ )
      end subroutine  read_basis_netcdf
#else 

      use netcdf

      implicit none


      integer ncid, iret

      integer nkbs, nkbs_id, ntb_id, proj_id,
     $     pjnl_l_id, pjnl_n_id, pjnl_ekb_id, kbdelta_id,
     $     kbcutoff_id
      integer norbs, norbs_id, orbnl_l_id, orbnl_n_id, orbnl_z_id,
     $     cutoff_id, delta_id, orb_id, orbnl_pop_id, orbnl_ispol_id
      integer vna_id, chlocal_id, core_id

      integer aux(maxnorbs)

      integer is, j, i, l, nrp_tables, core_flag, nor, nk, m

      character(len=20) filename
      character(len=20) dummy
                                !
      call read_chemical_types
      nspecies = number_of_species()

      allocate(species(nspecies))

      do is = 1, nspecies
         spp => species(is)
         spp%label = species_label(is)
         spp%read_from_file = .true.
         write(filename,'(a,a)') trim(spp%label), ".ion.nc"
         iret = nf90_open(trim(filename),NF90_NOWRITE,ncid)

         iret = nf90_inq_dimid(ncid,'norbs',norbs_id)
         iret = nf90_inquire_dimension(ncid,norbs_id,len=norbs)
         if (norbs .gt. maxnorbs)
     $  call die("read_user_basis: Increase maxnorbs in atm_types.f")

         spp%n_orbnl = norbs

         iret = nf90_inq_dimid(ncid,'nkbs',nkbs_id)
         iret = nf90_inquire_dimension(ncid,nkbs_id,len=nkbs)
         spp%n_pjnl = nkbs

!
!        For now, it is assumed that *all* the radial arrays have
!        the same length.
!
         iret = nf90_inq_dimid(ncid,'ntb',ntb_id)
         iret = nf90_inquire_dimension(ncid,ntb_id,len=nrp_tables)
         if (nrp_tables .ne. NTBMAX) call die("NTBMAX mismatch")

         allocate(spp%orbnl(norbs))
         allocate(spp%pjnl(nkbs))


         iret = nf90_get_att(ncid,nf90_global,'Element',spp%symbol)
!!!        iret = nf90_get_att(ncid,nf90_global,'Label',dummy)
                                !! Sanity check here??

         iret = nf90_get_att(ncid,nf90_global,'Atomic_number',spp%z)
         if (atomic_number(is) .ne. spp%z)
     $                            call die("Atomic number mismatch")

         iret = nf90_get_att(ncid,nf90_global,'Valence_charge',spp%zval)
         iret = nf90_get_att(ncid,nf90_global,'Mass',spp%mass)
         iret = nf90_get_att(ncid,nf90_global,'Self_energy',
     $                                         spp%self_energy)
         iret = nf90_get_att(ncid,nf90_global,
     $        'Number_of_orbitals',spp%norbs)
         iret = nf90_get_att(ncid,nf90_global,
     $        'L_max_basis',spp%lmax_basis)
         iret = nf90_get_att(ncid,nf90_global,
     $        'Number_of_projectors',spp%nprojs)
         iret = nf90_get_att(ncid,nf90_global,
     $        'L_max_projs',spp%lmax_projs)

                                !! Orbitals

         iret = nf90_inq_varid(ncid,'orbnl_l',orbnl_l_id)
         iret = nf90_inq_varid(ncid,'orbnl_n',orbnl_n_id)
         iret = nf90_inq_varid(ncid,'orbnl_z',orbnl_z_id)
         iret = nf90_inq_varid(ncid,'orbnl_ispol',orbnl_ispol_id)
         iret = nf90_inq_varid(ncid,'orbnl_pop',orbnl_pop_id)

         iret = nf90_inq_varid(ncid,'cutoff',cutoff_id)
         iret = nf90_inq_varid(ncid,'delta',delta_id)

                                !!      Projectors

         iret = nf90_inq_varid(ncid,'pjnl_l',pjnl_l_id)
         call check(iret)
         iret = nf90_inq_varid(ncid,'pjnl_n',pjnl_n_id)
         iret = nf90_inq_varid(ncid,'pjnl_ekb',pjnl_ekb_id)
         iret = nf90_inq_varid(ncid,'kbcutoff',kbcutoff_id)
         iret = nf90_inq_varid(ncid,'kbdelta',kbdelta_id)
         call check(iret)

         iret = nf90_inq_varid(ncid,'orb',orb_id)
         call check(iret)
                                !
                                !       Local potential
                                !
         iret = nf90_inq_varid(ncid,'vna',vna_id)
         iret = nf90_get_att(ncid,vna_id,
     $        'Vna_cutoff',spp%vna%cutoff)
         iret = nf90_get_att(ncid,vna_id,
     $        'Vna_delta',spp%vna%delta)
                                !
                                !       Local potential charge density
                                !
         iret = nf90_inq_varid(ncid,'chlocal',chlocal_id)
         iret = nf90_get_att(ncid,chlocal_id,
     $        'Chlocal_cutoff',spp%chlocal%cutoff)
         iret = nf90_get_att(ncid,chlocal_id,
     $        'Chlocal_delta',spp%chlocal%delta)
                                !
                                !       Core charge
                                !
         iret = nf90_get_att(ncid,nf90_global,'Core_flag',core_flag)
         spp%there_is_core = (core_flag .eq. 1)

         if (spp%there_is_core) then
            iret = nf90_inq_varid(ncid,'core',core_id)
            iret = nf90_get_att(ncid,core_id,
     $           'Core_cutoff',spp%core%cutoff)
            iret = nf90_get_att(ncid,core_id,
     $           'Core_delta',spp%core%delta)
         else
            call rad_zero(spp%core)
         endif

         call check(iret)
         iret = nf90_inq_varid(ncid,'proj',proj_id)
         call check(iret)

         iret = nf90_get_var(ncid,pjnl_l_id,spp%pjnl_l,count=(/nkbs/))
         call check(iret)
         iret = nf90_get_var(ncid,pjnl_n_id,spp%pjnl_n,count=(/nkbs/))
         call check(iret)
         iret = nf90_get_var(ncid,pjnl_ekb_id,spp%pjnl_ekb,
     $        count=(/nkbs/))
         call check(iret)

         iret=nf90_get_var(ncid,orbnl_l_id,spp%orbnl_l,count=(/norbs/))
         iret=nf90_get_var(ncid,orbnl_n_id,spp%orbnl_n,count=(/norbs/))
         iret=nf90_get_var(ncid,orbnl_z_id,spp%orbnl_z,count=(/norbs/))

         iret = nf90_get_var(ncid,orbnl_ispol_id,aux,count=(/norbs/))
         do i = 1, norbs
            spp%orbnl_ispol(i) =  aux(i) .eq. 1
         enddo
         call check(iret)
         iret = nf90_get_var(ncid,orbnl_pop_id,
     $        spp%orbnl_pop,count=(/norbs/))
         call check(iret)

c     
         nk = 0
         do i = 1, nkbs
            pp => spp%pjnl(i)
            call rad_alloc(pp,NTBMAX)
            iret = nf90_get_var(ncid,proj_id,pp%f(1:),
     $           start=(/1,i/),count=(/NTBMAX,1/))
            call check(iret)
            iret = nf90_get_var(ncid,kbcutoff_id,pp%cutoff,
     $           start=(/i/))
            call check(iret)
            iret = nf90_get_var(ncid,kbdelta_id,pp%delta,
     $           start=(/i/))
            call check(iret)
            call rad_setup_d2(pp)
            l = spp%pjnl_l(i)
            do m = -l,l
               nk = nk+1
               spp%pj_n(nk) = spp%pjnl_n(i)
               spp%pj_l(nk) = spp%pjnl_l(i)
               spp%pj_m(nk) = m
               spp%pj_index(nk) = i
            enddo
         enddo
         spp%nprojs = nk
                                !
                                !       Local potential
                                !
         call rad_alloc(spp%vna,NTBMAX)
         iret = nf90_get_var(ncid,vna_id,spp%vna%f(1:),
     $        start=(/1/),count=(/NTBMAX/))
         call check(iret)
         call rad_setup_d2(spp%vna)

                                !
                                !       Local potential charge density
                                !
         call rad_alloc(spp%chlocal,NTBMAX)
         iret = nf90_get_var(ncid,chlocal_id,spp%chlocal%f(1:),
     $        start=(/1/),count=(/NTBMAX/))
         call check(iret)
         call rad_setup_d2(spp%chlocal)

         if (spp%there_is_core) then
            call rad_alloc(spp%core,NTBMAX)
            iret = nf90_get_var(ncid,core_id,spp%core%f(1:),
     $           start=(/1/),count=(/NTBMAX/))
            call check(iret)
            call rad_setup_d2(spp%core)
         endif

         nor = 0
         do i = 1, norbs
            op => spp%orbnl(i)
            call rad_alloc(op,NTBMAX)
            iret = nf90_get_var(ncid,orb_id,op%f(1:),
     $           start=(/1,i/),count=(/NTBMAX,1/))
            call check(iret)
            iret = nf90_get_var(ncid,cutoff_id,op%cutoff,
     $           start=(/i/))
            call check(iret)
            iret = nf90_get_var(ncid,delta_id,op%delta,
     $           start=(/i/))
            call check(iret)
            call rad_setup_d2(op)
            l = spp%orbnl_l(i)
            do m = -l,l
               nor = nor+1
               spp%orb_n(nor) = spp%orbnl_n(i)
               spp%orb_l(nor) = spp%orbnl_l(i)
               spp%orb_m(nor) = m
               spp%orb_pop(nor) = spp%orbnl_pop(i) / (2*l+1)
               spp%orb_index(nor) = i
            enddo
         enddo
         spp%norbs = nor
         iret = nf90_close(ncid)
         call check(iret)

      enddo
!
      call elec_corr_setup()

      CONTAINS

      subroutine check(status)
      
      integer, intent(in):: status
      if (status .ne. nf90_noerr) then
         print  *, trim(nf90_strerror(status))
         stop 'Stopped'
      endif
      end subroutine check

      end subroutine read_basis_netcdf
#endif

!----------------------------------------------------------------
      subroutine elec_corr_setup

      use atmfuncs, only: floating

      integer is, is2, i
      real(dp) rchloc, rchloc2

      type(rad_func), pointer :: func

      npairs = ((nspecies+1)*nspecies)/2
      allocate(elec_corr(npairs))
      do is=1,npairs
        func=>elec_corr(is)
        call rad_zero(func)
      enddo

      do is=1,nspecies
         rchloc = species(is)%chlocal%cutoff
         do is2=is,1,-1
            rchloc2 = species(is2)%chlocal%cutoff
            i = ((is-1)*is)/2+is2
            func => elec_corr(i)
            if (floating(is) .or. floating(is2)) then
               call rad_zero(func)
            else
               call rad_alloc(func,NTBMAX)
               func%cutoff = rchloc + rchloc2 + 0.2_dp
               func%delta =  func%cutoff / (NTBMAX - 1)
               call ch_overlap(is,is2,func%cutoff,func%f)
               call rad_setup_d2(func)
            endif
         enddo
      enddo

      end subroutine elec_corr_setup
!----------------------------------------------------------------
#ifdef CDF
#endif
!=======================================================================

      subroutine read_basis_ascii

      implicit none

      integer lun
      character*20 filename


      integer is, j, i, l, ispol, nor, nk, m, npts
      real(dp) dummy

      call read_chemical_types
      nspecies = number_of_species()

      allocate(species(nspecies))

      do is = 1, nspecies
         spp => species(is)
         spp%label = species_label(is)
         spp%read_from_file = .true.
         call read_ion_ascii(spp)
      enddo
      call elec_corr_setup()

      end subroutine read_basis_ascii
!
!----------------------
      subroutine read_ion_ascii(spp)
      type(species_info), pointer  :: spp
      type(rad_func), pointer      :: op

      character(len=20) filename
      integer i, j, l, m, lun, nor, nk, ispol, npts
      real(dp) dummy

      write(filename,'(a,a)') trim(spp%label), ".ion"
      call io_assign(lun)
      open(lun,file=filename,status='old',form='formatted')
      rewind(lun)

      call read_header(spp,lun)
      read(lun,*) 
         do i=1,spp%n_orbnl
            read(lun,*)
     $           spp%orbnl_l(i), spp%orbnl_n(i), spp%orbnl_z(i),
     $           ispol, spp%orbnl_pop(i)
            spp%orbnl_ispol(i) =  ispol.eq.1
            call radial_read_ascii(spp%orbnl(i),lun)
         enddo
!
!        Update indexes
!
         nor = 0
         do i = 1, spp%n_orbnl
            l = spp%orbnl_l(i)
            do m = -l,l
               nor = nor+1
               spp%orb_n(nor) = spp%orbnl_n(i)
               spp%orb_l(nor) = spp%orbnl_l(i)
               spp%orb_m(nor) = m
               spp%orb_pop(nor) = spp%orbnl_pop(i) / (2*l+1)
               spp%orb_index(nor) = i
            enddo
         enddo
         spp%norbs = nor

! KBs
         read(lun,*)
         do i=1,spp%n_pjnl
            read(lun,*)
     $           spp%pjnl_l(i), spp%pjnl_n(i), spp%pjnl_ekb(i)
            call radial_read_ascii(spp%pjnl(i),lun)
         enddo
!
!        Update indexes
!

         nk = 0
         do i = 1, spp%n_pjnl
            l = spp%pjnl_l(i)
            do m = -l,l
               nk = nk+1
               spp%pj_n(nk) = spp%pjnl_n(i)
               spp%pj_l(nk) = spp%pjnl_l(i)
               spp%pj_m(nk) = m
               spp%pj_index(nk) = i
            enddo
         enddo
         spp%nprojs = nk
!
!Vna
         read(lun,*)
         call radial_read_ascii(spp%vna,lun)

!
!Chlocal
         read(lun,*)
         call radial_read_ascii(spp%chlocal,lun)
!
!Core
         read(lun,*,end=9999)
         call radial_read_ascii(spp%core,lun)

 9999    continue
         call io_close(lun)

      CONTAINS

      subroutine read_header(p,unit)

      type(species_info), pointer :: p
      integer, intent(in)         :: unit
      
      character(len=78) line

      read(unit,'(a)') line
      if (trim(line) .eq. '<preamble>') then
 1       continue
         read(unit,'(a)') line
         if (trim(line) .ne. '</preamble>') goto 1
      endif
         
      read(unit,'(a2)') p%symbol
      read(unit,'(a20)') p%label
      read(unit,*) p%z
      read(unit,*) p%zval
      read(unit,*) p%mass
      read(unit,*) p%self_energy
      read(unit,*) p%lmax_basis, p%n_orbnl
      read(unit,*) p%lmax_projs, p%n_pjnl

      allocate(p%orbnl(p%n_orbnl))
      allocate(p%pjnl(p%n_pjnl))

      end subroutine read_header

      end subroutine read_ion_ascii


#ifndef CDF
      subroutine dump_basis_netcdf
!     Do nothing
      end subroutine dump_basis_netcdf
#else
      subroutine dump_basis_netcdf

      use netcdf

      implicit none

      integer ncid, iret

      integer nkbs, nkbs_id, ntb_id, proj_id,
     $        pjnl_l_id, pjnl_n_id, pjnl_ekb_id, kbdelta_id,
     $        kbcutoff_id
      integer norbs, norbs_id, orbnl_l_id, orbnl_n_id, orbnl_z_id,
     $        cutoff_id, delta_id, orb_id, orbnl_pop_id, orbnl_ispol_id
      integer vna_id, chlocal_id, core_id

      integer aux(maxnorbs)

      integer is, j, i, l
      character*20 filename

      do is = 1, nspecies
         spp => species(is)
         if (spp%read_from_file) cycle   !! Do not dump

         write(filename,'(a,a)') trim(spp%label), ".ion.nc"
         write(6,'(2a)') 'Dumping basis to NetCDF file ',
     $                   trim(filename)

         iret = nf90_create(trim(filename),NF90_CLOBBER,ncid)

        nkbs =  spp%n_pjnl
        norbs = spp%n_orbnl

        iret = nf90_def_dim(ncid,'norbs',norbs,norbs_id)
        call check(iret)
        iret = nf90_def_dim(ncid,'nkbs',nkbs,nkbs_id)
        call check(iret)
        iret = nf90_def_dim(ncid,'ntb',NTBMAX,ntb_id)
        call check(iret)


!! Orbitals

        iret = nf90_put_att(ncid,nf90_global,'Element',spp%symbol)
        iret = nf90_put_att(ncid,nf90_global,'Label',spp%label)
        iret = nf90_put_att(ncid,nf90_global,'Atomic_number',spp%z)
        iret = nf90_put_att(ncid,nf90_global,'Valence_charge',spp%zval)
        iret = nf90_put_att(ncid,nf90_global,'Mass',spp%mass)
        iret = nf90_put_att(ncid,nf90_global,'Self_energy',
     $                                       spp%self_energy)
        iret = nf90_put_att(ncid,nf90_global,
     $                      'Number_of_orbitals',spp%norbs)
        iret = nf90_put_att(ncid,nf90_global,
     $                      'L_max_basis',spp%lmax_basis)
        iret = nf90_put_att(ncid,nf90_global,
     $                      'Number_of_projectors',spp%nprojs)
        iret = nf90_put_att(ncid,nf90_global,
     $                      'L_max_projs',spp%lmax_projs)

        iret = nf90_def_var(ncid,'orbnl_l',nf90_int,norbs_id,orbnl_l_id)
        iret = nf90_def_var(ncid,'orbnl_n',nf90_int,norbs_id,orbnl_n_id)
        iret = nf90_def_var(ncid,'orbnl_z',nf90_int,norbs_id,orbnl_z_id)
        iret = nf90_def_var(ncid,'orbnl_ispol',nf90_int,
     $                            norbs_id,orbnl_ispol_id)
        iret = nf90_def_var(ncid,'orbnl_pop',nf90_double,
     $                            norbs_id,orbnl_pop_id)

        iret = nf90_def_var(ncid,'cutoff',nf90_double,
     $                           norbs_id,cutoff_id)
        iret = nf90_def_var(ncid,'delta',nf90_double,
     $                           norbs_id,delta_id)

!!      Projectors

        iret = nf90_def_var(ncid,'pjnl_l',nf90_int,nkbs_id,pjnl_l_id)
        call check(iret)
        iret = nf90_def_var(ncid,'pjnl_n',nf90_int,nkbs_id,pjnl_n_id)
        iret = nf90_def_var(ncid,'pjnl_ekb',nf90_double,
     $                           nkbs_id,pjnl_ekb_id)
        iret = nf90_def_var(ncid,'kbcutoff',nf90_double,
     $                           nkbs_id,kbcutoff_id)
        iret = nf90_def_var(ncid,'kbdelta',nf90_double,
     $                           nkbs_id,kbdelta_id)
        call check(iret)

        iret = nf90_def_var(ncid,'orb',nf90_double,
     $                      (/ntb_id,norbs_id/),orb_id)
        call check(iret)
!
!       Local potential
!
        iret = nf90_def_var(ncid,'vna',nf90_double,
     $                      (/ntb_id/),vna_id)
        iret = nf90_put_att(ncid,vna_id,
     $                      'Vna_cutoff',spp%vna%cutoff)
        iret = nf90_put_att(ncid,vna_id,
     $                      'Vna_delta',spp%vna%delta)
!
!       Local potential charge density
!
        iret = nf90_def_var(ncid,'chlocal',nf90_double,
     $                      (/ntb_id/),chlocal_id)
        iret = nf90_put_att(ncid,chlocal_id,
     $                      'Chlocal_cutoff',spp%chlocal%cutoff)
        iret = nf90_put_att(ncid,chlocal_id,
     $                      'Chlocal_delta',spp%chlocal%delta)
!
!       Core charge
!
        if (spp%there_is_core) then
           iret = nf90_put_att(ncid,nf90_global,
     $                         'Core_flag',1)
           iret = nf90_def_var(ncid,'core',nf90_double,
     $                      (/ntb_id/),core_id)
           iret = nf90_put_att(ncid,core_id,
     $                      'Core_cutoff',spp%core%cutoff)
           iret = nf90_put_att(ncid,core_id,
     $                      'Core_delta',spp%core%delta)
        else
           iret = nf90_put_att(ncid,nf90_global,
     $                         'Core_flag',0)
        endif

        call check(iret)
        iret = nf90_def_var(ncid,'proj',nf90_double,
     $                      (/ntb_id,nkbs_id/),proj_id)
        call check(iret)

!!!!!!!!
        iret = nf90_enddef(ncid)
        call check(iret)

        iret = nf90_put_var(ncid,pjnl_l_id,spp%pjnl_l,count=(/nkbs/))
        call check(iret)
        iret = nf90_put_var(ncid,pjnl_n_id,spp%pjnl_n,count=(/nkbs/))
        call check(iret)
        iret = nf90_put_var(ncid,pjnl_ekb_id,spp%pjnl_ekb,
     $                           count=(/nkbs/))
        call check(iret)

        iret = nf90_put_var(ncid,orbnl_l_id,spp%orbnl_l,count=(/norbs/))
        iret = nf90_put_var(ncid,orbnl_n_id,spp%orbnl_n,count=(/norbs/))
        iret = nf90_put_var(ncid,orbnl_z_id,spp%orbnl_z,count=(/norbs/))

        if (norbs .gt. maxnorbs)
     $        call die("dump_basis_netcdf: Increase maxnorbs")
        aux(1:norbs) = 0
        do i = 1, norbs
           if (spp%orbnl_ispol(i)) aux(i)=1
        enddo
        iret = nf90_put_var(ncid,orbnl_ispol_id,aux,count=(/norbs/))
        call check(iret)
        iret = nf90_put_var(ncid,orbnl_pop_id,
     $                           spp%orbnl_pop,count=(/norbs/))
        call check(iret)

c
        do i = 1, nkbs
           pp => spp%pjnl(i)
           iret = nf90_put_var(ncid,proj_id,pp%f(1:),
     $                      start=(/1,i/),count=(/NTBMAX,1/))
           call check(iret)
           iret = nf90_put_var(ncid,kbcutoff_id,pp%cutoff,
     $                      start=(/i/))
           call check(iret)
           iret = nf90_put_var(ncid,kbdelta_id,pp%delta,
     $                      start=(/i/))
           call check(iret)
        enddo
!
!       Local potential
!
        iret = nf90_put_var(ncid,vna_id,spp%vna%f(1:),
     $                      start=(/1/),count=(/NTBMAX/))
        call check(iret)

!
!       Local potential charge density
!
        iret = nf90_put_var(ncid,chlocal_id,spp%chlocal%f(1:),
     $                      start=(/1/),count=(/NTBMAX/))
        call check(iret)

        if (spp%there_is_core) then
           iret = nf90_put_var(ncid,core_id,spp%core%f(1:),
     $                      start=(/1/),count=(/NTBMAX/))
           call check(iret)
        endif

        do i = 1, norbs
           op => spp%orbnl(i)
           iret = nf90_put_var(ncid,orb_id,op%f(1:),
     $                      start=(/1,i/),count=(/NTBMAX,1/))
           call check(iret)
           iret = nf90_put_var(ncid,cutoff_id,op%cutoff,
     $                      start=(/i/))
           call check(iret)
           iret = nf90_put_var(ncid,delta_id,op%delta,
     $                      start=(/i/))
           call check(iret)
        enddo
        iret = nf90_close(ncid)
        call check(iret)

      enddo

      contains
      subroutine check(status)
      
      integer, intent(in):: status
      if (status .ne. nf90_noerr) then
         print  *, trim(nf90_strerror(status))
         stop 'Stopped'
      endif
      end subroutine check

      end subroutine dump_basis_netcdf

#endif


      subroutine dump_basis_ascii

      integer is

      do is = 1, nspecies
         call dump_ion_ascii(is)
      enddo
      end subroutine dump_basis_ascii
!---
      subroutine dump_ion_ascii(is)
      integer, intent(in)      :: is

      type(species_info), pointer  ::   spp
      character*20 filename
      integer j, i, l, ispol, lun
      type(rad_func), pointer :: op

      spp => species(is)
      if (spp%read_from_file) return    !! Do not dump

      write(filename,'(a,a)') trim(spp%label), ".ion"
      call io_assign(lun)
      open(lun,file=filename,status='replace',form='formatted')

      write(lun,'(a)') '<preamble>'
      call write_basis_specs(lun,is)
      call pseudo_header_print(lun,
     $                   basis_parameters(is)%pseudopotential)
      write(lun,'(a)') '</preamble>'
      call write_header(spp,lun)
      write(lun,'(a)') "# PAOs:__________________________"
         do i=1,spp%n_orbnl
            ispol = 0
            if (spp%orbnl_ispol(i)) ispol = 1
            write(lun,'(4i3,f10.6,2x,a)')
     $           spp%orbnl_l(i), spp%orbnl_n(i), spp%orbnl_z(i),
     $           ispol, spp%orbnl_pop(i),
     $           " #orbital l, n, z, is_polarized, population"
            call radial_dump_ascii(spp%orbnl(i),lun)
         enddo
         
         write(lun,'(a)') "# KBs:__________________________"
         do i=1,spp%n_pjnl
            write(lun,'(2i3,f22.16,2x,a)')
     $           spp%pjnl_l(i), spp%pjnl_n(i), spp%pjnl_ekb(i),
     $           " #kb l, n (sequence number), Reference energy"
            call radial_dump_ascii(spp%pjnl(i),lun)
         enddo

         write(lun,'(a)') "# Vna:__________________________"
         call radial_dump_ascii(spp%vna,lun)

         write(lun,'(a)') "# Chlocal:__________________________"
         call radial_dump_ascii(spp%chlocal,lun)

         if (spp%there_is_core) then
            write(lun,'(a)') "# Core:__________________________"
            call radial_dump_ascii(spp%core,lun)
         endif
         call io_close(lun)

         CONTAINS

         subroutine write_header(p,unit)
         
         type(species_info), pointer :: p
         integer, intent(in)         :: unit

         write(unit,'(a2,28x,a)') p%symbol, "# Symbol"
         write(unit,'(a20,10x,a)') p%label, "# Label"
         write(unit,'(i5,25x,a)') p%z, "# Atomic number"
         write(unit,'(i5,25x,a)') p%zval, "# Valence charge"
         write(unit,'(g22.12,4x,a)') p%mass, "# Mass "
         write(unit,'(g22.12,4x,a)') p%self_energy, "# Self energy "
         write(unit,'(2i4,22x,a)') p%lmax_basis, p%n_orbnl,
     $        "# Lmax for basis, no. of nl orbitals "
         write(unit,'(2i4,22x,a)') p%lmax_projs, p%n_pjnl,
     $        "# Lmax for projectors, no. of nl KB projectors "

         end subroutine write_header

      end subroutine dump_ion_ascii
!
      subroutine dump_basis_xml

      integer is

      do is = 1, nspecies
         call dump_ion_xml(is)
      enddo
      end subroutine dump_basis_xml
!---
      subroutine dump_ion_xml(is)
      use xml

      integer, intent(in)      :: is

      type(species_info), pointer  ::   spp
      character*20 filename
      integer j, i, l, ispol, lun
      type(rad_func), pointer :: op

      spp => species(is)
      if (spp%read_from_file) return    !! Do not dump

      write(filename,'(a,a)') trim(spp%label), ".ion.xml"
      call io_assign(lun)
      open(lun,file=filename,status='replace',form='formatted')

      write(lun,'(a)') "<ion> "
      call xml_dump_element(lun,'symbol',str(spp%symbol))
      call xml_dump_element(lun,'label',str(spp%label))
      call xml_dump_element(lun,'z',str(spp%z))
      call xml_dump_element(lun,'valence',str(spp%zval))
      call xml_dump_element(lun,'mass',str(spp%mass))
      call xml_dump_element(lun,'self_energy',str(spp%self_energy))
      call xml_dump_element(lun,'lmax_basis',str(spp%lmax_basis))
      call xml_dump_element(lun,'norbs_nl',str(spp%n_orbnl))
      call xml_dump_element(lun,'lmax_projs',str(spp%lmax_projs))
      call xml_dump_element(lun,'nprojs_nl',str(spp%n_pjnl))

      write(lun,'(a)') '<preamble>'
      call write_basis_specs(lun,is)
      call pseudo_header_print(lun,
     $                   basis_parameters(is)%pseudopotential)
      write(lun,'(a)') '</preamble>'

      write(lun,'(a)') "<paos>"
         do i=1,spp%n_orbnl
            write(lun,'(a)') "<orbital "
            ispol = 0
            if (spp%orbnl_ispol(i)) ispol = 1
            call xml_dump_attribute(lun,'l',str(spp%orbnl_l(i)))
            call xml_dump_attribute(lun,'n',str(spp%orbnl_n(i)))
            call xml_dump_attribute(lun,'z',str(spp%orbnl_z(i)))
            call xml_dump_attribute(lun,'ispol',str(ispol))
            call xml_dump_attribute(lun,'population',
     $                            str(spp%orbnl_pop(i)))
            write(lun,'(a)') " >"
            call radial_dump_xml(spp%orbnl(i),lun)
            write(lun,'(a)') "</orbital>"
         enddo
      write(lun,'(a)') "</paos>"
         
      write(lun,'(a)') "<kbs>"
      do i=1,spp%n_pjnl
         write(lun,'(a)') "<projector "
         call xml_dump_attribute(lun,'l',str(spp%pjnl_l(i)))
         call xml_dump_attribute(lun,'n',str(spp%pjnl_n(i)))
         call xml_dump_attribute(lun,'ref_energy',
     $                                   str(spp%pjnl_ekb(i)))
         call radial_dump_xml(spp%pjnl(i),lun)
         write(lun,'(a)') " >"
       enddo
       write(lun,'(a)') "</kbs>"

       write(lun,'(a)') "<vna>"
       call radial_dump_xml(spp%vna,lun)
       write(lun,'(a)') "</vna>"

       write(lun,'(a)') "<chlocal>"
       call radial_dump_xml(spp%chlocal,lun)
       write(lun,'(a)') "</chlocal>"

       if (spp%there_is_core) then
          write(lun,'(a)') "<core>"
          call radial_dump_xml(spp%core,lun)
          write(lun,'(a)') "</core>"
       endif

       write(lun,'(a)') "</ion>"

       call io_close(lun)

      end subroutine dump_ion_xml
!
!
!======================================================================
!
!     This "atomic" routines use information **in the new structures**,
!     so they cannot be replaced by those in atom.f...
!
      SUBROUTINE CH_OVERLAP(IS1,IS2,RMX,CORR)
      use atmfuncs, only: izvalfis, psch
      use ionew, only: IOnode

      integer, intent(in)   :: is1, is2
      real(dp), intent(in )    :: rmx
      real(dp), intent(out)    :: corr(:)

C     Returns a table with the difference between the electrostatic energy 
C     of two spherical charge-densities and two punctual charges with the 
C     same total charge as a function of the distance between the centers 
C     of these charge densities. 
C     Written by D.Sanchez-Portal. March, 1997.(from routine MATEL, written 
C     by Jose M. Soler)

C     INTEGER IS1,IS2             :  Species indexes.
C     RMX                         :  Maximum range of the correction.
C     CORR(NTBMAX)                :  Electrostatic correction energy.

C     Distances in Bohr. Energy in Rydbergs.

C     Internal precision parameters  ------------------------------------
C     NQ is the number of radial points in reciprocal space.
C     Npoint , 2npoint+1 is the number of points used by RATINT in the 
C     interpolation.
C     Q2CUT is the required planewave cutoff for the expansion of
C     the 'local-pseudopotential atomic charge density'
C     (in Ry if lengths are in Bohr).
C     CHERR is a small number to check the precision of the charge density
C     integration.

      integer nq, npoint
      real(dp)            :: q2cut, cherr
      PARAMETER ( NQ     =  512  )
      PARAMETER ( NPOINT =  4     ) 
      PARAMETER ( Q2CUT  =  2.5e3_dp )
      PARAMETER ( CHERR   =  5.e-2_dp )

      real(dp)
     .     CH(0:NQ,2),VTB(NTBMAX,2),
     .     V(0:NQ,2),
     .     GRCH(3),RX(3),RAUX(2*NPOINT+1)


      REAL(DP) cons, qmax, rmax, delt, c, dlt, z1, z2, ch1, ch2, pi
      REAL(DP) r, vd, vv1, vv2, energ1, energ2, bessph
      integer iz1, iz2, itb, nr, nmin, nmax, nn, iq, ir

      REAL(DP) QTMP             

      PI= 4._DP * ATAN(1._DP)       
      CONS= 1.0_dp/(2.0_dp*PI)**1.5_DP
C     
C***  CUT-OFF IN REAL AND RECIPROCAL SPACE**
C     
      QMAX =  SQRT( Q2CUT )
      RMAX = PI * NQ / QMAX
      IF(RMX.GT.RMAX) THEN  
         if (IOnode) then
            WRITE(6,*) 'CH_OVERLAP: THE NUMBER OF INTEGRATION',
     .           ' POINTS MUST BE INCREASED'
            write(6,'(a,2f15.6)') 'ch_overlap: rmx,rmax =', rmx, rmax
         endif
         call die
      ENDIF 
      DELT=PI/QMAX
      C=4.0_DP*PI*DELT
      DLT=RMX/(NTBMAX-1)

      IZ1=IZVALFIS(IS1)
      IZ2=IZVALFIS(IS2)

      Z1=0.0_DP
      Z2=0.0_DP

      RX(2)=0.0_DP
      RX(3)=0.0_DP 

      DO IR=0,NQ
         R=IR*DELT
         RX(1)=R
             
         CALL PSCH(IS1,RX,CH1,GRCH)
         CALL PSCH(IS2,RX,CH2,GRCH)

         CH(IR,1)=-CH1
         CH(IR,2)=-CH2

         Z1=Z1-C*CH1*R*R    
         Z2=Z2-C*CH2*R*R

      ENDDO
           
      IF((ABS(Z1-IZ1).GT.CHERR).OR.
     .     (ABS(Z2-IZ2).GT.CHERR)) THEN 
         if (IOnode) then
            WRITE(6,*) 'CH_OVERLAP: THE NUMBER OF INTEGRATION',
     .           ' POINTS MUST BE INCREASED'
            WRITE(6,*) 'CH_OVERLAP: Z1=',Z1,' IZ1=',IZ1
            WRITE(6,*) 'CH_OVERLAP: Z2=',Z2,' IZ2=',IZ2
         endif
         call die
      ENDIF

      DO IR=0,NQ
         CH(IR,1)=real(IZ1,dp)*CH(IR,1)/Z1
         CH(IR,2)=real(IZ2,dp)*CH(IR,2)/Z2
      ENDDO 
C
C     REAL SPACE INTEGRATION OF POISSON EQUATION
C          
          
      CALL NUMEROV(NQ,DELT,CH(0,1),V(0,1))
      CALL NUMEROV(NQ,DELT,CH(0,2),V(0,2))
           
      DO ITB=1,NTBMAX
         R=DLT*(ITB-1)
         NR=NINT(R/DELT)
         NMIN=MAX(0,NR-NPOINT)
         NMAX=MIN(NQ,NR+NPOINT)
         NN=NMAX-NMIN+1
         DO IR=1,NN
            RAUX(IR)=DELT*(NMIN+IR-1) 
         ENDDO 
         CALL RATINT(RAUX,V(NMIN,1),NN,R,VV1,VD)
         CALL RATINT(RAUX,V(NMIN,2),NN,R,VV2,VD)
 
         VTB(ITB,1)=VV1
         VTB(ITB,2)=VV2
      ENDDO 
         
C****FOURIER-TRANSFORM OF RADIAL CHARGE DENSITY****
C
      CALL RADFFT( 0, NQ, RMAX, CH(0,1), CH(0,1) )
      CALL RADFFT( 0, NQ, RMAX, CH(0,2), CH(0,2) )
C

CNEUTRALIZE CHARGE DENSITY FOR FOURIER-SPACE CALCULATION
C
      DO IQ=0,NQ
         R=IQ*QMAX/NQ
         CH1 = (CH(IQ,1)-IZ1*CONS)*CH(IQ,2)
         CH2=  (CH(IQ,2)-IZ2*CONS)*CH(IQ,1)
         CH(IQ,1) = CH1
         CH(IQ,2) = CH2
      ENDDO
C
C     THE ELECTROSTATIC ENERGY CORRECTION IS STORED IN 'CORR'
C  
      DO IR=1,NTBMAX

         R=DLT*(IR-1)
         ENERG1=0.0_dp
         ENERG2=0.0_dp


         DO IQ=0,NQ
            QTMP=IQ*QMAX/NQ
            QTMP=QTMP*R 
            ENERG1=ENERG1+BESSPH(0,QTMP)*CH(IQ,1)
            ENERG2=ENERG2+BESSPH(0,QTMP)*CH(IQ,2)
         ENDDO 

         ENERG1=ENERG1*QMAX/NQ
         ENERG2=ENERG2*QMAX/NQ
   
         ENERG2=ENERG2*4.0_DP*(2.0_dp*PI)**2
         ENERG1=ENERG1*4.0_DP*(2.0_dp*PI)**2
              
         ENERG1=-(ENERG1*R)-(IZ2*(VTB(IR,1)*R-IZ1))
         ENERG2=-(ENERG2*R)-(IZ1*(VTB(IR,2)*R-IZ2))
  
         CORR(IR)=0.5_DP*(ENERG1+ENERG2)

      ENDDO 

      END subroutine ch_overlap


      SUBROUTINE NUMEROV(NR,DELT,Q,V)
      integer, intent(in)  :: nr
      REAL(DP), intent(in)   :: delt
      REAL(DP), intent(in)   :: q(0:nr)
      REAL(DP), intent(out)  :: v(0:nr)

C   Being Q(r) a spherical charge density in a homogeneus radial mesh
C   with distance DELT between consecutive points, this routine returns
C   the electrostatic potential generated by this charge distribution.
C   Written by D. Sanchez-Portal, March 1997.

C   INTEGER NR      :    Number of radial points.
C   REAL(DP)  DELT    :    Distance between consecutive points.
C   REAL(DP)  Q(0:NR) :    Spherical charge density.
C   REAL(DP)  V(0:NR) :    Electrostatic potential at mesh points.

C   Qtot/r asimptotic behaviour is imposed.


      integer ir
      REAL(DP) pi, fourpi, qtot, r, cons

      PI=4.0_DP*DATAN(1.0_DP)
      FOURPI=4.0_DP*PI

C     NUMEROV ALGORITHM* 
C
      V(0)=0.0_DP
      V(1)=1.0_DP

      DO IR=2,NR
         V(IR)=2.0_DP*V(IR-1)-V(IR-2) - FOURPI*DELT**3*
     .        ( Q(IR)*IR+10.0_DP*Q(IR-1)*(IR-1)+Q(IR-2)*(IR-2) )/12.0_DP
      ENDDO 

C***CALCULATE TOTAL CHARGE***
   
      QTOT=0.0_DP
      DO IR=1,NR
         R=IR*DELT
         QTOT=QTOT+R*R*Q(IR)
      ENDDO
      QTOT=4.0_DP*PI*QTOT*DELT

C** FIXING QTOT/R ASIMPTOTIC BEHAVIOUR*

      CONS=(QTOT-V(NR))/(NR*DELT)
             
      DO IR=1,NR
         R=IR*DELT
         V(IR)=V(IR)/(IR*DELT)+CONS
      ENDDO 
      V(0)=(4.0_DP*V(1)-V(2))/3.0_DP

      END subroutine numerov

      end module basis_io












