#ifndef CDF
      subroutine dump2cdf
      write(6,'(/,a,/)')
     $     '*** Install netCDF and get a nice dump of atomic info...'
      end subroutine dump2cdf
#else
      subroutine dump2cdf

      use types
      use atmfuncs
      use netcdf

      implicit none

      integer ncid, iret

      integer nkbs, nkbs_id, ntb_id, proj_id,
     $        nkbl_id, pjnl_l_id, pjnl_n_id, pjnl_ekb_id, kbdelta_id,
     $        kbcutoff_id
      integer norbs, norbs_id, orbnl_l_id, orbnl_n_id, orbnl_z_id,
     $        cutoff_id, delta_id, orb_id, orbnl_pop_id, orbnl_ispol_id
      integer vlocal_id, chlocal_id, core_id

      integer aux(maxnorbs)

      integer is, j, i, l
      character*20 filename

      write(6,'(/,a,/)')
     $     '*** Dumping atomic info in netCDF format (.nc files)'

      do is = 1, nspecies
         spp => species(is)
         write(filename,'(a,a)') trim(spp%label), ".dump.nc"
         write(6,*) 'Dumping to ', filename
         write(6,*) 'Trimmed: ', trim(filename)

         iret = nf90_create(trim(filename),NF90_CLOBBER,ncid)

        nkbs =  spp%n_pjnl
        norbs = spp%n_orbnl

        iret = nf90_def_dim(ncid,'norbs',norbs,norbs_id)
        call check(iret)
        iret = nf90_def_dim(ncid,'nkbs',nkbs,nkbs_id)
        call check(iret)
        iret = nf90_def_dim(ncid,'ntb',nrtmax,ntb_id)
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
        iret = nf90_def_var(ncid,'vlocal',nf90_double,
     $                      (/ntb_id/),vlocal_id)
        iret = nf90_put_att(ncid,vlocal_id,
     $                      'Vlocal_cutoff',spp%vlocal%cutoff)
        iret = nf90_put_att(ncid,vlocal_id,
     $                      'Vlocal_delta',spp%vlocal%delta)
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

        if (norbs .gt. maxnorbs) call die("dump2cdf: Increase maxnorbs")
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
     $                      start=(/1,i/),count=(/nrtmax,1/))
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
        iret = nf90_put_var(ncid,vlocal_id,spp%vlocal%f(1:),
     $                      start=(/1/),count=(/nrtmax/))
        call check(iret)

!
!       Local potential charge density
!
        iret = nf90_put_var(ncid,chlocal_id,spp%chlocal%f(1:),
     $                      start=(/1/),count=(/nrtmax/))
        call check(iret)

        if (spp%there_is_core) then
           iret = nf90_put_var(ncid,core_id,spp%core%f(1:),
     $                      start=(/1/),count=(/nrtmax/))
           call check(iret)
        endif

        do i = 1, norbs
           op => spp%orbnl(i)
           iret = nf90_put_var(ncid,orb_id,op%f(1:),
     $                      start=(/1,i/),count=(/nrtmax,1/))
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

      end subroutine dump2cdf

#endif







