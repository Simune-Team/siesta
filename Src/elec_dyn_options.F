! This is a legacy module from Daniel's old TDDFT code. 
! Most part of it is now redundant. 
! Needs cleaning up and a possible merger into siesta_options 
! Rafi Ullah February 2015.
 
       module elec_dyn_options 
       use files,   only : slabel
       use precision
       use fdf
       use m_fdf_global
       use m_steps 
#ifdef MPI
       use mpi
#endif

       implicit none 
CMAZ
	logical, save :: ClProj,stopacp,ewin_clproj,do_stps ! do_stps : do it after given num of steps ?
	logical :: doclproj,didclproj,chkmvdions
	integer :: iproj,nproj,projstp
	integer:: NStepProj ! NStepProjections; defualt=50
	logical, save :: wrtinstocc ! WriteInstOccupations; default=false
	double precision :: emin_clproj,emax_clproj  ! min and max energy for ClProj. 
	double precision :: dist
	double precision, save :: xa0,dxclproj ! xa0 :initial position of porjectile along vu axis
	integer, save :: vu ! unit vector of v of projectile, 1,2,3
	logical, save :: wrttdfiles ! to write dm,tdwf,xv& verlet_restart files at each iproj step 
                           ! (1 stp before Clproj so that it can be done in td-step).
	logical, save :: wrtSigmaZTflag	! whether to write SigmaZT is required by user.
	logical, save :: wrtSigmaZT ! whether to call SigmaZT to write Sigma(z,t), (if wrtSigmaZTflag =true).
	integer::	NStepSigmaZT
	integer::	iaProjectile
	logical ::	wrtxzSZmax
c maz ----

	double precision :: twall,twmax
	integer :: twmaxh,twmaxm

c maz. get pdos:
	logical, save :: getpdos
	integer, save :: pdos_strt_stp,pdosati,nstppdos
	logical :: writepdos

c maz. limit wall time 
	logical, save :: limwalltime

!	maz
	integer,save :: istep77
!	maz . write Force_vs_time file
	logical, save :: fvst
	real(dp),save :: kin_host,kin_proj


       logical electrondynamics
       logical onlyelectrons
       logical calculateEnew
       logical calculateforces
       logical movedions
       logical moveions
       logical response
       integer responsetype
       logical EigenvsTime
       logical EtotvsTime
       logical ElecDynPC
       logical usesaveTDWF, usesaveTDWF_default
       double precision timestep, restarttime
       integer itimestep, finalstep, initstep
       integer nstep_ion_move, nstep_ion_move_default
       character*70 dipolefile
       character*70 eigenfile
       character*70 etotfile
       integer ntddftwrite, ntddftwrite_default
       integer nelecsubsteps, nelecsubsteps_default
#ifdef MPI
      integer
     .  MPIerror
#endif
      
       CONTAINS

       subroutine read_dynoptions( )
       implicit none 
       integer node 
       character
     .   dyntyp*22, resp*22
      character
     .   dyntyp_default*22, resp_default*22
!Rafi
!      logical leqi
!End-Rafi    
   

CZEB 
      call fdf_global_get(wrtinstocc,"WriteInstOccupations",.false.)

	if(wrtinstocc) then !*************************

      call fdf_global_get(NStepProj,"NStepProjections",50)
c      call fdf_global_get(NenrgyDiv,'EnrgyDivisions'  , 3000)	! defined in siesta.f now.
      call fdf_global_get(stopacp,"StopAfterCalProj",.false.) ! stop after cal proj

      call fdf_global_get(ewin_clproj,"EWindowCalProj",.false.)

	if(ewin_clproj) then
      	call fdf_global_get(emin_clproj,"EminCalProj",0.0d0)
      	call fdf_global_get(emax_clproj,"EmaxCalProj",0.0d0)
	if(emin_clproj .ge. emax_clproj)
     .				Stop "emin_clproj .ge. emax_clproj"
	end if

      call fdf_global_get(do_stps,"DoNStepProj",.true.)

	if(.not. do_stps) then
	call fdf_global_get(xa0,"ProjInitPosition",0.0d0)
			 ! cartesian coor (x, y or z) along the velocity.
	call fdf_global_get(vu,"ProjVdirection",3)
			 ! 1,2,3 for x,y,z axis.
	if(vu .gt. 3 .or.  vu .lt. 1) then
	vu = 3
	write(6,*) 'ProjVdirection should be 1,2 or 3: set = 3'
	end if

	call fdf_global_get(dxclproj,"ProjDistanceInterval",1.0d0) 
			 ! cartesian coor (x, y or z) along the velocity.

	end if

	end if ! wrtinstocc !*************************




      call fdf_global_get(wrtSigmaZTflag,"WriteArealChargeDensity",
     .        .false.)

	if (wrtSigmaZTflag) then

      call fdf_global_get(NStepSigmaZT,"NStepArealChargeDensity",
     .        3)

      call fdf_global_get(iaProjectile,"iProjectileInListOfAtoms",
     .       1)

      call fdf_global_get(wrtxzSZmax,"WriteMaxACDensityPosition",
     .       .false.)

	end if


c maz.
      call fdf_global_get(getpdos,"GetPDOS",.false.)
	if (getpdos) then
      call fdf_global_get(pdos_strt_stp,'PDOSStartStep', 1) 
      call fdf_global_get(pdosati,"pdosAtomi",
     .       1)
      call fdf_global_get(nstppdos,"NStepPDOS", 10)

	end if


c maz
	call fdf_global_get(limwalltime,"WalltimeLimit",.false.)

	if (limwalltime) then
	 call fdf_global_get(twmaxh,"WalltimeHours", 0)
	 call fdf_global_get(twmaxm,"WalltimeMinutes", 0)
	 twmax = 3600.0d0*twmaxh + 60.0d0*twmaxm
	end if

! maz
!	write Forces_vs_time file
      call fdf_global_get(fvst,"WriteForcesFile",.false.)







       dyntyp_default = 'verlet'
       ntddftwrite_default=100
       Nelecsubsteps_default=1
       nstep_ion_move_default=20
c       call fdf_global_get(dyntyp,'MD.TypeOfRun',dyntyp_default)
c      if (leqi(dyntyp,'electrondynamics')) then
C         electrondynamics = .true.
C Number of electronic timesteps between ionic movements
          call fdf_global_get(nstep_ion_move,'MD.Nstep_ion_move',
     .      nstep_ion_move_default) 
          if(nstep_ion_move.lt.3) nstep_ion_move=3 
C Number of "substeps" larger or equal than 1
          call fdf_global_get(Nelecsubsteps,'MD.Nelecsubsteps', 
     .      Nelecsubsteps_default)
          if(Nelecsubsteps.lt.1) Nelecsubsteps=1
          call fdf_global_get(ElecDynPC,'MD.ElecDynPC',.false.)
          call fdf_global_get(EigenvsTime,'WriteEigenvsTime', .false. )
          call fdf_global_get(onlyelectrons,'MD.OnlyElectrons', .false.)
          usesaveTDWF_default=.false.
          call fdf_global_get(usesaveTDWF,'MD.UseSaveTDWF',
     .         usesaveTDWF_default)

          if(onlyelectrons) then      
           resp_default=' '
           response=.false. 
           if (Node.eq.0) then
             call fdf_global_get(resp,'MD.ResponseType',resp_default)
             if (leqi(resp,'molecule')) then
                 responsetype=1
                 response=.true. 
             else if (leqi(resp,'bulk')) then
                 responsetype=2
                 response=.true.
                write(6,*) 'Bulk optical response calculation'
                write(6,*) 'Not yet implemented, sorry'
                stop 
             else if (leqi(resp,' ')) then
                 responsetype=0
            endif
          endif
         endif

c      endif

      end subroutine read_dynoptions


C Returns total simulations time
C      double precision function tottime
C       function tottime(itimestep, initstep, timestep, restarttime)
C
C       double precision tottime, timestep, restarttime
C       integer itimestep, initstep
C
C       tottime= (itimestep-initstep)*timestep + restarttime 
C       return 
C       end function tottime



     
       subroutine  iodipole (dipole,extfield,frstme)
       double precision dipole(3), extfield(3) 
       integer iu
       character*15 fform
       logical frstme
       save iu 

       character*33 paste
      if(frstme) then
        dipolefile = paste (slabel, '.dipol_vs_time') 
        call io_assign( iu )
        fform='formatted'
        open( iu, file=dipolefile, form=fform, status='unknown' )
        write(iu,'(a,3f15.6)') '#',extfield(1), extfield(2), extfield(3)
      endif
      if(itimestep.gt.initstep) then 
       write(iu,'(4f15.6)') 
C     .    totime,
     .    dipole(1),dipole(2), dipole(3)
      endif
      if(itimestep.eq.finalstep) call io_close(iu)

      end subroutine iodipole

      subroutine ioetot(totime, etot, istp, itd, ntd,rstart_time)

      double precision totime, etot, eV, rstart_time
      logical frstme
      integer iu, istp, itd, ntd
       character*33 paste
      character*15 fform
      data frstme /.true./
      save frstme, iu, eV
       if(frstme) then
          etotfile = paste (slabel, '.etot_vs_time')
          call io_assign( iu )
          fform='formatted'
          open( iu, file=etotfile, form=fform,POSITION='APPEND',
     .         status='unknown' )
          frstme = .false. 
          eV=1.d0 / 13.60580d0
       end if
       if (totime .gt. rstart_time) then
          write(iu,'(2f15.6)') 
     .         totime, etot/eV
       end if

       if (istp .eq. fincoor .and. itd .eq. ntd) call io_close(iu)

      end subroutine ioetot 
 
      subroutine ioeigenvalues (eigen,maxo,maxspn,nocc,nspin,nk, totime,
     .                            rstart_time, istp, itd,ntd)   
 
      integer maxo, maxspn,nspin, nk, istp, itd, ntd
      double precision eigen(maxo,maxspn,nk)
      double precision eV , totime, rstart_time
      integer nocc(nk,nspin)
      integer ik,ispin,ie, iu
      logical frstme
      character*15 fform
       character*33 paste
      data frstme /.true./
      save frstme, iu, eV
       
        if(rstart_time .eq. 0.00_dp) then
          eigenfile = paste (slabel, '.eigen_vs_time')
          call io_assign( iu )
          fform='formatted'
          open( iu, file=eigenfile, form=fform, status='unknown' )
          write(iu,*) '# ',nspin,nk
          do ispin=1,nspin
             write(iu,*) '# ',(nocc(ik,ispin),ik=1,nk)
          enddo 
          eV=1.d0 / 13.60580d0
          frstme=.false.
        endif
        
        if (totime .gt. rstart_time) then
           write(iu,'(11f12.5,/,(5x,10f12.5))') 
     .        totime,
     .        (((eigen(ie,ispin,ik)/eV,ie=1,nocc(ik,ispin)),
     .        ispin=1,nspin),ik=1,nk)
        end if
       if (istp .eq. fincoor .and. itd .eq. ntd) call io_close(iu)
      end subroutine ioeigenvalues

! maz kin_host,kin_proj
      subroutine kinion(natoms,ma,va,cfa,dt, inst_kinion,
     .	kin_host1,kin_proj1)

      
      real(dp) ma(natoms), va(3,natoms), cfa(3,natoms)
      real(dp) kin, fovermp, Ang, eV, v, inst_kinion, dt,
     . 	kin_host1,kin_proj1 ! maz kin_host,kin_proj
      real(dp), pointer, save:: cfaold(:,:)
      logical first
      integer ia,ix, natoms
      data first /.true./

      Ang    = 1._dp / 0.529177_dp
      eV     = 1._dp / 13.60580_dp
      fovermp = 0.009579038_dp * Ang**2 / eV
   

	write(6,*)'£££ elec_dyn_options: kinion: dt = ',dt


      if (first) then 

!	write(6,*)"&& elec_dyn_options: kinion: "
!	write(6,*)"&&  line 342: first = ", first

         allocate(cfaold(3,natoms))
      kin=0.0d0
!      do ia=1, natoms
!         do ix=1,3
!          v=va(ix,ia)
!          kin= kin + 0.5_dp * ma(ia) * v**2 / fovermp
!         enddo
!      enddo

! maz
          kin_host1 = 0.0d0
      do ia=1, natoms-1 ! projectile is last atom
         do ix=1,3
          v=va(ix,ia)
          kin_host1= kin_host1 + 0.5_dp * ma(ia) * v**2 / fovermp
         enddo
      enddo

          kin_proj1 = 0.0d0
	! projectile is last atom
         do ix=1,3
          v=va(ix,natoms)
          kin_proj1= kin_proj1 + 0.5_dp * ma(natoms) * v**2 / fovermp
         enddo

	kin=kin_host1+kin_proj1

      inst_kinion=kin

!	maz
!	 	cfaold are not defined at the second step.
!		make them = cfa so, no need to change things there.
           do ia=1,natoms
             do ix=1,3
	       cfaold(ix,ia) = cfa(ix,ia)
	     end do
	   end do


      first=.false.
       return
      endif

C  With velocity verlet we do not need the Forces at
C  a given step to calculate velocities in such step.
C  However, we do not know such "new" forces until 
C  they are calculated.  Therefore, in order to print
c  the kinetic energy at each electron dynamic step, 
C   

      kin=0.0d0
      do ia=1, natoms
         do ix=1,3
          v=va(ix,ia)
     .   +0.5d0*dt*fovermp*(cfaold(ix,ia)+cfa(ix,ia))/ma(ia)
          kin= kin + 0.5_dp * ma(ia) * v**2 / fovermp
         enddo 
      enddo 
      inst_kinion=kin


          kin_host1 = 0.0d0
      do ia=1, natoms-1 ! projectile is last atom
         do ix=1,3
          v=va(ix,ia)
     .   +0.5d0*dt*fovermp*(cfaold(ix,ia)+cfa(ix,ia))/ma(ia)
          kin_host1= kin_host1 + 0.5_dp * ma(ia) * v**2 / fovermp
         enddo 
      enddo 

          kin_proj1 = 0.0d0
	! projectile is last atom
         do ix=1,3
          v=va(ix,natoms)
     .   +0.5d0*dt*fovermp*(cfaold(ix,natoms)+cfa(ix,natoms))/ma(natoms)
          kin_proj1=kin_proj1 + 0.5_dp *ma(natoms)* v**2 /fovermp
         enddo 
	     
	kin=kin_host1+kin_proj1

      inst_kinion=kin

      if(calculateforces) then 
!	write(6,*)'&& calculateforces = ', calculateforces
         do ia=1,natoms
             do ix=1,3
                 cfaold(ix,ia)=cfa(ix,ia)
             enddo 
         enddo 
      endif

     
      end subroutine kinion
        
      end module elec_dyn_options

