      subroutine read_user_basis

#ifndef CDF

      call die(
     $ '*** You need netCDF to read the new user-defined basis files...'
     $ )
      end subroutine  read_user_basis
#else 

      use precision
      use types
      use atmfuncs, only: nspecies, species, npairs, elec_corr
      use atmfuncs, only: floating, izvalfis, psch
      use radial
      use netcdf
      use fdf
      use ionew, only: IOnode

      implicit none

      type(species_info), pointer        :: spp
      type(rad_func), pointer            :: op
      type(rad_func), pointer            :: pp
      type(rad_func), pointer            :: func

      integer ncid, iret

      integer nkbs, nkbs_id, ntb_id, proj_id,
     $     nkbl_id, pjnl_l_id, pjnl_n_id, pjnl_ekb_id, kbdelta_id,
     $     kbcutoff_id
      integer norbs, norbs_id, orbnl_l_id, orbnl_n_id, orbnl_z_id,
     $     cutoff_id, delta_id, orb_id, orbnl_pop_id, orbnl_ispol_id
      integer vlocal_id, chlocal_id, core_id

      integer aux(maxnorbs)

      integer is, j, i, l, nrp_tables, core_flag, is2, nor, nk, m
      integer ia, noa, nkba
      real(dp) rchloc, rchloc2

      character(len=20) filename
      character(len=20) dummy
                                !
                                !     Automatic arrays
                                !
      character(len=20)              :: species_name(20)
      integer                        :: z(20)

      nspecies = fdf_integer('NumberOfSpecies',0)
      if (nspecies.eq.0) call die("Please specify number of species!!")
      if (nspecies.GT.20)
     $     call die("read_user_basis: Please increase bounds!!")

      allocate(species(nspecies))

      call get_species_names(nspecies,species_name,z)

      do is = 1, nspecies
         spp => species(is)
         spp%label = species_name(is)
         write(filename,'(a,a)') trim(species_name(is)), ".basis.nc"

         iret = nf90_open(trim(filename),NF90_NOWRITE,ncid)

         iret = nf90_inq_dimid(ncid,'norbs',norbs_id)
         iret = nf90_inquire_dimension(ncid,norbs_id,len=norbs)
         if (norbs .gt. maxnorbs)
     $        call die("read_user_basis: Increase maxnorbs in types.f")

         spp%n_orbnl = norbs

         iret = nf90_inq_dimid(ncid,'nkbs',nkbs_id)
         iret = nf90_inquire_dimension(ncid,nkbs_id,len=nkbs)
         spp%n_pjnl = nkbs

         iret = nf90_inq_dimid(ncid,'ntb',ntb_id)
         iret = nf90_inquire_dimension(ncid,ntb_id,len=nrp_tables)
         if (nrp_tables .ne. nrtmax) call die("nrtmax snafu")

         allocate(spp%orbnl(norbs))
         allocate(spp%pjnl(nkbs))


         iret = nf90_get_att(ncid,nf90_global,'Element',spp%symbol)
!!!        iret = nf90_get_att(ncid,nf90_global,'Label',dummy)
                                !! Sanity check here??

         iret = nf90_get_att(ncid,nf90_global,'Atomic_number',spp%z)
         if (z(is) .ne. spp%z) call die("Atomic number mismatch")

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
         iret = nf90_inq_varid(ncid,'vlocal',vlocal_id)
         iret = nf90_get_att(ncid,vlocal_id,
     $        'Vlocal_cutoff',spp%vlocal%cutoff)
         iret = nf90_get_att(ncid,vlocal_id,
     $        'Vlocal_delta',spp%vlocal%delta)
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
            spp%core%cutoff = 0._dp
            spp%core%delta = 0._dp
            spp%core%f = 0._dp
            spp%core%d2 = 0._dp
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
            iret = nf90_get_var(ncid,proj_id,pp%f(1:),
     $           start=(/1,i/),count=(/nrtmax,1/))
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
         iret = nf90_get_var(ncid,vlocal_id,spp%vlocal%f(1:),
     $        start=(/1/),count=(/nrtmax/))
         call check(iret)
         call rad_setup_d2(spp%vlocal)

                                !
                                !       Local potential charge density
                                !
         iret = nf90_get_var(ncid,chlocal_id,spp%chlocal%f(1:),
     $        start=(/1/),count=(/nrtmax/))
         call check(iret)
         call rad_setup_d2(spp%chlocal)

         if (spp%there_is_core) then
            iret = nf90_get_var(ncid,core_id,spp%core%f(1:),
     $           start=(/1/),count=(/nrtmax/))
            call check(iret)
            call rad_setup_d2(spp%core)
         endif

         nor = 0
         do i = 1, norbs
            op => spp%orbnl(i)
            iret = nf90_get_var(ncid,orb_id,op%f(1:),
     $           start=(/1,i/),count=(/nrtmax,1/))
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
!     Need to update the orb_XX and pj_XX lists and indexes...
!
      
      npairs = ((nspecies+1)*nspecies)/2
      allocate(elec_corr(npairs))

      do is=1,nspecies
         rchloc = species(is)%chlocal%cutoff
         do is2=is,1,-1
            rchloc2 = species(is2)%chlocal%cutoff
            i = ((is-1)*is)/2+is2
            func => elec_corr(i)
            if (floating(is) .or. floating(is2)) then
               call rad_zero(func)
            else
               func%cutoff = rchloc + rchloc2 + 0.2_dp
               func%delta =  func%cutoff / (nrtmax - 1)
               call ch_overlap(is,is2,func%cutoff,func%f)
               call rad_setup_d2(func)
            endif
         enddo
      enddo


      CONTAINS

      subroutine check(status)
      
      integer, intent(in):: status
      if (status .ne. nf90_noerr) then
         print  *, trim(nf90_strerror(status))
         stop 'Stopped'
      endif
      end subroutine check

      subroutine get_species_names(ns,atm_label,z)

      use fdf
      use parsing

c     Reads the 'long_names' for different species. 
c     The labels are neccesary to 
c     identify the files which contain the information about the 
c     basis.
c     
c     Reads fdf block. 
c     
      implicit          none

      integer, intent(in)  ::    ns ! Number of species
      integer, intent(out) ::    z(:) ! Atomic numbers
      character(len=*)     ::    atm_label(:) ! Names for the species


      character         line*130, names*80
      integer           ni, nn, nr, nv, ist, nst, iu, nsread
      integer           integs(4), lastc, lc(0:3)
      double precision  reals(4), values(4)

c     check for block and read

      if ( fdf_block('Chemical_species_label',iu) ) then

         nst = 0
         nsread=0
         do ist = 1, ns
            read(iu,'(a)', end=50) line
            lastc = index(line,'#') - 1
            if (lastc .le. 0) lastc = len(line)
            call parse( line(1:lastc), nn, lc, names, nv, values,
     .           ni, integs, nr, reals )
            if (ni .ge. 2)  then
               nst=nst+1
               if ((integs(1).gt.ns) .or. (integs(1).lt.1) ) then
                  call die("Bad spec. number in Chemical Species block")
               else
                  nsread=nsread+1
                  atm_label(integs(1)) = names(lc(0)+1:lc(1))
                  z(integs(1)) = integs(2)
               endif 
            else
               call die("Bad format in Chemical Species block")
            endif
         enddo
 50      continue
         if (nsread .ne. ns) call die("Species missing in block...")
      else
         call die("Cannot find Species block...")
      endif
      end subroutine get_species_names

      SUBROUTINE CH_OVERLAP(IS1,IS2,RMX,CORR)
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

      real(dp)            :: aux(size(corr))   ! automatic array

      integer nq, npoint
      real(dp)            :: q2cut, cherr
      PARAMETER ( NQ     =  512  )
      PARAMETER ( NPOINT =  4     ) 
      PARAMETER ( Q2CUT  =  2.5e3_dp )
      PARAMETER ( CHERR   =  5.e-2_dp )

      real(dp)
     .     CH(0:NQ,2),VTB(NRTMAX,2),
     .     V(0:NQ,2),
     .     GRCH(3),RX(3),RAUX(2*NPOINT+1)


      REAL(DP) cons, qmax, rmax, delt, c, dlt, z1, z2, ch1, ch2, pi
      REAL(DP) r, vd, vv1, vv2, energ1, energ2, bessph, dev1, devn
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
      DLT=RMX/(NRTMAX-1)

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
           
      DO ITB=1,NRTMAX
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
      DO IR=1,NRTMAX

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

      end subroutine read_user_basis

#endif







