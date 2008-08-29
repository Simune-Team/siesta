! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module atom

      use precision, only: dp
      use sys,       only: die
      use atmparams
      use m_recipes, only: spline, polint
      use basis_types, only: basis_parameters
      use basis_types, only: ground_state_t

!----------------------------------------------------------------
!     old_atmfuncs arrays
!
      use old_atmfuncs, only: tabpol, rcpoltb, table, tab2
      use old_atmfuncs, only: coretab, smasstb, tab2pol
      use old_atmfuncs, only: zvaltb, chargesave, nkbmax, nomax
      use old_atmfuncs, only: qtb, slfe
      use old_atmfuncs, only: chloctab, vlocaltab
      use old_atmfuncs, only: lmxosave, lmxkbsave, label_save
      use old_atmfuncs, only: semicsave, izsave, cnfigtb
      use old_atmfuncs, only: nzetasave, nsemicsave, nkblsave
      use old_atmfuncs, only: npolorbsave, basistype_save
      use old_atmfuncs, only: lambdatb, rcotb, rctb
!
!     old_atmfuncs procedures
!
      use old_atmfuncs, only: labelfis, izofis
!----------------------------------------------------------------
      use periodic_table
      use basis_types, only: basis_def_t
      use fdf
      use pseudopotential, only: pseudopotential_t

      implicit none      

      private

      public :: atom_main, prinput

!---------------------------------------------------------------------
! Some internal parameters

        real(dp), parameter             :: deltmax=0.05d0
! Maximum distance (in Bohrs) between points where function 
! is evaluated to generate the tables for interpolation
! If this distance is exceeded  a warning is printed
! In practice, this means that the cutoff can be as large as 25 bohrs
! with ntbmax=500

        real(dp), parameter             :: eshift_default=0.02d0
! Default energy-shift to define the cut off radius of orbitals
! In Rydbergs

        character(len=1)     :: sym(0:4) = (/ 's','p','d','f','g' /)

!---------------------------------------------------------------------
       
      CONTAINS

        SUBROUTINE ATOM_MAIN(IZIN,LMXKB_IN,NKBL,EREFKB,LMXO, NZETA, 
     .          RCO,LAMBDA_IN,ATM_LABEL,
     .          NPOLORB,SEMIC,NSEMIC,CNFIGMX,CHARGE_IN,SMASS,BASISTYPE, 
     .          ISIN,RINN,VCTE,split_norm,basp)

        implicit none

        type(basis_def_t), pointer   :: basp

      integer, intent(in) :: izin  
!       Atomic number (if IZ is negative it will be regarded as a set of 
!       floating orbitals rather than a real atom)

      integer, intent(IN)  :: lmxkb_in
!       Angular momentum cutoff for Kleinman-Bylander nonlocal pseudopotential

      integer, intent(in) :: nkbl(0:lmaxd) 
!       Number of Kleinman-Bylander projectors for each angular momentum

      real(dp), intent(in) :: erefkb(nkbmx,0:lmaxd)
!       Reference energies (in Ry) for the calculation  of the KB projectors

      integer, intent(IN) :: lmxo 
!       Angular momentum cutoff for basis orbitals 
!       (not including the polarization base functions)

      integer, intent(in) :: nzeta(0:lmaxd,nsemx)
!       Number of shells of PAOs with the same angular momentum. 
!       (i.e.not including the polarization base functions)

      real(dp), intent(IN) :: rco(nzetmx,0:lmaxd,nsemx)
!       Cutoff radius for Sankey-type basis orbitals

      real(dp), intent(IN) :: lambda_in(nzetmx,0:lmaxd,nsemx)
!       Scaling factor for the atomic orbital,  for the 'SPLIT' type 
!       of basis it is interpreted as the exponent (in Bohr-2)
!       of the gaussian orbitals for the double-Z, triple-Z, etc.. 

      character(len=20), intent(in) :: atm_label
!       Label used to name all the atom related files (i.e. pseudopotential
!       files, in/output files containing the basis set, etc..) 

      integer, intent(in) :: npolorb(0:lmaxd,nsemx)
!       Number of shells of polarization orbitals for pseudo-atomic
!       orbitals of momentum L.

      logical, intent(in) :: semic
!       If .true. there is a shell of semicore states with angular
!       momentum NSEMIC

      integer, intent(in) :: nsemic(0:lmaxd)
!       Number of semicore states with a given angular momentum

      integer, intent(in) :: cnfigmx(0:lmaxd)
!       Maximum principal quantum number present for a given angular momentum

      real(dp), intent(IN)  :: charge_in
!       Charge state of the atom, only for basis set generation purposes
!       AS SPECIFIED IN PAO.BASIS common block

      real(dp), intent(in)  :: smass
!       Atomic mass for the species

      character(len=10), intent(in)  :: basistype
!       Type of augmentation procedure for the construction of the basis set

      integer, intent(in)  :: isin
!       Species index 

      real(dp), intent(in) :: vcte(0:lmaxd,nsemx), rinn(0:lmaxd,nsemx)
!       Parameters for soft-confinement potential

      real(dp), intent(in) :: split_norm(0:lmaxd,nsemx)
!      User-defined split_norm values
!
!     Former arguments, no longer used
!
      integer no, nkb          ! total num of orbitals and KB projectors

!     Extra copies to avoid  problems with the old 'inout'
!
      integer lmxkb
      real(dp) charge
      real(dp) :: lambda(nzetmx,0:lmaxd,nsemx)

!-----------------------------------------------------------------------
!   Initializes Kleinman-Bylander pseudopotential and atomic orbitals.
!   Atomic orbitals basis sets are of several types.
!   **PAO's orbitals (Pseudo-atomic orbitals). We can form several
!    types of basis sets depending on how we double the basis set:
!     1) 'NODES': Double-Z, Triple-Z,... are orbitals of PAO 
!         type with one, two,... nodes
!     2) 'SPLIT': Double-Z, Triple-Z,... are generated from the 
!         original orbital, being equal to the PAO outside a 
!         radius Rm and given by (a*r^2+b)*r^l inside Rm.
!     3) 'NONODES': Double-Z, Triple-Z,... are orbitals of PAO
!         type with no-nodes, but different radii or different
!         scale-factors.
!     4) 'SPLITGAUSS':  Double-Z, Triple-Z,... are orbitals of GTO
!         (Gaussian Type Orbitals).
!     5) 'USER': Read from a file, provided by the user.
! 
!   Written by D.Sanchez-Portal.   1996-1998.
!   New version.   August, 1998. 
!   Modified by DSP to allow more than one KB projector per l, July 1999
!   Modified by DSP tp allow several shells of semicore states, July 1999
!   Rewritten by Alberto Garcia, October 1999
!                                       
!     Distances in Bohr.
!     Energies in Rydbergs.

! BEHAVIOUR
!   1) Several calls with the same input return the same output 
!    (in particular the same species index IS)
!   2) When called for the same atomic number but some other different
!    input, a different species index is assigned and a warning is 
!    printed.
!   3) IZ<0 is accepted, returning the orbitals of the corresponding 
!    atoms
!    but zero pseudopotential. If IZ=-100 square-box orbitals and zero 
!    pseudopot.
!   4) IZ=0 is also accepted, producing the re-initialization of tables,
!    i.e. setting to zero the number of tabulated tables. If this is 
!    done care must be taken to re-initialize also all routines which 
!    depend on RCUT, PHIATM, VLOCAL, etc... (like MATEL in O(N) program).
! 
! 


           
! Internal variables
 
       real(dp)
     .  rofi(nrmax), drdi(nrmax), s(nrmax),
     .  vps(nrmax,0:lmaxd), rphi(nrmax,0:lmaxd,nsemx),
     .  vlocal(nrmax), vxc(nrmax), ve(nrmax),
     .  rho(nrmax), chcore(nrmax), rho_PAO(nrmax), auxrho(nrmax),
     .  vePAO(nrmax), qPAO(0:lmaxd,nsemx), chlocal(nrmax),
     .  red_vlocal(nrmax)


       real(dp) 
     .  pi, a, b, zval, flting,
     .  ex, ec, dx, dc, r2, chgvps
     
      
       real(dp)
     .  Rgauss, Rgauss2, Rchloc

       character
     .   icorr*2, irel*3, nicore*4 

       integer 
     .  iz, nrval, ir , nrgauss, nchloc, nzcontr, l, nVna,
     .  irelt, nsm, nvlocal
    
       logical  new
   
         INTEGER is

         IS = ISIN
         pi=acos(-1.0d0)
! 
! Print some information about the atomic species
! and selected options
! 
      iz=izin
      lmxkb = lmxkb_in
      charge = charge_in
      lambda(:,:,:) = lambda_in(:,:,:)

!!** AG: Symbol is not adequate. Should use the label

             if (iz.gt.200) then  
               write(6,'(3a,i4,a)')
     .         'atom: Called for (synthetic) ', atm_label,
     $               '  (Z =', iz,')' 

             else if (iz.gt.0) then  
               write(6,'(3a,i4,a)')
     .         'atom: Called for ', atm_label, '  (Z =', iz,')' 

             elseif((iz.lt.0).and.(iz.ne.-100)) then  

               write(6,'(3a,i4,a,a)')
     .         'atom: Called for ', atm_label, '  (Z =', iz,')',
     .         ' ( Floating basis ) ' 

             elseif(iz.eq.-100) then
               write(6,'(a,i4,a)')
     .         'atom: Called for Z=',iz,'( Floating Bessel functions)'  

             endif
! 
! 
           
           call new_specie(iz,lmxkb, 
     .         nkbl, lmxo,
     .         nzeta, atm_label,
     .         npolorb, semic, nsemic, cnfigmx,
     .         is, new, no, nkb) 
! 
              basistype_save(is)=basistype
              smasstb(is)=smass

        if (iz.ne.-100) then

! Reading pseudopotentials*** 
! 
            call read_vps(lmxo, lmxkb,
     .        nrval,a,b,rofi,drdi,s,vps,
     .        rho, chcore, zval, chgvps, nicore, irel, icorr,basp) 

            do ir=1,nrval
              rho(ir)=chgvps*rho(ir)/zval
            enddo 
!
!           Rho now contains the 'true' charge in the pseudopotential
!           calculation, as rho in the VPS file is rescaled to the
!           charge of a neutral atom.
!
            if (abs(zval-chgvps).gt.1.0d-3) then 
              write(6,'(/,a,/,a,f5.2)') 
     .  'atom: Pseudopotential generated from an ionic configuration',
     .  'atom: with net charge', zval-chgvps
            endif 
!
!           AG: Note zval-chgvps = (Znuc-Zcore)-gen_zval
!           Example: Ba with 5s and 5p in the valence:
!           zval from atom (scaled with zratio): 10 (5s2 5p6 6s2)
!           chgvps (only 5s2 and 5p6 occupied in ref config): 8
!           gen_zval = chgvps
!           ==> net charge (zval-chgvps) = 2
!           Znuc = 56,  Zcore= Z(Xe)-2-6 = 54-2-6 = 46
!           Znuc-Zcore-true_val = 56 - 46 - 8 = 2
!
!           If we stop the practice of scaling the valence charge
!           in the ps file, we would need info about Zcore to setup
!           things correctly (both for Siesta and for the PW program...)
!           BUT actually, we have that information in the ps file!!
!           The potentials are stored as r*V, so at large r they all
!           go to -2*(Znuc-Zcore).....
!
!           Set 'charge':
!              1. If 'charge' is not set in the fdf file
!                 then set it to zval-chgvps.
!              2. If 'charge' is equal to zval-chgvps, set it to that.
!
            if ((abs(charge).eq.huge(1.0_dp)).or. 
     .          (abs(charge-(zval-chgvps)).lt.1.0d-3)) then   
c             write(6,'(/,2a)') 
c    .          'atom: The above configuration will be used ', 

              charge=zval-chgvps
            endif
! 
!  
!  Save read valence charge
!  This can be different from the standard one if we have included semicore
!  states in the valence shell.
!  For example: Ba with 5s and 5p in the valence: zval=10 (5s2 5p6 6s2)
! 
            zvaltb(is)= zval

! IF IZ IS NEGATIVE WE WILL HAVE FLOATING ORBITALS IN THE ATOMIC POSITIONS
! IF IZ POSITIVE WE WILL HAVE REAL ATOMS
! 
           flting=dsign(1.0d0,dble(iz))
           iz=abs(iz)
! 
! Common block with pseudocore information****
! used for non-linear-core-correction for exchange-correlation energy***
! 
           call comcore(is,a,b,rofi,chcore,nrval,nicore,flting)

! CALCULATION OF THE VALENCE SCREENING POTENTIAL FROM THE READ CHARGE
! DENSITY*
!  
! For Kleinman-Bylander projectors calculation
! We use the same valence charge as in ps generation
!
          call vhrtre(rho,ve,rofi,drdi,s,nrval,a)  

! For PAO basis functions calculations 
! We use the "scaled" charge density of an ion of total charge "charge"
! As seen above, this ion could be the one involved in ps generation,
! or another specified at the time of basis generation.
! Example: Ba: ps ionic charge: +2
!              basis gen specified charge: +0.7

          do ir=2,nrval 
              rho_PAO(ir)=(zval-charge)*rho(ir)/chgvps
          enddo 
          call vhrtre(rho_PAO,vePAO,rofi,drdi,s,nrval,a) 
          chargesave(is)=charge
! 
! Check the exchange correlation functionals
! 
          call xc_check(icorr)
! 
          if (irel.eq.'rel') irelt=1
          if (irel.ne.'rel') irelt=0 
!
!        NOTE that atomxc expects true rho(r), not 4pir^2*rho(r)
!        We use auxrho for that
!
          do ir=2,nrval
            r2=rofi(ir)**2
            r2=4.0d0*pi*r2
            dc=rho(ir)/r2
            if (nicore.ne.'nc  ')  dc=dc+chcore(ir)/r2
              auxrho(ir)=dc
          enddo

          r2=rofi(2)/(rofi(3)-rofi(2))
          auxrho(1)=auxrho(2) -(auxrho(3)-auxrho(2))*r2

          call atomxc(irelt,nrval,nrmax,rofi,
     .                1,auxrho,ex,ec,dx,dc,vxc)

          ve(1:nrval)=ve(1:nrval)+vxc(1:nrval)

          do ir=2,nrval
            r2=rofi(ir)**2
            r2=4.0d0*pi*r2
            dc= rho_PAO(ir)/r2
            if (nicore.ne.'nc  ')  dc=dc+chcore(ir)/r2
              auxrho(ir)=dc
          enddo

          r2=rofi(2)/(rofi(3)-rofi(2))
          auxrho(1)=auxrho(2) -(auxrho(3)-auxrho(2))*r2

          call atomxc(irelt,nrval,nrmax,rofi,
     .                1,auxrho,ex,ec,dx,dc,vxc)

          vePAO(1:nrval)=vePAO(1:nrval)+vxc(1:nrval)

! 
! Now, we are going  to calculate the local pseudopotential and the
! KB projectors, this is only necessary if we are dealing with real 
! atoms (not with floating orbitals), i.e. only if flting is greater
! than zero
! 
          if (flting.gt.0.0d0) then
 
! Rgauss is approximately the maximum cut-off radius used in the 
! pseudopotential generation. 
! Rgauss is determined by comparison of the pseudopot.   
! Corresponding to different l ( this is not possible if we have     
! just one pseudopotential)                                          
! Rgauss2 is the radius where the pseudopotentials reach  the        
! asymptotic behaviour 2*Zval/r.                                     
! For just one pseudopotential Rgauss is taken equal to Rgauss2       
! 
         call radii_ps(vps,rofi,Zval,nrval,lmxkb,
     .          nrgauss, rgauss, rgauss2)
! 
!  
! Calculate local pseudopotential
! 
            if (rgauss2.gt.1.3d0*rgauss) then

               write(6,'(a)') "Using large-core scheme for Vlocal"

! In this case the atom core is so big that we do not have an asymptotic
! of 2*Zval/r until Rgauss2 (> Rc) . To retain the same asymptotic
! behaviour as in the pseudopotentials we modify the definition
! of the local potential, making it join the Vps's smoothly at rgauss.
! 
        write(6,'(/,a,f10.5)') 'atom: Estimated core radius ',
     .           rgauss2
        if (nicore.eq.'nc ')
     .  write(6,'(/,2a)') 'atom: Including non-local core corrections',
     .  ' could be a good idea'
!
! As all potentials are equal beyond rgauss, we can just use the
! s-potential here.
!
         call vlocal2(Zval, nrval, a, rofi, drdi, s, vps(:,0),
     .              nrgauss,vlocal,nchloc,chlocal)
! 
              else
! 
! In this case the pseudopotential reach to an asymptotic behaviour 2*Zval/r  
! for a radius approximately equal to Rc. We build a generalized-gaussian
! "local charge density" and set Vlocal as the potential generated by
! it. Note that chlocal is negative.
! 
            call vlocal1(Zval, nrval, a, rofi, drdi, s, rgauss,
     .                     vlocal,nchloc,chlocal)

            endif 
! 
! Save local-pseudopotential charge
!  
          Rchloc=rofi(nchloc)
          write(6,'(2a,f10.5)') 'atom: Maximum radius for' ,
     .        ' 4*pi*r*r*local-pseudopot. charge ',Rchloc

            call  comlocal(is,a,b,rofi,chlocal,nchloc,flting)
!
!------------------------------------------------------------------
!new: Save vlocal (actually r*vlocal +2*Zval)
!
         red_vlocal(1:nrval) =
     $        rofi(1:nrval)*vlocal(1:nrval) + 2.0_dp*Zval
!
!  Use the same eps=1.0e-4 as in other Vlocal-related routines
!
         do ir=nrval,2,-1
            if (abs(red_vlocal(ir)) .gt. 1.0e-4_dp) then
               nvlocal = ir + 1
               exit
            endif
         enddo
         call  com_vlocal(is,a,b,rofi,red_vlocal,nvlocal,
     $                    zval,flting)
!
          write(6,'(2a,f10.5)') 'atom: Maximum radius for' ,
     .        ' r*vlocal+2*Zval: ', rofi(nvlocal)           
!
!-------------------------------------------------------------------
! 
! ARRAY S FOR THE SCHRODINGER EQ. INTEGRATION
! 
         s(2:nrval)=drdi(2:nrval)*drdi(2:nrval)
         s(1)=s(2)

! CALCULATION OF THE KLEINMAN-BYLANDER PROJECTOR FUNCTIONS
! 
            call KBgen(is, a,b,rofi,drdi,s,
     .         vps, vlocal, ve, nrval, Zval, lmxkb, 
     .         nkbl, erefkb, nkb) 
            nkbmax(is)=nkb


            elseif(flting.lt.0.0d0) then 

! No Kleinman-Bylander projectors if floating orbitals
!          
            nkb=0 
            nkbmax(is)=0


! Zero local pseudopotential if floating orbitals
! 
            nchloc=0
            call  comlocal(is,a,b,rofi,chlocal,nchloc,flting)
            nvlocal=0
            call  comlocal(is,a,b,rofi,red_vlocal,nvlocal,flting)
            call  com_vlocal(is,a,b,rofi,red_vlocal,nvlocal,
     $           zval,flting)
! 

! ARRAY S FOR THE SCHRODINGER EQ. INTEGRATION
!          (Needed for the calculation of the basis set )

            s(2:nrval)=drdi(2:nrval)*drdi(2:nrval)
            s(1)=s(2)

            endif  


! CONSTRUCTION OF THE BASIS ORBITALS

! Method for the augmentation of the basis set
!
         if (basistype.ne.'user') then
             write(6,'(a,73("-"))') 'atom: '
             write(6,'(/,a)') 'atom: SANKEY-TYPE ORBITALS:'
           nzcontr=0
           do l=0,lmxo
             do nsm=1,nsemic(l)+1
              if(nzeta(l,nsm).gt.1) nzcontr=1
             enddo
           enddo
           if (nzcontr.eq.1) then
               write(6,'(2a)')
     .        'atom: Selected multiple-zeta basis: ',basistype
           endif
         endif


! Generate the PAOs basis set **
! 
       if (charge.le.0.0d0) then

         if (charge.lt.0.0d0) then 
          write(6,'(/,a)') 
     .    'atom: basis set generated (by rescaling the valence charge)'
          write(6,'(a,f8.4)')
     .    'atom: for an anion of charge ',charge 
         endif

         call Basis_gen(Zval,is, iz, a,b,rofi,drdi,s,
     .                   vps, ve, vePAO, nrval, lmxo,nsemic,
     .                   nzeta, rco, lambda, npolorb,
     .                   basistype, rphi, no, rinn, vcte, split_norm,
     $                   atm_label)
 
      else
          if (abs(charge-zval+chgvps).gt.1.0d-3) then 
            write(6,'(/,a)')
     .  'atom: basis set generated (by rescaling the valence charge)'
            write(6,'(a,f8.4)')
     .  'atom: for a cation of charge ',charge 
          else
            write(6,'(/,a)')
     .  'atom: basis set generated from the ionic configuration used'
            write(6,'(a)')
     .  'atom: to generate the pseudopotential'
          endif

        call Basis_gen(Zval,is, iz, a,b,rofi,drdi,s,
     .                   vps, vePAO, vePAO, nrval, lmxo,nsemic,
     .                   nzeta, rco, lambda, npolorb,
     .                   basistype, rphi, no, rinn, vcte,
     $                   split_norm,atm_label)
      endif
        write(6,'(a,i3)')
     .      'atom: Total number of Sankey-type orbitals:', no
      nomax(is)=no
! 
      if(flting.gt.0.0d0) then 
 
! Calculate initial populations for the atomic orbitals*
! 
         call atm_pop(is,iz,qtb(1:,is),qPAO,lmxo,
     .        nzeta,semic,nsemic,npolorb,basp) 

! Screening of the atomic local pseudopotential
! 
         call  Vna(is,Zval,qPAO,rphi,rco,nsemic,vlocal,
     .        a,b,rofi,drdi,nrval,lmxo,nVna)
        
!  CALCULATION OF THE ELECTROSTATIC SELF-ENERGY Of THE *
!  *** LOCAL PSEUDOPOTENTIAL CHARGE DENSITY ** 
         call slfe_local(slfe(is),vlocal,rofi,a,nVna,drdi) 
 
      elseif(flting.lt.0.0d0) then 

! Populations are zero because we have floating orbitals,*
! not real atoms**
         qtb(1:maxos,is)=0.0d0
! Zero neutral-atom pseudopotential if floating orbitals**
         nVna=0
         call comVna(is,a,b,rofi,vlocal,nVna,flting)
! Zero self-energy for the local pseudopotential charge****
         slfe(is)=0.0d0
      endif 

      elseif(iz.eq.-100) then

! SET UP THE MESH POINTS AND ITS DERIVATIVE***
! ONLY FOR FLOATING ORBS
           call set_mesh(a,b,rofi,drdi,s)
           flting=dsign(1.0d0,dble(iz)) 
           zvaltb(is)=0.0_dp
           chargesave(is) = 0.d0
! No core charge
           call  comcore(is,a,b,rofi,chcore,
     .          nrval,nicore,flting) 
! No KB projectors
           nkb=0 
           nkbmax(is)=0
! No local potential
           nchloc=0
           call  comlocal(is,a,b,rofi,chlocal,nchloc,flting)    
           call  com_vlocal(is,a,b,rofi,red_vlocal,nvlocal,
     $          zval,flting)
! Zero neutral-atom pseudopotential if floating orbitals
           nVna=0
           call comVna(is,a,b,rofi,vlocal,nVna,flting)
! Zero self-energy for the local pseudopotential charge
           slfe(is)=0.0d0
! Calculate Bessel functions
           call Bessel (is,a,b,rofi,drdi,s,
     .          lmxo,nzeta,rco,lambda, no)
 
           write(6,'(/a,i3)')
     .          'atom: Total number of floating Bessel orbitals:', no
           nomax(is)=no
        endif 

        write(6,'(/,a,73("_"))') 'atom: '


        end subroutine atom_main


        subroutine rc_vs_e(a,b,r,vps,
     .      ve,nrval,l,el,nnodo,rnodo)
C**
C   Calculate the position, rnodo, of the node number nnode of the 
C   radial wavefunction of the pseudopotential Vps, with angular 
C   momentum  l, and energy el.
C   D. Sanchez-Portal, July 1997.
C   Modify by DSP, July 1999
C**

        integer nrval, ir
        real(dp) r(nrval),
     .   el,vps(nrval),g(nrmax),drdi(nrmax),h(nrmax),ve(nrval)
        
        integer  l, nnodo, nn
        real(dp) rnodo, a, b
       
        real(dp) dexpa, ab, hi, gold, r0, g0, r1, g1
c       if (nrval.gt.nrmax) then   
c        write(6,*) 'Rc_vs_E : Nrmx must be increased to at least',
c    .         nrval
c       endif

        dexpa=exp(a)
        ab=a*b
        do ir=1,nrval
           
           drdi(ir)=ab
c          r2=ab/a-b

c          if(abs(r2-r(ir)).gt.1.0d-6) then 
c                write(6,*) ir, abs(r2-r(ir)),r2,r(ir)
c          endif
           ab=dexpa*ab
        enddo            

 
          do ir=2,nrval
            hi=vps(ir)+ve(ir)+dble(l*(l+1))/r(ir)**2-el
            hi=hi*(drdi(ir)**2)
            hi=hi+0.25d0*a**2
            h(ir)=hi
          enddo 
          h(1)=h(2)

          
          g(1)=0.0d0
          g(2)=1.0d0
          gold=1.0d0
          rnodo=r(nrval)
          nn=0
          do ir=3,nrval

            hi=(10.0d0*h(ir-1)*g(ir-1)+h(ir-2)*g(ir-2))/12.0d0
 
            hi=hi+2.0d0*g(ir-1)-g(ir-2)

            g(ir)=hi/(1.0d0-h(ir)/12.0d0)


            if((g(ir).eq.0.0d0).or.(g(ir)*gold.lt.0.0d0)) then
                  r0=r(ir-1)
                  g0=gold
                  r1=r(ir)
                  g1=g(ir)
                  rnodo=r0-g0*(r1-r0)/(g1-g0)
                  nn=nn+1
                  if(nn.eq.nnodo) goto 50
            endif
            
            gold=g(ir)

          enddo 

50        continue 

          end subroutine rc_vs_e

          
        subroutine polarization(a,r,psi,vps,
     .      ve,drdi,nrc,l,el,psipol,nrval)

C*
C      This routine calculate the polarization (unoccupied) 
C      orbitals with angular momentum l from the atomic 
C      orbitals with l-1 angular momentum, using a perturbative
C      approach.
C      The routine solve and inhomogeneus version of 
C      Schrodinger eqn.  
C      It is not an optimized algorithm!!!!!!!!!!!!!!!!!
C
C
C       cons1 is a big number. If it is choosen too 
C       big, too many iterations will be needed.
C       If cons1 is too small it may happen that 
C       the routine never converges to the solution.
C
C       If Rc is too big the very simple (unoptimized)
C       algorithm used here cannot converge. That is 
C       why Rc's are required to be smaller than Rint,
C       where Rint should not be greater than aprox. 15 Bohr
C       Written by Daniel Sanchez-Portal, July 1997
C 
     
        integer nrmin, niter
        real(dp) cons1, rint
        parameter(nrmin=1,niter=1000,
     .                    cons1=1.0d5,rint=15.0d0)

        integer nrval, nrc
        real(dp) r(nrval),psi(nrval),psipol(nrval),
     .   el,vps(nrval),g(nrmax),drdi(nrmax),h(nrmax),ve(nrval)
        
        real(dp) a, rmax, reduc, dl, hi, rnd1, c1, c2, rnodo, cons, gold
        real(dp) gmax, r0, g0, r1, g1, grmx, dff1, dff2, savecons, dnrm
        integer index, nnodes, iter, nnd, ir
        integer  l

c       if ((nrval.gt.nrmax).or.(nrc.gt.nrmax)) then   
c        write(6,*) 'POLARIZATION: Nrmx must be increased to at least',
c    .        max(nrval,nrc)
c       endif

        rmax=r(nrval)

        if(rmax.gt.rint) then 
          write(6,*) 'POLARIZATION: Rc for the polarization orbitals'
          write(6,*) 'must be smaller than ',rint,' Bohr'
         call die()
        endif

        do ir=nrc+1,nrval 
          psi(ir)=0.0d0 
        enddo 
 
          reduc=-0.5d0
C**** We calculate the polarization function with angular** 
C momentum l=l+dl
          dl=1
C****
          do ir=2,nrval
            hi=vps(ir)+ve(ir)+(l+dl)*(l+dl+1)/r(ir)**2-el
            hi=hi*(drdi(ir)**2)
            hi=hi+0.25d0*a**2
            h(ir)=hi
          enddo 
          h(1)=h(2)

          rnd1=0.0d0
          index=1
          nnodes=1
CInitialized c1 and c2 to arbitrary values.............
          c1=0.0d0
          c2=0.0d0
          do iter=1,niter
           rnodo=0.0d0
           if(index.eq.1) then 
              cons=cons1
              index=2
           else
              cons=c2
           endif
          
          g(1)=0.0d0
          do ir=1,nrmin+1
            g(ir)=cons*(r(ir)**(l+dl+1))/sqrt(drdi(ir))
          enddo 
          gold=g(nrmin+1)

          nnd=0
          gmax=0.0d0
          do ir=nrmin+2,nrval
            hi=-((r(ir)*psi(ir)+10.0d0*r(ir-1)*psi(ir-1)
     .         +r(ir-2)*psi(ir-2))/12.0d0)

            hi=hi+(10.0d0*h(ir-1)*g(ir-1)+h(ir-2)*g(ir-2))/12.0d0
 
            hi=hi+2.0d0*g(ir-1)-g(ir-2)

            g(ir)=hi/(1.0d0-h(ir)/12.0d0)
            gmax=max(gmax,abs(g(ir)))
            if((g(ir).eq.0.0d0).or.(g(ir)*gold.lt.0.0d0)) then
              nnd=nnd+1
              if (nnd.eq.nnodes) then 
                  r0=r(ir-1)
                  g0=gold
                  r1=r(ir)
                  g1=g(ir)
                  rnodo=r0-g0*(r1-r0)/(g1-g0)
              endif 
            endif
           gold=g(ir)
          enddo 

          grmx=g(nrval)/gmax

          if(((abs(rnodo-rmax).lt.1.0d-3).and.
     .      (abs(grmx).lt.1.0d-7) )
     .        .or. 
     .        ((rnodo.eq.0.0d0).and.
     .        (abs(grmx).lt.1.0d-7) ) ) goto 100

*We begin by finding a node!!!!**

          if((rnd1.eq.0.0d0).and.(rnodo.eq.0.0d0)) then  
             c2=(-reduc)*cons
             if(abs(c2).le.1.0d0/abs(cons1)) then 
               index=1
               rnd1=0.0d0
               reduc=(1.0d0+reduc)/2.0d0
             endif  
          elseif((rnd1.eq.0.0d0).and.(rnodo.ne.0.0d0)) then
              rnd1=rnodo
              c1=cons
              c2=2.0d0*cons
          endif  
       
****

      
*Now we lead this node to Rc*
          if((rnd1.ne.0.0d0).and.(rnodo.eq.0.0d0)) then 
              c2=0.50d0*(c1+c2)
          elseif((rnd1.ne.0.0d0).and.(rnodo.ne.0.0d0)) then 
              if(abs(rnd1-rnodo).gt.1.0d-6)then 
                 dff1=abs(rnd1-rmax) 
                 dff2=abs(rnodo-rmax) 
                 if(dff1.gt.dff2) then 
                   savecons=c2
                   c2=(rmax-rnd1)*(c1-c2)/(rnd1-rnodo)+c1
                   c1=savecons
                   rnd1=rnodo
                 else
                   c2=1.10d0*c2
                 endif 
              else

               if(abs(cons).gt.1.0d15) then 
                  nnodes=nnodes+1
                  index=1
                  rnd1=0.0d0
               else
                 c2=1.1d0*c2
               endif  

              endif 
           endif 

          enddo 
            write(6,*)'POLARIZATION: Iteration to find the polarization'
            write(6,*)'orbital has failed !!!!!!!!!'
            write(6,*)'Please try with a Rc no bigger than ',rnd1,
     .        ' Bohr'
          call die()
                          
100       continue
          dnrm=0.0d0
          do ir=1,nrval
             g(ir)=g(ir)*sqrt(drdi(ir))
             dnrm=dnrm+drdi(ir)*(g(ir)**2)
          enddo 
          dnrm=sqrt(dnrm)
          do ir=1,nrval
               psipol(ir) = g(ir)/dnrm
          enddo 

          end subroutine polarization
!
!
         subroutine parabola(a,b,nrc,rphi,rnrm,l,
     .              splnorm,cons1,cons2,nm) 

C
C For a value of the SplitNorm parameter equal
C to splnorm, this routine returns 
C the parameters Cons1, Cons2 and Nm, such that 
C the doble-Z function would be defined as 
C
C     phi(r)=r^l (cons1*r^2+cons2) if r < Rm
C
C     phi(r)= rphi(ir)/r           if r > Rm
C    
C with  Rm= b*[exp(a(nm-1)) -1 ] 
C Continuity in the wavefunction and its derivative
C is imposed.   
C The arrays rphi and rnrm belong to the input
C rphi(nrmax): The original PAO function multiply
C   by the radius.
C rnrm(nrmax): PAO's norm as a function of the 
C   radius. 
C
C  Written by D. Sanchez-Portal, July 1997.
C  
C Algorithm based on routine for Golden Section Search
C from Numerical Recipes.
C

         real(dp), intent(in)    ::  a, b
         integer, intent(in)     ::  nrc, l
         real(dp), intent(in)    ::  rphi(nrc), rnrm(nrc)
         real(dp), intent(in)    ::  splnorm
         real(dp), intent(out)   ::  cons1, cons2
         integer, intent(out)    ::  nm
          
         real(dp),  parameter  :: Ratio=0.61803399D0       

         real(dp) slopold, slop, rmin, gmin, cmin, rnrmin
         real(dp) gmax, cmax, rmax, rnrmax, valmin, valmax, gmed
         real(dp) cmed, rmed, rnrmed, valmed, g1, c1, r, rn1, val1
         real(dp) g2, c2, rn2, val2
         integer n0, n1, n2, n3
         integer ir, nr_max, nmin, nmax, nmed, iter, nlast

!         do ir= 1, nrc
!            write(77,*) ir, b*(exp(a*(ir-1))-1.0d0), rphi(ir), rnrm(ir)
!         enddo
!         call pxfflush(77)
         
C         Find the last maximum of the wavefunction

          nlast=nrc-2
          slopold=0.0d0
          do ir=nlast,2,-1
             slop=rphi(ir)-rphi(ir-1)
             if(slop*slopold.lt.0.0d0) exit
             slopold=slop
          enddo

          nr_max=ir-1
          rmin=b*(exp(a*(nr_max-1))-1.0d0)
          rmin=1.01d0*rmin
          nmin=nint(dlog(rmin/b+1.0d0)/a)+1
          nmin=max(nmin,2)
          nmax=nrc-1

!
!         Initial brackets
!
          call findp(nrc,nmin,rphi,a,b,l,cmin,gmin)
          rmin=b*(exp(a*(nmin-1))-1.0d0)
          call nrmpal(cmin,gmin,rmin,l,rnrmin)
          rnrmin=1.0d0+rnrmin-rnrm(nmin)

          call findp(nrc,nmax,rphi,a,b,l,cmax,gmax)
          rmax=b*(exp(a*(nmax-1))-1.0d0)
          call nrmpal(cmax,gmax,rmax,l,rnrmax)
          rnrmax=1.0d0+rnrmax-rnrm(nmax)

!
!         Start the algorithm to find the matching point at
!         which the *total* norm of the parabola+tail = splitnorm
!         (compare with the JPC paper: there it appears that only
!         the tail should have a norm=splitnorm.
!

C Under certain circunstances the algorithm is not going to work
          if(rnrmin.gt.splnorm.and.rnrmax.gt.splnorm) then

!!             if(rnrmin.gt.rnrmax) then
!!               nm=nmax
!!               cons1=cmax
!!               cons2=gmax
!!               splnorm=rnrmax
!!             else
!!               nm=nmin
!!               cons1=cmin
!!               cons2=gmin
!!               splnorm=rnrmin
!!             endif
             write(6,'(/,A,/,A)')
     .    'parabola: The program failed in finding a SPLIT orbital ',
     .    'parabola: with the desired splitnorm'
             call die()
          endif


          valmin=(splnorm-rnrmin)**2
          valmax=(splnorm-rnrmax)**2
          nmed=(nmin+nmax)/2
          do iter=1,nrc
            call findp(nrc,nmed,rphi,a,b,l,cmed,gmed)
            rmed=b*(exp(a*(nmed-1))-1.0d0)
            call nrmpal(cmed,gmed,rmed,l,rnrmed)
            rnrmed=1.0d0+rnrmed-rnrm(nmed)

            valmed=(splnorm-rnrmed)**2

            if((valmed.lt.valmin).and.(valmed.lt.valmax)) goto 20
            nmed=nmed+1
            if(nmed.eq.nmax) goto 15
          enddo
 15               continue
          nmed=(nmin+nmax)/2
          do iter=1,nrc
             call findp(nrc,nmed,rphi,a,b,l,cmed,gmed)
             rmed=b*(exp(a*(nmed-1))-1.0d0)
             call nrmpal(cmed,gmed,rmed,l,rnrmed)
             rnrmed=1.0d0+rnrmed-rnrm(nmed)

             valmed=(splnorm-rnrmed)**2


             if((valmed.lt.valmin).and.(valmed.lt.valmax)) goto 20
             nmed=nmed-1
             if(nmed.eq.nmin) goto  20
          enddo
 20               continue

          if(nmed.eq.nmin) then
             if(valmin.lt.valmax) then
                nm=nmin
                cons1=cmin
                cons2=gmin
             elseif(valmax.le.valmin) then
                nm=nmax
                cons1=cmax
                cons2=gmax
             endif
             return
           endif

C    Ahora ya tenemos el minimo en un intervalo


            n0=nmin
            n3=nmax
            if(abs(nmed-nmax).gt.abs(nmed-nmin)) then
               n1=nmed
               n2=nmed+nint((1.0d0-ratio)*(nmax-nmed))
            else
               n2=nmed
               n1=nmed-nint((1.0d0-ratio)*(nmed-nmin))
            endif
            call findp(nrc,n1,rphi,a,b,l,c1,g1)
            r=b*(exp(a*(n1-1))-1.0d0)
            call nrmpal(c1,g1,r,l,rn1)
            rn1=1.0d0+rn1-rnrm(n1)
            val1=(splnorm-rn1)**2

            call findp(nrc,n2,rphi,a,b,l,c2,g2)
            r=b*(exp(a*(n2-1))-1.0d0)
            call nrmpal(c2,g2,r,l,rn2)
            rn2=1.0d0+rn2-rnrm(n2)
            val2=(splnorm-rn2)**2

1           if(abs(n3-n0).gt.1) then
              if(val2.lt.val1) then
               n0=n1
               n1=n2
               n2=nint(ratio*n1+(1-ratio)*n3)
c              val0=val1
               val1=val2
               call findp(nrc,n2,rphi,a,b,l,c2,g2)
               r=b*(exp(a*(n2-1))-1.0d0)
               call nrmpal(c2,g2,r,l,rn2)
               rn2=1.0d0+rn2-rnrm(n2)
               val2=(splnorm-rn2)**2
              else
               n3=n2
               n2=n1
               n1=nint(ratio*n2+(1-ratio)*n0)
c              val3=val2
               val2=val1
               call findp(nrc,n1,rphi,a,b,l,c1,g1)
               r=b*(exp(a*(n1-1))-1.0d0)
               call nrmpal(c1,g1,r,l,rn1)
               rn1=1.0d0+rn1-rnrm(n1)
               val1=(splnorm-rn1)**2
              endif
             goto1
             endif
             if(val1.lt.val2) then
                  cons2=g1
                  cons1=c1
                  nm=n1
             else
                 cons2=g2
                 cons1=c2
                 nm=n2
             endif

           end subroutine parabola


          subroutine findp(nrc,nm,rphi,a,b,l,cons1,cons2)
          integer, intent(in)   ::  nrc, nm, l
          real(dp), intent(in)  ::  a, b
          real(dp), intent(out) ::  cons1, cons2
          real(dp), intent(in)  ::  rphi(nrc)

C  This routine provides the constants Cons1 and 
C  Cons2 and described in subroutine 'parabola' 
C  for fitting at point of index nm

          logical  :: first_time = .true.   ! will be saved
          logical, save  :: keep_findp_bug  

          real(dp) rm, rm1, rm2, drdi_local, frsp
          real(dp) dfrsp

          if (first_time) then
             keep_findp_bug = fdf_boolean("PAO.Keep.Findp.Bug", .false.)
             first_time = .false.
          endif

          if (keep_findp_bug) then
             ! old, wrong code
             rm=b*(exp(a*(nm-1)) + 1.0d0) 
             rm1=b*(exp(a*(nm-2)) + 1.0d0)
             rm2=b*(exp(a*nm) + 1.0d0)
          else
             ! correct code
             rm=b*(exp(a*(nm-1)) - 1.0d0)     
             rm1=b*(exp(a*(nm-2)) -  1.0d0)   
             rm2=b*(exp(a*nm) - 1.0d0)        
          endif

          drdi_local =a*b*exp(a*(nm-1))

          frsp=rphi(nm)/rm
          dfrsp=0.5d0*(rphi(nm+1)/rm2
     .       -rphi(nm-1)/rm1)
          dfrsp=dfrsp/drdi_local

          cons1= 0.5d0*(dfrsp*rm-l*frsp)/(rm**(l+2))
          cons2= frsp/(rm**l)-cons1*(rm**2)

          end subroutine findp
!
!
          subroutine nrmpal(c1,c2,r,l,dnrm)
          real(dp), intent(in)  :: c1, c2, r
          real(dp), intent(out) :: dnrm
          integer, intent(in)   :: l

C returns the norm of a parabolic function
C    f(r')= r'^l (c1*r'^2 + c2)  r'< r
C           0 otherwise
C
!         The returned value is \int{ (rf)**2}
          
          dnrm=(c1**2)*r**(2*l+7)/(2*l+7) 
     .     + (2.0d0*c1*c2)*r**(2*l+5)/(2*l+5) 
     .     + (c2**2)*r**(2*l+3)/(2*l+3)

          end subroutine nrmpal

!---------------------

      subroutine radii_ps(vps,rofi,Zval,nrval,lmxkb,
     .     nrgauss, rgauss, rgauss2)

C     This routine returns the maximum radius for the
C     Kleinman-Bylander projectors with a standard choice
C     of the local potential.
C     Check also at which radius the asymptotic 2*Zval/r
C     behaviour is achieved.
C     D. Sanchez-Portal, Aug. 1998


      real(dp), intent(in)    ::  vps(:,0:), rofi(:)
      real(dp), intent(in)    ::  Zval
      integer, intent(in)   ::  nrval, lmxkb
      integer, intent(out)  ::  nrgauss
      real(dp), intent(out)   ::  rgauss, rgauss2

      real(dp) dincv, r
      integer ir, l, nrgauss2

      real(dp), parameter     ::  eps=1.0d-4

C**   Iterate over the possible local potentials**

      rgauss=0.0_dp
      rgauss2=0.0_dp
      nrgauss=0
      nrgauss2=0

      do l=0,lmxkb-1
         do ir=nrval,2,-1
            dincv=abs(vps(ir,l)-vps(ir,lmxkb))
            if(dincv.gt.eps) exit
         enddo
         rgauss=max(rofi(ir),rgauss)
         nrgauss=max(ir,nrgauss)
      enddo
!
!     New: Use all potentials, not just l=0, since
!     potentials with larger l can converge later...
!
      do l=0,lmxkb
         do ir=nrval,2,-1
            r=rofi(ir)
            dincv=abs(vps(ir,l)*r+2.0_dp*zval)
            if(dincv.gt.eps) exit
         enddo
         write(6,'(a,i1,a,f8.4)')
     $     'V l=', l,' = -2*Zval/r beyond r=', rofi(ir)
         rgauss2=max(rofi(ir),rgauss2)
         nrgauss2=max(ir,nrgauss2)
      enddo

      if(lmxkb.eq.0) then
         rgauss=rgauss2
         nrgauss=nrgauss2
      endif

      write(6,'(a,f8.4)') 'All V_l potentials equal beyond r=', rgauss
      write(6,'(a)')
     $     "This should be close to max(r_c) in ps generation"
      write(6,'(a,f8.4)')
     $     'All pots = -2*Zval/r beyond r=', rgauss2

      end subroutine radii_ps

            
!-----------------------------------------------------

      subroutine vlocal1(Zval, nrval, a, rofi, drdi, s, rgauss,
     .     vlocal, nchloc, chlocal)

C     This routine generates a smooth local pseudopotential.
C     Written by D. Sanchez-Portal, Aug. 1998

      real(dp), intent(in) :: Zval, a
      integer, intent(in) :: nrval
      real(dp), intent(in) :: rofi(:), drdi(:), s(:)
      real(dp), intent(out) :: vlocal(:)
      real(dp), intent(out) :: chlocal(:)
      integer, intent(out) :: nchloc
      real(dp), intent(inout) :: rgauss      !!???


C     *Internal variables* 

      real(dp) van, factor, alp, cutoff1, cutoff2,
     .         qtot, eps, pi, chc, r, Rchloc, rhor1, rhor
      integer ir
      character loctype*3

      parameter(eps=1.0d-4)

C**   Usual local potential 
C     (generated with an optimum Vandebilt function)**
      loctype='new' 

C***  The very first local potential used by SIESTA was 
C     the electrostatic potential generated by a gaussian 
C     distribution ===> loctype='old' 
C     loctype='old'
C***  
      pi=acos(-1.0_dp)

C     Local-potential size parameter 'rgauss'
C     We choose as a smooth pseudopotential the one generated 
C     by a 'Vanderbilt-function' charge distribution. We have to select 
C     the size of this distribution somehow.
C     'Vanderbilt-functions' are of the form :
C     p(r)=N*exp(-(sinh(van*r)/sinh(van))**2)
C     when van---> 0 we will obtain a 'gaussian'
C     when van---> Inf. we will obtain a step function
C     Some test has revealed that the best election to achieve 
C     a good convergence in real and reciprocal space is b in the 
C     range 0.5-1.0 .
C     *

C     So, the 'gaussian' charge distribution 
C     must go to zero at a distance 'rgauss'.

      if(loctype.eq.'new') then          

C     We take a 'Vanderbilt-function' as local potential
C     van=1.0_dp all the parameter have optimized for this value 

         van=1.0_dp
         cutoff1=3.63_dp
         cutoff2=5.48_dp
C**   99% of charge inside Rgauss**
c     factor=1.627_dp

C**   99.9% of charge inside Rgauss
         factor=1.815_dp
         
C     * Scaling factor for local-pseudopot. charge**
         alp=factor/rgauss

         write(6,'(/,a,f10.3,a)')
     .        'VLOCAL1: 99.0% of the norm of Vloc inside ',
     .        (alp*cutoff1)**2,' Ry'
         write(6,'(a,f10.3,a)')
     .        'VLOCAL1: 99.9% of the norm of Vloc inside ',
     .        (alp*cutoff2)**2,' Ry'

      else 

C     This is just a gaussian !!!!!!!!!!!!!!!!!                 
         van=0.00001_dp 
         rgauss=0.80_dp
         factor=2.0_dp
C     * Scaling factor for local-pseudopot. charge**
         alp=factor/rgauss  
      endif 
!--------------------

      qtot=0.0_dp 
      rhor1 = vander(van,alp*rofi(1))     ! This is 1...
      do ir=1,nrval
         r=rofi(ir) 
         rhor = vander(van,alp*r)
         chlocal(ir)=(-4.0_dp)*pi*rhor*r*r
         qtot=qtot+rhor*drdi(ir)*r*r
      enddo

      qtot=4.0_dp*pi*qtot 
      nchloc=0  
      do ir=nrval,1,-1
         chc=zval*chlocal(ir)/qtot
         chlocal(ir)=chc
         if((abs(chc).gt.eps).and.(nchloc.eq.0)) then    
            nchloc=ir+1
         endif
      enddo 
      Rchloc=rofi(nchloc)
!
!     Note that the above cutoff is for 4*pi*r*r*rho_local(r)...
!
      call vhrtre(chlocal,vlocal,rofi,drdi,s,nrval,a)

      do ir=2,nrval 
         r=rofi(ir)  
         chlocal(ir)=chlocal(ir)/(4.0_dp*pi*r*r)
!
!     Poor man's cutoff!! Largely irrelevant?
!      
         if (r.gt.1.1_dp*Rchloc) then
            vlocal(ir)=(-2.0_dp)*zval/rofi(ir)
         endif

      enddo 
      chlocal(1)= -rhor1* zval/qtot

      end subroutine vlocal1


      subroutine vlocal2(Zval, nrval, a, rofi, drdi, s, vps,
     .     nrgauss,vlocal,nchloc,chlocal) 

C     This routine generates the local pseudopotential appropiate 
C     for species with  a large core.
C     Written by D. Sanchez-Portal, Aug. 1998
      
      real(dp), intent(in) :: Zval, a
      integer, intent(in) :: nrval
      integer, intent(inout) :: nrgauss
      real(dp), intent(in) :: rofi(:), drdi(:), s(:), vps(:)
      real(dp), intent(out) :: vlocal(:), chlocal(:)
      integer, intent(out) :: nchloc


C     Internal variables****

      real(dp) 
     .     vlc, r, dev, dev2, dev3, var1, var2, var3, v1, v2, v3, v4,
     .     dm11, dm12, dm13, dm21, dm22, dm23, dm31, dm32, dm33, 
     .     g0, g1, g2, g3, g4, d2g, d2u, cons, a2b4, qtot, pi   
      integer 
     .     ndevfit, ir  

      real(dp), parameter  :: eps=1.0d-5
      
      pi=acos(-1.0_dp)        

C     Continuity up to second derivative***
      ndevfit=2
C     Continuity up to third derivative****
C     ndevfit=3

      nrgauss=nrgauss+3        !! For good measure...

      do ir=1,nrval
         vlocal(ir)=vps(ir)*rofi(ir)
      enddo 

      ir=nrgauss
      dev=(vlocal(ir+1)-vlocal(ir-1))*0.5_dp
      dev2=(vlocal(ir+1)+vlocal(ir-1)-2.0_dp*vlocal(ir))
      dev3=(vlocal(ir+2)-2.0_dp*vlocal(ir+1)
     .     +2.0_dp*vlocal(ir-1)-vlocal(ir-2))*0.5_dp
      dev3=(dev3-3.0_dp*a*dev2+2.0_dp*(a**2)*dev)
     .     /(drdi(ir)**3)
      dev2=(dev2-a*dev)/(drdi(ir)**2)
      dev=dev/drdi(ir)

C     Local potential is Vloc(r)=v3*exp(v1*r^2+v2*r^3) 
C     inside Rgauss and equals the 
C     all-electron atomic potential outside Rgauss
C     We impose the continuity up to second derivative
      
      if(ndevfit.eq.2) then               
         vlc=vlocal(nrgauss)
         r=rofi(nrgauss)

         var1=dev/vlc-1.0_dp/r
         var2=dev2/vlc-2.0_dp*var1/r -(var1**2)

         dm11=2.0_dp*r
         dm12=3.0_dp*r*r
         dm21=2.0_dp
         dm22=6.0_dp*r

         v1=(dm22*var1-dm12*var2)/(6.0_dp*r*r)
         v2=(dm11*var2-dm21*var1)/(6.0_dp*r*r)
         v3=vlc/(r*exp((v1+v2*r)*r*r))


c     elseif(ndevfit.eq.3) then 
      else

C     We can also construct a local potential 
C     Vloc(r)=v4*exp(v1*r^2+v2*r^3+v3*r^4),
C     this new coefficient allows us to impose the continuity 
C     of the potential up  to the third derivative.

         vlc=vlocal(nrgauss)
         r=rofi(nrgauss)
         
         var1=dev/vlc-1.0_dp/r
         var2=dev2/vlc-2.0_dp*var1/r-(var1**2)
         var3=dev3/vlc-3.0_dp*var1*var2-(var1**3)
     .        -3.0_dp*(var1**2+var2)/r

         dm11=2.0_dp*r
         dm12=3.0_dp*r*r
         dm13=4.0_dp*r*r*r
         dm21=2.0_dp
         dm22=6.0_dp*r
         dm23=12.0_dp*r*r
         dm31=0.0_dp
         dm32=6.0_dp
         dm33=24.0_dp*r

         v1=((var1*dm22*dm33+var2*dm13*dm32+var3*dm12*dm23)
     . -(var3*dm22*dm13+var1*dm32*dm23+var2*dm12*dm33))/(48.0_dp*r*r*r)
         v2=((var2*dm11*dm33+var3*dm21*dm13+var1*dm23*dm31)
     . -(var2*dm31*dm13+var3*dm23*dm11+var1*dm21*dm33))/(48.0_dp*r*r*r)
         v3=((var3*dm11*dm22+var2*dm12*dm31+var1*dm32*dm21)
     . -(var1*dm22*dm31+var3*dm21*dm12+var2*dm11*dm32))/(48.0_dp*r*r*r)
         v4=vlc/(r*exp((v1+v2*r+v3*r*r)*r*r))
         
      endif 
      
      do ir=1,nrval
         r=rofi(ir)
         if(ir.le.nrgauss) then 
            
C**   If second derivative fit***
            if(ndevfit.eq.2) then 
               vlocal(ir)=v3*exp((v1+v2*r)*r*r)

C**   If third derivative fit****
            elseif(ndevfit.eq.3) then 

               vlocal(ir)=v4*exp((v1+v2*r+v3*r*r)*r*r)
c**** 
            endif 

         else
            vlocal(ir)=vps(ir)
         endif 

      enddo 

C     Once we have the local potential we define the 'local-pseudopotential 
C     charge' which help us to calculate the electrostatic interation 
C     between the ions
!
!     Poisson's eq.:
!
!           1/r* d2(rV)/dr2 = -8*pi*rho
!
      a2b4=0.25_dp*a*a 
      qtot=0._dp 
      do ir=1,nrval-1
         
         g2=vlocal(ir)*rofi(ir)
!
!        To determine the chlocal cutoff, use the reduced_vlocal cutoff
!
         if(abs(g2+2.0_dp*zval).lt.eps) exit   !exit loop

         if(ir.gt.nrgauss) then  

            if((ir.gt.2).and.(ir.lt.(nrval-1))) then 
               g0=vlocal(ir-2)*rofi(ir-2)/s(ir-2)
               g1=vlocal(ir-1)*rofi(ir-1)/s(ir-1)
               g2=vlocal(ir)*rofi(ir)/s(ir)
               g3=vlocal(ir+1)*rofi(ir+1)/s(ir+1)
               g4=vlocal(ir+2)*rofi(ir+2)/s(ir+2)

               d2g=(16.0_dp*(g1+g3)-(g0+g4)-30.0_dp*g2)/12.0_dp
               
            else
               g1=vlocal(ir-1)*rofi(ir-1)/s(ir-1)
               g2=vlocal(ir)*rofi(ir)/s(ir) 
               g3=vlocal(ir+1)*rofi(ir+1)/s(ir+1)  
               
               d2g=g1+g3-2.0_dp*g2

            endif  

            d2u=d2g-a2b4*g2

            r=rofi(ir)
            cons=8.0_dp*pi*r*drdi(ir)*s(ir)
            chlocal(ir)=(-d2u)/cons
            qtot= qtot + 0.5_dp*d2u*r/s(ir)

         else

C     If second derivative fit
            if(ndevfit.eq.2)  then
               r=rofi(ir)

               g0=v3*exp((v1+v2*r)*r**2)
               g1=(2.0_dp*v1+3.0_dp*v2*r)
               g2=2.0_dp*v1+6.0_dp*v2*r
               g3=(g2+g1*g1*r*r+2.0_dp*g1)*g0
               
               cons=8.0_dp*pi
               chlocal(ir)= (-g3)/cons
               qtot= qtot  + 0.5_dp*g3*r*r*drdi(ir)
               
C**** If third derivative fit
               
            elseif(ndevfit.eq.3)  then

               r=rofi(ir)
               
               g0=v4*exp((v1+v2*r+v3*r*r)*r*r)     
               g1=(2.0_dp*v1+3.0_dp*v2*r+4.0_dp*v3*r*r)
               g2=(2.0_dp*v1+6.0_dp*v2*r+12.0_dp*v3*r*r)    
               g3=(g2+g1*g1*r*r+2.0_dp*g1)*g0   

               cons=8.0_dp*pi
               chlocal(ir)= -g3/cons
               qtot= qtot  + 0.5_dp*g3*r*r*drdi(ir)
            endif 

         endif
      enddo              
!
!     This sets the cutoff point for chlocal in a rather
!     arbitrary way, as that in which Vlocal "equals" 2Z/r
!
      nchloc=ir          

      do ir=1,nchloc-1
         chlocal(ir)=zval*chlocal(ir)/qtot
      enddo  
      do ir=nchloc,nrval
         chlocal(ir)=0.0_dp
      enddo 

      end subroutine vlocal2

!
!----------------------------------------------------------------
!
        subroutine schro_eq(Zval,rofi,vps,ve,s,drdi,
     .        nrc,l,a,b,nnodes,nprin,
     .        e,g) 
        
          implicit none   

          integer  
     .        nrc,l,nnodes,nprin
          real(dp)
     .        Zval, rofi(*),vps(*),ve(*),s(nrc),drdi(*),a,b,e,g(*)
     
C* Internal variables***  

           real(dp)
     .       a2b4, h(nrmax), r2, vtot, rmax, dr,
     .        y(nrmax), dnrm, phi, dsq
           integer
     .        ir

         a2b4=a*a*0.25d0

            do ir=2,nrc 
               g(ir)=0.0d0
               r2=(rofi(ir)**2)
               vtot=vps(ir)+ve(ir)+dble(l*(l+1))/r2
               h(ir)=vtot*s(ir)+a2b4
            enddo
            h(1)=h(2)
            g(1)=0.0d0 

            e=-((zval/dble(nprin))**2)
            dr=-1.0d6
            rmax=rofi(nrc)
            call egofv(h,s,nrc,e,g,y,l,zval,a,b,rmax,
     .        nprin,nnodes,dr)

            
            do ir=2,nrc
              phi=g(ir)
              dsq=sqrt(drdi(ir))
              phi=phi*dsq
              g(ir)=phi
            enddo 
            g(1)=0.0d0 

            dnrm=0.0d0
            do ir=2,nrc
               phi=g(ir)
               dnrm=dnrm+phi*phi*drdi(ir) 
            enddo
            dnrm=sqrt(dnrm)

            do ir=2,nrc
               g(ir)=g(ir)/dnrm
            enddo 
            
            end subroutine schro_eq
             

      subroutine ghost(Zval,rofi,vps,vlocal,
     .     ve,s,drdi,nrval,l,a,b,
     .     nrc, eigenl,rphi,ighost)

C     This routine checks the possible existence of ghost states. 
C     Input:
C     vps(*)    : pseudopotential for angular momentum l
C     vlocal(*) : local pseudopotential
C     rphi (*)  : first radial pseudowavefunctions for Vps.
C     eigenl    : eigenvalue 
C     
C     Output:
C     Sets ighost to 1 if it finds a ghost state.
!     ighost is a saved variable in the caller.
C     
C     Written by D. Sanchez-Portal, Aug. 1998
C     

      integer, intent(in)   :: nrval,l
      integer, intent(inout)   :: nrc
      integer, intent(inout)  :: ighost

      real(dp), intent(in) :: zval, a, b, eigenl
      real(dp), intent(in) :: rofi(:), vps(:), vlocal(:), ve(:)
      real(dp), intent(in) :: s(:), drdi(:), rphi(:)

C     * Internal variables***

      real(dp)   dnrm, avgv, phi,
     .     elocal1, elocal2, g(nrmax), vl, vphi, dkbcos

      integer  ir, nprin, nnode

!     Compares the reference energy with the eigenvalues of the
!     local potential, for ghost analysis.

! ATTENTION , 'Ve' is the screenig potential generated from valence
! pseudo-charge given by the pseudopotential generation program

!!!      write(6,*) 'Entering ghost: l, eigenl', l, eigenl

C     CALCULATE KB-COSINE

      nrc=min(nrc,nrval)
      dnrm=0.0d0
      avgv=0.0d0
      do 30 ir=2,nrc
         vl=(vps(ir)-vlocal(ir))
         phi=rphi(ir)
         vphi=vl*phi
         dnrm=dnrm+vphi*vphi*drdi(ir)
         avgv=avgv+vphi*phi*drdi(ir)
 30   continue
      dkbcos=dnrm/(avgv+1.0d-20)

C     GHOST ANALYSIS

!
!     Do not compute the second eigenvalue it it is not necessary
!

      if (dkbcos.lt.0d0) then

         nprin=l+1
         nnode = 1
         call schro_eq(Zval,rofi,vlocal,ve,s,drdi,
     .        nrval,l,a,b,nnode,nprin,
     .        elocal1,g)

         if(eigenl.gt.elocal1) then
            write(6,"(a,i3)")
     .           'GHOST: WARNING: Ghost state for L =', l
            ighost=1
         else
            write(6,'(a,i3)') 'GHOST: No ghost state for L =',l
         endif

      else if (dkbcos.gt.0.0d0) then

         nprin=l+1
         nnode = 2
         call schro_eq(Zval,rofi,vlocal,ve,s,drdi,
     .        nrval,l,a,b,nnode,nprin,
     .        elocal2,g)

         if(eigenl.gt.elocal2) then
            write(6,"(a,i3)")
     .           'GHOST: WARNING: Ghost state for L =', l
            ighost=1
         else
            write(6,'(a,i3)') 'GHOST: No ghost state for L =',l
         endif

      else if (dkbcos.eq.0.0d0) then

         write(6,"('GHOST: vps = vlocal, no ghost for L =',i3)") l

      endif

c     write(6,*) 'GHOST: Ground state vlocal for L=',l,elocal1
c     write(6,*) 'GHOST: First excited state for L=',l,elocal2

      end subroutine ghost


        subroutine KBproj(rofi,drdi,vps,vlocal,nrwf,l,rphi,
     .    dkbcos,ekb,proj,nrc)  
C****
C    This routine calculates the Kleinman-Bylander projector
C    with angular momentum l.
C
C  Written by D. Sanchez-Portal, Aug. 1998
C  Modified by DSP to allow more than one projector per l, July 1999.
C****

          implicit none

          real(dp)
     .        rofi(*),vps(*),drdi(*),proj(*),
     .        rphi(*), vlocal(*), dkbcos,ekb
          integer
     .        nrwf, l, nrc

C* Internal variables***

          real(dp) 
     .        eps, dnrm, vl, vphi, avgv, r, phi, dknrm,
     .        dincv, rc, rphi2(nrmax,nkbmx), vii(nkbmx),
     .        sum, vij(nkbmx)
          integer
     .        ir, l_last, nkb_last, jkb, ikb

          logical  called
             
          parameter (eps=1.0d-6)


          save called, l_last, nkb_last, rphi2, vii
          data called /.false./
   
              
             if(called) then 
               if(l.ne.l_last) then 
                  l_last=l 
                  nkb_last=0
               endif
             else
               called=.true.
               l_last=l
               nkb_last=0
             endif


C We need to orthogonalize to the previous projectors. 
C We follow the scheme proposed by Blochl, PRB 41, 5414 (1990)
             ikb=nkb_last+1
             if(nkb_last.eq.0) then 
               do ir=1,nrwf
                  rphi2(ir,1)=rphi(ir)
               enddo 
             else
               do jkb=1,nkb_last
                 sum=0.0d0
                 do ir=1,nrwf
                    vl=vps(ir)-vlocal(ir)
                    sum=sum+
     .              rphi(ir)*vl*rphi2(ir,jkb)*drdi(ir)
                 enddo
                 vij(jkb)=sum
               enddo
               do ir=1,nrwf
                  sum=0.0d0
                  do jkb=1,nkb_last
                     sum=sum+
     .               vij(jkb)*rphi2(ir,jkb)/(vii(jkb)+1.0d-20)
                  enddo
                  rphi2(ir,ikb)=rphi(ir)-sum
               enddo
             endif
             nkb_last=ikb
C Normalize the new function            
             dnrm=0.0d0
             do ir=1,nrwf
               dnrm=dnrm
     .          + drdi(ir)*(rphi2(ir,ikb)**2)
             enddo 
             dnrm=sqrt(dnrm)
             if(dnrm.lt.eps) then 
               do ir=1,nrwf
                 rphi2(ir,ikb)=0.0d0
               enddo 
             else
               do ir=1,nrwf
                 rphi2(ir,ikb)=rphi2(ir,ikb)/dnrm
               enddo 
             endif

             dnrm=0.0d0
             avgv=0.0d0
             do 10 ir=2,nrwf
               r=rofi(ir)
               vl=(vps(ir)-vlocal(ir))
               phi=rphi2(ir,ikb)
               vphi=vl*phi
               dnrm=dnrm+vphi*vphi*drdi(ir)
               avgv=avgv+vphi*phi*drdi(ir)
  10         continue
             vii(ikb)=avgv

             ekb=dnrm/(avgv+1.0d-20)
             dknrm=1.0d0/(sqrt(dnrm)+1.0d-20)
             dkbcos=avgv*dknrm
             
            

CDEFINE THE CUT-OFF RADII FOR THE KB PROJECTORS**
C Warning these radii should be quite short, if it is not the case
C something is probably wrong in this part of the program.
C It will display a warning if Rc>4.5 a.u.or Rc < 0.5a.u.!!!!!!!!!!!!

              do 20 ir=nrwf,2,-1 
                 phi=(rphi2(ir,ikb)/rofi(ir))*dknrm
                 dincv=abs((vps(ir)-vlocal(ir))*phi)
                 if (dincv.gt.eps) then 
                   if (ir.ge.nrwf-1) then
                     write(6,"(2a,/,2a)") 'KBproj: WARNING: ',
     .               'KB projector does not decay to zero',
     .               'KBproj: WARNING: ',       
     .        'parameter Rmax in routine KBgen should be increased'
                   endif
                   goto 21
                 endif 
20            continue

21            nrc=ir+1
              rc=rofi(nrc)
      
              if(rc.lt.0.5d0) then
                write(6,"('KBproj: WARNING: Rc(',i2,')=',f12.4)")l,rc
                write(6,"(2a)") 'KBproj: WARNING: ',
     .            'Cut of radius for the KB projector too small' 
              elseif(rc.gt.4.5d0) then
               write(6,"('KBproj: WARNING: Rc(',i2,')=',f12.4)")l,rc
              write(6,"(2a)") 'KBproj: WARNING: ',
     .            'Cut of radius for the KB projector too big'
               write(6,"(2a)") 'KBproj: WARNING: ',
     .            'Increasing the tolerance parameter eps'
               write(6,"(a)") 'KBproj: WARNING: might be a good idea'
              endif

              do 30 ir=2,nrc
                r=rofi(ir)
                vl=(vps(ir)-vlocal(ir))
                phi=rphi2(ir,ikb)/r
                vphi=vl*phi*dknrm
                proj(ir)=vphi/r**l
30           continue
             proj(1)= ( proj(2)*rofi(3)**2 - proj(3)*rofi(2)**2 ) /
     .             (      rofi(3)**2 -      rofi(2)**2 )

             end subroutine KBproj

        subroutine xc_check(icorr)

C Checking the functional used for exchange-correlation energy.
C Written by D. Sanchez-Portal, Aug. 1998

        use xcmod, only : nXCfunc, XCfunc, XCauth

C Passed variable
        character(len=2), intent(in) :: icorr

C Local variables
        integer                      :: nf
        character(len=40)            :: ps_string

        ps_string = "Unknown atomic XC code"
        if (icorr .eq. "ca") ps_string ="Ceperley-Alder"
        if (icorr .eq. "pw") ps_string ="Perdew-Wang 1992"
        if (icorr .eq. "pb")
     $         ps_string ="GGA Perdew, Burke & Ernzerhof 1996"
        if (icorr .eq. "rp") ps_string ="GGA RPBE"
        if (icorr .eq. "rv") ps_string ="GGA revPBE"
        if (icorr .eq. "wc") ps_string ="GGA Wu-Cohen"
        if (icorr .eq. "bl") ps_string ="GGA Becke-Lee-Yang-Parr"

C Loop over functionals
        do nf = 1,nXCfunc

          write(6,'(/a)') 'xc_check: Exchange-correlation functional:'
          if (((XCauth(nf).eq.'CA').or.(XCauth(nf).eq.'PZ')).and.
     .       ((XCfunc(nf).eq.'LDA').or.(XCfunc(nf).eq.'LSD'))) then

            write(6,'(a)') 'xc_check: Ceperley-Alder'
            if (icorr.ne.'ca'.and.nXCfunc.eq.1)
     $          write(6,'(a,1x,2a)')
     .          'xc_check: WARNING: Pseudopotential generated with',
     $           trim(ps_string), " functional"


          elseif((XCauth(nf).eq.'PW92').and.
     .      ((XCfunc(nf).eq.'LDA').or.(XCfunc(nf).eq.'LSD'))) then

            write(6,'(a)') 'xc_check: Perdew-Wang 1992'
            if (icorr.ne.'pw'.and.nXCfunc.eq.1) 
     $          write(6,'(a,1x,2a)')
     .          'xc_check: WARNING: Pseudopotential generated with',
     $           trim(ps_string), " functional"


          elseif((XCauth(nf).eq.'PBE').and.(XCfunc(nf).eq.'GGA')) then

            write(6,'(a)')
     .       'xc_check: GGA Perdew, Burke & Ernzerhof 1996'
            if (icorr.ne.'pb'.and.nXCfunc.eq.1) 
     $          write(6,'(a,1x,2a)')
     .          'xc_check: WARNING: Pseudopotential generated with',
     $           trim(ps_string), " functional"

          elseif ((XCauth(nf).eq.'LYP').and.(XCfunc(nf).eq.'GGA')) then

            write(6,'(a)')  'xc_check: GGA Becke Lee Yang Parr'
            if (icorr.ne.'bl'.and.nXCfunc.eq.1) 
     $          write(6,'(a,1x,2a)')
     .          'xc_check: WARNING: Pseudopotential generated with',
     $           trim(ps_string), " functional"

          elseif ((XCauth(nf).eq.'RPBE').and.(XCfunc(nf).eq.'GGA')) then

            write(6,'(a)')  'xc_check: GGA RPBE'
            if (icorr.ne.'rp'.and.nXCfunc.eq.1) 
     $          write(6,'(a,1x,2a)')
     .          'xc_check: WARNING: Pseudopotential generated with',
     $           trim(ps_string), " functional"

          elseif ((XCauth(nf).eq.'revPBE')
     $                  .and.(XCfunc(nf).eq.'GGA')) then

            write(6,'(a)')  'xc_check: GGA revPBE'
            if (icorr.ne.'rp'.and.nXCfunc.eq.1) 
     $          write(6,'(a,1x,2a)')
     .          'xc_check: WARNING: Pseudopotential generated with',
     $           trim(ps_string), " functional"

          elseif ((XCauth(nf).eq.'WC')
     $                  .and.(XCfunc(nf).eq.'GGA')) then

            write(6,'(a)')  'xc_check: GGA Wu-Cohen'
            if (icorr.ne.'wc'.and.nXCfunc.eq.1) 
     $          write(6,'(a,1x,2a)')
     .          'xc_check: WARNING: Pseudopotential generated with',
     $           trim(ps_string), " functional"

         else

           write(6,'(a)')
     . 'xc_check: ERROR: Exchange-correlation functional not allowed'
           write(6,'(a)') 'xc_check: ERROR: xc.functional= ',XCfunc(nf)
           write(6,'(a)') 'xc_check: ERROR: xc.authors= ',XCauth(nf)
           call die()

         endif

       enddo

       end subroutine xc_check
!
       subroutine comcore(is,a,b,rofi,chcore,nrval,nicore,flting)

C Generates the common blocks with the pseudo-core information.
C  D. Sanchez-Portal, Aug. 1998

           implicit none

           integer nrval, is
 
           real(dp) 
     .        chcore(nrmax), rofi(nrmax), flting, a ,b

           character nicore*4
          
CInternal variables
C
           integer nrcore,ir, nr, nmin, nmax, nn, itb
           real(dp) r2, chc, r, pi, delt, Rcore, dy,
     .     yp1, ypn
  
      logical   :: filterPCC
      real(dp)  :: kmaxpcc,grid,filterFactor,minnorm
      
           real(dp) eps
            parameter (eps=1.0d-6)  
          
C****NUMBER OF POINTS USED BY POLINT FOR THE INTERPOLATION***
C
          integer npoint
          parameter(npoint=4)
C
          pi=acos(-1.0d0) 
C
         coretab(1,2,is)=0

         if (flting.lt.0.0d0) then
            do itb=1,ntbmax+1
               coretab(itb,1,is)=0.0d0
               coretab(itb,2,is)=0.0d0
            enddo

            RETURN

         endif

         if(nicore.eq.'nc  ') then
            do itb=1,ntbmax+1
               coretab(itb,1,is)=0.0d0
               coretab(itb,2,is)=0.0d0
            enddo
            RETURN
         endif

            coretab(1,2,is)=1
            nrcore=0 
            chcore(1)=0.d0
            do ir=nrval,2,-1
              r=rofi(ir)
              r2=4.0d0*pi*r*r
              chc=chcore(ir)/r2

              if((chc.gt.eps).and.(nrcore.eq.0)) then
                  nrcore=ir+1
                  Rcore=rofi(nrcore)
                  goto 10
              endif
            enddo
10         continue

           write(6,'(a,f10.6)') 
     .         'comcore: Pseudo-core radius Rcore=',Rcore

C*TABLE WITH THE PSEUDO_CORE DATA

            delt=Rcore/(dble(ntbmax-1)+1.0d-20)  

            if(delt.gt.deltmax) then 
              write(6,'(a)') 
     .    'comcore: WARNING It might be a good idea to increase' 
              write(6,'(a)')
     .    'comcore: WARNING parameter ntbmax (in file atmparams.f) '
              write(6,'(a,i6)')
     .    'comcore: WARNING to at least ntbmax = ', 
     .            nint(Rcore/deltmax)+2 
            endif 

      filterPCC = fdf_boolean("PCC.Filter",.false.)

      if (filterPCC) then
         
         !Filter the real chcore
         do ir=2,nrval
            r=rofi(ir)
            r2=4.0d0*pi*r*r
            chcore(ir)=chcore(ir)/r2
         enddo

         minnorm = fdf_double("PCC.minnorm",0.999_dp)
         kmaxpcc = fdf_double("PCC.FilterCutoff",-10.0_dp)
         if (kmaxpcc .eq. -10.0)then
            filterFactor = fdf_double("PCC.FilterFactor",1.0_dp)
            grid = fdf_physical('MeshCutoff',100.0_dp,'Ry')
            kmaxpcc = filterFactor*sqrt(grid)
         endif
         write(6,"(A,f8.3)")"Corecharge: Filtered PCC cutoff" //
     $    " (Bohr^-1):",kmaxpcc
         call filter(0,nrcore,rofi(1:nrcore), 
     $        chcore(1:nrcore),kmaxpcc,0,minnorm)   
         
         !Store chcore*r2
         do ir=2,nrval
            r=rofi(ir)
            r2=4.0d0*pi*r*r
            chcore(ir)=chcore(ir)*r2
         enddo

      endif
   
            coretab(1,1,is)=delt
            do itb=2,ntbmax
              r=delt*(itb-1)
              nr=nint(log(r/b+1.0d0)/a)+1
              nmin=max(1,nr-npoint)
              nmax=min(nrcore,nr+npoint)
              nn=nmax-nmin+1
              call polint(rofi(nmin),chcore(nmin),nn,r,chc,dy)
              r2=4.0d0*pi*r*r

              coretab(itb+1,1,is)=chc/r2
            enddo
            coretab(2,1,is)=coretab(3,1,is)

C****TABLE WITH THE SECOND DERIVATIVE OF THE PSEUDO_CORE***

            yp1=huge(1.d0)
            ypn=huge(1.d0)
            call SPLINE( delt, coretab(2:1+ntbmax,1,is), ntbmax,
     &                   yp1,ypn,coretab(2:1+ntbmax,2,is))

           end  subroutine comcore
!
           subroutine comlocal(is,a,b,rofi,chlocal,nchloc,flting)

C Generates the common blocks with the local-pseudopotential charge
C  D. Sanchez-Portal, Aug. 1998

           implicit none

           integer nchloc, is
 
           real(dp) 
     .        chlocal(nrmax), rofi(nrmax), a ,b, flting

          
CInternal variables
C
           integer nr, nmin, nmax, nn, itb
           real(dp) chc, r, delt, Rchloc, dy,
     .     yp1, ypn
  
C****NUMBER OF POINTS USED BY POLINT FOR THE INTERPOLATION***
C
          integer npoint
          parameter(npoint=4)
C
        if(flting.gt.0.0d0) then 
  
          Rchloc=rofi(nchloc)


          delt=Rchloc/(dble(ntbmax-1)+1.0d-20) 
          
          if(delt.gt.deltmax) then
              write(6,'(a)')
     .    'comlocal: WARNING It might be a good idea to increase'
              write(6,'(a)')
     .    'comlocal: WARNING parameter ntbmax (in file atmparams.f) '
              write(6,'(a,i6)')
     .    'comlocal: WARNING to at least ntbmax = ', 
     .        nint(Rchloc/deltmax)+2 
          endif 

          chloctab(1,1,is)=delt
          chloctab(1,2,is)=1.0d0
          do itb=1,ntbmax
             r=delt*(itb-1)
             nr=nint(log(r/b+1.0d0)/a)+1
             nmin=max(1,nr-npoint)
             nmax=min(nchloc,nr+npoint)
             nn=nmax-nmin+1
             call polint(rofi(nmin),chlocal(nmin),nn,r,chc,dy)

             chloctab(itb+1,1,is)=chc
          enddo

          chloctab(2,1,is)=chloctab(3,1,is)

C****TABLE WITH THE SECOND DERIVATIVE OF THE LOCAL-PSEUDOTENTIAL***
C***CHARGE DENSITY****

         yp1=huge(1.d0)
         ypn=huge(1.d0)

         call SPLINE( delt, chloctab(2:1+ntbmax,1,is), ntbmax,
     &                yp1, ypn, chloctab(2:1+ntbmax,2,is) )

        elseif( flting.lt.0.0d0) then 
 

            do itb=1,ntbmax+1
               chloctab(itb,1,is)=0.0d0
               chloctab(itb,2,is)=0.0d0
            enddo

        endif 
        end subroutine comlocal

!
!       New routine to deal with Vlocal
!
        subroutine com_vlocal(is,a,b,rofi,red_vlocal,nvlocal,
     $                        zval,flting)

! Save vlocal (actually r*vlocal +2*Zval)

        implicit none

        integer nvlocal, is
 
        real(dp)  red_vlocal(nrmax), rofi(nrmax),
     $                    a ,b, zval, flting

           integer nr, nmin, nmax, nn, itb
           real(dp) chc, r, delt, Rvlocal, dy,
     .     yp1, ypn
  
          integer npoint
          parameter(npoint=4)
C
        if(flting.gt.0.0d0) then 
  
          Rvlocal=rofi(nvlocal)

          delt=Rvlocal/(dble(ntbmax-1)+1.0d-20) 
          
          if(delt.gt.deltmax) then
              write(6,'(a)')
     .    'comlocal: WARNING It might be a good idea to increase'
              write(6,'(a)')
     .    'comlocal: WARNING parameter ntbmax (in file atmparams.f) '
              write(6,'(a,i6)')
     .    'comlocal: WARNING to at least ntbmax = ', 
     .        nint(Rvlocal/deltmax)+2 
          endif 

          vlocaltab(1,1,is)=delt
          vlocaltab(1,2,is)=1.0d0
          do itb=1,ntbmax
             r=delt*(itb-1)
             nr=nint(log(r/b+1.0d0)/a)+1
             nmin=max(1,nr-npoint)
             nmax=min(nvlocal,nr+npoint)
             nn=nmax-nmin+1
             call polint(rofi(nmin),red_vlocal(nmin),nn,r,chc,dy)
             vlocaltab(itb+1,1,is)=chc
          enddo
!
!         Note rigorous limit at r=0...
!
          vlocaltab(2,1,is)= 2.0d0*zval
!
         yp1=huge(1.d0)
         ypn=huge(1.d0)

         call SPLINE( delt, vlocaltab(2:1+ntbmax,1,is), ntbmax,
     &                yp1, ypn, vlocaltab(2:1+ntbmax,2,is) )

        elseif( flting.lt.0.0d0) then 
 
            do itb=1,ntbmax+1
               vlocaltab(itb,1,is)=0.0d0
               vlocaltab(itb,2,is)=0.0d0
            enddo

        endif 

        end subroutine com_vlocal

!
       subroutine new_specie(iz,lmxkb, 
     .  nkbl, lmxo,
     .  nzeta, atm_label,
     .  npolorb, semic, nsemic, cnfigmx,
     .  is, new, no, nkb)

       implicit none

       integer
     .  iz, lmxkb, lmxo, nzeta(0:lmaxd,nsemx), npolorb(0:lmaxd,nsemx),
     .  nsemic(0:lmaxd), is, no, nkb, nkbl(0:lmaxd), cnfigmx(0:lmaxd)
      
       character 
     .   atm_label*20 

       logical new, semic
     
C      Internal and common variables

       integer  l, lmax, nzetamax, nkblmx,nsm, nsm_max
 
C**ADDING A NEW SPECIES TO THE LIST**
           
          new=.true.  
          if(iz.lt.0) lmxkb=0  
          lmax=max(lmxo,lmxkb) 

          if (lmax.gt.lmaxd) then 
              write(6,"(2a,i4)") 
     .        'new_specie: ERROR: Parameter lmaxd must be increased ',
     .        'to at least ', lmax 
              call die()
          endif                

          nzetamax=0
          nsm_max=0
          do l=0,lmxo
            nsm_max=max(nsemic(l)+1,nsm_max)
            do nsm=1, nsemic(l)+1
              nzetamax=max(nzeta(l,nsm),nzetamax) 
              nzetamax=max(npolorb(l,nsm),nzetamax) 
            enddo
          enddo

          if (nsm_max.gt.nsemx) then 
             write(6,"(2a,i4)")
     .        'new_specie: ERROR: Parameter nsmx must be increased ',
     .        'to at least ', nsm_max-1
             call die
          endif

          if (nzetamax.gt.nzetmx) then 
             write(6,"(2a,i4)") 
     .        'new_specie: ERROR: Parameter nzetmx must be increased ',
     .        'to at least ', nzetamax
             call die
          endif                
!
          nkblmx= maxval(nkbl(0:lmxkb))
           
          if (nkblmx.gt.nkbmx) then
             write(6,"(2a,i4)")
     .        'new_specie: ERROR: Parameter nkbmx must be increased ',
     .        'to at least ', nkblmx
             call die
          endif


          izsave(is)=iz
          lmxosave(is)=lmxo
          lmxkbsave(is)=lmxkb
          label_save(is)=atm_label
          semicsave(is)=semic   
            
          cnfigtb(:,:,is) = 0

           if (iz.ne.-100) then  
            do l=0,lmxo
              do nsm=1, nsemic(l)+1
                 cnfigtb(l,nsm,is)=cnfigmx(l)-(nsemic(l)+1)+nsm
              enddo
            enddo 
            do l=min(lmaxd,lmxo+1),lmaxd
               cnfigtb(l,1,is)=l+1         ! AG*** Why??
                                           ! To deal with polarization
                                           ! orbitals beyond lmxo?
! BUT, what happens to "in-between" polarization orbitals?
! 
!
            enddo
           else
!
!          Arbitrary "n" for "semicore" bessel states...
!
            do l=0,lmaxd
               do nsm=1,nsemic(l)+1
                  cnfigtb(l,nsm,is)=nsm
               enddo
            enddo 
           endif   ! iz.ne.-100

          no=0
          do l=0,lmxo
             nsemicsave(l,is)=nsemic(l)
             do nsm=1,nsemic(l)+1
               nzetasave(l,nsm,is)=nzeta(l,nsm)
               no=no+(2*l+1)*nzeta(l,nsm)
             enddo 
           enddo   

           nkb=0
           do l=0,lmxkb
              nkblsave(l,is)=nkbl(l)
              nkb=nkb+(2*l+1)*nkbl(l)
           enddo  

           do l=0,lmxo
             do nsm=1,nsemic(l)+1
               npolorbsave(l,nsm,is)=npolorb(l,nsm)
               no=no+(2*(l+1)+1)*npolorb(l,nsm)
             enddo
           enddo 

           nomax(is)=no
           if(iz.lt.0) nkb=0
           nkbmax(is)=nkb  

        end subroutine new_specie

!
        subroutine read_vps(lmxo, lmxkb,
     .             nrval,a,b,rofi,drdi,s,vps,
     .             rho, chcore, zval, chgvps,
     .             nicore, irel, icorr,basp)

        use basis_specs, only: restricted_grid
        use basis_specs, only: rmax_radial_grid

C   Read the file generated by the pseudopotential generation code.
C   Written by D. Sanchez-Portal, Aug. 1998
C**

           implicit none 

       type(basis_def_t), pointer   :: basp
      
           real(dp)
     .        rofi(nrmax), drdi(nrmax), s(nrmax), vps(nrmax,0:lmaxd),
     .        rho(nrmax), chcore(nrmax) 

           real(dp)
     .         a, b, zval
           
           integer  nrval, lmxo, lmxkb 

           character nicore*4, irel*3, icorr*2

C*Internal variables ****
           
           real(dp) 
     .     ve(nrmax)

           real(dp) 
     .        ea, rpb, chgvps

           integer  
     .        nr, nodd, lmax, linput, npotd, npotu,
     .        ndown, l, ir, i, itext
           character 
     .         orb*2,method(6)*10,text*70

           type(pseudopotential_t), pointer :: vp

           vp => basp%pseudopotential

             icorr = vp%icorr
             irel = vp%irel
             nicore = vp%nicore
             method = vp%method
             text = vp%text
             npotd = vp%npotd
             npotu = vp%npotu
             nr = vp%nr
             b = vp%b
             a = vp%a
             zval = vp%zval


           linput=max(lmxo,lmxkb)
           lmax=min(npotd-1,linput)


           if (lmax.lt.linput) then
               write(6,'(a)') 
     .          'read_vps: ERROR: You must generate a pseudopotential'
               write(6,'(a,i4)') 
     .          'read_vps: ERROR: for each L up to ',linput
             call die
           endif

           nrval=nr+1 
           if (rmax_radial_grid /= 0.0_dp) then
              nrval = nint(log(rmax_radial_grid/b+1.0d0)/a)+1
              write(6,"(a,f10.5,i5)")
     $             "Maximum radius (at nrval) set to ",
     $             rmax_radial_grid, nrval
           endif

           if (restricted_grid) then
              nodd=mod(nrval,2)
              nrval=nrval-1+nodd ! Will be less than or equal to vp%nrval
           endif

           if(nrval.gt.nrmax) then
               write(6,'(a,i4)')
     .     'read_vps: ERROR: Nrmax must be increased to at least',nrval
             call die
           endif

!=======
             write(6,'(/,a)')
     .           'read_vps: Pseudopotential generation method:'
             write(6,'(7a)')   
     .            'read_vps: ',method(1),(method(i),i=3,6)
          
C We are going to find the charge configuration
C used for the pseudopotential generation using the information given in
C the 'text' variable.

           chgvps = vp%gen_zval

           write(6,'(a,f10.5)') 'Total valence charge: ', chgvps

           if (nicore.ne.'nc  ') then
           write(6,'(/,a)')
     .      'read_vps: Pseudopotential includes a core correction:'

             if(nicore.eq.'pcec') then
               write(6,'(a)') 'read_vps: Pseudo-core for xc-correction'
             elseif(nicore.eq.'pche') then
               write(6,'(a)') 
     .            'read_vps: Pseudo-core for hartree and xc-correction'
           write(6,'(a)') 'Siesta cannot use this pseudopotential'
           write(6,'(a)') 'Use option pe instead of ph in ATOM program'
           call die()
             elseif(nicore.eq.'fcec') then
               write(6,'(a)') 'read_vps: Full-core for xc-correction'
             elseif(nicore.eq.'fche') then
               write(6,'(a)') 
     .            'read_vps: Full-core for hartree and xc-correction'
           write(6,'(a)') 'Siesta cannot use this pseudopotential'
           write(6,'(a)') 'Use option pe instead of ph in ATOM program'
           call die()
             endif

           endif

!          Radial mesh

           rofi(1:nrval) = vp%r(1:nrval)

C    Calculate drdi and s 
C    drdi is the derivative of the radial distance respect to the mesh index
C    i.e. rofi(ir)= b*[ exp( a*(i-1) ) - 1 ] and therefore 
C    drdi=dr/di =a*b*exp(a*(i-1))= a*[rofi(ir)+b] 

           rpb=b
           ea=exp(a)
           do ir=1,nrval
             drdi(ir)=a*rpb
             s(ir)=sqrt(a*rpb)
             rpb=rpb*ea
           enddo 


!          Ionic pseudopotentials   (Only 'down' used)

           do 20 ndown=1,lmax+1
               l = vp%ldown(ndown) 
               if(l.ne.ndown-1) then
                  write(6,'(a)')
     . 'atom: Unexpected angular momentum  for pseudopotential'
                  write(6,'(a)')
     . 'atom: Pseudopotential should be ordered by increasing l'
               endif
               vps(1:nrval,l) = vp%vdown(ndown,1:nrval)
               do ir=2,nrval
                  vps(ir,l)=vps(ir,l)/rofi(ir)
               enddo
               vps(1,l) = vps(2,l)     ! AG

  20       continue


!          Core and valence charge density

            chcore(1:nrval) = vp%chcore(1:nrval)
            rho(1:nrval) = vp%chval(1:nrval)


C OBTAIN AN IONIC-PSEUDOPOTENTIAL IF CORE CORRECTION FOR HARTREE
C POTENTIAL
!
! AG: OBSOLETE, as the program will stop in this case.

        if((nicore.eq.'pche').or.(nicore.eq.'fche')) then
            call vhrtre(chcore,ve,rofi,drdi,s,nrval,a)
            do l=0,lmax
              do ir=2,nrval
                vps(ir,l)=vps(ir,l)+ve(ir)
              enddo
              vps(1,l) = vps(2,l)    ! AG
            enddo
         endif

         return

 5000    continue
         write(6,*)
     $      'ERROR: You are using an old pseudopotential file.',
     $      ' Siesta needs a newer version.'
         call die
          
         end subroutine read_vps
!
               subroutine comKB(is,a,b,rofi,proj,
     .            l,ikb,rc,ekb,nrc)
C*
C  Creates the common block with all the information about the 
C  Kleinman-Bylander projectors.
C  Written by D. Sanchez-Portal, Aug. 1998.
C  Modified by DSP to allow more than one projector per l, July 1999.
C*
              
               implicit none

               integer l, nrc,is, ikb 

               real(dp) rc, ekb, proj(nrmax), a, b,
     .            rofi(nrmax)  

C****Internal variables
C
             integer indx, itb, nr, nmax, nmin, nn, il
             real(dp) delt, r, vphi, dy, yp1, ypn
             
****NUMBER OF POINTS USED BY POLINT FOR THE INTERPOLATION***
C
          integer npoint
          parameter(npoint=4)
C
             rctb(ikb,l,is)=rc 
 
CINTERPOLATION TO GENERATE TABLES WITH KB PROJECTORS**
C
             indx=0
             do il=0,l-1
               indx=indx+nkblsave(il,is)
             enddo 
             indx=indx+ikb
             if(ikb.gt.nkblsave(l,is)) then
                write(6,'(/,2a,i3,a,i3)')
     .         'comKB: ERROR: Maximum number of KB projectors',
     .         ' for l=',l,' must be', nkblsave(l,is)
              call die 
             endif 

             delt=rc/(dble(ntbmax-1)+1.0d-20) 
          if(delt.gt.deltmax) then
              write(6,'(a)')
     .    'comKB: WARNING It might be a good idea to increase'
              write(6,'(a)')
     .    'comKB: WARNING parameter ntbmax (in file atmparams.f) '
              write(6,'(a,i6)')
     .    'comKB: WARNING to at least ntbmax = ', 
     .        nint(Rc/deltmax)+2
           endif

             table(1,-indx,is)=delt
             table(2,-indx,is)=ekb

             do itb=1,ntbmax-1
                r=delt*(itb-1)
                nr=nint(log(r/b+1.0d0)/a)+1
                nmin=max(1,nr-npoint)
                nmax=min(nrc,nr+npoint)
                nn=nmax-nmin+1
                call polint(rofi(nmin),proj(nmin),nn,r,vphi,dy)
                table(itb+2,-indx,is)=vphi
            enddo 
            table(ntbmax+2,-indx,is)=0.0d0
C
C****

C****TABLE WITH THE SECOND DERIVATIVE ***
C
            yp1=huge(1.d0)
            ypn=huge(1.d0)

            call SPLINE( delt, table(3:2+ntbmax,-indx,is), ntbmax,
     &                   yp1, ypn, tab2(1:ntbmax,-indx,is) )

            end subroutine comkb
!
               subroutine KBgen(is, a,b,rofi,drdi,s, 
     .         vps, vlocal, ve, nrval, Zval, lmxkb, 
     .         nkbl, erefkb, nkb)

               use basis_specs, only: restricted_grid

C****
C Call routines for 1) the generation of the Kleinman-Bylander projectors,
C 2) Cheking for the presence of ghost states and 3), the storage of 
C all the information in the corresponding common blocks.
C
C  Written D. Sanchez-Portal, Aug. 1998.
C  Modified by DSP to allow more than one projector per l, July 1999.
C****

               implicit none

               real(dp) 
     .            a, b, rofi(nrmax), vps(nrmax,0:lmaxd),
     .            drdi(nrmax), s(nrmax), ve(nrmax),vlocal(nrmax),
     .            Zval, erefkb(nkbmx,0:lmaxd)

               integer
     .           nrval, lmxkb, nkb, is, nkbl(0:lmaxd)

C***Internal variables**

               integer 
     .           l,nprin, nnodes, ighost, nrwf, ikb, ir,
     .           nrc
               real(dp)
     .           rc(nkbmx,0:lmaxd), dkbcos(nkbmx,0:lmaxd),
     .           ekb(nkbmx,0:lmaxd)
               
               real(dp)
     .           rphi(nrmax,nkbmx), rmax, dnrm, 
     .           proj(nrmax)
                 
C  The atomic wavefunctions and/or its energy derivatives are
C  calculated only inside a sphere of radius Rmax. To define the
C  KB projectors they will not be needed very far from the nucleus,
C  and this limitation simplifies the handling of not bound states
C 
         parameter (Rmax=6.0d0)
C
         save ighost
         data ighost / 0 /

         nrwf=nint(log(Rmax/b+1.0d0)/a)+1
         nrwf=min(nrwf,nrval)
         if (restricted_grid)  nrwf=nrwf+1-mod(nrwf,2)
         
         do l=0,lmxkb
            do ir=1,nrmax
               do ikb=1,nkbmx
                  rphi(ir,ikb)=0.0d0
               enddo 
               proj(ir)=0.0d0
            enddo
            do ikb= 1, nkbl(l)
CAtomic wavefunctions and eigenvalues*
C****for the construction of the KB projectors**
 
C If the reference energies have not been specifed, the eigenstates
C with the condition of being zero at r(nrval) will be used.
C       
           if(erefkb(ikb,l).ge.1.0d3) then             
              nnodes=ikb
              nprin=l+1
              call schro_eq(Zval,rofi,vps(1,l),ve,s,drdi,
     .                      nrval,l,a,b,nnodes,nprin,
     .                      erefkb(ikb,l),rphi(1,ikb)) 
C Normalization of the eigenstates inside a sphere of radius Rmax
              dnrm=0.0d0
              do ir=1,nrwf
                dnrm=dnrm+drdi(ir)*rphi(ir,ikb)**2
              enddo 
              dnrm=sqrt(dnrm)
              do ir=1,nrwf
                 rphi(ir,ikb)=rphi(ir,ikb)/dnrm
              enddo 
C
           elseif((erefkb(ikb,l).le.-1.0d3).and.
     .       (ikb.gt.1) ) then 
C If the energy is specified to be 1000 Ry, the energy derivative
C of the previous wavefunction will be used
C
              call energ_deriv(a,rofi,rphi(1,ikb-1),vps(1,l),
     .                         ve,drdi,nrwf,l,erefkb(ikb-1,l),
     .                         rphi(1,ikb),nrval)
              erefkb(ikb,l)=0.0d0
C
           else 
C If the reference energies have been specified, we just use them
C
              call rphi_vs_e(a,b,rofi,vps(1,l),
     .                       ve,nrval,l,erefkb(ikb,l),
     .                       rphi(1,ikb),Rmax)
C
           endif 
C
C*

C***GHOST ANALYSIS****
C 
            if(nkbl(l).eq.1) then
             call ghost(Zval,rofi,vps(:,l),vlocal,
     .        ve,s,drdi,nrval,l,a,b,nrwf,
     .        erefkb(ikb,l),rphi(:,ikb),ighost)
            else 
             if (ikb.eq.1)
     .        write(6,'(a,i3,/a)') 
     .         'KBgen: More than one KB projector for l=',l,
     .         'KBgen: ghost states analysis will be not performed'
            endif

C***KB Projectors
C
          call KBproj(rofi,drdi,vps(1,l),vlocal,nrwf,l,rphi(1,ikb),
     .                 dkbcos(ikb,l),ekb(ikb,l),proj,nrc)

C
C*

          rc(ikb,l)=rofi(nrc)
 
CCommon block with the information about the  KB projectors***
C
          call comKB(is,a,b,rofi,proj,
     .                 l,ikb,rc(ikb,l),ekb(ikb,l),nrc)
C
C***

          enddo 
         enddo   
        
         if (ighost.eq.1) then
            write(6,"(2a)")'KBgen: WARNING: ',
     .            'Ghost states have been detected'
            write(6,"(2a)")'KBgen: WARNING: ',
     .            'Some parameter should be changed in the '
            write(6,"(2a)")'KBgen: WARNING: ',
     .            'pseudopotential generation procedure.'
          call die
         endif


         write(6,'(/,a)')'KBgen: Kleinman-Bylander projectors: '
         do l=0,lmxkb
           do ikb=1, nkbl(l)
              write(6,'(3x,a,i2,4(3x,a,f10.6))')
     .        'l=',l, 'rc=',rc(ikb,l), 'el=',erefkb(ikb,l), 
     .        'Ekb=',ekb(ikb,l),'kbcos=',dkbcos(ikb,l)
           enddo
         enddo

CTOTAL NUMBER OF KLEINMAN-BYLANDER PROJECTORS****
C
           nkb=0
           do l=0,lmxkb
              do ikb=1,nkbl(l)
                 nkb=nkb+(2*l+1) 
              enddo 
           enddo 
           write(6,'(/,a, i4)')
     .'KBgen: Total number of  Kleinman-Bylander projectors: ', nkb
C
        end  subroutine KBgen
!
              subroutine Basis_gen(Zval,is, iz, a,b,rofi,drdi,s,
     .                   vps, ve, vePAO, nrval, lmxo,nsemic, 
     .                   nzeta, rco, lambda, polorb,
     .                   basis_type, rphi, notot,
     $                   rinn, vcte, split_norm,atm_label)

C****
C Generates the basis set and stores all the information in the 
C correspoinding common blocks.
C
C Written by D. Sanchez-Portal, Aug. 1998.
C Modify by DSP, July 1999
C****


               implicit none

               real(dp)
     .            a, b, rofi(nrmax), vps(nrmax,0:lmaxd),
     .            drdi(nrmax), s(nrmax), ve(nrmax),
     .            rphi(nrmax,0:lmaxd,nsemx), rco(nzetmx,0:lmaxd,nsemx),
     .            lambda(nzetmx, 0:lmaxd,nsemx), vePAO(nrmax),
     .            Zval,rinn(0:lmaxd,nsemx),vcte(0:lmaxd,nsemx),
     $            split_norm(0:lmaxd,nsemx)
               character(len=*) atm_label
               integer
     .           nrval, lmxo, notot, is, iz, nzeta(0:lmaxd,nsemx),
     .           polorb(0:lmaxd,nsemx),nsemic(0:lmaxd)

               character
     .           basis_type*10

C***Internal variables**

                integer noPAO, noPOL

                real(dp) ePAO(0:lmaxd,nsemx)

             
                noPAO=0
                noPOL=0
                if(basis_type.eq.'split') then  

                 call SPLIT(Zval,is, iz, a,b,rofi,drdi,s,
     .                   vps, ve, vePAO, 
     .                   nrval, lmxo, nsemic,
     .                   nzeta, rco, lambda,
     .                   rphi, ePAO, noPAO,
     $                   rinn, vcte, split_norm,atm_label)
         
                elseif(basis_type.eq.'nodes') then 
 
                 call NODES(Zval,is, a,b,rofi,drdi,s,
     .                   vps, ve, vePAO,
     .                   nrval, lmxo, nsemic,
     .                   nzeta, rco, lambda,
     .                   rphi, ePAO, noPAO)

                elseif(basis_type.eq.'nonodes') then 
 
                 call NONODES(Zval,is, a,b,rofi,drdi,s,
     .                   vps, ve, vePAO,
     .                   nrval, lmxo, nsemic,
     .                   nzeta, rco, lambda,
     .                   rphi, ePAO, noPAO)
          
                elseif(basis_type.eq.'splitgauss') then 
      
                 call SPLITGAUSS(Zval,is, a,b,rofi,drdi,s,
     .                   vps, ve, vePAO,
     .                   nrval, lmxo, nsemic,
     .                   nzeta, rco, lambda,
     .                   rphi, ePAO, noPAO)

                endif            
 

C***Polarization orbitals*
C 
                  call POLgen(is, iz, a,b,rofi,drdi,
     .               ePAO,rphi,rco,vps,vePAO,
     .               polorb,lmxo,nsemic,noPOL,
     $               rinn, vcte,nrval, split_norm,atm_label) 
C
C

C*Total number of orbitals
C 
                   notot=noPAO+noPOL

          end subroutine basis_gen
!

           subroutine SPLIT(Zval,is, iz, a,b,rofi,drdi,s,
     .             vps,ve,vePAO,
     .             nrval,lmxo, nsemic,
     .             nzeta,rco,lambda, rphi, ePAO, norb,
     $             rinn,vcte,split_norm,atm_label) 

           use basis_specs, only:  restricted_grid
C****
C Calculates the atomic orbitals basis set, using the option SPLIT 
C for the generation of the augmentation orbitals.
C  Written by D. Sanchez-Portal, Aug. 1998
C  Modified by DSP, July 1999
C****

               implicit none

               real(dp)
     .         a, b, rofi(nrmax), vps(nrmax,0:lmaxd),
     .         drdi(nrmax), s(nrmax), ve(nrmax),
     .         rphi(nrmax,0:lmaxd,nsemx), rco(nzetmx,0:lmaxd,nsemx),
     .         lambda(nzetmx,0:lmaxd,nsemx), Zval,vePAO(nrmax),
     .         ePAO(0:lmaxd,nsemx),
     .         rinn(0:lmaxd,nsemx),vcte(0:lmaxd,nsemx),
     .         split_norm(0:lmaxd,nsemx)
               character(len=*) atm_label


               integer
     .           nrval, lmxo, is, iz, nzeta(0:lmaxd,nsemx),
     .           norb, nsemic(0:lmaxd)


C***Internal variables**

               integer
     .           l,nprin, nnodes, nodd, nrc, nsp, i, ir,indx,
     .           izeta, nmax, nmin, nn, nr, nrcomp, nsm, nrc1, 
     .           nrc2, nrc3, nrc4, ism

               real(dp)
     .           eigen(0:lmaxd), rc, split_table(nrmax),
     .           rnrm(nrmax), dnrm, phi,
     .           cons1, cons2, rnp, spln, eshift, 
     .           g(nrmax), r, el, ekin, rdummy,
     .           r1, r2, dfdi, d2fdi2, d2fdr2, dr,
     .           epot, epot2, rh, dy, eorb, eps, 
     .           over(nsemx), vsoft(nrmax), vePAOsoft(nrmax),
     $           exponent, dlt, d, dn, norm(nsemx), rcsan,
     $     kmax,grid,filterFactor,minnorm,
     $           spln_min
      real(dp),allocatable,dimension(:) :: forb !filtered orbital
      logical :: filterOrbitals
      logical :: new_split_code, fix_split_table, split_tail_norm

               parameter (dlt=0.60d0)

               character  filename*80, paste*80


C****NUMBER OF POINTS USED BY POLINT FOR THE INTERPOLATION***
C 
               integer  npoint 
               parameter(npoint=4)
C
C***READING THE ENERGY-SHIFT TO DEFINE THE CUT-OFF RADIUS OF ORBITALS***

           eshift=fdf_physical('PAO.EnergyShift',eshift_default,'Ry')

           split_tail_norm=fdf_boolean('PAO.SplitTailNorm',.false.)
           fix_split_table=fdf_boolean('PAO.FixSplitTable',.false.)
           if (split_tail_norm .or. fix_split_table) then
              new_split_code = .true.
           else
              new_split_code=fdf_boolean('PAO.NewSplitCode',.false.)
           endif

             norb=0 
             indx=0
!
!            LOOP over angular momenta
!
             do l=0,lmxo 
              
              do nsm=1,nsemic(l)+1 
                if(nzeta(l,nsm).gt.0) then
                    write(6,'(/A,I2)')
     .               'SPLIT: Orbitals with angular momentum L=',l 
                  goto 50
                endif
              enddo 

50            continue
              
              do nsm=1,nsemic(l)+1

                if(nzeta(l,nsm).gt.0) then
 
                    write(6,'(/A,I1,A)')
     .               'SPLIT: Basis orbitals for state ',
     .                cnfigtb(l,nsm,is), sym(l)

                  if(rco(1,l,nsm).lt.1.0d-5) then    
C**Automatic determination of the cut off radius for the PAOs***
CAtomic eigenvalues
C
                      nnodes=nsm
                      nprin=l+nsm
                      call schro_eq(Zval,rofi,vps(1,l),ve,s,drdi,
     .                  nrval,l,a,b,nnodes,nprin,
     .                  eigen(l),rphi(1,l,nsm))
C
C*Rc given by eshift****   
C
                       if(eigen(l).gt.0.0d0) then 
                          write(6,'(/A,I2,A)')
     .       'SPLIT: ERROR Orbital with angular momentum L=',l,
     .       ' not bound in the atom'
                         write(6,'(A)')
     .       'SPLIT: ERROR a cut off radius must be explicitely given' 
           call die
                       endif 
                       if(abs(eshift).gt.1.0d-5) then
                          el=eigen(l)+eshift
                          call rc_vs_e(a,b,rofi,vps(1,l),
     .                     ve,nrval,l,el,nnodes,rco(1,l,nsm))
                       else
                          rco(1,l,nsm)=rofi(nrval-2)
                       endif
                  write(6,'(/,A,/,A,f10.6,A)')
     .        'SPLIT: PAO cut-off radius determined from an',
     .        'SPLIT: energy shift=',eshift,' Ry'

                 endif  
C

C IF THE COMPRESSION FACTOR IS NEGATIVE OR ZERO THE ORBITALS ARE
C LEFT UNTOUCHED
             if(lambda(1,l,nsm).le.0.0d0) lambda(1,l,nsm)=1.0d0
!
! SOFT-CONFINEMENT 
!
C Calculate the soft confinement potential for a given atomic species
C and angular momentum 
!!AG** rco might not be the PAO cutoff radius! See below

              call build_vsoft(is,l,nsm,rinn(l,nsm),vcte(l,nsm),
     $                       a, b, rco(1,l,nsm),rofi,nrval,
     $                       vsoft,plot=.true.)

              do ir = 1, nrval
                vePAOsoft(ir) = vePAO(ir) + vsoft(ir)
              enddo

! END SOFT-CONFINEMENT 

              do izeta=1, nzeta(l,nsm)

C COMPRESSION FACTOR IS ONLY ACTIVE FOR THE INITIAL PAO WHEN USING
C SPLIT OPTION FOR THE GENERATION OF THE BASIS SET

                 lambda(izeta,l,nsm)=lambda(1,l,nsm)

!
!                 If rc is negative, treat it as a fractional value
!
                  if (rco(izeta,l,nsm) < 0.0_dp) then
                     if (izeta == 1) then
                        call die("rc < 0 for first-zeta orbital")
                     else
                        rco(izeta,l,nsm) = rco(1,l,nsm) *
     $                       (-rco(izeta,l,nsm))
                     endif
                  endif

                  rc=rco(izeta,l,nsm)/lambda(1,l,nsm)
                  nrc=nint(log(rc/b+1.0d0)/a)+1

! Note that rco is redefined here, to make it fall on an odd-numbered
! grid point.
!
                  if (restricted_grid) then
                     nodd=mod(nrc,2)
                     if(nodd.eq.0) then
                        nrc=nrc+1
                     endif
                  endif
                  rc=b*(exp(a*(nrc-1))-1.0d0)
                  rco(izeta,l,nsm)=rc*lambda(1,l,nsm)
                  
                if(izeta.eq.1) then 

!
!AG**    Soft-confinement setup might need to go here
!
C****Generate PAO orbitals for the first shell of the basis set***
C 
                   nnodes=nsm
                   nprin=l+nsm
                   call schro_eq(Zval,rofi,vps(1,l),vePAOsoft,s,drdi,
     .                  nrc,l,a,b,nnodes,nprin,
     .                  eorb,rphi(1,l,nsm)) 
                   dnrm=0.0d0
                   do ir=2,nrc
                      phi=rphi(ir,l,nsm)
                      dnrm=dnrm+drdi(ir)*phi*phi
                      rnrm(ir)=dnrm 
                      g(ir)=rphi(ir,l,nsm)/(rofi(ir)**(l+1))
                      write(99,*) ir, rofi(ir), phi, rnrm(ir)
                   enddo 
                   g(1)=g(2)         

                   if (split_tail_norm) then
                      write(*,"(a)") "Split based on tail norm"
                      ! Follows the criterion of the JPC paper, but
                      ! with more contrast 
                      !(use the actual math norm (sqrt(int(f^2)))

                      split_table(1:nrc) =
     $                     sqrt(max(1.0_dp-rnrm(1:nrc),0.0_dp))

                   else
                   !  Do a full scan of the old method
                   !  (norm of tail+parabola)
                      call split_scan(nrc, rofi, drdi, l,
     $                        rphi(1:,l,nsm),rnrm,atm_label,
     $                        split_table,fix_split_table)
                   endif

                else              ! izeta > 1...
       
C***Cut-off radius for double-Z, triple-Z,..., if it is set to**** 
C***zero in the input then it is calculated from the splitnorm**** 
C*** parameter *

                if(rco(izeta,l,nsm).gt.rco(1,l,nsm)) then
                    write(6,'(/,A)') 
     . 'SPLIT: ERROR: SPLIT OPTION FOR BASIS SET '
                    write(6,'(A)')  
     . 'SPLIT: ERROR: Rc FOR DOUBLE-Z, TRIPLE-Z,... SHOULD BE SMALLER '
                    write(6,'(A)') 
     . 'SPLIT: ERROR:  THAN THAT OF THE INITIAL PAO !!!!!'
                  call die
                endif
                                  
                
            if(rco(izeta,l,nsm).gt.1.0d-5) then

               rc=rco(izeta,l,nsm)/lambda(1,l,nsm)
               nrc=nint(log(rc/b+1.0d0)/a)+1

               call fit_parabola(nrc,rofi,drdi,rphi(1:,l,nsm),
     $              l, cons1, cons2, rnp)
               spln=split_table(nrc)
C**
                do i=1,izeta-1
                 if(abs(rco(izeta,l,nsm)-rco(i,l,nsm)).lt.1.0d-5) then
                   write(6,'(/,A,I2,A,I2,A,I2)')
     .            'SPLIT: WARNING: Split-orbital with zeta=',izeta,
     .            ' and zeta=',i,' are identical for l=',l
                   call die()
                 endif
                enddo

            else        ! Generate multiple zeta
                        ! using split-norm parameters

            rc=rco(1,l,nsm)/lambda(1,l,nsm)
            nrc=nint(log(rc/b+1.0d0)/a)+1

            spln = split_norm(l,nsm)
            if(izeta.gt.2) then
              spln=split_norm(l,nsm)/(2.0d0*(izeta-2) )
            endif
!
!!            Maybe fix on the fly in the future?
!!            if (spln < minval(split_table(1:nrc))
!!                ...  ! maybe fix split_table

            if (new_split_code) then
               call find_split_location(nrc,rofi,drdi,
     $              rphi(1:,l,nsm),split_table,
     $              l,spln,cons1,cons2,nsp)
            else
               spln_min = minval(split_table(1:nrc))
               if (spln < spln_min) then
                  write(6,"(a,f8.5,a,f8.5)")
     $            "WARNING: Minimum split_norm parameter: ",
     $            spln_min, ". Will not be able to generate "
     $            // "orbital with split_norm = ", spln
                  call die("See manual for new split options")
               endif
               call parabola(a,b,nrc,rphi(1,l,nsm),rnrm,
     .              l,spln,cons1,cons2,nsp)
            endif

C***Cut-off radius for the split orbital with a desired norm*
         nrc=nsp
         rco(izeta,l,nsm)= rofi(nsp) *lambda(izeta,l,nsm)
C

             do i=1,izeta-1
                if(abs(rco(izeta,l,nsm)-rco(i,l,nsm))
     .                                   .lt.1.0d-5) then
                   write(6,'(/,A,I2,A,I2,A,I2)')
     .            'SPLIT: WARNING: Split-orbital with zeta=',izeta,
     .            ' and zeta=',i,' are identical for l=',l
                endif
             enddo

            endif

            do ir=1,nrval
C***parabolic split***
               r=rofi(ir)
               if (ir.ge.nrc) then
                 g(ir)=0.0d0
               else
                 ! Store first-zeta minus second-zeta in g
                 g(ir)=-(cons1*r**2+cons2)*r**(l+1)+rphi(ir,l,nsm) 
               endif
C*
            enddo 

C***Orthogonalize to the inner shells if present***
C***In some cases we have to use a "smooth" version of the orginal
C***PAO orbitals in this orthogonalization process***************
             if(nsm.gt.1) then
              do ism=1,nsm-1
               rc=rco(1,l,ism)/lambda(1,l,ism)
!
               nrc4=nint(dlog(rc/b+1.0d0)/a)+1
               rc=rc-dlt
               nrc3=nint(dlog(rc/b+1.0d0)/a)+1
               rc=rofi(nrc3)
               cons1=rphi(nrc3,l,ism)
               d=dsqrt(cons1**2+dlt**2)
               rc=rc+d+dlt
               nrc2=nint(dlog(rc/b+1.0d0)/a)+1
               d=rofi(nrc2)-rofi(nrc3)
               rc=rofi(nrc2)
               cons2=0.5d0*(cons1**2+(d)**2)/cons1
               dnrm=0.0d0
               dn=0.0d0
!
               do ir=1, nrc2
                  r=rofi(ir)
                  if(ir.gt.nrc3.and.nrc.gt.nrc4) then
                    rh=cons2-
     .              cons1*dsqrt(cons2**2-(r-rc)**2)/dabs(cons1)
                  else
                    rh=rphi(ir,l,ism)
                  endif
c                  write(68,*) r,rh,rphi(ir,l,ism)
                  dnrm=dnrm+drdi(ir)*rh*g(ir)
                  dn=dn+drdi(ir)*rh*rh
               enddo
               over(ism)=dnrm
               norm(ism)=dn
              enddo
!
              nrc1=nrc
              do ism=1,nsm-1
               rc=rco(1,l,ism)/lambda(1,l,ism)
!
               nrc4=nint(dlog(rc/b+1.0d0)/a)+1
               rc=rc-dlt
               rc=rco(1,l,ism)/lambda(1,l,ism)-dlt
               nrc3=nint(dlog(rc/b+1.0d0)/a)+1
               rc=rofi(nrc3)
               cons1=rphi(nrc3,l,ism)
               d=dsqrt(cons1**2+dlt**2)
               rc=rc+d+dlt
               nrc2=nint(dlog(rc/b+1.0d0)/a)+1
               d=rofi(nrc2)-rofi(nrc3)
               rc=rofi(nrc2)
               cons2=0.5d0*(cons1**2+(d)**2)/cons1
               if(nrc4.lt.nrc) nrc1=max(nrc1,nrc2)

               do ir=1, nrc2
                  r=rofi(ir)
                  if(ir.gt.nrc3.and.nrc.gt.nrc4) then
                    rh=cons2-
     .              cons1*dsqrt(cons2**2-(r-rc)**2)/dabs(cons1)
                  else
                    rh=rphi(ir,l,ism)
                  endif
                  g(ir)=g(ir)-over(ism)*rh/(norm(ism)+1.0d-20)
               enddo
              enddo


              if(nrc.ne.nrc1) then 
                 nrc=nrc1
                 rco(izeta,l,nsm)=
     .           b*(exp(a*(nrc1-1))-1.0d0)*lambda(izeta,l,nsm)
              endif 
            endif
              
            dnrm=0.0d0
            do ir=2,nrc-1
              r=rofi(ir)
C***parabolic split***
              phi=g(ir)/(r**(l+1))
C*

              dnrm=dnrm+drdi(ir)*(phi*r**(l+1))**2
              g(ir)=phi
            enddo
            g(1)=g(2)
            g(nrc)=0.0d0

            endif
              
C****Normalization of basis functions***
            eps=1.0d-4
            if(abs(dnrm-1.0d0).gt.eps) then
               do ir=1,nrc
                 g(ir)=g(ir)/sqrt(dnrm)
                 if(izeta.eq.1) then
                    rphi(ir,l,nsm)=rphi(ir,l,nsm)/sqrt(dnrm)
                    rnrm(ir)=rnrm(ir)/dnrm
                 endif
               enddo
            endif
C****  

C*Calculation of the mean value of kinetic and potential energy**
C    Potential and kinetic energy of the orbital before compression



             ekin=0.0d0
             do ir=2,nrc-1
                r=rofi(ir)
                r1=rofi(ir-1)
                r2=rofi(ir+1)
                d2fdi2=(g(ir-1)*r1**(l+1)+g(ir+1)*r2**(l+1)
     .                       -2.0d0*g(ir)*r**(l+1))
                dfdi=0.5d0*(g(ir+1)*r2**(l+1)-g(ir-1)*r1**(l+1))
                dr=drdi(ir)
                d2fdr2= ((-a)*dfdi +  d2fdi2)/dr**2
                ekin=ekin+
     .              dr*g(ir)*r**(l+1)*(-d2fdr2)
     .             +dr*l*(l+1)*(g(ir)*r**l)**2
             enddo


C Kinetic energy after compression

             ekin=ekin/(lambda(izeta,l,nsm)**2)

C Potential energy after compression

             nrcomp=nint(log(rco(izeta,l,nsm)/b+1.0d0)/a)+1
             epot=0.0d0
             epot2=0.0d0
             do ir=1,nrcomp
                r=rofi(ir)
                r2=r/lambda(izeta,l,nsm)
                nr=nint(log(r2/b+1.0d0)/a)+1
                nmin=max(1,nr-npoint)
                nmax=min(nrc,nr+npoint)
                nn=nmax-nmin+1
                call polint(rofi(nmin),g(nmin),nn,r2,rh,dy)
                rh=rh/sqrt(lambda(izeta,l,nsm)**(2*l+3))
                epot=epot+
     .          drdi(ir)*(ve(ir)+vps(ir,l))*(rh*r**(l+1))**2
                epot2=epot2+
     .          drdi(ir)*vps(ir,l)*(rh*r**(l+1))**2
             enddo
             eorb=ekin+epot



            if(izeta.eq.1) then  

             write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')
     .          'izeta =',izeta,
     .          'lambda =',lambda(izeta,l,nsm),
     .          'rc =',rco(izeta,l,nsm),
     .          'energy =',eorb  

                ePAO(l,nsm)=eorb

            elseif(izeta.gt.1) then 

            write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')
     .         'izeta =',izeta,
     .         'rmatch =',rco(izeta,l,nsm),
     .         'splitnorm =',spln,
     .         'energy =',eorb 
            if (spln < 0.05_dp) then
               write(6,"(a)")
     $           "* WARNING: effective split_norm is quite small."
     $              // " Orbitals will be very similar."
            endif

            endif 

          write(6,'(a25,f12.6)') 'kinetic =',ekin
          write(6,'(a25,f12.6)') 'potential(screened) =',epot
          write(6,'(a25,f12.6)') 'potential(ionic) =',epot2 

            norb=norb+(2*l+1)
            indx=indx+1
!Filter
                  filterOrbitals = fdf_boolean("PAO.Filter",.false.)

                  kmax = fdf_double("PAO.FilterCutoff",-10.0_dp)
                  if (kmax .eq. -10.0)then
                     filterFactor = fdf_double("PAO.FilterFactor",
     $                    0.7_dp)
                     grid = fdf_physical('MeshCutoff',100.0_dp,'Ry')
                     kmax = filterFactor*sqrt(grid)
                  endif          

                  if (filterOrbitals)then
                     
                     allocate(forb(1:nrc))
                     
                     do ir=1,nrc
                        if(l .eq. 0) then
                           forb(ir) = g(ir)
                        else
                           forb(ir) = g(ir)*rofi(ir)**(l)
                        endif
                     enddo

          
                     minnorm = fdf_double("PAO.minnorm",0.999_dp)
                 write(6,"(A,A,f8.3)")"paogen: Filtered orbital cutoff",
     $                    " (Bohr^-1):",kmax
                     call filter(l,nrc,rofi(1:nrc),forb(1:nrc),kmax,2,
     $                    minnorm)
             

c                    Store the filtered orbital in g
                     g = 0.0_dp                   
                     do ir=2,nrc
                        if(l .eq. 0)then
                           g(ir) = forb(ir) 
                        else
                           g(ir) = forb(ir)/rofi(ir)**l
                        endif
                     enddo
                     g(1) = g(2)    
                     deallocate(forb)
                     
                  endif !End of filtering

            call comBasis(is,a,b,rofi,g,l,
     .              rco(izeta,l,nsm),lambda(izeta,l,nsm),izeta,
     .              nsm,nrc,indx)

               
              enddo 
        
            call compress_PAO(a,b,rofi,rphi(1,l,nsm),
     .              rco(1,l,nsm),lambda(1,l,nsm))  
                       

             endif    
           enddo 
          enddo 
          
          end subroutine SPLIT
!
           subroutine NODES(Zval,is,a,b,rofi,drdi,s,
     .             vps,ve,vePAO,
     .             nrval,lmxo,nsemic,
     .             nzeta,rco,lambda, rphi, ePAO, norb) 

           use basis_specs, only: restricted_grid
C****
C Calculates the atomic orbitals basis set, using the option NODES
C for the generation of the augmentation orbitals.
C  Written by D. Sanchez-Portal, Aug. 1998
C   Modified by DSP, July 1999
C****

               implicit none

               real(dp)
     .         a, b, rofi(nrmax), vps(nrmax,0:lmaxd),
     .         drdi(nrmax), s(nrmax), ve(nrmax),
     .         rphi(nrmax,0:lmaxd,nsemx), 
     .         rco(nzetmx,0:lmaxd,nsemx),
     .         lambda(nzetmx, 0:lmaxd,nsemx), Zval,vePAO(nrmax),
     .         ePAO(0:lmaxd,nsemx)

               integer
     .           nrval, lmxo, is, nzeta(0:lmaxd,nsemx),
     .           norb,nsemic(0:lmaxd)


C***Internal variables**

               integer
     .           l,nprin, nnodes, nodd, nrc, ir, indx,
     .           izeta, nmax, nmin, nn, nr, nrcomp, nsm

               real(dp)
     .           eigen(0:lmaxd), rc,
     .           dnrm, phi, eshift,
     .           g(nrmax), r, el, ekin, 
     .           r1, r2, dfdi, d2fdi2, d2fdr2, dr,
     .           epot, epot2, rh, dy, eorb, eps


C****NUMBER OF POINTS USED BY POLINT FOR THE INTERPOLATION***
C 
               integer  npoint 
               parameter(npoint=4)
C
C***READING THE ENERGY-SHIFT TO DEFINE THE CUT-OFF RADIUS OF ORBITALS***

           eshift=fdf_physical('PAO.EnergyShift',eshift_default,'Ry')

             norb=0 
             indx=0
             do l=0,lmxo


              do nsm=1,nsemic(l)+1
                if(nzeta(l,nsm).gt.0) then
                    write(6,'(/A,I2)')
     .               'NODES: Orbitals with angular momentum L=',l
                  goto 50
                endif
              enddo

50            continue

              do nsm=1,nsemic(l)+1

                if(nzeta(l,nsm).gt.0) then

                    write(6,'(/A,I1,A)')
     .               'NODES: Basis orbitals for state ',
     .                cnfigtb(l,nsm,is), sym(l)

                  if(rco(1,l,nsm).lt.1.0d-5) then    
C**Automatic determination of the cut off radius for the PAOs***
CAtomic eigenvalues
C
                      nnodes=nsm
                      nprin=l+nsm
                      call schro_eq(Zval,rofi,vps(1,l),ve,s,drdi,
     .                  nrval,l,a,b,nnodes,nprin,
     .                  eigen(l),rphi(1,l,nsm))
C
C*Rc given by eshift****   
C                
                       if(eigen(l).gt.0.0d0) then
                          write(6,'(/A,I2,A)')
     .       'NODES: ERROR Orbital with angular momentum L=',l,
     .       ' not bound in the atom'
                         write(6,'(A)')
     .       'NODES: ERROR a cut off radius must be explicitely given' 
          call die
                       endif 
 
                       if(abs(eshift).gt.1.0d-5) then
                          el=eigen(l)+eshift
                          call rc_vs_e(a,b,rofi,vps(1,l),
     .                      ve,nrval,l,el,nnodes,rco(1,l,nsm))
                       else
                          rco(1,l,nsm)=rofi(nrval-2)
                       endif 

                  write(6,'(/,A,/,A,f10.6,A)')
     .         'NODES: PAO cut-off radius determinated from an',
     .         'NODES: energy shift=',eshift,' Ry'

                 endif  
C
C***


              do izeta=1, nzeta(l,nsm)

CIF THE COMPRESSION FACTOR IS NEGATIVE OR ZERO THE ORBITALS ARE***
CUNTOUCHED**
          if(lambda(izeta,l,nsm).le.0.0d0) lambda(izeta,l,nsm)=1.0d0 
C*
             if(abs(rco(izeta,l,nsm)).le.1.0d-5) then 
                 rco(izeta,l,nsm)=rco(1,l,nsm)
             endif 

                  rc=rco(izeta,l,nsm)/lambda(izeta,l,nsm)
                  nrc=nint(log(rc/b+1.0d0)/a)+1
                  if (restricted_grid) then
                     nodd=mod(nrc,2)
                     if(nodd.eq.0) then
                        nrc=nrc+1
                     endif
                  endif
                  rc=b*(exp(a*(nrc-1))-1.0d0)
                  rco(izeta,l,nsm)=rc*lambda(izeta,l,nsm)
C****Generate PAO orbitals with increasing number of nodes ***
C*for the different shells*
C 
                      nnodes=izeta+nsm
                      nprin=l+izeta+nsm
                      eorb=0.0d0 
                   call schro_eq(Zval,rofi,vps(1,l),vePAO,s,drdi,
     .              nrc,l,a,b,nnodes,nprin,
     .                eorb,g)  

                       dnrm=0.0d0
                       do ir=2,nrc
                          phi=g(ir) 
                          if(izeta.eq.1) rphi(ir,l,nsm)=phi
                          dnrm=dnrm+drdi(ir)*phi*phi
                          g(ir)=g(ir)/(rofi(ir)**(l+1))
                       enddo 
                       g(1)=g(2)         

C****Normalization of basis functions***
            eps=1.0d-4
            if(abs(dnrm-1.0d0).gt.eps) then
               do ir=1,nrc
                 g(ir)=g(ir)/sqrt(dnrm)
                 if(izeta.eq.1) then
                    rphi(ir,l,nsm)=rphi(ir,l,nsm)/sqrt(dnrm)
                 endif
               enddo
            endif
C****


C*Calculation of the mean value of kinetic and potential energy**
C    Potential and kinetic energy of the orbital before compression



             ekin=0.0d0
             do ir=2,nrc-1
                r=rofi(ir)
                r1=rofi(ir-1)
                r2=rofi(ir+1)
                d2fdi2=(g(ir-1)*r1**(l+1)+g(ir+1)*r2**(l+1)
     .                       -2.0d0*g(ir)*r**(l+1))
                dfdi=0.5d0*(g(ir+1)*r2**(l+1)-g(ir-1)*r1**(l+1))
                dr=drdi(ir)
                d2fdr2= ((-a)*dfdi +  d2fdi2)/dr**2
                ekin=ekin+
     .              dr*g(ir)*r**(l+1)*(-d2fdr2)
     .             +dr*l*(l+1)*(g(ir)*r**l)**2
             enddo


C Kinetic energy after compression

             ekin=ekin/(lambda(izeta,l,nsm)**2)

C Potential energy after compression

             nrcomp=nint(log(rco(izeta,l,nsm)/b+1.0d0)/a)+1
             epot=0.0d0
             epot2=0.0d0
             do ir=1,nrcomp
                r=rofi(ir)
                r2=r/lambda(izeta,l,nsm)
                nr=nint(log(r2/b+1.0d0)/a)+1
                nmin=max(1,nr-npoint)
                nmax=min(nrc,nr+npoint)
                nn=nmax-nmin+1
                call polint(rofi(nmin),g(nmin),nn,r2,rh,dy)
                rh=rh/sqrt(lambda(izeta,l,nsm)**(2*l+3))
                epot=epot+
     .          drdi(ir)*(ve(ir)+vps(ir,l))*(rh*r**(l+1))**2
                epot2=epot2+
     .          drdi(ir)*vps(ir,l)*(rh*r**(l+1))**2
             enddo
             eorb=ekin+epot



               write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')
     .          'izeta =',izeta,
     .          'lambda =',lambda(izeta,l,nsm),
     .          'rc =',rco(izeta,l,nsm),
     .          'energy =',eorb  


                 if(izeta.eq.1) ePAO(l,nsm)=eorb

          write(6,'(a25,f12.6)') 'kinetic =',ekin
          write(6,'(a25,f12.6)') 'potential(screened) =',epot
          write(6,'(a25,f12.6)') 'potential(ionic) =',epot2 

            norb=norb+(2*l+1)
            indx=indx+1
            call comBasis(is,a,b,rofi,g,l,
     .              rco(izeta,l,nsm),lambda(izeta,l,nsm),izeta,
     .              nsm,nrc,indx)


              enddo 
        
            call compress_PAO(a,b,rofi,rphi(1,l,nsm),
     .              rco(1,l,nsm),lambda(1,l,nsm))


              endif
             enddo     
           enddo 

           end subroutine nodes
!

           subroutine NONODES(Zval,is,a,b,rofi,drdi,s,
     .             vps,ve,vePAO,
     .             nrval,lmxo,nsemic,
     .             nzeta,rco,lambda, rphi, ePAO, norb) 

           use basis_specs, only: restricted_grid

C****
C Calculates the atomic orbitals basis set, using the option NONODES
C for the generation of the augmentation orbitals.
C  Written by D. Sanchez-Portal, Aug. 1998
C****

               implicit none

               real(dp)
     .         a, b, rofi(nrmax), vps(nrmax,0:lmaxd),
     .         drdi(nrmax), s(nrmax), ve(nrmax),
     .         rphi(nrmax,0:lmaxd,nsemx), 
     .         rco(nzetmx,0:lmaxd,nsemx),
     .         lambda(nzetmx, 0:lmaxd,nsemx), Zval,vePAO(nrmax),
     .         ePAO(0:lmaxd,nsemx)

               integer
     .           nrval, lmxo,  is, nzeta(0:lmaxd,nsemx),
     .           norb, nsemic(0:lmaxd)


               integer
     .           l,nprin, nnodes, nodd, nrc, i, ir, indx,
     .           izeta, nmax, nmin, nn, nr, nrcomp, nsm

               real(dp)
     .           eigen(0:lmaxd), rc,
     .           dnrm, phi, eshift,
     .           g(nrmax), r, el, ekin, 
     .           r1, r2, dfdi, d2fdi2, d2fdr2, dr,
     .           epot, epot2, rh, dy, eorb, eps


C****NUMBER OF POINTS USED BY POLINT FOR THE INTERPOLATION***
C 
               integer  npoint 
               parameter(npoint=4)
C
C***READING THE ENERGY-SHIFT TO DEFINE THE CUT-OFF RADIUS OF ORBITALS***

           eshift=fdf_physical('PAO.EnergyShift',eshift_default,'Ry')


             norb=0 
             indx=0
             do l=0,lmxo


              do nsm=1,nsemic(l)+1
                if(nzeta(l,nsm).gt.0) then
                  write(6,'(/A,I2)')
     .             'NONODES: Orbitals with angular momentum L=',l
                  goto 50
                endif
              enddo

50            continue

              do nsm=1,nsemic(l)+1

                if(nzeta(l,nsm).gt.0) then

                    write(6,'(/A,I1,A)')
     .                'NONODES: Basis orbitals for state ',
     .                cnfigtb(l,nsm,is), sym(l)


                  if(rco(1,l,nsm).lt.1.0d-5) then    
C**Automatic determination of the cut off radius for the PAOs***
CAtomic eigenvalues
C
                      nnodes=nsm
                      nprin=l+nsm
                      call schro_eq(Zval,rofi,vps(1,l),ve,s,drdi,
     .                  nrval,l,a,b,nnodes,nprin,
     .                  eigen(l),rphi(1,l,nsm))
C
C*Rc given by eshift****   
C                
                       if(eigen(l).gt.0.0d0) then
                          write(6,'(/A,I2,A)')
     .       'NONODES: ERROR Orbital with angular momentum L=',l,
     .       ' not bound in the atom'
                         write(6,'(A)')
     .       'NONODES: ERROR a cut off radius must be explicitely given'
                         call die
                       endif 
 
                       if(abs(eshift).gt.1.0d-5) then
                          el=eigen(l)+eshift
                          call rc_vs_e(a,b,rofi,vps(1,l),
     .                       ve,nrval,l,el,nnodes,rco(1,l,nsm))
                       else
                          rco(1,l,nsm)=rofi(nrval-2)
                       endif 

                  write(6,'(/,A,/,A,f10.6,A)')
     .         'NONODES: PAO cut-off radius determinated from an',
     .         'NONODES: energy shift=',eshift,' Ry'

                 endif  
C
C***


              do izeta=1, nzeta(l,nsm)

CIF THE COMPRESSION FACTOR IS NEGATIVE OR ZERO THE ORBITALS ARE***
CUNTOUCHED**
            if(lambda(izeta,l,nsm).le.0.0d0) lambda(izeta,l,nsm)=1.0d0 
C* 
            if(abs(rco(izeta,l,nsm)).lt.1.0d-5) 
     .                           rco(izeta,l,nsm)=rco(1,l,nsm)
            do i=1,izeta-1
             if((abs(rco(izeta,l,nsm)-rco(i,l,nsm)).lt.1.0d-5).and.
     . (abs(lambda(izeta,l,nsm)-lambda(i,l,nsm)).lt.1.0d-5)) then
                 write(6,'(/,A,I2,A,I2,A,I2,2A)')
     .  'NONODES: WARNING: PAO base function with zeta=',izeta,
     .  ' and zeta=',i,' are identical for ',cnfigtb(l,nsm,is),
     .       sym(l),' state'
               call die
             endif
            enddo

                  rc=rco(izeta,l,nsm)/lambda(izeta,l,nsm)
                  nrc=nint(log(rc/b+1.0d0)/a)+1
                  if (restricted_grid) then
                     nodd=mod(nrc,2)
                     if(nodd.eq.0) then
                        nrc=nrc+1
                     endif
                  endif
                  rc=b*(exp(a*(nrc-1))-1.0d0)
                  rco(izeta,l,nsm)=rc*lambda(izeta,l,nsm)
C****Generate PAO orbitals with increasing number of nodes ***
C*for the different shells*
C 
                     nnodes=nsm
                     nprin=l+nsm
                   call schro_eq(Zval,rofi,vps(1,l),vePAO,s,drdi,
     .              nrc,l,a,b,nnodes,nprin,
     .               eorb,g)  


                       dnrm=0.0d0
                       do ir=2,nrc
                          phi=g(ir) 
                          if(izeta.eq.1) rphi(ir,l,nsm)=phi
                          dnrm=dnrm+drdi(ir)*phi*phi
                          g(ir)=g(ir)/(rofi(ir)**(l+1))
                       enddo 
                       g(1)=g(2)         

C****Normalization of basis functions***
            eps=1.0d-4
            if(abs(dnrm-1.0d0).gt.eps) then
               do ir=1,nrc
                 g(ir)=g(ir)/sqrt(dnrm)
                 if(izeta.eq.1) then
                    rphi(ir,l,nsm)=rphi(ir,l,nsm)/sqrt(dnrm)
                 endif
               enddo
            endif
C****


C*Calculation of the mean value of kinetic and potential energy**
C    Potential and kinetic energy of the orbital before compression



             ekin=0.0d0
             do ir=2,nrc-1
                r=rofi(ir)
                r1=rofi(ir-1)
                r2=rofi(ir+1)
                d2fdi2=(g(ir-1)*r1**(l+1)+g(ir+1)*r2**(l+1)
     .                       -2.0d0*g(ir)*r**(l+1))
                dfdi=0.5d0*(g(ir+1)*r2**(l+1)-g(ir-1)*r1**(l+1))
                dr=drdi(ir)
                d2fdr2= ((-a)*dfdi +  d2fdi2)/dr**2
                ekin=ekin+
     .              dr*g(ir)*r**(l+1)*(-d2fdr2)
     .             +dr*l*(l+1)*(g(ir)*r**l)**2
             enddo


C Kinetic energy after compression

             ekin=ekin/(lambda(izeta,l,nsm)**2)

C Potential energy after compression

             nrcomp=nint(log(rco(izeta,l,nsm)/b+1.0d0)/a)+1
             epot=0.0d0
             epot2=0.0d0
             do ir=1,nrcomp
                r=rofi(ir)
                r2=r/lambda(izeta,l,nsm)
                nr=nint(log(r2/b+1.0d0)/a)+1
                nmin=max(1,nr-npoint)
                nmax=min(nrc,nr+npoint)
                nn=nmax-nmin+1
                call polint(rofi(nmin),g(nmin),nn,r2,rh,dy)
                rh=rh/sqrt(lambda(izeta,l,nsm)**(2*l+3))
                epot=epot+
     .          drdi(ir)*(ve(ir)+vps(ir,l))*(rh*r**(l+1))**2
                epot2=epot2+
     .          drdi(ir)*vps(ir,l)*(rh*r**(l+1))**2
             enddo
             eorb=ekin+epot




             write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')
     .          'izeta =',izeta,
     .          'lambda =',lambda(izeta,l,nsm),
     .          'rc =',rco(izeta,l,nsm),
     .          'energy =',eorb  

                if(izeta.eq.1) ePAO(l,nsm)=eorb

          write(6,'(a25,f12.6)') 'kinetic =',ekin
          write(6,'(a25,f12.6)') 'potential(screened) =',epot
          write(6,'(a25,f12.6)') 'potential(ionic) =',epot2 

            norb=norb+(2*l+1)
            indx=indx+1
            call comBasis(is,a,b,rofi,g,l,
     .              rco(izeta,l,nsm),lambda(izeta,l,nsm),izeta,
     .              nsm,nrc,indx)



              enddo 
               

            call compress_PAO(a,b,rofi,rphi(1,l,nsm),
     .              rco(1,l,nsm),lambda(1,l,nsm))



              endif    
             enddo 
            enddo 
      
            end subroutine nonodes

!
               subroutine comBasis(is,a,b,rofi,rphi,
     .            l,rc,lambda,nzeta,nsemic,nrc,norb)
C*
C  Generates the common blocks for the storage of the information 
C  about the basis set orbitals.
C
C  Written by D. Sanchez-Portal, Aug. 1998.
C*

               implicit none

               integer l, is, norb, nrc, nzeta, nsemic

               real(dp) rc, lambda, rphi(nrmax), a, b,
     .            rofi(nrmax)  

C****Internal variables
C
             integer itb, nr, nmax, nmin, nn
        
             real(dp) delt, r, phi, dy, yp1, ypn
             
 
****NUMBER OF POINTS USED BY POLINT FOR THE INTERPOLATION***
C
          integer npoint
          parameter(npoint=4)
C
              rcotb(nzeta,l,nsemic,is)=rc
              lambdatb(nzeta,l,nsemic,is)=lambda


CINTERPOLATION TO GENERATE TABLES
C
            delt=rc/(dble(ntbmax-1)+1.0d-20) 
  
            if(delt.gt.deltmax) then
              write(6,'(a)')
     .    'comBasis: WARNING It might be a good idea to increase'
              write(6,'(a)')
     .    'comBasis: WARNING parameter ntbmax (in file atmparams.f) '
              write(6,'(a,i6)')
     .    'comBasis: WARNING to at least ntbmax = ',
     .        nint(Rc/deltmax)+2
            endif
!
!           First two entries are used for other purposes...
!
            table(1,norb,is)=delt
            table(2,norb,is)=dble(l)

            do itb=1,ntbmax-1
                 r=delt*(itb-1)
                 r=r/lambda
                 nr=nint(log(r/b+1.0d0)/a)+1
                 nmin=max(1,nr-npoint)
                 nmax=min(nrc,nr+npoint)
                 nn=nmax-nmin+1
                 call polint(rofi(nmin),rphi(nmin),nn,r,phi,dy)
                 phi=phi/sqrt(lambda**(2*l+3))
                 table(itb+2,norb,is)=phi
            enddo
!           Why do we stop in ntbmax-1 ??
!!          If we use the next statement, the last few elements of
!!          tab2 will be *very* wrong!!! (use dump_atom)
!!          In other places of the program, it is set to zero, as in (2)
!!!!!            table(ntbmax+2,norb,is) = table(ntbmax+1,norb,is)  !AG
            table(ntbmax+2,norb,is) = 0.d0
!!          Either way, it looks quite arbitrary
!!
C****TABLE WITH THE SECOND DERIVATIVE ***

            yp1=huge(1.d0)
            ypn=huge(1.d0)
            call SPLINE( delt, table(3:2+ntbmax,norb,is), ntbmax,
     &                   yp1, ypn, tab2(1:ntbmax,norb,is) )

C
      end subroutine combasis
!
      subroutine SPLITGAUSS(Zval,is,a,b,rofi,drdi,s,
     .             vps,ve,vePAO,
     .             nrval,lmxo,nsemic,
     .             nzeta,rco,lambda, rphi, ePAO, norb) 

           use basis_specs, only: restricted_grid

C****
C Calculates the atomic orbitals basis set, using the option SPLITGAUSS
C for the generation of the augmentation orbitals.
C  Written by D. Sanchez-Portal, Aug. 1998
C****

               implicit none

               real(dp)
     .         a, b, rofi(nrmax), vps(nrmax,0:lmaxd),
     .         drdi(nrmax), s(nrmax), ve(nrmax),
     .         rphi(nrmax,0:lmaxd,nsemx), rco(nzetmx,0:lmaxd,nsemx),
     .         lambda(nzetmx, 0:lmaxd,nsemx), Zval,vePAO(nrmax),
     .         ePAO(0:lmaxd,nsemx)

               integer
     .           nrval, lmxo, is, nzeta(0:lmaxd,nsemx),
     .           norb, nsemic(0:lmaxd)


               integer
     .           l,nprin, nnodes, nodd, nrc, i, ir, indx,
     .           izeta, nmax, nmin, nn, nr, nrcomp, nsm

               real(dp)
     .           eigen(0:lmaxd), rc,
     .           dnrm, phi, 
     .           cons, fac, eshift, pi, gexp,
     .           g(nrmax), r, el, ekin, 
     .           r1, r2, dfdi, d2fdi2, d2fdr2, dr,
     .           epot, epot2, rh, dy, eorb, eps, dlapl


C****NUMBER OF POINTS USED BY POLINT FOR THE INTERPOLATION***
C 
               integer  npoint 
               parameter(npoint=4)

C***READING THE ENERGY-SHIFT TO DEFINE THE CUT-OFF RADIUS OF ORBITALS***

           eshift=fdf_physical('PAO.EnergyShift',eshift_default,'Ry')

                   pi=acos(-1.0d0) 
 
             norb=0 
             indx=0
             do l=0,lmxo

              do nsm=1,nsemic(l)+1
                if(nzeta(l,nsm).gt.0) then
                    write(6,'(/A,I2)')
     .               'SPLITGAUSS: Orbitals with angular momentum L=',l

                  goto 50

                endif
              enddo

50            continue

              do nsm=1,nsemic(l)+1

                if(nzeta(l,nsm).gt.0) then

                    write(6,'(/A,I1,A)')
     .               'SPLITGAUSS: Basis orbitals for state ',
     .                cnfigtb(l,nsm,is), sym(l)

                  if(rco(1,l,nsm).lt.1.0d-5) then    
C**Automatic determination of the cut off radius for the PAOs***
CAtomic eigenvalues
C
                      nnodes=nsm
                      nprin=l+nsm
                      call schro_eq(Zval,rofi,vps(1,l),ve,s,drdi,
     .                  nrval,l,a,b,nnodes,nprin,
     .                  eigen(l),rphi(1,l,nsm))
C
C*Rc given by eshift****   
C                
                       if(eigen(l).gt.0.0d0) then
                          write(6,'(/A,I2,A)')
     .  'SPLITGAUSS: ERROR Orbital with angular momentum L=',l,
     .       ' not bound in the atom'
                         write(6,'(A)')
     .  'SPLITGAUSS: ERROR a cut off radius must be explicitely given' 
          call die
                       endif 
 
                       if(abs(eshift).gt.1.0d-5) then
                          el=eigen(l)+eshift
                          call rc_vs_e(a,b,rofi,vps(1,l),
     .                           ve,nrval,l,el,nnodes,rco(1,l,nsm))
                       else
                          rco(1,l,nsm)=rofi(nrval-2)
                       endif 

                  write(6,'(/,A,/,A,f10.6,A)')
     .   'SPLITGAUSS: PAO cut-off radius determinated from an',
     .   'SPLITGAUSS: energy shift=',eshift,' Ry'

                 endif  
C
C***

CIF THE COMPRESSION FACTOR IS NEGATIVE OR ZERO THE ORBITALS ARE***
CUNTOUCHED**
               if(lambda(1,l,nsm).le.0.0d0) lambda(1,l,nsm)=1.0d0
C*

              do izeta=1, nzeta(l,nsm)

             if(abs(rco(izeta,l,nsm)).lt.1.0d-5) 
     .                      rco(izeta,l,nsm)=rco(1,l,nsm)
CWith spligauss option, compression factor must be taken****
Cas the gaussian exponent
              if(izeta.gt.1) then 
                  if(lambda(izeta,l,nsm).le.0.0d0) then 
                    write(6,'(/a,/a,a)')
     .'SPLITGAUSS: ERROR: with SPLITGAUSS option the compression ',
     .'SPLITGAUSS: ERROR: factors for all the augmentation functions',
     .   ' must be explicitely specified' 
                    call die
                  endif
                  gexp=abs(lambda(izeta,l,nsm))
                  gexp=1.0d0/(gexp**2)
                  lambda(izeta,l,nsm)=1.0d0
              endif
C****
                  
                  rc=rco(izeta,l,nsm)/lambda(izeta,l,nsm)
                  nrc=nint(log(rc/b+1.0d0)/a)+1
                  if (restricted_grid) then
                     nodd=mod(nrc,2)
                     if(nodd.eq.0) then
                        nrc=nrc+1
                     endif
                  endif
                  rc=b*(exp(a*(nrc-1))-1.0d0)
                  rco(izeta,l,nsm)=rc*lambda(izeta,l,nsm)  


                  if(izeta.eq.1) then 
C****Generate a PAO orbital for the first shell of basis functions****
C 
                      nnodes=nsm
                      nprin=l+nsm
                      call schro_eq(Zval,rofi,vps(1,l),vePAO,s,drdi,
     .                  nrc,l,a,b,nnodes,nprin,
     .                  eorb,rphi(1,l,nsm)) 
                       dnrm=0.0d0
                       do ir=2,nrc
                          phi=rphi(ir,l,nsm) 
                          dnrm=dnrm+drdi(ir)*phi*phi
                          g(ir)=phi/(rofi(ir)**(l+1))
                       enddo 
                       g(1)=g(2)         
c                   elseif(izeta.gt.1) then 
                    else
                     fac=1.0d0
                     do i=0,l
                       fac=(2*i+1)*fac
                     enddo

                     cons=sqrt(pi)*fac/(2.0d0**(l+2))
                     cons=cons/((2.0d0*gexp)**(l+1.5d0))
                     cons=1.0d0/sqrt(cons)

                     dnrm=0.0d0
                     do ir=1,nrc
                       r=rofi(ir)
                       phi=cons*exp((-gexp)*r**2)
                       dnrm=dnrm+drdi(ir)*(phi*r**(l+1))**2
                       g(ir)=phi
                     enddo
                    endif  



C****Normalization of basis functions***
            eps=1.0d-4
            if(abs(dnrm-1.0d0).gt.eps) then
               do ir=1,nrc
                 g(ir)=g(ir)/sqrt(dnrm)
                 if(izeta.eq.1) then
                    rphi(ir,l,nsm)=rphi(ir,l,nsm)/sqrt(dnrm)
                 endif
               enddo
            endif
C****


C*Calculation of the mean value of kinetic and potential energy**
C    Potential and kinetic energy of the orbital before compression


           if(izeta.eq.1) then 

             ekin=0.0d0
             do ir=2,nrc-1
                r=rofi(ir)
                r1=rofi(ir-1)
                r2=rofi(ir+1)
                d2fdi2=(g(ir-1)*r1**(l+1)+g(ir+1)*r2**(l+1)
     .                       -2.0d0*g(ir)*r**(l+1))
                dfdi=0.5d0*(g(ir+1)*r2**(l+1)-g(ir-1)*r1**(l+1))
                dr=drdi(ir)
                d2fdr2= ((-a)*dfdi +  d2fdi2)/dr**2
                ekin=ekin+
     .              dr*g(ir)*r**(l+1)*(-d2fdr2)
     .             +dr*l*(l+1)*(g(ir)*r**l)**2
             enddo


C Kinetic energy after compression

             ekin=ekin/(lambda(izeta,l,nsm)**2)

C Potential energy after compression

             nrcomp=nint(log(rco(izeta,l,nsm)/b+1.0d0)/a)+1
             epot=0.0d0
             epot2=0.0d0
             do ir=1,nrcomp
                r=rofi(ir)
                r2=r/lambda(izeta,l,nsm)
                nr=nint(log(r2/b+1.0d0)/a)+1
                nmin=max(1,nr-npoint)
                nmax=min(nrc,nr+npoint)
                nn=nmax-nmin+1
                call polint(rofi(nmin),g(nmin),nn,r2,rh,dy)
                rh=rh/sqrt(lambda(izeta,l,nsm)**(2*l+3))
                epot=epot+
     .          drdi(ir)*(ve(ir)+vps(ir,l))*(rh*r**(l+1))**2
                epot2=epot2+
     .          drdi(ir)*vps(ir,l)*(rh*r**(l+1))**2
             enddo
             eorb=ekin+epot
        
            elseif(izeta.gt.1) then  

             epot=0.0d0
             epot2=0.0d0
             ekin=0.0d0
             do ir=2,nrc
               r=rofi(ir)
               phi=g(ir)*r**l
               epot=epot+
     .         drdi(ir)*(ve(ir)+vps(ir,l))*(phi*r)**2
               epot2=epot2+
     .         drdi(ir)*vps(ir,l)*(phi*r)**2
               dlapl=
     .        -((l*(l+1)-2*gexp*(2*l+3)*r**2+4*(gexp*r**2)**2)*phi)
               dlapl=(dlapl+l*(l+1)*phi)
               ekin=ekin +
     .         drdi(ir)*dlapl*phi
             enddo
             eorb=ekin+epot

            endif 

            if(izeta.eq.1) then  

             write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')
     .          'izeta =',izeta,
     .          'lambda =',lambda(izeta,l,nsm),
     .          'rc =',rco(izeta,l,nsm),
     .          'energy =',eorb 

                 ePAO(l,nsm)=eorb

            elseif(izeta.gt.1) then  

                write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')
     .            'izeta=',izeta,'gaussian exponent=',gexp,
     .            'rc=',rco(izeta,l,nsm),'energy=',eorb

            endif 

          write(6,'(a25,f12.6)') 'kinetic =',ekin
          write(6,'(a25,f12.6)') 'potential(screened) =',epot
          write(6,'(a25,f12.6)') 'potential(ionic) =',epot2 

            norb=norb+(2*l+1)
            indx=indx+1
            call comBasis(is,a,b,rofi,g,l,
     .              rco(izeta,l,nsm),lambda(izeta,l,nsm),izeta,
     .              nsm,nrc,indx)

              enddo 
               
            call compress_PAO(a,b,rofi,rphi(1,l,nsm),
     .              rco(1,l,nsm),lambda(1,l,nsm))


              endif    
              
              enddo
             enddo 

            end subroutine splitgauss
!
               subroutine compress_PAO(a,b,rofi,rphi,
     .            rc,lambda)
C***
C   Compression of a PAO orbital according to the compression
C   factor lambda. Input and outputare stored in the same array, 
C   and the radial grid is identical in input and output. 
C  Written by D. Sanchez-Portal, Aug. 1998
C***
               implicit none

               real(dp) rc, lambda, rphi(nrmax), a, b,
     .            rofi(nrmax)  

C****Internal variables
C
             integer nr, nmax, nmin, nn, maxpoint, nrc, ir  
        
             real(dp)  r, phi, dy, rmax
             
             real(dp)
     .          aux(nrmax)

****NUMBER OF POINTS USED BY POLINT FOR THE INTERPOLATION***
C
          integer npoint
          parameter(npoint=4)
C

C*INTERPOLATION TO CALCULATE THE VALUE OF THE FIRST-SHELL PAO ****
C***BASIS FUNCTIONS IN THE LOGARITHMIC MESH AFTER COMPRESSION ****
C
              nrc=nint(log(rc/b+1.0d0)/a)+1
              rmax=rc/lambda
              maxpoint=nint(log(rmax/b+1.0d0)/a)+1
              
              do ir=2,nrc
                r=rofi(ir)/lambda
                nr=nint(log(r/b+1.0d0)/a)+1
                nmin=max(1,nr-npoint)
                nmax=min(maxpoint,nr+npoint)
                nn=nmax-nmin+1
                call polint(rofi(nmin),rphi(nmin),nn,r,phi,dy) 
                aux(ir)=phi/sqrt(lambda)
              enddo 
              rphi(1)=0.0d0 
              do ir=2,nrc
               rphi(ir)=aux(ir)
              enddo 
              do ir=nrc+1,nrmax
                rphi(ir)=0.0d0 
              enddo 

         end subroutine compress_PAO
!
       subroutine atm_pop(is,iz,q,qPAO,lmxo,
     .       nzeta,semic,nsemic,polorb,basp) 
C
C Returns the ground states atomic population for each species.
C This information is required for the screening of the local 
C pseudopotential.
!! Even though qtb is passed, atom itself only uses qPAO...
C Written by D. Sanchez-Portal, Aug. 1998
C

         implicit none

         real(dp), intent(out) ::   q(1:), qPAO(0:lmaxd,nsemx) 

         integer, intent(in)   ::
     .     nzeta(0:lmaxd,nsemx),polorb(0:lmaxd,nsemx), 
     .     nsemic(0:lmaxd),iz,
     .     lmxo ,is 

         logical, intent(in)   ::    semic
         type(basis_def_t), pointer   :: basp


C***Internal variables*

        real(dp)  qatm(0:3)
          
        integer noPAO, l, izeta, m, norb, noPol, lpop,
     .     nsm, nvalence, config(0:lmaxd), i, j
        character*70  line

        qatm(0:3) = basp%ground_state%occupation(0:3)
        config(:) = 0
        config(0:3) = basp%ground_state%n(0:3)

        qPAO(0:lmaxd,1:nsemx) = 0.d0  ! AG
!
!       What is this? semic says whether nsemic is defined....
!
        do l=0,lmxo 
          nvalence=nsemic(l)+1
     $        -(cnfigtb(l,nsemic(l)+1,is)-config(l))
!!          print *, "l, config(l), nvalence:", l, config(l), nvalence
          nsm=nvalence
C Add a check on whether nsemx is large enough
          if (nsm.gt.nsemx) then
            call die('need to increase nsemx')
          endif
          qPAO(l,nsm)=0.0d0
          if(l.le.3) qPAO(l,nsm)=qatm(l) 
        enddo 

        if(semic) then 
        do l=0,lmxo
           do nsm=1,nsemic(l)+1
           nvalence=nsemic(l)+1
     .        -(cnfigtb(l,nsemic(l)+1,is)-config(l))
              if(nsm.lt.nvalence) qPAO(l,nsm)=2*(2*l+1)
              if(nsm.gt.nvalence) qPAO(l,nsm)=0.0d0
           enddo 
        enddo 
        endif 

        noPAO=0
        do l=0,lmxo 
          do nsm=1,nsemic(l)+1
           if(nzeta(l,nsm).gt.0) then 
            do  izeta=1,nzeta(l,nsm) 
                  do m=1,2*l+1 
                   q(noPAO+m)=0.0d0
                   if(izeta.eq.1) 
     .                 q(noPAO+m)=qPAO(l,nsm)/(2*l+1)
                  enddo  
                  noPAO=noPAO+2*l+1
            enddo
           endif 
          enddo 
        enddo 

        noPol=0
        do l=0,lmxo
          do nsm=1,nsemic(l)+1
            if(polorb(l,nsm).gt.0) then 
             do  izeta=1,polorb(l,nsm)
                do m=1,2*(l+1)+1 
                  q(noPAO+noPol+m)=0.0d0
                enddo  
                noPol=noPol+2*(l+1)+1
             enddo    
            endif 
          enddo
        enddo      
        norb=noPAO+noPol 

        lpop=min(3,lmxo)  
          write(6,'(/,2a)') 'atm_pop: Valence configuration',
     .                      ' (for local Pseudopot. screening):' 
          do l=0,lpop 
            write(line,'(7(1x,i1,a1,a1,f5.2,a1))')
     .          (cnfigtb(l,nsm,is),sym(l),"(",qPAO(l,nsm),")",
     .       nsm=1,nsemic(l)+1-(cnfigtb(l,nsemic(l)+1,is)-config(l)))
            write(6,'(a)') line
 
c           write(6,*) 
c    .         (cnfigtb(l,nsm,is),sym(l),'(',qPAO(l,nsm),')',
c    .                                      nsm=1,nsemic(l)+1) 
          enddo

        end subroutine atm_pop
!

        subroutine Vna(is,Zval,qPAO,rphi,rco,nsemic,vloc,
     .        a,b,rofi,drdi,nrval,lmxo,nVna) 
C*
C  Generates the neutral-atom pseudopotential.
C  D. Sanchez-Portal, Aug. 1998.
C  Modify by DSP, July 1999
C*
          
          implicit none
 
          real(dp)
     .    Zval, qPAO(0:lmaxd,nsemx), rofi(nrmax), 
     .    rphi(nrmax,0:lmaxd,nsemx),
     .    drdi(nrmax),rco(nzetmx,0:lmaxd,nsemx),a,b, vloc(nrmax) 
      
          integer
     .    nrval, lmxo,is, nVna, nsemic(0:lmaxd),nsm


C***Internal variables****

         real(dp)
     .    rho(nrmax), chval, ve(nrmax), s(nrmax),eps, phi,
     .    rcocc, dincv, rVna,kmax,tail,filterfactor,grid,minnorm
       
         logical filterVna
  
         integer
     .    nrc,ir, l, ncocc

          do ir=1,nrval
             rho(ir)=0.0d0
             s(ir)=sqrt(drdi(ir))
          enddo

          chval=0.0d0 
          ncocc=0 
          do l=0,lmxo
            do nsm=1,nsemic(l)+1
             if(qPAO(l,nsm).gt.0.0d0) then
              nrc=nint(log(rco(1,l,nsm)/b+1.0d0)/a)+1 
              ncocc=max(ncocc,nrc)
              do ir=2,nrc 
                phi=rphi(ir,l,nsm)
                rho(ir)=rho(ir)+qPAO(l,nsm)*phi**2 
                chval=chval+drdi(ir)*qPAO(l,nsm)*phi**2
              enddo
            endif
          enddo 
         enddo
         rho(1)=0.0d0 
!
!        Now chval contains the total charge
!        in the occupied basis states.
!
!        It should be equal to the (nominal) 
!        valence charge (i.e. Znuc-Zcore).
!
         write(6,'(a,2f10.5)') 'Vna: chval, zval: ', chval, zval
!
         eps=1.0d-4
         if(abs(chval-zval).gt.eps) then
           do ir=2,nrval
              rho(ir)=zval*rho(ir)/chval
           enddo
         endif

C**CALCULATION OF THE HARTREE POTENTIAL DUE TO THE NEW VALENCE CHARGE**
C
          call vhrtre(rho,ve,rofi,drdi,s,nrval,a)
C
C****LOCAL NEUTRAL-ATOM PSEUDOPOTENTIAL*
          eps=1.0d-5
          nVna=0
          do ir=nrval,2,-1
               dincv=vloc(ir)+ve(ir)
               if((abs(dincv).gt.eps).and.(nVna.eq.0)) nVna=ir+1
               ve(ir)=dincv
          enddo 
          nVna=max(nVna,ncocc)

C****CUT-OFF RADIUS FOR THE LOCAL NEUTRAL-ATOM PSEUDOPOTENTIAL 

          rcocc=b*(exp(a*(ncocc-1))-1.0d0)

          if(nVna.eq.ncocc) then
            rVna=rcocc
          else
            rVna=b*(exp(a*(nVna-1))-1.0d0)
          endif 

          write(6,'(/,a,f10.6)')
     .  'Vna:  Cut-off radius for the neutral-atom potential: ', 
     .  rVna

         if(rVna.gt.(rcocc+0.5d0)) then
           write(6,"(2a,f12.5)")'Vna: WARNING: ',
     .  'Cut-off radius for the neutral-atom potential, rVna =', 
     .      rVna
           write(6,"(2a,f12.5)")'Vna: WARNING: ',
     .        'Cut-off radius for charge density =', rcocc
           write(6,"(2a)")'Vna: WARNING: ',
     .        'Check atom: Look for the sentence:'
           write(6,"(2a)")'Vna: WARNING: ',
     .        'LOCAL NEUTRAL-ATOM PSEUDOPOTENTIAL'
           write(6,"(2a)")'Vna: WARNING: ',
     .        'Increasing the tolerance parameter EPS'
           write(6,"(2a)")'Vna: WARNING: ',
     .        'might be a good idea'
         endif

          ve(1)= ( ve(2)*rofi(3)**2 - ve(3)*rofi(2)**2 ) /
     .          (      rofi(3)**2 -      rofi(2)**2 )
         
!Filter
               !Filter the high-k components of Vna
      filterVna = fdf_boolean("Vna.Filter",.false.)

      if (filterVna)then
         kmax = fdf_double("Vna.FilterCutoff",-10.0_dp)
         minnorm = fdf_double("VNA.minnorm",0.999_dp)
         if (kmax .eq. -10.0)then
            filterFactor = fdf_double("Vna.FilterFactor",1.0_dp)
            grid = fdf_physical('MeshCutoff',100.0_dp,'Ry')
            kmax = filterFactor*sqrt(grid)
         endif

         write(6,"(A,f8.3)")"na: Filtered Vna cutoff (Bohr^-1):",
     $        kmax       
         
         call filter(0,nVna,rofi,Ve,kmax,0,minnorm)      

         tail = 0.0_dp
         
      endif
      

C**Construct the common block with the neutral-atom potential**** 
C
          call comVna(is,a,b,rofi,Ve,nVna,1.0d0)
C
          end subroutine vna
!
!
           subroutine comVna(is,a,b,rofi,Vna,nVna,flting)

C  Creates the common block with the information about the neutral atom
C  pseudoptential.
C   D. Sanchez-Portal, Aug. 1998.

           implicit none

           integer nVna, is
           real(dp) 
     .        Vna(nrmax), rofi(nrmax), a ,b, flting

           integer nr, nmin, nmax, nn, itb
           real(dp)  yp1, ypn, dy, v, rVna, delt, r

C****NUMBER OF POINTS USED BY POLINT FOR THE INTERPOLATION***
C
          integer npoint
          parameter(npoint=4)

       if(flting.lt.0.0d0) then    !! Floating orbital
          table(1,0,is)=0.0d0
          table(2,0,is)=0.0d0 
          do itb=1,ntbmax-1
             tab2(itb,0,is)=0.0d0
             table(itb+2,0,is)=0.0d0
          enddo
          table(ntbmax+2,0,is)=0.0d0
          return                   !! Return
       endif

       rVna=b*(exp(a*(nVna-1))-1.0d0)
       delt=rVna/(dble(ntbmax-1)+1.0d-20) 
 
       if(delt.gt.deltmax) then
              write(6,'(a)')
     .    'comVna: WARNING It might be a good idea to increase'
              write(6,'(a)')
     .    'comVna: WARNING parameter ntbmax (in file atmparams.f) '
              write(6,'(a,i6)')
     .    'comVna: WARNING to at least ntbmax = ',
     .        nint(rVna/deltmax)+2
       endif

       table(1,0,is)=delt
       table(2,0,is)=rVna

       do itb=1,ntbmax-1
          r=delt*(itb-1)
          nr=nint(log(r/b+1.0d0)/a)+1
          nmin=max(1,nr-npoint)
          nmax=min(nVna,nr+npoint)
          nn=nmax-nmin+1 
          call polint(rofi(nmin),Vna(nmin),nn,r,v,dy) 
          table(itb+2,0,is)=v 
       enddo 
       table(ntbmax+2,0,is)=0.0d0

C****TABLE WITH THE SECOND DERIVATIVE ***
C 
       yp1=0.d0
       ypn=huge(1.d0)

       call SPLINE( delt, table(3:2+ntbmax,0,is), ntbmax,
     &              yp1, ypn, tab2(1:ntbmax,0,is))
       end  subroutine comVna
!
       subroutine slfe_local(slfe,vlocal,rofi,a,nVna,drdi)

C Calculates the self-energy associated to the local-pseudopotential
C charge density.
C Written by D. Sanchez-Portal, Aug. 1998.

         implicit none
           
         real(dp) slfe, vlocal(nrmax),rofi(nrmax),a,
     .       drdi(nrmax) 

         integer nVna 
           
C*Internal variables**
          
         real(dp) slf, a2b4, s(nrmax), g0, g1, g2,
     .       g3, g4, d2g, d2u

         integer ir

          do ir=1,nVna
             s(ir)=sqrt(drdi(ir))
          enddo

           a2b4=0.25d0*a*a
           slf=0.0d0
           do ir=2,nVna-1
              if((ir.gt.2).and.(ir.lt.(nVna-1))) then
                g0=vlocal(ir-2)*rofi(ir-2)/s(ir-2)
                g1=vlocal(ir-1)*rofi(ir-1)/s(ir-1)
                g2=vlocal(ir)*rofi(ir)/s(ir)
                g3=vlocal(ir+1)*rofi(ir+1)/s(ir+1)
                g4=vlocal(ir+2)*rofi(ir+2)/s(ir+2)

                d2g=(16.0d0*(g1+g3)-(g0+g4)-30.0d0*g2)/12.0d0

              else
                g1=vlocal(ir-1)*rofi(ir-1)/s(ir-1)
                g2=vlocal(ir)*rofi(ir)/s(ir)
                g3=vlocal(ir+1)*rofi(ir+1)/s(ir+1)

                d2g=g1+g3-2.0d0*g2

              endif

              d2u=d2g-a2b4*g2

              slf=slf - g2*d2u*0.25d0

           enddo

           slfe=slf
   
           end subroutine slfe_local
!
            subroutine POLgen(is,iz,a,b, rofi, drdi,
     .          ePAO,rphi,rco,vps,ve,
     .          polorb,lmxo, nsemic,norb,
     $          rinn,vcte,nrval, split_norm,atm_label)
C****
C Calculates the polarization  orbitals for the basis set augmentation.
C Written by D. Sanchez-Portal, Aug. 1998.
C Modify by DSP, July 1999
C****
C
               implicit none

               real(dp)
     .         a, b, rofi(nrmax), vps(nrmax,0:lmaxd),
     .         ve(nrmax), drdi(nrmax),
     .         rphi(nrmax,0:lmaxd,nsemx), 
     .         rco(nzetmx,0:lmaxd,nsemx),
     .         ePAO(0:lmaxd,nsemx),
     .         rinn(0:lmaxd,nsemx), vcte(0:lmaxd,nsemx),
     $         split_norm(0:lmaxd,nsemx)

               character(len=*) atm_label

               integer
     .           lmxo, is, iz, norb, polorb(0:lmaxd,nsemx), 
     .           nsemic(0:lmaxd), nrval

               integer
     .           l, nrc, nsp, ir,indx,
     .           ipol, nsm, nrcomp

               real(dp)
     .           rc, rcpol(nzetmx,0:lmaxd,nsemx),
     .           phipol(nrmax), split_table(nrmax),
     .           rnrm(nrmax), dnrm, phi,
     .           cons1, cons2, spln, 
     .           g(nrmax), r, ekin, 
     .           r1, r2, dfdi, d2fdi2, d2fdr2, dr,
     .           epot, epot2, eorb, eps,
     $           rcsan, exponent, vsoft(nrmax), vePAOsoft(nrmax)

      logical::filterorbitals, new_split_code, split_tail_norm
      logical :: fix_split_table
      real(dp)::kmax,grid,filterFactor,minnorm, spln_min
      
           split_tail_norm=fdf_boolean('PAO.SplitTailNorm',.false.)
           fix_split_table=fdf_boolean('PAO.FixSplitTable',.false.)
           if (split_tail_norm .or. fix_split_table) then
              new_split_code = .true.
           else
              new_split_code=fdf_boolean('PAO.NewSplitCode',.false.)
           endif

             norb=0 
             indx=0
             
             do l=0,lmxo
               do nsm=1,nsemic(l)+1
                if(polorb(l,nsm).gt.0) then
                    write(6,'(/A,I2)')
     .    'POLgen: Perturbative polarization orbital with L= ',l+1
                  goto 50
                endif
               enddo

50            continue


              do nsm=1,nsemic(l)+1

               if (polorb(l,nsm).gt.0) then 

C Calculate the soft-confinement potential for the polarization orbitals

              call build_vsoft(is,l,nsm,rinn(l,nsm),vcte(l,nsm),
     $                       a, b, rco(1,l,nsm),rofi,nrval,vsoft)

              do ir = 1, nrval
                vePAOsoft(ir) = ve(ir) + vsoft(ir)
              enddo

               do ipol=1,polorb(l,nsm)

                   write(6,'(/A,I1,A)')
     .              'POLgen: Polarization orbital for state ',
     .               cnfigtb(l,nsm,is), sym(l)
                  
                if (ipol.eq.1) then  
                  rc=rco(1,l,nsm) 
                  rcpol(ipol,l,nsm)=rc
                  nrc=nint(log(rc/b+1.0d0)/a)+1

!                 Generate the polarization function perturbatively 
!                 from the original PAO**

                  call polarization(a,rofi,rphi(1,l,nsm),vps(1,l),
     .                 vePAOsoft,drdi,nrc,l,ePAO(l,nsm),g,nrc)

                       dnrm=0.0d0
                       do ir=2,nrc-1
                          phi=g(ir)  
                          phipol(ir)=phi
                          dnrm=dnrm+drdi(ir)*phi*phi
                          rnrm(ir)=dnrm 
                          g(ir)=g(ir)/(rofi(ir)**(l+2))
                       enddo   
                       g(1)=g(2)
                       g(nrc)=0.0d0
                       phipol(nrc)=0.0d0

                   if (split_tail_norm) then
                      write(*,"(a)") "Split based on tail norm"
                      ! Follows the criterion of the JPC paper, but
                      ! with more contrast 
                      !(use the actual math norm (sqrt(int(f^2)))

                      split_table(1:nrc) =
     $                     sqrt(max(1.0_dp-rnrm(1:nrc),0.0_dp))

                   else
                   !  Do a full scan of the old method
                   !  (norm of tail+parabola)
                      call split_scan(nrc, rofi, drdi, l+1,
     $                        phipol,rnrm,atm_label,
     $                        split_table,fix_split_table)
                   endif

                else

!                  Multiple shells can be generated using the split scheme

                   rc=rco(1,l,nsm) 
                   nrc=nint(log(rc/b+1.0d0)/a)+1
       
            spln=split_norm(l,nsm)
            if(ipol.gt.2) then
              spln=split_norm(l,nsm)/(2.0d0*(ipol-2) )
            endif

            if (new_split_code) then
               call find_split_location(nrc,rofi,drdi,
     $              phipol,split_table,
     $              l+1,spln,cons1,cons2,nsp)

            else
               spln_min = minval(split_table(1:nrc))
               if (spln < spln_min) then
                  write(6,"(a,f8.5,a,f8.5)")
     $            "WARNING: Minimum split_norm parameter: ",
     $            spln_min, ". Will not be able to generate "
     $            // "orbital with split_norm = ", spln
                  call die("See manual for new split options")
               endif
               call parabola(a,b,nrc,phipol,rnrm,
     .                   l+1,spln,cons1,cons2,nsp)
            endif



C***Cut-off radius for the split orbital with a desired norm* 
         nrc=nsp
         rcpol(ipol,l,nsm)=rofi(nrc)

              
            dnrm=0.0d0
            do ir=2,nrc-1
              r=rofi(ir)
C***parabolic split***
              phi=-(cons1*r**2+cons2)+phipol(ir)/(r**(l+2))
C*

              dnrm=dnrm+drdi(ir)*(phi*r**(l+2))**2
              g(ir)=phi
            enddo
            g(1)=g(2)
            g(nrc)=0.0d0
               
            
              endif

C****Normalization of basis functions***
            eps=1.0d-4
            if(abs(dnrm-1.0d0).gt.eps) then
               do ir=1,nrc
                 g(ir)=g(ir)/sqrt(dnrm)
                 if(ipol.eq.1) then
                    phipol(ir)=phipol(ir)/sqrt(dnrm)
                    rnrm(ir)=rnrm(ir)/dnrm
                 endif
               enddo
            endif
C****

C*Calculation of the mean value of kinetic and potential energy**
C    Potential and kinetic energy of the orbital

             ekin=0.0d0 
             epot=0.0d0
             epot2=0.0d0
             do ir=2,nrc-1
                r=rofi(ir)
                r1=rofi(ir-1)
                r2=rofi(ir+1)
                d2fdi2=(g(ir-1)*r1**(l+2)+g(ir+1)*r2**(l+2)
     .                       -2.0d0*g(ir)*r**(l+2))
                dfdi=0.5d0*(g(ir+1)*r2**(l+2)-g(ir-1)*r1**(l+2))
                dr=drdi(ir)
                d2fdr2= ((-a)*dfdi +  d2fdi2)/dr**2
                ekin=ekin+
     .              dr*g(ir)*r**(l+2)*(-d2fdr2)
     .             +dr*(l+1)*(l+2)*(g(ir)*r**(l+1))**2  

                epot=epot+ 
     .          drdi(ir)*(vePAOsoft(ir)+vps(ir,l))*(g(ir)*r**(l+2))**2 
                epot2=epot2+
     .          drdi(ir)*vps(ir,l)*(g(ir)*r**(l+2))**2

             enddo
             eorb=ekin+epot
       
            if(ipol.eq.1) then  

             write(6,'(/,(3x,a,i2),2(/,a25,f12.6))')
     .          'izeta =',ipol,
     .          'rc =',rcpol(ipol,l,nsm),
     .          'energy =',eorb 

            elseif(ipol.gt.1) then 

            write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')
     .         'izeta =',ipol,
     .         'rmatch =',rcpol(ipol,l,nsm),
     .         'splitnorm =',spln,
     .         'energy =',eorb 

            endif 

          write(6,'(a25,f12.6)') 'kinetic =',ekin
          write(6,'(a25,f12.6)') 'potential(screened) =',epot
          write(6,'(a25,f12.6)') 'potential(ionic) =',epot2 

            norb=norb+(2*(l+1)+1)
            indx=indx+1

!Filter
                  !Filter the high-k components of the orbital
                  filterOrbitals = fdf_boolean("PAO.Filter",.false.)
                  filterOrbitals = .false.

                  kmax = fdf_double("PAO.FilterCutoff",-10.0_dp)
                  if (kmax .eq. -10.0)then
                     filterFactor = fdf_double("PAO.FilterFactor",
     $                    0.7_dp)
                     grid = fdf_physical('MeshCutoff',100.0_dp,'Ry')
                     kmax = filterFactor*sqrt(grid)
                  endif          

                  if (filterOrbitals)then

                     do ir=1,nrc
                        if(l .eq. 0) then
                           g(ir) = g(ir)
                        else
                           g(ir) = g(ir)*rofi(ir)**(l)
                        endif
                     enddo

                     minnorm = fdf_double("PAO.minnorm",0.999_dp)
                    write(6,"(a,f8.3)")
     $                 "paogen: Filtered orbital cutoff" //
     $                    " (Bohr^-1):",kmax
                     call filter(l,nrc,rofi(1:nrc),g(1:nrc),kmax,2,
     $                    minnorm)
                     
!Store the filtered orbital in g
                     g = 0.0_dp                   
                     do ir=2,nrc
                        if(l .eq. 0)then
                           g(ir) = g(ir) 
                        else
                           g(ir) = g(ir)/rofi(ir)**l
                        endif
                     enddo
                     g(1) = g(2)                   

!tail = 0.0_dp
!if(g(nrc) .ne. 0)then
!   tail = -g(nrc)  
!   g    = g + tail
!endif
          
                  endif

            call comPOL(is,a,b,rofi,g,l,
     .             nsm,rcpol(ipol,l,nsm),ipol,nrc,indx)

              enddo 
        

              endif    
            enddo  
          enddo 
      
          end subroutine polgen
!

               subroutine comPOL(is,a,b,rofi,rphi,
     .            l,nsm,rc,ipol,nrc,norb)
C****
C Generates the common block with the information about the polarization
C orbitals.
C  Written by D. Sanchez-Portal, Aug. 1998.
C  Modify by DSP, July 1999

               implicit none

               integer l, is, norb, nrc, ipol, nsm

               real(dp) rc, rphi(nrmax), a, b,
     .            rofi(nrmax)  

C****Internal variables
C
             integer itb, nr, nmax, nmin, nn
        
             real(dp) delt, r, phi, dy, yp1, ypn
             
****NUMBER OF POINTS USED BY POLINT FOR THE INTERPOLATION***
C
          integer npoint
          parameter(npoint=4)
C
          rcpoltb(ipol,l,nsm,is)=rc

CINTERPOLATION TO GENERATE TABLES WITH KB PROJECTORS**
C
            delt=rc/(dble(ntbmax-1)+1.0d-20) 

          if(delt.gt.deltmax) then
              write(6,'(a)')
     .    'comPOL: WARNING It might be a good idea to increase'
              write(6,'(a)')
     .    'comPOL: WARNING parameter ntbmax (in file atmparams.f) '
              write(6,'(a,i6)')
     .    'comPOL: WARNING to at least ntbmax = ',
     .        nint(Rc/deltmax)+2
          endif

            tabpol(1,norb,is)=delt

            tabpol(2,norb,is)=dble(l)


            do itb=1,ntbmax-1
                 r=delt*(itb-1)
                 nr=nint(log(r/b+1.0d0)/a)+1
                 nmin=max(1,nr-npoint)
                 nmax=min(nrc,nr+npoint)
                 nn=nmax-nmin+1
                 call polint(rofi(nmin),rphi(nmin),nn,r,phi,dy)
                 tabpol(itb+2,norb,is)=phi 
            enddo
            tabpol(ntbmax+2,norb,is)=0.0d0  
                       
C
C****

C****TABLE WITH THE SECOND DERIVATIVE ***
C

            yp1=huge(1.d0)
            ypn=huge(1.d0)

            call SPLINE( delt, tabpol(3:2+ntbmax,norb,is), ntbmax,
     &                   yp1, ypn, tab2pol(1:ntbmax,norb,is) )
            end subroutine compol
!
        subroutine set_mesh(a,b,rofi,drdi,s)

C**
C    Setting up mesh points an its derivatives from standard
C    values
C    D. Sanchez-Portal, Aug. 98
 
        implicit none

        real(dp) 
     .     rofi(nrmax), drdi(nrmax), s(nrmax), a, b


C**** Internal variables**
        real(dp) 
     .     aa, bb, zt, rpb, ea, ea2
        integer ir

        parameter(zt=1.0d0)
        parameter(aa=80.0d0)
        parameter(bb=6.0d0) 


C*STANDART VALUES FOR MESH PARAMETERS

          b=exp(-bb)/zt
          a=1.0d0/aa

C*SET UP THE MESH POINTS AND ITS DERIVATIVE***

          rpb=b
          ea=exp(a)
          ea2=1.0d0
          do ir=1,nrmax
            drdi(ir)=a*rpb
            rofi(ir)=b*(ea2-1.0d0)
            s(ir)=(a*rpb)**2
            rpb=rpb*ea
            ea2=ea2*ea
          enddo
          end subroutine set_mesh
!
          subroutine BESSEL(is,a,b,rofi,drdi,s,
     .             lmxo,
     .             nzeta,rco,lambda, norb) 

           use basis_specs, only: restricted_grid

C****
C  Caculates Bessel functions as a floating basis
C  Written by D. Sanchez-Portal, Aug. 1998.
C  Modify by DSP, July 1999

               implicit none

               real(dp)
     .         a, b, rofi(nrmax),
     .         drdi(nrmax), s(nrmax), 
     .         rco(nzetmx,0:lmaxd,nsemx),
     .         lambda(nzetmx, 0:lmaxd,nsemx)


               integer
     .           lmxo, is, nzeta(0:lmaxd,nsemx),
     .           norb


C***Internal variables**

               integer
     .           l,nprin, nnodes, nodd, nrc, ir,indx,
     .           izeta

               real(dp)
     .           rc, 
     .           dnrm, phi,
     .           g(nrmax), eorb, eps

               real(dp) :: v(nrmax)=0.0d0

             norb=0 
             indx=0
             do l=0,lmxo 

               if(nzeta(l,1).gt.0) then

            write(6,'(/2A,I2)')
     .       'Bessel: floating Bessel functions ',
     .           'with angular momentum L=',l

              do izeta=1, nzeta(l,1) 

C**Cut-off radius for Bessel functions must be an explicit input***
C
                if (rco(izeta,l,1).lt.1.0d-5) then 
                    write(6,'(a)')
     .     'Bessel: ERROR Zero cut-off radius with Z=-100 option'
                    write(6,'(a)')
     .     'Bessel: ERROR Cut-off radius must be explicitely specified'
                    write(6,'(a)')
     .     'Bessel: ERROR using Z=-100 (Floating Bessel functions)'
                  call die
 
                endif
C
C***

          if(abs(lambda(izeta,l,1)).lt.1.0d-3) lambda(izeta,l,1)=1.0d0
C*
           if(abs(lambda(izeta,l,1)-1.0d0).gt.1.0d-3) 
     .         then
             write(6,'(/,a)')
     . 'Bessel: WARNING Scale factor is not active with Z=-100 option' 
           endif 
           lambda(izeta,l,1)=1.0d0
C*

                  rc=rco(izeta,l,1)
                  nrc=nint(log(rc/b+1.0d0)/a)+1
                  if (restricted_grid) then
                     nodd=mod(nrc,2)
                     if(nodd.eq.0) then
                        nrc=nrc+1
                     endif
                  endif
                  rc=b*(exp(a*(nrc-1))-1.0d0)
                  rco(izeta,l,1)=rc
                  

                      nnodes=izeta
                      nprin=l+1
                      call schro_eq(1.0d0,rofi,v,v,s,drdi,
     .                  nrc,l,a,b,nnodes,nprin,
     .                  eorb,g) 
                       dnrm=0.0d0
                       do ir=2,nrc
                          phi=g(ir)
                          dnrm=dnrm+drdi(ir)*phi*phi
                          g(ir)=phi/(rofi(ir)**(l+1))
                       enddo 
                       g(1)=g(2)        

C****Checking normalization of the wavefunctions**
                 eps=1.0d-4
                 if(abs(dnrm-1.0d0).gt.eps) then
                   do ir=1,nrc
                   g(ir)=g(ir)/sqrt(dnrm)
                  enddo
                 endif
C****  



             write(6,'(/,(3x,a,i2),2(/,a25,f12.6))')
     .          'izeta =',izeta,
     .          'rc =',rco(izeta,l,1),
     .          'energy =',eorb  

            norb=norb+(2*l+1)
            indx=indx+1
            call comBasis(is,a,b,rofi,g,l,
     .           rco(izeta,l,1),lambda(izeta,l,1),izeta,
     .           1,nrc,indx)


              enddo 
        
              endif  

             enddo 

            end subroutine bessel
!
!
        subroutine energ_deriv(a,r,psi,vps,
     .      ve,drdi,nrc,l,el,psidev,nrval)
      
        implicit none

C*
C      This routine calculate the energy derivative of 
C      a given wavefunction.
C      The routine solve and inhomogeneus version of 
C      Schrodinger eqn.  
C      It is not an optimized algorithm!!!!!!!!!!!!!!!!!
C       Written by Daniel Sanchez-Portal, July 1999
C 
        integer  l, nrmin, nrval, ir, nrc   

        real(dp)  r(nrval),psi(nrval),psidev(nrval),
     .   el,vps(nrval),g(nrmax),drdi(nrmax),h(nrmax),ve(nrval), 
     .   hi, dnrm, cons, a, ortog, dnrm2
        
         parameter(nrmin=1) 

C

          
          nrc=min(nrc,nrval)
 
C Solving the inhomogeneus Schrodinger equation
          do ir=2,nrc
            hi=vps(ir)+ve(ir)+l*(l+1)/r(ir)**2-el
            hi=hi*(drdi(ir)**2)
            hi=hi+0.25d0*a**2
            h(ir)=hi
          enddo 
          h(1)=h(2)
          
          cons=psi(nrmin+1)/(vps(nrmin+1)+ve(nrmin+1)-el)
          cons=cons/r(nrmin+1)**(l+1) 
          g(1)=0.0d0
          do ir=1,nrmin+1
            g(ir)=cons*(r(ir)**(l+1))/sqrt(drdi(ir))
          enddo 

          do ir=nrmin+2,nrc
            hi=-((psi(ir)+10.0d0*psi(ir-1)
     .         +psi(ir-2))/12.0d0)

            hi=hi+(10.0d0*h(ir-1)*g(ir-1)+h(ir-2)*g(ir-2))/12.0d0
 
            hi=hi+2.0d0*g(ir-1)-g(ir-2)

            g(ir)=hi/(1.0d0-h(ir)/12.0d0)


          enddo 

C Orthogonalize the energy derivative to the original wavefunction
C and normalize
          dnrm2=0.0d0
          ortog=0.0d0
          do ir=1, nrc
            g(ir)=g(ir)*sqrt(drdi(ir))
            dnrm2=dnrm2+drdi(ir)*(psi(ir)**2)
            ortog=ortog+drdi(ir)*g(ir)*psi(ir)
          enddo
          dnrm=0.0d0
          do ir=1, nrc
             g(ir)=g(ir)-ortog*psi(ir)/dnrm2
             dnrm=dnrm+drdi(ir)*(g(ir)**2)
          enddo 
          dnrm=sqrt(dnrm)
          do ir=1,nrc
             psidev(ir)=g(ir)/dnrm
          enddo

          end subroutine energ_deriv
!
        subroutine rphi_vs_e(a,b,r,vps,
     .      ve,nrval,l,el,rphi,rmax)

           use basis_specs, only: restricted_grid

C**
C   Calculate the atomic 
C   radial wavefunction of the pseudopotential Vps, with angular
C   momentum  l, and energy el, inside r<Rmax
C   The Schrodinger equation is solved using a simple Numerov 
C   scheme. Rmax should not be taken too big. 
C   D. Sanchez-Portal, July 1999.
C**

        real(dp) a, b
        integer nrval
        real(dp) r(nrval),
     .   el,vps(nrval),g(nrmax),drdi(nrmax),h(nrmax),ve(nrval),
     .   rphi(nrval), rmax, dnrm

        real(dp) big, dexpa, ab, hi
        parameter (big=1.0d6)
        integer  l, nrc, jr, ir


        dexpa=exp(a)
        ab=a*b
        do ir=1,nrval
           drdi(ir)=ab
           ab=dexpa*ab
        enddo

        
          do ir=2,nrval
            hi=vps(ir)+ve(ir)+dble(l*(l+1))/r(ir)**2-el
            hi=hi*(drdi(ir)**2)
            hi=hi+0.25d0*a**2
            h(ir)=hi
          enddo
          h(1)=h(2)

         
          g(1)=0.0d0
          g(2)=1.0d0
          nrc=nint(log(rmax/b+1.0d0)/a)+1
          nrc=min(nrc,nrval)
          if (restricted_grid) nrc=nrc+1-mod(nrc,2)
          do ir=3,nrc

            hi=(10.0d0*h(ir-1)*g(ir-1)+h(ir-2)*g(ir-2))/12.0d0

            hi=hi+2.0d0*g(ir-1)-g(ir-2)

            g(ir)=hi/(1.0d0-h(ir)/12.0d0)
            
            if(abs(g(ir)).gt.big) then 
             dnrm=0.0d0
             do jr=1,ir
               dnrm=dnrm+drdi(jr)*(g(jr)*sqrt(drdi(jr)))**2
             enddo 
             dnrm=sqrt(dnrm)
             do jr=1,ir
               g(jr)=g(jr)/dnrm
             enddo 
            endif 
          enddo


C Normalize the wavefunction
          dnrm=0.0d0
          do ir=1, nrc
             g(ir)=g(ir)*sqrt(drdi(ir))
             dnrm=dnrm+drdi(ir)*(g(ir)**2)
          enddo
          dnrm=sqrt(dnrm)
          do ir=1, nrc
             rphi(ir)=g(ir)/dnrm
          enddo
          
          end subroutine rphi_vs_e

!
!   This is duplicate with redbasis.F...
!
        subroutine type_name(basistype)

C Written by D. Sanchez-Portal, Aug. 1998.

         character basistype*(*)

         if(basistype.eq.'NODES') then
               basistype='nodes'
         elseif(basistype.eq.'nodes') then
         elseif(basistype.eq.'NONODES') then
              basistype='nonodes'
         elseif(basistype.eq.'nonodes') then
         elseif(basistype.eq.'SPLIT') then
              basistype='split'
         elseif(basistype.eq.'split') then
         elseif(basistype.eq.'SPLITGAUSS') then
              basistype='splitgauss'
         elseif(basistype.eq.'splitgauss') then
         elseif(basistype.eq.'USER') then
              basistype='user'
         elseif(basistype.eq.'user') then
         else

              write(6,'(/,2a,(/,5(3x,a)),(/,2(3x,a)))')
     .        'type_name: Incorrect basis-type option specified,',
     .        ' active options are:',
     .        'NODES','SPLIT','USER','SPLITGAUSS','NONODES'

              call die

         endif

        end subroutine type_name
!
        subroutine prinput(ntotsp)

C**
C Prints the values of the parameters which have been actually 
C used in the generation of the basis (cut-off radius, contraction
C factors, and BasisType option to augment the basis set).
C The information is written in the same format as required for the 
C input file.
C Written by D. Sanchez-Portal, Oct. 1998.
C**

          implicit none 
          integer ntotsp

         character(len=10), parameter  :: basistype_default='split'

C***Internal variables
          integer is, nshell, l, lo, nzt, izt, nsm, i
          logical :: synthetic_atoms = .false.
          
          character basistype*10
          character(len=11) rcchar(nzetmx), lambdachar(nzetmx)

          type(ground_state_t), pointer ::  gs

             write(6,'(/a,58("-"))')
     .            'prinput: Basis input '

             basistype=fdf_string('PAO.BasisType',basistype_default)
             call type_name(basistype) 

             write(6,'(/2a)')'PAO.BasisType ',basistype

             write(6,'(/a)')
     .                   '%block ChemicalSpeciesLabel' 
             do is=1,ntotsp
                write(6,'(2(1x,i4),1x,2a)')
     .                  is,izofis(is),labelfis(is), 
     .              '    # Species index, atomic number, species label' 
                if (izofis(is) > 200) synthetic_atoms = .true.
             enddo 
             write(6,'(a)')
     .                  '%endblock ChemicalSpeciesLabel' 

!
!            Possible SyntheticAtoms block here
!
             if (synthetic_atoms) then
                write(6,'(/a)')        
     .               '%block SyntheticAtoms   # Valence config'

                do is=1,ntotsp
                   if (izofis(is) < 200) cycle
                   gs => basis_parameters(is)%ground_state
                   write(6,"(i3)") is
                   do i = 0, 3
                      write(6,fmt="(i3)",advance="no") gs%n(i)
                   enddo
                   write(6,*)
                   do i = 0, 3
                      write(6,fmt="(f9.5)",advance="no")
     $                      gs%occupation(i)
                   enddo
                   write(6,*)
                enddo 
                write(6,'(a)') "%endblock SyntheticAtoms"
             endif     ! synthetic_atoms
!
!
                
             write(6,'(/a)')        
     .   '%block PAO.Basis                 # Define Basis set'
             do is=1, ntotsp  
                
                nshell=0 
                lo=lmxosave(is)
                do l=0,lo
                   do nsm=1,nsemicsave(l,is)+1
                    if(nzetasave(l,nsm,is).ne.0) nshell=nshell+1
                   enddo 
                enddo 
               if(basistype_save(is).eq.basistype) then  

                if(abs(chargesave(is)).lt.1.0d-4) then 
                     write(6,'(a,1x,i2,20x,a)')
     .                labelfis(is), nshell,
     .                 '# Species label, number of l-shells'
                else
                    write(6,'(a,1x,i2,1x,f7.3,12x,2a)')
     .                labelfis(is), nshell, chargesave(is),
     .                  '# Label, l-shells,',
     .              ' ionic net charge'
                endif 
         
               else 

               if(abs(chargesave(is)).lt.1.0d-4) then
                     write(6,'(a,1x,i2,1x,a,10x,2a)')
     .             labelfis(is), nshell, basistype_save(is), 
     .                  '# Species label, l-shells,',
     .                  ' basis type '
                else
                    write(6,'(a,1x,i2,1x,a,1x,f7.3,1x,2a)')
     .              labelfis(is), nshell, basistype_save(is),
     .              chargesave(is),
     .                  '# Label, l-shells, type,',
     .              ' ionic net charge'
                endif 
           
               endif 
 
                   do l=0,lo
                     do nsm=1,nsemicsave(l,is)+1
                      nzt=nzetasave(l,nsm,is)
                      if(nzt.ne.0) then  
                        if(npolorbsave(l,nsm,is).gt.0)then
                          write(6,'(1x,a,i1,2(1x,i3),a,i3,19x,2a)') 
     .                       'n=',cnfigtb(l,nsm,is),
     .                       l, nzt, ' P ',npolorbsave(l,nsm,is),
     .                       '# n, l, Nzeta, ','Polarization, NzetaPol'
                        else
                          write(6,'(1x,a,i1,2(1x,i3),25x,a)')
     .                           'n=',cnfigtb(l,nsm,is),
     .                           l, nzt, '# n, l, Nzeta '
                        endif  
                        do izt=1, nzt
                           write(rcchar(izt),'(1x,f7.3)') 
     .                                          rcotb(izt,l,nsm,is)
                           write(lambdachar(izt),'(1x,f7.3)') 
     .                                       lambdatb(izt,l,nsm,is)
                        enddo
                        write(6,'(20a)')
     .                               (rcchar(izt), izt=1,nzt)
c    .                     ,'        # rc(izeta=1,Nzeta)(Bohr)'
                        write(6,'(20a)') 
     .                               (lambdachar(izt), izt=1,nzt)
c    .                     ,'        # scaleFactor(izeta=1,Nzeta)'
                     endif 
                    enddo
                   enddo  
                  enddo 
             write(6,'(a)')
     .                       '%endblock PAO.Basis' 


             write(6,'(/a,70("-")/)') 'prinput: '

             end subroutine prinput

!--------------------------------------------------------------
!
!      The famous "Vanderbilt generalized cutoff"
!
       function vander(a,x) result(f)
       real(dp), intent(in) :: a    ! Generalized gaussian shape
       real(dp), intent(in) :: x    
       real(dp)             :: f

       real(dp), parameter :: log10_e = 0.4343
       real(dp), parameter :: exp_range = (range(1.0_dp)-1)/log10_e

!!       real(dp), parameter :: exp_range = 40.0_dp
       real(dp)   :: gexp

       gexp = sinh(a*x)/sinh(a)
       gexp = gexp*gexp

       if (gexp .lt. exp_range) then
          f=exp(-gexp)  
       else
          f = 0.0_dp
       endif

       end function vander
!----------------------------------------------------- Soft Confinement
       subroutine build_vsoft(is,l,nsm,rinn,vcte,a,b,
     $                        rc,rofi,nrval,vsoft,plot)

       integer, intent(in)             :: is
       integer, intent(in)             :: l, nsm
       integer, intent(in)             :: nrval
       real(dp), intent(inout)         :: rinn
       real(dp), intent(in)            :: vcte
       real(dp), intent(in)            :: a, b
       real(dp), intent(in)            :: rc
       real(dp), intent(in)            :: rofi(*)
       real(dp), intent(out)           :: vsoft(*)
       logical, intent(in), optional   :: plot

       character(len=80)   :: filename
       integer             :: nrcomp, iu, ir
       real(dp)            :: rcsan, exponent
       logical             :: write_to_file

       write_to_file = .false.
       if (present(plot)) then
          write_to_file = plot
       endif
       

       nrcomp = nint(log(rc/b+1.0_dp)/a)+1
       rcsan = rofi(nrcomp+1) + 1.0e-6_dp
!
!      Fractional rinn indicated by a negative value
!
       if (rinn < 0.0_dp) rinn = -rinn*rcsan
            
       if (write_to_file) then
          write(filename,"(a,i1,a,i1,a,a)")
     $         trim(labelfis(is)) // ".L", l, ".", nsm,
     $         ".", "confpot"
          call io_assign(iu)
          open(unit=iu,file=filename,status='unknown')
          write(iu,'(2a)')
     .         '# Soft confinement potential for ', labelfis(is)
          write(iu,'(a,2i3)')
     .         '#    Soft confinement for shell l, nsm = ',l, nsm
          write(iu,'(a,f10.4,a)')
     .         '#        Inner radius    (r_inn) = ', rinn,' bohrs'
          write(iu,'(a,f10.4,a)')
     .         '#        External radius (r_ext) = ',   rcsan,' bohrs'
          write(iu,'(a,f10.4,a)')
     .         '#        Prefactor       (V_0)   = ', vcte,' Ry'
       endif

       do ir = 1, nrval
          if(rofi(ir) .le. rinn) then
             ! This avoids the singularity
             vsoft(ir) = 0.0_dp
          elseif (rofi(ir) .lt. rcsan) then
             exponent = -(rcsan-rinn)/(rofi(ir)-rinn)
             ! The maximum exponent is machine-dependent
             if( abs(exponent) .gt. 500.0_dp) then
                vsoft(ir) = 0.0_dp
             else
                vsoft(ir) = vcte / (rcsan-rofi(ir)) *
     .               exp(exponent)
             endif
          else
             vsoft(ir) = 0.0_dp
          endif
       enddo

       if (write_to_file) then
          do ir = 1, nrcomp + 1
             write(iu,'(2f20.7)') rofi(ir),vsoft(ir)
          enddo
          call io_close(iu)
       endif

       end subroutine build_vsoft

       subroutine find_split_location(nrc,rofi,
     $     drdi,rphi, split_table, l,spln,cons1,cons2,nsp)

       !
       ! Searches for the right point to fit, and
       ! fits a parabola.

       integer, intent(in)   :: nrc  ! Index of last point of orbital
       real(dp), intent(in)  :: rofi(*), drdi(*)
       real(dp), intent(in)  :: rphi(*) ! orb
       integer, intent(in)   :: l
       real(dp), intent(in)  :: spln
       real(dp), intent(in)  :: split_table(*)  ! Precomputed table
       real(dp), intent(out) :: cons1, cons2
       integer, intent(out)  :: nsp

       integer :: ir
       real(dp) :: parab_norm

       nsp =nrc
       do ir = nrc, 2, -1
          if (split_table(ir) > spln) then
             if (ir == nrc) then ! borderline case
                nsp = ir
             else
                ! Choose closest point
                if ( (split_table(ir)-spln) >
     $               (spln-split_table(ir+1)) ) then
                   nsp = ir + 1
                else
                   nsp = ir   
                endif
             endif
             exit
          endif
       enddo
       if (nsp == nrc) then
          write(6,"(a,f10.6)")
     $         "Split-norm parameter is too small, "
     $         // "(degenerate 2nd zeta): ",
     $          spln
          call die()
       endif
       if (nsp <= 2) then
          call die("Cannot find split_valence match point")
       endif

       call fit_parabola(nsp,rofi, drdi,rphi, l,
     $      cons1,cons2,parab_norm)

       end subroutine find_split_location

       subroutine fit_parabola(nsp,rofi,
     $     drdi,rphi, l, cons1,cons2,parab_norm)

       ! Fits C1 and C2 to match a parabola to rphi at
       ! point of index nsp

       integer, intent(in)   :: nsp
       real(dp), intent(in)  :: rofi(*), drdi(*)
       real(dp), intent(in)  :: rphi(*)
       integer, intent(in)   :: l
       real(dp), intent(out) :: cons1, cons2
       real(dp), intent(out) :: parab_norm

       real(dp) :: rsp, frsp, dfrsp

       rsp = rofi(nsp)
       frsp=rphi(nsp)/rsp
       dfrsp=0.5d0*(rphi(nsp+1)/rofi(nsp+1)
     .      -rphi(nsp-1)/rofi(nsp-1))
       dfrsp=dfrsp/drdi(nsp)
       cons1= 0.5d0*(dfrsp*rsp-l*frsp)/(rsp**(l+2))
       cons2= frsp/(rsp**l)-cons1*rsp**2

       call nrmpal(cons1,cons2,rsp,l,parab_norm)

       end subroutine fit_parabola
!---------------------------------
       subroutine split_scan(nrc,rofi, drdi, l, rphi, rnrm, label,
     $                       split_table,fix_split_table)
       !
       ! Scans the whole orbital, fitting parabolas and
       ! reporting effective norms to file

       integer, intent(in)   :: nrc  ! Index of last point of orbital
       real(dp), intent(in)  :: rofi(*), drdi(*)
       real(dp), intent(in)  :: rphi(*), rnrm(*) ! orb, norm(r)
       integer, intent(in)   :: l
       character(len=*), intent(in) :: label
       real(dp), intent(out) :: split_table(*)
       logical, intent(in)   :: fix_split_table

       integer :: ir, iu, jr
       character(len=50) :: fname
       real(dp) :: cons1, cons2, parab_norm, spln, rmin, rc, factor

       real(dp) :: split_table_raw(nrc)   ! Automatic array

       do ir = 3, nrc-1          ! Have to avoid /0 
          call fit_parabola(ir, rofi, drdi,rphi, l,
     $      cons1,cons2,parab_norm)
          spln = 1.0_dp - rnrm(ir)
          split_table_raw(ir) =  (spln+parab_norm)
       enddo
       split_table_raw(2) = split_table_raw(3)
       split_table_raw(1) = split_table_raw(2)
       split_table_raw(nrc) = split_table_raw(nrc-1)

       if (fix_split_table) then
          jr = nrc - 20         ! heuristic
          rmin = rofi(jr)
          rc = rofi(nrc)
          split_table(1:jr) = split_table_raw(1:jr)
          do ir = jr+1, nrc
             factor = dampen(rmin,rc,rofi(ir))
             split_table(ir) = factor*split_table_raw(ir)
          enddo
       else
          split_table(1:nrc) = split_table_raw(1:nrc)
       endif

       write(fname,"(3a,i1)") "SPLIT_SCAN.", trim(label), ".", l
       call io_assign(iu)
       open(iu,file=trim(fname),form="formatted",
     $      status="unknown",action="write",
     $      position="rewind")
       do ir = 1, nrc
          write(iu,"(i4,5f14.8)") ir, rofi(ir), rphi(ir),
     $         (1.0_dp - rnrm(ir)), split_table_raw(ir),
     $         split_table(ir)
       enddo
       call io_close(iu)

       end subroutine split_scan

      function dampen(a,b,r) result (y)
      real(dp), intent(in) :: a, b, r
      real(dp)             :: y

      real(dp)             :: x
      real(dp), parameter  :: tiny = 1.0e-12_dp

      x = (r-a)/(b-a)
      y = tanh(1.0_dp/(x+tiny) - 1.0_dp)
      end function dampen

      end module atom
