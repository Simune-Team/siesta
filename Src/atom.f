C $Id: atom.f,v 1.39.2.1 1999/06/03 10:03:58 wdpgaara Exp $

       SUBROUTINE ATOM (NTOTSP,IZIN,LMXKB,LMXO,
     .          NZETA,RCO,LAMBDA,ATM_LABEL,
     .          NPOLORB,SEMIC,LSEMIC,CHARGE,SMASS,BASISTYPE, 
     .          IS,NKB,NO,Q)

C***********************************************************************
C  Initializes Kleinman-Bylander pseudopotential and atomic orbitals.
C  Atomic orbitals basis sets are of several types.
C  **PAO's orbitals (Pseudo-atomic orbitals). We can form several
C   types of basis sets depending on how we double the basis set:
C    1) 'NODES': Double-Z, Triple-Z,... are orbitals of PAO 
C        type with one, two,... nodes
C    2) 'SPLIT': Double-Z, Triple-Z,... are generated from the 
C        original orbital, being equal to the PAO outside a 
C        radius Rm and given by (a*r^2+b)*r^l inside Rm.
C    3) 'NONODES': Double-Z, Triple-Z,... are orbitals of PAO
C        type with no-nodes, but different radii or different
C        scale-factors.
C    4) 'SPLITGAUSS':  Double-Z, Triple-Z,... are orbitals of GTO
C        (Gaussian Type Orbitals).
C    5) 'USER': Read from a file, provided by the user.
C
C  Written by D.Sanchez-Portal.   1996-1998.
C  New version.   August, 1998. 
C************************INPUT******************************************
C   INTEGER NTOTSP        : Expected total number of different chemical 
C                           species 
C   INTEGER IZIN            : Atomic number (if IZ is negative it will be 
C                          regarded as a set of floating orbitals 
C                          rather than a real atom) 
C   INTEGER LMXKB          : Angular momentum cutoff for
C                            Kleinman-Bylander nonlocal pseudopotential
C   INTEGER LMXO           : Angular momentum cutoff for basis orbitals 
C                            (not including the polarization base functions)
C   CHARACTER*20 ATM_LABEL: Label used to name all the atom related 
C                           files (i.e. pseudopotential files, in/output 
C                           files containing the basis set, etc..) 
C   INTEGER NZETA(L)      : Number of shells of PAOs with the same
C                            angular momentum. (i.e.not including 
C                            the polarization base functions)
C   REAL*8  RCO(NZETMX,L)    : Cutoff radius for Sankey-type basis orbitals
C   REAL*8  LAMBDA(NZETMX,L) : Scaling factor for the atomic orbital,
C                           for the 'SPLIT' type of basis it is
C                           interpreted as the exponent (in Bohr-2)
C                           of the gaussian orbitals for the
C                           the double-Z, triple-Z, etc.. 
C   INTEGER NPOLORB(L)      : Number of shells of polarization orbitals 
C                             for pseudo-atomic orbitals of momentum L.  
C   LOGICAL  SEMIC          : If .true. there is a shell of semicore states
C                             with angular momentum LSEMIC
C   INTEGER LSEMIC          : Angular momentum for the semicore states.          
C   REAL*8  CHARGE          : Charge state of the atom, only for basis set 
C                             generation purposes.            
C   REAL*8  SMASS           : Atomic mass for the species
C   CHARACTER*10 BASISTYPE     : Type of augmentation procedure for the 
C                             construction of the basis set
C***********************OUTPUT******************************************
C
C      INTEGER IS       : Species index assigned to atom
C      INTEGER LMXKB    : Angular momentum cutoff for
C                         Kleinman-Bylander nonlocal pseudopotential
C      INTEGER LMXO     : Angular momentum cutoff for basis orbitals 
C                             (including the polarization base functions)
C   REAL*8  RCO(NZETMX,L)    : Cutoff radius for Sankey-type basis orbitals
C   REAL*8  LAMBDA(NZETMX,L) : Scaling factor for the atomic orbital,
C                             for the 'SPLIT' type of basis it is
C                           interpreted as the exponent (in Bohr-2)
C                           of the gaussian orbitals for the
C                           the double-Z, triple-Z, etc..
C      INTEGER NKB      : Total number of Kleinman-Bylander projector 
C                         functions
C      INTEGER NO       : Total number of basis orbitals
C      INTEGER NZTOUT   : Actual value for NZT (if insufficient memory).
C      REAL*8  Q(NO)    : Neutral-atom occupations of basis orbitals, 
C                         with which pseudopotential is screened
C***********************RELATED ROUTINES********************************
C  EXTERNAL REAL*8  RCUT(IS,IO)      : Function which returns cutoff
C                                        radius of:
C                                        a) basis orbitals (IO > 0)
C                                        b) KB projectors  (IO < 0)
C                                        c) Local pseudopotential(IO=0)
C  EXTERNAL PHIATM(IS,IO,R,PHI,GRPHI): Routine which returns values and 
C                                        gradients of:
C                                        a) basis orbitals (IO > 0)
C                                        b) KB projectors  (IO < 0)
C                                        c) Local pseudopotential(IO=0)
C  EXTERNAL VLOCAL(IS,R,V,GRV)       : Routine which returns local part 
C                                        of neutral-atom pseudopotential
C                                        and its gradient
C  EXTERNAL RCORE(IS)                : Function which returns cutoff
C                                         radius of the pseudo-core
C  EXTERNAL CHCORE(IS,R,CH,GRCH)     : Routine which returns the 
C                                        pseudo-core for non-linear
C                                        core corrections in the
C                                        XC potential
C  EXTERNAL PSCH(IS,R,CH,GRCH)       : Routine which returns the
C                                        local-pseudotential charge
C                                        density (laplacian of the
C                                        local-pseudotential)
C  EXTERNAL PSOVER(IS1,IS2,R,E,DEDR) : Routine which returns the
C                                       correction to the electrostatic
C                                       interaction energy of the ions
C                                       due to the overlap of the
C                                       local-pseudotential charge
C                                       densities.
C  EXTERNAL REAL*8 EPSKB(IS,IO)      : Function which returns the KB
C                                        projector energies
C  EXTERNAL REAL*8 UION(IS)          : Function which returns the  
C                                      electrostatic self-energy of 
C                                      the charge density of the ion
C  EXTERNAL INTEGER IZOFIS(IS)       : Function which returns atomic
C                                       number 
C  EXTERNAL INTEGER IZVALFIS(IS)     : Function which returns the 
C                                      valence charge
C  EXTERNAL INTEGER LOFIO(IS,IO)     : Function which returns total 
C                                      angular momentum quantum number
C                                      of orbitals 
C  EXTERNAL INTEGER LOFIO(IS,IO)     : Function which returns magnetic
C                                      quantum number of the orbitals
C  EXTERNAL INTEGER LOMAXFIS(IS)     : Function which returns the maximum 
C                                      angular momentum of the basis functions
C                                      for a given specie
C  EXTERNAL INTEGER LMXKBFIS(IS)     : Function which returns the maximum 
C                                      angular momentum of the KB projectors
C                                      for a given specie.
C  EXTERNAL INTEGER NZTFL(IS,L)      : Function which returns the number of 
C                                      basis functions with a given angular 
C                                      momentum and for a given species
C  EXTERNAL INTEGER NOFIS(IS)        : Functions wich returns the total number 
C                                      of basis functions for the species is
C  EXTERNAL INTEGER NKBFIS(IS)       : Functions wich returns the total number
C                                      of KB projectors for the species is.
C  EXTERNAL REAL*8  MASSFIS(IS)      : Functions wich returns the atomic mass
C
C  EXTERNAL CHARACTER LABELFIS(IS)  : Functions wich returns the atomic label 
C  EXTERNAL CHARACTER SYMFIO(IS,IO) : Symmetry (s,px,py,pz,dz2,....) for 
C                                     the base functions. 
C                                      
C***********************UNITS*******************************************
C    Distances in Bohr.
C    Energies in Rydbergs.
C**************************BEHAVIOUR************************************
C  1) Several calls with the same input return the same output 
C   (in particular the same species index IS)
C  2) When called for the same atomic number but some other different
C   input, a different species index is assigned and a warning is 
C   printed.
C  3) IZ<0 is accepted, returning the orbitals of the corresponding 
C   atoms
C   but zero pseudopotential. If IZ=-100 square-box orbitals and zero 
C   pseudopot.
C  4) IZ=0 is also accepted, producing the re-initialization of tables,
C   i.e. setting to zero the number of tabulated tables. If this is 
C   done care must be taken to re-initialize also all routines which 
C   depend on RCUT, PHIATM, VLOCAL, etc... (like MATEL in O(N) program).
C***********************************************************************
C
C
C***************************PRECISION***********************************

                   implicit none      

C***********************************************************************
C
C
C**************************PARAMETERS (in file atom.h)******************
C INTEGER NSMAX     : Maximum number of species.
C INTEGER NTBMAX    : Maximum number of points in the tables defining
C                      orbitals,projectors and local neutral-atom 
C                      pseudopotential.
C INTEGER LMAXD     : Maximum angular momentum for both orbitals and 
C                      projectors.
C INTEGER  NZETMX   : Maximum number of PAO orbitals or polarization 
C                      orbitals with the same angular momentum and for
C                      the same species.
C INTEGER  NRMAX    : Maximum number of points in the functions read 
C                      from file '.vps' or '.psatom.data' (this number is
C                      determined by the parameter nrmax in the 
C                      program atm, which generates the files with 
C                      the pseudopotential information).
C**********************************************************************
C
         include 'atom.h'
C
C**********************************************************************

       character atm_label*20 , basistype*10
       
       integer izin, lmxo, lmxkb, lsemic, ntotsp, is, nkb, no
  
       integer 
     .    nzeta(0:lmaxd), npolorb(lmaxd) 

       integer maxos
           parameter (maxos=2*nzetmx*lmx2)
 
       double precision
     .  rco(nzetmx,0:lmaxd), lambda(nzetmx,0:lmaxd),  
     .  q(maxos), charge,smass

       logical semic
         
           
C**************Some internal parameters********************************

C******Maximum distance (in Bohrs) between points where function ******
C******is evaluated to generate the tables for intepolation************
C******If this distance is exceeded  a warning is printed**************  
         double precision  deltmax, deltmx 
         parameter (deltmx=0.05d0)
         common/cmdelt/deltmax


C****Default energy-shift to define the cut off radius of orbitals*****
C****In Rydbergs*******************************************************
        double precision eshift_default,eshift_deflt
        parameter(eshift_deflt=0.02d0)


C****Default norm-percentage for the automatic definition of **********
C***********multiple-zeta orbitals with the 'SPLIT' option*************
        double precision splnorm_default,splnorm_deflt
        parameter(splnorm_deflt=0.15d0)
       
C***We put these parameters in a common block to do it accesible to 
C*****other subroutines********************************************
C
        common/cmeshdefault/eshift_default
        common/cmspldefault/splnorm_default
C**********************************************************************
          integer ns2
          parameter(ns2=((nsmax+1)*nsmax)/2)

C**********************************************************************
C
C
C**********************************************************************



C*******************Variables in common blocks*************************


       double precision
     .  rcotb(nzetmx,0:lmaxd,nsmax), 
     .  rcpoltb(nzetmx,lmaxd,nsmax),
     .  lambdatb(nzetmx,0:lmaxd,nsmax),
     .  qtb(maxos,nsmax), slfe(nsmax),
     .  rctb(0:lmaxd,nsmax),smasstb(nsmax),
     .  chargesave(nsmax),qltb(0:3,nsmax)

      double precision
     .  table((ntbmax+2),-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tab2(ntbmax,-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tabpol((ntbmax+2),nzetmx*lmaxd,nsmax),
     .  tab2pol(ntbmax,nzetmx*lmaxd,nsmax),
     .  coretab(ntbmax+1,2,nsmax), 
     .  corrtab((ntbmax+1),2,ns2),
     .  chloctab((ntbmax+1),2,nsmax)
      
      

      integer
     .  izsave(nsmax),lmxkbsave(nsmax),
     .  lmxosave(nsmax),
     .  npolorbsave(lmaxd,nsmax),
     .  lsemicsave(nsmax), nzetasave(0:lmaxd,nsmax),
     .  nomax(nsmax), nkbmax(nsmax), ismax,
     .  izvaltb(nsmax), cnfigtb(0:3,nsmax) 

        
       character
     .   label_save(nsmax)*20, basistype_save(nsmax)*10 

       logical
     .    semicsave(nsmax)

       

C********Internal variables*********************************************
 

       double precision
     .  rofi(nrmax), drdi(nrmax), s(nrmax),
     .  vps(nrmax,0:lmaxd), rphi(nrmax,0:lmaxd),
     .  vlocal(nrmax), vxc(nrmax), ve(nrmax),
     .  rho(nrmax), chcore(nrmax), auxrho(nrmax),
     .  vePAO(nrmax), qPAO(0:lmaxd)


       double precision 
     .  pi, a, b, zval, flting,
     .  ex, ec, dx, dc, r2
     
      
       double precision
     .  Rgauss, Rgauss2, Rchloc

       character
     .   icorr*2, irel*3, nicore*4, 
     .   symbol*2,paste*50,
     .   xcfunc*3, xcauth*4



       integer 
     .  iz, nrval, ir , nrgauss, nchloc, nzcontr, l, nVna,
     .  irelt, lun
    
      
       logical
     .  called, new
   
C***********Declarations for FDF package *******************************
 
         include 'fdf/fdfdefs.h'
   
C*********************COMMON DECLARATIONS******************************

             common/cmlmxo/lmxosave
             common/cmlmxkb/lmxkbsave
             common/cmlabel/label_save
             common/cmsemic/semicsave,lsemicsave
             common/cmzeta/nzetasave
             common/cmiz/izsave
             common/cmpolorb/npolorbsave
             common/cmradorb/rcotb
             common/cmlambda/lambdatb 
             common/cmradpol/rcpoltb
             common/control/ismax,nomax,nkbmax 
             common/cmtab/table,tabpol
             common/cmspline/tab2,tab2pol
             common/cmradkb/rctb
             common/cizvaltb/izvaltb
             common/cmslfe/slfe
             common/cmcore/coretab
             common/cmchloc/chloctab
             common/cmcorr/corrtab 
             common/cmmass/smasstb
             common/cmq/qtb 
             common/cmql/qltb
             common/cmcnfig/cnfigtb
             common/cmcharge/chargesave
             common/cmbasistype/basistype_save 

C**********************SAVE FOR NEXT CALL******************************

                  save called

C**********************************************************************

                  data  called /.false./
 

C********************      PI     *************************************
                  pi=dacos(-1.0d0)

C**********************************************************************





C*******************INITIALITATION IN FIRST CALL*********************
C     
           if (.not.called) then 

C****SIZE OF THE BIGGER ARRAYS IN COMMON BLOCKS******************
                  call prmem( 0, 'atom', 'table', 'D',
     .                (ntbmax+2)*((lmaxd+1)*(nzetmx+1)+1)*nsmax )
                  call prmem( 0, 'atom', 'tab2',  'D',
     .                 ntbmax*((lmaxd+1)*(nzetmx+1)+1)*nsmax )
                  call prmem( 0, 'atom', 'tab2pol',  'D',
     .                (ntbmax+2)*(lmaxd+1)*nzetmx*nsmax )
                  call prmem( 0, 'atom', 'tab2pol',  'D',
     .                ntbmax*(lmaxd+1)*nzetmx*nsmax )
                  call prmem( 0, 'atom', 'coretab',  'D',
     .                 2*(ntbmax+1)*nsmax )
                  call prmem( 0, 'atom', 'corrtab',  'D',
     .                 2*(ntbmax+1)*ns2 ) 
                  call prmem( 0, 'atom', 'chloctab',  'D',
     .                 2*(ntbmax+1)*nsmax ) 

                  call prmem( 0, 'atom', ' ',' ',0) 
C**********************************************************************

                  called=.true.
           endif                
C
C**********************************************************************


C****************Print some information about the atomic species*******
C****************and selected options**********************************
C
           iz=izin

           write(6,'(/a,73(1h*),/)') 'ATOM: '
           if (iz.gt.0) then  
             write(6,'(3a,i4,a)')
     .       'ATOM: CALLED FOR ', symbol(iz), '  (Z =', iz,')' 

           elseif((iz.lt.0).and.(iz.ne.-100)) then  

             write(6,'(3a,i4,a,a)')
     .       'ATOM: Called for ', symbol(abs(iz)), '  (Z =', iz,')',
     .       ' ( Floating basis ) ' 

           elseif(iz.eq.-100) then
             write(6,'(a,i4,a)')
     .       'ATOM: Called for Z=',iz,'( Floating Bessel functions)'  

           elseif(iz.eq.0) then 

              write(6,'(a)') 'ATOM: re-initialization of the tables'
              write(6,'(a)') 'ATOM: all data will be set to zero'
 
          
           endif  
C
C**********************************************************************


C***Checking if the has been called previously for the same specie*****
C
           call new_specie(iz,lmxkb, lmxo,
     .         nzeta, rco, lambda, atm_label,
     .         npolorb, semic, lsemic,
     .         is, new, no, nkb, q) 
C
C***************************************************************************


C****If iz=0 the tables haves been initializated and we exit****************
C
              if(iz.eq.0) then 
                 called=.false.
                 write(6,'(a,73(1h*))') 'ATOM: '
                 return
              endif 
C
C***************************************************************************
 

C******If it was not a new species there is nothing more to do, so we exit** 
C
              if(.not.new) then 
                  write(6,'(a,73(1h*))') 'ATOM: ' 
                  return
              endif                 
C
C***************************************************************************

 

C*****IF IT WAS A NEW SPECIES WE START HERE ALL THE CALCULATIONS************
     
C***********Save the type of basis set generation procedure*****************
C
           basistype_save(is)=basistype
C
C***************************************************************************



C***********Save atomic mass in a common block****************************** 
C
          smasstb(is)=smass
C
C***************************************************************************


C**Maximum distance between consecutives points in the interpolation tables**
C
            deltmax=deltmx
C
C***************************************************************************

        if (iz.ne.-100) then

C Default values of some parameter for the common block**********************
C
            eshift_default=eshift_deflt
            splnorm_default=splnorm_deflt         
C
C***************************************************************************
     
C******************Reading pseudopotentials********************************* 
C
               call read_vps(atm_label, lmxo, lmxkb,
     .             nrval,a,b,rofi,drdi,s,vps,
     .             rho, chcore, zval, nicore, irel, icorr) 
C
C*************************************************************************** 


C*************** Save readed valence charge*********************************
C This can be different from the standard one if we have included semicore
C states in the valence shell.
C
            izvaltb(is)=nint(zval)
C
C***************************************************************************



C**IF IZ IS NEGATIVE WE WILL HAVE FLOATING ORBITALS IN THE ATOMIC POSITIONS*
C****************IF IZ POSITIVE WE WILL HAVE REAL ATOMS*********************
C
           flting=dsign(1.0d0,dble(iz))
           iz=abs(iz)
C
C***************************************************************************




C****************Common block with pseudocore information*******************
C*****used for non-linear-core-correction for exchange-correlation energy***
C
               call  comcore(is,a,b,rofi,chcore,
     .                       nrval,nicore,flting)
C
C***************************************************************************




C*CALCULATION OF THE VALENCE SCREENING POTENTIAL FROM THE READ CHARGE
C***********************DENSITY*****************************************
C 
CFor Kleinman-Bylander projectors calculation***************************
          call vhrtre(rho,ve,rofi,drdi,s,nrval,a)  
CFor PAO basis functions calculations*********************************** 
          do ir=2,nrval 
              auxrho(ir)=(1.0d0-charge/zval)*rho(ir) 
          enddo 
          call vhrtre(auxrho,vePAO,rofi,drdi,s,nrval,a) 
          chargesave(is)=charge
C
C***************************************************************************
C*************ADD THE EXCHANGE-CORRELATION POTENTIAL*******************
C*****Choosing the adecuate functional for the exchange-correlation****
C
       xcfunc = fdf_string('xc.functional','LDA')
       xcauth = fdf_string('xc.authors','PZ')
        
            call xc_check(xcfunc,xcauth,icorr)
C
C***************************************************************************
C
            if(irel.eq.'rel') irelt=1
            if(irel.ne.'rel') irelt=0 

          do ir=2,nrval
            r2=rofi(ir)**2
            r2=4.0d0*pi*r2
            dc=rho(ir)/r2
            if(nicore.ne.'nc  ')  dc=dc+chcore(ir)/r2
              auxrho(ir)=dc
          enddo

          r2=rofi(2)/(rofi(3)-rofi(2))
          auxrho(1)=auxrho(2) -(auxrho(3)-auxrho(2))*r2


          call atomxc(xcfunc,xcauth,irelt,
     .             nrval,nrmax,rofi,1,auxrho,
     .             ex,ec,dx,dc,vxc)


          do ir=1,nrval
             ve(ir)=ve(ir)+vxc(ir) 
          enddo  


          do ir=2,nrval
            r2=rofi(ir)**2
            r2=4.0d0*pi*r2
            dc=(1.0d0-charge/zval)*rho(ir)/r2
            if(nicore.ne.'nc  ')  dc=dc+chcore(ir)/r2
              auxrho(ir)=dc
          enddo

          r2=rofi(2)/(rofi(3)-rofi(2))
          auxrho(1)=auxrho(2) -(auxrho(3)-auxrho(2))*r2


          call atomxc(xcfunc,xcauth,irelt,
     .             nrval,nrmax,rofi,1,auxrho,
     .             ex,ec,dx,dc,vxc)


          do ir=1,nrval
             vePAO(ir)=vePAO(ir)+vxc(ir) 
          enddo 

C
C***************************************************************************







C*************************************************************************
C**Now, we are going  to calculate the local pseudopotential and the  ****
C**KB projectors, this is only necessary if we are dealing with real  ****
C**atoms (not with floating orbitals), i.e. only if flting is greater ****
C**than zero**************************************************************
C
          if(flting.gt.0.0d0) then
 


C***Rgauss is the maximum cut-off radius used in the pseudopotential   *****
C***generation, Rgauss is determined by comparison of the pseudopot.   *****
C***Corresponding to different l ( this is not possible is we have     *****
C***just one pseudopotential)                                          *****
C***Rgauss2 is the radius where the pseudopotentials reach  the        *****
C***asymptotic behaviour 2*Zval/r.                                     *****
C***For just one pseudopotential Rgauss is taken equal to Rgauss2      ***** 
C
         call  radii_ps(vps,rofi,Zval,nrval,lmxkb,
     .          nrgauss, rgauss, rgauss2)
C
C*************************************************************************** 
C******************Calculate local pseudopotential**************************
C
            if (rgauss2.gt.1.30d0*rgauss) then


CIn this case the atom core is so big that we do not have an asymptotic
Cof 2*Zval/r until Rgauss2 > Rc . To retain the same asymptotic
Cbehaviour as in the pseudopotentials we modified the definition
Cof the local potential
C
         write(6,'(/,a,f10.5)') 'ATOM: Estimated core radius ',
     .           rgauss2
         if (nicore.eq.'nc ')
     .   write(6,'(/,2a)') 'ATOM: Include non-local core corrections',
     .   ' could be a good idea'
           call vlocal2(Zval, nrval, a, rofi, drdi, s, vps(1,0),
     .              nrgauss,vlocal,nchloc,chcore)
C
C***************************************************************************
C
              else
C
CIn this case the pseudopotential reach to an asymptotic behaviour 2*Zval/r  
Cfor a radius approximately equal to Rc. 
C
            call vlocal1(Zval, nrval, a, rofi, drdi, s, rgauss,
     .                     vlocal,nchloc,chcore)
C
C***************************************************************************

            endif 
C
C***************************************************************************


C***************************************************************************
C*****Common block with local-pseudopotential charge************************
C*****used for non-linear-core-correction for exchange-correlation energy***
C 
          Rchloc=rofi(nchloc)
          write(6,'(2a,f10.5)') 'ATOM: Maximum radius for' ,
     .      ' local-pseudopot. charge ',Rchloc

            call  comlocal(is,a,b,rofi,chcore,nchloc,flting)
C
C***************************************************************************



C**************Write a file to plot the local potential**********************
C
            call io_assign(lun)
         open(lun,file=paste(atm_label,'.vlocal'),
     .       status='unknown')
          do ir=1,nrval
          write(lun,*) rofi(ir),vlocal(ir),chcore(ir)
          enddo
          call io_close(lun) 
C
C********************************************************************

C********ARRAY S FOR THE SCHRODINGER EQ. INTEGRATION***************** 
C
         do ir=2,nrval
             s(ir)=drdi(ir)*drdi(ir)
         enddo 
         s(1)=s(2)
C
C********************************************************************





C******CALCULATION OF THE KLEINMAN-BYLANDER PROYECTOR FUNCTIONS*******
C
            call KBgen(is, a,b,rofi,drdi,s,
     .         vps, vlocal, ve, nrval, Zval, lmxkb, nkb) 
            nkbmax(is)=nkb
C
C********************************************************************



            elseif(flting.lt.0.0d0) then 

C*******None Kleinman-Bylander projectors if floating orbitals*******
C         
            nkb=0 
            nkbmax(is)=0
C
C********************************************************************


C*********Zero local pseudopotential if floating orbitals*******
C
            nchloc=0
            call  comlocal(is,a,b,rofi,chcore,nchloc,flting)
C
C********************************************************************



C********ARRAY S FOR THE SCHRODINGER EQ. INTEGRATION*****************
C         (Needed for the calculation of the basis set )
         do ir=2,nrval
             s(ir)=drdi(ir)*drdi(ir)
         enddo
         s(1)=s(2)
C
C********************************************************************

            endif  





C*********CONSTRUCTION OF THE BASIS ORBITALS*************************


C********Method for the augmentation of the basis set****************
C
         if(basistype.ne.'user') then
           write(6,'(a,73(1h-))') 'ATOM: '
           write(6,'(/,a)') 'ATOM: SANKEY-TYPE ORBITALS:'
           nzcontr=0
           do l=0,lmxo
              if(nzeta(l).gt.1) nzcontr=1
           enddo
           if( nzcontr.eq.1) then
               write(6,'(2a)')
     .        'ATOM: Selected multiple-zeta basis: ',basistype
           endif
         endif

C
C********************************************************************

C*************Generate the PAOs basis set ***************************
C
         call Basis_gen(Zval,is, a,b,rofi,drdi,s,
     .                   vps, ve, vePAO, nrval, lmxo,
     .                   nzeta, rco, lambda, npolorb,
     .                   basistype, rphi, no) 
 
         write(6,'(a,i3)')
     .      'ATOM: Total number of Sankey-type orbitals:', no
           nomax(is)=no
C
C********************************************************************



           if(flting.gt.0.0d0) then 
 
C***Calculate initial populations for the atomic orbitals***********
C
           call atm_pop(is,iz,q,qPAO,lmxo,
     .       nzeta,semic,lsemic,npolorb) 
C
C********************************************************************

C***Screening of the atomic local pseudopotential********************
C
           call  Vna(is,Zval,qPAO,rphi,rco,vlocal,
     .        a,b,rofi,drdi,nrval,lmxo,nVna)
C
C********************************************************************
        
C CALCULATION OF THE ELECTROSTATIC SELF-ENERGY Of THE ****************
C ******** LOCAL PSEUDOPOTENTIAL CHARGE DENSITY ********************** 
C
          call slfe_local(slfe(is),vlocal,rofi,a,nVna,drdi) 
C
C********************************************************************
 
           elseif(flting.lt.0.0d0) then 

C*******Populations are zero because we have floating orbitals,******
C*******not real atoms***********************************************
C   
               do ir=1,no
                  q(ir)=0.0d0
               enddo 
C
C********************************************************************


C*********Zero neutral-atom pseudopotential if floating orbitals*******
C 
              nVna=0
              call comVna(is,a,b,rofi,vlocal,nVna,flting)
C
C********************************************************************

C**********Zero self-energy for the local pseudopotential charge****
C
              slfe(is)=0.0d0
C
C********************************************************************  
            endif 

          
                  
C********************************************************************
C********************************************************************
C********************************************************************
C********************************************************************
C********************************************************************
C********************************************************************
C********************************************************************
C********************************************************************
C********************************************************************
C********************************************************************

        elseif(iz.eq.-100) then


C***********SET UP THE MESH POINTS AND ITS DERIVATIVE******************
C
               call set_mesh(a,b,rofi,drdi,s)
C
C**********************************************************************
      

C***Local pseudopotential, KB projector, pseudocore charge density, etc...
C****blocks do not contain any information in the case of floating Bessel
C****functions but they have to be intialize in anycase
C
           flting=dsign(1.0d0,dble(iz)) 
C
C**** Zero valence charge for floating Bessel functions***********
C
            izvaltb(is)=0
C
C****************Zero pseudocore charge density******************
C
           call  comcore(is,a,b,rofi,chcore,
     .                       nrval,nicore,flting) 
C
C*******None Kleinman-Bylander projectors if floating orbitals*******
C         
            nkb=0 
            nkbmax(is)=0
C
C*********Zero local pseudopotential if floating orbitals*******
C
            nchloc=0
           call  comlocal(is,a,b,rofi,chcore,nchloc,flting)    
C
C*******Populations are zero because we have floating orbitals,******
C*******not real atoms***********************************************
C   
               do ir=1,no
                  q(ir)=0.0d0
               enddo 
C
C*********Zero neutral-atom pseudopotential if floating orbitals*******
C 
              nVna=0
              call comVna(is,a,b,rofi,vlocal,nVna,flting)
C
C**********Zero self-energy for the local pseudopotential charge****
C
              slfe(is)=0.0d0
C
C***************************************************************************


C******Calculate Bessel functions******************************************
C
             call Bessel (is,a,b,rofi,drdi,s,
     .             lmxo,nzeta,rco,lambda, no)

 
         write(6,'(/a,i3)')
     .      'ATOM: Total number of floating Bessel orbitals:', no
           nomax(is)=no

C
C*************************************************************************


           endif 



C*****Now atom has been called for all the species present in the 
C*****calculation
         if(is.eq.ntotsp) then 
C**********************************************************************
C     Writing basis information into a file 
C**********************************************************************
           call draw_basis(ntotsp)
           call prinput(ntotsp)   
           write(6,'(a,73(1h*))') 'ATOM: '
        endif  

        return

        end




