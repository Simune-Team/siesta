       SUBROUTINE ATOM (NTOTSP,IZIN,LMXKB,LMXO,MAXL,
     .          NZETA,NZT,NZTOUT,
     .          RCO,LAMBDA,
     .          IS,ATM_LABEL,NKB,NO,Q,MAXOS)

C***********************************************************************
C  Initializes Kleinman-Bylander pseudopotential and atomic orbitals.
C  Atomic orbitals basis sets are of seven types :
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
C    5) 'USER': Readed from a file, provided by the user.
C
C  Written by D.Sanchez-Portal.   Oct, 1996.
C  Strongly modified by D.Sanchez-Portal.   March, 1997
C  Strongly modified again by D.Sanchez-Portal.  July-August, 1997
C  the generation procedure of the SPLIT basis set has change 
C  and it has been included a new procedure to generate polarization 
C  orbitals.
C************************INPUT******************************************
C   INTEGER NTOTSP        : Expected total number of different chemical 
C                           species 
C   INTEGER IZ            : Atomic number (if IZ is negative it will be 
C                          regarded as a set of floating orbitals 
C                          rather than a real atom)
C   INTEGER MAXL          : Maximum angular momentum for atomic orbitals
C                            as defined in 'siesta.h' parameter file.
C   INTEGER NZT           : Maximum NZETA in the main program 
C                            (needed to define the RCO matrix)
C   CHARACTER*20 ATM_LABEL: Label used to name all the atom related 
C                           files (i.e. pseudopotential files, in/output 
C                           files containing the basis set, etc..) 
C   INTEGER MAXOS         : Maximum number of orbitals per atom in 
C                           calling routine
C***********************OUTPUT******************************************
C
C      INTEGER IS       : Species index assigned to atom
C      INTEGER LMXKB    : Angular momentum cutoff for
C                         Kleinman-Bylander nonlocal pseudopotential
C      INTEGER LMXO     : Angular momentum cutoff for basis orbitals
C   INTEGER NZETA(L)      : Number of different basis orbitals for each
C                            angular momentum and magnetic number
C   REAL*8  RCO(NZT,L)    : Cutoff radius for Sankey-type basis orbitals
C   REAL*8  LAMBDA(NZT,L) : Scaling factor for the atomic orbital,
C                           for the 'SPLIT' type of basis it is
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
C  EXTERNAL IZOFIS(IS)               : Function which returns atomic
C                                       number
C  EXTERNAL LOFIO(IS,IO)             : Function which returns total 
C                                      angular momentum quantum number
C                                      of orbitals
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

       implicit double precision (a-h,o-z)

C***********************************************************************
C
C
C**************************PARAMETERS***********************************
C INTEGER NSMAX     : Maximum number of species.
C INTEGER NTBMAX    : Maximum number of points in the tables defining
C                      orbitals,projectors and local neutral-atom 
C                      pseudopotential.
C INTEGER LMAXD     : Maximum angular momentum for both orbitals and 
C                      projectors.
C INTEGER  NZETMX   : Maximum number of orbitals with the same angular
C                      momentum and for the same species.       
C INTEGER  NRMAX    : Maximum number of points in the functions read 
C                      from file '.vps' or '.psatom.data' (this number is
C                      determined by the parameter nrmax in the 
C                      program atm, which generates the files with 
C                      the pseudopotential information).
C**********************************************************************

        include 'atom.h'
        parameter(ns2=((nsmax+1)*nsmax)/2)



C****NUMBER OF POINTS USED BY RATINT FOR THE INTERPOLATION************* 

        parameter(npoint=4)


C*****NUMBER OF CONTINOUS DERIVATIVES FOR VLOCAL, IN CASE IT IS******** 
C*****GENERATED MATCHING THE ORIGINALS PSEUDOPOTENTIALS AT Rc(ps)******
C******************ndevfit can value 2 or 3****************************
C**********with ndevfit=2 usually we obtain a smoother vlocal**********

c       parameter(ndevfit=2)

C**********************************************************************



C****Default energy-shift to define the cut off radius of orbitals*****
C****In Rydbergs*******************************************************

        parameter(eshift_default=0.02d0)

C****Default norm-percentage for the automatic definition of **********
C***********multiple-zeta orbitals with the 'SPLIT' option*************

        parameter(splnorm_default=0.15d0)

C**********************************************************************

C********Default number of polarization orbitals***********************
        parameter(numb_pol_default=1) 
C**********************************************************************
C
C
C**********************************************************************

       character
     .   namatm*2, icorr*2, irel*3, nicore*4, 
     .   method(6)*10,text(10)*10,symbol*2,paste*50,fname*50,
     .   xcfunc*3, xcauth*4,basistype*15,basistype_default*15,
     .   ipol*5,
c    .   loctype*3,
     .   atm_label*20,label_old(nsmax)*20,
     .   basis_size*15,basis_size_default*15,
     .   gen_basis*15

       double precision
     .  rofi(nrmax),vps(nrmax,0:lmaxd),rphi(nrmax,0:lmaxd),
     .  ve(nrmax),rc(0:lmaxd),rco(nzt,0:lmaxd),q(maxos),
     .  ql(0:lmaxd),
     .  drdi(nrmax),elocal(2,0:lmaxd),s(nrmax),h(nrmax),
     .  g(nrmax),y(nrmax),rho(nrmax),eigen(0:lmaxd),dkbcos(0:lmaxd),
     .  table((ntbmax+2),-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  rcoold(nzetmx,0:lmaxd,nsmax),
     .  rcoin(nzetmx,0:lmaxd,nsmax),rctb(0:lmaxd,nsmax),
     .  rcotb(nzetmx,0:lmaxd,nsmax),vlocal(nrmax), 
     .  tab2(ntbmax,-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),aux(ntbmax),
     .  qold(lmx2*nzetmx,nsmax),slfe(nsmax),lambda(nzt,0:lmaxd),
     .  lambdain(nzetmx,0:lmaxd,nsmax),lambdaold(nzetmx,0:lmaxd,nsmax),
     .  chcore(nrmax),dknrm(0:lmaxd),coretab(ntbmax+1,2,nsmax),
     .  vxc(nrmax),auxrho(nrmax),chloctab((ntbmax+1),2,nsmax),
     .  corrtab((ntbmax+1),2,ns2),rnrm(nrmax,0:lmaxd),
     .  saveeorb(0:lmaxd,nzetmx)



       integer 
     .  iz,nkb,no,lmxkb,lmxo,nzeta(0:lmaxd),izold(nsmax),
     .  lkbold(nsmax),loold(nsmax),nzetold(0:lmaxd,nsmax),
     .  lkbin(nsmax),loin(nsmax),nzetin(0:lmaxd,nsmax),
     .  isold,nzettb(0:lmaxd,nsmax),iztb(nsmax),ismax,
     .  nomax(nsmax),nkbmax(nsmax),config(0:4),loctab(nsmax),
     .  lmxPAOold(nsmax),two
      
       logical
     .  called,found,polorb,polorb_default
   

       include 'fdf/fdfdefs.h'

C******LOCAL POTENTIAL TYPE (OLD (GAUSSIAN), NEW (SMOOTHER) )**********


c           parameter(loctype='new')
c           parameter(loctype='old')

C****Default size of the automatic basis set***************************
         parameter(basis_size_default='standard')
C*********************************************************************

C****Default polarization orbitals************************************
         parameter(polorb_default=.false.)
C*********************************************************************

C****Default type of basis set*****************************************
         parameter(basistype_default='split')
C**********************************************************************
   
C*********************COMMON DECLARATIONS******************************

        common/cmtab/table
        common/cmradkb/rctb
        common/cmradorb/rcotb
        common/cmzeta/nzettb
        common/cmspline/tab2
        common/ciztb/iztb
        common/control/ismax,nomax,nkbmax
        common/cmslfe/slfe
        common/cmloc/loctab
        common/cmcore/coretab
        common/cmchloc/chloctab
        common/cmcorr/corrtab
 
C**********************************************************************
C
       data  called /.false./
C
C**********************SAVE FOR NEXT CALL******************************

        save called,isold,izold,lkbold,loold,nzetold,rcoold,qold
        save label_old,lkbin,loin,nzetin,rcoin,lambdain,lambdaold
        save lmxPAOold  

C**********INITIALITATION IN FIRST CALL********************************
     
        if (.not.called) write(6,'(/,a,73(1h*))') 'ATOM: '
           iz=izin
           if (iz.gt.0) then  
             write(6,'(3a,i4,a,/)')
     .       'ATOM: Called for ', symbol(iz), '  (Z =', iz,')'
           elseif((iz.lt.0).and.(iz.ne.-100)) then 
             write(6,'(3a,i4,a,a,/)')
     .       'ATOM: Called for ', symbol(abs(iz)), '  (Z =', iz,')',
     .       ' ( Floating basis ) '
           elseif(iz.eq.-100) then
             write(6,'(a,i4,a,/)')
     .       'ATOM: Called for Z=',iz,'( Floating Bessel functions)'
           endif  
        if (.not.called) then

          isold=0
          called=.true.
          goto 9

        endif


C***** IS THIS A NEW SPECIES ? ****************************************
C***** We first compare the inputs, if all the inputs are equal then***
C***** this is not a new species **************************************  

          do 1 ns=1,isold    
             if((izold(ns).eq.iz).and.
     .               (label_old(ns).eq.atm_label)) then
               nsold=ns
               goto 5
             endif 

  1       continue


          goto 9

   
  5       if(lmxkb.ne.lkbin(nsold)) goto 8
          if(lmxo.ne.loin(nsold)) goto 8

C If lmxo=-1 then the block with the basis-set information has not been 
C specified, the program will use an automatic basis set.
    
          if (lmxo.ne.-1) then 
          do 6 l=0,lmxo
            if(nzeta(l).ne.nzetin(l,nsold)) goto 8

            if(nzeta(l).gt.nzt) then 
             write(6,"('ATOM: ERROR: Bad dimensions for matrix Rco')")
             write(6,*) 'Nzt=',nzt,' but nzeta(',l,')=',nzeta(l)
             stop
            endif
           
   
            do izeta=1,nzeta(l)  
              if(rco(izeta,l).ne.rcoin(izeta,l,nsold)) goto 8   
              if(lambda(izeta,l).ne.lambdain(izeta,l,nsold)) goto 8   
            enddo 

  6       continue 
          endif 
          
c         is=nsold
c         no=nomax(is)
c         nkb=nkbmax(is)
c         lmxo=loold(is)
c         lmxkb=lkbold(is)
c         do l=0,lmxo
c           nzeta(l)=nzetold(l,is)
c           do izeta=1,nzeta(l)
c               rco(izeta,l)=rcoold(izeta,l,is)
c               lambda(izeta,l)=lambdaold(izeta,l,is)
c           enddo 
c         enddo 
c         do 7 i=1,no
c       
c             if (i .le. maxos) then
c                q(i)=qold(i,is)
c      else
c          write(6,'(/,a)')'ATOM: WARNING: Number of orbitals per atom'
c          write(6,'(a)') 'ATOM: WARNING: too small in calling routine'
c          write(6,'(a)') i,maxos
c             endif
c 7       continue 
c         goto 1000


          write(6,'(/,a,i2,a,i2)') 
     .     'ATOM: WARNING: Data for species', is+1, 
     .     ' are identical to those of the previous species',nsold

          goto 9 

  8       write(6,'(/,a)')
     .     'ATOM: WARNING: There are previous data for the same species'
          write(6,'(a)') 
     .     'ATOM: WARNING: Some of the arguments have been changed'
 
  9       continue
        
C*******ADDING A NEW SPECIES TO THE LIST********************************
          
          do ns=1,isold
              if(atm_label.eq.label_old(ns) ) then  
                  write(6,'(/,2a)') 
     .             'ATOM: WARNING: Two species with the same label  ', 
     .                atm_label
              endif
          enddo 

          is=isold+1
          isold=is
          ismax=is
          if(is.gt.nsmax) then
             write(6,"(a)")
     .           'ATOM: ERROR: Parameter nsmax must be increased'
             stop
          endif 
          if(is.gt.ntotsp) then 
              write(6,'(/,a,i3,2a,i3)')
     .     'ATOM: WARNING: Species ',is,' expected total number ',
     .      'of chemical species ',ntotsp

          endif 


          izold(is)=iz
          iztb(is)=iz
          loin(is)=lmxo
          lkbin(is)=lmxkb
          label_old(is)=atm_label

          
          if(lmxo.ne.-1) then  
           do l=0,lmxo

             nzetin(l,is)=nzeta(l)
           
             if(nzeta(l).gt.nzetmx) then
               write(6,"(a)")
     .           'ATOM: ERROR: Parameter nzetmx must be increased'
               stop
             endif

             do izeta=1,nzeta(l)
                rcoin(izeta,l,is)=rco(izeta,l) 
             enddo 

           enddo 
          endif 

          do i=1,lmx2*nzetmx
             qold(i,is)=0.0d0
          enddo


C********************      PI     *************************************
         
          pi=dacos(-1.0d0)

C**********************************************************************

        if ((iz.ne.0).and.(iz.ne.-100)) then 

C*********IF IZ NEGATIVE WE WILL HAVE FLOATING ORBITALS****************

       flting=dsign(1.0d0,dble(iz)) 
       iz=abs(iz)
       
         
       basistype=fdf_string('PAO.BasisType',basistype_default)
       if(basistype.eq.'SPLITGAUSS') basistype='splitgauss'
       if(basistype.eq.'NONODES')  basistype='nonodes'

       polorb=fdf_boolean('PAO.PolarizationOrbitals',polorb_default)
       if(polorb) then 
          numb_pol=fdf_integer('PAO.SplitPolarizationOrbitals',
     .         numb_pol_default)
       endif 

       nzeta_overflow=0
       norb_overflow=0
       nkb_overflow=0
       if (lmxo.eq.-1) then 
           if((basistype.eq.'nonodes')
     .         .or.(basistype.eq.'splitgauss')) then 
            write(6,*)
     .        'ATOM: BasisType equal to SPLIGAUSS or NONODES'
            write(6,*)
     .        'ATOM: needs the specification of the block'
            write(6,*) 'ATOM: PAO_basis_and_PS_lmax '
            stop
           endif 
           gen_basis='automatic'

C*IF THE BASIS SET IS AUTOMATIC, WE HAVE TO DEFINE SEVERAL PARAMETERS**
C******** FIRST WE READ DESIRED SIZE FOR THE BASIS SET*****************
    
           basis_size=fdf_string('PAO.BasisSize',basis_size_default)

           if(basis_size.eq.'STANDARD') basis_size='dzp'
           if(basis_size.eq.'standard') basis_size='dzp'
           if(basis_size.eq.'DZP')  basis_size='dzp'
 
           if(basis_size.eq.'DZ') basis_size='dz'

           if(basis_size.eq.'MINIMAL') basis_size='sz'
           if(basis_size.eq.'minimal')  basis_size='sz'
           if(basis_size.eq.'SZ')  basis_size='sz'

           
           if(basis_size.eq.'dzp') polorb=.true. 
               
           call lmxofz(iz,lmxPAO,latm)
           lmxPAOold(is)=lmxPAO

           if(polorb) then 
                 lmxo=lmxPAO+1
           else
                 lmxo=lmxPAO
           endif 
 
           if(lmxo+1.gt.lmaxd) then 
              write(6,"(2a,i3)")'ATOM: ERROR: Parameter lmaxd ',
     .           'must be increased to at least ', lmxo+1
              stop
           endif

              
           if((basis_size.eq.'dzp').
     .          or.(basis_size.eq.'dz')) then 
                two = 2
                if(two.gt.nzetmx) then
                  write(6,"(2a)")'ATOM: ERROR: Parameter nzetmx ',
     .              'must be increased to at least 2'
                  stop
                endif

                do l=0,lmxPAO
                  nzeta_overflow=max(nzeta_overflow,2)
                  norb_overflow=norb_overflow+2*(2*l+1)
                  nzeta(l)=2
                enddo 
           elseif(basis_size.eq.'sz') then 
                do l=0,lmxPAO
                  nzeta_overflow=max(nzeta_overflow,1)
                  norb_overflow=norb_overflow+(2*l+1)
                  nzeta(l)=1
                enddo 
           endif   
           if(polorb) then
            numb_pol=fdf_integer('PAO.SplitPolarizationOrbitals',
     .         numb_pol_default)
               if(basis_size.eq.'sz') numb_pol=1
               norb_overflow=norb_overflow+numb_pol*(2*lmxo+1)
               nzeta_overflow=max(nzeta_overflow,numb_pol)
           endif

           if(polorb) nzeta(lmxo)=numb_pol              
 
           lmxkb=lmxo+1

       else
           gen_basis='block'
           lmxPAO=lmxo
           lmxPAOold(is)=lmxo
           do l=0,lmxPAO
              nzeta_overflow=max(nzeta_overflow,nzeta(l))
              norb_overflow=norb_overflow+nzeta(l)*(2*l+1)
           enddo 
   
           if(polorb) then 
                lmxo=lmxPAO+1
                nzeta(lmxo)=numb_pol
                norb_overflow=norb_overflow+numb_pol*(2*lmxo+1)
                nzeta_overflow=max(nzeta_overflow,numb_pol)
           endif 
           
           
           if((lmxkb.lt.lmxo+1).and.(flting.gt.0.0d0)) then 
       write(6,'(/,2a,i2)') 'ATOM: WARNING: Maximum angular momentum ',
     .        'for KB projectors',lmxkb
       write(6,'(a,i2,a)') 'ATOM: WARNING: Increase that to ',lmxo+1,
     .        ' could be a good idea'
           endif 
       
       endif 
       
C**********If dimension are not enough, return to recompile*************

       nkb_overflow=(lmxkb+1)**2

       no=norb_overflow
       nkb=nkb_overflow
       nztout=nzeta_overflow

       if(nzt.lt.nzeta_overflow) then 
c        write(6,*) 'no',norb_overflow
c        write(6,*) 'nkb',nkb_overflow
c        write(6,*) 'nzeta',nzeta_overflow
c        write(6,*) 'lmxo',lmxo
c        write(6,*) 'lmxkb',lmxkb 
         return
       endif
       if(maxl.lt.lmxo) then
c        write(6,*) 'no',norb_overflow
c        write(6,*) 'nkb',nkb_overflow
c        write(6,*) 'nzeta',nzeta_overflow
c        write(6,*) 'lmxo',lmxo
c        write(6,*) 'lmxkb',lmxkb
         return
       endif

            
C**********************************************************************
             do l=0,lmxo
               nzetold(l,is)=nzeta(l)
               if(nzeta(l).gt.nzetmx) then
                 write(6,"(a)")
     .             'ATOM: ERROR: Parameter nzetmx must be increased'
                 stop
               endif
             enddo 
             lkbold(is)=lmxkb
             loold(is)=lmxo
C**********************************************************************




   
C*************READING INFORMATION ABOUT THE PSEUDOPOTENTIAL************


           fname = paste(atm_label,'.vps')
           inquire(file=fname, exist=found)
           if (.not.found) then
             write(6,'(/,2a,a20)') 'ATOM: WARNING: ',
     .         'Pseudopotential file not found: ', fname
             fname = paste(atm_label,'.psatom.data')
             write(6,'(2a)') 'ATOM: WARNING: Looking for ', fname
             inquire(file=fname, exist=found)
           endif
           if (.not.found) then
             write(6,'(/,2a,a20,/)') 'ATOM: ERROR: ', 
     .         'Pseudopotential file not found: ', fname
             stop 'ATOM: ERROR: Pseudopotential file not found'
           endif

           open(unit=1,file=fname,form='unformatted',status='unknown')

           read(1) namatm, icorr, irel, nicore,
     .     (method(i),i=1,6), (text(i),i=1,7),
     .     npotd, npotu, nr, b, a, zval
 
          
           write(6,'(/,a)')  'ATOM: Pseudopotential generation method:'
           write(6,'(7a)') 'ATOM: ',method(1),(method(i),i=3,6)



           if(irel.eq.'rel') then 
         write(6,'(/,2a)') 'ATOM: Pseudopotential generated from an ',
     .   'atomic-relativistic calculation'
         write(6,*) 'There are availables spin-orbital pseudopotentials'
         write(6,*) 'Spin-orbital interaction will not be included in ',
     .     'solid-state calculation'
           elseif (irel.eq.'isp') then
         write(6,'(/,2a)') 'ATOM: Pseudopotential generated from an ',
     .  'atomic spin-polarized calculation'
           endif 
 

           if (nicore.ne.'nc  ') then 
           write(6,'(/,a)') 
     .      'ATOM: Pseudopotential includes a core correction'

             if(nicore.eq.'pcec') then
               write(6,*) 'Pseudo-core for xc-correction'
             elseif(nicore.eq.'pche') then
               write(6,*) 'Pseudo-core for hartree and xc-correction'
             elseif(nicore.eq.'fcec') then
               write(6,*) 'Full-core for xc-correction'
             elseif(nicore.eq.'fche') then
               write(6,*) 'Full-core for hartree and xc-correction'
             endif

           endif 

           linp=max(lmxo,lmxkb)
           lmax=min(npotd-1,linp)


           if (lmax.gt.lmaxd) then 
             write(6,*) 'ATOM: Parameter lmaxd must be increased'
             write(6,*) 'to at least ',lmax
             stop
           endif

           if (lmax.lt.linp) then
             write(6,*) 'ATOM:  '
             write(6,*) 'You must generate a pseudopotential'
             write(6,*) 'for each L up to ',linp
             stop
           endif

           if(flting.gt.0.0d0) then 
            call lmxofz(iz,lval,latm)
            if(lval.gt.lmxPAO) then 
             write(6,*) 'ATOM:  '
             write(6,*) 'For Z=',IZ,' Lmax must be at least ',lval
             write(6,*) 'due to the ground state atomic configuration'
             stop
            endif
            do l=0,lmaxd
               ql(l)=0.0d0 
            enddo 
            call qvlofz(iz,ql)
           endif         
 
           
            
           nrval=nr+1
           if(nrval.gt.nrmax) then
              write(6,*) 'ATOM:'
              write(6,*) 'Nrmax must be increased to at least',nrval
              stop
           endif

           read(1) (rofi(ir),ir=2,nrval)
           rofi(1)=0.0d0

           do 20 ndown=1,lmax+1
               read(1) l,(vps(ir,l),ir=2,nrval)
               if(l.ne.ndown-1) then 
                  write(6,'(a)')
     . 'ATOM: Unexpected angular momentum  for pseudopotential'
                  write(6,'(a)')
     . 'ATOM: Pseudopotential should be ordered by incrising l'
               endif 
               do 19 ir=2,nrval
                    vps(ir,l)=vps(ir,l)/rofi(ir)
  19           continue 
  20       continue
           if(lmax+2.le.npotd)then 
           do ndown=lmax+2,npotd
              read(1) l
           enddo 
           endif 
           do 22 nup=1,npotu
              read(1) l
  22       continue 
           
 

C******* READ THE CORE CORRECTION CHARGE DENSITY *********************

          r2=rofi(2)/(rofi(3)-rofi(2))

          read(1) (chcore(ir),ir=2,nrval)
          chcore(1)=chcore(2)-(chcore(3)-chcore(2))*r2
 

C******** READ THE PSEUDO VALENCE DENSITY ****************************
   
          read(1) (rho(ir),ir=2,nrval)  
          rho(1)=rho(2)-(rho(3)-rho(2))*r2

          close(1)

C****CONSTRUCTION OF THE KLEINMAN-BYLANDER PROYECTOR FUNCTIONS**********
      
        nodd=mod(nrval,2)
        nrval=nrval-1+nodd

C********Calculate drdi*************************************************

          rpb=b
          ea=dexp(a)
          rofi(1)=0.0d0
          do 40 ir=1,nrval
             drdi(ir)=a*rpb
             s(ir)=dsqrt(a*rpb)
             rpb=rpb*ea
  40      continue


C****************PSEUDO-CORE RADIUS*************************************

         coretab(1,2,is)=0
         if(nicore.ne.'nc  ') then



           if(flting.gt.0.0d0) then
            
            coretab(1,2,is)=1
            nrcore=0
            eps=1.0d-6
            do ir=nrval,2,-1
              r=rofi(ir)
              r2=4.0d0*pi*r*r
              chc=chcore(ir)/r2

              if((chc.gt.eps).and.(nrcore.eq.0)) then 
                  nrcore=ir+1
                  Rcore=rofi(nrcore)
                  goto 30
              endif
            enddo
30         continue
 
           write(6,*) 'Pseudo-core radius Rcore=',Rcore

C****************TABLE WITH THE PSEUDO_CORE DATA******************** 

            delt=Rcore/(dble(ntbmax-1)+1.0d-20)
            coretab(1,1,is)=delt
            do itb=2,ntbmax
              r=delt*(itb-1)
              nr=nint(dlog(r/b+1.0d0)/a)+1
              nmin=max(1,nr-npoint)
              nmax=min(nrcore,nr+npoint)
              nn=nmax-nmin+1
              call ratint(rofi(nmin),chcore(nmin),nn,r,chc,dy)
              r2=4.0d0*pi*r*r 

              coretab(itb+1,1,is)=chc/r2
            enddo 
            coretab(2,1,is)=coretab(3,1,is) 

C*********TABLE WITH THE SECOND DERIVATIVE OF THE PSEUDO_CORE********
            yp1=0.0d0 
            ypn=1.0d50
            
          call spline(delt,coretab(2,1,is),ntbmax,     
     .      yp1,ypn,coretab(2,2,is),aux)   
                       
            elseif(flting.lt.0.0d0) then 
            
             do itb=1,ntbmax+1
                 coretab(itb,1,is)=0.0d0 
                 coretab(itb,2,is)=0.0d0 
             enddo 
        
            endif
 
           elseif(nicore.eq.'nc  ') then 

             do itb=1,ntbmax+1
               coretab(itb,1,is)=0.0d0 
               coretab(itb,2,is)=0.0d0 
             enddo 

           endif         
 
C********More information about pseudopot. generation*****************
          write(6,'(/,a)') 'ATOM: Valence configuration '//
     .                             '(pseudopotential generation):'
            do 10 l=0,lmax
              write(6,'(7x,a)')  text(2*l+1)
   10       continue
 
C***OBTAIN AN IONIC-PSEUDOPOTENTIAL IF CORE CORRECTION FOR HARTREE**** 
C*************************POTENTIAL***********************************

        if((nicore.eq.'pche').or.(nicore.eq.'fche')) then
            call vhrtre(chcore,ve,rofi,drdi,s,nrval,a)
            do l=0,lmax
              do ir=2,nrval
                vps(ir,l)=vps(ir,l)+ve(ir)
              enddo
            enddo
         endif



C*CALCULATION OF THE VALENCE SCREENING POTENTIAL FROM THE READED CHARGE 
C***********************DENSITY*****************************************
            

          call vhrtre(rho,ve,rofi,drdi,s,nrval,a)


C*************ADD THE EXCHANGE-CORRELATION POTENTIAL*******************

C*****Choosing the adecuate functional for the exchange-correlation****
   
       xcfunc = fdf_string('xc.functional','LDA')
       xcauth = fdf_string('xc.authors','PZ')


       write(6,'(/a)') 'ATOM: Exchange-correlation functional:'
       if(((xcauth.eq.'CA').or.(xcauth.eq.'PZ')).and.
     .    ((xcfunc.eq.'LDA').or.(xcfunc.eq.'LSD'))) then

          write(6,'(a)') 'ATOM: Ceperley-Alder'
          if(icorr.ne.'ca') then
           write(6,'(a)') 
     .      'ATOM: WARNING: Pseudopotential generated with'
           if(icorr.eq.'pw') write(6,'(a)') 
     .      'ATOM: WARNING: Perdew-Wang 1992 functional'
           if(icorr.eq.'pb') write(6,'(a)') 
     .  'ATOM: WARNING: GGA Perdew, Burke & Ernzerhof 1996 functional'
          endif

       elseif((xcauth.eq.'PW92').and.
     .    ((xcfunc.eq.'LDA').or.(xcfunc.eq.'LSD'))) then

         write(6,'(a)') 'ATOM: Perdew-Wang 1992'
         if(icorr.ne.'pw') then 
           write(6,'(a)') 
     .       'ATOM: WARNING: Pseudopotential generated with'
           if(icorr.eq.'ca') 
     .        write(6,'(a)') 'ATOM: WARNING: Ceperly-Alder functional'
           if(icorr.eq.'pb') write(6,'(a)') 
     .  'ATOM: WARNING: GGA Perdew, Burke & Ernzerhof 1996 functional'
         endif
       elseif((xcauth.eq.'PBE').and.(xcfunc.eq.'GGA')) then

         write(6,'(a)') 
     .     'ATOM: GGA Perdew, Burke & Ernzerhof 1996'
         if(icorr.ne.'pb') then
           write(6,'(a)') 
     .       'ATOM: WARNING: Pseudopotential generated with'
           if(icorr.eq.'ca') 
     .       write(6,'(a)') 'ATOM: WARNING: Ceperly-Alder functional'
           if(icorr.eq.'pw') 
     .       write(6,'(a)') 'ATOM: WARNING: Perdew-Wang 1992 functional'
         endif

       else  
       
          write(6,'(a)')
     .      'ATOM: Exchange-correlation functional not allowed'
          write(6,*) 'xc.functional= ',xcfunc
          write(6,*) 'xc.authors= ',xcauth
          stop

       endif

        if(irel.eq.'rel') irelt=1
        if(irel.ne.'rel') irelt=0

C*********************************************************************

 

          do ir=2,nrval
            r2=(rofi(ir))**2
            r2=4.0d0*pi*r2
            chgd=rho(ir)/r2
            if(nicore.ne.'nc  ')  chgd=chgd+chcore(ir)/r2
              auxrho(ir)=chgd
          enddo

          r2=rofi(2)/(rofi(3)-rofi(2))
          auxrho(1)=auxrho(2) -(auxrho(3)-auxrho(2))*r2
       

          call atomxc(xcfunc,xcauth,irelt,
     .             nrval,nrmax,rofi,1,auxrho,
     .             ex,ec,dx,dc,vxc)



          do ir=2,nrval
             ve(ir)=ve(ir)+vxc(ir)
          enddo

C******************************************************************



          if(flting.gt.0.0d0) then 

C***********WE HAVE TO CHOOSE A LOCAL POTENTIAL********************
C           In this case we will choose an 'gaussian' potential.
C******************************************************************

C******************************************************************
C All the pseudopotential will have a Kleinmann-Bylander projector
C
             lloc=lmxkb+1  
             loctab(is)=lloc        
C******************************************************************



C*******************************************************************
C   As a first step be check the maximum radius for the 
C   Kleinman-Bylander projectors with a standart choice
C   of the local potential
C
C*******Iterate over the possible local potentials*****************

          eps=1.0d-4
          rgauss=0.0d0
          nrgauss=0
          do l=0,lmxkb-1
               
             do ir=nrval,2,-1
                 dincv=dabs(vps(ir,l)-vps(ir,lmxkb))
                 if(dincv.gt.eps) goto 148
             enddo 
148          rgauss=max(rofi(ir),rgauss)
             nrgauss=max(ir,nrgauss)
          enddo
            dmax=0.0
            do ir=nrval,2,-1
              r=rofi(ir)
              dincv=dabs(vps(ir,0)*r+2.0d0*zval)
              if(dincv.gt.eps) goto 149
            enddo 
            
149         rgauss2=rofi(ir)
            

            if (rgauss2.gt.1.30d0*rgauss) then 

CIn this case the atom core is so big that we do not have an asymptotic 
Cof 2*Zval/r until Rgauss2 > Rc . To retain the same asymptotic 
Cbehaviour as in the pseudopotentials we modified the definition 
Cof the local potential

 
         write(6,'(/,a,f10.5)') 'ATOM: Estimated core radius ',
     .           rgauss2
         if (nicore.eq.'nc ') 
     .   write(6,'(/,2a)') 'ATOM: Include non-local core corrections',
     .   ' could be a good idea'
                nrgauss=nrgauss+3

                do ir=1,nrval
                   vlocal(ir)=vps(ir,0)*rofi(ir)
                enddo 

                   ir=nrgauss
                   dev=(vlocal(ir+1)-vlocal(ir-1))*0.5d0
                   dev2=(vlocal(ir+1)+vlocal(ir-1)-2.0d0*vlocal(ir))
                   dev3=(vlocal(ir+2)-2.0d0*vlocal(ir+1)
     .                 +2.0d0*vlocal(ir-1)-vlocal(ir-2))*0.5d0
                   dev3=(dev3-3.0d0*a*dev2+2.0d0*(a**2)*dev)
     .               /(drdi(ir)**3)
                   dev2=(dev2-a*dev)/(drdi(ir)**2)
                   dev=dev/drdi(ir)

C Local potential is Vloc(r)=v3*exp(v1*r^2+v2*r^3) inside Rgauss and equals the 
C all-electron atomic potential outside Rgauss
C We impose the continuity up to second            
c            if(ndevfit.eq.2) then               
               vlc=vlocal(nrgauss)
               r=rofi(nrgauss)

               var1=dev/vlc-1.0d0/r
               var2=dev2/vlc-2.0d0*var1/r -(var1)**2

               dm11=2.0d0*r
               dm12=3.0d0*r*r
               dm21=2.0d0
               dm22=6.0d0*r

               v1=(dm22*var1-dm12*var2)/(6.0d0*r*r)
               v2=(dm11*var2-dm21*var1)/(6.0d0*r*r)
               v3=vlc/(r*dexp((v1+v2*r)*r*r))


c            elseif(ndevfit.eq.3) then 

C*********************************************************************************
C  We can also construct a local potential Vloc(r)=v4*exp(v1*r^2+v2*r^3+v3*r^4),
C  this new coefficient allows us to impose the continuity of the potential up
C  to the third derivative.
C*******************************************************************************
 
c           vlc=vlocal(nrgauss)
c           r=rofi(nrgauss)
c           
c           var1=dev/vlc-1.0d0/r
c           var2=dev2/vlc-2.0d0*var1/r-(var1)**2
c           var3=dev3/vlc-3.0d0*var1*var2-(var1**3)
c    .                           -3.0d0*(var1**2+var2)/r

c           dm11=2.0d0*r
c           dm12=3.0d0*r*r
c           dm13=4.0d0*r*r*r
c           dm21=2.0d0
c           dm22=6.0d0*r
c           dm23=12.0d0*r*r
c           dm31=0.0d0
c           dm32=6.0d0
c           dm33=24.0d0*r

c           v1=((var1*dm22*dm33+var2*dm13*dm32+var3*dm12*dm23)
c    .   -(var3*dm22*dm13+var1*dm32*dm23+var2*dm12*dm33))/(48.0d0*r*r*r) 
c           v2=((var2*dm11*dm33+var3*dm21*dm13+var1*dm23*dm31)
c    .   -(var2*dm31*dm13+var3*dm23*dm11+var1*dm21*dm33))/(48.0d0*r*r*r)
c           v3=((var3*dm11*dm22+var2*dm12*dm31+var1*dm32*dm21)
c    .   -(var1*dm22*dm31+var3*dm21*dm12+var2*dm11*dm32))/(48.0d0*r*r*r)
c           v4=vlc/(r*dexp((v1+v2*r+v3*r*r)*r*r))
            
c         endif 


 
            open(unit=12,file=paste(atm_label,'.vlocal'),
     .         status='unknown')
             do ir=1,nrval
               r=rofi(ir)
               if(ir.le.nrgauss) then 
 
C************If second derivative fit**************************************
c               if(ndevfit.eq.2) then 
                    vlocal(ir)=v3*dexp((v1+v2*r)*r*r)
C**************************************************************************

C************If third derivative fit***************************************
c               elseif(ndevfit.eq.3) then 

c                   vlocal(ir)=v4*dexp((v1+v2*r+v3*r*r)*r*r)
c**************************************************************************
c               endif 

               else
                    vlocal(ir)=vps(ir,0)
               endif 
c              write(12,*) r,vlocal(ir)*r,(vps(ir,l)*r,l=0,lmxkb)
               write(12,*) r,vlocal(ir)*r

             enddo 
             close(12)



C Once we have the local potential we define the 'local-pseudopotential 
C charge' which help us to calculate the electrostatic interation 
C between the ions


          a2b4=0.25d0*a*a 
          qtot=0.d0 
          eps=1.0d-4
          do ir=1,nrval-1
             
            g2=vlocal(ir)*rofi(ir)
            if(dabs(g2+2.0d0*zval).lt.eps) goto 150

             if(ir.gt.nrgauss) then  

              if((ir.gt.2).and.(ir.lt.(nrval-1))) then 
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

              r=rofi(ir)
              cons=8.0d0*pi*r*drdi(ir)*s(ir)
              chcore(ir)=(-d2u)/cons
              qtot= qtot + 0.5d0*d2u*r/s(ir)




             else

C***************If second derivative fit************************************** 
c            if(ndevfit.eq.2)  then
              r=rofi(ir)

              g0=v3*dexp((v1+v2*r)*r**2)
              g1=(2.0d0*v1+3.0d0*v2*r)
              g2=2.0d0*v1+6.0d0*v2*r
              g3=(g2+g1*g1*r*r+2.0d0*g1)*g0
             
              cons=8.0d0*pi
              chcore(ir)= (-g3)/cons
              qtot= qtot  + 0.5d0*g3*r*r*drdi(ir)
 
C**************If third derivative fit****************************************
     
c             elseif(ndevfit.eq.3)  then

c             r=rofi(ir)
c          
c             g0=v4*dexp((v1+v2*r+v3*r*r)*r*r)     
c             g1=(2.0d0*v1+3.0d0*v2*r+4.0d0*v3*r*r)
c             g2=(2.0d0*v1+6.0d0*v2*r+12.0d0*v3*r*r)    
c             g3=(g2+g1*g1*r*r+2.0d0*g1)*g0   

c             cons=8.0d0*pi
c             chcore(ir)= -g3/cons
c             qtot= qtot  + 0.5d0*g3*r*r*drdi(ir)
c            endif 
C****************************************************************************



             endif



          enddo              


150       continue
          nchloc=ir          
          Rchloc=rofi(ir)

          write(6,'(2a,f10.5)') 'ATOM: Maximum radius for' ,
     .      ' local-pseudopot. charge ',Rchloc

          do ir=1,nchloc
            r=rofi(ir)
             chc=zval*chcore(ir)/qtot
             chcore(ir)=chc
          enddo 

          delt=Rchloc/(dble(ntbmax-1)+1.0d-20)
  
          chloctab(1,1,is)=delt
          chloctab(1,2,is)=1.0d0
          do itb=1,ntbmax
             r=delt*(itb-1)
             nr=nint(dlog(r/b+1.0d0)/a)+1
             nmin=max(1,nr-npoint)
             nmax=min(nchloc,nr+npoint)
             nn=nmax-nmin+1
             call ratint(rofi(nmin),chcore(nmin),nn,r,chc,dy)

             chloctab(itb+1,1,is)=chc
          enddo

          chloctab(2,1,is)=chloctab(3,1,is)
 


C****TABLE WITH THE SECOND DERIVATIVE OF THE LOCAL-PSEUDOTENTIAL********
C***********************CHARGE DENSITY********************************** 

         yp1=1.0d50
         ypn=1.0d50

         call spline(delt,chloctab(2,1,is),ntbmax,
     .      yp1,ypn,chloctab(2,2,is),aux)
        
C*****CALCULATION OF THE ELECTROSTATIC CORRECTION***********************

         do is2=is,1,-1
            rmax=chloctab(1,1,is2)*(ntbmax-1)+Rchloc+0.2d0
            indx=((is-1)*is)/2+is2
            corrtab(1,1,indx)=rmax/(ntbmax-1)
            corrtab(1,2,indx)=1.0d0
            call choverlp(is,is2,rmax,corrtab(2,1,indx),
     .       corrtab(2,2,indx),aux)
     
            if(dabs(rmax).lt.1.0d-8) corrtab(1,2,indx)=0.0d0

         enddo 
         
         



          



           else 






C******************************************************************
C           Local-potential size parameter 'rgauss'
C   We choose as a smooth pseudopotential the one generated 
C   by a 'Vanderbilt-function' charge distribution. We have to select 
C   the size of this distribution somehow.
C   'Vandebilt-functions' are of the form :
C    p(r)=N*exp(-(sinh(van*r)/sinh(van)**2)
C    when van---> 0 we will obtain a 'gaussian'
C    when van---> Inf. we will obtain a step function
C    Some test has revealed that the best election to achieve 
C    a good converge in real and reciprocal space is b in the 
C    range 0.5-1.0 .
C******************************************************************

C  So, the 'gaussian' charge distribution 
C  must go to zero at a distance 'rgauss'.


         open(unit=12,file=paste(atm_label,'.vlocal'),
     .       status='unknown')


c       if(loctype.eq.'new') then          


C***********************************************************************
C     We take a 'Vanderbilt-function' as local potential
C     van=1.0d0 all the parameter have optimized for this value 
C***********************************************************************
                van=1.0d0
                cutoff1=3.63d0
                cutoff2=5.48d0
C**********************99% of charge inside Rgauss**********************
c               factor=1.627d0
C***********************************************************************

C**********************99.9% of charge inside Rgauss********************
                factor=1.815
C***********************************************************************
              
C*********** Scaling factor for local-pseudopot. charge*****************
                  alp=factor/rgauss
C***********************************************************************


       write(6,'(/,a,f10.3,a)')
     . 'ATOM: 99.0% of the norm of Vloc inside ',(alp*cutoff1)**2,' Ry'
       write(6,'(a,f10.3,a)')
     . 'ATOM: 99.9% of the norm of Vloc inside ',(alp*cutoff2)**2,' Ry'


c        elseif(loctype.eq.'old') then 
                 
c                van=0.00001d0 
c                rgauss=0.80d0
c                factor=2.0d0


C*********** Scaling factor for local-pseudopot. charge*****************
c               alp=factor/rgauss  
C***********************************************************************
                
c        endif 

          qtot=0.0d0 
          do ir=1,nrval
             r=rofi(ir) 
             gexp=dsinh(van*alp*r)/dsinh(van)
             gexp=gexp*gexp
             rhor=dexp(-gexp)
             chcore(ir)=(-4.0d0)*pi*rhor*r*r
             qtot=qtot+rhor*drdi(ir)*r*r
          enddo
          qtot=4.0d0*pi*qtot 
          eps=1.0d-4
          nchloc=0 
          do ir=nrval,1,-1
             chc=zval*chcore(ir)/qtot
             chcore(ir)=chc   
             if((dabs(chc).gt.eps).and.(nchloc.eq.0)) then    
                 nchloc=ir+1
                 Rchloc=rofi(nchloc)
             endif
          enddo 
          
          write(6,'(/,2a,f10.5)') 'ATOM: Maximum radius for' ,
     .      ' local-pseudopot. charge ',Rchloc

          delt=Rchloc/(dble(ntbmax-1)+1.0d-20)
  
          chloctab(1,1,is)=delt
          chloctab(1,2,is)=1.0d0
          do itb=2,ntbmax
             r=delt*(itb-1)
             nr=nint(dlog(r/b+1.0d0)/a)+1
             nmin=max(1,nr-npoint)
             nmax=min(nchloc,nr+npoint)
             nn=nmax-nmin+1
             call ratint(rofi(nmin),chcore(nmin),nn,r,chc,dy)
             r2=4.0d0*pi*r*r

             chloctab(itb+1,1,is)=chc/r2
          enddo

          chloctab(2,1,is)=chloctab(3,1,is)

C****TABLE WITH THE SECOND DERIVATIVE OF THE LOCAL-PSEUDOTENTIAL********
C***********************CHARGE DENSITY********************************** 

         yp1=1.0d50
         ypn=1.0d50

         call spline(delt,chloctab(2,1,is),ntbmax,
     .      yp1,ypn,chloctab(2,2,is),aux)
        
C*****CALCULATION OF THE ELECTROSTATIC CORRECTION***********************

         do is2=is,1,-1
            rmax=chloctab(1,1,is2)*(ntbmax-1)+Rchloc+0.2d0
            indx=((is-1)*is)/2+is2
            corrtab(1,1,indx)=rmax/(ntbmax-1)
            corrtab(1,2,indx)=1.0d0
            call choverlp(is,is2,rmax,corrtab(2,1,indx),
     .       corrtab(2,2,indx),aux)

            if(dabs(rmax).lt.1.0d-8) corrtab(1,2,indx)=0.0d0
         enddo 
         


 
C***** Calculation of the local potential generated by the ************ 
C********************* charge distribution*****************************


          call vhrtre(chcore,vlocal,rofi,drdi,s,nrval,a)
  
          
          do ir=2,nrval 
             r=rofi(ir)
             if (r.gt.1.1d0*Rchloc) then
                 vlocal(ir)=(-2.0d0)*zval/rofi(ir)
             endif

             write(12,*) r,vlocal(ir)*rofi(ir)
c    .          ,(vps(ir,l)*r,l=0,lmxkb)
          enddo 
          close(12) 

          r2=rofi(2)/(rofi(3)-rofi(2))
          vlocal(1)=vlocal(2)-(vlocal(3)-vlocal(2))*r2
   



          endif 



C********************************************************************

         elseif(flting.lt.0.0d0) then

            do itb=1,ntbmax+1
               chloctab(itb,1,is)=0.0d0
               chloctab(itb,2,is)=0.0d0
            enddo 
            do is2=is,1,-1
              indx=((is-1)*is)/2+is2
              do itb=1,ntbmax+1
                  corrtab(itb,1,indx)=0.0d0
                  corrtab(itb,2,indx)=0.0d0
              enddo 
            enddo

         endif


        
 


C********ARRAY S FOR THE SCHRODINGER EQ. INTEGRATION*****************
         do 42  ir=2,nrval
             s(ir)=drdi(ir)*drdi(ir)
  42     continue
         s(1)=s(2)


C***AND ONLY IF IT IS NEEDED****************************************** 
c          if(lmxkb.gt.0)then

C******CALCULATION OF THE KLEINMAN-BYLANDER PROYECTOR FUNCTIONS*******


         a2b4=a*a*0.25d0
         do 45 l=0,lmxkb

            do 43 ir=2,nrval
               r2=(rofi(ir))**2
               vtot=vps(ir,l)+ve(ir)+dble(l*(l+1))/r2
               h(ir)=vtot*s(ir)+a2b4
  43        continue
            h(1)=h(2)


            nnode=1
            nprin=l+1
            e=-((zval/dble(nprin))**2)
            z=zval
            dr=-1.0d5
            rmax=rofi(nrval)
            call egofv(h,s,nrval,e,g,y,l,z,a,b,rmax,
     .        nprin,nnode,dr)


            eigen(l)=e
            
            if(flting.gt.0.0d0) then 
              do 44 ir=2,nrval
                 r=rofi(ir)
                 f=g(ir)
                 dsq=dsqrt(drdi(ir))
                 f=f*dsq
                 rphi(ir,l)=f
   44         continue
              rphi(1,l)=rphi(2,l)
            endif   
   45     continue

C***ONLY CALCULATE THE PROJECTORS IF IT IS A REAL ATOM****************
         if(flting.gt.0.0d0) then
         
C       checking normalization of the calculated wave functions

          do 47 l=0,lmxkb
            dnrm=0.0d0
            do 46 ir=2,nrval
               phi=rphi(ir,l)
               dnrm=dnrm+phi*phi*drdi(ir)
   46       continue
               dnrm=dsqrt(dnrm)
               ddnrm=dabs(dnrm-1.0d0)

            if (ddnrm.gt.1.0d-5) then
               write(6,'(/,a,i2,a)') 
     .    'ATOM: WARNING: Eigenstate for l=',l,'is not normalized'
               write(6,'(a,f12.6)') 'ATOM: WARNING: Norm=',dnrm
            endif
   47     continue


C******************GHOST ANALYSIS**********************************


C*******Iterate over the possible local potentials*****************
c         do 53 loc=lmxkb,0,-1
           ighost=0

c          write(6,'(/,a,i2,a)')'ATOM: Testing pseudopot. for l=',loc,
c    .       ' as local'
            
C***CALCULATE EIGENVALUES OF LOCAL POTENTIAL FOR GHOST ANALYSIS****** 
   
          do 50 l=0,lmxkb
             do 48 ir=2,nrval
                r2=rofi(ir)**2
c               vtot=vps(ir,loc)+ve(ir)+dble(l*(l+1))/r2

C* ATTENTION , 'Ve' is the screenig potential generated from valence*  
C* pseudo-charge given by the pseudopotential generation program ****
C********************************************************************
                vtot=vlocal(ir)+ve(ir)+dble(l*(l+1))/r2
C********************************************************************
                h(ir)=vtot*s(ir)+a2b4
   48       continue
            do 49 nnode=1,2
             nprin=l+1
             e=-((zval/dble(nprin))**2)
             z=zval
             dr=-1.0d5
             rmax=rofi(nrval)
             call egofv(h,s,nrval,e,g,y,l,z,a,b,rmax,
     .        nprin,nnode,dr)
               elocal(nnode,l)=e
   49       continue
c        write(6,*) 'ATOM: Ground state vlocal for L=',l,elocal(1,l)  
c        write(6,*) 'ATOM: First excited state for L=',l,elocal(2,l) 
   50     continue 
           
C****************CALCULATE KB-COSINE**********************************

           do 52 l=0,lmxkb

             dnrm=0.0d0
             avgv=0.0d0
             do 51 ir=2,nrval
               r=rofi(ir)
c              vl=(vps(ir,l)-vps(ir,loc))
               vl=(vps(ir,l)-vlocal(ir))
               phi=rphi(ir,l)
               vphi=vl*phi
               dnrm=dnrm+vphi*vphi*drdi(ir)
               avgv=avgv+vphi*phi*drdi(ir)
  51         continue


             dkbcos(l)=dnrm/(avgv+1.0d-20) 
             dknrm(l)=1.0d0/(dsqrt(dnrm)+1.0d-20)

C***************GHOST ANALYSIS*****************************************


          if(dkbcos(l).gt.0.0d0) then

              if(eigen(l).gt.elocal(2,l)) then
                 write(6,"(a,i3)")
     .            'ATOM: WARNING: Ghost state for L =', l
                 ighost=1
              else
                 write(6,'(a,i3)') 'ATOM: No ghost state for L =',l
              endif

           elseif(dkbcos(l).lt.0d0) then

              if(eigen(l).gt.elocal(1,l)) then
                 write(6,"(a,i3)")
     .            'ATOM: WARNING: Ghost state for L =', l
                 ighost=1
              else
                 write(6,'(a,i3)') 'ATOM: No ghost state for L =',l
              endif

           elseif(dkbcos(l).eq.0.0d0) then

               write(6,"('ATOM: vps = vlocal, no ghost for L =',i3)") l

           endif
           if(ighost.eq.1) then 
c            write(6,"('ATOM: WARNING: Ghost state for L =',i3)") L
c            write(6,"('Trying with the next possible local potential')")
             goto 53

           endif
 
 52      continue
         goto 54 
C***********************************************************************
           
       
 53       continue



 54       if (ighost.eq.1) then
            write(6,"(2a)")'ATOM: WARNING: ',
     .            'No pseudopotential free of ghost states'
            write(6,"(2a)")'ATOM: WARNING: ',
     .            'Some parameter must be changed in the '
            write(6,"(2a)")'ATOM: WARNING: ',
     .            'pseudopotential generation procedure.'
            stop
          else
c           lloc=loc
c           loctab(is)=loc
          endif



C*******DEFINE THE CUT-OFF RADII************************************
C Warning these radii should be quite short, if it is not the case 
C something is probably wrong in this part of the program.
C It will display a warning if Rc>4.5 a.u.or Rc < 0.5a.u.!!!!!!!!!!!!
            
            eps=1.0d-6

            do 55 l=0,lmxkb
              if(l.eq.lloc) goto 55
              do 56 ir=nrval,2,-1
                 phi=(rphi(ir,l)/rofi(ir))*dknrm(l)
c                dincv=dabs(vps(ir,l)-vps(ir,lloc))*phi
                 dincv=dabs(vps(ir,l)-vlocal(ir))*phi
                 if(dincv.gt.eps) goto 57
56            continue
57            rc(l)=rofi(ir+1) 
              rctb(l,is)=rofi(ir+1)

              if(rc(l).lt.0.5d0) then
                write(6,"('ATOM: WARNING: Rc(',i2,')=',f12.4)")l,rc(l)
                write(6,"(2a)") 'ATOM: WARNING: ',
     .            'Check ATOM: look for the sentence:'
                write(6,"(2a)") 'ATOM: WARNING: ',
     .            'DEFINE THE CUT-OFF RADII'       
              elseif(rc(l).gt.4.5d0) then
               write(6,"('ATOM: WARNING: Rc(',i2,')=',f12.4)")l,rc(l)
                write(6,"(2a)") 'ATOM: WARNING: ',
     .            'Check ATOM: look for the sentence:'
                write(6,"(2a)") 'ATOM: WARNING: ',
     .            'DEFINE THE CUT-OFF RADII'       
               write(6,"(2a)") 'ATOM: WARNING: ',
     .            'Increasing the tolerance parameter eps'
               write(6,"(a)") 'ATOM: WARNING: might be a good idea'
              endif
55           continue
c            rc(lloc)=0.0d0
c            rctb(lloc,is)=0.0d0 
             

C*******KLEINMAN-BYLANDER PROJECTION FUNCTIONS******************      
   
          indx=0
          do 58 l=0,lmxkb
             if(l.eq.lloc) goto 58 
             nrc=nint(dlog(rc(l)/b+1.0d0)/a)+1
             indx=indx+1
             do 59 ir=2,nrc
               r=rofi(ir)
c              vl=(vps(ir,l)-vps(ir,lloc))
               vl=(vps(ir,l)-vlocal(ir))
               phi=rphi(ir,l)/r
               vphi=vl*phi*dknrm(l)
               h(ir)=vphi
               h(ir)=h(ir)/r**l
  59         continue

             h(1)= ( h(2)*rofi(3)**2 - h(3)*rofi(2)**2 ) /
     .             (      rofi(3)**2 -      rofi(2)**2 )


C**********INTERPOLATION TO GENERATE TABLES WITH KB PROJECTORS*******

              
             delt=rc(l)/(dble(ntbmax-1)+1.0d-20)
             table(1,-indx,is)=delt
          
             table(2,-indx,is)=dkbcos(l)
 
             do 60 itb=1,ntbmax
                r=delt*(itb-1)
                nr=nint(dlog(r/b+1.0d0)/a)+1
                nmin=max(1,nr-npoint)
                nmax=min(nrc,nr+npoint)
                nn=nmax-nmin+1
                call ratint(rofi(nmin),h(nmin),nn,r,vphi,dy)
                table(itb+2,-indx,is)=vphi
  60         continue        
             
  58          continue 

           write(6,'(/,a)')'ATOM: Kleinman-Bylander projectors: '
           do l=0,lmxkb
             if(l.ne.lloc) then 
             write(6,'(3x,a,i2,3(3x,a,f10.6))')
     .        'l=',l, 'rc=',rc(l), 'el=',eigen(l), 'kbcos=',dkbcos(l)
             else
              write(6,'(3x,a,i2,a,f10.6)') 'l=',l,
     .         '   Local pseudopotential, el=',eigen(l)
             endif
           enddo
         
             
C***********TOTAL NUMBER OF KLEINMAN-BYLANDER PROJECTORS************** 

c          nkb=(lmxkb+1)*(lmxkb+1)-(2*lloc+1)
           nkb=(lmxkb+1)*(lmxkb+1)

c         elseif(lmxkb.eq.0) then 
c          lloc=0
c          loctab(is)=0
c          nkb=0
c         endif
          
         elseif(flting.lt.0.0d0) then

          nkb=0
          loctab(is)=0

         endif


C********************************************************************** 
         nkbmax(is)=nkb
C**********************************************************************



         write(6,'(a,i3)')
     .     'ATOM: Total number of KB projectors:',nkb
 
C***CONSTRUCTION OF THE SANKEY-TYPE ORBITALS FOR EACH ANGULAR MOMENTUM**
         
         basistype=fdf_string('PAO.BasisType',basistype_default) 


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
     .        'ATOM: Incorrect basis-type option specified,',
     .        ' active options are:',
     .        'NODES','SPLIT','USER','SPLITGAUSS','NONODES'
              stop
         endif 
           
         if(basistype.ne.'user') then  
           write(6,'(/,a)') 'ATOM: Sankey-type orbitals:' 
           nzcontr=0
           do l=0,lmxo
              if(nzeta(l).gt.1) nzcontr=1
           enddo  
           if( nzcontr.eq.1) then 
               write(6,'(2a)')
     .        'ATOM: Selected multiple-zeta basis: ',basistype
           endif
         endif 

C***************Open a file to write the basis*************************
           fname= paste(atm_label,'.PAO.basis')
           open(unit=2,file=fname,status='unknown',form='formatted')
           write(2,*) lmxo
C***************User basis ********************************************
         if(basistype.eq.'user') then

           fname = paste(atm_label,'.user.basis')
           inquire(file=fname, exist=found)
           if (.not.found) then
              write(6,'(/,3a,/)') 'ATOM: ', symbol(iz),
     .         ' user-basis file ', fname, ' not found.'
              stop 'ATOM: user-basis file not found'
           endif


           write(6,'(/,2a,/a)') 
     .     'ATOM: User-basis. The basis orbitals will be readed ', 
     .     'from the file: ',fname

           open(unit=1,file=fname,
     .       status='old',form='formatted')
 
            read(1,*) lmxfile
            if(lmxfile.lt.lmxo) then 
              write(6,'(/,a,a,/,2a,i2,/,a,i2)')
     .       'ATOM: ERROR: ',
     .       'maximum angular momentum of the orbitals in file',
     .        paste(symbol(iz),'.user.basis'), 'is ',lmxfile,
     .        'it should be at least', lmxo
              stop
            endif 


         endif  


C**********************************************************************
         if(flting.lt.0.0d0) write(6,'(2a)')
     .       'ATOM: Floating orbitals corresponding to ',symbol(iz)




C***********************************************************************
         norb=0
C***********************************************************************



C***READING THE ENERGY-SHIFT TO DEFINE THE CUT-OFF RADIUS OF ORBITALS***

         eshift=fdf_physical('PAO.EnergyShift',eshift_default,'Ry')
        
C***********************************************************************


C***READING SPLNORM TO GENERATE THE SPLIT IF Rmatch IS ZERO IN INPUT****
     
         splnorm=fdf_double('PAO.SplitNorm',splnorm_default)

C*********************************************************************** 



         do  100 l=0,lmxo

            write(6,'(/A,I2)') 
     .       'ATOM: Orbitals with angular momentum L=',l


C**************************  write the basis  **************************
              write(2,*) l,nzeta(l)
C***********************************************************************



C***********************************************************************
          if(basistype.eq.'user') then 
 
               read(1,*) luser,nzetl

            if(luser.ne.l) then
               write(6,'(/,2a,/,2(a,i2),/,a)')
     .      'ATOM: Reading user-basis orbitals from file:',fname,
     .      'expected l=',l,' readed l=',luser,
     .      'check order!!!!!!!!!!!!!!!!!!!!!'
               stop
            endif
              
             if (nzetl.lt.nzeta(l)) then 
                  write(6,'(/,2a,i2,a,/,2a,i2,/,a,i2)')
     .           'ATOM: ERROR: ',
     .           'number of orbitals with l=', l,' in file',
     .            paste(symbol(iz),'.user.basis'), ' is ',nzetl,
     .            'it should be at least', nzeta(l)
                stop
             endif 

          endif 
C*********************************************************************** 



C********TABLE WITH THE NUMBER OF ORBITALS FOR EACH ANG. MOMENTUM*******
            nzettb(l,is)=nzeta(l)
C***********************************************************************


C********Calculating hamiltonian for solving Schrodinger eqn.***********       
          if(basistype.ne.'user') then      
            a2b4=a*a*0.25d0
        
            do 70 ir=2,nrval
              vtot=vps(ir,l)+ve(ir)+dble(l*(l+1))/(rofi(ir)**2)
              h(ir)=vtot*s(ir)+a2b4
  70        continue
            h(1)=h(2)
          endif 
C************************************************************************




C*********Are there enough orbitals ??????******************************
           if(nzeta(l).lt.1) then
            write(6,"(2a,i3)")'ATOM: WARNING: ',
     .        'No orbital has been calculated for l =', l
            write(6,"(2a,i2,a,i3)")'ATOM: WARNING: ',
     .        'Nzeta(',l,') =',nzeta(l)
            if((ql(l).ne.0.0d0).and.(flting.gt.0.0d0)) then
              write(6,"(2a,i3,a)")'ATOM: ERROR: ',
     .          'Orbital with l =',l,' must exist'
              write(6,"(2a)")'ATOM: ERROR: ',
     .          'due to the atomic ground state configuration'
              stop
            endif
           endif
C***********************************************************************



C****ESTIMATED ORBITAL CUT-OFF RADIUS FROM ENERGY SHIFT***************      
C****THIS IS ONLY TRUE BEFORE COMPRESSION IF THE COMPRESSION FACTOR***
C****IS DIFFERENT FROM ONE THE ENERGY SHIFT FOR THE ORBITAL WILL BE***
C***************************DIFFERENT*********************************   

          if(dabs(eshift).gt.1.0d-5) then
             el=eigen(l)+eshift
             call rc_vs_e(a,b,rofi,vps(1,l),ve,nrval,l,el,rnodo)
c            write(6,'(/,A,I2,A,f10.6,A)')
c    .           'L= ',l,' Energy shift= ',eshift,' Ry'
c            write(6,'(A,f10.6)') 
c    .            'Aprox. Rc for this energy shift= ',rnodo
          else
             rnodo=rofi(nrval-2) 
          endif 
           
           
C*******************************************************************
           
C*We decide if an orbital is calculated as a PAO, or perturbatively*
C*******from the orbital with lower angular momentum****************
C*******************************************************************
C***To generate a double-zeta orbital from a polarization orbital
C** the only active option is 'split'****************************
            ipol='no'
            if(l.gt.lmxPAO) then 
                 ipol='yes'
                 basistype='split'
            endif 


                          
            do 80 izeta=1,nzeta(l)

C**********With spligauss option, compression factor must be taken****
C**********as the gaussian exponent***********************************
              if((basistype.eq.'splitgauss').and.(izeta.gt.1)) then
                  gexp=dabs(lambda(izeta,l))
                  gexp=1.0d0/(gexp**2)
                  lambda(izeta,l)=1.0d0
              endif
C*********************************************************************


C*****IF THE COMPRESSION FACTOR IS NEGATIVE OR ZERO THE ORBITALS ARE***
C********************UNTOUCHED******************************************
             if(lambda(izeta,l).le.0.0d0) lambda(izeta,l)=1.0d0
C***********************************************************************


C*********Comparing requested cut-off radius with that define using*****
C*******the energy-shift parameter**************************************
            if ((izeta.eq.1).and.(ipol.eq.'no')) then
               if(rco(1,l).le.1.0d-5) then
                 if(basistype.ne.'user') then
                  write(6,'(/,A,/,A,f10.6,A)')
     .             'ATOM: PAO cut-off radius determinated from an',
     .             'ATOM: energy shift=',eshift,' Ry'
                  endif

                   rco(1,l) = rnodo
               endif 
            endif  
             
            if((basistype.eq.'nodes').and.(rco(izeta,l).le.1.0d-5))then
               rco(izeta,l)=rco(1,l)
             endif
C************************************************************************


C********Final cut-off radius for basis orbitals*************************
              if((basistype.eq.'split').and.(izeta.gt.1)) then
                  Rsp2=rco(izeta,l)/lambda(1,l)
                  nsp2=nint(dlog(Rsp2/b+1.0d0)/a)+1
                  Rsp2=b*(dexp(a*(nsp2-1))-1.0d0)
C****COMPRESSION FACTOR IS ONLY ACTIVE FOR THE INITIAL PAO WHEN USING**** 
C**** SPLIT OPTION FOR THE GENERATION OF THE BASIS SET*******************
                 lambda(izeta,l)=lambda(1,l)
C************************************************************************ 
                if(Rsp2.gt.rco(1,l)) then 
                  write(6,'(/,2A)') 'ATOM: SPLIT OPTION FOR BASIS SET ',
     .             'Rc FOR DOUBLE-Z, TRIPLE-Z,... SHOULD BE SMALLER '
                  write(6,'(A)') 'THAN THAT OF THE INITIAL PAO !!!!!'
                  STOP
                endif 
                rco(izeta,l)=Rsp2
              elseif(ipol.eq.'no') then 
                  rco(izeta,l)=rco(izeta,l)/lambda(izeta,l)
                  nrc=nint(dlog(rco(izeta,l)/b+1.0d0)/a)+1
                  nodd=mod(nrc,2)
                  if(nodd.eq.0) then 
                     nrc=nrc+1
                  endif 
                  rmax=b*(dexp(a*(nrc-1))-1.0d0) 
                  rco(izeta,l)=rmax*lambda(izeta,l)
              elseif(ipol.eq.'yes') then 
                      lambda(1,l)=lambda(1,l-1)
                      rpol=rco(1,l-1)/lambda(1,l)
                      nrc=nint(dlog(rpol/b+1.0d0)/a)+1
                      nodd=mod(nrc,2)
                      if(nodd.eq.0) then
                          nrc=nrc+1
                      endif
                      rpol=b*(dexp(a*(nrc-1))-1.0d0)
                      rco(1,l)=rpol*lambda(izeta,l)
              endif 

C****************************************************************


C********************************************************************
          norb=norb+1
C********************************************************************


C****************Calculation of the basis functions**************
          if((basistype.eq.'splitgauss').and.(izeta.gt.1)) then 
             
             fac=1.0
             do i=0,l
                fac=(2*i+1)*fac
             enddo

            cons=sqrt(pi)*fac/(2.0d0**(l+2))
            cons=cons/((2.0d0*gexp)**(l+1.5d0))
            cons=1.0d0/sqrt(cons)

            dnrm=0.0d0
            do ir=1,nrc
              r=rofi(ir)
              f=cons*dexp((-gexp)*r**2)
              dnrm=dnrm+drdi(ir)*(f*r**(l+1))**2
              g(ir)=f
            enddo 

          elseif((basistype.eq.'split').and.(izeta.gt.1)) then

             if(Rsp2.gt.1.0d-5) then 

               frsp=rphi(nsp2,l)/Rsp2
               dfrsp=0.5d0*(rphi(nsp2+1,l)/rofi(nsp2+1)
     .             -rphi(nsp2-1,l)/rofi(nsp2-1))
               dfrsp=dfrsp/drdi(nsp2)

C*****************gaussian split************************************
c           dfrsp=dfrsp/frsp
c           gexp=0.5d0*(dble(l)-Rsp2*dfrsp)/(Rsp2)**2
c           cons=frsp/((Rsp2**l)*dexp(-gexp*(Rsp2)**2))
c           call gauss(a,b,nrc,rphi(1,l),l,splnorm,cons3,gexp3,nsp3)
C*******************************************************************


C**********************parabolic split******************************
            cons1= 0.5d0*(dfrsp*Rsp2-l*frsp)/(Rsp2**(l+2))
            cons2= frsp/(Rsp2**l)-cons1*Rsp2**2
            call nrmpal(cons1,cons2,Rsp2,l,rnp)
            spln=1.0d0-rnrm(nsp2,l)+rnp
C*******************************************************************

                do i=1,izeta-1
                 if(dabs(rco(izeta,l)-rco(i,l)).lt.1.0d-5) then 
                   write(6,'(/,A,I2,A,I2,A,I2)') 
     .            'ATOM: WARNING: Split-orbital with zeta=',izeta,
     .            ' and zeta=',i,' are identicals for l=',l
                 endif  
                enddo 
            else
            rmax=rco(1,l)/lambda(1,l)
            nrc1=nint(dlog(rmax/b+1.0d0)/a)+1
            spln=splnorm
            if(izeta.gt.2) then 
              spln=spln/(2.0d0*(izeta-2) )
            endif 

            call parabola(a,b,nrc1,rphi(1,l),rnrm(1,l)
     .                  ,l,spln,cons1,cons2,nsp2)
             

C***Cut-off radius for the split orbital with a desired norm******
       nrc=nsp2
       rco(izeta,l)=b*(dexp(a*(nsp2-1))-1.0d0)*lambda(izeta,l)
C*****************************************************************


             do i=1,izeta-1
                if(dabs(rco(izeta,l)-rco(i,l)).lt.1.0d-5) then
                   write(6,'(/,A,I2,A,I2,A,I2)')
     .            'ATOM: WARNING: Split-orbital with zeta=',izeta,
     .            ' and zeta=',i,' are identicals for l=',l
                endif 
             enddo 
 
            endif 

            dnrm=0.0d0
            do ir=2,nsp2-1
              r=rofi(ir)
C***********************parabolic split****************************
              f=-(cons1*r**2+cons2)+rphi(ir,l)/(r)**(l+1)
C******************************************************************

C**********************gaussian split******************************
c             f=-cons*dexp(-gexp*r**2)+rphi(ir,l)/(r)**(l+1)
C******************************************************************

              dnrm=dnrm+drdi(ir)*(f*r**(l+1))**2
              g(ir)=f
            enddo
            g(1)=g(2)
            g(nsp2)=0.0d0

 

          elseif(basistype.ne.'user') then 

            if(basistype.eq.'nodes') then
             nnode=izeta
            elseif(basistype.eq.'nonodes') then 
             nnode=1
            elseif((basistype.eq.'splitgauss').and.(izeta.eq.1)) then
             nnode=1 
            elseif((basistype.eq.'split').and.(izeta.eq.1)) then
             nnode=1
            endif

            if(ipol.eq.'no') then 

            z=zval
            dr=-1.0d6
            nprin=l+izeta
            eorb=-((zval/dble(nprin))**2)

            call egofv(h,s,nrc,eorb,g,y,l,z,a,b,rmax,nprin,
     .       nnode,dr)
           
  
            dnrm=0.0d0
            do 90 ir=2,nrc
               r=rofi(ir)
               f=g(ir)
               dsq=dsqrt(drdi(ir))
               f=f*dsq
               dnrm=dnrm+drdi(ir)*f*f
               g(ir)=f/r
               if(izeta.eq.1) rphi(ir,l)=f
               if(izeta.eq.1) rnrm(ir,l)=dnrm
               g(ir)=g(ir)/r**l

 90         continue

            g(1)= ( g(2)*rofi(3)**2 - g(3)*rofi(2)**2 ) /
     .            (      rofi(3)**2 -      rofi(2)**2 )

           elseif(ipol.eq.'yes') then 

             rmax=rco(1,l-1)/lambda(1,l-1)
             nrc1=nint(dlog(rmax/b+1.0d0)/a)+1

             eref=saveeorb(l-1,1)

             call polarization(a,b,rofi,rphi(1,l-1),vps(1,l-1),
     .            ve,nrc1,l-1,eref,g,nrc)
             
              dnrm=0.0d0
              do ir=2,nrc
                r=rofi(ir)
                f=g(ir)
                dnrm=dnrm+drdi(ir)*f*f
                g(ir)=f/r
                if(izeta.eq.1) rphi(ir,l)=f
                if(izeta.eq.1) rnrm(ir,l)=dnrm

                g(ir)=g(ir)/r**l

              enddo 
             
             endif 

         elseif(basistype.eq.'user') then 

C******Reading basis functions from the file provided by the user*****
          read(1,*) izetaread,npread,rcread

          rcread=min(rofi(nrval),rcread)
 
          if(izetaread.ne.izeta) then 
             write(6,'(/,2a,/,2(a,i2),/,a)') 
     .      'ATOM: Reading user-basis orbitals from file:',fname,
     .      'expected zeta=',izeta,' readed zeta=',izetaread,
     .      'check order!'
             stop
          endif      
          if ((npread+1).gt.nrmax) then 
             write(6,'(/,a,2a,/,a,i6)') 
     .       'ATOM: ERROR: ',
     .       'Too many grid points required to read functions in file ',
     .        fname,  
     .       'Parameter nrmax in atom.h must be increased to at least ',
     .        npread
              stop
          endif 

          do ir=1,npread
              read(1,*) r,y(ir)
              h(ir)=r
              if((r.gt.rcread).or.(dabs(rcread-r).lt.1.0d-5)) goto 88
          enddo  
          write(6,'(/,2a,/,a)') 
     .     'ATOM: WARNING: Basis orbitals read from file ',fname,
     .     'ATOM: WARNING: The required Rc is larger ',
     .     'than the maximum radial'
          write(6,'(a,f12.6)')
     .     'ATOM: WARNING: grid point specified in the file ',r
88        nrcread=ir 
          do ir=nrcread+1,npread
             read(1,*)h(ir),y(ir)
          enddo 
     
 
C**************Definition of the cut-off radius*************************
               nrc=nint(dlog(rcread/b+1.0d0)/a)+1
               nodd=mod(nrc,2)
               if(nodd.eq.0) then
                   nrc=nrc+1
               endif
               rco(izeta,l)=b*(dexp(a*(nrc-1))-1.0d0)*lambda(izeta,l)
C**********************************************************************     

C****Interpolation in the logaritmic mesh where pseudopotentials*******
C****are defined*******************************************************
              if(h(1).ne.0.0d0) then 
                do ir=1,npread
                   g(ir)=h(ir)
                   auxrho(ir)=y(ir) 
                enddo 
                do ir=1,npread
                   h(ir+1)=g(ir) 
                   y(ir+1)=auxrho(ir)
                enddo 
                h(1)=0.0d0 
                y(1)=0.0d0 
                npread=npread+1
                nrcread=nrcread+1
              endif 

              dnrm=0.0d0
              nr=1     
              do ir=2,nrc  
                 r=rofi(ir)
                 do jr=nr,npread
                    if(h(jr).ge.rofi(ir)) goto 87
                 enddo 
87               nr=jr

                 nmin=max(1,nr-npoint)
                 nmax=min(npread,nr+npoint)
                 nn=nmax-nmin+1
                 call ratint(h(nmin),y(nmin),nn,r,phi,dy)
                 dnrm=dnrm+drdi(ir)*(phi)**2
                 if(izeta.eq.1) rphi(ir,l)=phi 
                 if(izeta.eq.1) rnrm(ir,l)=dnrm
                 g(ir)=phi/(r**(l+1))
              enddo  
              g(1)=g(2)
          

 
C********************************************************************* 
          endif 
 

C**************Normalization of basis functions***********************
            eps=1.0d-4
            if(dabs(dnrm-1.0d0).gt.eps) then 
               do ir=1,nrc
                 g(ir)=g(ir)/dsqrt(dnrm)
                 if(izeta.eq.1) then
                    rphi(ir,l)=rphi(ir,l)/dsqrt(dnrm)
                    rnrm(ir,l)=rnrm(ir,l)/dnrm
                 endif
               enddo 
            endif
C*********************************************************************




C*Calculation of the mean value of kinetic and potential energy******* 
   
           
C    Potential and kinetic energy of the orbital before compression

          if((basistype.eq.'splitgauss').and.(izeta.gt.1)) then

             epot=0.0d0
             epot2=0.0d0
             ekin=0.0d0
             do ir=2,nrc
               r=rofi(ir)
               f=g(ir)*r**l
               epot=epot+
     .         drdi(ir)*(ve(ir)+vps(ir,l))*(f*r)**2
               epot2=epot2+
     .         drdi(ir)*vps(ir,l)*(f*r)**2
               dlapl=
     .           -((l*(l+1)-2*gexp*(2*l+3)*r**2+4*(gexp*r**2)**2)*f)
               dlapl=(dlapl+l*(l+1)*f)
               ekin=ekin +
     .         drdi(ir)*dlapl*f
               sum=sum+drdi(ir)*(f*r)**2
             enddo
             eorb=ekin+epot

          else

             ekin=0.0d0
             do ir=2,nrc 
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
             
             ekin=ekin/(lambda(izeta,l)**2)

C Potential energy after compression

             nrcomp=nint(dlog(rco(izeta,l)/b+1.0d0)/a)+1
             nmaxpoint=nrc 
             epot=0.0d0
             epot2=0.0d0
             do ir=1,nrcomp
                r=rofi(ir)
                r2=r/lambda(izeta,l)
                nr=nint(dlog(r2/b+1.0d0)/a)+1
                nmin=max(1,nr-npoint)
                nmax=min(nmaxpoint,nr+npoint)
                nn=nmax-nmin+1
                call ratint(rofi(nmin),g(nmin),nn,r2,rh,dy)
                rh=rh/dsqrt(lambda(izeta,l)**(2*l+3))
                epot=epot+
     .          drdi(ir)*(ve(ir)+vps(ir,l))*(rh*r**(l+1))**2
                epot2=epot2+
     .          drdi(ir)*vps(ir,l)*(rh*r**(l+1))**2
             enddo 
             eorb=ekin+epot

          endif 


           if(ipol.eq.'yes') then 
             write(6,'(/,(3x,a,i2),(3x,a,f12.6))')
     .       'izeta=',izeta,'Perturbative polarization orbital'
           endif 

          if((basistype.eq.'splitgauss').and.(izeta.gt.1)) then
 
           write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')        
     .         'izeta=',izeta,'gaussian exponent=',gexp,
     .         'rc=',rco(izeta,l),'energy=',eorb
 
         elseif((basistype.eq.'split').and.(izeta.gt.1)) then

            write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')
     .         'izeta =',izeta,
     .         'rmatch =',rco(izeta,l),
     .         'splitnorm =',spln,
     .         'energy =',eorb

          else

             write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')
     .          'izeta =',izeta,
     .          'lambda =',lambda(izeta,l),
     .          'rc =',rco(izeta,l),
     .          'energy =',eorb             

          endif

           saveeorb(l,izeta)=eorb

          write(6,'(a25,f12.6)') 'kinetic =',ekin
          write(6,'(a25,f12.6)') 'potential(screened) =',epot
          write(6,'(a25,f12.6)') 'potential(ionic) =',epot2


           
C*********TABLE WITH CUTOFF RADIUS FOR EACH ORBITAL*********************
               rcotb(izeta,l,is)=rco(izeta,l)
C***********************************************************************



C***********************************************************************
               rcoold(izeta,l,is)=rco(izeta,l)
               lambdaold(izeta,l,is)=lambda(izeta,l)
C***********************************************************************


C********************** write the basis *********************************
                   write(2,*) izeta, nrc, rco(izeta,l)
                   ndraw=1000
                   sum=0.0d0 
                   delt=rco(izeta,l)/(dble(ndraw-1)+1.0d-20)
                   sum=0.0d0 
                   do ir=1,ndraw
                       r=delt*(ir-1)
                       r=r/lambda(izeta,l)
                       nr=nint(dlog(r/b+1.0d0)/a)+1
                       nmin=max(1,nr-npoint)
                       nmax=min(nrc,nr+npoint)
                       nn=nmax-nmin+1
                       call ratint(rofi(nmin),g(nmin),nn,r,phi,dy)
                       phi=phi/dsqrt(lambda(izeta,l)**(2*l+3))
                       r=delt*(ir-1)
                       phi=phi*(r**dble(l+1)) 
                       sum=sum+delt*phi**2
                       write(2,*) r,phi,sum
                   enddo 


C***INTERPOLATION FOR THE TABLES WITH THE ATOMIC ORBITALS***************
 
     
            delt=rco(izeta,l)/(dble(ntbmax-1)+1.0d-20)
            table(1,norb,is)=delt
            
            table(2,norb,is)=dble(l)

            do 89 itb=1,ntbmax
                 r=delt*(itb-1)
                 r=r/lambda(izeta,l) 
                 nr=nint(dlog(r/b+1.0d0)/a)+1
                 nmin=max(1,nr-npoint)
                 nmax=min(nrc,nr+npoint)
                 nn=nmax-nmin+1
                 call ratint(rofi(nmin),g(nmin),nn,r,phi,dy)
                 phi=phi/dsqrt(lambda(izeta,l)**(2*l+3))
                 table(itb+2,norb,is)=phi
  89        continue
 80        continue

           if(basistype.eq.'user') then 
             do izeta=nzeta(l)+1,nzetl
                  read(1,*) izetaread,npread,rcread
                  do ir=1,npread
                    read(1,*) 
                  enddo 
             enddo 
           endif      

 
100        continue 
   
C************ close the file where basis has been written  *************
                 close(2)
C***********************************************************************



C***TOTAL NUMBER OF CALCULATED ORBITALS AND ITS INITIAL OCCUPATIONS*****
           nonew=0
           rcocc=0.0d0
           do 120 l=0,lmxo
             if((ql(l).gt.0.0d0).and.(rco(1,l).gt.rcocc))rcocc=rco(1,l)
              do  110 izeta=1,nzeta(l)
                  if(izeta.eq.1) then
                   if(flting.gt.0.0d0) then
                    do 105 m=1,2*l+1
                      if (nonew+m .le. maxos) then
                        q(nonew+m)=ql(l)/dble(2*l+1)
	              else
              write(6,"(2a)")'ATOM: WARNING: Number of orbitals ', 
     .         'per atom too small in calling routine'
	         write(6,*) nonew+m,maxos
                      endif
                      qold(nonew+m,is)=ql(l)/dble(2*l+1)
105                 continue
                   endif
                  elseif(izeta.gt.1) then 
                    do 106 m=1,2*l+1
                      if (nonew+m .le. maxos) then
                        q(nonew+m)=0.0d0
	              else
              write(6,"(2a)")'ATOM: WARNING: Number of orbitals ',
     .          'per atom too small in calling routine'
	         write(6,*) nonew+m,maxos
                      endif
                      qold(nonew+m,is)=0.0d0 
106                 continue
                  endif
                  nonew=nonew+2*l+1
110           continue
120        continue
           
C***********************************************************************
           
           write(6,'(a,i3)')
     .      'ATOM: Total number of Sankey-type orbitals:', nonew

           nomax(is)=nonew




C******* CONSTRUCTION OF THE NEW ELECTRONIC CHARGE DENSITY**************

         if (flting.gt.0.0d0) then

            write(6,'(/,2a)') 'ATOM: Valence configuration',
     .                      '(local Pseudopot. screening):'
            call cnfig(iz,config)
            linp=min(3,lmxo)
            do l=0,linp
               if(l.eq.0)
     .       write(6,'(7x,i2,a,f5.2,a)') config(0),'s(',ql(0),')'
               if(l.eq.1)
     .       write(6,'(7x,i2,a,f5.2,a)') config(1),'p(',ql(1),')'
               if(l.eq.2)
     .       write(6,'(7x,i2,a,f5.2,a)') config(2),'d(',ql(2),')'
               if(l.eq.3)
     .       write(6,'(7x,i2,a,f5.2,a)') config(3),'f(',ql(3),')'
           enddo

C**WE HAVE TO INTERPOLATE BECAUSE OF THE DIFFERENT COMPRESSION FACTORS**
            
          do ir=1,nrval
             rho(ir)=0.0d0 
             s(ir)=dsqrt(drdi(ir))
          enddo
 
          chval=0.0d0 
          zval=0.0d0
          do l=0,lmxo
             if(ql(l).gt.0.0d0) then 
              zval=zval+ql(l)
              nrc=nint(dlog(rco(1,l)/b+1.0d0)/a)+1
              rmax=rco(1,l)/lambda(1,l)
              maxpoint=nint(dlog(rmax/b+1.0d0)/a)+1

              do ir=2,nrc
                r=rofi(ir)/lambda(1,l)
                nr=nint(dlog(r/b+1.0d0)/a)+1
                nmin=max(1,nr-npoint)
                nmax=min(maxpoint,nr+npoint)
                nn=nmax-nmin+1
                call ratint(rofi(nmin),rphi(nmin,l),nn,r,rh,dy)
                rho(ir)=rho(ir)+ql(l)*rh**2/lambda(1,l)
                chval=chval+drdi(ir)*ql(l)*rh**2/lambda(1,l)
              enddo 
            endif
         enddo 
         rho(1)=0.0d0
         
         eps=1.0d-4
         if(dabs(chval-zval).gt.eps) then      
           do ir=2,nrval     
              rho(ir)=zval*rho(ir)/chval 
           enddo 
         endif     
    
C**CALCULATION OF THE HARTREE POTENTIAL DUE TO THE NEW VALENCE CHARGE**



          call vhrtre(rho,ve,rofi,drdi,s,nrval,a)
      

C*********LOCAL NEUTRAL-ATOM PSEUDOPOTENTIAL**************************
          nocc=nint(dlog(rcocc/b+1.0d0)/a)+2
          eps=1.0d-6
          nvlocal=0 
          do 280 ir=nrval,2,-1
c              dincv=vps(ir,lloc)+ve(ir)
               dincv=vlocal(ir)+ve(ir)
               if((abs(dincv).gt.eps).and.(nvlocal.eq.0)) nvlocal=ir+1
               h(ir)=dincv
 280      continue
          nvlocal=max(nvlocal,nocc)
          
C*********CUT-OFF RADIUS FOR THE LOCAL NEUTRAL-ATOM PSEUDOPOTENTIAL*****

          if(nvlocal.eq.nocc) then 
            rvlmx=rcocc
          else
            rvlmx=b*(dexp(a*(nvlocal-1))-1.0d0) 

          endif
  
          write(6,'(/,a,f10.6)') 
     .      'ATOM:  Cut-off radius for vlocal:',rvlmx

         if(rvlmx.gt.(rcocc+0.5d0)) then 
           write(6,"(2a,f12.5)")'ATOM: WARNING: ',
     .        'Cut-off radius for vlocal, Rvl =', rvlmx
           write(6,"(2a,f12.5)")'ATOM: WARNING: ',
     .        'Cut-off radius for charge density =', rcocc
           write(6,"(2a)")'ATOM: WARNING: ',
     .        'Check ATOM: Look for the sentence:'
           write(6,"(2a)")'ATOM: WARNING: ',
     .        'LOCAL NEUTRAL-ATOM PSEUDOPOTENTIAL'
           write(6,"(2a)")'ATOM: WARNING: ',
     .        'Increasing the tolerance parameter EPS'
           write(6,"(2a)")'ATOM: WARNING: ',
     .        'might be a good idea'
         endif

          h(1)= ( h(2)*rofi(3)**2 - h(3)*rofi(2)**2 ) /
     .          (      rofi(3)**2 -      rofi(2)**2 )


C CALCULATION OF THE ELECTROSTATIC SELF-ENERGY FOR THE 'CORE'*********
C ******************  CHARGE DENSITY ********************************


           a2b4=0.25d0*a*a 
           slf=0.0d0
           do ir=2,nvlocal-1
              if((ir.gt.2).and.(ir.lt.(nvlocal-1))) then 
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

           slfe(is)=slf


           elseif (flting.lt.0.0d0) then 
           
              slfe(is)=0.0d0 
            
           endif   

C***********************************************************************
 
C INTERPOLATION FOR THE TABLES OF THE LOCAL NEUTRAL-ATOM PSEUDOPOTENTIAL

           if (flting.gt.0.0d0) then  
              delt=rvlmx/(dble(ntbmax-1)+1.0d-20)
              table(1,0,is)=delt
              table(2,0,is)=rvlmx
           elseif(flting.lt.0.0d0) then 
              table(1,0,is)=0.0d0 
              table(2,0,is)=0.0d0 
           endif 
            
           do 300 itb=1,ntbmax
               r=delt*(itb-1)
               nr=nint(dlog(r/b+1.0d0)/a)+1
               nmin=max(1,nr-npoint)
               nmax=min(nvlocal,nr+npoint)
               nn=nmax-nmin+1
               if(flting.gt.0.0d0) then 
                 call ratint(rofi(nmin),h(nmin),nn,r,vloc,dy)
               else  
                 vloc=0.0d0
               endif
               table(itb+2,0,is)=vloc
300        continue
       


C**********************************************************************
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**********************************************************************



         elseif (iz.eq.-100) then


C AT THE MOMENT YOU HAVE TO SPECIFY A RADIUS FOR THE BESSEL 
C FUNCTIONS, SO NO AUTOMATIC BASIS IS AVAILABLE. 
       if (lmxo.eq.-1) then 
         write(6,'(2a)')
     .     'ATOM: ERROR: A radius must be specified for the',
     .     'floating Bessel functions.'
         write(6,'(a)')
     .     'ATOM: ERROR: it must be introduced using' 
         write(6,'(a)')
     .     'ATOM: ERROR: the PAO_basis_and_PS_lmax block format.'
           stop
       endif


       nzeta_overflow=0
       norb_overflow=0
       if (lmxo.eq.-1) then
           gen_basis='automatic'

C*IF THE BASIS SET IS AUTOMATIC, WE HAVE TO DEFINE SEVERAL PARAMETERS**
C******** FIRST WE READ DESIRED SIZE FOR THE BASIS SET*****************

           basis_size=fdf_string('PAO.BasisSize',basis_size_default)

           if(basis_size.eq.'STANDARD') basis_size='dzp'
           if(basis_size.eq.'standard') basis_size='dzp'
           if(basis_size.eq.'DZP')  basis_size='dzp'

           if(basis_size.eq.'DZ') basis_size='dz'

           if(basis_size.eq.'MINIMAL') basis_size='sz'
           if(basis_size.eq.'minimal')  basis_size='sz'
           if(basis_size.eq.'SZ')  basis_size='sz'


           if(basis_size.eq.'dzp') polorb=.true.
         
C******************Default L=0, L=1 if polarization is included*************
           lmxPAO=0
C***************************************************************************
           if(polorb) then
                 lmxo=lmxPAO+1
           else
                 lmxo=lmxPAO
           endif
           if(lmxo.gt.lmaxd) then 
              write(6,"(2a,i3)") 'ATOM: Parameter lmaxd ',
     .          'must be increased to at least ',lmax
              stop
           endif


           if((basis_size.eq.'dzp').
     .          or.(basis_size.eq.'dz')) then
                two = 2
                if(two.gt.nzetmx) then
                  write(6,"(2a)")'ATOM: ERROR: ',
     .            'Parameter nzetmx must be increased to at least 2'
                  stop
                endif

                do l=0,lmxPAO
                  nzeta_overflow=max(nzeta_overflow,2)
                  norb_overflow=norb_overflow+2*(2*l+1)
                  nzeta(l)=2
                enddo
           elseif(basis_size.eq.'sz') then
                do l=0,lmxPAO
                  nzeta_overflow=max(nzeta_overflow,1)
                  norb_overflow=norb_overflow+(2*l+1)
                  nzeta(l)=1
                enddo
           endif
           if(polorb) then
            numb_pol=fdf_integer('PAO.SplitPolarizationOrbitals',
     .         numb_pol_default)
               if(basis_size.eq.'sz') numb_pol=1
               norb_overflow=norb_overflow+numb_pol*(2*lmxo+1)
               nzeta_overflow=max(nzeta_overflow,numb_pol)
           endif

           if(polorb)then 
              nzeta(lmxo)=numb_pol
              do i=1,numb_pol
                 rco(i,lmxo)=rco(i,lmxPAO)
                 lambda(i,lmxo)=1.0d0
              enddo 
           endif 
       else
           gen_basis='block'
           lmxPAO=lmxo
           do l=0,lmxPAO
              nzeta_overflow=max(nzeta_overflow,nzeta(l))
              norb_overflow=norb_overflow+nzeta(l)*(2*l+1)
           enddo
  
           if(polorb) then
                lmxo=lmxPAO+1
                nzeta(lmxo)=numb_pol
                norb_overflow=norb_overflow+numb_pol*(2*lmxo+1)
                nzeta_overflow=max(nzeta_overflow,numb_pol)
                do i=1,numb_pol
                    rco(i,lmxo)=rco(i,lmxPAO)
                    lambda(i,lmxo)=1.0d0
                enddo 
           endif

       endif
  
C**********If dimension are not enough, return to recompile*************


       no=norb_overflow
       nztout=nzeta_overflow

       if(nzt.lt.nzeta_overflow) then
c        write(6,*) 'no',norb_overflow
c        write(6,*) 'nzeta',nzeta_overflow
c        write(6,*) 'lmxo',lmxo
c        write(6,*) 'lmxkb',lmxkb
         return
       endif
       if(maxl.lt.lmxo) then
c        write(6,*) 'no',norb_overflow
c        write(6,*) 'nzeta',nzeta_overflow
c        write(6,*) 'lmxo',lmxo
c        write(6,*) 'lmxkb',lmxkb
         return
       endif


C**********************************************************************
             do l=0,lmxo
               nzetold(l,is)=nzeta(l)
               if(nzeta(l).gt.nzetmx) then
                 write(6,"(2a)")'ATOM: ERROR: ',
     .             'Parameter nzetmx must be increased'
                 stop
               endif
             enddo
             lkbold(is)=0
             loold(is)=lmxo
             lmxPAOold(is)=lmxo
C**********************************************************************




C***********THE CASE OF FLOATING ORBITAL********************************       

       write(6,'(/,2a)') 'ATOM: Bases augmented with floating orbitals',
     .     ' consisting in Bessel functions' 
          

          if (lmxo.gt.lmaxd) then
            write(6,"(2a,i3)")'ATOM: ERROR: Paramater lmaxd ',
     .        'must be increased to at least',lmxo
            stop
          endif

C*****FOR THIS CASE ALL THE PSEUDOPOTENCIALS ARE ZERO*******************

         lmxkb=0
         nkb=0
         nkbmax(is)=0
         loctab(is)=0
         table(1,0,is)=0.0d0
         table(2,0,is)=0.0d0

         do 400 itb=1,ntbmax
            table(itb+2,0,is)=0.0d0
400      continue   
         slfe(is)=0.0d0   

C*********NO PSEUDO_CORE***********************************************
   
          do itb=1,ntbmax
             coretab(itb,1,is)=0.0d0 
             coretab(itb,2,is)=0.0d0 
          enddo  

C******INITIAL OCCUPATIONS ARE ALSO ZERO*******************************

          no=0
          rmax=0.0d0
          do 500 l=0,lmxo
             do 490 izeta=1,nzeta(l)
                 if(rco(izeta,l).gt.rmax) rmax=rco(izeta,l)
                 do 480 m=1,2*l+1
                    q(no+m)=0.0d0
                    qold(no+m,is)=0.0d0
480              enddo
                 no=no+2*l+1
490          continue
500       continue

          write(6,'(/,a,i3)')
     .     'ATOM: Total number of floating-type orbitals',no
          nomax(is)=no

C******CALCULATING THE FLOATING ORBITALS****************************** 
C*****SOLUTION OF A SPHERICAL POTENTIAL WELL**************************



C***********STANDART VALUES FOR MESH PARAMETERS*************************

          zt=1.0d0
          bb=6.0d0
          aa=80.0d0
          
          b=dexp(-bb)/zt
          a=1.0d0/aa

C***********SET UP THE MESH POINTS AND ITS DERIVATIVE******************

          nrcmx=nint(dlog(rmax/b+1.0d0)/a)+1
          nrcmx=nrcmx+1
          rpb=b
          ea=dexp(a)
          ea2=1.0d0
          do 550 ir=1,nrcmx
            drdi(ir)=a*rpb
            rofi(ir)=b*(ea2-1.0d0)
            rpb=rpb*ea
            ea2=ea2*ea 
550       continue

C***************Open a file to write the basis*************************
           fname= paste(atm_label,'.PAO.basis')
           open(unit=2,file=fname,status='unknown',form='formatted')
           write(2,*) lmxo
C**********************************************************************

          norb=0
          do  600 l=0,lmxo

C**************************  write the basis  **************************
              write(2,*) l,nzeta(l)
C***********************************************************************



C*******TABLE WITH THE NUMBER OF ORBITALS FOR EACH ANG. MOMENTUM*******
             nzettb(l,is)=nzeta(l)
C***********************************************************************


             a2b4=a*a*0.25d0

             do 570 ir=2,nrcmx
               s(ir)=drdi(ir)*drdi(ir)
               vtot=dble(l*(l+1))/(rofi(ir)**2)
               h(ir)=vtot*s(ir)+a2b4
570         continue
            h(1)=h(2)
            s(1)=s(2)




         
           do 580 izeta=1,nzeta(l)

           nrc=nint(dlog(rco(izeta,l)/b+1.0d0)/a)+1
           nodd=mod(nrc,2)
           nrc=nrc-1+nodd
           rmax=b*(dexp(a*(nrc-1))-1.0d0) 
           rco(izeta,l)=rmax
           

          z=zval
          dr=-1.0d5

          if(nzeta(l).lt.1) then 
            write(6,"(2a,i3)")'ATOM: WARNING: ',
     .          'No orbital has been calculated for l =', l
          endif


C********TABLE WITH THE CUTOFF RADIUS OF EACH ORBITAL*******************

           rcotb(izeta,l,is)=rco(izeta,l)
C***********************************************************************
           if(dabs(lambda(izeta,l)-1.0d0).gt.1.0d-3) then  
               write(6,'(/,a)') 
     .  'ATOM: WARNING: Scale factor is not active with Z=-100 option'
           endif
C***********************************************************************
           rcoold(izeta,l,is)=rco(izeta,l)
           lambdaold(izeta,l,is)=1.0d0
           lambda(izeta,l)=1.0d0
C***********************************************************************

            norb=norb+1
            nnode=izeta
            nprin=l+izeta
            e=-((zval/dble(nprin))**2)
c           e=0.0d0
            call egofv(h,s,nrc,eorb,g,y,l,z,a,b,rmax,nprin,
     .       nnode,dr)

           write(6,'(2(3x,a,i2),2(3x,a,f10.6))')
     .       'l=',l,'nzeta',izeta,'rc',rco(izeta,l),'energy=',eorb

            dnrm=0.0d0
            do 590 ir=2,nrc
               r=rofi(ir)
               f=g(ir)
               dsq=dsqrt(drdi(ir))
               f=f*dsq
               dnrm=dnrm+drdi(ir)*f*f
               g(ir)=f/r
             
               g(ir)=g(ir)/r**l

 590        continue
        
             eps=1.0d-4
             if(dabs(dnrm-1.0d0).gt.eps) then
                do ir=2,nrc
                    g(ir)=g(ir)/dsqrt(dnrm)
                enddo 
             endif

            g(1)= ( g(2)*rofi(3)**2 - g(3)*rofi(2)**2 ) /
     .            (      rofi(3)**2 -      rofi(2)**2 )
             

C********************** write the basis *********************************
                   write(2,*) izeta, nrc, rco(izeta,l)
                   ndraw=1000
                   sum=0.0d0
                   delt=rco(izeta,l)/(dble(ndraw-1)+1.0d-20)
                   sum=0.0d0
                   do ir=1,ndraw
                       r=delt*(ir-1)
                       r=r
                       nr=nint(dlog(r/b+1.0d0)/a)+1
                       nmin=max(1,nr-npoint)
                       nmax=min(nrc,nr+npoint)
                       nn=nmax-nmin+1
                       call ratint(rofi(nmin),g(nmin),nn,r,phi,dy)
                       phi=phi
                       r=delt*(ir-1)
                       phi=phi*(r**dble(l+1))
                       sum=sum+delt*phi**2
                       write(2,*) r,phi,sum
                   enddo


C****INTERPOLATION FOR THE TABLES WITH THE ATOMIC ORBITALS*************
 
     
            delt=rco(izeta,l)/(dble(ntbmax-1)+1.0d-20)
            table(1,norb,is)=delt
            
            table(2,norb,is)=dble(l)

            do 589 itb=1,ntbmax
                 r=delt*(itb-1)
                 nr=nint(dlog(r/b+1.0d0)/a)+1
                 nmin=max(1,nr-npoint)
                 nmax=min(nrc,nr+npoint)
                 nn=nmax-nmin+1
                 call ratint(rofi(nmin),g(nmin),nn,r,phi,dy)
                 table(itb+2,norb,is)=phi
  589       continue

                       
 580       continue

 
600        continue 
           write(6,*) ' '
 




         elseif (iz.eq.0) then



            write(6,*) 'ATOM: re-initialization of the tables'
            write(6,*) 'ATOM: all data will be set to zero'
  
            called=.false.
            
            ismax=0
            isold=0
            do i=1,nsmax
              nomax(i)=0 
              nkbmax(i)=0
              iztb(i)=0
              izold(i)=0
              lkbold(i)=0 
              lkbin(i)=0
              loold(i)=0
              loin(i)=0
              loctab(i)=0
              lmxPAOold(i)=0
              slfe(i)=0.0d0
              do itb=1,ntbmax+1
                 coretab(itb,1,i)=0.0d0 
                 coretab(itb,2,i)=0.0d0 
              enddo 

              do l=0,lmaxd
                do izeta=1,nzetmx
                 rcoold(izeta,l,i)=0.0d0
                 rcotb(izeta,l,i)=0.0d0
                 rcoin(izeta,l,i)=0.0d0
                 lambdaold(izeta,l,i)=0.0d0 
                 lambdain(izeta,l,i)=0.0d0
                enddo 
                 nzetold(l,i)=0
                 rctb(l,i)=0.0d0
                 nzettb(l,i)=0
                 nzetin(l,i)=0
              enddo
              do l=1,lmx2*nzetmx
                 qold(l,i)=0.0d0
              enddo
              do l=-(lmaxd+1),nzetmx*(lmaxd+1)
                 do ix=1,ntbmax
                    table(ix,l,i)=0.0d0
                    tab2(ix,l,i)=0.0d0
                 enddo
                 table(ntbmax+1,l,i)=0.0d0
                 table(ntbmax+2,l,i)=0.0d0
              enddo
            
             enddo
             is=0
             no=0
             nkb=0
             
            goto 1000
            
         endif
 


C***TABLES WITH THE SECOND DERIVATIVE FOR CUBIC SPLINE INTERPOLATION***

         
C*********FOR K-B PROJECTORS*******************************************

        if ((flting.gt.0.0d0).and.(iz.ne.-100)) then 
        indx=0
        do  700 l=0,lmxkb
c           if(l.eq.lloc) goto 700
            indx=indx+1
            delt=table(1,-indx,is)

            yp1=0.d0
            ypn=1.0d50

            call spline(delt,table(3,-indx,is),ntbmax,
     .        yp1,ypn,tab2(1,-indx,is),aux)

 

700     continue

        endif 


C*******FOR LOCAL PSEUDOPOTENTIAL*************************************


        if ((flting.gt.0.0d0).and.(iz.ne.-100)) then
           delt=table(1,0,is)
           yp1=0.d0
           ypn=1.0d50
 
             call spline(delt,table(3,0,is),ntbmax,
     .        yp1,ypn,tab2(1,0,is),aux)
        endif 


C*******FOR THE ATOMIC ORBITALS***************************************



           norb=0
           do 800 l=0,lmxo
              
             do 850 izeta=1,nzeta(l)
   
               norb=norb+1
      
               delt=table(1,norb,is)

               yp1=0.d0
               ypn=1.0d50 

               call spline(delt,table(3,norb,is),ntbmax,
     .           yp1,ypn,tab2(1,norb,is),aux)

 
850          continue


800       continue

C**Print an alternative input if the basis has been generated in******** 
C****************an automatic way***************************************

         if(is.eq.ntotsp)then 
            if(gen_basis.eq.'automatic') then 
             write(6,'(/,a,73(1h*))') 'ATOM: '
             write(6,'(/,a)')
     .         'ATOM: Basis sets have been generated automatically'
             write(6,'(a)') 'ATOM: using the options:'
             write(6,'(2a)') 'ATOM: PAO.BasisSize ',basis_size
             write(6,'(2a)') 'ATOM: PAO.BasisType  ', basistype
             write(6,'(a,f12.6,a)')
     .             'ATOM: PAO.EnergyShift', eshift ,' Ry '
             write(6,'(a)') 
     .  'ATOM: An alternative input to obtain the same basis set:'
             write(6,'(/,a,72(1h*))') 'INPUT: '

             write(6,'(2a)')'PAO.BasisType  ', basistype
             write(6,'(a,f12.6,a)')'PAO.EnergyShift', eshift ,' Ry '
             
             if((basistype.eq.'split').or.
     .        (numb_pol.gt.1)) 
     .            write(6,'(a,f12.6)') 'PAO.SplitNorm',splnorm
             if(polorb) then 
                write(6,'(a)') 'PAO.PolarizationOrbitals true'
             else
                write(6,'(a)') 'PAO.PolarizationOrbitals false'
             endif 

             if(polorb)
     .        write(6,'(a,i2)')'PAO.SplitPolarizationOrbitals  ',
     .              numb_pol
             write(6,'(a)')
     .         '%block PAO_basis_and_PS_lmax        # Define Basis set'
                do i=1,ntotsp
            write(6,'(/,2i3,2i2,a)') i,izold(i),loold(i),lkbold(i),
     .              '                   # Index, Z, LmaxPAO, LmaxPS'    
                   do l=0,lmxPAOold(i)
                     write(6,*) l,nzetold(l,i) ,
     .               '              # l, Nzeta  '
                    write(6,*) 
     .                 (rcoold(izeta,l,i),izeta=1,nzetold(l,i)),
     .                         '  # rc(izeta=1,Nzeta)(Bohr)'
                     write(6,*) 
     .                 (lambdaold(izeta,l,i),izeta=1,nzetold(l,i)),
     .                        '  # scaleFactor(izeta=1,Nzeta)'
                   enddo 
           
             enddo 
            write(6,'(a)') 
     .        '%endblock PAO_basis_and_PS_lmax     '

         else
            write(6,'(/,a,73(1h*))') 'ATOM: '
           write(6,'(/,a)')
     .         'ATOM: Basis set has been generated using block data'        
           write(6,'(a)')
     .   'ATOM: Active additional parameters' 
           write(6,'(2a)')'PAO.BasisType  ', basistype
          if(polorb) then
            write(6,'(a)') 'PAO.PolarizationOrbitals true'
          else
            write(6,'(a)') 'PAO.PolarizationOrbitals false'
          endif
           if(polorb )
     .        write(6,'(a,i2)')'PAO.SplitPolarizationOrbitals  ',
     .          numb_pol
              do l=0,lmxPAO
                if(rcoin(1,l,is).lt.1.0d-5) then 
                   write(6,'(a,f12.6,a)')
     .                  'PAO.EnergyShift', eshift,' Ry'
                endif 
             enddo 
             if((basistype.eq.'split').or.
     .        (numb_pol.gt.1))
     .        write(6,'(a,f12.6)') 'PAO.SplitNorm',splnorm
 
               
            endif  
          
           endif 
              
1000     continue

         write(6,'(a,73(1h*))') 'ATOM: '
        
        return

        end




