C $Id: atom_functions.f,v 1.7 1999/02/26 21:18:09 daniel Exp $

C************************************************************************
C  This file contains a set of routines which provide all the information
C  about the basis set, pseudopotential, atomic mass, etc... of all the
C  chemical species present in the calculation.
C
C   Written by D. Sanchez-Portal, Sept. 1998
C
C************************************************************************
C  The routines contained in this file can only be called after they
C  are initialized by calling the subroutine 'atom' for all the 
C  different chemical species in the calculation:
C    
C   REAL*8  RCUT(IS,IO)          : Function which returns cutoff
C                                        radius of:
C                                        a) basis orbitals (IO > 0)
C                                        b) KB projectors  (IO < 0)
C                                        c) Local pseudopotential(IO=0)
C  PHIATM(IS,IO,R,PHI,GRPHI)     : Routine which returns values and 
C                                        gradients of:
C                                        a) basis orbitals (IO > 0)
C                                        b) KB projectors  (IO < 0)
C                                        c) Local pseudopotential(IO=0) 
C  RPHIATM(IS,IO,R,PHI,DPHIDR)   : Like PHIATM, just for the radial 
C                                        component of the orbitals 
C  VLOCAL(IS,R,V,GRV)            : Routine which returns local part 
C                                        of neutral-atom pseudopotential
C                                        and its gradient
C  REAL*8  RCORE(IS)             : Function which returns cutoff
C                                         radius of the pseudo-core
C  CHCORE(IS,R,CH,GRCH)          : Routine which returns the 
C                                        pseudo-core for non-linear
C                                        core corrections in the
C                                        XC potential
C  PSCH(IS,R,CH,GRCH)            : Routine which returns the
C                                        local-pseudotential charge
C                                        density (laplacian of the
C                                        local-pseudotential)
C  PSOVER(IS1,IS2,R,E,DEDR)      : Routine which returns the
C                                       correction to the electrostatic
C                                       interaction energy of the ions
C                                       due to the overlap of the
C                                       local-pseudotential charge
C                                       densities.
C  REAL*8 EPSKB(IS,IO)           : Function which returns the KB
C                                       projector energies
C  REAL*8 UION(IS)               : Function which returns the  
C                                      electrostatic self-energy of 
C                                      the charge density of the ion
C  INTEGER IZOFIS(IS)            : Function which returns atomic
C                                       number 
C  INTEGER IZVALFIS(IS)          : Function which returns the 
C                                      valence charge
C  INTEGER LOFIO(IS,IO)          : Function which returns total 
C                                      angular momentum quantum number
C                                      of orbitals 
C  INTEGER MOFIO(IS,IO)          : Function which returns magnetic
C                                      quantum number of the orbitals
C  INTEGER LOMAXFIS(IS)          : Function which returns the maximum 
C                                      angular momentum of the basis functions
C                                      for a given specie
C  INTEGER LMXKBFIS(IS)          : Function which returns the maximum 
C                                      angular momentum of the KB projectors
C                                      for a given specie.
C  INTEGER NZTFL(IS,L)           : Function which returns the number of 
C                                     basis functions with a given angular 
C                                     momentum and for a given species
C  INTEGER NOFIS(IS)             : Functions wich returns the total number 
C                                      of basis functions for the species is
C  INTEGER NKBFIS(IS)            : Functions wich returns the total number
C                                     of KB projectors for the species is.
C  REAL*8  MASSFIS(IS)           : Functions wich returns the atomic mass
C
C  CHARACTER LABELFIS(IS)        : Functions wich returns the atomic label 
C  CHARACTER SYMFIO(IS,IO)       : Symmetry (s,px,py,pz,dz2,....) for 
C                                     the base functions. 
C  LOGICAL  POL(IS,IO)           : If POL=.true. the orbital is a 
C                                   perturbative polarization one.
C  INTEGER  CNFIGFL(IS,L)       : Function which returns the atomic 
C                                  configuration.
C  REAL*8   ATMPOPFL(IS,L)      : Function which returns the population
C                                  of the different l-shells in the 
C                                  atomic ground state configuration.       
C  REAL*8   ATMPOPFIO(IS,IO)    : Function which returns the population
C                                  of the basis orbitals in the 
C                                  atomic ground state configuration.
C***********************UNITS*******************************************
C    Distances in Bohr.
C    Energies in Rydbergs.
C***********************************************************************







       INTEGER FUNCTION IZOFIS( IS )

C***********************************************************************
C  Returns atomic number of a given species, assigned by routine ATOM
C  Written by D. Sanchez-Portal. Oct, 1996 
C**************************INPUT****************************************
C  INTEGER IS     : Species index
C**************************OUTPUT***************************************
C  INTEGER IZOFIS : Atomic number
C**************************BEHAVIOUR************************************
C  0) Before using IZOFIS, the pseudopotential must be initialized 
C     by calling ATOM for each atomic species required.
C  1) Prints a message and stops when no data exists for IS.
C**********************************************************************
         
         implicit none 

         include 'atom.h'

 
         integer
     .    is,izsave(nsmax),ismax, nomax(nsmax), nkbmax(nsmax)

          

         common/cmiz/izsave
         common/control/ismax,nomax,nkbmax

         if ((is.lt.1).or.(is.gt.ismax)) then 
            write(6,*) 'IZOFIS: THERE ARE NO DATA FOR IS=',IS
            write(6,*) 'IZOFIS: ISMIN= 1, ISMAX= ',ismax
            STOP
         endif


 
         izofis=izsave(is)


         return
         end





           INTEGER FUNCTION IZVALFIS( IS )

C***********************************************************************
C  Returns valence charge of a given species, assigned by routine ATOM
C  Written by D. Sanchez-Portal. Aug. 1998
C**************************INPUT****************************************
C  INTEGER IS     : Species index
C**************************OUTPUT***************************************
C  INTEGER IZVALFIS : Valence charge
C**************************BEHAVIOUR************************************
C  0) Before using IZALFIS, the pseudopotential must be initialized 
C     by calling ATOM for each atomic species required.
C  1) Prints a message and stops when no data exists for IS.
C********************************************************************** 


         implicit  none 

         include 'atom.h'

 
         integer
     .    is,izvaltb(nsmax),ismax, nomax(nsmax), nkbmax(nsmax)

          

         common/cizvaltb/izvaltb
         common/control/ismax, nomax, nkbmax

         if ((is.lt.1).or.(is.gt.ismax)) then 
            write(6,*) 'IZALFIS: THERE ARE NO DATA FOR IS=',IS
            write(6,*) 'IZALFIS: ISMIN= 1, ISMAX= ',ismax
            STOP
         endif


 
         izvalfis=izvaltb(is)


         return

         end





          subroutine psch(is,r,ch,grch)


C***********************************************************************
C Returns 'local-pseudotential charge density'.
C Written by D.Sanchez-Portal. March, 1997
C************************INPUT******************************************
C INTEGER IS     :  Species index
C REAL*8  R(3)   :  Point vector, relative to atom
C ***********************OUTPUT*****************************************
C REAL*8  CH     :  Value of local-pseudotential charge density.
C REAL*8  GRCH(3):  Gradient of local-pseudotential charge density.
C************************UNITS******************************************
C Distances in Bohr
C Energies in Rydbergs
C Density in electrons/Bohr**3
C************************BEHAVIOUR**************************************
C  0) Before using PSCH, the pseudopotential must be initialized
C     by calling ATOM for each atomic species required
C  1) Prints a message and stops when no data exits for IS.
C  2) Returns exactly zero when |R| > Rchloc
C***********************************************************************
 
          implicit  none
          
          include 'atom.h'

          double precision 
     .     r(3),grch(3)

          integer
     .     is

C*****************Internal variables***************************************
C
          integer i
          double precision 
     .      rmod, dchdr, dloc, delt, rcmx, ch
C
C***********************************************************************
 
C****************Variables in common blocks******************************
C
         double precision
     .     chloctab((ntbmax+1),2,nsmax)

          integer   ismax, nomax(nsmax), nkbmax(nsmax) 

          common/control/ismax, nomax, nkbmax
          common/cmchloc/chloctab
C
C***********************************************************************


          if ((is.lt.1).or.(is.gt.ismax)) then
            write(6,*) 'PSCH: THERE ARE NO DATA FOR IS=',IS
            write(6,*) 'PSCH: ISMIN= 1, ISMAX= ',ismax
            STOP
          endif

          dloc=chloctab(1,2,is)
          if(dabs(dloc).lt.1.0d-8) then 
            ch=0.0d0 
            grch(1)=0.0d0 
            grch(2)=0.0d0
            grch(3)=0.0d0 
            return
          endif

          delt=chloctab(1,1,is) 
          rcmx=delt*(ntbmax-1)
          rmod=0.0d0
          do i=1,3
            rmod=rmod+r(i)*r(i)
          enddo
          rmod=dsqrt(rmod)

          if(rmod.gt.rcmx-1.d-12) then
             ch=0.0d0
             grch(1)=0.0d0
             grch(2)=0.0d0
             grch(3)=0.0d0
        
          else
            call splint(delt,chloctab(2,1,is),chloctab(2,2,is),ntbmax,
     .        rmod,ch,dchdr)

             rmod=rmod+1.0d-20
           
             grch(1)=dchdr*r(1)/rmod
             grch(2)=dchdr*r(2)/rmod
             grch(3)=dchdr*r(3)/rmod
       

          endif

 
          return
 
          end


       REAL*8 FUNCTION ATMPOPFIO (IS,IO)
C**********************************************************************
C Returns the population of the atomic basis orbitals in the atomic 
C ground state configuration.
C Written by D.Sanchez-Portal. Oct, 1998
C************************INPUT******************************************
C    INTEGER  IS   : Species index
C    INTEGER  IO   : Orbital index (within atom)
C************************OUTPUT*****************************************
C    REAL*8 ATMPOPFIO: Population of the atomic basis orbitals in the  
C                      atomic ground state configuration.
C************************BEHAVIOUR**************************************
C  0) Before using ATMPOPFIO, the pseudopotential must be initialized
C      by calling ATOM for each atomic species required.
C  1) Prints a message and stops when no data exist for IS and/or IO
C  2) Returns zero for IO = 0
C***********************************************************************


       implicit none

       include 'atom.h'

       integer  is,io


C
C********************************************************************
C************Variables in common blocks******************************
C
       integer
     .  ismax,nomax(nsmax),
     .  nkbmax(nsmax) 

       integer maxos
       parameter (maxos=2*nzetmx*lmx2) 

       double precision  
     .    qtb(maxos,nsmax)
C
C*******************************************************************

C*******************************************************************
C
       common/control/ismax,nomax,nkbmax 
       common/cmq/qtb
C
C*******************************************************************
 
        if ((is.lt.1).or.(is.gt.ismax)) then
          write(6,*) 'ATMPOPFIO: THERE ARE NO DATA FOR IS=',IS
          write(6,*) 'ATMPOPFIO: ISMIN= 1, ISMAX= ',ismax
          STOP
        endif
        if((io.gt.nomax(is)).or.(io.lt.1)) then
          write(6,*) 'ATMPOPFIO: THERE ARE NO DATA FOR IO=',IO
          write(6,*) 'ATMPOPFIO: IOMIN= 1',
     .     ' IOMAX= ',nomax(is)
          STOP
        endif
 


 
        atmpopfio=qtb(io,is) 

 

        return 
        end
 
  

       INTEGER FUNCTION LOFIO (IS,IO)
C**********************************************************************
C Returns total angular momentum quantum number of a given atomic basis
C   basis orbital or Kleynman-Bylander projector.
C Written by D.Sanchez-Portal. Oct, 1996
C************************INPUT******************************************
C    INTEGER  IS   : Species index
C    INTEGER  IO   : Orbital index (within atom)
C                    IO > 0 => Basis orbitals
C                    IO < 0 => Kleynman-Bylander projectors
C************************OUTPUT*****************************************
C   INTEGER LOFIO  : Quantum number L of orbital or KB projector
C************************BEHAVIOUR**************************************
C  0) Before using LOFIO, the pseudopotential must be initialized
C      by calling ATOM for each atomic species required.
C  1) Prints a message and stops when no data exist for IS and/or IO
C  2) Returns zero for IO = 0
C***********************************************************************
       

       implicit none

       include 'atom.h'

       integer  is,io

C*************Internal variables*************************************
C
         integer l, norb, izeta, ipol, nkb



C
C********************************************************************
C************Variables in common blocks******************************
C
       integer
     .  nzetasave(0:lmaxd,nsmax),ismax,nomax(nsmax),
     .  nkbmax(nsmax),npolorbsave(lmaxd,nsmax),
     .  lmxosave(nsmax), lmxkbsave(nsmax)
C
C*******************************************************************

C*******************************************************************
C
       common/cmzeta/nzetasave
       common/control/ismax,nomax,nkbmax
       common/cmpolorb/npolorbsave
       common/cmlmxo/lmxosave
       common/cmlmxkb/lmxkbsave
C
C*******************************************************************
  
        if ((is.lt.1).or.(is.gt.ismax)) then 
          write(6,*) 'LOFIO: THERE ARE NO DATA FOR IS=',IS
          write(6,*) 'LOFIO: ISMIN= 1, ISMAX= ',ismax
          STOP
        endif
        if((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then 
          write(6,*) 'LOFIO: THERE ARE NO DATA FOR IO=',IO
          write(6,*) 'LOFIO: IOMIN= ',-nkbmax(is),
     .     ' IOMAX= ',nomax(is)
          STOP
        endif
 
       if (io.gt.0) then

        norb=0
        do 10 l=0,lmxosave(is)
          do 5 izeta=1,nzetasave(l,is)
            norb=norb+(2*l+1)
            if(norb.ge.io) goto 30
 5        continue
10      continue

        do  20 l=1,min(lmxosave(is)+1,lmaxd)
            do 15 ipol=1, npolorbsave(l,is)
              norb=norb+(2*l+1)
              if(norb.ge.io) goto 40
15          continue
20      continue
        write(6,*) 'LOFIO: ERROR: ORBITAL INDEX IO=',IO
        write(6,*) 'LOFIO: ERROR: NOT FOUND'
        stop

30      lofio=l
        return

40      lofio=l
        return


       elseif(io.lt.0) then


        nkb=0
        do 50 l=0,lmxkbsave(is)
          nkb=nkb-(2*l+1)
          if(nkb.le.io) goto 60
50      continue 

60      lofio=l       

        elseif (io.eq.0) then

        lofio=0

        endif
 
         
        return
         
          end 










       INTEGER FUNCTION MOFIO (IS,IO)
C**********************************************************************
C Returns magnetic quantum number of a given atomic basis
C   basis orbital or Kleynman-Bylander projector.
C Written by D.Sanchez-Portal. Oct, 1996
C************************INPUT******************************************
C    INTEGER  IS   : Species index
C    INTEGER  IO   : Orbital index (within atom)
C                    IO > 0 => Basis orbitals
C                    IO < 0 => Kleynman-Bylander projectors
C************************OUTPUT*****************************************
C   INTEGER MOFIO  : Quantum number M of orbital or KB projector
C************************BEHAVIOUR**************************************
C  0) Before using MOFIO, the pseudopotential must be initialized
C      by calling ATOM for each atomic species required.
C  1) Prints a message and stops when no data exist for IS and/or IO
C  2) Returns zero for IO = 0
C***********************************************************************
       

        implicit none

       include 'atom.h'

       integer  is,io

C*************Internal variables*************************************
C
         integer l, norb, izeta, ipol, nkb, lorb, lkb



C
C********************************************************************
C************Variables in common blocks******************************
C
       integer
     .  nzetasave(0:lmaxd,nsmax),ismax,nomax(nsmax),
     .  nkbmax(nsmax),npolorbsave(lmaxd,nsmax),
     .  lmxosave(nsmax), lmxkbsave(nsmax)
C
C*******************************************************************

C*******************************************************************
C
       common/cmzeta/nzetasave
       common/control/ismax,nomax,nkbmax
       common/cmpolorb/npolorbsave
       common/cmlmxo/lmxosave
       common/cmlmxkb/lmxkbsave
C
C*******************************************************************

        if ((is.lt.1).or.(is.gt.ismax)) then
          write(6,*) 'MOFIO: THERE ARE NO DATA FOR IS=',IS
          write(6,*) 'MOFIO: ISMIN= 1, ISMAX= ',ismax
          STOP
        endif
        if((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
          write(6,*) 'MOFIO: THERE ARE NO DATA FOR IO=',IO
          write(6,*) 'MOFIO: IOMIN= ',-nkbmax(is),
     .     ' IOMAX= ',nomax(is)
          STOP
        endif
      if (io.gt.0) then

        norb=0
        do 10 l=0,lmxosave(is)
          do 5 izeta=1,nzetasave(l,is)
            norb=norb+(2*l+1)
            if(norb.ge.io) goto 30
 5        continue
10      continue

        do  20 l=1,min(lmxosave(is)+1,lmaxd)
            do 15 ipol=1, npolorbsave(l,is)
              norb=norb+(2*l+1)
              if(norb.ge.io) goto 40
15          continue
20      continue
        write(6,*) 'MOFIO: ERROR: ORBITAL INDEX IO=',IO
        write(6,*) 'MOFIO: ERROR: NOT FOUND'
        stop

30      lorb=l 
        mofio=io-norb+lorb
        return

40      lorb=l 
        mofio=io-norb+lorb
        return


       elseif(io.lt.0) then


        nkb=0
        do 50 l=0,lmxkbsave(is)
          nkb=nkb-(2*l+1)
          if(nkb.le.io) goto 60
50      continue

60      lkb=l
        mofio=-io+nkb+lkb
        elseif (io.eq.0) then

        mofio=0

        endif
        
        end 














          subroutine vlocal(is,r,v,grv)


C***********************************************************************
C Returns local part of neutral-atom Kleynman-Bylander pseudopotential.
C Written by D.Sanchez-Portal. Oct, 1996
C************************INPUT******************************************
C INTEGER IS     :  Species index
C REAL*8  R(3)   :  Point vector, relative to atom
C ***********************OUTPUT*****************************************
C REAL*8  V      :  Value of local pseudopotential
C REAL*8  GRV(3) :  Gradient of local pseudopotential
C************************UNITS******************************************
C Distances in Bohr
C Energies in Rydbergs
C************************BEHAVIOUR**************************************
C  0) Before using VLOCAL, the pseudopotential must be initialized
C     by calling ATOM for each atomic species required
C  1) Prints a message and stops when no data exits for IS.
C  2) Returns exactly zero when |R| > RCUT(IS,0)
C***********************************************************************
 
          implicit none
          
          include 'atom.h'
        
        double precision
     .      r(3), v, grv(3)

        integer
     .      is

C*************Internal variables*************************************
C
         integer i

         double precision  rvmx, rmod, dvdr, delt

C
C********************************************************************


C************Variables in common blocks******************************
C
       double precision
     .  table((ntbmax+2),-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tab2(ntbmax,-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tabpol((ntbmax+2),nzetmx*lmaxd,nsmax),
     .  tab2pol(ntbmax,nzetmx*lmaxd,nsmax)


       integer
     .  ismax,nomax(nsmax),
     .  nkbmax(nsmax)
C
C*******************************************************************

C*******************************************************************
C
       common/cmtab/table,tabpol 
       common/cmspline/tab2, tab2pol
       common/control/ismax,nomax,nkbmax
C
C*******************************************************************


          if ((is.lt.1).or.(is.gt.ismax)) then
            write(6,*) 'VLOCAL: THERE ARE NO DATA FOR IS=',IS
            write(6,*) 'VLOCAL: ISMIN= 1, ISMAX= ',ismax
            STOP
          endif

 
          delt=table(1,0,is)
          rvmx=table(2,0,is)
         
          rmod=0.0d0
          do i=1,3
            rmod=rmod+r(i)*r(i)
          enddo
          rmod=dsqrt(rmod)
           
          if(rmod.gt.rvmx) then
             v=0.0d0
             grv(1)=0.0d0
             grv(2)=0.0d0
             grv(3)=0.0d0
        
          else
            call splint(delt,table(3,0,is),tab2(1,0,is),ntbmax,
     .        rmod,v,dvdr)

             rmod=rmod+1.0d-20
           
             grv(1)=dvdr*r(1)/rmod
             grv(2)=dvdr*r(2)/rmod
             grv(3)=dvdr*r(3)/rmod
       

          endif

 
          return
 
          end


       INTEGER FUNCTION LOMAXFIS (IS)
C**********************************************************************
C Returns maximum angular momentum of the basis functions of a atomic
C specie
C Written by D.Sanchez-Portal. Aug., 1998
C************************INPUT******************************************
C    INTEGER  IS   : Species index
C************************OUTPUT*****************************************
C   INTEGER LOMAXFIS  : Angular momentum L of the orbital
C************************BEHAVIOUR**************************************
C  0) Before using LOMAXFIS, the pseudopotential must be initialized
C      by calling ATOM for each atomic species required.
C  1) Prints a message and stops when no data exist for IS 
C***********************************************************************
       
         implicit none
         include 'atom.h'
        
         integer
     .      is,ismax,
     .      nomax(nsmax),nkbmax(nsmax),lmxosave(nsmax),
     .      npolorbsave(lmaxd,nsmax)


         common/control/ismax,nomax,nkbmax
         common/cmlmxo/lmxosave  
         common/cmpolorb/npolorbsave 



        if ((is.lt.1).or.(is.gt.ismax)) then 
          write(6,*) 'LOMAXFIS: THERE ARE NO DATA FOR IS=',IS
          write(6,*) 'LOMAXFIS: ISMIN= 1, ISMAX= ',ismax
          STOP
        endif
 
          
          if((lmaxd.ge.lmxosave(is)+1).and.
     .            (npolorbsave(lmxosave(is)+1,is).gt.0)) then 
           lomaxfis=lmxosave(is)+1
          else
           lomaxfis=lmxosave(is)
          endif 
         
          return
         
          end 










       INTEGER FUNCTION LMXKBFIS (IS)
C**********************************************************************
C Returns maximum angular momentum of the Kleinman-Bylander projectors
C for an atomic specie
C Written by D.Sanchez-Portal. Aug., 1998
C************************INPUT******************************************
C    INTEGER  IS   : Species index
C************************OUTPUT*****************************************
C   INTEGER LMXKBFIS  : Angular momentum L of the projector
C************************BEHAVIOUR**************************************
C  0) Before using LMXKBFIS, the pseudopotential must be initialized
C      by calling ATOM for each atomic species required.
C  1) Prints a message and stops when no data exist for IS 
C***********************************************************************
       

         include 'atom.h'
        
         integer
     .      is,ismax,
     .      nomax(nsmax),nkbmax(nsmax),lmxkbsave(nsmax)


         common/control/ismax,nomax,nkbmax
         common/cmlmxkb/lmxkbsave
  
        if ((is.lt.1).or.(is.gt.ismax)) then 
          write(6,*) 'LMXKBFIS: THERE ARE NO DATA FOR IS=',IS
          write(6,*) 'LMXKBFIS: ISMIN= 1, ISMAX= ',ismax
          STOP
        endif
 
 

          lmxkbfis=lmxkbsave(is)
         
          return
         
          end 






       INTEGER FUNCTION CNFIGFL(IS,L)
C**********************************************************************
C Returns the valence-shell configuration in the atomic ground state
C (i.e. the principal quatum number for orbitals of angular momentum l
C in the valence shell)
C Written by D.Sanchez-Portal. Aug, 1998
C************************INPUT******************************************
C    INTEGER  IS   : Species index
C    INTEGER  L    : Angular momentum of the shell
C************************OUTPUT*****************************************
C   INTEGER CNFIGFL: Valence-shell configuration in the atomic 
C                     ground state (Principal quantum number).
C************************BEHAVIOUR**************************************
C  0) Before using CNFIGFL, the pseudopotential must be initialized
C      by calling ATOM for each atomic species required.
C  1) Prints a message and stops when no data exist for IS and/or L
C***********************************************************************

         implicit none
         include 'atom.h'

         integer
     .      is,l


C************Variables in common blocks******************************
C

       integer
     .  ismax,nomax(nsmax),
     .  nkbmax(nsmax), cnfigtb(0:3,nsmax)
C
C*******************************************************************

C*******************************************************************
C
       common/control/ismax,nomax,nkbmax 
       common/cmcnfig/cnfigtb
C
C*******************************************************************

        if ((is.lt.1).or.(is.gt.ismax)) then
          write(6,*) 'CNFIGFL: THERE ARE NO DATA FOR IS=',IS
          write(6,*) 'CNFIGFL: ISMIN= 1, ISMAX= ',ismax
          STOP
        endif

          if(l.gt.3) then 
             cnfigfl=l+1
          else 
            cnfigfl=cnfigtb(l,is) 
          endif 


          return

          end




       DOUBLE PRECISION FUNCTION ATMPOPFL(IS,L)
C**********************************************************************
C Returns the populations of the valence l-shell 
C in the atomic ground state configuration
C Written by D.Sanchez-Portal. Aug, 1998
C************************INPUT******************************************
C    INTEGER  IS   : Species index
C    INTEGER  L    : Angular momentum of the shell
C************************OUTPUT*****************************************
C    REAL*8 ATMPOPFL: Population of the valence l-shell in the
C                     ground state 
C************************BEHAVIOUR**************************************
C  0) Before using ATMPOPFL, the pseudopotential must be initialized
C      by calling ATOM for each atomic species required.
C  1) Prints a message and stops when no data exist for IS and/or L
C***********************************************************************

         implicit none
         include 'atom.h'

         integer
     .      is,l


C************Variables in common blocks******************************
C

       integer
     .  ismax,nomax(nsmax),
     .  nkbmax(nsmax)
   
       double precision 
     .  qltb(0:3, nsmax)
C
C********************************************************************

C*******************************************************************
C
       common/control/ismax,nomax,nkbmax
       common/cmql/qltb
C
C*******************************************************************

        if ((is.lt.1).or.(is.gt.ismax)) then
          write(6,*) 'ATMPOPFL: THERE ARE NO DATA FOR IS=',IS
          write(6,*) 'ATMPOPFL: ISMIN= 1, ISMAX= ',ismax
          STOP
        endif 

        if(l.gt.3) then  
            atmpopfl=0.0d0
        else
            atmpopfl=qltb(l,is)
        endif

        return
        end
        


       INTEGER FUNCTION NZTFL (IS,L)
C**********************************************************************
C Returns the number of different basis functions
C with the same angular momentum and for a given specie
C Written by D.Sanchez-Portal. Aug, 1998
C************************INPUT******************************************
C    INTEGER  IS   : Species index
C    INTEGER  L    : Angular momentum of the basis functions
C************************OUTPUT*****************************************
C   INTEGER NZTFL  : Number of different basis functions with the same
C                    angular momentum.
C************************BEHAVIOUR**************************************
C  0) Before using NZTL, the pseudopotential must be initialized
C      by calling ATOM for each atomic species required.
C  1) Prints a message and stops when no data exist for IS and/or L
C***********************************************************************
       
         implicit none
         include 'atom.h'
        
         integer
     .      is,l


C************Variables in common blocks******************************
C

       integer
     .  nzetasave(0:lmaxd,nsmax),ismax,nomax(nsmax),
     .  nkbmax(nsmax),npolorbsave(lmaxd,nsmax),
     .  lmxosave(nsmax)
C
C*******************************************************************

C*******************************************************************
C
       common/cmzeta/nzetasave
       common/control/ismax,nomax,nkbmax
       common/cmpolorb/npolorbsave
       common/cmlmxo/lmxosave
C
C*******************************************************************


        if ((is.lt.1).or.(is.gt.ismax)) then 
          write(6,*) 'NZTFL: THERE ARE NO DATA FOR IS=',IS
          write(6,*) 'NZTLF: ISMIN= 1, ISMAX= ',ismax
          STOP
        endif
         
          if(l.gt.lmxosave(is)+1) then 
            nztfl=0
          elseif(l.eq.lmxosave(is)+1) then 
            if(lmaxd.ge.lmxosave(is)+1) then 
               nztfl=npolorbsave(l,is)  
            else
               nztfl=0
            endif 
          elseif(l.eq.0) then 
            nztfl=nzetasave(l,is)
          else
            nztfl=nzetasave(l,is)+npolorbsave(l,is)
          endif 
 
          return
         
          end 







          CHARACTER*(*) FUNCTION LABELFIS (IS)
C**********************************************************************
C Returns the atomic label for species is
C Written by D.Sanchez-Portal. Aug, 1998
C************************INPUT******************************************
C    INTEGER  IS   : Species index
C************************OUTPUT*****************************************
C    CHARACTER*(*) LABELFIS: Atomic label for species is
C************************BEHAVIOUR**************************************
C  0) Before using LABELFIS subroutine ATOM should be called 
C      for each atomic species required.
C  1) Prints a message and stops when no data exist for IS 
C***********************************************************************

         implicit none
         include 'atom.h'

         integer
     .      is,ismax, nomax(nsmax), nkbmax(nsmax) 
         character*20
     .      label_save(nsmax)

         common/control/ismax,nomax,nkbmax
         common/cmlabel/label_save

        if ((is.lt.1).or.(is.gt.ismax)) then
          write(6,*) 'LABELFIS: THERE ARE NO DATA FOR IS=',IS
          write(6,*) 'LABELFIS: ISMIN= 1, ISMAX= ',ismax
          STOP
        endif


          labelfis=label_save(is)


          return

          end




         INTEGER FUNCTION NOFIS(IS)
C**********************************************************************
C Returns the total number of basis functions for the species is
C Written by D.Sanchez-Portal. Aug, 1998
C************************INPUT******************************************
C    INTEGER  IS   : Species index
C************************OUTPUT*****************************************
C    INTEGER  NOFIS: Number of basis funxtions for species is
C************************BEHAVIOUR**************************************
C  0) Before using NOFIS subroutine ATOM should be called
C      for each atomic species required.
C  1) Prints a message and stops when no data exist for IS
C***********************************************************************

         implicit none
         include 'atom.h'

         integer
     .      is,ismax, nomax(nsmax), nkbmax(nsmax)

         common/control/ismax,nomax,nkbmax

        if ((is.lt.1).or.(is.gt.ismax)) then
          write(6,*) 'NOFIS: THERE ARE NO DATA FOR IS=',IS
          write(6,*) 'NOFIS: ISMIN= 1, ISMAX= ',ismax
          STOP
        endif


          nofis=nomax(is)


          return

          end



         INTEGER FUNCTION NKBFIS(IS)
C**********************************************************************
C Returns the total number of KB projectors for the species is
C Written by D.Sanchez-Portal. Aug, 1998
C************************INPUT******************************************
C    INTEGER  IS   : Species index
C************************OUTPUT*****************************************
C    INTEGER  NKBFIS: Number of KB projectors for species is
C************************BEHAVIOUR**************************************
C  0) Before using NKBFIS subroutine ATOM should be called
C      for each atomic species required.
C  1) Prints a message and stops when no data exist for IS
C***********************************************************************

         implicit none
         include 'atom.h'

         integer
     .      is,ismax, nomax(nsmax), nkbmax(nsmax)

         common/control/ismax,nomax,nkbmax

        if ((is.lt.1).or.(is.gt.ismax)) then
          write(6,*) 'NKBFIS: THERE ARE NO DATA FOR IS=',IS
          write(6,*) 'NKBFIS: ISMIN= 1, ISMAX= ',ismax
          STOP
        endif


          nkbfis=nkbmax(is)


          return

          end



         DOUBLE PRECISION FUNCTION MASSFIS(IS)
C**********************************************************************
C Returns the atomic mass for the species is
C Written by D.Sanchez-Portal. Aug, 1998
C************************INPUT******************************************
C    INTEGER  IS   : Species index
C************************OUTPUT*****************************************
C    REAL*8  MASSFIS: Atomic mass for species is
C************************BEHAVIOUR**************************************
C  0) Before using MASSFIS subroutine ATOM should be called
C      for each atomic species required.
C  1) Prints a message and stops when no data exist for IS
C***********************************************************************

         implicit none
         include 'atom.h'

         integer
     .      is,ismax, nomax(nsmax), nkbmax(nsmax)
         double precision 
     .      smasstb(nsmax) 
     
         common/control/ismax,nomax,nkbmax 
         common/cmmass/smasstb

        if ((is.lt.1).or.(is.gt.ismax)) then
          write(6,*) 'MASSFIS: THERE ARE NO DATA FOR IS=',IS
          write(6,*) 'MASSFIS: ISMIN= 1, ISMAX= ',ismax
          STOP
        endif


          massfis=smasstb(is)


          return

          end





          CHARACTER*(*) FUNCTION SYMFIO (IS,IO)
C********************************************************************** 
C Returns a label describing the symmetry of the
C   basis orbital or Kleynman-Bylander projector.
C Written by D.Sanchez-Portal. Aug., 1996
C************************INPUT******************************************
C    INTEGER  IS   : Species index
C    INTEGER  IO   : Orbital index (within atom)
C                    IO > 0 => Basis orbitals
C                    IO < 0 => Kleynman-Bylander projectors
C************************OUTPUT*****************************************
C   INTEGER SYMFIO  : Symmetry of the orbital or KB projector
C************************BEHAVIOUR**************************************
C  0) Before using SYMFIO, the pseudopotential must be initialized
C      by calling ATOM for each atomic species required.
C  1) Prints a message and stops when no data exist for IS and/or IO
C  2) Returns 's' for IO = 0
C***********************************************************************

        implicit none

       include 'atom.h'

       integer  is,io

C*************Internal variables*************************************
C
         integer lofio, mofio, lmax_sym, ilm, i, lorb, morb
         parameter(lmax_sym=3)
   
         character  sym_label((lmax_sym+1)*(lmax_sym+1))*6, 
     .        paste*7
         
         logical     pol      
C
C********************************************************************
C************Variables in common blocks******************************
C
       integer
     .  ismax, nomax(nsmax), nkbmax(nsmax)
C
C*******************************************************************

C*******************************************************************
C
       common/control/ismax,nomax,nkbmax
C
C*******************************************************************

C*******************External declarations***************************
C
              external lofio, mofio, pol, paste
C
C*******************************************************************

C*******************DATA*********************************************
C
        data  sym_label(1)          / 's' /
        data (sym_label(i),i=2,4)   / 'py', 'pz', 'px' /
        data (sym_label(i),i=5,9)   / 'dxy', 'dyz', 'dz2',
     .                                'dxz', 'dx2-y2' / 
        data (sym_label(i),i=10,16) / 'f', 'f', 'f', 'f', 
     .                                'f', 'f', 'f' /
C
C*******************************************************************




        if ((is.lt.1).or.(is.gt.ismax)) then
          write(6,*) 'SYMFIO: THERE ARE NO DATA FOR IS=',IS
          write(6,*) 'SYMFIO: ISMIN= 1, ISMAX= ',ismax
          STOP
        endif
        if((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
          write(6,*) 'SYMFIO: THERE ARE NO DATA FOR IO=',IO
          write(6,*) 'SYMFIO: IOMIN= ',-nkbmax(is),
     .     ' IOMAX= ',nomax(is)
          STOP
        endif

        lorb=lofio(is,io)
        morb=mofio(is,io)

        if(lorb.gt.lmax_sym ) then 
            symfio=' '
        else
            ilm=lorb*lorb+lorb+morb+1  
            if(pol(is,io)) then 
              symfio=paste('P',sym_label(ilm))
            else
              symfio=sym_label(ilm) 
            endif 
        endif         



        return
        end









        subroutine phiatm(is,io,r,phi,grphi)
C***********************************************************************
C  Returns Kleynman-Bylander local pseudopotential, nonlocal projectors,
C  and atomic basis orbitals (and their gradients).
C  Written by D.Sanchez-Portal. May, 1996. 
C  Modified August, 1998.
C*************************INPUT*****************************************
C  INTEGER IS : Species index
C  INTEGER IO : Orbital index (within atom):
C              IO > 0 =>  Basis orbitals
C              IO = 0 =>  Local screened pseudopotential
C              IO < 0 =>  Kleynman-Bylander projectors
C  REAL*8 R(3): Point vector, relative to atom 
C*************************OUTPUT****************************************
C  REAL*8 PHI     : Value of the basis orbital, KB-projector or local
C                    pseudopotential
C  REAL*8 GRPHI(3): Gradient of the basis orbital, KB-projector or local
C                    pseudopotential
C*************************UNITS*****************************************
C Distances in Bohr
C************************BEHAVIOUR**************************************
C 0) Before using PHIATM, the pseudopotential must be initialized 
C    by calling ATOM for each atomic species required
C 1) Each projector and basis function has a well defined total
C    angular momentum (quantum number l).
C 2) Basis functions are normalized and mutually orthogonal
C 3) Projection functions are normalized and mutually orthogonal
C 4) Normalization of KB projectors |Phi_lm> is such that 
C     <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                   Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C    where epsKB_l is returned by function EPSKB
C 5) Prints a message and stops when no data exits for IS and/or IO
C 6) Returns exactly zero when |R| > RCUT(IS,IO)
C 7) PHIATM with IO = 0 is strictly equivalent to VLOCAL
C*********************************************************************** 

        implicit none
        include 'atom.h'
        
        double precision
     .      r(3), phi,grphi(3)

        integer
     .      is,io

C*************Internal variables*************************************
C
         integer l, norb, lorb, izeta, ipol, nkb,
     .    indx, morb, ilm, i

         double precision  rly(lmx2), grly(3,lmx2), rmax, rmod,
     .      phir, dphidr, delt
         
         logical pol
C
C********************************************************************


C************Variables in common blocks******************************
C
       double precision
     .  table((ntbmax+2),-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tabpol((ntbmax+2),nzetmx*lmaxd,nsmax),
     .  tab2(ntbmax,-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tab2pol(ntbmax,nzetmx*lmaxd,nsmax)


       integer
     .  nzetasave(0:lmaxd,nsmax),ismax,nomax(nsmax),
     .  nkbmax(nsmax),npolorbsave(lmaxd,nsmax),
     .  lmxosave(nsmax), lmxkbsave(nsmax)
C
C*******************************************************************

C*******************************************************************
C
       common/cmtab/table,tabpol 
       common/cmspline/tab2,tab2pol
       common/cmzeta/nzetasave
       common/control/ismax,nomax,nkbmax
       common/cmpolorb/npolorbsave
       common/cmlmxo/lmxosave
       common/cmlmxkb/lmxkbsave
C
C*******************************************************************


        if ((is.lt.1).or.(is.gt.ismax)) then
           write(6,*) 'PHIATM: THERE ARE NO DATA FOR IS=',IS
           write(6,*) 'PHIATM: ISMIN= 1, ISMAX= ',ismax
           STOP
        endif
        if((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
          write(6,*) 'PHIATM: THERE ARE NO DATA FOR IO=',IO
          write(6,*) 'PHIATM: IOMIN= ',-nkbmax(is),
     .     ' IOMAX= ',nomax(is)
          STOP
        endif

       pol=.false.
       if (io.gt.0) then

          norb=0 
          indx=0
          do  l=0,lmxosave(is)
            do izeta=1,nzetasave(l,is)
               norb=norb+(2*l+1)
               indx=indx+1
               if(norb.ge.io) then 
                   lorb=l
                   morb=io-norb+lorb
                   goto 20 
               endif 
            enddo 
          enddo 

          indx=0
          do  l=1,min(lmxosave(is)+1,lmaxd)
            do ipol=1, npolorbsave(l,is)
              norb=norb+(2*l+1)
              indx=indx+1
              if(norb.ge.io) then  
                    lorb=l
                    morb=io-norb+lorb
                    pol=.true.
                    goto 20 
             endif   
           enddo 
          enddo  

20       continue

         elseif(io.lt.0) then
         nkb=0
         indx=0
         do l=0,lmxkbsave(is)
            indx=indx+1
            nkb=nkb-(2*l+1)
            if(nkb.le.io) goto 30
         enddo 

 30      lorb=l
         morb=-io+nkb+lorb 
         indx=-indx
      
         else

             indx=0 
C            Next two lines introduced by J.Soler. 2/7/96.
             lorb=0
             morb=0

         endif 



 




        if(.not.pol) then 
           delt=table(1,indx,is)  
        else
           delt=tabpol(1,indx,is)
        endif 

        rmax=delt*(ntbmax-1)
  
        rmod=0.0d0
        do i=1,3
           rmod=rmod+r(i)*r(i)
        enddo
        rmod=dsqrt(rmod)+1.0d-20

        if(rmod.gt.rmax-1.d-12) then

           phi=0.0d0
           grphi(1)=0.0d0
           grphi(2)=0.0d0
           grphi(3)=0.0d0

        else
        
           if(.not.pol) then
               call splint(delt,table(3,indx,is),
     .               tab2(1,indx,is),ntbmax,
     .               rmod,phir,dphidr)
           else
               call splint(delt,tabpol(3,indx,is),
     .               tab2pol(1,indx,is),ntbmax,
     .               rmod,phir,dphidr)
           endif

           if (indx.eq.0) then
              phi=phir
              do i=1,3 
                 grphi(i)=dphidr*r(i)/rmod
              enddo
           else

              ilm = lorb*lorb + lorb + morb + 1
              call rlylm( lorb, r, rly, grly )

              phi = phir * rly(ilm)
              
              do i = 1,3
                grphi(i)=dphidr*rly(ilm)*r(i)/rmod+phir*grly(i,ilm)
              enddo


*             write(6,'(a,i4,2f12.6)')
*    .         'phiatm: ilm,phi/rl,rl*ylm=', ilm, phi, rly(ilm)

           endif
                   
        endif


  
        return
        
        end







       double precision function rcut(is,io)




C**************************************************************************
C  Returns cutoff radius of Kleynman-Bylander projectors and
C  atomic basis orbitals.
C  To be written by D.Sanchez-Portal. May 1996.
C***********************INPUT**********************************************
C  INTEGER IS : Species index
C  INTEGER IO : Orbital index (within atom):
C                IO > 0 =>  Basis orbitals
C                IO < 0 =>  Kleynman-Bylander projectors
C                IO = 0 =>  Local screened pseudopotential
C***********************OUTPUT********************************************* 
C  REAL*8  RCUT : Cutoff radius
C***********************UNITS**********************************************
C  Distances in Bohr
C***********************BEHAVIOUR******************************************
C   0) Before using RCUT, the pseudopotential must be initialized
C      by calling ATOM for each atomic species required.
C   1) Prints a message and stops when no data exits for IS and/or IO
C**************************************************************************



       implicit none

       include 'atom.h'
       
       integer  is,io

C*************Internal variables*************************************
C
         integer l, norb, lorb, nzetorb, izeta, ipol, nkb,lkb
         
         integer  indx        


C
C******************************************************************** 
C************Variables in common blocks******************************
C
       double precision
     .  rctb(0:lmaxd,nsmax),rcotb(nzetmx,0:lmaxd,nsmax),
     .  table((ntbmax+2),-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tabpol((ntbmax+2),nzetmx*lmaxd,nsmax),
     .  rcpoltb(nzetmx,lmaxd,nsmax)

       integer 
     .  nzetasave(0:lmaxd,nsmax),ismax,nomax(nsmax),
     .  nkbmax(nsmax),npolorbsave(lmaxd,nsmax),
     .  lmxosave(nsmax), lmxkbsave(nsmax)
C
C*******************************************************************

C*******************************************************************
C
       common/cmradkb/rctb
       common/cmradorb/rcotb  
       common/cmradpol/rcpoltb
       common/cmtab/table,tabpol
       common/cmzeta/nzetasave
       common/control/ismax,nomax,nkbmax 
       common/cmpolorb/npolorbsave
       common/cmlmxo/lmxosave 
       common/cmlmxkb/lmxkbsave
C
C*******************************************************************

       if ((is.lt.1).or.(is.gt.ismax)) then
         write(6,*) 'RCUT: THERE ARE NO DATA FOR IS=',IS
         write(6,*) 'RCUT: ISMIN= 1, ISMAX= ',ismax
         STOP
       endif
       if((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
         write(6,*) 'RCUT: THERE ARE NO DATA FOR IO=',IO
         write(6,*) 'RCUT: IOMIN= ',-nkbmax(is),
     .     ' IOMAX= ',nomax(is)
         STOP
       endif


       if (io.gt.0) then

        norb=0 
        indx=0
        do 10 l=0,lmxosave(is)
          do 5 izeta=1,nzetasave(l,is)
            norb=norb+(2*l+1)
            indx=indx+1 
            if(norb.ge.io) goto 30 
 5        continue       
10      continue   

        indx=0
        do  20 l=1,lmxosave(is)+1
            do 15 ipol=1, npolorbsave(l,is)
              norb=norb+(2*l+1)  
              indx=indx+1  
              if(norb.ge.io) goto 40
15          continue
20      continue     
        write(6,*) 'RCUT: ERROR: ORBITAL INDEX IO=',IO
        write(6,*) 'RCUT: ERROR: NOT FOUND'
        stop

30      lorb=l
        nzetorb=izeta
        rcut=rcotb(nzetorb,lorb,is)  
C       rcut=table(1,indx,is)*(ntbmax-1)
        return  

40      lorb=l
        nzetorb=ipol
        rcut=rcpoltb(nzetorb,lorb,is) 
C       rcut=tabpol(1,indx,is)*(ntbmax-1)
        return 


       elseif(io.lt.0) then 


        nkb=0
        do 50 l=0,lmxkbsave(is)
          nkb=nkb-(2*l+1)
          if(nkb.le.io) goto 60
50      continue 

60      lkb=l 
        indx=-(lkb+1)
        rcut=rctb(lkb,is) 
c       rcut=table(1,indx,is)*(ntbmax-1)
        return 

        elseif (io.eq.0) then


        rcut=table(2,0,is)

        endif

        return 

        end









       double precision function rcore(is)




C**************************************************************************
C  Returns cutoff radius of the pseudo-core charge density for the non-linear
C   core corrections for xc potential.
C  To be written by D.Sanchez-Portal. Oct 1996.
C***********************INPUT**********************************************
C  INTEGER IS : Species index
C***********************OUTPUT********************************************* 
C  REAL*8  RCORE  : Cutoff radius
C***********************UNITS**********************************************
C  Distances in Bohr
C***********************BEHAVIOUR******************************************
C   0) Before using RCORE, the pseudopotential must be initialized
C      by calling ATOM for each atomic species required.
C   1) Prints a message and stops when no data exits for IS 
C**************************************************************************



       implicit none
       include 'atom.h'

       double precision
     .  coretab(ntbmax+1,2,nsmax)

       integer 
     .  is,ismax, nomax(nsmax), nkbmax(nsmax)

       common/control/ismax, nomax, nkbmax
       common/cmcore/coretab

       if ((is.lt.1).or.(is.gt.ismax)) then
         write(6,*) 'RCORE: THERE ARE NO DATA FOR IS=',IS
         write(6,*) 'RCORE: ISMIN= 1, ISMAX= ',ismax
         STOP
       endif

        rcore=coretab(1,1,is)*dble(ntbmax-1)

        return 

        end








          subroutine chcore(is,r,ch,grch)


C***********************************************************************
C Returns returns pseudo-core charge density for non-linear core correction
C in the xc potential.
C Written by D.Sanchez-Portal. May 1996
C************************INPUT******************************************
C INTEGER IS     :  Species index
C REAL*8  R(3)   :  Point vector, relative to atom
C ***********************OUTPUT*****************************************
C REAL*8  CHC    :  Value of pseudo-core charge density.
C REAL*8  GRV(3) :  Gradient of pseudo-core charge density.
C************************UNITS******************************************
C Distances in Bohr
C Energies in Rydbergs
C Density in electrons/Bohr**3
C************************BEHAVIOUR**************************************
C  0) Before using CHCORE, the pseudopotential must be initialized
C     by calling ATOM for each atomic species required
C  1) Prints a message and stops when no data exits for IS.
C  2) Returns exactly zero when |R| > Rcore
C***********************************************************************
 
          implicit double precision (a-h,o-z)
          
          include 'atom.h'

          double precision 
     .     r(3),ch,grch(3),coretab(ntbmax+1,2,nsmax)
 

          integer
     .     is,ismax, nomax(nsmax), nkbmax(nsmax)
 
          common/control/ismax, nomax, nkbmax
          common/cmcore/coretab



          if ((is.lt.1).or.(is.gt.ismax)) then
            write(6,*) 'CHCORE: THERE ARE NO DATA FOR IS=',IS
            write(6,*) 'CHCORE: ISMIN= 1, ISMAX= ',ismax
            STOP
          endif

          core=coretab(1,2,is)
          if(core.eq.0) then 
            ch=0.0d0 
            grch(1)=0.0d0 
            grch(2)=0.0d0
            grch(3)=0.0d0 
            return
          endif

          delt=coretab(1,1,is) 
          rcmx=delt*(ntbmax-1)
          
          rmod=0.0d0
          do i=1,3
            rmod=rmod+r(i)*r(i)
          enddo
          rmod=dsqrt(rmod)
           
          if(rmod.gt.rcmx-1.d-12) then
             ch=0.0d0
             grch(1)=0.0d0
             grch(2)=0.0d0
             grch(3)=0.0d0
        
          else
            call splint(delt,coretab(2,1,is),coretab(2,2,is),ntbmax,
     .        rmod,ch,dchdr)

             rmod=rmod+1.0d-20
           
             grch(1)=dchdr*r(1)/rmod
             grch(2)=dchdr*r(2)/rmod
             grch(3)=dchdr*r(3)/rmod
       

          endif

 
          return
 
          end 


          subroutine psover(is1,is2,r,energ,dedr)

C***********************************************************************
C Returns electrostatic correction to the ions interaction energy
C due to the overlap of the two 'local pseudopotential charge densities'
C Written by D.Sanchez-Portal. March, 1997
C************************INPUT******************************************
C INTEGER IS     :  Species index
C REAL*8  R      :  Distance between atoms.
C ***********************OUTPUT*****************************************
C REAL*8  ENERG  :  Value of the correction interaction energy.
C REAL*8  DEDR   :  Radial derivative of the correction.
C************************UNITS******************************************
C Distances in Bohr
C Energies in Rydbergs
C************************BEHAVIOUR**************************************
C  0) Before using PSOVER, the pseudopotential must be initialized
C     by calling ATOM for each atomic species required
C  1) Prints a message and stops when no data exits for IS.
C  2) Returns exactly zero when |R| > Rchloc
C***********************************************************************
 
          implicit double precision (a-h,o-z)
          
          include 'atom.h'
          parameter (ns2=((nsmax+1)*nsmax)/2)

          double precision 
     .     r,energ,dedr,corrtab((ntbmax+1),2,ns2)
 

          integer
     .     ismax, nomax(nsmax), nkbmax(nsmax)
 
          common/control/ismax, nomax, nkbmax
          common/cmcorr/corrtab



          if ((is1.lt.1).or.(is1.gt.ismax)) then
            write(6,*) 'PSOVER: THERE ARE NO DATA FOR IS1=',IS1
            write(6,*) 'PSOVER: ISMIN= 1, ISMAX= ',ismax
            STOP
          endif
          if ((is2.lt.1).or.(is2.gt.ismax)) then
            write(6,*) 'PSOVER: THERE ARE NO DATA FOR IS2=',IS2
            write(6,*) 'PSOVER: ISMIN= 1, ISMAX= ',ismax
            STOP
          endif
         
          ismx=max(is1,is2)
          ismn=min(is1,is2)
    
          indx=((ismx-1)*ismx)/2 + ismn

          dloc=corrtab(1,2,indx)
          if(dabs(dloc).lt.1.0d-8) then 
            energ=0.0d0 
            dedr=0.0d0 
            return
          endif

          delt=corrtab(1,1,indx) 
          rcmx=delt*(ntbmax-1)
          if(r.gt.rcmx-1.d-12) then

             energ=0.0d0
             dedr=0.0d0
        
          else


            call splint(delt,corrtab(2,1,indx),corrtab(2,2,indx),ntbmax,
     .        r,energ,dedr)
         
             r=r+1.0d-20
             
             energ=2.0d0*energ/r
             dedr=(-energ + 2.0d0*dedr)/r 

          endif

 
          return
 
          end




         DOUBLE PRECISION FUNCTION UION ( IS )

C**************************************************************************
C  Returns electrostatic self-energy of the 'ions', assigned by routine PSEUDO
C  Written by D. Sanchez-Portal. May, 1996 
C**************************INPUT*******************************************
C  INTEGER IS     : Species index
C**************************OUTPUT******************************************
C  INTEGER UION   : Self-energy of the 'ion'
C***************************UNITS*******************************************
C   Energy in Rydbergs
C**************************BEHAVIOUR***************************************
C  0) Before using UION, the pseudopotential must be initialized 
C     by calling ATOM for each atomic species required.
C  1) Prints a message and stops when no data exists for IS.
C**************************************************************************

         include 'atom.h'

         double precision 
     .    slfe(nsmax)
 
         integer
     .    is,ismax, nomax(nsmax), nkbmax(nsmax)

          

         common/cmslfe/slfe
         common/control/ismax, nomax, nkbmax

         if ((is.lt.1).or.(is.gt.ismax)) then 
            write(6,*) 'UION: THERE ARE NO DATA FOR IS=',IS
            write(6,*) 'UION: ISMIN= 1, ISMAX= ',ismax
            STOP
         endif


 
         uion=slfe(is)


         return

         end








       DOUBLE PRECISION FUNCTION EPSKB (IS,IO)



C**************************************************************************
C  Returns the energies epsKB_l of the Kleynman-Bylander projectors:
C       <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                 Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C  where Phi_lm is returned by subroutine PHIATM.
C  Written by D.Sanchez-Portal.  May,  1996.
C************************INPUT*********************************************
C  INTEGER IS    : Species index
C  INTEGER IO    : Kleynman-Bylander projector index (within atom).
C                  May be positive or negative (only ABS(IO) is used).
C***********************OUTPUT*********************************************
C  REAL*8 EPSKB  : Kleynman-Bylander projector energy
C***********************UNITS**********************************************
C  Energy in Rydbergs.
C***********************BEHAVIOUR******************************************
C  0) Before using EPSKB, the pseudopotential must be initialized by 
C    calling ATOM for each atomic species required.
C  1) Prints a message and stops when no data exist for IS,IO.
C**************************************************************************
 
        implicit none
        include 'atom.h'

        integer is, io

C*****************Internal variables******************************
C
        integer  ionew, nkb, indx, l
C
C******************************************************************
C*****************Variables in the common blocks*******************
C
        double precision
     .   table((ntbmax+2),-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .   tabpol((ntbmax+2),nzetmx*lmaxd,nsmax)


        integer
     .   ismax,nomax(nsmax),nkbmax(nsmax),lmxkbsave(nsmax)

         common/cmtab/table,tabpol
         common/control/ismax,nomax,nkbmax
         common/cmlmxkb/lmxkbsave
C
C****************************************************************

         if ((is.lt.1).or.(is.gt.ismax)) then
            write(6,*) 'EPSKB: THERE ARE NO DATA FOR IS=',IS
            write(6,*) 'EPSKB: ISMIN= 1, ISMAX= ',ismax
            STOP
         endif
 
         ionew=-abs(io)
         if (ionew.eq.0) then 
            write(6,*) 'EPSKB: FUNCTION CANNOT BE CALLED WITH'
     .       ,' ARGUMENT EQUAL TO ZERO' 
              STOP 
         endif

         if(ionew.lt.-nkbmax(is)) then
            write(6,*) 'EPSKB: THERE ARE NO DATA FOR IO=',IONEW
            write(6,*) 'EPSKB: IOMIN= ',-nkbmax(is)
            STOP
         endif

         nkb=0
         indx=0
         do 10  l=0,lmxkbsave(is)
             indx=indx+1
             nkb=nkb-(2*l+1)
             if(nkb.le.ionew) goto 20
10       continue 

20        indx=-indx
 

 
        epskb=table(2,indx,is)


        return 
 
        end









        subroutine rphiatm(is,io,r,phi,dphidr)
C***********************************************************************
C  Returns the radial component of 
C  Kleynman-Bylander local pseudopotential, nonlocal projectors,
C  and atomic basis orbitals (and their radial drivatives)
C  
C*************************INPUT*****************************************
C  INTEGER IS : Species index
C  INTEGER IO : Type of radial function to plot:
C              IO > 0 =>  Basis orbitals
C              IO = 0 =>  Local screened pseudopotential
C              IO < 0 =>  Kleynman-Bylander projectors 
C  REAL*8 R   : Radial distance, relative to atom 
C*************************OUTPUT****************************************
C  REAL*8 PHI     : Value of the basis orbital, KB-projector or local
C                    pseudopotential
C  REAL*8 DPHIDR:  Radial derivative of the basis orbital, 
C                    KB-projector or local  pseudopotential
C*************************UNITS*****************************************
C Distances in Bohr
C************************BEHAVIOUR**************************************
C 0) Before using RPHIATM, the pseudopotential must be initialized 
C    by calling ATOM for each atomic species required
C 1) Each projector and basis function has a well defined total
C    angular momentum (quantum number l).
C 2) Basis functions are normalized and mutually orthogonal
C 3) Projection functions are normalized and mutually orthogonal
C 4) Normalization of KB projectors |Phi_lm> is such that 
C     <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                   Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C    where epsKB_l is returned by function EPSKB
C 5) Prints a message and stops when no data exits for IS and/or IO
C 6) Returns exactly zero when |R| > RCUT(IS,IO)
C 7) RPHIATM with ITYPE = 0 is strictly equivalent to VLOCAL
C*********************************************************************** 

        implicit none
        include 'atom.h'
        
        double precision
     .      r, phi,dphidr

        integer
     .      is,io

C*************Internal variables*************************************
C
         integer l, norb, lorb, izeta, ipol, nkb,
     .    indx

         double precision  rmax, rmod, phir, delt
         
         logical pol
C
C********************************************************************


C************Variables in common blocks******************************
C
       double precision
     .  table((ntbmax+2),-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tabpol((ntbmax+2),nzetmx*lmaxd,nsmax),
     .  tab2(ntbmax,-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tab2pol(ntbmax,nzetmx*lmaxd,nsmax)


       integer
     .  nzetasave(0:lmaxd,nsmax),ismax,nomax(nsmax),
     .  nkbmax(nsmax),npolorbsave(lmaxd,nsmax),
     .  lmxosave(nsmax), lmxkbsave(nsmax)
C
C*******************************************************************

C*******************************************************************
C
       common/cmtab/table,tabpol 
       common/cmspline/tab2,tab2pol
       common/cmzeta/nzetasave
       common/control/ismax,nomax,nkbmax
       common/cmpolorb/npolorbsave
       common/cmlmxo/lmxosave
       common/cmlmxkb/lmxkbsave
C
C*******************************************************************


        if ((is.lt.1).or.(is.gt.ismax)) then
           write(6,*) 'RPHIATM: THERE ARE NO DATA FOR IS=',IS
           write(6,*) 'RPHIATM: ISMIN= 1, ISMAX= ',ismax
           STOP
        endif
         
        if((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
          write(6,*) 'RPHIATM: THERE ARE NO DATA FOR IO=',IO
          write(6,*) 'RPHIATM: IOMIN= ',-nkbmax(is),
     .     ' IOMAX= ',nomax(is)
          STOP
        endif

       pol=.false.
       if (io.gt.0) then

          norb=0 
          indx=0
          do  l=0,lmxosave(is)
            do izeta=1,nzetasave(l,is)
               norb=norb+(2*l+1)
               indx=indx+1
               if(norb.ge.io) then 
                   lorb=l
                   goto 20 
               endif 
            enddo
          enddo 

          indx=0
          do  l=1,min(lmxosave(is)+1,lmaxd)
            do ipol=1, npolorbsave(l,is)
              norb=norb+(2*l+1)
              indx=indx+1
              if(norb.ge.io) then  
                    lorb=l
                    pol=.true.
                    goto 20 
             endif   
           enddo 
          enddo  

20       continue

         elseif(io.lt.0) then
         nkb=0
         indx=0
         do l=0,lmxkbsave(is)
            indx=indx+1
            nkb=nkb-(2*l+1)
            if(nkb.le.io) goto 30
         enddo 

 30      lorb=l
         indx=-indx
      
         elseif (io.eq.0) then

             indx=0 
C            Next two lines introduced by J.Soler. 2/7/96.
             lorb=0
         endif 

        if(.not.pol) then 
           delt=table(1,indx,is)  
        else
           delt=tabpol(1,indx,is)
        endif 

        rmax=delt*(ntbmax-1)
  
        rmod=r+1.0d-20

        if(rmod.gt.rmax-1.d-12) then

           phi=0.0d0
           dphidr=0.0d0
        else
        
           if(.not.pol) then
               call splint(delt,table(3,indx,is),
     .               tab2(1,indx,is),ntbmax,
     .               rmod,phir,dphidr)
           else
               call splint(delt,tabpol(3,indx,is),
     .               tab2pol(1,indx,is),ntbmax,
     .               rmod,phir,dphidr)
           endif

            if (l.eq.0) then
               phi=phir
            elseif (l.eq.1) then
               phi=phir*r
               dphidr=dphidr*r
               dphidr=dphidr+phir 
            else
               phi=phir*r**l 
               dphidr=dphidr*r**l
               dphidr=dphidr+l*r**(l-1)*phir 
            endif
 
                   
        endif


  
        return
        
        end 



       LOGICAL FUNCTION POL (IS,IO)
C********************************************************************** 
C If true, the orbital IO is a perturbative polarization orbital
C Written by D.Sanchez-Portal. Oct, 1996
C************************INPUT******************************************
C    INTEGER  IS   : Species index
C    INTEGER  IO   : Orbital index (within atom)
C                    IO > 0 => Basis orbitals
C************************OUTPUT*****************************************
C    LOGICAL POL   : if true, this is a perturbative polarization orbital
C************************BEHAVIOUR**************************************
C  0) Before using POL, the basis set must be initialized
C      by calling ATOM for each atomic species required.
C  1) Prints a message and stops when no data exist for IS and/or IO
C  2) Returns zero for IO = 0
C***********************************************************************
       

       implicit none

       include 'atom.h'

       integer  is,io

C*************Internal variables*************************************
C
         integer l, norb, izeta, ipol



C
C********************************************************************
C************Variables in common blocks******************************
C
       integer
     .  nzetasave(0:lmaxd,nsmax),ismax,nomax(nsmax),
     .  nkbmax(nsmax),npolorbsave(lmaxd,nsmax),
     .  lmxosave(nsmax)
C
C*******************************************************************

C*******************************************************************
C
       common/cmzeta/nzetasave
       common/control/ismax,nomax,nkbmax
       common/cmpolorb/npolorbsave
       common/cmlmxo/lmxosave
C
C*******************************************************************
  
        if ((is.lt.1).or.(is.gt.ismax)) then 
          write(6,*) 'POL: THERE ARE NO DATA FOR IS=',IS
          write(6,*) 'POL: ISMIN= 1, ISMAX= ',ismax
          STOP
        endif
        if(io.gt.nomax(is)) then 
          write(6,*) 'POL: THERE ARE NO DATA FOR IO=',IO
          write(6,*) 'POL: IOMIN= 0',
     .     ' IOMAX= ',nomax(is)
          STOP
        endif
 

        norb=0
        do 10 l=0,lmxosave(is)
          do 5 izeta=1,nzetasave(l,is)
            norb=norb+(2*l+1)
            if(norb.ge.io) goto 30
 5        continue
10      continue

        do  20 l=1,min(lmxosave(is)+1,lmaxd)
            do 15 ipol=1, npolorbsave(l,is)
              norb=norb+(2*l+1)
              if(norb.ge.io) goto 40
15          continue
20      continue
        write(6,*) 'POL: ERROR: ORBITAL INDEX IO=',IO
        write(6,*) 'POL: ERROR: NOT FOUND'
        stop

30      pol=.false.
        return

40      pol=.true.
        return
         
          end 



