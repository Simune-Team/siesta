






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



       implicit double precision(a-h,o-z)

       include 'atom.h'

       double precision
     .  rctb(0:lmaxd,nsmax),rcotb(nzetmx,0:lmaxd,nsmax),
     .  table((ntbmax+2),-(lmaxd+1):nzetmx*(lmaxd+1),nsmax)

       integer 
     .  is,io,nzettb(0:lmaxd,nsmax),ismax,nomax(nsmax),
     .  nkbmax(nsmax),loctab(nsmax)

       common/cmradkb/rctb
       common/cmradorb/rcotb
       common/cmzeta/nzettb
       common/cmtab/table
       common/control/ismax,nomax,nkbmax
       common/cmloc/loctab

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
        do 10 l=0,lmaxd
          do 5 izeta=1,nzettb(l,is)
            norb=norb+(2*l+1)
            if(norb.ge.io) goto 20
 5        continue       
10      continue 

20      lorb=l
        nzetorb=izeta

        rcut=rcotb(nzetorb,lorb,is)
       
       elseif(io.lt.0) then
        lloc=loctab(is)
        nkb=0
        do 30 l=0,lmaxd
          if(lloc.eq.l) goto 30
          nkb=nkb-(2*l+1)
          if(nkb.le.io) goto 40
30      continue 

40      lkb=l
     
        rcut=rctb(lkb,is)

        elseif (io.eq.0) then


        rcut=table(2,0,is)

        endif

        return 

        end


