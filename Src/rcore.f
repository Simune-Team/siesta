






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



       implicit double precision(a-h,o-z)

       include 'atom.h'

       double precision
     .  coretab(ntbmax+1,2,nsmax)

       integer 
     .  is,ismax

       common/control/ismax
       common/cmcore/coretab

       if ((is.lt.1).or.(is.gt.ismax)) then
         write(6,*) 'RCORE: THERE ARE NO DATA FOR IS=',IS
         write(6,*) 'RCORE: ISMIN= 1, ISMAX= ',ismax
         STOP
       endif

        rcore=coretab(1,1,is)*dble(ntbmax-1)

        return 

        end


