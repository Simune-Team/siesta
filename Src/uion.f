



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
     .    is,ismax

          

         common/cmslfe/slfe
         common/control/ismax

         if ((is.lt.1).or.(is.gt.ismax)) then 
            write(6,*) 'UION: THERE ARE NO DATA FOR IS=',IS
            write(6,*) 'UION: ISMIN= 1, ISMAX= ',ismax
            STOP
         endif


 
         uion=slfe(is)


         return

         end



