





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
     .     is,ismax
 
          common/control/ismax
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
