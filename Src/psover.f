C $Id: psover.f,v 1.3 1999/01/31 11:20:15 emilio Exp $

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
     .     ismax
 
          common/control/ismax
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
