




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
 
        implicit double precision (a-h,o-z)


        include 'atom.h'


        double precision
     .   table((ntbmax+2),-(lmaxd+1):nzetmx*(lmaxd+1),nsmax)

        integer
     .   is,io,ismax,nomax(nsmax),nkbmax(nsmax),ionew,
     .   loctab(nsmax) 

         common/cmtab/table
         common/control/ismax,nomax,nkbmax
         common/cmloc/loctab


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

         lloc=loctab(is) 
         nkb=0
         indx=0
         do 10  l=0,lmaxd
             if(l.eq.lloc) goto 10
             indx=indx+1
             nkb=nkb-(2*l+1)
             if(nkb.le.ionew) goto 20
10       continue 

20      lorb=l
        morb=-ionew+nkb+lorb
        indx=-indx
 

 
        epskb=table(2,indx,is)


        return 
 
        end
