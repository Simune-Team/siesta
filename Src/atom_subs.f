C $Id: atom_subs.f,v 1.23 1999/02/26 21:13:41 daniel Exp $

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    This file contains all the routines called by subroutine 'atom'
C    for the generation of the basis set orbitals, Kleinman-Bylander 
C    projectors, local pseudopotential, neutral-atom pseudopotential,
C    etc... and also used to create the common blocks for the storage 
C    and later use by SIESTA of all this information. 
C
C    D. Sanchez-Portal  1996-1998
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C * The routines contained in this file are:
C
c     subroutine rc_vs_e
c     subroutine polarization
c     subroutine parabola
c     subroutine findp
c     subroutine nrmpal
c     subroutine radii_ps
c     subroutine vlocal2
c     subroutine vlocal1
c     subroutine schro_eq
c     subroutine ghost
c     subroutine KBproj
c     subroutine xc_check
c     subroutine comcore
c     subroutine comlocal
c     subroutine new_specie
c     subroutine read_vps
c     subroutine comKB
c     subroutine KBgen
c     subroutine Basis_gen
c     subroutine SPLIT
c     subroutine NODES
c     subroutine NONODES
c     subroutine comBasis
c     subroutine SPLITGAUSS
c     subroutine compress_PAO
c     subroutine atm_pop
c     subroutine Vna
c     subroutine comVna
c     subroutine slfe_local
c     subroutine POLgen
c     subroutine comPOL
c     subroutine set_mesh
c     subroutine BESSEL
c     subroutine draw_basis
c     subroutine USER
c     subroutine prinput
c     SUBROUTINE CHOVERLP
c     SUBROUTINE NUMEROV
c     CHARACTER*2 FUNCTION SYMBOL
c     CHARACTER*(*) FUNCTION PASTE
c     CHARACTER*(*) FUNCTION PASTEB
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




      CHARACTER*2 FUNCTION SYMBOL( IZ )
C RETURNS THE SYMBOL OF THE ELEMENT OF ATOMIC NUMBER IZ
C Written by J. Soler

      PARAMETER (NZ=103)
      CHARACTER*2 NAME(NZ)
      DATA NAME /'H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne',
     .           'Na','Mg','Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca',
     .           'Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn',
     .           'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y' ,'Zr',
     .           'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     .           'Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd',
     .           'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     .           'Lu','Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg',
     .           'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     .           'Pa','U' ,'Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     .           'Md','No','Lr'/
      IF (IZ.LT.1 .OR. IZ.GT.NZ) THEN
         WRITE(6,*) ' SYMBOL: OUT OF RANGE IZ =',IZ
         SYMBOL = ' '
      ELSE
         SYMBOL = NAME(IZ)
      ENDIF
      END





        SUBROUTINE CHOVERLP(IS1,IS2,RMX,CORR,CORR2,AUX)



C***********************************************************************
C  Returns a table with the difference between the electrostatic energy 
C  of two spherical charge-densities and two punctual charges with the 
C  same total charge as a function of the distance between the centers 
C  of these charge densities. 
C  Written by D.Sanchez-Portal. March, 1997.(from routine MATEL, written 
C  by Jose M. Soler)
C*************************INPUT*****************************************
C  INTEGER IS1,IS2             :  Species indexes.
C  EXTERNAL PSCH(IS,R,CH,GRCH) :  This functions provides the
C                                 'local-pseudopotential charge density'
C                                 for each species and its gradient
C                                 at point R.
C  RMX                         :  Maximum range of the correction.
C*************************OUTPUT****************************************
C  CORR(NTBMAX)                :  Electrostatic correction energy.
C  CORR2(NTBMAX)               :  Table with the second derivative 
C                                 of CORR for spline interpolation.
C  RMX                         :  Rmx is zero in output is one of 
C                                 the 'atoms' is not an atom but 
C                                 just a floating basis set. 
C*************************UNITS*****************************************
C Distances in Bohr.
C Energy in Rydbergs.
C***********************************************************************



C****************PRECISION**********************************************
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
C
C***********************************************************************


C******************INTERNAL PARAMETERS**********************************
C
C Internal precision parameters  ------------------------------------
C NQ is the number of radial points in reciprocal space.
C Npoint , 2npoint+1 is the number of points used by RATINT in the 
C interpolation.
C Q2CUT is the required planewave cutoff for the expansion of
C the 'local-pseudopotential atomic charge density'
C  (in Ry if lengths are in Bohr).
C TINY is a small number to check the precision of the charge density
C integration.
C
C***********************************************************************

        INCLUDE    'atom.h'
        PARAMETER ( NQ     =  512  )
        PARAMETER ( NPOINT =  4     ) 
        PARAMETER ( Q2CUT  =  2.5D3 )
        PARAMETER ( CHERR   =  5.D-2 )

C
C**********************************************************************

C********************ARRAYS DECLARATION********************************
C
        DIMENSION 
     .    CH(0:NQ,2),VTB(NTBMAX,2),CORR(NTBMAX),
     .    CORR2(NTBMAX),AUX(NTBMAX),V(0:NQ,2),
     .    GRCH(3),RX(3),RAUX(2*NPOINT+1)

        EXTERNAL PSCH
C
C**********************************************************************

C***********************PI********************************************* 
C
          PI= 4.D0 * ATAN(1.D0)       
          CONS= 1.0d0/(2.0d0*PI)**1.5D0
C
C***********************************************************************

C******************CUT-OFF IN REAL AND RECIPROCAL SPACE*****************
C
           QMAX =  DSQRT( Q2CUT )
           RMAX = PI * NQ / QMAX
           IF(RMX.GT.RMAX) THEN  
              WRITE(6,*) 'CHOVERLP: THE NUMBER OF INTEGRATION',
     .       ' POINTS MUST BE INCREASED'
              write(6,'(a,2f15.6)') 'chovrlap: rmx,rmax =', rmx, rmax
              STOP
           ENDIF 
           DELT=PI/QMAX
           C=4.0D0*PI*DELT
           DLT=RMX/(NTBMAX-1)
C
C*********************************************************************** 
     

C*********RADIAL CHARGE DENSITIES(CHECKING TOTAL CHARGE)****************
C
          IZ1=IZOFIS(IS1)
          IZ2=IZOFIS(IS2)

          IF((IZ1.LT.0.0D0).OR.(IZ2.LT.0.0D0)) THEN 
              DO ITB=1,NTBMAX
                 CORR(ITB)=0.0D0
                 CORR2(ITB)=0.0D0
              ENDDO 
              RMX=0.0D0
              RETURN
          ENDIF

          IZ1=IZVALFIS(IS1)
          IZ2=IZVALFIS(IS2)



          Z1=0.0D0
          Z2=0.0D0

          RX(2)=0.0D0
          RX(3)=0.0D0 

          DO IR=0,NQ
             R=IR*DELT
        
             RX(1)=R
             
             CALL PSCH(IS1,RX,CH1,GRCH)
             CALL PSCH(IS2,RX,CH2,GRCH)

             CH(IR,1)=-CH1
             CH(IR,2)=-CH2

             Z1=Z1-C*CH1*R*R    
             Z2=Z2-C*CH2*R*R

           ENDDO
           
           IF((DABS(Z1-IZ1).GT.CHERR).OR.
     .        (DABS(Z2-IZ2).GT.CHERR)) THEN 
              WRITE(6,*) 'CHOVERLP: THE NUMBER OF INTEGRATION',
     .       ' POINTS MUST BE INCREASED'
              WRITE(6,*) 'CHOVERLP: Z1=',Z1,' IZ1=',IZ1
              WRITE(6,*) 'CHOVERLP: Z2=',Z2,' IZ2=',IZ2
             STOP
           ENDIF

           DO IR=0,NQ
             CH(IR,1)=DBLE(IZ1)*CH(IR,1)/Z1
             CH(IR,2)=DBLE(IZ2)*CH(IR,2)/Z2
           ENDDO 
C
C**********************************************************************





C*****REAL SPACE INTEGRATION OF POISSON EQUATION***********************
C          
          
           CALL NUMEROV(NQ,DELT,CH(0,1),V(0,1))
           CALL NUMEROV(NQ,DELT,CH(0,2),V(0,2))
           
           DO ITB=1,NTBMAX
              R=DLT*(ITB-1)
              NR=NINT(R/DELT)
              NMIN=MAX(0,NR-NPOINT)
              NMAX=MIN(NQ,NR+NPOINT)
              NN=NMAX-NMIN+1
              DO IR=1,NN
                 RAUX(IR)=DELT*(NMIN+IR-1) 
              ENDDO 
              CALL RATINT(RAUX,V(NMIN,1),NN,R,VV1,VD)
              CALL RATINT(RAUX,V(NMIN,2),NN,R,VV2,VD)
 
              VTB(ITB,1)=VV1
              VTB(ITB,2)=VV2
           ENDDO 
         

C
C**********************************************************************


C****FOURIER-TRANSFORM OF RADIAL CHARGE DENSITY************************
C
           CALL RADFFT( 0, NQ, RMAX, CH(0,1), CH(0,1) )
           CALL RADFFT( 0, NQ, RMAX, CH(0,2), CH(0,2) )
C
C**********************************************************************


C*****NEUTRALIZE CHARGE DENSITY FOR FOURIER-SPACE CALCULATION**********
C
           DO IQ=0,NQ
              R=IQ*QMAX/NQ

              CH1 = (CH(IQ,1)-IZ1*CONS)*CH(IQ,2)
              CH2=  (CH(IQ,2)-IZ2*CONS)*CH(IQ,1)
              
              CH(IQ,1) = CH1
              CH(IQ,2) = CH2

           ENDDO
C
C**********************************************************************



C*********THE ELECTROSTATIC ENERGY CORRECTION IS STORED IN 'CORR'******
C  
            DO IR=1,NTBMAX

               R=DLT*(IR-1)
               ENERG1=0.0d0
               ENERG2=0.0d0


               DO IQ=0,NQ
                  Q=IQ*QMAX/NQ
                  Q=Q*R 
                  ENERG1=ENERG1+BESSPH(0,Q)*CH(IQ,1)
                  ENERG2=ENERG2+BESSPH(0,Q)*CH(IQ,2)
               ENDDO 

               ENERG1=ENERG1*QMAX/NQ
               ENERG2=ENERG2*QMAX/NQ
   
               ENERG2=ENERG2*4.0D0*(2.0d0*PI)**2
               ENERG1=ENERG1*4.0D0*(2.0d0*PI)**2
              
               ENERG1=-(ENERG1*R)-(IZ2*(VTB(IR,1)*R-IZ1))
               ENERG2=-(ENERG2*R)-(IZ1*(VTB(IR,2)*R-IZ2))
  
               CORR(IR)=0.5D0*(ENERG1+ENERG2)

            ENDDO 

C
C*******************************************************************


C********CREATING A TABLE WITH SECOND DERIVATIVES FOR SPLINES*******
C
            DEV1=1.0D50
            DEVN=1.0D50
            CALL SPLINE(DLT,CORR,NTBMAX,DEV1,DEVN,CORR2,AUX)

C
C******************************************************************
      
          RETURN
          END




           SUBROUTINE NUMEROV(NR,DELT,Q,V)


C*********************************************************************
C   Being Q(r) a spherical charge density in a homogeneus radial mesh
C   with distance DELT between consecutive points, this routine returns
C   the electrostatic potential generated by this charge distribution.
C   Written by D. Sanchez-Portal, March 1997.
C********************INPUT********************************************
C   INTEGER NR      :    Number of radial points.
C   REAL*8  DELT    :    Distance between consecutive points.
C   REAL*8  Q(0:NR) :    Spherical charge density.
C********************OUTPUT*******************************************
C   REAL*8  V(0:NR) :    Electrostatic potential at mesh points.
C********************BEHAVIOUR****************************************
C   Qtot/r asimptotic behaviour is imposed.
C*********************************************************************

           IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
           DIMENSION Q(0:*),V(0:*) 



            PI=4.0D0*DATAN(1.0D0)
            FOURPI=4.0D0*PI

C********************NUMEROV ALGORITHM******************************* 
C
             V(0)=0.0D0
             V(1)=1.0D0

             DO IR=2,NR
              V(IR)=2.0D0*V(IR-1)-V(IR-2) - FOURPI*DELT**3*
     .      ( Q(IR)*IR+10.0D0*Q(IR-1)*(IR-1)+Q(IR-2)*(IR-2) )/12.0D0
             ENDDO 
C
C********************************************************************

C******************CALCULATE TOTAL CHARGE****************************
C
   
             QTOT=0.0D0
             DO IR=1,NR
               R=IR*DELT
               QTOT=QTOT+R*R*Q(IR)
             ENDDO
             QTOT=4.0D0*PI*QTOT*DELT
C
C********************************************************************

C******* FIXING QTOT/R ASIMPTOTIC BEHAVIOUR**************************
C

             CONS=(QTOT-V(NR))/(NR*DELT)
             
             DO IR=1,NR
                R=IR*DELT
                V(IR)=V(IR)/(IR*DELT)+CONS
             ENDDO 
             V(0)=(4.0D0*V(1)-V(2))/3.0D0
C
C********************************************************************

             RETURN 
             END
              
       



        subroutine rc_vs_e(a,b,r,vps,
     .      ve,nrval,l,el,rnodo)
     
C**************************************************************
C   Calculate the position, rnodo, of the first node of the 
C   radial wavefunction of the pseudopotential Vps, with angular 
C   momentum  l, and energy el.
C   D. Sanchez-Portal, July 1997.
C**************************************************************

        implicit double precision (a-h,o-z) 

     
        include 'atom.h'

        double precision r(nrval),
     .   el,vps(nrval),g(nrmax),drdi(nrmax),h(nrmax),ve(nrval)
        
         integer  l
       
c       if (nrval.gt.nrmax) then   
c        write(6,*) 'Rc_vs_E : Nrmx must be increased to at least',
c    .         nrval
c       endif

        dexpa=dexp(a)
        ab=a*b
        do ir=1,nrval
           
           drdi(ir)=ab
c          r2=ab/a-b

c          if(dabs(r2-r(ir)).gt.1.0d-6) then 
c                write(6,*) ir, dabs(r2-r(ir)),r2,r(ir)
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
                  goto 50
            endif
            
            gold=g(ir)

          enddo 

50        continue 


        
          return 
          end

          
             





        subroutine polarization(a,r,psi,vps,
     .      ve,drdi,nrc,l,el,psipol,nrval)
      
        implicit double precision (a-h,o-z) 
        include 'atom.h'
C*************************************************************
C      This routine calculate the polarization (unoccupied) 
C      orbitals with angular momentum l from the atomic 
C      orbitals with l-1 angular momentum, using a perturbative
C      approach.
C      The routine solve and inhomogeneus version of 
C      Schrodinger eqn.  
C      It is not an optimized algorithm!!!!!!!!!!!!!!!!!
C************************************************************
C************************************************************
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
C************************************************************ 
     
        parameter(nrmin=1,niter=1000,
     .                    cons1=1.0d5,rint=15.0d0)

        double precision r(nrval),psi(nrval),psipol(nrval),
     .   el,vps(nrval),g(nrmax),drdi(nrmax),h(nrmax),ve(nrval)
        
         integer  l

C************************************************************

       
c       if ((nrval.gt.nrmax).or.(nrc.gt.nrmax)) then   
c        write(6,*) 'POLARIZATION: Nrmx must be increased to at least',
c    .        max(nrval,nrc)
c       endif

        rmax=r(nrval)

        if(rmax.gt.rint) then 
          write(6,*) 'POLARIZATION: Rc for the polarization orbitals'
          write(6,*) 'must be smaller than ',rint,' Bohr'
          stop
        endif

        do ir=nrc+1,nrval 
          psi(ir)=0.0d0 
        enddo 
 
          reduc=-0.5d0
          dl=1

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
            g(ir)=cons*(r(ir)**(l+dl+1))/dsqrt(drdi(ir))
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
            gmax=max(gmax,dabs(g(ir)))
            if((g(ir).eq.0.0d0).or.(g(ir)*gold.lt.0.0d0)) then

              nnd=nnd+1
                
              if (nnd.eq.nnodes) then 
                  r0=r(ir-1)
                  g0=gold
                  r1=r(ir)
                  g1=g(ir)
                  rnodo=r0-g0*(r1-r0)/(g1-g0)
                  goto 50
              endif 
           endif
            
           gold=g(ir)


          enddo 

 50        continue




          grmx=g(nrval)/gmax

          if(((dabs(rnodo-rmax).lt.1.0d-3).and.
     .      (dabs(grmx).lt.1.0d-7) )
     .        .or. 
     .        ((rnodo.eq.0.0d0).and.
     .        (dabs(grmx).lt.1.0d-7) ) ) goto 100

*********************We begin by finding a node!!!!**********************

          if((rnd1.eq.0.0d0).and.(rnodo.eq.0.0d0)) then  
             c2=(-reduc)*cons
             if(dabs(c2).le.1.0d0/dabs(cons1)) then 
               index=1
               rnd1=0.0d0
               reduc=(1.0d0+reduc)/2.0d0
             endif  
          elseif((rnd1.eq.0.0d0).and.(rnodo.ne.0.0d0)) then
              rnd1=rnodo
              c1=cons
              c2=2.0d0*cons
          endif  
       
**************************************************************************

      
****************Now we lead this node to Rc*******************************
          if((rnd1.ne.0.0d0).and.(rnodo.eq.0.0d0)) then 
              c2=0.50d0*(c1+c2)
          elseif((rnd1.ne.0.0d0).and.(rnodo.ne.0.0d0)) then 
              if(dabs(rnd1-rnodo).gt.1.0d-6)then 
                 dff1=dabs(rnd1-rmax) 
                 dff2=dabs(rnodo-rmax) 
                 if(dff1.gt.dff2) then 
                   savecons=c2
                   c2=(rmax-rnd1)*(c1-c2)/(rnd1-rnodo)+c1
                   c1=savecons
                   rnd1=rnodo
                 else
                   c2=1.10d0*c2
                 endif 
              else

               if(dabs(cons).gt.1.0d15) then 
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
          write(6,*)'Please try with a Rc no bigger than ',rnd1,' Bohr'
          stop
                          
100       continue
          dnrm=0.0d0
          do ir=1,nrval
             g(ir)=g(ir)*dsqrt(drdi(ir))
             dnrm=dnrm+drdi(ir)*(g(ir)**2)
          enddo 
          dnrm=dsqrt(dnrm)
          do ir=1,nrval
               psipol(ir) = g(ir)/dnrm
          enddo 

 

          end

          
             


         subroutine parabola(a,b,nrc,rphi,rnrm,l,
     .              splnorm,cons1,cons2,nm) 

C*******************************************************
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
C*******************************************************  
C Algorithm based on routine for Golden Section Search
C from Numerical Recipes.
C*******************************************************

          implicit double precision (a-h,o-z) 
          
          parameter (Ratio=0.61803399D0)       
          double precision rphi(nrc),rnrm(nrc) 


C         Hallar el maximo de la funcion de onda
          rfirst=0.05d0
          nfirst=nint(dlog(rfirst/b+1.0d0)/a)+1
          slopold=0.0d0
          do ir=nfirst,nrc
             slop=rphi(ir)-rphi(ir-1) 
             if(slop*slopold.lt.0.0d0) goto 10
             slopold=slop
          enddo 
10        continue
          nrmax=ir-1
          rmin=b*(dexp(a*(nrmax-1))-1.0d0) 
          rmin=1.01d0*rmin
          nmin=nint(dlog(rmin/b+1.0d0)/a)+1   
          nmin=max(nmin,2)
          nmax=nrc-1 
          
          call findp(nrc,nmin,rphi,a,b,l,cmin,gmin) 
          rmin=b*(dexp(a*(nmin-1))-1.0d0) 
          call nrmpal(cmin,gmin,rmin,l,rnrmin)
          rnrmin=1.0d0+rnrmin-rnrm(nmin)
 
          call findp(nrc,nmax,rphi,a,b,l,cmax,gmax)
          rmax=b*(dexp(a*(nmax-1))-1.0d0) 
          call nrmpal(cmax,gmax,rmax,l,rnrmax) 
          rnrmax=1.0d0+rnrmax-rnrm(nmax)         
          
          valmin=(splnorm-rnrmin)**2
          valmax=(splnorm-rnrmax)**2

          nmed=(nmin+nmax)/2 
          do iter=1,nrc
            call findp(nrc,nmed,rphi,a,b,l,cmed,gmed)
            rmed=b*(dexp(a*(nmed-1))-1.0d0) 
            call nrmpal(cmed,gmed,rmed,l,rnrmed) 
            rnrmed=1.0d0+rnrmed-rnrm(nmed)
              
            valmed=(splnorm-rnrmed)**2

            if((valmed.lt.valmin).and.(valmed.lt.valmax)) goto 20
            nmed=nmed+1
            if(nmed.eq.nmax) goto 15
          enddo 
15        continue
          nmed=(nmin+nmax)/2
          do iter=1,nrc
             call findp(nrc,nmed,rphi,a,b,l,cmed,gmed)
             rmed=b*(dexp(a*(nmed-1))-1.0d0)
             call nrmpal(cmed,gmed,rmed,l,rnrmed)
             rnrmed=1.0d0+rnrmed-rnrm(nmed)

             valmed=(splnorm-rnrmed)**2


             if((valmed.lt.valmin).and.(valmed.lt.valmax)) goto 20
             nmed=nmed-1
             if(nmed.eq.nmin) goto  20
          enddo 
20        continue

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
            r=b*(dexp(a*(n1-1))-1.0d0)
            call nrmpal(c1,g1,r,l,rn1)
            rn1=1.0d0+rn1-rnrm(n1)
            val1=(splnorm-rn1)**2

            call findp(nrc,n2,rphi,a,b,l,c2,g2)
            r=b*(dexp(a*(n2-1))-1.0d0)
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
               r=b*(dexp(a*(n2-1))-1.0d0)
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
               r=b*(dexp(a*(n1-1))-1.0d0)
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

             return


            end




          subroutine findp(nrc,nm,rphi,a,b,l,cons1,cons2)

C******************************************************** 
C  This routine provides the constants Cons1 and 
C  Cons2 and described in subroutine 'parabola' 
C
C Written by D. Sanchez-Portal, July 1997.
C********************************************************

          implicit double precision (a-h,o-z) 
          double precision rphi(nrc)
          
          rm=b*(dexp(a*(nm-1)) + 1.0d0) 
          rm1=b*(dexp(a*(nm-2)) + 1.0d0)
          rm2=b*(dexp(a*nm) + 1.0d0)
          drdi=a*b*dexp(a*(nm-1))

          frsp=rphi(nm)/rm
          dfrsp=0.5d0*(rphi(nm+1)/rm2
     .       -rphi(nm-1)/rm1)
          dfrsp=dfrsp/drdi

          cons1= 0.5d0*(dfrsp*rm-l*frsp)/(rm**(l+2))
          cons2= frsp/(rm**l)-cons1*(rm**2)

          
 
          end
        


 
          subroutine nrmpal(c1,c2,r,l,dnrm)

C returns the norm of a parabolic function
C    f(r')= r'^l (c1*r'^2 + c2)  r'< r
C           0 otherwise
C
C D. Sanchez-Portal, July 1997.

          implicit double precision (a-h,o-z)
          
          dnrm=(c1**2)*r**(2*l+7)/(2*l+7) 
     .     + (2.0d0*c1*c2)*r**(2*l+5)/(2*l+5) 
     .     + (c2**2)*r**(2*l+3)/(2*l+3)


          end 




        subroutine radii_ps(vps,rofi,Zval,nrval,lmxkb,
     .          nrgauss, rgauss, rgauss2)

C*******************************************************************
C   This rouitne returns the maximum radius for the
C   Kleinman-Bylander projectors with a standart choice
C   of the local potential.
C   Check also at which radius the asymptotic 2*Zval/r
C   behaviour is achieved.  
C   D. Sanchez-Portal, Aug. 1998
C*******************************************************************


        implicit none
        include 'atom.h'

        double precision 
     .     vps(nrmax,0:lmaxd), rofi(nrmax), 
     .     rgauss, rgauss2, Zval
        integer nrval, lmxkb, nrgauss


C***************internal variables*****************************
        double precision eps, dincv, r
        integer ir, l
        parameter (eps=1.0d-4)

           

C*******Iterate over the possible local potentials*****************

          rgauss=0.0d0
          nrgauss=0
          do l=0,lmxkb-1
               
             do ir=nrval,2,-1
                 dincv=dabs(vps(ir,l)-vps(ir,lmxkb))
                 if(dincv.gt.eps) goto 10
             enddo 
10           rgauss=max(rofi(ir),rgauss)
             nrgauss=max(ir,nrgauss)
          enddo

          do ir=nrval,2,-1
              r=rofi(ir)
              dincv=dabs(vps(ir,0)*r+2.0d0*zval)
              if(dincv.gt.eps) goto 20
          enddo 
            
20        rgauss2=rofi(ir)
          if(lmxkb.eq.0) then 
             rgauss=rgauss2 
             nrgauss=ir
          endif 



             return

             end
            
           



       subroutine vlocal2(Zval, nrval, a, rofi, drdi, s, vps,
     .              nrgauss,vlocal,nchloc,chcore) 

C**************************************************************************** 
C     This routine generates the local pseudopotential appropiate 
C     for species with  a large core.
C     
C     Written by D. Sanchez-Portal, Aug. 1998
C**************************************************************************** 
   

       implicit none 
       include 'atom.h'
        
       double precision rofi(nrmax), drdi(nrmax), s(nrmax),
     .     vps(nrmax), vlocal(nrmax), a, Zval, chcore(nrmax)
       integer nrgauss, nrval, nchloc


C********************Internal variables***************************************
       double precision 
     .      vlc, r, dev, dev2, dev3, var1, var2, var3, v1, v2, v3, v4,
     .      dm11, dm12, dm13, dm21, dm22, dm23, dm31, dm32, dm33, eps,
     .      g0, g1, g2, g3, g4, d2g, d2u, cons, a2b4, qtot, pi   
       integer 
     .      ndevfit, ir  

       parameter (eps=1.0d-5)
                     
 
C***********************PI*****************************************************
                     pi=dacos(-1.0d0)        





C**********Continuity up to second derivative*********************************
                ndevfit=2
C**********Continuity up to third derivative**********************************
C               ndevfit=3



                nrgauss=nrgauss+3

                do ir=1,nrval
                   vlocal(ir)=vps(ir)*rofi(ir)
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

CLocal potential is Vloc(r)=v3*exp(v1*r^2+v2*r^3) inside Rgauss and equals the 
Call-electron atomic potential outside Rgauss
CWe impose the continuity up to second        
    
             if(ndevfit.eq.2) then               
               vlc=vlocal(nrgauss)
               r=rofi(nrgauss)

               var1=dev/vlc-1.0d0/r
               var2=dev2/vlc-2.0d0*var1/r -(var1**2)

               dm11=2.0d0*r
               dm12=3.0d0*r*r
               dm21=2.0d0
               dm22=6.0d0*r

               v1=(dm22*var1-dm12*var2)/(6.0d0*r*r)
               v2=(dm11*var2-dm21*var1)/(6.0d0*r*r)
               v3=vlc/(r*dexp((v1+v2*r)*r*r))


             elseif(ndevfit.eq.3) then 

C******************************************************************************
C We can also construct a local potential Vloc(r)=v4*exp(v1*r^2+v2*r^3+v3*r^4),
C this new coefficient allows us to impose the continuity of the potential up
C to the third derivative.
C******************************************************************************
 
            vlc=vlocal(nrgauss)
            r=rofi(nrgauss)
            
            var1=dev/vlc-1.0d0/r
            var2=dev2/vlc-2.0d0*var1/r-(var1**2)
            var3=dev3/vlc-3.0d0*var1*var2-(var1**3)
     .                           -3.0d0*(var1**2+var2)/r

            dm11=2.0d0*r
            dm12=3.0d0*r*r
            dm13=4.0d0*r*r*r
            dm21=2.0d0
            dm22=6.0d0*r
            dm23=12.0d0*r*r
            dm31=0.0d0
            dm32=6.0d0
            dm33=24.0d0*r

            v1=((var1*dm22*dm33+var2*dm13*dm32+var3*dm12*dm23)
     .   -(var3*dm22*dm13+var1*dm32*dm23+var2*dm12*dm33))/(48.0d0*r*r*r) 
            v2=((var2*dm11*dm33+var3*dm21*dm13+var1*dm23*dm31)
     .   -(var2*dm31*dm13+var3*dm23*dm11+var1*dm21*dm33))/(48.0d0*r*r*r)
            v3=((var3*dm11*dm22+var2*dm12*dm31+var1*dm32*dm21)
     .   -(var1*dm22*dm31+var3*dm21*dm12+var2*dm11*dm32))/(48.0d0*r*r*r)
            v4=vlc/(r*dexp((v1+v2*r+v3*r*r)*r*r))
            
          endif 


 
             do ir=1,nrval
               r=rofi(ir)
               if(ir.le.nrgauss) then 
 
C************If second derivative fit**************************************
                if(ndevfit.eq.2) then 
                    vlocal(ir)=v3*dexp((v1+v2*r)*r*r)
C**************************************************************************

C************If third derivative fit***************************************
                elseif(ndevfit.eq.3) then 

                    vlocal(ir)=v4*dexp((v1+v2*r+v3*r*r)*r*r)
c**************************************************************************
                endif 

               else
                    vlocal(ir)=vps(ir)
               endif 

             enddo 



C Once we have the local potential we define the 'local-pseudopotential 
C charge' which help us to calculate the electrostatic interation 
C between the ions


          a2b4=0.25d0*a*a 
          qtot=0.d0 
          do ir=1,nrval-1
             
            g2=vlocal(ir)*rofi(ir)
            if(dabs(g2+2.0d0*zval).lt.eps) goto 10

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
             if(ndevfit.eq.2)  then
              r=rofi(ir)

              g0=v3*dexp((v1+v2*r)*r**2)
              g1=(2.0d0*v1+3.0d0*v2*r)
              g2=2.0d0*v1+6.0d0*v2*r
              g3=(g2+g1*g1*r*r+2.0d0*g1)*g0
             
              cons=8.0d0*pi
              chcore(ir)= (-g3)/cons
              qtot= qtot  + 0.5d0*g3*r*r*drdi(ir)
 
C**************If third derivative fit****************************************
     
              elseif(ndevfit.eq.3)  then

              r=rofi(ir)
           
              g0=v4*dexp((v1+v2*r+v3*r*r)*r*r)     
              g1=(2.0d0*v1+3.0d0*v2*r+4.0d0*v3*r*r)
              g2=(2.0d0*v1+6.0d0*v2*r+12.0d0*v3*r*r)    
              g3=(g2+g1*g1*r*r+2.0d0*g1)*g0   

              cons=8.0d0*pi
              chcore(ir)= -g3/cons
              qtot= qtot  + 0.5d0*g3*r*r*drdi(ir)
             endif 
C****************************************************************************



             endif



          enddo              


10        continue
          nchloc=ir          

          do ir=1,nchloc
            r=rofi(ir)
            chcore(ir)=zval*chcore(ir)/qtot
          enddo  





           end



       subroutine vlocal1(Zval, nrval, a, rofi, drdi, s, rgauss,
     .          vlocal, nchloc, chcore)


C****************************************************************************
C     This routine generates a smooth local pseudopotential.
C     Written by D. Sanchez-Portal, Aug. 1998
C****************************************************************************

         implicit none
         include 'atom.h'
         double precision 
     .     Zval, rofi(nrmax), vlocal(nrmax), rgauss, chcore(nrmax),
     .     rhor, drdi(nrmax), s(nrmax), a
         integer
     .     nrval, nchloc



C*********************Internal variables******************************* 

          double precision van, factor, alp, cutoff1, cutoff2,
     .        gexp, qtot, eps, pi, chc, r, Rchloc, rhor1
          integer ir 
          character loctype*3

          parameter(eps=1.0d-4)


C**Usual local potential (generated with an optimum Vandebilt function)**
                   loctype='new' 
C************************************************************************

C***The very first local potential used by SIESTA was a the electrostatic
C***potential generated by a gaussian distribution ===> loctype='old' 
C                  loctype='old'
C*************************************************************************


C*******************************PI****************************************
                        pi=dacos(-1.0d0)


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




        if(loctype.eq.'new') then          


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
                factor=1.815d0
C***********************************************************************
              
C*********** Scaling factor for local-pseudopot. charge*****************
                  alp=factor/rgauss
C***********************************************************************


       write(6,'(/,a,f10.3,a)')
     .             'VLOCAL1: 99.0% of the norm of Vloc inside ',
     .                 (alp*cutoff1)**2,' Ry'
       write(6,'(a,f10.3,a)')
     .             'VLOCAL1: 99.9% of the norm of Vloc inside ',
     .                 (alp*cutoff2)**2,' Ry'


         elseif(loctype.eq.'old') then 

C This is just a gaussian !!!!!!!!!!!!!!!!!                 
                 van=0.00001d0 
                 rgauss=0.80d0
                 factor=2.0d0


C*********** Scaling factor for local-pseudopot. charge*****************
                alp=factor/rgauss  
C***********************************************************************
                
         endif 

          qtot=0.0d0 
          do ir=1,nrval
             r=rofi(ir) 
             gexp=dsinh(van*alp*r)/dsinh(van)
             gexp=gexp*gexp
             rhor=dexp(-gexp)  
             if(ir.eq.1) rhor1=-rhor
             chcore(ir)=(-4.0d0)*pi*rhor*r*r
             qtot=qtot+rhor*drdi(ir)*r*r
          enddo

          qtot=4.0d0*pi*qtot 
          nchloc=0  
          do ir=nrval,1,-1
             chc=zval*chcore(ir)/qtot
             chcore(ir)=chc
             if((dabs(chc).gt.eps).and.(nchloc.eq.0)) then    
                 nchloc=ir+1
             endif
          enddo 
          Rchloc=rofi(nchloc)


         call vhrtre(chcore,vlocal,rofi,drdi,s,nrval,a)
   
          do ir=2,nrval 
             r=rofi(ir)  
             chcore(ir)=chcore(ir)/(4.0d0*pi*r*r)
             if (r.gt.1.1d0*Rchloc) then
                 vlocal(ir)=(-2.0d0)*zval/rofi(ir)
             endif

          enddo 
          chcore(1)=rhor1*zval/qtot


          end





        subroutine schro_eq(Zval,rofi,vps,ve,s,drdi,
     .        nrc,l,a,b,nnodes,nprin,
     .        e,g) 
        
          implicit none   
          include 'atom.h'   
          double precision
     .        Zval, rofi(*),vps(*),ve(*),s(*),drdi(*),a,b,e,g(*)
          integer  
     .        nrc,l,nnodes,nprin
     
C********************* Internal variables*********************************  

           double precision
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
              dsq=dsqrt(drdi(ir))
              phi=phi*dsq
              g(ir)=phi
            enddo 
            g(1)=0.0d0 

            dnrm=0.0d0
            do ir=2,nrc
               phi=g(ir)
               dnrm=dnrm+phi*phi*drdi(ir) 
            enddo
            dnrm=dsqrt(dnrm)

            do ir=2,nrc
               g(ir)=g(ir)/dnrm
            enddo 
            
            end
             







         subroutine ghost(Zval,rofi,vps,vlocal,
     .        ve,s,drdi,nrval,l,a,b,
     .        eigenl,rphi,ighost)

C**********************************************************************
C This routine checks the possible existence of ghost states. 
C   Input:
C        vps(*)    : pseudopotential for angular momentum l
C        vlocal(*) : local pseudopotential
C        rphi (*)  : first radial pseudowavefunctions for Vps.
C        eigenl    : eigenvalue 
C     
C  Output:
C        ighost:   if ighost=0 no ghost states
C                  if ighost=1, ghost states exist
C
C   Written by D. Sanchez-Portal, Aug. 1998
C**********************************************************************

         implicit none
          include 'atom.h'
          double precision
     .        Zval, rofi(*),vps(*),ve(*),s(*),drdi(*),a,b,eigenl,
     .        rphi(*), vlocal(*)
          integer
     .        nrval,l,ighost

C********************* Internal variables*********************************

           double precision
     .        dnrm, avgv, phi, e,
     .        elocal(2), g(nrmax), vl, vphi, dkbcos
C    .        , dknrm
           integer
     .        ir, nprin, nnode
           logical  called
C**********************SAVE FOR NEXT CALL******************************
C
              save called
C
C**********************************************************************

C**********************************************************************
C
       data  called /.false./
C
C**********************************************************************

              if(.not.called) ighost=0 



C******************GHOST ANALYSIS**********************************


C***CALCULATE EIGENVALUES OF LOCAL POTENTIAL FOR GHOST ANALYSIS******
              

C* ATTENTION , 'Ve' is the screenig potential generated from valence*
C* pseudo-charge given by the pseudopotential generation program ****
C******************************************************************** 
          do 20 nnode=1,2
              nprin=l+1

              call schro_eq(Zval,rofi,vlocal,ve,s,drdi,
     .        nrval,l,a,b,nnode,nprin,
     .        e,g)

              elocal(nnode)=e
              
  20       continue

c        write(6,*) 'GHOST: Ground state vlocal for L=',l,elocal(1)
c        write(6,*) 'GHOST: First excited state for L=',l,elocal(2)

C****************CALCULATE KB-COSINE**********************************


             dnrm=0.0d0
             avgv=0.0d0
             do 30 ir=2,nrval
               vl=(vps(ir)-vlocal(ir))
               phi=rphi(ir)
               vphi=vl*phi
               dnrm=dnrm+vphi*vphi*drdi(ir)
               avgv=avgv+vphi*phi*drdi(ir)
  30         continue


             dkbcos=dnrm/(avgv+1.0d-20)
c            dknrm=1.0d0/(dsqrt(dnrm)+1.0d-20)

C***************GHOST ANALYSIS*****************************************


          if(dkbcos.gt.0.0d0) then

              if(eigenl.gt.elocal(2)) then
                 write(6,"(a,i3)")
     .            'GHOST: WARNING: Ghost state for L =', l
                 ighost=1
              else
                 write(6,'(a,i3)') 'GHOST: No ghost state for L =',l
              endif

           elseif(dkbcos.lt.0d0) then

              if(eigenl.gt.elocal(1)) then
                 write(6,"(a,i3)")
     .            'GHOST: WARNING: Ghost state for L =', l
                 ighost=1
              else
                 write(6,'(a,i3)') 'GHOST: No ghost state for L =',l
              endif

           elseif(dkbcos.eq.0.0d0) then

               write(6,"('GHOST: vps = vlocal, no ghost for L =',i3)") l

           endif



 
           end




        subroutine KBproj(rofi,drdi,vps,vlocal,nrval,l,
     .    dkbcos,rphi,nrc)  

C****************************************************************
C    This routine calculates the Kleinman-Bylander projector
C    with angular momentum l.
C
C  Written by D. Sanchez-Portal, Aug. 1998
C****************************************************************

         implicit none
          include 'atom.h'
          double precision
     .        rofi(*),vps(*),drdi(*),
     .        rphi(*), vlocal(*), dkbcos
          integer
     .        nrval, l, nrc


C********************* Internal variables*********************************

          double precision 
     .        eps, dnrm, vl, vphi, avgv, r, phi, dknrm,
     .        dincv, rc
          integer
     .        ir


              parameter (eps=1.0d-6)

             dnrm=0.0d0
             avgv=0.0d0
             do 10 ir=2,nrval
               r=rofi(ir)
               vl=(vps(ir)-vlocal(ir))
               phi=rphi(ir)
               vphi=vl*phi
               dnrm=dnrm+vphi*vphi*drdi(ir)
               avgv=avgv+vphi*phi*drdi(ir)
  10         continue


             dkbcos=dnrm/(avgv+1.0d-20)
             dknrm=1.0d0/(dsqrt(dnrm)+1.0d-20)

             
            

C**********DEFINE THE CUT-OFF RADII FOR THE KB PROJECTORS************
C Warning these radii should be quite short, if it is not the case
C something is probably wrong in this part of the program.
C It will display a warning if Rc>4.5 a.u.or Rc < 0.5a.u.!!!!!!!!!!!!

              do 20 ir=nrval,2,-1 
                 phi=(rphi(ir)/rofi(ir))*dknrm
                 dincv=dabs(vps(ir)-vlocal(ir))*phi
                 if(dincv.gt.eps) goto 21
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
                phi=rphi(ir)/r
                vphi=vl*phi*dknrm
                rphi(ir)=vphi/r**l
30           continue
             rphi(1)= ( rphi(2)*rofi(3)**2 - rphi(3)*rofi(2)**2 ) /
     .             (      rofi(3)**2 -      rofi(2)**2 )




              end




               subroutine xc_check(xcfunc,xcauth,icorr)

C**************************************************************
C   Checking the functional used for exchange-correlation energy.
C
C Written by D. Sanchez-Portal, Aug. 1998
C**************************************************************


                character
     .             xcfunc*3, xcauth*4,icorr*2





       write(6,'(/a)') 'xc_check: Exchange-correlation functional:'
       if(((xcauth.eq.'CA').or.(xcauth.eq.'PZ')).and.
     .    ((xcfunc.eq.'LDA').or.(xcfunc.eq.'LSD'))) then

          write(6,'(a)') 'xc_check: Ceperley-Alder'
          if(icorr.ne.'ca') then
           write(6,'(a)')
     .      'xc_check: WARNING: Pseudopotential generated with'
           if(icorr.eq.'pw') write(6,'(a)')
     .      'xc_check: WARNING: Perdew-Wang 1992 functional'
           if(icorr.eq.'pb') write(6,'(a)')
     .'xc_check: WARNING: GGA Perdew, Burke & Ernzerhof 1996 functional'
          endif

       elseif((xcauth.eq.'PW92').and.
     .    ((xcfunc.eq.'LDA').or.(xcfunc.eq.'LSD'))) then

         write(6,'(a)') 'xc_check: Perdew-Wang 1992'
         if(icorr.ne.'pw') then
           write(6,'(a)')
     .       'xc_check: WARNING: Pseudopotential generated with'
           if(icorr.eq.'ca')
     .   write(6,'(a)') 'xc_check: WARNING: Ceperly-Alder functional'
           if(icorr.eq.'pb') write(6,'(a)')
     .'xc_check:WARNING: GGA Perdew, Burke & Ernzerhof 1996 functional'
         endif
       elseif((xcauth.eq.'PBE').and.(xcfunc.eq.'GGA')) then

         write(6,'(a)')
     .     'xc_check: GGA Perdew, Burke & Ernzerhof 1996'
         if(icorr.ne.'pb') then
           write(6,'(a)')
     .       'xc_check: WARNING: Pseudopotential generated with'
           if(icorr.eq.'ca')
     .    write(6,'(a)') 'xc_check: WARNING: Ceperly-Alder functional'
           if(icorr.eq.'pw')
     . write(6,'(a)') 'xc_check: WARNING: Perdew-Wang 1992 functional'
         endif

       else

          write(6,'(a)')
     . 'xc_check: ERROR: Exchange-correlation functional not allowed'
          write(6,'(a)') 'xc_check: ERROR: xc.functional= ',xcfunc
          write(6,'(a)') 'xc_check: ERROR: xc.authors= ',xcauth
          stop

       endif

        return
       end



           subroutine comcore(is,a,b,rofi,chcore,
     .                               nrval,nicore,flting)

C****************************************************************
C Generates the common blocks with the pseudo-core information.
C  D. Sanchez-Portal, Aug. 1998
C****************************************************************


           implicit none
           include 'atom.h'

           integer nrval, is
 
           double precision 
     .        chcore(nrmax), rofi(nrmax), flting, a ,b

           character nicore*4
          
C**********Internal variables*********************************************
C
           integer nrcore,ir, nr, nmin, nmax, nn, itb

           double precision r2, chc, r, pi, delt, Rcore, dy,
     .     yp1, ypn
  
           double precision 
     .          aux(ntbmax)

           double precision eps
            parameter (eps=1.0d-6)  
          
           double precision deltmax
           common/cmdelt/deltmax
C
C***********************************************************************

C****NUMBER OF POINTS USED BY RATINT FOR THE INTERPOLATION*************
C
          integer npoint
          parameter(npoint=4)
C
C***********************************************************************

C**********Common block variables****************************************
C  
          double precision 
     .        coretab(ntbmax+1,2,nsmax)

          common/cmcore/coretab
C
C***********************************************************************



C*******************PI****************************************************
C
                 pi=dacos(-1.0d0) 
C
C***************************************************************************

         coretab(1,2,is)=0
         if(nicore.ne.'nc  ') then



           if(flting.gt.0.0d0) then

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
     .       'comcore: Pseudo-core radius Rcore=',Rcore

C****************TABLE WITH THE PSEUDO_CORE DATA********************

            delt=Rcore/(dble(ntbmax-1)+1.0d-20)  

            if(delt.gt.deltmax) then 
              write(6,'(a)') 
     .    'comcore: WARNING It might be a good idea to increase' 
              write(6,'(a)')
     .    'comcore: WARNING parameter ntbmax (in file atom.h) '
              write(6,'(a,i6)')
     .    'comcore: WARNING to at least ntbmax = ', 
     .            nint(Rcore/deltmax)+2 
            endif 
   
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

            yp1=1.0d50
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



           return 
 
           end                  



           subroutine comlocal(is,a,b,rofi,
     .            chlocal,nchloc,flting)

C****************************************************************
C Generates the common blocks with the local-pseudopotential charge
C  D. Sanchez-Portal, Aug. 1998
C****************************************************************

           implicit none
           include 'atom.h'

           integer nchloc, is
 
           double precision 
     .        chlocal(nrmax), rofi(nrmax), a ,b, flting

          
C**********Internal variables*********************************************
C
           integer nr, nmin, nmax, nn, itb, indx, is2

           double precision chc, r, delt, Rchloc, dy,
     .     yp1, ypn, rmax
  
           double precision 
     .          aux(ntbmax) 
           
           double precision 
     .          deltmax
           common/cmdelt/deltmax

C
C***********************************************************************

C****NUMBER OF POINTS USED BY RATINT FOR THE INTERPOLATION*************
C
          integer npoint
          parameter(npoint=4)
C
C***********************************************************************

C**********Common block variables****************************************
C  
          integer ns2
           parameter(ns2=((nsmax+1)*nsmax)/2)

          double precision 
     .        chloctab(ntbmax+1,2,nsmax),
     .        corrtab(ntbmax+1,2,ns2)

          common/cmchloc/chloctab 
          common/cmcorr/corrtab
  
C
C***********************************************************************



        if(flting.gt.0.0d0) then 
  
          Rchloc=rofi(nchloc)


          delt=Rchloc/(dble(ntbmax-1)+1.0d-20) 
          
          if(delt.gt.deltmax) then
              write(6,'(a)')
     .    'comlocal: WARNING It might be a good idea to increase'
              write(6,'(a)')
     .    'comlocal: WARNING parameter ntbmax (in file atom.h) '
              write(6,'(a,i6)')
     .    'comlocal: WARNING to at least ntbmax = ', 
     .        nint(Rchloc/deltmax)+2 
           endif 

          chloctab(1,1,is)=delt
          chloctab(1,2,is)=1.0d0
          do itb=1,ntbmax
             r=delt*(itb-1)
             nr=nint(dlog(r/b+1.0d0)/a)+1
             nmin=max(1,nr-npoint)
             nmax=min(nchloc,nr+npoint)
             nn=nmax-nmin+1
             call ratint(rofi(nmin),chlocal(nmin),nn,r,chc,dy)

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

        elseif( flting.lt.0.0d0) then 
 

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
           return 
 
           end                  




       subroutine new_specie(iz,lmxkb, lmxo,
     .  nzeta, rco, lambda, atm_label,
     .  npolorb, semic, lsemic, 
     .  isnew, new, no, nkb, q)

C*************************************************************
C  Checks if 'atom' has been called for a new species or if 
C  the requested information was previously calculted and, 
C  therefore, it is yet available.
C  Output: 
C         logical new
C         integer isnew: species index for the new species
C
C  Written D. Sanchez-Portal, Aug. 1998
C*************************************************************

       implicit none
       include 'atom.h' 
       integer  maxos
       parameter (maxos=2*nzetmx*lmx2)


       double precision 
     .   rco(nzetmx,0:lmaxd), lambda(nzetmx,0:lmaxd),
     .   q(maxos)

       integer
     .  iz, lmxkb, lmxo, nzeta(0:lmaxd), npolorb(lmaxd),
     .  lsemic, isnew, no, nkb 
      
       character 
     .   atm_label*20 

       logical
     .     new, semic
     
     

C      Internal and common variables


       double precision 
     .  rcosave(nzetmx,0:lmaxd,nsmax),
     .  lambdasave(nzetmx,0:lmaxd,nsmax),
     .  rcotb(nzetmx,0:lmaxd,nsmax), 
     .  rcpoltb(nzetmx,lmaxd,nsmax),
     .  lambdatb(nzetmx,0:lmaxd,nsmax),
     .  qtb(maxos,nsmax)

       double precision
     .  table((ntbmax+2),-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tab2(ntbmax,-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tabpol((ntbmax+2),nzetmx*lmaxd,nsmax),
     .  tab2pol(ntbmax,nzetmx*lmaxd,nsmax)

       integer
     .  ns, nsold, izeta, l, ismax,ix, isold,
     .  lmax, nzetamax
 

       integer
     .  izsave(nsmax),lmxkbsave(nsmax), 
     .  lmxosave(nsmax),  
     .  npolorbsave(lmaxd,nsmax), 
     .  lsemicsave(nsmax), nzetasave(0:lmaxd,nsmax),
     .  nomax(nsmax), nkbmax(nsmax)
       
       character
     .   label_save(nsmax)*20

       logical 
     .    semicsave(nsmax), overflow
     

C**************common declarations*********************************

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
             common/cmq/qtb
             common/cmtab/table,tabpol
             common/cmspline/tab2,tab2pol


C***********Save for the next call********************************

         save rcosave, lambdasave, isold
         
C*****************************************************************
              data isold / 0 /
C*****************************************************************   


          if(iz.ne.0) then 



C***** IS THIS A NEW SPECIES ? ****************************************
C***** We first compare the inputs, if all the inputs are equal then***
C***** this is not a new species **************************************
          
          if(isold.eq.0) goto 9

 
          do 1 ns=1,isold
             if((izsave(ns).eq.iz).and.
     .               (label_save(ns).eq.atm_label)) then
               nsold=ns
               goto 5
             endif

  1       continue


          goto 9

  
  5       if(  (lmxkb.ne.lmxkbsave(nsold)).and.
     .              (iz.gt.0)) goto 8
          if(lmxo.ne.lmxosave(nsold)) goto 8

          do 6 l=0,lmxo
            if(nzeta(l).ne.nzetasave(l,nsold)) goto 8

  
            do izeta=1,nzeta(l)
              if(rco(izeta,l).ne.rcosave(izeta,l,nsold)) goto 8 
              rco(izeta,l)=rcotb(izeta,l,nsold)
              if(lambda(izeta,l).ne.lambdasave(izeta,l,nsold)) goto 8 
              lambda(izeta,l)=lambdatb(izeta,l,nsold)
            enddo  
  6       continue
          do l=1,min(lmxo+1,lmaxd)
              if(npolorb(l).ne.npolorbsave(l,nsold)) goto 8  
          enddo  


          if(semic.neqv.semicsave(nsold)) goto 8
          if (lsemic.ne.lsemicsave(nsold)) goto 8
           
          no=nomax(nsold)
          nkb=nkbmax(nsold)  
          if(iz.lt.0) then 
             if(nkb.ne.0) then
               write(6,'(2a)')
     .          'new_specie: ERROR: Specie: ', label_save(nsold) 
               write(6,'(a)')
     .          'new_specie: ERROR: For negative atomic number'  
               write(6,'(a)')
     .          'new_specie: ERROR: we should have no KB projectors' 
                 stop
             endif 
          endif 

          isnew=nsold
          new=.false.  
          do izeta=1,no
              q(izeta)=qtb(izeta,nsold)
          enddo 
       

          write(6,'(/,a,i2,a,i2)')
     .     'new_specie: WARNING: Data for species', isold+1,
     .     ' are identical to those of the previous species',nsold

          return 

  8       write(6,'(/,a)')
     .'new_specie:WARNING: There are previous data for the same species'
          write(6,'(a)')
     .'new_specie: WARNING: Some of the arguments have been changed'

  9       continue


C*******ADDING A NEW SPECIES TO THE LIST********************************
           
          overflow=.false.
          do ns=1,isold
              if(atm_label.eq.label_save(ns) ) then
                  write(6,'(/,2a)')
     .        'new_specie: WARNING: Two species with the same label  ',
     .                atm_label
              endif
          enddo

          isnew=isold+1
          isold=isold+1
          ismax=isnew
          new=.true.  
          if(iz.lt.0) lmxkb=0
          lmax=max(lmxo,lmxkb) 
          nzetamax=0
          do l=0,lmaxd
              nzetamax=max(nzeta(l),nzetamax) 
              if(l.gt.0) nzetamax=max(npolorb(l),nzetamax) 
          enddo

          if(isnew.gt.nsmax) then
             write(6,"(2a,i4)")
     .        'new_specie: ERROR: Parameter nsmax must be increased ',
     .        'to at least ', isnew 
               overflow=.true.
          endif 

          if(lmax.gt.lmaxd) then 
              write(6,"(2a,i4)") 
     .        'new_specie: ERROR: Parameter lmaxd must be increased ',
     .        'to at least ', lmax 
               overflow=.true.
          endif                
           
          if(nzetamax.gt.nzetmx) then 
             write(6,"(2a,i4)") 
     .        'new_specie: ERROR: Parameter nzetmx must be increased ',
     .        'to at least ', nzetamax
               overflow=.true.
          endif                

          if(.not.overflow) then  

            izsave(isnew)=iz
            lmxosave(isnew)=lmxo
            lmxkbsave(isnew)=lmxkb
            label_save(isnew)=atm_label
            nkb=(lmxkb+1)**2 
            if(iz.lt.0) nkb=0
            nkbmax(isnew)=nkb 


          no=0
          do l=0,lmxo
             
             nzetasave(l,isnew)=nzeta(l)

             do izeta=1,nzeta(l)
                rcosave(izeta,l,isnew)=rco(izeta,l)
                lambdasave(izeta,l,isnew)=lambda(izeta,l)
             enddo
             
             no=no+(2*l+1)*nzeta(l)

           enddo   

           do l=1,min(lmxo+1,lmaxd)
            npolorbsave(l,isnew)=npolorb(l)
            no=no+(2*l+1)*npolorb(l)
           enddo 


           nomax(isnew)=no
          
           
          do izeta=1,maxos
             q(izeta)=0.0d0
          enddo


          elseif(overflow) then 
             write(6,'(a)') 
     .        'new_specie: ERROR: Check dimensions in file atom.h'
               stop 'new_specie: ERROR: Check dimensions in file atom.h'  
          endif 

 





          else
C****If iz=0 all the tables are set to zero, everything is reinitialized** 


           new=.false. 
           isnew=0
           ismax=0
           isold=0
           do ns=1,nsmax 
              izsave(ns)=0
              lmxosave(ns)=0
              lmxkbsave(ns)=0
              label_save(ns)='  '
              nkbmax(ns)=0
              nomax(ns)=0  
              semicsave(ns)=.false.
              lsemicsave(ns)=0 
              do l=0,lmaxd 
                  nzetasave(l,ns)=0
                  do izeta=1,nzetmx 
                     rcotb(izeta,l,ns)=0.0d0 
                     lambdatb(izeta,l,ns)=0.0d0  
          if(l.gt.0) rcpoltb(izeta,l,ns)=0.0d0
                     rcosave(izeta,l,ns)=0.0d0 
                     lambdasave(izeta,l,ns)=0.0d0 
                  enddo    
              enddo  
              do l=-(lmaxd+1),nzetmx*(lmaxd+1)
                 do ix=1,ntbmax 
                    table(ix,l,ns)=0.0d0
                    tab2(ix,l,ns)=0.0d0
                 enddo
                 table(ntbmax+1,l,ns)=0.0d0
                 table(ntbmax+2,l,ns)=0.0d0 
             enddo  

              do l=1,nzetmx*lmaxd
                 do ix=1,ntbmax
                    tabpol(ix,l,ns)=0.0d0
                    tab2pol(ix,l,ns)=0.0d0
                 enddo
                 tabpol(ntbmax+1,l,ns)=0.0d0
                 tabpol(ntbmax+2,l,ns)=0.0d0
             enddo  


             do izeta=1, maxos
                qtb(izeta,ns)=0.0d0
             enddo 
           enddo 
           

          endif 



        end








            subroutine read_vps(atm_label, lmxo, lmxkb,
     .             nrval,a,b,rofi,drdi,s,vps,
     .             rho, chcore, zval, nicore, irel, icorr)

C*******************************************************************
C   Read the file generated by the pseudopotential generation code.
C   Written by D. Sanchez-Portal, Aug. 1998
C*******************************************************************

           implicit none 
      
           include 'atom.h'  
           double precision
     .        rofi(nrmax), drdi(nrmax), s(nrmax), vps(nrmax,0:lmaxd),
     .        rho(nrmax), chcore(nrmax) 

           double precision
     .         a, b, zval
           
           integer  nrval, lmxo, lmxkb 

           character atm_label*20, nicore*4, irel*3, icorr*2


C***********Internal variables **************
           
           double precision 
     .     ve(nrmax)

           double precision 
     .     r2, ea, rpb 

           integer  
     .        nr, nodd, lmax, linput, npotd, npotu,
     .        ndown, nup, l, ir, i, itext 

           character 
     .         fname*50, namatm*2, 
     .         method(6)*10,text*70,paste*50

           logical 
     .         found          
       
           fname = paste(atm_label,'.vps')
           inquire(file=fname, exist=found)
           if (.not.found) then
             write(6,'(/,2a,a20)') 'read_vps: WARNING: ',
     .         'Pseudopotential file not found: ', fname
             fname = paste(atm_label,'.psatom.data')
             write(6,'(2a)') 'read_vps: WARNING: Looking for ', fname
             inquire(file=fname, exist=found)
           endif
           if (.not.found) then
             write(6,'(/,2a,a20,/)') 'read_vps: ERROR: ',
     .         'Pseudopotential file not found: ', fname
             stop 'read_vps: ERROR: Pseudopotential file not found'
           endif

           open(unit=1,file=fname,form='unformatted',status='unknown')

           read(1) namatm, icorr, irel, nicore,
     .     (method(i),i=1,6), text,
     .     npotd, npotu, nr, b, a, zval


           linput=max(lmxo,lmxkb)
           lmax=min(npotd-1,linput)


           if (lmax.lt.linput) then
             write(6,'(a)')
     .        'read_vps: ERROR: You must generate a pseudopotential'
             write(6,'(a)')
     .        'read_vps: ERROR: for each L up to ',linput
             stop
           endif

           nrval=nr+1
           if(nrval.gt.nrmax) then
              write(6,'(a, i4)')
     .       'read_vps: ERROR: Nrmax must be increased to at least',
     .          nrval
              stop
           endif


C********Write information about pseudopotential generation*******************


           write(6,'(/,a)')
     .           'read_vps: Pseudopotential generation method:'
           write(6,'(7a)')   
     .            'read_vps: ',method(1),(method(i),i=3,6)
          

 
          write(6,'(/,a)') 'read_vps: Valence configuration '//
     .                             '(pseudopotential generation):' 

            do 10 l=0,lmax
              itext=l*17
              write(6,*)  text(1+itext:2+itext),
     .               '(',text(3+itext:7+itext),')'
c    .        ,'   rc=', text(12+itext:16+itext)
   10       continue



           if(irel.eq.'rel') then
              write(6,'(/,2a)') 
     .          'read_vps: Pseudopotential generated from a ',
     .          'relativistic atomic calculation'
              write(6,'(2a)') 
     .          'read_vps: There are spin-orbit pseudopotentials ',
     .          'available'
              write(6,'(2a)') 
     .          'read_vps: Spin-orbit interaction is not included in ',
     .          'this calculation'
           elseif (irel.eq.'isp') then
              write(6,'(/,2a)')
     .          'read_vps: Pseudopotential generated from an ',
     .          'atomic spin-polarized calculation'
           endif


           if (nicore.ne.'nc  ') then
           write(6,'(/,a)')
     .      'read_vps: Pseudopotential includes a core correction:'

             if(nicore.eq.'pcec') then
               write(6,'(a)') 'read_vps: Pseudo-core for xc-correction'
             elseif(nicore.eq.'pche') then
               write(6,'(a)') 
     .            'read_vps: Pseudo-core for hartree and xc-correction'
             elseif(nicore.eq.'fcec') then
               write(6,'(a)') 'read_vps: Full-core for xc-correction'
             elseif(nicore.eq.'fche') then
               write(6,'(a)') 
     .            'read_vps: Full-core for hartree and xc-correction'
             endif

           endif

           linput=max(lmxo,lmxkb)
           lmax=min(npotd-1,linput)


           if (lmax.lt.linput) then
             write(6,'(a)') 
     .        'read_vps: ERROR: You must generate a pseudopotential'
             write(6,'(a,i4)') 
     .        'read_vps: ERROR: for each L up to ',linput
             stop
           endif

           nrval=nr+1 
           nodd=mod(nrval,2)
           nrval=nrval-1+nodd

           if(nrval.gt.nrmax) then
              write(6,'(a,i4)')
     .     'read_vps: ERROR: Nrmax must be increased to at least',nrval
              stop
           endif


C*******************Radial mesh*****************************************

           read(1) (rofi(ir),ir=2,nrval)
           rofi(1)=0.0d0 



C********Calculate drdi and s ***********************************************
C****drdi is the derivative of the radial distance respect to the mesh index
C*****i.e. rofi(ir)= b*[ exp( a*(i-1) ) - 1 ] and therefore *****************
C*****drdi=dr/di =a*b*exp(a*(i-1))= a*[rofi(ir)+b] **************************

           rpb=b
           ea=dexp(a)
           do ir=1,nrval
             drdi(ir)=a*rpb
             s(ir)=dsqrt(a*rpb)
             rpb=rpb*ea
           enddo 

C*****************************************************************************




C********Reading ionic pseudopotentials*******************
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








          return 
          
         end





               subroutine comKB(is,a,b,rofi,rphi,
     .            l,rc,dkbcos,nrc)

C***********************************************************************
C  Creates the common block with all the information about the 
C  Kleinman-Bylander projectors.
C  Written by D. Sanchez-Portal, Aug. 1998.
C***********************************************************************
              
               implicit none
               include 'atom.h'

               integer l, nrc,is

               double precision rc, dkbcos, rphi(nrmax), a, b,
     .            rofi(nrmax)  

C*********Internal variables*********************************************
C
             integer indx, itb, nr, nmax, nmin, nn  
        
             double precision delt, r, vphi, dy, yp1, ypn
             
             double precision
     .          aux(ntbmax)
 
             double precision  deltmax 
             common/cmdelt/deltmax

****NUMBER OF POINTS USED BY RATINT FOR THE INTERPOLATION*************
C
          integer npoint
          parameter(npoint=4)
C
C***********************************************************************

C***********Variables in common blocks**********************************
C
          double precision
     .  table((ntbmax+2),-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tab2(ntbmax,-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tabpol((ntbmax+2),nzetmx*lmaxd,nsmax),
     .  tab2pol(ntbmax,nzetmx*lmaxd,nsmax),
     .  rctb(0:lmaxd,nsmax)
     

C***********************************************************************
C
              common/cmtab/table,tabpol
              common/cmspline/tab2,tab2pol
              common/cmradkb/rctb
C
C***********************************************************************

             rctb(l,is)=rc 
   
C**********INTERPOLATION TO GENERATE TABLES WITH KB PROJECTORS*******
C
             indx=l+1

             delt=rc/(dble(ntbmax-1)+1.0d-20) 
          if(delt.gt.deltmax) then
              write(6,'(a)')
     .    'comKB: WARNING It might be a good idea to increase'
              write(6,'(a)')
     .    'comKB: WARNING parameter ntbmax (in file atom.h) '
              write(6,'(a,i6)')
     .    'comKB: WARNING to at least ntbmax = ', 
     .        nint(Rc/deltmax)+2
           endif

             table(1,-indx,is)=delt
             table(2,-indx,is)=dkbcos

             do itb=1,ntbmax-1
                r=delt*(itb-1)
                nr=nint(dlog(r/b+1.0d0)/a)+1
                nmin=max(1,nr-npoint)
                nmax=min(nrc,nr+npoint)
                nn=nmax-nmin+1
                call ratint(rofi(nmin),rphi(nmin),nn,r,vphi,dy)
                table(itb+2,-indx,is)=vphi
            enddo 
            table(ntbmax+2,-indx,is)=0.0d0
C
C*********************************************************************

C*********TABLE WITH THE SECOND DERIVATIVE ****************************
C
            yp1=1.0d50
            ypn=1.0d50

            call spline(delt,table(3,-indx,is),ntbmax,
     .        yp1,ypn,tab2(1,-indx,is),aux)

C
C*********************************************************************



             return

              end




               subroutine KBgen(is, a,b,rofi,drdi,s, 
     .         vps, vlocal, ve, nrval, Zval, lmxkb, nkb)

C*********************************************************************
C Call routines for 1) the generation of the Kleinman-Bylander projectors,
C 2) Cheking for the presence of ghost states and 3), the storage of 
C all the inforamtion in the corresponding common blocks.
C
C  Written D. Sanchez-Portal, Aug. 1998.
C*********************************************************************


               implicit none
               include 'atom.h'

               double precision 
     .            a, b, rofi(nrmax), vps(nrmax,0:lmaxd),
     .            drdi(nrmax), s(nrmax), ve(nrmax),vlocal(nrmax),
     .            Zval

               integer
     .           nrval, lmxkb, nkb, is



C********Internal variables*************************************

               integer 
     .           l,nprin, nnodes, ighost, nrc 
     
               double precision
     .           eigen(0:lmaxd), rc(0:lmaxd), dkbcos(0:lmaxd) 
               
               double precision
     .           rphi(nrmax,0:lmaxd)



         do l=0,lmxkb


C***************Atomic wavefunctions and eigenvalues****************
C
            nnodes=1
            nprin=l+1
            call schro_eq(Zval,rofi,vps(1,l),ve,s,drdi,
     .        nrval,l,a,b,nnodes,nprin,
     .        eigen(l),rphi(1,l))
C
C******************************************************************

C******************GHOST ANALYSIS**********************************
C
             call ghost(Zval,rofi,vps(1,l),vlocal,
     .        ve,s,drdi,nrval,l,a,b,
     .        eigen(l),rphi(1,l),ighost)


C******************KB Projectors***********************************
C
          call KBproj(rofi,drdi,vps(1,l),vlocal,nrval,l,
     .    dkbcos(l),rphi(1,l),nrc)
C
C******************************************************************

          rc(l)=rofi(nrc)
 
C*****Common block with the information about the  KB projectors********
C
          call comKB(is,a,b,rofi,rphi(1,l),
     .            l,rc(l),dkbcos(l),nrc)
C
C********************************************************************


         enddo   


         if (ighost.eq.1) then
            write(6,"(2a)")'KBgen: WARNING: ',
     .            'Ghost states have been detected'
            write(6,"(2a)")'KBgen: WARNING: ',
     .            'Some parameter should be changed in the '
            write(6,"(2a)")'KBgen: WARNING: ',
     .            'pseudopotential generation procedure.'
            stop
          endif


         write(6,'(/,a)')'KBgen: Kleinman-Bylander projectors: '
         do l=0,lmxkb
              write(6,'(3x,a,i2,3(3x,a,f10.6))')
     .        'l=',l, 'rc=',rc(l), 'el=',eigen(l), 'kbcos=',dkbcos(l)
         enddo


C**********TOTAL NUMBER OF KLEINMAN-BYLANDER PROJECTORS**************
C
           nkb=(lmxkb+1)*(lmxkb+1) 
           write(6,'(/,a, i4)')
     .'KBgen: Total number of  Kleinman-Bylander projectors: ', nkb
C
C********************************************************************


        return
   
        end 








              subroutine Basis_gen(Zval,is, a,b,rofi,drdi,s,
     .                   vps, ve, vePAO, nrval, lmxo, 
     .                   nzeta, rco, lambda, polorb,
     .                   basis_type, rphi, notot)

C*********************************************************************
C Generates the basis set and stores all the information in the 
C correspoinding common blocks.
C
C Written by D. Sanchez-Portal, Aug. 1998.
C*********************************************************************


               implicit none
               include 'atom.h'

               double precision
     .            a, b, rofi(nrmax), vps(nrmax,0:lmaxd),
     .            drdi(nrmax), s(nrmax), ve(nrmax),
     .            rphi(nrmax,0:lmaxd), rco(nzetmx,0:lmaxd),
     .            lambda(nzetmx, 0:lmaxd), vePAO(nrmax),
     .            Zval

               integer
     .           nrval, lmxo, notot, is, nzeta(0:lmaxd),
     .           polorb(lmaxd)

               character
     .           basis_type*10

C********Internal variables*************************************

                integer noPAO, noPOL

                double precision ePAO(0:lmaxd)

             
                noPAO=0
                noPOL=0
                if(basis_type.eq.'split') then  

                 call SPLIT(Zval,is, a,b,rofi,drdi,s,
     .                   vps, ve, vePAO, 
     .                   nrval, lmxo,
     .                   nzeta, rco, lambda,
     .                   rphi, ePAO, noPAO)
         
                elseif(basis_type.eq.'nodes') then 
 
                 call NODES(Zval,is, a,b,rofi,drdi,s,
     .                   vps, ve, vePAO,
     .                   nrval, lmxo,
     .                   nzeta, rco, lambda,
     .                   rphi, ePAO, noPAO)

                elseif(basis_type.eq.'nonodes') then 
 
                 call NONODES(Zval,is, a,b,rofi,drdi,s,
     .                   vps, ve, vePAO,
     .                   nrval, lmxo,
     .                   nzeta, rco, lambda,
     .                   rphi, ePAO, noPAO)
          
                elseif(basis_type.eq.'splitgauss') then 
      
                 call SPLITGAUSS(Zval,is, a,b,rofi,drdi,s,
     .                   vps, ve, vePAO,
     .                   nrval, lmxo,
     .                   nzeta, rco, lambda,
     .                   rphi, ePAO, noPAO)

                elseif(basis_type.eq.'user') then 
 
                 call USER(is, a, b, rofi, drdi,
     .             vps, ve, lmxo, nzeta,
     .             rco, lambda, rphi, ePAO, noPAO)

                endif            
 

C********Polarization orbitals*******************************
C 
                  call POLgen(is,a,b,rofi,drdi,
     .               ePAO,rphi,rco,vps,vePAO,
     .               polorb,lmxo,noPOL) 
C
C************************************************************

C******Total number of orbitals******************************
C 
                   notot=noPAO+noPOL
C
C************************************************************

        

              return
              end






           subroutine SPLIT(Zval,is,a,b,rofi,drdi,s,
     .             vps,ve,vePAO,
     .             nrval,lmxo,
     .             nzeta,rco,lambda, rphi, ePAO, norb) 

C*********************************************************************
C Calculates the atomic orbitals basis set, using the option SPLIT 
C for the generation of the augmentation orbitals.
C  Written by D. Sanchez-Portal, Aug. 1998
C*********************************************************************

               implicit none
               include 'atom.h'

               double precision
     .         a, b, rofi(nrmax), vps(nrmax,0:lmaxd),
     .         drdi(nrmax), s(nrmax), ve(nrmax),
     .         rphi(nrmax,0:lmaxd), rco(nzetmx,0:lmaxd),
     .         lambda(nzetmx, 0:lmaxd), Zval,vePAO(nrmax),
     .         ePAO(0:lmaxd)


               integer
     .           nrval, lmxo, is, nzeta(0:lmaxd),
     .           norb


C********Internal variables*************************************

               integer
     .           l,nprin, nnodes, nodd, nrc, nsp, i, ir,indx,
     .           izeta, nmax, nmin, nn, nr, nrcomp

               double precision
     .           eigen(0:lmaxd), rc,
     .           eshift_default, splnorm_default,
     .           rnrm(nrmax), dnrm, phi, frsp, dfrsp,
     .           cons1, cons2, rnp, spln, eshift, 
     .           splnorm, g(nrmax), r, el, ekin, 
     .           r1, r2, dfdi, d2fdi2, d2fdr2, dr,
     .           epot, epot2, rh, dy, eorb, eps



C****NUMBER OF POINTS USED BY RATINT FOR THE INTERPOLATION*************
C 
               integer  npoint 
               parameter(npoint=4)
C
C***********************************************************************


C***********Declarations for FDF package *******************************
C
         include 'fdf/fdfdefs.h'
C
C***********************************************************************

C******Common blocks with the values of some parameters*********
C
        common/cmeshdefault/eshift_default  
        common/cmspldefault/splnorm_default
C
C***************************************************************       

      


C***READING THE ENERGY-SHIFT TO DEFINE THE CUT-OFF RADIUS OF ORBITALS***

         eshift=fdf_physical('PAO.EnergyShift',eshift_default,'Ry')

C***********************************************************************


C***READING SPLNORM TO GENERATE THE SPLIT IF Rmatch IS ZERO IN INPUT****

         splnorm=fdf_double('PAO.SplitNorm',splnorm_default)

C***********************************************************************





 
             norb=0 
             indx=0
             do l=0,lmxo 

               if(nzeta(l).gt.0) then

            write(6,'(/A,I2)')
     .       'SPLIT: Orbitals with angular momentum L=',l



                  if(rco(1,l).lt.1.0d-5) then    
C**Automatic determination of the cut off radius for the PAOs********
C***************Atomic eigenvalues***********************************
C
                      nnodes=1
                      nprin=l+1
                      call schro_eq(Zval,rofi,vps(1,l),ve,s,drdi,
     .                  nrval,l,a,b,nnodes,nprin,
     .                  eigen(l),rphi(1,l))
C
C****************Rc given by eshift**********************************   
C
                       if(eigen(l).gt.0.0d0) then 
                          write(6,'(/A,I2,A)')
     .       'SPLIT: ERROR Orbital with angular momentum L=',l,
     .       ' not bound in the atom'
                         write(6,'(A)')
     .       'SPLIT: ERROR a cut off radius must be explicitely given' 
                          stop
                       endif 
                       if(dabs(eshift).gt.1.0d-5) then
                          el=eigen(l)+eshift
                          call rc_vs_e(a,b,rofi,vps(1,l),
     .                                ve,nrval,l,el,rco(1,l))
c                               write(6,'(/,A,I2,A,f10.6,A)')
                       else
                          rco(1,l)=rofi(nrval-2)
                       endif
                  write(6,'(/A/A,f10.6,A)')
     .        'SPLIT: PAO cut-off radius determinated from an',
     .        'SPLIT: energy shift=',eshift,' Ry'

                 endif  
C
C********************************************************************


C*****IF THE COMPRESSION FACTOR IS NEGATIVE OR ZERO THE ORBITALS ARE***
C**************LEFT UNTOUCHED******************************************
             if(lambda(1,l).le.0.0d0) lambda(1,l)=1.0d0
C*********************************************************************** 




              do izeta=1, nzeta(l)

C****COMPRESSION FACTOR IS ONLY ACTIVE FOR THE INITIAL PAO WHEN USING****
C**** SPLIT OPTION FOR THE GENERATION OF THE BASIS SET*******************
                 lambda(izeta,l)=lambda(1,l)
C************************************************************************

                  rc=rco(izeta,l)/lambda(1,l)
                  nrc=nint(dlog(rc/b+1.0d0)/a)+1
                  nodd=mod(nrc,2)
                  if(nodd.eq.0) then
                     nrc=nrc+1
                  endif
                  rc=b*(dexp(a*(nrc-1))-1.0d0)
                  rco(izeta,l)=rc*lambda(1,l)
                  
                if(izeta.eq.1) then 
C****Generate PAO orbitals for the first shell of the basis set********
C 
                      nnodes=1
                      nprin=l+1
                      call schro_eq(Zval,rofi,vps(1,l),vePAO,s,drdi,
     .                  nrc,l,a,b,nnodes,nprin,
     .                  eorb,rphi(1,l)) 
                       dnrm=0.0d0
                       do ir=2,nrc
                          phi=rphi(ir,l)
                          dnrm=dnrm+drdi(ir)*phi*phi
                          rnrm(ir)=dnrm 
                          g(ir)=rphi(ir,l)/(rofi(ir)**(l+1))
                       enddo 
                       g(1)=g(2)         
                elseif(izeta.gt.1) then 

       
C********Cut-off radius for double-Z, triple-Z,..., if it is set to**** 
C********zero in the input then it is calculated from the splitnorm**** 
C******** parameter ***************************************************

                if(rco(izeta,l).gt.rco(1,l)) then
                  write(6,'(/,A)') 
     . 'SPLIT: ERROR: SPLIT OPTION FOR BASIS SET '
                  write(6,'(A)')  
     . 'SPLIT: ERROR: Rc FOR DOUBLE-Z, TRIPLE-Z,... SHOULD BE SMALLER '
                  write(6,'(A)') 
     . 'SPLIT: ERROR:  THAN THAT OF THE INITIAL PAO !!!!!'
                  STOP
                endif
                                  
                
            if(rco(izeta,l).gt.1.0d-5) then
               rc=rco(izeta,l)/lambda(1,l)
               frsp=rphi(nrc,l)/rc
               dfrsp=0.5d0*(rphi(nrc+1,l)/rofi(nrc+1)
     .             -rphi(nrc-1,l)/rofi(nrc-1))
               dfrsp=dfrsp/drdi(nrc)


C**********************parabolic split******************************
            cons1= 0.5d0*(dfrsp*rc-l*frsp)/(rc**(l+2))
            cons2= frsp/(rc**l)-cons1*rc**2
            call nrmpal(cons1,cons2,rc,l,rnp)
            spln=1.0d0-rnrm(nrc)+rnp
C*******************************************************************

                do i=1,izeta-1
                 if(dabs(rco(izeta,l)-rco(i,l)).lt.1.0d-5) then
                   write(6,'(/,A,I2,A,I2,A,I2)')
     .            'SPLIT: WARNING: Split-orbital with zeta=',izeta,
     .            ' and zeta=',i,' are identical for l=',l
                 endif
                enddo
            else 

            rc=rco(1,l)/lambda(1,l)
            nrc=nint(dlog(rc/b+1.0d0)/a)+1
            spln=splnorm
            if(izeta.gt.2) then
              spln=spln/(2.0d0*(izeta-2) )
            endif

            call parabola(a,b,nrc,rphi(1,l),rnrm,
     .                   l,spln,cons1,cons2,nsp)


C***Cut-off radius for the split orbital with a desired norm******
         nrc=nsp
         rco(izeta,l)=b*(dexp(a*(nsp-1))-1.0d0)*lambda(izeta,l)
C*****************************************************************


             do i=1,izeta-1
                if(dabs(rco(izeta,l)-rco(i,l)).lt.1.0d-5) then
                   write(6,'(/,A,I2,A,I2,A,I2)')
     .            'SPLIT: WARNING: Split-orbital with zeta=',izeta,
     .            ' and zeta=',i,' are identicals for l=',l
                endif
             enddo

            endif
              
            dnrm=0.0d0
            do ir=2,nrc-1
              r=rofi(ir)
C***********************parabolic split****************************
              phi=-(cons1*r**2+cons2)+rphi(ir,l)/(r**(l+1))
C******************************************************************

              dnrm=dnrm+drdi(ir)*(phi*r**(l+1))**2
              g(ir)=phi 
             
            enddo
            g(1)=g(2)
            g(nrc)=0.0d0
            
              endif

C**************Normalization of basis functions***********************
            eps=1.0d-4
            if(dabs(dnrm-1.0d0).gt.eps) then
               do ir=1,nrc
                 g(ir)=g(ir)/dsqrt(dnrm)
                 if(izeta.eq.1) then
                    rphi(ir,l)=rphi(ir,l)/dsqrt(dnrm)
                    rnrm(ir)=rnrm(ir)/dnrm
                 endif
               enddo
            endif
C*********************************************************************  

C*Calculation of the mean value of kinetic and potential energy*******
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

             ekin=ekin/(lambda(izeta,l)**2)

C Potential energy after compression

             nrcomp=nint(dlog(rco(izeta,l)/b+1.0d0)/a)+1
             epot=0.0d0
             epot2=0.0d0
             do ir=1,nrcomp
                r=rofi(ir)
                r2=r/lambda(izeta,l)
                nr=nint(dlog(r2/b+1.0d0)/a)+1
                nmin=max(1,nr-npoint)
                nmax=min(nrc,nr+npoint)
                nn=nmax-nmin+1
                call ratint(rofi(nmin),g(nmin),nn,r2,rh,dy)
                rh=rh/dsqrt(lambda(izeta,l)**(2*l+3))
                epot=epot+
     .          drdi(ir)*(ve(ir)+vps(ir,l))*(rh*r**(l+1))**2
                epot2=epot2+
     .          drdi(ir)*vps(ir,l)*(rh*r**(l+1))**2
             enddo
             eorb=ekin+epot



            if(izeta.eq.1) then  

             write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')
     .          'izeta =',izeta,
     .          'lambda =',lambda(izeta,l),
     .          'rc =',rco(izeta,l),
     .          'energy =',eorb  

                ePAO(l)=eorb

            elseif(izeta.gt.1) then 

            write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')
     .         'izeta =',izeta,
     .         'rmatch =',rco(izeta,l),
     .         'splitnorm =',spln,
     .         'energy =',eorb 

            endif 

          write(6,'(a25,f12.6)') 'kinetic =',ekin
          write(6,'(a25,f12.6)') 'potential(screened) =',epot
          write(6,'(a25,f12.6)') 'potential(ionic) =',epot2 

            norb=norb+(2*l+1)
            indx=indx+1
            call comBasis(is,a,b,rofi,g,l,
     .              rco(izeta,l),lambda(izeta,l),izeta,
     .              nrc,indx)


              enddo 
        
            call compress_PAO(a,b,rofi,rphi(1,l),
     .              rco(1,l),lambda(1,l))  
                       

             endif     
             enddo 

 
            return
      
            end




                






           subroutine NODES(Zval,is,a,b,rofi,drdi,s,
     .             vps,ve,vePAO,
     .             nrval,lmxo,
     .             nzeta,rco,lambda, rphi, ePAO, norb) 

C*********************************************************************
C Calculates the atomic orbitals basis set, using the option NODES
C for the generation of the augmentation orbitals.
C  Written by D. Sanchez-Portal, Aug. 1998
C*********************************************************************

               implicit none
               include 'atom.h'

               double precision
     .         a, b, rofi(nrmax), vps(nrmax,0:lmaxd),
     .         drdi(nrmax), s(nrmax), ve(nrmax),
     .         rphi(nrmax,0:lmaxd), rco(nzetmx,0:lmaxd),
     .         lambda(nzetmx, 0:lmaxd), Zval,vePAO(nrmax),
     .         ePAO(0:lmaxd)

               integer
     .           nrval, lmxo, is, nzeta(0:lmaxd),
     .           norb


C********Internal variables*************************************

               integer
     .           l,nprin, nnodes, nodd, nrc, ir, indx,
     .           izeta, nmax, nmin, nn, nr, nrcomp

               double precision
     .           eigen(0:lmaxd), rc,
     .           eshift_default, 
     .           dnrm, phi, eshift,
     .           g(nrmax), r, el, ekin, 
     .           r1, r2, dfdi, d2fdi2, d2fdr2, dr,
     .           epot, epot2, rh, dy, eorb, eps




C****NUMBER OF POINTS USED BY RATINT FOR THE INTERPOLATION*************
C 
               integer  npoint 
               parameter(npoint=4)
C
C***********************************************************************


C***********Declarations for FDF package *******************************
C
         include 'fdf/fdfdefs.h'
C
C***********************************************************************

C******Common blocks with the values of some parameters*********
C
        common/cmeshdefault/eshift_default  
C
C***************************************************************       

      


C***READING THE ENERGY-SHIFT TO DEFINE THE CUT-OFF RADIUS OF ORBITALS***

         eshift=fdf_physical('PAO.EnergyShift',eshift_default,'Ry')

C***********************************************************************


 
             norb=0 
             indx=0
             do l=0,lmxo

               if(nzeta(l).gt.0) then 

            write(6,'(/A,I2)')
     .       'NODES: Orbitals with angular momentum L=',l






                  if(rco(1,l).lt.1.0d-5) then    
C**Automatic determination of the cut off radius for the PAOs********
C***************Atomic eigenvalues***********************************
C
                      nnodes=1
                      nprin=l+1
                      call schro_eq(Zval,rofi,vps(1,l),ve,s,drdi,
     .                  nrval,l,a,b,nnodes,nprin,
     .                  eigen(l),rphi(1,l))
C
C****************Rc given by eshift**********************************   
C                
                       if(eigen(l).gt.0.0d0) then
                          write(6,'(/A,I2,A)')
     .       'NODES: ERROR Orbital with angular momentum L=',l,
     .       ' not bound in the atom'
                         write(6,'(A)')
     .       'NODES: ERROR a cut off radius must be explicitely given' 
                          stop
                       endif 
 
                       if(dabs(eshift).gt.1.0d-5) then
                          el=eigen(l)+eshift
                          call rc_vs_e(a,b,rofi,vps(1,l),
     .                                ve,nrval,l,el,rco(1,l))
                       else
                          rco(1,l)=rofi(nrval-2)
                       endif 

                  write(6,'(/,A,/,A,f10.6,A)')
     .         'NODES: PAO cut-off radius determinated from an',
     .         'NODES: energy shift=',eshift,' Ry'

                 endif  
C
C********************************************************************


              do izeta=1, nzeta(l)

C*****IF THE COMPRESSION FACTOR IS NEGATIVE OR ZERO THE ORBITALS ARE***
C********************UNTOUCHED******************************************
             if(lambda(izeta,l).le.0.0d0) lambda(izeta,l)=1.0d0 
C***********************************************************************
             if(dabs(rco(izeta,l)).le.1.0d-5) then 
                 rco(izeta,l)=rco(1,l)
             endif 

                  rc=rco(izeta,l)/lambda(izeta,l)
                  nrc=nint(dlog(rc/b+1.0d0)/a)+1
                  nodd=mod(nrc,2)
                  if(nodd.eq.0) then
                     nrc=nrc+1
                  endif
                  rc=b*(dexp(a*(nrc-1))-1.0d0)
                  rco(izeta,l)=rc*lambda(izeta,l)
C****Generate PAO orbitals with increasing number of nodes *************
C*********************for the different shells**************************
C 
                      nnodes=izeta
                      nprin=l+izeta 
                      eorb=0.0d0 
                   call schro_eq(Zval,rofi,vps(1,l),vePAO,s,drdi,
     .              nrc,l,a,b,nnodes,nprin,
     .                eorb,g)  

                       dnrm=0.0d0
                       do ir=2,nrc
                          phi=g(ir) 
                          if(izeta.eq.1) rphi(ir,l)=phi
                          dnrm=dnrm+drdi(ir)*phi*phi
                          g(ir)=g(ir)/(rofi(ir)**(l+1))
                       enddo 
                       g(1)=g(2)         

C**************Normalization of basis functions***********************
            eps=1.0d-4
            if(dabs(dnrm-1.0d0).gt.eps) then
               do ir=1,nrc
                 g(ir)=g(ir)/dsqrt(dnrm)
                 if(izeta.eq.1) then
                    rphi(ir,l)=rphi(ir,l)/dsqrt(dnrm)
                 endif
               enddo
            endif
C*********************************************************************


C*Calculation of the mean value of kinetic and potential energy*******
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

             ekin=ekin/(lambda(izeta,l)**2)

C Potential energy after compression

             nrcomp=nint(dlog(rco(izeta,l)/b+1.0d0)/a)+1
             epot=0.0d0
             epot2=0.0d0
             do ir=1,nrcomp
                r=rofi(ir)
                r2=r/lambda(izeta,l)
                nr=nint(dlog(r2/b+1.0d0)/a)+1
                nmin=max(1,nr-npoint)
                nmax=min(nrc,nr+npoint)
                nn=nmax-nmin+1
                call ratint(rofi(nmin),g(nmin),nn,r2,rh,dy)
                rh=rh/dsqrt(lambda(izeta,l)**(2*l+3))
                epot=epot+
     .          drdi(ir)*(ve(ir)+vps(ir,l))*(rh*r**(l+1))**2
                epot2=epot2+
     .          drdi(ir)*vps(ir,l)*(rh*r**(l+1))**2
             enddo
             eorb=ekin+epot



             write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')
     .          'izeta =',izeta,
     .          'lambda =',lambda(izeta,l),
     .          'rc =',rco(izeta,l),
     .          'energy =',eorb  


                 if(izeta.eq.1) ePAO(l)=eorb

          write(6,'(a25,f12.6)') 'kinetic =',ekin
          write(6,'(a25,f12.6)') 'potential(screened) =',epot
          write(6,'(a25,f12.6)') 'potential(ionic) =',epot2 

            norb=norb+(2*l+1)
            indx=indx+1
            call comBasis(is,a,b,rofi,g,l,
     .              rco(izeta,l),lambda(izeta,l),izeta,
     .              nrc,indx)




              enddo 
        
            call compress_PAO(a,b,rofi,rphi(1,l),
     .              rco(1,l),lambda(1,l))




              endif    
             enddo 

 
            return
      
            end




                






           subroutine NONODES(Zval,is,a,b,rofi,drdi,s,
     .             vps,ve,vePAO,
     .             nrval,lmxo,
     .             nzeta,rco,lambda, rphi, ePAO, norb) 

C*********************************************************************
C Calculates the atomic orbitals basis set, using the option NONODES
C for the generation of the augmentation orbitals.
C  Written by D. Sanchez-Portal, Aug. 1998
C*********************************************************************

               implicit none
               include 'atom.h'

               double precision
     .         a, b, rofi(nrmax), vps(nrmax,0:lmaxd),
     .         drdi(nrmax), s(nrmax), ve(nrmax),
     .         rphi(nrmax,0:lmaxd), rco(nzetmx,0:lmaxd),
     .         lambda(nzetmx, 0:lmaxd), Zval,vePAO(nrmax),
     .         ePAO(0:lmaxd)

               integer
     .           nrval, lmxo,  is, nzeta(0:lmaxd),
     .           norb


C********Internal variables*************************************

               integer
     .           l,nprin, nnodes, nodd, nrc, i, ir, indx,
     .           izeta, nmax, nmin, nn, nr, nrcomp

               double precision
     .           eigen(0:lmaxd), rc,
     .           eshift_default, 
     .           dnrm, phi, eshift,
     .           g(nrmax), r, el, ekin, 
     .           r1, r2, dfdi, d2fdi2, d2fdr2, dr,
     .           epot, epot2, rh, dy, eorb, eps




C****NUMBER OF POINTS USED BY RATINT FOR THE INTERPOLATION*************
C 
               integer  npoint 
               parameter(npoint=4)
C
C***********************************************************************


C***********Declarations for FDF package *******************************
C
         include 'fdf/fdfdefs.h'
C
C***********************************************************************

C******Common blocks with the values of some parameters*********
C
        common/cmeshdefault/eshift_default  
C
C***************************************************************       

      


C***READING THE ENERGY-SHIFT TO DEFINE THE CUT-OFF RADIUS OF ORBITALS***

         eshift=fdf_physical('PAO.EnergyShift',eshift_default,'Ry')

C***********************************************************************



 
             norb=0 
             indx=0
             do l=0,lmxo
               if(nzeta(l).gt.0) then 
               
            write(6,'(/A,I2)')
     .       'NONODES: Orbitals with angular momentum L=',l






                  if(rco(1,l).lt.1.0d-5) then    
C**Automatic determination of the cut off radius for the PAOs********
C***************Atomic eigenvalues***********************************
C
                      nnodes=1
                      nprin=l+1
                      call schro_eq(Zval,rofi,vps(1,l),ve,s,drdi,
     .                  nrval,l,a,b,nnodes,nprin,
     .                  eigen(l),rphi(1,l))
C
C****************Rc given by eshift**********************************   
C                
                       if(eigen(l).gt.0.0d0) then
                          write(6,'(/A,I2,A)')
     .       'NONODES: ERROR Orbital with angular momentum L=',l,
     .       ' not bound in the atom'
                         write(6,'(A)')
     .       'NONODES: ERROR a cut off radius must be explicitely given' 
                          stop
                       endif 
 
                       if(dabs(eshift).gt.1.0d-5) then
                          el=eigen(l)+eshift
                          call rc_vs_e(a,b,rofi,vps(1,l),
     .                                ve,nrval,l,el,rco(1,l))
                       else
                          rco(1,l)=rofi(nrval-2)
                       endif 

                  write(6,'(/,A,/,A,f10.6,A)')
     .         'NONODES: PAO cut-off radius determinated from an',
     .         'NONODES: energy shift=',eshift,' Ry'

                 endif  
C
C********************************************************************


              do izeta=1, nzeta(l)

C*****IF THE COMPRESSION FACTOR IS NEGATIVE OR ZERO THE ORBITALS ARE***
C********************UNTOUCHED******************************************
             if(lambda(izeta,l).le.0.0d0) lambda(izeta,l)=1.0d0 
C*********************************************************************** 
            if(dabs(rco(izeta,l)).lt.1.0d-5) rco(izeta,l)=rco(1,l)
            do i=1,izeta-1
             if((dabs(rco(izeta,l)-rco(i,l)).lt.1.0d-5).and.
     . (dabs(lambda(izeta,l)-lambda(i,l)).lt.1.0d-5)) then
               write(6,'(/,A,I2,A,I2,A,I2)')
     .  'NONODES: WARNING: PAO base function with zeta=',izeta,
     .  ' and zeta=',i,' are identical for l=',l   
             stop   
             endif
            enddo

                  rc=rco(izeta,l)/lambda(izeta,l)
                  nrc=nint(dlog(rc/b+1.0d0)/a)+1
                  nodd=mod(nrc,2)
                  if(nodd.eq.0) then
                     nrc=nrc+1
                  endif
                  rc=b*(dexp(a*(nrc-1))-1.0d0)
                  rco(izeta,l)=rc*lambda(izeta,l)
C****Generate PAO orbitals with increasing number of nodes *************
C*********************for the different shells**************************
C 
                     nnodes=1
                     nprin=l+1
                   call schro_eq(Zval,rofi,vps(1,l),vePAO,s,drdi,
     .              nrc,l,a,b,nnodes,nprin,
     .               eorb,g)  


                       dnrm=0.0d0
                       do ir=2,nrc
                          phi=g(ir) 
                          if(izeta.eq.1) rphi(ir,l)=phi
                          dnrm=dnrm+drdi(ir)*phi*phi
                          g(ir)=g(ir)/(rofi(ir)**(l+1))
                       enddo 
                       g(1)=g(2)         

C**************Normalization of basis functions***********************
            eps=1.0d-4
            if(dabs(dnrm-1.0d0).gt.eps) then
               do ir=1,nrc
                 g(ir)=g(ir)/dsqrt(dnrm)
                 if(izeta.eq.1) then
                    rphi(ir,l)=rphi(ir,l)/dsqrt(dnrm)
                 endif
               enddo
            endif
C*********************************************************************


C*Calculation of the mean value of kinetic and potential energy*******
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

             ekin=ekin/(lambda(izeta,l)**2)

C Potential energy after compression

             nrcomp=nint(dlog(rco(izeta,l)/b+1.0d0)/a)+1
             epot=0.0d0
             epot2=0.0d0
             do ir=1,nrcomp
                r=rofi(ir)
                r2=r/lambda(izeta,l)
                nr=nint(dlog(r2/b+1.0d0)/a)+1
                nmin=max(1,nr-npoint)
                nmax=min(nrc,nr+npoint)
                nn=nmax-nmin+1
                call ratint(rofi(nmin),g(nmin),nn,r2,rh,dy)
                rh=rh/dsqrt(lambda(izeta,l)**(2*l+3))
                epot=epot+
     .          drdi(ir)*(ve(ir)+vps(ir,l))*(rh*r**(l+1))**2
                epot2=epot2+
     .          drdi(ir)*vps(ir,l)*(rh*r**(l+1))**2
             enddo
             eorb=ekin+epot




             write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')
     .          'izeta =',izeta,
     .          'lambda =',lambda(izeta,l),
     .          'rc =',rco(izeta,l),
     .          'energy =',eorb  

                if(izeta.eq.1) ePAO(l)=eorb

          write(6,'(a25,f12.6)') 'kinetic =',ekin
          write(6,'(a25,f12.6)') 'potential(screened) =',epot
          write(6,'(a25,f12.6)') 'potential(ionic) =',epot2 

            norb=norb+(2*l+1)
            indx=indx+1
            call comBasis(is,a,b,rofi,g,l,
     .              rco(izeta,l),lambda(izeta,l),izeta,
     .              nrc,indx)


50           continue

              enddo 
               

            call compress_PAO(a,b,rofi,rphi(1,l),
     .              rco(1,l),lambda(1,l))



              endif    
             enddo 

 
            return
      
            end




                





               subroutine comBasis(is,a,b,rofi,rphi,
     .            l,rc,lambda,nzeta,nrc,norb)

C***********************************************************************
C  Generates the common blocks for the storage of the information 
C  about the basis set orbitals.
C
C  Written by D. Sanchez-Portal, Aug. 1998.
C***********************************************************************

               implicit none
               include 'atom.h'

               integer l, is, norb, nrc, nzeta

               double precision rc, lambda, rphi(nrmax), a, b,
     .            rofi(nrmax)  

C*********Internal variables*********************************************
C
             integer itb, nr, nmax, nmin, nn  
        
             double precision delt, r, phi, dy, yp1, ypn
             
             double precision
     .          aux(ntbmax) 
 
             double precision deltmax
             common/cmdelt/deltmax

****NUMBER OF POINTS USED BY RATINT FOR THE INTERPOLATION*************
C
          integer npoint
          parameter(npoint=4)
C
C***********************************************************************

C***********Variables in common blocks**********************************
C
          double precision
     .  table((ntbmax+2),-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tab2(ntbmax,-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tabpol((ntbmax+2),nzetmx*lmaxd,nsmax),
     .  tab2pol(ntbmax,nzetmx*lmaxd,nsmax),
     .  rcotb(nzetmx,0:lmaxd,nsmax),
     .  lambdatb(nzetmx,0:lmaxd,nsmax)

C***********************************************************************
C
              common/cmtab/table,tabpol
              common/cmspline/tab2,tab2pol
              common/cmradorb/rcotb
              common/cmlambda/lambdatb
C
C***********************************************************************    
              
              rcotb(nzeta,l,is)=rc
              lambdatb(nzeta,l,is)=lambda


C**********INTERPOLATION TO GENERATE TABLES WITH KB PROJECTORS*******
C
            delt=rc/(dble(ntbmax-1)+1.0d-20) 
  
            if(delt.gt.deltmax) then
              write(6,'(a)')
     .    'comBasis: WARNING It might be a good idea to increase'
              write(6,'(a)')
     .    'comBasis: WARNING parameter ntbmax (in file atom.h) '
              write(6,'(a,i6)')
     .    'comBasis: WARNING to at least ntbmax = ',
     .        nint(Rc/deltmax)+2
            endif

            table(1,norb,is)=delt

            table(2,norb,is)=dble(l)

            do itb=1,ntbmax-1
                 r=delt*(itb-1)
                 r=r/lambda
                 nr=nint(dlog(r/b+1.0d0)/a)+1
                 nmin=max(1,nr-npoint)
                 nmax=min(nrc,nr+npoint)
                 nn=nmax-nmin+1
                 call ratint(rofi(nmin),rphi(nmin),nn,r,phi,dy)
                 phi=phi/dsqrt(lambda**(2*l+3))
                 table(itb+2,norb,is)=phi
            enddo
C
C*********************************************************************

C*********TABLE WITH THE SECOND DERIVATIVE ****************************
C
            yp1=1.d50
            ypn=1.0d50

            call spline(delt,table(3,norb,is),ntbmax,
     .        yp1,ypn,tab2(1,norb,is),aux)

C
C*********************************************************************



             return

              end






           subroutine SPLITGAUSS(Zval,is,a,b,rofi,drdi,s,
     .             vps,ve,vePAO,
     .             nrval,lmxo,
     .             nzeta,rco,lambda, rphi, ePAO, norb) 

C*********************************************************************
C Calculates the atomic orbitals basis set, using the option SPLITGAUSS
C for the generation of the augmentation orbitals.
C  Written by D. Sanchez-Portal, Aug. 1998
C*********************************************************************

               implicit none
               include 'atom.h'

               double precision
     .         a, b, rofi(nrmax), vps(nrmax,0:lmaxd),
     .         drdi(nrmax), s(nrmax), ve(nrmax),
     .         rphi(nrmax,0:lmaxd), rco(nzetmx,0:lmaxd),
     .         lambda(nzetmx, 0:lmaxd), Zval,vePAO(nrmax),
     .         ePAO(0:lmaxd)

               integer
     .           nrval, lmxo, is, nzeta(0:lmaxd),
     .           norb


C********Internal variables*************************************

               integer
     .           l,nprin, nnodes, nodd, nrc, i, ir, indx,
     .           izeta, nmax, nmin, nn, nr, nrcomp

               double precision
     .           eigen(0:lmaxd), rc,
     .           eshift_default,
     .           dnrm, phi, 
     .           cons, fac, eshift, pi, gexp,
     .           g(nrmax), r, el, ekin, 
     .           r1, r2, dfdi, d2fdi2, d2fdr2, dr,
     .           epot, epot2, rh, dy, eorb, eps, dlapl




C****NUMBER OF POINTS USED BY RATINT FOR THE INTERPOLATION*************
C 
               integer  npoint 
               parameter(npoint=4)
C
C***********************************************************************


C***********Declarations for FDF package *******************************
C
         include 'fdf/fdfdefs.h'
C
C***********************************************************************

C******Common blocks with the values of some parameters*********
C
        common/cmeshdefault/eshift_default  
C
C***************************************************************       

      


C***READING THE ENERGY-SHIFT TO DEFINE THE CUT-OFF RADIUS OF ORBITALS***

         eshift=fdf_physical('PAO.EnergyShift',eshift_default,'Ry')

C***********************************************************************

C****************************PI*****************************************
                   pi=dacos(-1.0d0) 
C***********************************************************************


 
             norb=0 
             indx=0
             do l=0,lmxo
               if(nzeta(l).gt.0) then 
               
            write(6,'(/A,I2)')
     .       'SPLITGAUSS: Orbitals with angular momentum L=',l






                  if(rco(1,l).lt.1.0d-5) then    
C**Automatic determination of the cut off radius for the PAOs********
C***************Atomic eigenvalues***********************************
C
                      nnodes=1
                      nprin=l+1
                      call schro_eq(Zval,rofi,vps(1,l),ve,s,drdi,
     .                  nrval,l,a,b,nnodes,nprin,
     .                  eigen(l),rphi(1,l))
C
C****************Rc given by eshift**********************************   
C                
                       if(eigen(l).gt.0.0d0) then
                          write(6,'(/A,I2,A)')
     .  'SPLITGAUSS: ERROR Orbital with angular momentum L=',l,
     .       ' not bound in the atom'
                         write(6,'(A)')
     .  'SPLITGAUSS: ERROR a cut off radius must be explicitely given' 
                          stop
                       endif 
 
                       if(dabs(eshift).gt.1.0d-5) then
                          el=eigen(l)+eshift
                          call rc_vs_e(a,b,rofi,vps(1,l),
     .                                ve,nrval,l,el,rco(1,l))
                       else
                          rco(1,l)=rofi(nrval-2)
                       endif 

                  write(6,'(/,A,/,A,f10.6,A)')
     .   'SPLITGAUSS: PAO cut-off radius determinated from an',
     .   'SPLITGAUSS: energy shift=',eshift,' Ry'

                 endif  
C
C********************************************************************

C*****IF THE COMPRESSION FACTOR IS NEGATIVE OR ZERO THE ORBITALS ARE***
C********************UNTOUCHED******************************************
               if(lambda(1,l).le.0.0d0) lambda(1,l)=1.0d0
C***********************************************************************

              do izeta=1, nzeta(l)

             if(dabs(rco(izeta,l)).lt.1.0d-5) rco(izeta,l)=rco(1,l)
C**********With spligauss option, compression factor must be taken****
C**********as the gaussian exponent***********************************
              if(izeta.gt.1) then 
                  if(lambda(izeta,l).le.0.0d0) then 
                  write(6,'(/a,/a,a)')
     .'SPLITGAUSS: ERROR: with SPLITGAUSS option the compression ',
     .'SPLITGAUSS: ERROR: factors for all the augmentation functions',
     .   ' must be explicitely specified' 
                  stop
                  endif
                  gexp=dabs(lambda(izeta,l))
                  gexp=1.0d0/(gexp**2)
                  lambda(izeta,l)=1.0d0
              endif
C*********************************************************************
                  
                  rc=rco(izeta,l)/lambda(izeta,l)
                  nrc=nint(dlog(rc/b+1.0d0)/a)+1
                  nodd=mod(nrc,2)
                  if(nodd.eq.0) then
                     nrc=nrc+1
                  endif
                  rc=b*(dexp(a*(nrc-1))-1.0d0)
                  rco(izeta,l)=rc*lambda(izeta,l)  


                  if(izeta.eq.1) then 
C****Generate a PAO orbital for the first shell of basis functions****
C 
                      nnodes=1
                      nprin=l+1
                      call schro_eq(Zval,rofi,vps(1,l),vePAO,s,drdi,
     .                  nrc,l,a,b,nnodes,nprin,
     .                  eorb,rphi(1,l)) 
                       dnrm=0.0d0
                       do ir=2,nrc
                          phi=rphi(ir,l) 
                          dnrm=dnrm+drdi(ir)*phi*phi
                          g(ir)=phi/(rofi(ir)**(l+1))
                       enddo 
                       g(1)=g(2)         
                    elseif(izeta.gt.1) then 
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
                       phi=cons*dexp((-gexp)*r**2)
                       dnrm=dnrm+drdi(ir)*(phi*r**(l+1))**2
                       g(ir)=phi
                     enddo
                    endif  



C**************Normalization of basis functions***********************
            eps=1.0d-4
            if(dabs(dnrm-1.0d0).gt.eps) then
               do ir=1,nrc
                 g(ir)=g(ir)/dsqrt(dnrm)
                 if(izeta.eq.1) then
                    rphi(ir,l)=rphi(ir,l)/dsqrt(dnrm)
                 endif
               enddo
            endif
C*********************************************************************


C*Calculation of the mean value of kinetic and potential energy*******
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

             ekin=ekin/(lambda(izeta,l)**2)

C Potential energy after compression

             nrcomp=nint(dlog(rco(izeta,l)/b+1.0d0)/a)+1
             epot=0.0d0
             epot2=0.0d0
             do ir=1,nrcomp
                r=rofi(ir)
                r2=r/lambda(izeta,l)
                nr=nint(dlog(r2/b+1.0d0)/a)+1
                nmin=max(1,nr-npoint)
                nmax=min(nrc,nr+npoint)
                nn=nmax-nmin+1
                call ratint(rofi(nmin),g(nmin),nn,r2,rh,dy)
                rh=rh/dsqrt(lambda(izeta,l)**(2*l+3))
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
     .          'lambda =',lambda(izeta,l),
     .          'rc =',rco(izeta,l),
     .          'energy =',eorb 

                 ePAO(l)=eorb

            elseif(izeta.gt.1) then  

           write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')
     .         'izeta=',izeta,'gaussian exponent=',gexp,
     .         'rc=',rco(izeta,l),'energy=',eorb

            endif 

          write(6,'(a25,f12.6)') 'kinetic =',ekin
          write(6,'(a25,f12.6)') 'potential(screened) =',epot
          write(6,'(a25,f12.6)') 'potential(ionic) =',epot2 

            norb=norb+(2*l+1)
            indx=indx+1
            call comBasis(is,a,b,rofi,g,l,
     .              rco(izeta,l),lambda(izeta,l),izeta,
     .              nrc,norb)



              enddo 
               

            call compress_PAO(a,b,rofi,rphi(1,l),
     .              rco(1,l),lambda(1,l))


              endif    
             enddo 

 
            return
      
            end




                





               subroutine compress_PAO(a,b,rofi,rphi,
     .            rc,lambda)

C**********************************************************
C   Compression of a PAO orbital according to the compression
C   factor lambda. Input and outputare stored in the same array, 
C   and the radial grid is identical in input and output. 
C  Written by D. Sanchez-Portal, Aug. 1998
C**********************************************************
               implicit none
               include 'atom.h'


               double precision rc, lambda, rphi(nrmax), a, b,
     .            rofi(nrmax)  

C*********Internal variables*********************************************
C
             integer nr, nmax, nmin, nn, maxpoint, nrc, ir  
        
             double precision  r, phi, dy, rmax
             
             double precision
     .          aux(nrmax)

****NUMBER OF POINTS USED BY RATINT FOR THE INTERPOLATION*************
C
          integer npoint
          parameter(npoint=4)
C
C***********************************************************************



C*INTERPOLATION TO CALCULATE THE VALUE OF THE FIRST-SHELL PAO *********
C***BASIS FUNCTIONS IN THE LOGARITHMIC MESH AFTER COMPRESSION *********
C
              nrc=nint(dlog(rc/b+1.0d0)/a)+1
              rmax=rc/lambda
              maxpoint=nint(dlog(rmax/b+1.0d0)/a)+1
              
              do ir=2,nrc
                r=rofi(ir)/lambda
                nr=nint(dlog(r/b+1.0d0)/a)+1
                nmin=max(1,nr-npoint)
                nmax=min(maxpoint,nr+npoint)
                nn=nmax-nmin+1
                call ratint(rofi(nmin),rphi(nmin),nn,r,phi,dy) 
                aux(ir)=phi/dsqrt(lambda)
              enddo 
              rphi(1)=0.0d0 
              do ir=2,nrc
               rphi(ir)=aux(ir)
              enddo 
              do ir=nrc+1,nrmax
                rphi(ir)=0.0d0 
              enddo 

C
C*********************************************************************



             return

              end
 

       subroutine atm_pop(is,iz,q,qPAO,lmxo,
     .       nzeta,semic,lsemic,polorb) 


C*****************************************************************
C Returns the ground states atomic population for each species.
C This information is required for the screening of the local 
C pseudopotential.
C Written by D. Sanchez-Portal, Aug. 1998
C*****************************************************************

         implicit none
         include 'atom.h'

         integer maxos
         parameter (maxos=2*nzetmx*lmx2)

         double precision  q(maxos), qPAO(0:lmaxd) 

         integer 
     .     nzeta(0:lmaxd),polorb(lmaxd), lsemic,iz,
     .     lmxo ,is 

         logical semic


C********Internal variables************************************

        double precision  qatm(0:3)
          
        integer noPAO, l, izeta, m, norb, noPol, iorb, lpop,
     .     config(0:4)         
            
C*************Common block**************************************

        double precision qtb(maxos,nsmax), qltb(0:3,nsmax)
 
        integer cnfigtb(0:3,nsmax)

        common/cmq/qtb
        common/cmcnfig/cnfigtb
        common/cmql/qltb         

C***************************************************************
        do l=0,3
           qatm(l)=0.0d0
        enddo

        call qvlofz(iz,qatm)  

        do l=0,lmxo 
          qPAO(l)=0.0d0
          if(l.le.3) qPAO(l)=qatm(l) 
        enddo 

        if(semic) then 
          qPAO(lsemic)=qPAO(lsemic)+2*(2*lsemic+1)  
          if(lsemic.le.3) qatm(lsemic)=qatm(lsemic)+2*(2*lsemic+1)
        endif 

        noPAO=0
        do l=0,lmxo 
            if(nzeta(l).gt.0) then 
            do  izeta=1,nzeta(l) 
                  do m=1,2*l+1 
                   q(noPAO+m)=0.0d0
                   if(izeta.eq.1) q(noPAO+m)=qPAO(l)/(2*l+1)
                  enddo  
                  noPAO=noPAO+2*l+1
            enddo
            endif 
        enddo 

        noPol=0
        do l=1,min(lmxo+1, lmaxd)
           if(polorb(l).gt.0) then 
           do  izeta=1,polorb(l)
              do m=1,2*l+1 
                q(noPAO+noPol+m)=0.0d0
              enddo  
              noPol=noPol+2*l+1
           enddo    
           endif
        enddo      
        norb=noPAO+noPol 
        
        do iorb=1,norb 
           qtb(iorb,is)=q(iorb)
        enddo 

            write(6,'(/,2a)') 'atm_pop: Valence configuration',
     .                      '(local Pseudopot. screening):' 
           call cnfig(iz,config)
           lpop=min(3,lmxo)  
           if((semic).and.(lsemic.le.3)) config(lsemic)=config(lsemic)-1
           do l=0,lpop 
               if(l.eq.0)
     .    write(6,'(7x,i2,a,f5.2,a)') config(0),'s(',qatm(0),')'
               if(l.eq.1)
     .    write(6,'(7x,i2,a,f5.2,a)') config(1),'p(',qatm(1),')'
               if(l.eq.2)
     .    write(6,'(7x,i2,a,f5.2,a)') config(2),'d(',qatm(2),')'
               if(l.eq.3)
     .    write(6,'(7x,i2,a,f5.2,a)') config(3),'f(',qatm(3),')'
           enddo

          do l=0,3
              cnfigtb(l,is)=config(l)
              qltb(l,is)=qatm(l)
          enddo
 
 
       return 
       end








        subroutine Vna(is,Zval,qPAO,rphi,rco,vloc,
     .        a,b,rofi,drdi,nrval,lmxo,nVna) 

C******************************************************************
C  Generates the neutral-atom pseudopotential.
C  D. Sanchez-Portal, Aug. 1998.
C******************************************************************
          
          implicit none
          include 'atom.h'
 
          double precision
     .    Zval, qPAO(0:lmaxd), rofi(nrmax), rphi(nrmax,0:lmaxd),
     .    drdi(nrmax),rco(nzetmx,0:lmaxd),a,b, vloc(nrmax) 
      
          integer
     .    nrval, lmxo,is, nVna


C********Internal variables***************************************

         double precision
     .    rho(nrmax), chval, ve(nrmax), s(nrmax),eps, phi,
     .    rcocc, dincv, rVna
         
         integer
     .    nrc,ir, l, ncocc 
 



          do ir=1,nrval
             rho(ir)=0.0d0
             s(ir)=dsqrt(drdi(ir))
          enddo

          chval=0.0d0 
          ncocc=0 
          do l=0,lmxo
             if(qPAO(l).gt.0.0d0) then
              nrc=nint(dlog(rco(1,l)/b+1.0d0)/a)+1 
              ncocc=max(ncocc,nrc)
              do ir=2,nrc 
                phi=rphi(ir,l)
                rho(ir)=rho(ir)+qPAO(l)*phi**2 
                chval=chval+drdi(ir)*qPAO(l)*phi**2
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
C
          call vhrtre(rho,ve,rofi,drdi,s,nrval,a)
C
C*********************************************************************

 
C*********LOCAL NEUTRAL-ATOM PSEUDOPOTENTIAL**************************
          eps=1.0d-5
          nVna=0
          do ir=nrval,2,-1
               dincv=vloc(ir)+ve(ir)
               if((dabs(dincv).gt.eps).and.(nVna.eq.0)) nVna=ir+1
               ve(ir)=dincv
          enddo 
          nVna=max(nVna,ncocc)

C*********CUT-OFF RADIUS FOR THE LOCAL NEUTRAL-ATOM PSEUDOPOTENTIAL***** 

          rcocc=b*(dexp(a*(ncocc-1))-1.0d0)

          if(nVna.eq.ncocc) then
            rVna=rcocc
          else
            rVna=b*(dexp(a*(nVna-1))-1.0d0)
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
     .        'Check ATOM: Look for the sentence:'
           write(6,"(2a)")'Vna: WARNING: ',
     .        'LOCAL NEUTRAL-ATOM PSEUDOPOTENTIAL'
           write(6,"(2a)")'Vna: WARNING: ',
     .        'Increasing the tolerance parameter EPS'
           write(6,"(2a)")'Vna: WARNING: ',
     .        'might be a good idea'
         endif

          ve(1)= ( ve(2)*rofi(3)**2 - ve(3)*rofi(2)**2 ) /
     .          (      rofi(3)**2 -      rofi(2)**2 )
         

C*******Construct the common block with the neutral-atom potential**** 
C
          call comVna(is,a,b,rofi,Ve,nVna,1.0d0)
C
C**********************************************************************

           return
           end        



           subroutine comVna(is,a,b,rofi,Vna,nVna,flting)

C***********************************************************************
C  Creates the common block with the information about the neutral atom
C  pseudoptential.
C   D. Sanchez-Portal, Aug. 1998.
C***********************************************************************

           implicit none
           include 'atom.h'

           integer nVna, is
 
           double precision 
     .        Vna(nrmax), rofi(nrmax), a ,b, flting

          
C**********Internal variables*********************************************
C
           integer nr, nmin, nmax, nn, itb

           double precision 
     .     yp1, ypn, dy, v, rVna, delt, r
  
           double precision 
     .          aux(ntbmax)

           double precision deltmax
           common/cmdelt/deltmax
C
C***********************************************************************

C****NUMBER OF POINTS USED BY RATINT FOR THE INTERPOLATION*************
C
          integer npoint
          parameter(npoint=4)
C
C***********************************************************************

C***********Variables in common blocks**********************************
C
          double precision
     .  table((ntbmax+2),-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tab2(ntbmax,-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tabpol((ntbmax+2),nzetmx*lmaxd,nsmax),
     .  tab2pol(ntbmax,nzetmx*lmaxd,nsmax)



              common/cmtab/table,tabpol
              common/cmspline/tab2,tab2pol
  
C
C***********************************************************************



       if (flting.gt.0.0d0) then 
         rVna=b*(dexp(a*(nVna-1))-1.0d0)
         delt=rVna/(dble(ntbmax-1)+1.0d-20) 
 
         if(delt.gt.deltmax) then
              write(6,'(a)')
     .    'comVna: WARNING It might be a good idea to increase'
              write(6,'(a)')
     .    'comVna: WARNING parameter ntbmax (in file atom.h) '
              write(6,'(a,i6)')
     .    'comVna: WARNING to at least ntbmax = ',
     .        nint(rVna/deltmax)+2
          endif

         table(1,0,is)=delt
         table(2,0,is)=rVna
       elseif(flting.lt.0.0d0) then
          table(1,0,is)=0.0d0
          table(2,0,is)=0.0d0 
       endif

           do itb=1,ntbmax-1
               r=delt*(itb-1)
               nr=nint(dlog(r/b+1.0d0)/a)+1
               nmin=max(1,nr-npoint)
               nmax=min(nVna,nr+npoint)
               nn=nmax-nmin+1 
               if(flting.gt.0.0d0) then
                 call ratint(rofi(nmin),Vna(nmin),nn,r,v,dy) 
               else
                 v=0.0d0  
                 tab2(itb,0,is)=0.0d0
               endif 
               table(itb+2,0,is)=v 
          enddo 
          table(ntbmax+2,0,is)=0.0d0

C*********TABLE WITH THE SECOND DERIVATIVE ****************************
C 
          if (flting.gt.0.0d0) then
            yp1=0.d0
            ypn=1.0d50

            call spline(delt,table(3,0,is),ntbmax,
     .        yp1,ypn,tab2(1,0,is),aux)
          endif
C
C*********************************************************************

           return 
 
           end     







       subroutine slfe_local(slfe,vlocal,rofi,a,nVna,drdi)
C*********************************************************************
C Calculates the self-energy associated to the local-pseudopotential
C charge density.
C Written by D. Sanchez-Portal, Aug. 1998.
C*********************************************************************


         implicit none
         include 'atom.h'
           
         double precision slfe, vlocal(nrmax),rofi(nrmax),a,
     .       drdi(nrmax) 
 

         integer nVna 
           
C***********Internal variables********************************
          
         double precision slf, a2b4, s(nrmax), g0, g1, g2,
     .       g3, g4, d2g, d2u

         integer ir

          do ir=1,nVna
             s(ir)=dsqrt(drdi(ir))
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
   

           return 
           end


       
           subroutine POLgen(is,a,b, rofi, drdi,
     .          ePAO,rphi,rco,vps,ve,
     .          polorb,lmxo, norb)

C*********************************************************************
C Calculates the polarization  orbitals for the basis set augmentation.
C Written by D. Sanchez-Portal, Aug. 1998.
C*********************************************************************
               implicit none
               include 'atom.h'

               double precision
     .         a, b, rofi(nrmax), vps(nrmax,0:lmaxd),
     .         ve(nrmax), drdi(nrmax),
     .         rphi(nrmax,0:lmaxd), rco(nzetmx,0:lmaxd),
     .         ePAO(0:lmaxd)

               integer
     .           lmxo, is, norb, polorb(lmaxd)


C********Internal variables*************************************

               integer
     .           l, nrc, nsp, ir,indx,
     .           ipol

               double precision
     .           rc, rcpol(nzetmx,lmaxd),
     .           splnorm_default, phipol(nrmax),
     .           rnrm(nrmax), dnrm, phi,
     .           cons1, cons2, spln, 
     .           splnorm, g(nrmax), r, ekin, 
     .           r1, r2, dfdi, d2fdi2, d2fdr2, dr,
     .           epot, epot2, eorb, eps




C***********Declarations for FDF package *******************************
C
         include 'fdf/fdfdefs.h'
C
C***********************************************************************

C******Common blocks with the values of some parameters*********
C
        common/cmspldefault/splnorm_default
C
C***************************************************************       

      

C***READING SPLNORM TO GENERATE THE SPLIT IF Rmatch IS ZERO IN INPUT****

         splnorm=fdf_double('PAO.SplitNorm',splnorm_default)

C***********************************************************************





 
             norb=0 
             indx=0
             do l=1,lmxo+1


               if(polorb(l).gt.0) then 

            write(6,'(/,(a,i2))')
     .    'POLgen: Perturbative polarization orbital with L= ',l

              do ipol=1,polorb(l)

                  
                if(ipol.eq.1) then  
                  rc=rco(1,l-1) 
                  rcpol(ipol,l)=rc
                  nrc=nint(dlog(rc/b+1.0d0)/a)+1
C**Generate the polarization function perturbatively from the original PAO**
C 

            call polarization(a,rofi,rphi(1,l-1),vps(1,l-1),
     .            ve,drdi,nrc,l-1,ePAO(l-1),g,nrc)

                       dnrm=0.0d0
                       do ir=2,nrc-1
                          phi=g(ir)  
                          phipol(ir)=phi
                          dnrm=dnrm+drdi(ir)*phi*phi
                          rnrm(ir)=dnrm 
                          g(ir)=g(ir)/(rofi(ir)**(l+1))
                       enddo   
                       g(1)=g(2)
                       g(nrc)=0.0d0
                       phipol(nrc)=0.0d0
                elseif(ipol.gt.1) then  

                   rc=rco(1,l-1) 
                   nrc=nint(dlog(rc/b+1.0d0)/a)+1
C***Multiple shells can be generated using the split scheme*************
       
            spln=splnorm
            if(ipol.gt.2) then
              spln=spln/(2.0d0*(ipol-2) )
            endif

            call parabola(a,b,nrc,phipol,rnrm,
     .                   l,spln,cons1,cons2,nsp)


C***Cut-off radius for the split orbital with a desired norm****** 
         nrc=nsp
         rcpol(ipol,l)=b*(dexp(a*(nrc-1))-1.0d0)
C*****************************************************************

              
            dnrm=0.0d0
            do ir=2,nrc-1
              r=rofi(ir)
C***********************parabolic split****************************
              phi=-(cons1*r**2+cons2)+phipol(ir)/(r**(l+1))
C******************************************************************

              dnrm=dnrm+drdi(ir)*(phi*r**(l+1))**2
              g(ir)=phi
            enddo
            g(1)=g(2)
            g(nrc)=0.0d0
               
            
              endif

C**************Normalization of basis functions***********************
            eps=1.0d-4
            if(dabs(dnrm-1.0d0).gt.eps) then
               do ir=1,nrc
                 g(ir)=g(ir)/dsqrt(dnrm)
                 if(ipol.eq.1) then
                    phipol(ir)=phipol(ir)/dsqrt(dnrm)
                    rnrm(ir)=rnrm(ir)/dnrm
                 endif
               enddo
            endif
C*********************************************************************


C*Calculation of the mean value of kinetic and potential energy*******
C    Potential and kinetic energy of the orbital



             ekin=0.0d0 
             epot=0.0d0
             epot2=0.0d0
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

                epot=epot+ 
     .          drdi(ir)*(ve(ir)+vps(ir,l))*(g(ir)*r**(l+1))**2 
                epot2=epot2+
     .          drdi(ir)*vps(ir,l)*(g(ir)*r**(l+1))**2

             enddo
             eorb=ekin+epot
       



            if(ipol.eq.1) then  

             write(6,'(/,(3x,a,i2),2(/,a25,f12.6))')
     .          'izeta =',ipol,
     .          'rc =',rcpol(ipol,l),
     .          'energy =',eorb 

            elseif(ipol.gt.1) then 

            write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')
     .         'izeta =',ipol,
     .         'rmatch =',rcpol(ipol,l),
     .         'splitnorm =',spln,
     .         'energy =',eorb 

            endif 

          write(6,'(a25,f12.6)') 'kinetic =',ekin
          write(6,'(a25,f12.6)') 'potential(screened) =',epot
          write(6,'(a25,f12.6)') 'potential(ionic) =',epot2 

            norb=norb+(2*l+1)
            indx=indx+1
            call comPOL(is,a,b,rofi,g,l,
     .              rcpol(ipol,l),ipol,nrc,indx)




              enddo 
        




              endif    
             enddo 

 
            return
      
            end




                





               subroutine comPOL(is,a,b,rofi,rphi,
     .            l,rc,ipol,nrc,norb)

C*********************************************************************
C Generates the common block with the information about the polarization
C orbitals.
C  Written by D. Sanchez-Portal, Aug. 1998.
C*********************************************************************

               implicit none
               include 'atom.h'

               integer l, is, norb, nrc, ipol

               double precision rc, rphi(nrmax), a, b,
     .            rofi(nrmax)  

C*********Internal variables*********************************************
C
             integer itb, nr, nmax, nmin, nn  
        
             double precision delt, r, phi, dy, yp1, ypn
             
             double precision
     .          aux(ntbmax)

             double precision deltmax 
             common/cmdelt/deltmax

****NUMBER OF POINTS USED BY RATINT FOR THE INTERPOLATION*************
C
          integer npoint
          parameter(npoint=4)
C
C***********************************************************************

C***********Variables in common blocks**********************************
C
          double precision
     .  table((ntbmax+2),-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tab2(ntbmax,-(lmaxd+1):nzetmx*(lmaxd+1),nsmax),
     .  tabpol((ntbmax+2),nzetmx*lmaxd,nsmax),
     .  tab2pol(ntbmax,nzetmx*lmaxd,nsmax),
     .  rcpoltb(nzetmx,lmaxd,nsmax)
     

C***********************************************************************
C
              common/cmtab/table,tabpol
              common/cmspline/tab2,tab2pol 
              common/cmradpol/rcpoltb
C
C***********************************************************************

            rcpoltb(ipol,l,is)=rc


C**********INTERPOLATION TO GENERATE TABLES WITH KB PROJECTORS*******
C
            delt=rc/(dble(ntbmax-1)+1.0d-20) 

          if(delt.gt.deltmax) then
              write(6,'(a)')
     .    'comPOL: WARNING It might be a good idea to increase'
              write(6,'(a)')
     .    'comPOL: WARNING parameter ntbmax (in file atom.h) '
              write(6,'(a,i6)')
     .    'comPOL: WARNING to at least ntbmax = ',
     .        nint(Rc/deltmax)+2
           endif

            tabpol(1,norb,is)=delt

            tabpol(2,norb,is)=dble(l)


            do itb=1,ntbmax-1
                 r=delt*(itb-1)
                 nr=nint(dlog(r/b+1.0d0)/a)+1
                 nmin=max(1,nr-npoint)
                 nmax=min(nrc,nr+npoint)
                 nn=nmax-nmin+1
                 call ratint(rofi(nmin),rphi(nmin),nn,r,phi,dy)
                 tabpol(itb+2,norb,is)=phi 
            enddo
            tabpol(ntbmax+2,norb,is)=0.0d0  
                       
C
C*********************************************************************

C*********TABLE WITH THE SECOND DERIVATIVE ****************************
C

            yp1=1.0d50
            ypn=1.0d50

            call spline(delt,tabpol(3,norb,is),ntbmax,
     .        yp1,ypn,tab2pol(1,norb,is),aux)
        
C
C*********************************************************************



             return

              end



        subroutine set_mesh(a,b,rofi,drdi,s)

C**************************************************************
C    Setting up mesh points an its derivatives from standard
C    values
C    D. Sanchez-Portal, Aug. 98
C**************************************************************
 
        implicit none
        include 'atom.h'

        double precision 
     .     rofi(nrmax), drdi(nrmax), s(nrmax), a, b


C**** Internal variables******************************************
        double precision 
     .     aa, bb, zt, rpb, ea, ea2
        integer ir

        parameter(zt=1.0d0)
        parameter(aa=80.0d0)
        parameter(bb=6.0d0) 


C***********STANDART VALUES FOR MESH PARAMETERS*************************

          b=dexp(-bb)/zt
          a=1.0d0/aa

C***********SET UP THE MESH POINTS AND ITS DERIVATIVE******************

          rpb=b
          ea=dexp(a)
          ea2=1.0d0
          do ir=1,nrmax
            drdi(ir)=a*rpb
            rofi(ir)=b*(ea2-1.0d0)
            s(ir)=(a*rpb)**2
            rpb=rpb*ea
            ea2=ea2*ea
          enddo



          return
          end




           subroutine BESSEL(is,a,b,rofi,drdi,s,
     .             lmxo,
     .             nzeta,rco,lambda, norb) 

C*********************************************************************
C  Caculates Bessel functions as a floating basis
C  Written by D. Sanchez-Portal, Aug. 1998.
C*********************************************************************
               implicit none
               include 'atom.h'

               double precision
     .         a, b, rofi(nrmax),
     .         drdi(nrmax), s(nrmax), 
     .         rco(nzetmx,0:lmaxd),
     .         lambda(nzetmx, 0:lmaxd)


               integer
     .           lmxo, is, nzeta(0:lmaxd),
     .           norb


C********Internal variables*************************************

               integer
     .           l,nprin, nnodes, nodd, nrc, ir,indx,
     .           izeta

               double precision
     .           rc, v(nrmax),
     .           dnrm, phi,
     .           g(nrmax), eorb, eps

C*************************************************************
               data  v/ nrmax*0.0d0 /


 
             norb=0 
             indx=0
             do l=0,lmxo 

               if(nzeta(l).gt.0) then

            write(6,'(/2A,I2)')
     .       'Bessel: floating Bessel functions ',
     .           'with angular momentum L=',l





              do izeta=1, nzeta(l) 

C*******Cut-off radius for Bessel functions must be an explicit input***
C
                if(rco(izeta,l).lt.1.0d-5) then 
                   write(6,'(a)')
     .     'Bessel: ERROR Zero cut-off radius with Z=-100 option'
                   write(6,'(a)')
     .     'Bessel: ERROR Cut-off radius must be explicitely specified'
                   write(6,'(a)')
     .     'Bessel: ERROR using Z=-100 (Floating Bessel functions)'
                    stop
 
                endif
C
C********************************************************************


C***********************************************************************
           if(dabs(lambda(izeta,l)-1.0d0).gt.1.0d-3) then
               write(6,'(/,a)')
     . 'Bessel: WARNING Scale factor is not active with Z=-100 option' 
           endif 
           lambda(izeta,l)=1.0d0
C***********************************************************************

                  rc=rco(izeta,l)
                  nrc=nint(dlog(rc/b+1.0d0)/a)+1
                  nodd=mod(nrc,2)
                  if(nodd.eq.0) then
                     nrc=nrc+1
                  endif
                  rc=b*(dexp(a*(nrc-1))-1.0d0)
                  rco(izeta,l)=rc
                  

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

C*********Checking normalization of the wavefunctions*****************
                 eps=1.0d-4
                 if(dabs(dnrm-1.0d0).gt.eps) then
                   do ir=1,nrc
                   g(ir)=g(ir)/dsqrt(dnrm)
                  enddo
                 endif
C*********************************************************************  



             write(6,'(/,(3x,a,i2),2(/,a25,f12.6))')
     .          'izeta =',izeta,
     .          'rc =',rco(izeta,l),
     .          'energy =',eorb  

            norb=norb+(2*l+1)
            indx=indx+1
            call comBasis(is,a,b,rofi,g,l,
     .              rco(izeta,l),lambda(izeta,l),izeta,
     .              nrc,indx)


              enddo 
        
              endif  

             enddo 

 
            return
      
            end




                






         subroutine draw_basis(ntotsp)

C*******************************************************************
C This routine prints a file with information about the basis set. 
C The format used is compatible with that required by USER option       
C for basis generation.
C 
C  INPUT: 
C       ntotsp    : Total number of different chemical species
C Written by D. Sanchez-Portal, Sept. 98
C******************************************************************

      implicit none
      include 'atom.h'
 
      integer  ntotsp  


C**********************External declarations***********************
C   This routine, which is called in the last call of atom, uses 
C   the interface functions provided in the file atom_functons.f
C   This functions has been created to provide all the information 
C   about the chemical species and the basis sets, and designed to be 
C   called out of routine 'atom'. Here, they can be used 
C   because all the common blocks have been initialized previously.
C
        integer
     .   nztfl, lomaxfis, nofis, lofio 
        double precision
     .   rcut, phiatm
        character 
     .   labelfis*20, paste*80

        external 
     .   nztfl, lomaxfis, rcut, phiatm, labelfis, nofis, lofio,
     .   paste 
C
C*****************************************************************


C*********Internal variables**************************************
C

C
C      Number of points in the homogeneus radial mesh 
C      for plotting purposes
C
       integer ndraw
       parameter(ndraw=200) 

       integer l, ir, is, lmax, nztmx, izt, io, iomax, izeta

       double precision r, phi, dphidr, delt, sum, rc
       
        character filename*80
C
C******************************************************************
C

             do  is=1,ntotsp 

                filename=paste(labelfis(is),'.PAO.basis')
                open(unit=1, file=filename, status='unknown',
     .                      form='formatted')        
        
                lmax=lomaxfis(is) 
                write(1,*) lmax  

                do l=0,lmax 
 
                   nztmx=nztfl(is,l) 

                   write(1,*) l, nztmx 
                   
                   do izt=1, nztmx

C**********************************************************************
C The only problem now is to translate the l,izt index into the orbital 
C index io (index inside each atom). The problem is that index for
C PAOs and polarization orbitals of each symmetry are not consecutive. 
C
                         iomax=nofis(is)
                         izeta=0
                         do io=1,iomax
                            if(l.eq.lofio(is,io)) then 
                              izeta=izeta+1
                              if(izeta.gt.(2*l+1)*(izt-1)) goto 10 
                            endif 
                         enddo 
                                 
 10                      continue 
C 
C   Now we have the needed orbital index, io   
C*********************************************************************

                         rc=rcut(is,io)
                         write(1,*) izt, ndraw, rc      


                         delt=rc/(ndraw-1)
                         sum=0.0d0
                         do ir=1, ndraw

                          r=delt*(ir-1)
                          call rphiatm(is,io,r,phi, dphidr)  
                          sum=sum+(r*phi)**2*delt 
  
c                         write(1,*) r, r*phi, sum 
                          write(1,*) r, r*phi


 
                         enddo  


                      enddo  

                  enddo 

               enddo  



               return
               end



           subroutine USER(is,a,b,rofi,drdi,
     .             vps,ve,
     .             lmxo,
     .             nzeta,rco,lambda, rphi, ePAO, norb) 

C*****************************************************************
C Read the basis orbitals provided by the user and calculates its
C energy.
C Written by D. Sanchez-Portal, Aug. 1998
C*****************************************************************
               implicit none
               include 'atom.h'

               double precision
     .         a, b, rofi(nrmax), vps(nrmax,0:lmaxd),
     .         drdi(nrmax), ve(nrmax),
     .         rphi(nrmax,0:lmaxd), rco(nzetmx,0:lmaxd),
     .         lambda(nzetmx, 0:lmaxd),
     .         ePAO(0:lmaxd)


               integer
     .           lmxo, is, nzeta(0:lmaxd),
     .           norb


C********Internal variables*************************************

               integer
     .           l,nodd, nrc, ir,indx,nrcfile, jr, nrcomp,
     .           izeta, nmax, nmin, nn, nr,
     .           lmxfile, luser, nzetfile, npfile, izetafile

               double precision
     .           rc,rcfile,
     .           dnrm, phi,
     .           g(nrmax), r, ekin, rfile(nrmax),
     .           r1, r2, dfdi, d2fdi2, d2fdr2, dr,
     .           epot, epot2, rh, dy, eorb, eps,
     .           swap(nrmax), rf(nrmax) 
 
               character filename*80, paste*80 
               external paste


C******Common block with the species labels*****************************
                 character  label_save(nsmax)*20
                 common/cmlabel/label_save


C****NUMBER OF POINTS USED BY RATINT FOR THE INTERPOLATION*************
C 
               integer  npoint 
               parameter(npoint=4)
C
C***********************************************************************


             filename=paste(label_save(is),'.user.basis')

             open(unit=1, file=filename, status='old',
     .                                  form='formatted')
             
            write(6,'(/,2a,/2a)')
     .     'USER: User-basis. The basis orbitals will be read ',
     .     'from the file: ','USER: ',filename 
             
            read(1,*) lmxfile 

            if(lmxfile.lt.lmxo) then
              write(6,'(/a)')
     .       'USER: ERROR maximum angular momentum of the orbitals'
              write(6,'(2a)')
     .       'USER: ERROR in file ', filename 
              write(6,'(a,i2,a,i2)')
     .       'USER: ERROR is ',lmxfile,
     .        'while it should be at least', lmxo
              stop
            endif


             norb=0 
             indx=0
             do l=0,lmxo 
               if(nzeta(l).gt.0) then
                 read(1,*) luser, nzetfile 


            if(luser.ne.l) then
               write(6,'(/,2a, /2(a,i2),/a)')
     .      'USER: ERROR Reading user-basis orbitals from file:',
     .       filename,
     .      'USER: ERROR expected l=',l,' read l=',luser,
     .      'USER: ERROR check order'
               stop
            endif

             if (nzetfile.lt.nzeta(l)) then
                  write(6,'(/,2a,i2,/2a,/2(a,i2))')
     .           'USER: ERROR ',
     .           'number of orbitals with l=', l,'USER: ERROR in file',
     .            filename, 'USER: ERROR is ',nzetfile,
     .            'it should be at least', nzeta(l)
                stop
             endif

            write(6,'(/A,I2)')
     .       'USER: Orbitals with angular momentum L=',l


              do izeta=1,nzeta(l)

C*****IF THE COMPRESSION FACTOR IS NEGATIVE OR ZERO THE ORBITALS ARE***
C**************LEFT UNTOUCHED******************************************
             if(lambda(izeta,l).le.0.0d0) lambda(izeta,l)=1.0d0
C*********************************************************************** 

                  read(1,*) izetafile,npfile,rcfile

                  if(izetafile.ne.izeta) then
                    write(6,'(/,2a,/2(a,i2),/a)')
     .      'USER: ERROR Reading user-basis orbitals from file:',
     .       filename,
     .      'USER: ERROR expected zeta=',izeta,
     .      'USER: ERROR read zeta=',izetafile,
     .      'USER: ERROR check order'
             stop
                  endif
                  if ((npfile+1).gt.nrmax) then
                     write(6,'(/,2a,/2a,/2a,i6)')
     .       'USER: ERROR ',
     .       'Too many grid points required to read functions', 
     .       'USER: ERROR in file: ',filename,
     .       'USER: ERROR Parameter nrmax in atom.h must be increased',
     .        'to at least',npfile+1
                     stop
                  endif


              if(rco(izeta,l).lt.1.0d-5) 
     .            rco(izeta,l)=rcfile*lambda(izeta,l)
              if(rco(izeta,l).lt.1.0d-5) then  
                write(6,'(/,a,/,2a,/,a)')
     .      'USER: ERROR A non-zero cut-off radius should be given',
     .      'USER: ERROR either in the input file or in the file',
     .       filename, 'USER: ERROR while using USER basistype option' 
              stop
              endif
           
                    


                  rc=rco(izeta,l)/lambda(izeta,l) 
                  do ir=1,npfile
                     read(1,*) r, rf(ir)
                     rfile(ir)=r
                     if((r.gt.rc).or.
     .                   (dabs(rc-r).lt.1.0d-5)) goto 10
                  enddo 
                  rc=r
                     write(6,'(/2a,/2a)')
     .     'USER: WARNING Basis orbitals read from file ',filename,
     .     'USER: WARNING The required Rc is larger ',
     .     'than the maximum radial'
          write(6,'((a,f12.6),(/2a,f12.6))')
     .     'USER: WARNING grid point specified in the file ',r,
     .     'USER: WARNING New cut-off radius for this orbital ',
     .      '(after compression) ', r*lambda(izeta,l)  

10              nrcfile=ir
                do ir=nrcfile+1,npfile
                      read(1,*) rfile(ir),rf(ir)
                enddo
           
                nrc=nint(dlog(rc/b+1.0d0)/a)+1
                nodd=mod(nrc,2)
                if(nodd.eq.0) then
                    nrc=nrc+1
                endif
                rc=b*(dexp(a*(nrc-1))-1.0d0)
                rco(izeta,l)=rc*lambda(izeta,l)


C*********First point should be zero***********************************
C
              if(rfile(1).ne.0.0d0) then  
                do ir=1,npfile
                   swap(ir)=rfile(ir) 
                   g(ir)=rf(ir)  
                enddo
                do ir=1,npfile
                   rfile(ir+1)=swap(ir)
                   rf(ir+1)=g(ir)
                enddo
                rfile(1)=0.0d0
                g(1)=0.0d0
                npfile=npfile+1
                nrcfile=nrcfile+1
              endif
C
C********************************************************************** 


C****Interpolation in the logaritmic mesh where pseudopotentials*******
C****are defined******************************************************* 
C
              dnrm=0.0d0
              nr=1    
              do ir=2,nrc
                 r=rofi(ir)
                 do jr=nr,npfile
                    if(rfile(jr).ge.rofi(ir)) goto 20
                 enddo
20               nr=jr
                 nmin=max(1,nr-npoint)
                 nmax=min(npfile,nr+npoint) 
                 nn=nmax-nmin+1
                 call ratint(rfile(nmin),rf(nmin),nn,r,phi,dy) 
c                call polint(rfile(nmin),rf(nmin),nn,r,phi,dy)
                 dnrm=dnrm+drdi(ir)*(phi**2)
                 if(izeta.eq.1) rphi(ir,l)=phi   
                 g(ir)=phi/(r**(l+1))  
              enddo 
              g(1)=g(2) 
C
C*********************************************************************

    
C**************Normalization of basis functions***********************
            eps=1.0d-4
            if(dabs(dnrm-1.0d0).gt.eps) then
               do ir=1,nrc
                 g(ir)=g(ir)/dsqrt(dnrm)
                 if(izeta.eq.1) then
                    rphi(ir,l)=rphi(ir,l)/dsqrt(dnrm)
                 endif
               enddo
            endif
C*********************************************************************  

C*Calculation of the mean value of kinetic and potential energy*******
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

             ekin=ekin/(lambda(izeta,l)**2)

C Potential energy after compression

             nrcomp=nint(dlog(rco(izeta,l)/b+1.0d0)/a)+1
             epot=0.0d0
             epot2=0.0d0
             do ir=1,nrcomp
                r=rofi(ir)
                r2=r/lambda(izeta,l)
                nr=nint(dlog(r2/b+1.0d0)/a)+1
                nmin=max(1,nr-npoint)
                nmax=min(nrc,nr+npoint)
                nn=nmax-nmin+1
                call ratint(rofi(nmin),g(nmin),nn,r2,rh,dy)
                rh=rh/dsqrt(lambda(izeta,l)**(2*l+3))
                epot=epot+
     .          drdi(ir)*(ve(ir)+vps(ir,l))*(rh*r**(l+1))**2
                epot2=epot2+
     .          drdi(ir)*vps(ir,l)*(rh*r**(l+1))**2
             enddo
             eorb=ekin+epot



             write(6,'(/,(3x,a,i2),3(/,a25,f12.6))')
     .          'izeta =',izeta,
     .          'lambda =',lambda(izeta,l),
     .          'rc =',rco(izeta,l),
     .          'energy =',eorb  

              if(izeta.eq.1) ePAO(l)=eorb


          write(6,'(a25,f12.6)') 'kinetic =',ekin
          write(6,'(a25,f12.6)') 'potential(screened) =',epot
          write(6,'(a25,f12.6)') 'potential(ionic) =',epot2 

            norb=norb+(2*l+1)
            indx=indx+1
            call comBasis(is,a,b,rofi,g,l,
     .              rco(izeta,l),lambda(izeta,l),izeta,
     .              nrc,indx)


              enddo 
              do izeta=nzeta(l)+1,nzetfile
                  read(1,*) izetafile,npfile
                  do ir=1,npfile
                    read(1,*)
                  enddo
             enddo
 
            call compress_PAO(a,b,rofi,rphi(1,l),
     .              rco(1,l),lambda(1,l))  
                       

             endif     
             enddo 

            return
            end

                
            subroutine prinput(ntotsp)

C*******************************************************************
C Prints the values of the parameter which have been actually 
C used in the generation of the basis (cut-off radius, contraction
C factors, and BasisType option to augment the basis set).
C The information is written in the same format as required for the 
C input file.
C Written by D. Sanchez-Portal, Oct. 1998.
C*******************************************************************

            implicit none 
            include 'atom.h'

            integer ntotsp

C******************Internal variables***********************************
           integer is, izofis, nshell, l, lo, nzt, izt
          
           character labelfis*20, basistype*10,
     .       basistype_default*10, rcchar(nzetmx)*11,
     .       lambdachar(nzetmx)*11


C**********************External declarations***********************
C   This routine, which is called in the last call of atom, uses
C   the interface functions provided in the file atom_functions.f
C   This functions has been created to provide all the information
C   about the chemical species and the basis sets, and designed to be
C   called out of routine 'atom'.  Here, they can be used
C   because all the common blocks have been initialized previously.
C
           external 
     .      izofis, labelfis
C
C**********************************************************************
           include 'fdf/fdfdefs.h'
           parameter(basistype_default='split')
C**********************************************************************


C********************Variables in common blocks************************ 
           integer
     .         lmxosave(nsmax),
     .         npolorbsave(lmaxd,nsmax),
     .         nzetasave(0:lmaxd,nsmax)

           double precision
     .          rcotb(nzetmx,0:lmaxd,nsmax),
     .          lambdatb(nzetmx,0:lmaxd,nsmax),
     .          chargesave(nsmax)

           character*10  
     .          basistype_save(nsmax)

              common/cmradorb/rcotb
              common/cmlambda/lambdatb 
              common/cmlmxo/lmxosave
              common/cmzeta/nzetasave
              common/cmpolorb/npolorbsave
              common/cmcharge/chargesave 
              common/cmbasistype/basistype_save
C**********************************************************************



             write(6,'(/2a)')
     .            'prinput: ************************* Basis input ',
     .            '********************************'

             basistype=fdf_string('PAO.BasisType',basistype_default)
             call type_name(basistype) 

             write(6,'(/2a)')'PAO.BasisType ',basistype

             write(6,'(/a)')
     .                   '%block ChemicalSpeciesLabel' 
             do is=1,ntotsp
                write(6,'(2(1x,i3),1x,2a)')
     .                  is,izofis(is),labelfis(is), 
     .              '    # Species index, atomic number, species label' 
             enddo 
             write(6,'(a)')
     .                  '%endblock ChemicalSpeciesLabel' 

             write(6,'(/a)')        
     .   '%block PAO.Basis                 # Define Basis set'
             do is=1, ntotsp  
                
                nshell=0 
                lo=lmxosave(is)
                do l=0,lo
                    if(nzetasave(l,is).ne.0) nshell=nshell+1
                enddo 
                if(basistype_save(is).eq.basistype) then  

                if(dabs(chargesave(is)).lt.1.0d-4) then 
                     write(6,'(a10,1x,i2,20x,a)')
     .                labelfis(is), nshell,
     .                 '# Species label, number of l-shells'
                else
                    write(6,'(a10,1x,i2,1x,f7.3,12x,2a)')
     .                labelfis(is), nshell, chargesave(is),
     .                  '# Label, l-shells,',
     .              ' ionic net charge'
                endif 
         
               else 

               if(dabs(chargesave(is)).lt.1.0d-4) then
                     write(6,'(a10,1x,i2,1x,a,10x,2a)')
     .             labelfis(is), nshell, basistype_save(is), 
     .                  '# Species label, l-shells,',
     .                  ' basis type '
                else
                    write(6,'(a10,1x,i2,1x,a,1x,f7.3,1x,2a)')
     .              labelfis(is), nshell, basistype_save(is),
     .              chargesave(is),
     .                  '# Label, l-shells, type,',
     .              ' ionic net charge'
                endif 
           
                endif 
 
                   do l=0,lo
                      nzt=nzetasave(l,is)
                      if(nzt.ne.0) then  
                        if(((l+1).le.lmaxd).and.
     .                     (npolorbsave(l+1,is).gt.0)) then 
                           write(6,'(2(1x,i3),a,i3,19x,2a)') 
     .                            l, nzt, ' P ',npolorbsave(l+1,is),
     .                          '# l, Nzeta, ','Polarization, NzetaPol'
                        else
                           write(6,'(2(1x,i3),25x,a)')
     .                           l, nzt, '# l, Nzeta '
                        endif  
                        do izt=1, nzt
                           write(rcchar(izt),'(1x,f7.3)') 
     .                                          rcotb(izt,l,is)
                           write(lambdachar(izt),'(1x,f7.3)') 
     .                                       lambdatb(izt,l,is)
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
             write(6,'(a)')
     .                       '%endblock PAO.Basis' 


             write(6,'(/2a)')
     .             'prinput: **********************************',
     .             '************************************'

             
               return 
               end 
