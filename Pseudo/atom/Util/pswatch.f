c
      program pswatch
c
C Program pswatch to analize the existence of ghost states 
C in the pseudopotentials generated using Troullier-Martins code
C Written by Daniel Sanchez Portal.Oct. 1995. 

       implicit double precision (a-h,o-z)

       parameter(nrmax=2000,lmaxd=4)

       character
     .   NAMATM*2, ICORR*2, IREL*3, NICORE*4, NAMAUX*2,
     .   method(6)*10,text(10)*10,SYMBOL*2,PASTE*23,
     .   functl*3 ,author*4

       double precision
     .  rofi(nrmax),vps(nrmax,0:lmaxd),rphi(nrmax,0:lmaxd),
     .  ve(nrmax),chcore(nrmax),vxc(nrmax),auxrho(nrmax),
     .  drdi(nrmax),elocal(2,0:lmaxd),s(nrmax),h(nrmax),
     .  g(nrmax),y(nrmax),rho(nrmax),eigen(0:lmaxd),dkbcos(0:lmaxd)

      
   
C****************     PI     ***********************************
               pi=dacos(-1.0d0)
C***************************************************************      



        
           open(unit=1,file='VPSOUT',form='unformatted')

           read(1) namatm, icorr, irel, nicore,
     .     (method(i),i=1,6), (text(i),i=1,7),
     .     npotd, npotu, nr, b, a, zval
 
          
           if(irel.eq.'rel') then 
         write(6,*) 'Pseudopotential generated from an ',
     .   'atomic-relativistic calculation'
         write(6,*) 'There are availables spin-orbital pseudopotentials'
         write(6,*) 'Spin-orbital interaction will not be included in ',
     .     'solid-state calculation'
           elseif (irel.eq.'isp') then
         write(6,*) 'Pseudopotential generated from an ',
     .  'atomic spin-polarized calculation'
           endif 
 

           if (nicore.ne.'nc') then 
       write(6,*) 'Pseudopotential includes a core correction',
     .    'which will be not included in the solid-state calculation'
           endif 
           write(6,*) ' '
           write(6,*) 'Pseudopotential generation method:'
           write(6,*)  method(1),(method(i),i=3,6)



          
          
            
           nrval=nr+1
           if(nrval.gt.nrmax) then
              write(6,*) 'ERROR!!!!'
              write(6,*) 'Nrmax must be increased to at least',nrval
              stop
           endif

           read(1) (rofi(ir),ir=2,nrval)
           rofi(1)=0.0d0

           lmxq=-1
           lmax=-1
           do 20 ndown=1,npotd
               read(1) l,(vps(ir,l),ir=2,nrval)
               do 19 ir=2,nrval
                    vps(ir,l)=vps(ir,l)/rofi(ir)
  19           continue 
               if(l.gt.lmax) lmax=l
  20       continue
           write(6,*) 'Lmax=',l
           write(6,*) 'Valence configuration:'
           do l=0,lmax
              write(6,*) text(2*l+1)
           enddo
           write(6,*) 
           ispin=0
           do 22 nup=1,npotu
              read(1) l  
  22       continue 
           

C******* WE READ THE CORE CORRECTION CHARGE DENSITY ****************


        r2=rofi(2)/(rofi(3)-rofi(2))

        read(1) (chcore(ir),ir=2,nrval)
        chcore(1)=chcore(2)-(chcore(3)-chcore(2))*r2


C******** WE READ THE PSEUDO VALENCE DENSITY ***************************
   
          read(1) (rho(ir),ir=2,nrval)  
          rho(1)=rho(2)-(rho(3)-rho(2))*r2




          close(1)


          nodd=mod(nrval,2)
          nrval=nrval-1+nodd
          

C****CONSTRUCTION OF THE KLEINMAN-BYLANDER PROYECTOR FUNCTIONS**************

C       calculate drdi

          rpb=b
          ea=dexp(a)
          rofi(1)=0.0d0
          do 40 ir=1,nrval
             drdi(ir)=a*rpb
             s(ir)=dsqrt(a*rpb)
             rpb=rpb*ea
  40      continue

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
   

C    CALCULATION OF THE VALENCE SCREENING POTENTIAL 
C             FROM THE READED CHARGE DENSITY

          call vhrtre(rho,ve,rofi,drdi,s,nrval,a)

C   ADD THE EXCHANGE-CORRELATION POTENTIAL

C*****Choosing the adecuate functional for the exchange-correlation****

        if(icorr.eq.'ca') then
           write(6,*) 'Ceperly-Alder for exchange-correlation '
           functl='LDA'
           author='ca'
        elseif(icorr.eq.'pw') then   
           write(6,*) 'PW for exchange-correlation '
           functl='LDA'
           author='pw92'
         elseif(icorr.eq.'pb') then 
           write(6,*) 'GGA-PB for exchange-correlation '
           functl='GGA'
           author='pbe'
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
       

          call atomxc(functl,author,irelt,
     .             nrval,nrmax,rofi,1,auxrho,
     .             ex,ec,dx,dc,vxc)



          do ir=2,nrval
             ve(ir)=ve(ir)+vxc(ir)
          enddo
           
C CALCULATION OF THE KLEINMAN-BYLANDER PROYECTOR FUNCTIONS

         do 42  ir=2,nrval
             s(ir)=drdi(ir)*drdi(ir)
  42     continue
         s(1)=s(2)

         write(6,*) 'Calculation of the Kleinman-Bylander projectors'
         write(6,*) 'Eigenvalues:'

         a2b4=a*a*0.25d0
         do 45 l=0,lmax

            do 43 ir=2,nrval
               r2=(rofi(ir))**2
               vtot=vps(ir,l)+ve(ir)+dble(l*(l+1))/r2
               h(ir)=vtot*s(ir)+a2b4
  43        continue
            h(1)=h(2)


            nnode=1
            nprin=l+1
            e=-(zval/dble(nprin))**2
            z=zval
            dr=-1.0d5
            rmax=rofi(nrval)
            call egofv(h,s,nrval,e,g,y,l,z,a,b,rmax,
     .        nprin,nnode,dr)

            write(6,*) 'L=',l,'EL=',e

            eigen(l)=e

            do 44 ir=2,nrval
                 r=rofi(ir)
                 f=g(ir)
                 dsq=dsqrt(drdi(ir))
                 f=f*dsq
                 rphi(ir,l)=f
   44       continue
            rphi(1,l)=rphi(2,l)
 
   45     continue
           
C       checking normalization of the calculated wave functions 

           do 50 l=0,lmax
             dnrm=0.0d0   
             do 49 ir=2,nrval
                 phi=rphi(ir,l)
                 dnrm=dnrm+phi*phi*drdi(ir)
  49         continue
             dnrm=dsqrt(dnrm)
             ddnrm=dabs(dnrm-1.0d0)
         
             if (ddnrm.gt.1.0d-5) then 
                write(6,*) 'WARNING!!!'
                write(6,*) 'Eigenstate for l=',l,'is not normalized'
                write(6,*) 'Norm=',dnrm
             endif
  50       continue


C CALCULATE EIGENVALUES OF LOCAL POTENTIAL FOR GHOST ANALYSIS
           do lloc=0,lmax
           write(6,*) ' '
           write(6,*)'**********************************************'
           write(6,*) ' '
           write(6,*) 'We take the pseudopotential with L=',lloc
     .      ,'as local'
           do 48 l=0,lmax
              do 46 ir=2,nrval
                 r2=rofi(ir)**2
                 vtot=vps(ir,lloc)+ve(ir)+dble(l*(l+1))/r2
                 h(ir)=vtot*s(ir)+a2b4
 46       continue
          do 47 nnode=1,2
              nprin=l+1
              e=-(zval/dble(nprin))**2
              z=zval
              dr=-1.0d5
              rmax=rofi(nrval)
              call egofv(h,s,nrval,e,g,y,l,z,a,b,rmax,
     .        nprin,nnode,dr)
               elocal(nnode,l)=e
 47       continue
 48     continue

           write(6,*) ' '
           write(6,*) 'KB-projectors cosine:'
           do 60 l=0,lmax

             dnrm=0.0d0
             avgv=0.0d0
             do 59 ir=2,nrval
               r=rofi(ir)
               vl=(vps(ir,l)-vps(ir,lloc))
               phi=rphi(ir,l)/r
               vphi=vl*phi
               dnrm=dnrm+vphi*vphi*r*r*drdi(ir)
               avgv=avgv+vphi*phi*r*r*drdi(ir)
  59         continue


             dkbcos(l)=dnrm/(avgv+1.0d-20) 
             write(6,*) 'kbcos(',l,')',dkbcos(l)
             dnrm=1.0d0/(dsqrt(dnrm)+1.0d-20)



             do 58 ir=2,nrval
               r=rofi(ir)
               vl=(vps(ir,l)-vps(ir,lloc))
               phi=rphi(ir,l)/r
               vphi=vl*phi*dnrm
               h(ir)=vphi
  58         continue
             h(1)=h(2)
             
  60       continue 
           close(21)
           write(6,*) ' '

C*******************GHOST ANALYSIS*****************************************
           do 65 l=0,lmax
             if(dkbcos(l).gt.0.0d0) then 

               if(eigen(l).gt.elocal(2,l)) then
                 write(6,*) 'ATTENTION!!!!!!!'
                 write(6,*) 'GHOST STATE FOR L=',L
               else 
                 write(6,*) 'No ghost state for l=',l
               endif

             elseif(dkbcos(l).lt.0d0) then

               if(eigen(l).gt.elocal(1,l)) then
                  write(6,*) 'ATTENTION!!!!!!!'
                  write(6,*) 'GHOST STATE FOR L=',L
               else
                  write(6,*) 'No ghost state for l=',l
               endif

             elseif(dkbcos(l).eq.0.0d0) then

               write(6,*) 'vps = vlocal, no ghost for l=',l
         
             endif

   65      continue
           enddo        


           end 
        
      SUBROUTINE EGOFV(H,S,N,E,G,Y,L,Z,A,B,RMAX,NPRIN,NNODE,DR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  EIOFV DETERMINES THE EIGENENERGY AND WAVEFUNCTION CORRESPONDING
C  TO A PARTICULAR L, PRINCIPAL QUANTUM NUMBER AND BOUNDARY CONDITION.
C
C  TWO FUNDAMENTAL TECHNIQUES ARE USED TO LOCATE THE SOLUTION:
C       1) NODE COUNTING AND BISECTION
C       2) VARIATIONAL ESTIMATE BASED ON A SLOPE DISCONTINUITY IN PSI
C  THE ARGUMENTS ARE DEFINED AS FOLLOWS:
C       H,S: G" = (H-E*S)*G
C       NR: MAXIMUM ALLOWED NUMBER OF RADIAL POINTS
C       E: E(I) IS THE I-TH ENERGY FOUND
C       NE: NUMBER OF ENERGIES FOUND
C       L: THE ANGULAR MOMENTUM
C       NCOR: THE NUMBER OF LOWER-ENERGY STATES
C
C  THE INDIVIDUAL ENERGIES ARE RESOLVED BY PERFORMING A FIXED NUMBER
C  OF BISECTIONS AFTER A GIVEN EIGENVALUE HAS BEEN ISOLATED
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(N),S(N),G(N),Y(*)
      DATA TOL   /1.D-5/
      NCOR=NPRIN-L-1
      N1=NNODE
      N2=NNODE-1
      E1=E
      E2=E
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  THE LABELS 1 AND 2 REFER TO THE BISECTION PROCESS, DEFINING THE
C  RANGE IN WHICH THE DESIRED SOLUTION IS LOCATED.  THE INITIAL
C  SETTINGS OF N1, N2, E1 AND E2 ARE NOT CONSISTENT WITH THE BISECTION
C  ALGORITHM; THEY ARE SET TO CONSISTENT VALUES WHEN THE DESIRED
C  ENERGY INTERVAL HAS BEEN LOCATED.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DEL=5.D-1
      DE=0.D0
      NITER = 0
 1    NITER = NITER + 1
      IF(NITER.GT.40) GO TO 3
      ET=E+DE
C  THE FOLLOWING LINE IS THE FUNDAMENTAL "BISECTION"
      E=0.5*(E1+E2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  THE FOLLOWING CONCATENATION OF LOGICAL ORS ENSURES THAT NODE
C  COUNTING IS USED UNLESS THE PREVIOUS INTEGRATION OF THE RADIAL
C  EQ PRODUCED BOTH THE CORRECT NUMBER OF NODES AND A SENSIBLE
C  PREDICTION FOR THE ENERGY.
C
C     SENSIBLE MEANS THAT ET MUST BE GREATER THAN E1 AND LESS THAN E2
C     CORRECT NUMBER OF NODES MEANS THAT NT = NNODE OR NNODE-1.
C
C     LEAVING E SET TO ITS BISECTION VALUE, AND TRANSFERING TO
C     THE CALL TO YOFE MEANS THAT WE ARE PERFORMING BISECTION,
C     WHEREAS SETTING E TO ET IS USE OF THE VARIATIONAL ESTIMATE.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(ET.LE.E1 .OR. ET.GE.E2 .OR.
     1       NT.LT.NNODE-1 .OR. NT.GT.NNODE) GO TO 2
      E=ET
      IF(DABS(DE).LT.TOL) GO TO 6
 2    CALL YOFE(E,DE,DR,RMAX,H,S,Y,N,L,NCOR,NT,Z,A,B)
C     WRITE(6,101) L,DR,N1,NT,NNODE,N2,E1,E,E2,DE
C101  FORMAT('  L     DR     N1  NT   N  N2       E1           E',
C    1       '          E2          DE'/I3,D10.3,4I4,4F12.5)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  YOFE INTEGRATES THE SCHRO EQ.; NOW THE BISECTION LOGIC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(NT.GE.NNODE) GO TO 5
C  TOO FEW NODES; SET E1 AND N1
      E1=E
      N1=NT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  AT THIS POINT, WE HAVE JUST SET THE BOTTOM OF THE BISECTION RANGE;
C  IF THE TOP IS ALSO SET, WE PROCEDE.  IF THE TOP OF THE RANGE HAS NOT
C  BEEN SET, IT MEANS THAT WE HAVE YET TO FIND AN E GREATER THAN THE
C  DESIRED ENERGY.  THE UPPER END OF THE RANGE IS EXTENDED.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(N2.GE.NNODE) GO TO 1
      DEL=DEL*2.D0
      E2=E1+DEL
      GO TO 1
C  TOO MANY NODES; SET E2 AND N2
 5    E2=E
      N2=NT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  AT THIS POINT, WE HAVE JUST SET THE TOP OF THE BISECTION RANGE;
C  IF THE TOP IS ALSO SET, WE PROCEDE.  IF THE TOP OF THE RANGE HAS
C  NOT BEEN SET, IT MEANS THAT WE HAVE YET TO FIND AN E LESS THAN THE
C  DESIRED ENERGY.  THE LOWER END OF THE RANGE IS EXTENDED.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(N1.LT.NNODE) GO TO 1
      DEL=DEL*2.D0
      E1=E2-DEL
      GO TO 1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  THE NUMEROV METHOD USES A TRANSFORMATION OF THE RADIAL WAVE FCN.
C  THAT WE CALL "Y".  HAVING LOCATED THE EIGENENERGY, WE TRANSFORM
C  Y TO "G", FROM WHICH THE DENSITY IS EASILY CONSTRUCTED.
C  FINALLY, THE CALL TO "NRMLZG" NORMALIZES G TO ONE ELECTRON.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 6    G(1) = 0.D0
      DO 7 I=2,N
           T=H(I)-E*S(I)
           G(I)=Y(I)/(1.D0-T/12.D0)
 7    CONTINUE
      CALL NRMLZG(G,S,N,N)
      RETURN
 3    WRITE(6,4) Z,L,NNODE,E,DE
 4    FORMAT(' EGOFV: TOO MANY ITERATIONS; EXECUTION STOPPING'/
     1       ' Z=',F3.0,'  L=',I2,'  NNODE=',I2,'  E=',F12.5,
     2       '  DE=',F12.5)
      STOP 8
      END





      SUBROUTINE YOFE(E,DE,DR,RMAX,H,S,Y,NMAX,L,NCOR,NNODE,Z,A,B)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   YOFE INTEGRATES THE RADIAL SCHRODINGER EQN USING THE NUMEROV
C   METHOD.
C
C   THE ARGUMENTS ARE DEFINED AS FOLLOWS:
C       E IS THE OLD ENERGY(OVERWRITTEN) BY THE NEW ENERGY
C       DE IS THE E CHANGE PREDICTED TO ELIM THE KINK IN PSI
C       DR IS THE LOG DERIV (THE BOUNDARY CONDITION)
C       G" = (H-ES)G (ALL DIAGONAL IN I (RADIUS) )
C       Y IS THE NUMEROV INDEPENDENT VARIABLE Y = G - G"/12
C       N IS THE NUMBER OF RADIAL MESH POINTS
C       L IS THE ANGULAR MOMENTUM
C       NCOR IS THE NUMBER OF STATES OF LOWER ENERGY
C       NNODE IS 1 + THE NUMBER OF INTERIOR NODES IN PSI
C       Z IS THE ATOMIC NUMBER
C       A AND B SPECIFY THE RADIAL MESH R(I)=(EXP(A*(I-1))-1)*B
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION H(NMAX),S(NMAX),Y(NMAX)
      ZDR = Z*A*B
      N=NMAX
 8    IF( H(N)-E*S(N) .LT. 1.D0 ) GO TO 9
      Y(N)=0.D0
      N=N-1
      GO TO 8
 9    CONTINUE
      CALL BCORGN(E,H,S,L,ZDR,Y2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  BCORGN COMPUTES Y2, WHICH EMBODIES THE BOUNDARY CONDITION
C  SATISFIED BY THE RADIAL WAVE FUNCTION AT THE ORIGIN
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
      KNK=N
      CALL NUMOUT(E,H,S,Y,NCOR,KNK,NNODE,Y2,G,GSG,X)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  THE OUTWARD INTEGRATION IS NOW COMPLETE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     WE FIRST DECIDE IF THE KINETIC ENERGY IS SUFFICIENTLY NON
C     NEGATIVE TO PERMIT USE OF THE NUMEROV EQ AT RMAX.  IF
C     IT IS NOT, THEN ZERO-VALUE BOUNDARY CONDITION IS USED
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      YN=0.D0
      IF(N.LT.NMAX .OR. DABS(DR).GT.1.D3) GO TO 7
      CALL BCRMAX(E,DR,RMAX,H,S,N,YN,A,B)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  BCRMAX COMPUTES YN, WHICH EMBODIES THE BOUNDARY CONDITION
C  SATISFIED BY THE RADIAL WAVE FUNCTION AT RMAX
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
 7    CALL NUMIN(E,H,S,Y,N,NNDIN,YN,GIN,GSGIN,XIN,KNK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  NUMIN PERFORMS THE INWARD INTEGRATION BY THE NUMEROV METHOD
C
C  THE ENERGY INCREMENT IS NOW EVALUATED FROM THE KINK IN PSI
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
      RATIO = G/GIN
      XIN=XIN*RATIO
      GSG=GSG+GSGIN*RATIO*RATIO
      T=H(KNK)-E*S(KNK)
      DE=G*(X+XIN+T*G)/GSG
      NNODE=NNODE+NNDIN
      IF(DE.LT.0.D0) NNODE=NNODE+1
      DO 6 I=KNK,N
         Y(I) = Y(I)*RATIO
 6    CONTINUE
      RETURN
      END
      SUBROUTINE NRMLZG(G,S,N)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   NRMLZG NORMALIZES THE SUPPLIED RADIAL WAVE FUNCTION
C
C   THE ARGUMENTS ARE DEFINED AS FOLLOWS:
C       G IS THE RADIAL WAVE FUNCTION APPROPRIATE TO THE NUMEROV
C             REPRESENTATION OF THE RADIAL SCHRODINGER EQUATION
C             THAT IS, THE RADIAL FCN R(R) = (DRDI)**1/2 G(I) / R(I)
C       G" = (H-ES)G (ALL DIAGONAL IN I (RADIUS) )
C       N1 IS THE NUMBER OF RADIAL MESH POINTS CORRESPONDING TO
C             THE PORTION OF THE RADIAL MESH ON WHICH THE NORM
C             IS DEFINED
C       N2 IS THE NUMBER OF RADIAL MESH POINTS CORRESPONDING TO
C             THE PORTION OF THE RADIAL MESH ON WHICH THE WAVE
C             FUNCTION IS DEFINED.  FOR THE INTENDED USE OF THIS
C             ROUTINE, N1 = NRVAL AND N2 = NRCOR
C       A AND B ARE THE RADIAL MESH PARAMETERS
C             R(I) = ( EXP(A*(I-1)) - 1 ) * B
C             (DR/DI = A*B AT THE ORIGIN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION S(*),G(*),NORM
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  DETERMINE THE NORM OF G(I) USING SIMPSON'S RULE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (MOD(N,2).NE.1) WRITE(6,*) ' NRMLZG: N SHOULD BE ODD. N =',N
      NORM = 0.D0
      NM1 = N - 1
      DO 2 I = 2,NM1,2
         NORM=NORM + G(I)*S(I)*G(I)
 2    CONTINUE
      NORM = NORM * 2.D0
      NM2  = N - 2
      DO 3 I = 3,NM2,2
         NORM=NORM + G(I)*S(I)*G(I)
 3    CONTINUE
      NORM = NORM * 2.D0
      NM2  = N - 2
      DO 4 I = 1,N,NM1
         NORM=NORM + G(I)*S(I)*G(I)
 4    CONTINUE
      NORM = NORM/3.D0
      SRNRM = DSQRT(NORM)
      DO 5 I=1,N
         G(I) = G(I)/SRNRM
 5    CONTINUE
      RETURN
      END






      SUBROUTINE BCORGN(E,H,S,L,ZDR,Y2)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   YOFE INTEGRATES THE RADIAL SCHRODINGER EQN USING THE NUMEROV
C   METHOD.
C
C   THE ARGUMENTS ARE DEFINED AS FOLLOWS:
C       E IS THE OLD ENERGY(OVERWRITTEN) BY THE NEW ENERGY
C       DE IS THE E CHANGE PREDICTED TO ELIM THE KINK IN PSI
C       DR IS THE LOG DERIV (THE BOUNDARY CONDITION)
C       G" = (H-ES)G (ALL DIAGONAL IN I (RADIUS) )
C       Y IS THE NUMEROV INDEPENDENT VARIABLE Y = G - G"/12
C       N IS THE NUMBER OF RADIAL MESH POINTS
C       L IS THE ANGULAR MOMENTUM
C       NNODE IS 1 + THE NUMBER OF INTERIOR NODES IN PSI
C       Z IS THE ATOMIC NUMBER
C       A AND B SPECIFY THE RADIAL MESH R(I)=(EXP(A*(I-1))-1)*B
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION H(*),S(*)
C
C   THE QUANTITY CALLED D(I) IN THE PROGRAM IS ACTUALLY THE INVERSE
C   OF THE DIAGONAL OF THE TRI-DIAGONAL NUMEROV MATRIX
C
      T2=H(2)-E*S(2)
      D2=-(24.D0+10.D0*T2)/(12.D0-T2)
C
C=================================================================
C  THE FOLLOWING SECTION DEALS WITH THE FACT THAT THE INDEPENDENT
C  VARIABLE "Y" IN THE NUMEROV EQUATION IS NOT ZERO AT THE ORIGIN
C  FOR L LESS THAN 2
C  THE L=0 SOLUTION G VANISHES, BUT THE FIRST AND SECOND
C  DERIVATIVES ARE FINITE, MAKING THE NUMEROV VARIABLE Y FINITE
C  THE L=1 SOLUTION G VANISHES, AND G' ALSO VANISHES, BUT
C  THE SECOND DERIVATIVE G" IS FINITE MAKING Y FINITE.  FOR L > 1,
C  G AND ITS FIRST TWO DERIVATIVES VANISH, MAKING Y ZERO.
C=================================================================
      IF(L.GE.2) GOTO 3
      IF(L.GT.0) GOTO 1
      C0=ZDR/6.D0
      C0=C0/(1.D0-0.75*ZDR)
      GO TO 2
 1    C0=1.D0/12.D0
      C0=-C0*8.D0/3.D0
 2    C1=C0*(12.D0+13.D0*T2)/(12.D0-T2)
      T3=H(3)-E*S(3)
      C2=-5.D-1*C0*(24.D0-T3)/(12.D0-T3)
      D2=(D2-C1)/(1.D0-C2)
 3    Y2=-1.D0/D2
      RETURN
      END



      SUBROUTINE BCRMAX(E,DR,RMAX,H,S,N,YN,A,B)
C
C 22.7.85
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION  H(*),S(*),
     .   E,DR,RMAX,YN,A,B,TNM1,TN,TNP1,BETA,DG,C1,C2,C3,DN
C
C     WRITE(6,*) 'BCRMAX:',DR
      TNM1=H(N-1)-E*S(N-1)
      TN  =H(N  )-E*S(N  )
      TNP1=H(N+1)-E*S(N+1)
      BETA=1.D0+B/RMAX
      DG=A*BETA*(DR+1.D0-5.D-1/BETA)
C
C
      C2=24.D0*DG/(12.D0-TN)
      DN=-(24.D0+10.D0*TN)/(12.D0-TN)
C
      C1= (1.D0-TNM1/6.D0)/(1.D0-TNM1/12.D0)
      C3=-(1.D0-TNP1/6.D0)/(1.D0-TNP1/12.D0)
      YN=-(1.D0-C1/C3)/(DN-C2/C3)
C
C
      RETURN
      END



      SUBROUTINE NUMIN(E,H,S,Y,N,NNODE,YN,G,GSG,X,KNK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   YOFE INTEGRATES THE RADIAL SCHRODINGER EQN USING THE NUMEROV
C   METHOD.
C
C   THE ARGUMENTS ARE DEFINED AS FOLLOWS:
C       E IS THE OLD ENERGY(OVERWRITTEN) BY THE NEW ENERGY
C       DE IS THE E CHANGE PREDICTED TO ELIM THE KINK IN PSI
C       DR IS THE LOG DERIV (THE BOUNDARY CONDITION)
C       G" = (H-ES)G (ALL DIAGONAL IN I (RADIUS) )
C       Y IS THE NUMEROV INDEPENDENT VARIABLE Y = G - G"/12
C       N IS THE NUMBER OF RADIAL MESH POINTS
C       L IS THE ANGULAR MOMENTUM
C       NNODE IS 1 + THE NUMBER OF INTERIOR NODES IN PSI
C       Z IS THE ATOMIC NUMBER
C       A AND B SPECIFY THE RADIAL MESH R(I)=(EXP(A*(I-1))-1)*B
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION H(N),S(N),Y(N)
      Y(N)=YN
      T=H(N)-E*S(N)
      G=Y(N)/(1.D0-T/12.D0)
      GSG=G*S(N)*G
      I=N-1
      Y(I)=1.D0
      T=H(I)-E*S(I)
      G=Y(I)/(1.D0-T/12.D0)
      GSG=GSG+G*S(I)*G
      X=Y(I)-Y(N)
      NNODE=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  BEGIN THE INWARD INTEGRATIONBY THE NUMEROV METHOD
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
 1    X=X+T*G
      I=I-1
      Y(I)=Y(I+1)+X
      IF( Y(I)*Y(I+1) .LT. 0.D0) NNODE=NNODE+1
      T=H(I)-E*S(I)
      G=Y(I)/(1.D0-T/12.D0)
      GSG=GSG+G*S(I)*G
      IF(I.GT.KNK) GO TO 1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  THE LAST STATEMENT DEFINES THE KINK RADIUS AS THE POINT WHERE
C  PSI FIRST TURNS DOWNWARD.  THIS USUALLY MEANS AT THE OUTERMOST
C  MAXIMUM
C
C  THE INWARD INTEGRATION IS NOW COMPLETE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
      KNK=I
      RETURN
      END




      SUBROUTINE NUMOUT(E,H,S,Y,NCOR,KNK,NNODE,Y2,G,GSG,X)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   YOFE INTEGRATES THE RADIAL SCHRODINGER EQN USING THE NUMEROV
C   METHOD.
C
C   THE ARGUMENTS ARE DEFINED AS FOLLOWS:
C       E IS THE OLD ENERGY(OVERWRITTEN) BY THE NEW ENERGY
C       DE IS THE E CHANGE PREDICTED TO ELIM THE KINK IN PSI
C       DR IS THE LOG DERIV (THE BOUNDARY CONDITION)
C       G" = (H-ES)G (ALL DIAGONAL IN I (RADIUS) )
C       Y IS THE NUMEROV INDEPENDENT VARIABLE Y = G - G"/12
C       N IS THE NUMBER OF RADIAL MESH POINTS
C       L IS THE ANGULAR MOMENTUM
C       NNODE IS 1 + THE NUMBER OF INTERIOR NODES IN PSI
C       Z IS THE ATOMIC NUMBER
C       A AND B SPECIFY THE RADIAL MESH R(I)=(EXP(A*(I-1))-1)*B
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION H(KNK),S(KNK),Y(KNK)
      Y(1)=0.D0
      Y(2)=Y2
      T=H(2)-E*S(2)
      G=Y(2)/(1.D0-T/12.D0)
      GSG=G*S(2)*G
      Y(3)=1.D0
      T=H(3)-E*S(3)
      G=Y(3)/(1.D0-T/12.D0)
      GSG=GSG+G*S(3)*G
      X=Y(3)-Y(2)
      I=3
      NNODE=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  BEGIN THE OUTWARD INTEGRATIONBY THE NUMEROV METHOD
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
      NM4=KNK-4
 1    XL=X
      X=X+T*G
      I=I+1
      Y(I)=Y(I-1)+X
C     WRITE(6,300) I,Y(I),X,T,H(I),S(I)
C300  FORMAT(I5,5D14.5)
      IF( Y(I)*Y(I-1) .LT. 0.D0) NNODE=NNODE+1
      T=H(I)-E*S(I)
      G=Y(I)/(1.D0-T/12.D0)
      GSG=GSG+G*S(I)*G
      IF(I.EQ.NM4) GO TO 2
      IF(NNODE.LT.NCOR) GO TO 1
      IF(XL*X.GT.0.D0) GO TO 1
 2    KNK=I
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  THE OUTWARD INTEGRATION IS NOW COMPLETE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
      RETURN
      END



      SUBROUTINE VHRTRE(RHO,V,R,DRDI,SRDRDI,NR,A)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   VHRTRE CONSTRUCTS THE ELECTROSTATIC POTENTIAL DUE TO A SUPPLIED
C   ELECTRON DENSITY.  THE NUMEROV METHOD IS USED TO INTEGRATE
C   POISSON'S EQN.
C
C   DESCRIPTION OF ARGUMENTS:
C      RHO....4*PI*R**2 * THE ELECTRON DENSITY FOR WHICH WE CALCULATING
C             THE ELECTROSTATIC POTENTIAL
C      V......THE ELECTROSTATIC POTENTIAL DUE TO THE ELECTRON DENSITY
C             RHO.  THE CONSTANTS OF INTEGRATION ARE FIXED SO THAT THE
C             POTENTIAL TENDS TO A CONSTANT AT THE ORIGIN AND TO
C             2*Q/R AT R=R(NR), WHERE Q IS THE INTEGRATED CHARGE
C             CONTAINED IN RHO(R)
C      R......THE RADIAL MESH R(I) = B*(EXP(A(I-1))-1)
C      NR.....THE NUMBER OF RADIAL MESH POINTS
C      DRDI...DR(I)/DI
C      SRDRDI.SQRT(DR/DI)
C      A......THE PARAMETER APPEARING IN R(I) = B*(EXP(A(I-1))-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION RHO(*),V(*),R(*),DRDI(*),SRDRDI(*)
      NRM1=NR-1
      NRM2=NR-2
      A2BY4=A*A/4.D0
      YBYQ=1.D0-A*A/48.D0
      QBYY=1.D0/YBYQ
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  SIMPSON'S RULE IS USED TO PERFORM TWO INTEGRALS OVER THE ELECTRON
C  DENSITY.  THE TOTAL CHARGE QT IS USED TO FIX THE POTENTIAL AT R=R(NR)
C  AND V0 (THE INTEGRAL OF THE ELECTRON DENSITY DIVIDED BY R) FIXES
C  THE ELECTROSTATIC POTENTIAL AT THE ORIGIN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      V0=0.D0
      QT=0.D0
      DO 1 IR=2,NRM1,2
      DZ=DRDI(IR)*RHO(IR)
      QT=QT+DZ
 1    V0=V0+DZ/R(IR)
      V0=V0+V0
      QT=QT+QT
      DO 2 IR=3,NRM2,2
      DZ=DRDI(IR)*RHO(IR)
      QT=QT+DZ
 2    V0=V0+DZ/R(IR)
      DZ=DRDI(NR)*RHO(NR)
      QT=(QT+QT+DZ)/3.D0
      V0=(V0+V0+DZ/R(NR))/3.D0
      V(1)=2.D0*V0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  THE ELECTROSTATIC POTENTIAL AT R=0 IS SET EQUAL TO
C                       THE AVERAGE VALUE OF RHO(R)/R
C  BEGIN CONSTRUCTION OF THE POTENTIAL AT FINITE R
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IR=2
      T=SRDRDI(IR)/R(IR)
      BETA=DRDI(IR)*T*RHO(IR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  THE NEXT 4 STATEMENTS INDICATE THAT WE FIRST FIND THE PARTICULAR
C  SOLUTION TO THE INHOMOGENEOUS EQN. FOR WHICH Q(2)=0, WE THEN
C  ADD TO THIS PARTICULAR SOLUTION A SOLUTION OF THE HOMOGENEOUS EQN.
C  (A CONSTANT IN V OR A Q PROPORTIONAL TO R)
C  WHICH WHEN DIVIDED BY R IN GOING FROM Q TO V GIVES
C  THE POTENTIAL THE DESIRED COULOMB TAIL OUTSIDE THE ELECTRON DENSITY.
C  THE SIGNIFICANCE OF THE SOLUTION VANISHING AT THE SECOND RADIAL
C  MESH POINT IS THAT, SINCE ALL REGULAR SOLUTIONS OF THE EQUATION
C  FOR Q=R*V VANISH AT THE ORIGIN, THE KNOWLEDGE OF THE SOLUTION
C  VALUE AT THE SECOND MESH POINT PROVIDES THE TWO SOLUTION VALUES
C  REQUIRED TO START THE NUMEROV PROCEDURE.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      X=0.D0
      Y=0.D0
      Q=(Y-BETA/12.D0)*QBYY
      V(IR)=2.D0*T*Q
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  BEGINNING OF THE NUMEROV ALGORITHM
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 3    X=X+A2BY4*Q-BETA
      Y=Y+X
      IR=IR+1
      T=SRDRDI(IR)/R(IR)
      BETA=T*DRDI(IR)*RHO(IR)
      Q=(Y-BETA/12.D0)*QBYY
      V(IR)=2.D0*T*Q
      IF(IR.LT.NR) GO TO 3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  END OF THE NUMEROV ALGORITHM
C
C  WE HAVE NOW FOUND A PARTICULAR SOLUTION TO THE INHOMOGENEOUS EQN.
C  FOR WHICH Q(R) AT THE SECOND RADIAL MESH POINT EQUALS ZERO.
C  NOTE THAT ALL REGULAR SOLUTIONS TO THE EQUATION FOR Q=R*V
C  VANISH AT THE ORIGIN.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      QPARTC = R(NR)*V(NR)/2.D0
      DZ=QT-QPARTC
      DV=2.D0*DZ/R(NR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  THE LOOP FOLLOWING ADDS THE CONSTANT SOLUTION OF THE HOMOGENEOUS
C  EQN TO THE PARTICULAR SOLUTION OF THE INHOMOGENEOUS EQN.
C  NOTE THAT V(1) IS CONSTRUCTED INDEPENDENTLY
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 4 IR=2,NR
         V(IR)=V(IR)+DV
 4    CONTINUE
      RETURN
      END





      CHARACTER*2 FUNCTION SYMBOL( IZ )

C RETURNS THE SYMBOL OF THE ELEMENT OF ATOMIC NUMBER IZ

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





      CHARACTER*(*) FUNCTION PASTE( STR1, STR2 )

C CONCATENATES THE STRINGS STR1 AND STR2 REMOVING BLANKS IN BETWEEN

      CHARACTER*(*) STR1, STR2
      DO 10 L = LEN( STR1 ), 1, -1
         IF (STR1(L:L) .NE. ' ') GOTO 20
   10 CONTINUE
   20 PASTE = STR1(1:L)//STR2
      END





      SUBROUTINE RATINT(XA,YA,N,X,Y,DY)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NMAX=20,TINY=1.D-15)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      HH=ABS(X-XA(1))
      DO 11 I=1,N
        H=ABS(X-XA(I))
        IF (H.EQ.0.0D0)THEN
          Y=YA(I)
          DY=0.0
          RETURN
        ELSE IF (H.LT.HH) THEN
          NS=I
          HH=H
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)+TINY
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          W=C(I+1)-D(I)
          H=XA(I+M)-X
          T=(XA(I)-X)*D(I)/H
          DD=T-C(I+1)
          IF(DD.EQ.0.0D0)PAUSE
          DD=W/DD
          D(I)=C(I+1)*DD
          C(I)=T*DD
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END





      SUBROUTINE SPLINE(DELT,Y,N,YP1,YPN,Y2,U)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(N),Y2(N),U(N)

      IF (YP1.GT..99E30) THEN
        Y2(1)=0.0D0
        U(1)=0.0D0
      ELSE
        Y2(1)=-0.5D0
        U(1)=(3.0D0/DELT)*((Y(2)-Y(1))/DELT-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=0.5D0
        P=SIG*Y2(I-1)+2.0D0
        Y2(I)=(SIG-1.0D0)/P
        U(I)=(3.0D0*( Y(I+1)+Y(I-1)-2.0D0*Y(I) )/(DELT*DELT)
     *      -SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.0D0
        UN=0.0D0
      ELSE
        QN=0.5D0
        UN=(3.0D0/DELT)*(YPN-(Y(N)-Y(N-1))/DELT)
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END






      SUBROUTINE SPLINT(DELT,YA,Y2A,N,X,Y,DYDX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION YA(N),Y2A(N)
      
      NLO=INT(X/DELT)+1
      NHI=NLO+1

 
      IF (DELT.EQ.0.) PAUSE 'Bad DELT input.'

      A=NHI-X/DELT-1
      B=1.0D0-A
      Y=A*YA(NLO)+B*YA(NHI)+
     *      ((A**3-A)*Y2A(NLO)+(B**3-B)*Y2A(NHI))*(DELT**2)/6.

      DYDX=(YA(NHI)-YA(NLO))/DELT +
     * (-(3*(A**2)-1.0)*Y2A(NLO)+(3*(B**2)-1.0)*Y2A(NHI))*DELT/6.
      RETURN
      END







        SUBROUTINE YLMR (R,LI,MI,YLM,GRYLM)

C COMPUTES REAL SPHERICAL HARMONICS YLM IN THE DIRECTION OF VECTOR R:
C    YLM = C * PLM( COS(THETA) ) * SIN(M*PHI)   FOR   M <  0
C    YLM = C * PLM( COS(THETA) ) * COS(M*PHI)   FOR   M >= 0
C WITH (THETA,PHI) THE POLAR ANGLES OF R, C A POSITIVE NORMALIZATION
C CONSTANT AND PLM ASSOCIATED LEGENDRE POLYNOMIALS.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LMAXD=10,NLD=LMAXD+1,TINY=1.D-30,
     .   ZERO=0.D0,HALF=0.5D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,SIX=6.D0)
      DIMENSION R(3),C(0:NLD*NLD),GRYLM(3)
      SAVE LMAX,C
      DATA LMAX /-1/

C EVALUATE NORMALIZATION CONSTANTS ONCE AND FOR ALL
      IF (LI.GT.LMAXD) THEN
         WRITE(6,*) 'YLM: ERROR PARAMETER LMAXD MUST BE INCREASED'
         STOP
      ENDIF

      IF (LI.GT.LMAX) THEN
         FOURPI=TWO**4*ATAN(ONE)
         DO 20 L=0,LI
            ILM0=L*L+L
            DO 15 M=0,L
               FAC=(2*L+1)/FOURPI
               DO 10 I=L-M+1,L+M
                  FAC=FAC/I
   10          CONTINUE
               C(ILM0+M)=SQRT(FAC)
C              NEXT LINE BECAUSE YLM'S ARE REAL COMBINATIONS OF M AND -M
               IF (M.NE.0) C(ILM0+M)=C(ILM0+M)*SQRT(TWO)
               C(ILM0-M)=C(ILM0+M)
   15       CONTINUE
   20    CONTINUE
         LMAX=LI
      ENDIF

C IF L=0, NO CALCULATIONS ARE REQUIRED

      IF (LI.EQ.0) THEN

      YLM=C(0)
      GRYLM(1)=ZERO
      GRYLM(2)=ZERO
      GRYLM(3)=ZERO

       GOTO 999

      ENDIF

C IF R=0, DIRECTION IS UNDEFINED => MAKE YLM=0 EXCEPT FOR L=0
      R2=R(1)*R(1)+R(2)*R(2)+R(3)*R(3)
      IF (R2.LT.TINY) THEN
            YLM=ZERO
         GOTO 999
      ENDIF
      RSIZE=SQRT(R2)
        
      X=R(1)/RSIZE
      Y=R(2)/RSIZE
      Z=R(3)/RSIZE

C EXPLICIT FORMULAS FOR L=1 AND L=2
      IF (LI.EQ.1) THEN
         IF(MI.EQ.-1) THEN 
           YLM=-C(1)*Y
           GRYLM(1)= C(1)*X*Y/RSIZE
           GRYLM(2)=-C(1)*(ONE-Y*Y)/RSIZE
           GRYLM(3)= C(1)*Z*Y/RSIZE 
         ELSEIF(MI.EQ.0) THEN
           YLM= C(2)*Z
           GRYLM(1)=-C(2)*X*Z/RSIZE
           GRYLM(2)=-C(2)*Y*Z/RSIZE
           GRYLM(3)= C(2)*(ONE-Z*Z)/RSIZE
         ELSEIF(MI.EQ.1) THEN
           YLM=-C(3)*X
           GRYLM(1)=-C(3)*(ONE-X*X)/RSIZE
           GRYLM(2)= C(3)*Y*X/RSIZE
           GRYLM(3)= C(3)*Z*X/RSIZE
         ENDIF
         GOTO 999
      ENDIF
         
       IF (LI.EQ.2) THEN
         IF(MI.EQ.-2)THEN
           YLM= C(4)*SIX*X*Y
           GRYLM(1)=-C(4)*SIX*(TWO*X*X*Y-Y)/RSIZE
           GRYLM(2)=-C(4)*SIX*(TWO*Y*X*Y-X)/RSIZE
           GRYLM(3)=-C(4)*SIX*(TWO*Z*X*Y)/RSIZE
         ELSEIF(MI.EQ.-1)THEN
           YLM=-C(5)*THREE*Y*Z
           GRYLM(1)= C(5)*THREE*(TWO*X*Y*Z)/RSIZE
           GRYLM(2)= C(5)*THREE*(TWO*Y*Y*Z-Z)/RSIZE
           GRYLM(3)= C(5)*THREE*(TWO*Z*Y*Z-Y)/RSIZE
         ELSEIF(MI.EQ.0) THEN
           YLM= C(6)*HALF*(THREE*Z*Z-ONE)
           GRYLM(1)=-C(6)*THREE*(X*Z*Z)/RSIZE
           GRYLM(2)=-C(6)*THREE*(Y*Z*Z)/RSIZE
           GRYLM(3)=-C(6)*THREE*(Z*Z-ONE)*Z/RSIZE
         ELSEIF(MI.EQ.1) THEN
           YLM=-C(7)*THREE*X*Z
           GRYLM(1)= C(7)*THREE*(TWO*X*X*Z-Z)/RSIZE
           GRYLM(2)= C(7)*THREE*(TWO*Y*X*Z)/RSIZE
           GRYLM(3)= C(7)*THREE*(TWO*Z*X*Z-X)/RSIZE
         ELSEIF(MI.EQ.2) THEN
           YLM= C(8)*THREE*(X*X-Y*Y)
           GRYLM(1)=-C(8)*SIX*(X*X-Y*Y-ONE)*X/RSIZE
           GRYLM(2)=-C(8)*SIX*(X*X-Y*Y+ONE)*Y/RSIZE
           GRYLM(3)=-C(8)*SIX*(X*X-Y*Y)*Z/RSIZE
         ENDIF

         GOTO 999
      ENDIF

C GENERAL ALGORITHM BASED ON ROUTINE PLGNDR OF 'NUMERICAL RECIPES'

      MABS=ABS(MI)
      XYSIZE=SQRT(MAX(X*X+Y*Y,TINY))
      COSPHI=X/XYSIZE
      SINPHI=Y/XYSIZE
      COSM=ONE
      SINM=ZERO
      DO M=1,MABS
        COSMM1=COSM
        SINMM1=SINM
        COSM=COSMM1*COSPHI-SINMM1*SINPHI
        SINM=COSMM1*SINPHI+SINMM1*COSPHI
      ENDDO

     
      IF(MI.LT.0) THEN
        PHASE=SINM
        DPHASE=MABS*COSM
      ELSE
        PHASE=COSM
        DPHASE=-MABS*SINM
      ENDIF

         PMM=ONE
         FACT=ONE
       IF(MABS.GT.ZERO) THEN
         DO 30 I=1,MABS
            PMM=-PMM*FACT*XYSIZE
            FACT=FACT+TWO
   30    CONTINUE
       ENDIF

       IF(LI.EQ.MABS) THEN
         PLGNDR=PMM
         DPLG=-LI*Z*PMM/(XYSIZE**2)
       ELSE   
         PMMP1=Z*(2*MABS+1)*PMM
         IF(LI.EQ.MABS+1) THEN
           PLGNDR=PMMP1
           DPLG=-(LI*Z*PMMP1-(MABS+LI)*PMM)/(XYSIZE**2)
         ELSE 

         DO 40 L=MABS+2,LI
            PLL=(Z*(2*L-1)*PMMP1-(L+MABS-1)*PMM)/(L-MABS)
            PMM=PMMP1
            PMMP1=PLL
   40    CONTINUE
         PLGNDR=PLL
         DPLG=-(LI*Z*PLL-(L+MABS-1)*PMM)/(XYSIZE**2)
         ENDIF
       ENDIF         
        

         ILM0=LI*LI+LI
         CMI=C(ILM0+MI)
         YLM=CMI*PLGNDR*PHASE
         GRYLM(1)=-CMI*DPLG*X*Z*PHASE/RSIZE
     .     -CMI*PLGNDR*DPHASE*Y/(RSIZE*XYSIZE**2)

         GRYLM(2)=-CMI*DPLG*Y*Z*PHASE/RSIZE
     .     +CMI*PLGNDR*DPHASE*X/(RSIZE*XYSIZE**2)

         GRYLM(3)= CMI*DPLG*(ONE-Z*Z)*PHASE/RSIZE
   

  999 CONTINUE

      END

      SUBROUTINE ATOMXC( FUNCTL, AUTHOR, IREL,
     .                   NR, MAXR, RMESH, NSPIN, DENS,
     .                   EX, EC, DX, DC, VXC )

C *******************************************************************
C Finds total exchange-correlation energy and potential for a
C spherical electron density distribution.
C This version implements the Local (spin) Density Approximation and
C the Generalized-Gradient-Aproximation with the 'explicit mesh 
C functional' method of White & Bird, PRB 50, 4954 (1994).
C Gradients are 'defined' by numerical derivatives, using 2*NN+1 mesh
C   points, where NN is a parameter defined below
C Coded by L.C.Balbas and J.M.Soler. December 1996. Version 0.5.
C ************************* INPUT ***********************************
C CHARACTER*(*) FUNCTL : Functional to be used:
C              'LDA' or 'LSD' => Local (spin) Density Approximation
C                       'GGA' => Generalized Gradient Corrections
C                                Uppercase is optional
C CHARACTER*(*) AUTHOR : Parametrization desired:
C     'CA' or 'PZ' => LSD Perdew & Zunger, PRB 23, 5075 (1981)
C           'PW92' => LSD Perdew & Wang, PRB, 45, 13244 (1992). This is
C                     the local density limit of the next:
C            'PBE' => GGA Perdew, Burke & Ernzerhof, PRL 77, 3865 (1996)
C                     Uppercase is optional
C INTEGER IREL         : Relativistic exchange? (0=>no, 1=>yes)
C INTEGER NR           : Number of radial mesh points
C INTEGER MAXR         : Physical first dimension of RMESH, DENS and VXC
C REAL*8  RMESH(MAXR)  : Radial mesh points
C INTEGER NSPIN        : NSPIN=1 => unpolarized; NSPIN=2 => polarized
C REAL*8  DENS(MAXR,NSPIN) : Total (NSPIN=1) or spin (NSPIN=2) electron
C                            density at mesh points
C ************************* OUTPUT **********************************
C REAL*8  EX              : Total exchange energy
C REAL*8  EC              : Total correlation energy
C REAL*8  DX              : IntegralOf( rho * (eps_x - v_x) )
C REAL*8  DC              : IntegralOf( rho * (eps_c - v_c) )
C REAL*8  VXC(MAXR,NSPIN) : (Spin) exch-corr potential
C ************************ UNITS ************************************
C Distances in atomic units (Bohr).
C Densities in atomic units (electrons/Bohr**3)
C Energy unit depending of parameter EUNIT below
C ********* ROUTINES CALLED *****************************************
C GGAXC, LDAXC
C *******************************************************************

C Next line is nonstandard but may be suppressed
      IMPLICIT NONE

C Argument types and dimensions
      CHARACTER*(*)     FUNCTL, AUTHOR
      INTEGER           IREL, MAXR, NR, NSPIN
      DOUBLE PRECISION  DENS(MAXR,NSPIN), RMESH(MAXR), VXC(MAXR,NSPIN)
      DOUBLE PRECISION  DC, DX, EC, EX

C Fix the order of the numerical derivatives: the number of radial 
C points used is 2*NN+1
C MAXAUX must be larger than the number of radial mesh points
      INTEGER NN, MAXAUX
      PARAMETER ( NN     =    5 )
      PARAMETER ( MAXAUX = 2000 )

C Fix energy unit:  EUNIT=1.0 => Hartrees,
C                   EUNIT=0.5 => Rydbergs,
C                   EUNIT=0.03674903 => eV
      DOUBLE PRECISION EUNIT
      PARAMETER ( EUNIT = 0.5D0 )

C DVMIN is added to differential of volume to avoid division by zero
      DOUBLE PRECISION DVMIN
      PARAMETER ( DVMIN = 1.D-12 )

C Local variables and arrays
      LOGICAL
     .  GGA
      INTEGER
     .  IN, IN1, IN2, IR, IS, JN
      DOUBLE PRECISION
     .  AUX(MAXAUX), D(2), DECDD(2), DECDGD(3,2), DEXDD(2), DEXDGD(3,2),
     .  DGDM(-NN:NN), DGIDFJ(-NN:NN), DRDM, DVOL, 
     .  EPSC, EPSX, F1, F2, GD(3,2), PI
      EXTERNAL
     .  GGAXC, LDAXC

C Set GGA switch
      IF ( FUNCTL.EQ.'LDA' .OR. FUNCTL.EQ.'lda' .OR.
     .     FUNCTL.EQ.'LSD' .OR. FUNCTL.EQ.'lsd' ) THEN
        GGA = .FALSE.
      ELSEIF ( FUNCTL.EQ.'GGA' .OR. FUNCTL.EQ.'gga') THEN
        GGA = .TRUE.
      ELSE
        WRITE(6,*) 'CELLXC: Unknown functional ', FUNCTL
        STOP
      ENDIF

C Check size of auxiliary array
      IF (GGA .AND. MAXAUX.LT.NR)
     .  STOP 'ATOMXC: Parameter MAXAUX too small'

C Initialize output
      EX = 0
      EC = 0
      DX = 0
      DC = 0
      DO 20 IS = 1,NSPIN
        DO 10 IR = 1,NR
          VXC(IR,IS) = 0
   10   CONTINUE
   20 CONTINUE

C Get number pi
      PI = 4 * ATAN(1.D0)

C Loop on mesh points
      DO 140 IR = 1,NR

C       Find interval of neighbour points to calculate derivatives
        IN1 = MAX(  1, IR-NN ) - IR
        IN2 = MIN( NR, IR+NN ) - IR

C       Find weights of numerical derivation from Lagrange
C       interpolation formula
        DO 50 IN = IN1,IN2
          IF (IN .EQ. 0) THEN
            DGDM(IN) = 0
            DO 30 JN = IN1,IN2
              IF (JN.NE.0) DGDM(IN) = DGDM(IN) + 1.D0 / (0 - JN)
   30       CONTINUE
          ELSE
            F1 = 1
            F2 = 1
            DO 40 JN = IN1,IN2
              IF (JN.NE.IN .AND. JN.NE.0) F1 = F1 * (0  - JN)
              IF (JN.NE.IN)               F2 = F2 * (IN - JN)
   40       CONTINUE
            DGDM(IN) = F1 / F2
          ENDIF
   50   CONTINUE

C       Find dr/dmesh
        DRDM = 0
        DO 60 IN = IN1,IN2
          DRDM = DRDM + RMESH(IR+IN) * DGDM(IN)
   60   CONTINUE

C       Find differential of volume. Use trapezoidal integration rule
        DVOL = 4 * PI * RMESH(IR)**2 * DRDM
C       DVMIN is a small number added to avoid a division by zero
        DVOL = DVOL + DVMIN
        IF (IR.EQ.1 .OR. IR.EQ.NR) DVOL = DVOL / 2
        IF (GGA) AUX(IR) = DVOL

C       Find the weights for the derivative d(gradF(i))/d(F(j)), of
C       the gradient at point i with respect to the value at point j
        IF (GGA) THEN
          DO 80 IN = IN1,IN2
            DGIDFJ(IN) = DGDM(IN) / DRDM
   80     CONTINUE
        ENDIF

C       Find density and gradient of density at this point
        DO 90 IS = 1,NSPIN
          D(IS) = DENS(IR,IS)
   90   CONTINUE
        IF (GGA) THEN
          DO 110 IS = 1,NSPIN
            GD(1,IS) = 0
            GD(2,IS) = 0
            GD(3,IS) = 0
            DO 100 IN = IN1,IN2
              GD(3,IS) = GD(3,IS) + DGIDFJ(IN) * DENS(IR+IN,IS)
  100       CONTINUE
  110     CONTINUE
        ENDIF

C       Find exchange and correlation energy densities and their 
C       derivatives with respect to density and density gradient
        IF (GGA) THEN
          CALL GGAXC( AUTHOR, IREL, NSPIN, D, GD,
     .                EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
        ELSE
          CALL LDAXC( AUTHOR, IREL, NSPIN, D, EPSX, EPSC, DEXDD, DECDD )
        ENDIF

C       Add contributions to exchange-correlation energy and its
C       derivatives with respect to density at all points
        DO 130 IS = 1,NSPIN
          EX = EX + DVOL * D(IS) * EPSX
          EC = EC + DVOL * D(IS) * EPSC
          DX = DX + DVOL * D(IS) * (EPSX - DEXDD(IS))
          DC = DC + DVOL * D(IS) * (EPSC - DECDD(IS))
          IF (GGA) THEN
            VXC(IR,IS) = VXC(IR,IS) + DVOL * ( DEXDD(IS) + DECDD(IS) )
            DO 120 IN = IN1,IN2
              DX= DX - DVOL * DENS(IR+IN,IS) * DEXDGD(3,IS) * DGIDFJ(IN)
              DC= DC - DVOL * DENS(IR+IN,IS) * DECDGD(3,IS) * DGIDFJ(IN)
              VXC(IR+IN,IS) = VXC(IR+IN,IS) + DVOL *
     .               (DEXDGD(3,IS) + DECDGD(3,IS)) * DGIDFJ(IN)
  120       CONTINUE
          ELSE
            VXC(IR,IS) = DEXDD(IS) + DECDD(IS)
          ENDIF
  130   CONTINUE

  140 CONTINUE

C Divide by volume element to obtain the potential (per electron)
      IF (GGA) THEN
        DO 160 IS = 1,NSPIN
          DO 150 IR = 1,NR
            DVOL = AUX(IR)
            VXC(IR,IS) = VXC(IR,IS) / DVOL
  150     CONTINUE
  160   CONTINUE
      ENDIF

C Divide by energy unit
      EX = EX / EUNIT
      EC = EC / EUNIT
      DX = DX / EUNIT
      DC = DC / EUNIT
      DO 180 IS = 1,NSPIN
        DO 170 IR = 1,NR
          VXC(IR,IS) = VXC(IR,IS) / EUNIT
  170   CONTINUE
  180 CONTINUE

      END


      SUBROUTINE EXCHNG( IREL, NSP, DS, EX, VX )

C *****************************************************************
C  Finds local exchange energy density and potential
C  Adapted by J.M.Soler from routine velect of Froyen's 
C    pseudopotential generation program. Madrid, Jan'97. Version 0.5.
C **** Input ******************************************************
C INTEGER IREL    : relativistic-exchange switch (0=no, 1=yes)
C INTEGER NSP     : spin-polarizations (1=>unpolarized, 2=>polarized)
C REAL*8  DS(NSP) : total (nsp=1) or spin (nsp=2) electron density
C **** Output *****************************************************
C REAL*8  EX      : exchange energy density
C REAL*8  VX(NSP) : (spin-dependent) exchange potential
C **** Units ******************************************************
C Densities in electrons/Bohr**3
C Energies in Hartrees
C *****************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.D0,ONE=1.D0,PFIVE=.5D0,OPF=1.5D0,C014=0.014D0)
      DIMENSION DS(NSP), VX(NSP)

       PI=4*ATAN(ONE)
       TRD = ONE/3
       FTRD = 4*TRD
       TFTM = 2**FTRD-2
       A0 = (4/(9*PI))**TRD

C      X-alpha parameter:       
       ALP = 2 * TRD

       IF (NSP .EQ. 2) THEN
         D = DS(1) + DS(2)
         IF (D .LE. ZERO) THEN
           EX = ZERO
           VX(1) = ZERO
           VX(2) = ZERO
           RETURN
         ENDIF
         Z = (DS(1) - DS(2)) / D
         FZ = ((1+Z)**FTRD+(1-Z)**FTRD-2)/TFTM
         FZP = FTRD*((1+Z)**TRD-(1-Z)**TRD)/TFTM 
       ELSE
         D = DS(1)
         IF (D .LE. ZERO) THEN
           EX = ZERO
           VX(1) = ZERO
           RETURN
         ENDIF
         Z = ZERO
         FZ = ZERO
         FZP = ZERO
       ENDIF
       RS = (3 / (4*PI*D) )**TRD
       VXP = -3*ALP/(2*PI*A0*RS)
       EXP = 3*VXP/4
       IF (IREL .EQ. 1) THEN
         BETA = C014/RS
         SB = SQRT(1+BETA*BETA)
         ALB = LOG(BETA+SB)
         VXP = VXP * (-PFIVE + OPF * ALB / (BETA*SB))
         EXP = EXP * (ONE-OPF*((BETA*SB-ALB)/BETA**2)**2) 
       ENDIF
       VXF = 2**TRD*VXP
       EXF = 2**TRD*EXP
       IF (NSP .EQ. 2) THEN
         VX(1) = VXP + FZ*(VXF-VXP) + (1-Z)*FZP*(EXF-EXP)
         VX(2) = VXP + FZ*(VXF-VXP) - (1+Z)*FZP*(EXF-EXP)
         EX    = EXP + FZ*(EXF-EXP)
       ELSE
         VX(1) = VXP
         EX    = EXP
       ENDIF
      END




      SUBROUTINE GGAXC( AUTHOR, IREL, NSPIN, D, GD,
     .                  EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )

C Finds the exchange and correlation energies at a point, and their
C derivatives with respect to density and density gradient, in the
C Generalized Gradient Correction approximation.
C Lengths in Bohr, energies in Hartrees
C Written by L.C.Balbas and J.M.Soler, Dec'96. Version 0.5.

      IMPLICIT          NONE
      CHARACTER*(*)     AUTHOR
      INTEGER           IREL, NSPIN
      DOUBLE PRECISION  D(NSPIN), DECDD, DECDGD, DEXDD, DEXDGD,
     .                  EPSC, EPSX, GD(3,NSPIN)

      IF (AUTHOR.EQ.'PBE' .OR. AUTHOR.EQ.'pbe') THEN
        CALL PBEXC( IREL, NSPIN, D, GD,
     .              EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
      ELSE
        WRITE(6,*) 'GGAXC: Unknown author ', AUTHOR
        STOP
      ENDIF
      END




      SUBROUTINE LDAXC( AUTHOR, IREL, NSPIN, D, EPSX, EPSC, VX, VC )

C Finds the exchange and correlation energies and potentials, in the
C Local (spin) Density Approximation.
C Lengths in Bohr, energies in Hartrees
C Written by L.C.Balbas and J.M.Soler, Dec'96. Version 0.5.

      IMPLICIT          NONE
      CHARACTER*(*)     AUTHOR
      INTEGER           IREL, NSPIN
      DOUBLE PRECISION  D(NSPIN), EPSC, EPSX, VX(NSPIN), VC(NSPIN)

      IF ( AUTHOR.EQ.'CA' .OR. AUTHOR.EQ.'ca' .OR.
     .     AUTHOR.EQ.'PZ' .OR. AUTHOR.EQ.'pz') THEN
        CALL PZXC( IREL, NSPIN, D, EPSX, EPSC, VX, VC )
      ELSEIF ( AUTHOR.EQ.'PW92' .OR. AUTHOR.EQ.'pw92' ) THEN
        CALL PW92XC( IREL, NSPIN, D, EPSX, EPSC, VX, VC )
      ELSE
        WRITE(6,*) 'GGAXC: Unknown author ', AUTHOR
        STOP
      ENDIF
      END



      SUBROUTINE PBEXC( IREL, NSPIN, DENS, GDENS,
     .                  EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Perdew-Burke-Ernzerhof Generalized-Gradient-Approximation.
C Ref: J.P.Perdew, K.Burke & M.Ernzerhof, PRL 77, 3865 (1996)
C Written by L.C.Balbas and J.M.Soler. December 1996. Version 0.5.
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER NSPIN          : Number of spin polarizations (1 or 2)
C REAL*8  DENS(NSPIN)    : Total electron density (if NSPIN=1) or
C                           spin electron density (if NSPIN=2)
C REAL*8  GDENS(3,NSPIN) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(NSPIN)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( DENS(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(NSPIN)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( DENS(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,NSPIN): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,NSPIN): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C EXCHNG, PW92C
C ********************************************************************

      IMPLICIT          NONE
      INTEGER           IREL, NSPIN
      DOUBLE PRECISION  DENS(NSPIN), DECDD(NSPIN), DECDGD(3,NSPIN),
     .                  DEXDD(NSPIN), DEXDGD(3,NSPIN), GDENS(3,NSPIN)

C Internal variables
      INTEGER
     .  IS, IX
      DOUBLE PRECISION
     .  A, BETA, D(2), DADD, DECUDD, DENMIN, 
     .  DF1DD, DF2DD, DF3DD, DF4DD, DF1DGD, DF3DGD, DF4DGD,
     .  DFCDD(2), DFCDGD(3,2), DFDD, DFDGD, DFXDD(2), DFXDGD(3,2),
     .  DHDD, DHDGD, DKFDD, DKSDD, DPDD, DPDZ, DRSDD, 
     .  DS, DSDD, DSDGD, DT, DTDD, DTDGD, DZDD(2), 
     .  EC, ECUNIF, EX, EXUNIF,
     .  F, F1, F2, F3, F4, FC, FX, FOUTHD,
     .  GAMMA, GD(3,2), GDM(2), GDMIN, GDMS, GDMT, GDS, GDT(3),
     .  H, HALF, KAPPA, KF, KFS, KS, MU, PHI, PI, RS, S,
     .  T, THD, THRHLF, TWO, TWOTHD, VCUNIF(2), VXUNIF, ZETA

C Lower bounds of density and its gradient to avoid divisions by zero
      PARAMETER ( DENMIN = 1.D-12 )
      PARAMETER ( GDMIN  = 1.D-12 )

C Fix some numerical parameters
      PARAMETER ( FOUTHD=4.D0/3.D0, HALF=0.5D0,
     .            THD=1.D0/3.D0, THRHLF=1.5D0,
     .            TWO=2.D0, TWOTHD=2.D0/3.D0 )

C Fix some more numerical constants
      PI = 4 * ATAN(1.D0)
      BETA = 0.066725D0
      GAMMA = (1 - LOG(TWO)) / PI**2
      MU = BETA * PI**2 / 3
      KAPPA = 0.804D0

C Translate density and its gradient to new variables
      IF (NSPIN .EQ. 1) THEN
        D(1) = HALF * DENS(1)
        D(2) = D(1)
        DT = MAX( DENMIN, DENS(1) )
        DO 10 IX = 1,3
          GD(IX,1) = HALF * GDENS(IX,1)
          GD(IX,2) = GD(IX,1)
          GDT(IX) = GDENS(IX,1)
   10   CONTINUE
      ELSE
        D(1) = DENS(1)
        D(2) = DENS(2)
        DT = MAX( DENMIN, DENS(1)+DENS(2) )
        DO 20 IX = 1,3
          GD(IX,1) = GDENS(IX,1)
          GD(IX,2) = GDENS(IX,2)
          GDT(IX) = GDENS(IX,1) + GDENS(IX,2)
   20   CONTINUE
      ENDIF
      GDM(1) = SQRT( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
      GDM(2) = SQRT( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
      GDMT   = SQRT( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
      GDMT = MAX( GDMIN, GDMT )

C Find local correlation energy and potential
      CALL PW92C( 2, D, ECUNIF, VCUNIF )

C Find total correlation energy
      RS = ( 3 / (4*PI*DT) )**THD
      KF = (3 * PI**2 * DT)**THD
      KS = SQRT( 4 * KF / PI )
      ZETA = ( D(1) - D(2) ) / DT
      ZETA = MAX( -1.D0+DENMIN, ZETA )
      ZETA = MIN(  1.D0-DENMIN, ZETA )
      PHI = HALF * ( (1+ZETA)**TWOTHD + (1-ZETA)**TWOTHD )
      T = GDMT / (2 * PHI * KS * DT)
      F1 = ECUNIF / GAMMA / PHI**3
      F2 = EXP(-F1)
      A = BETA / GAMMA / (F2-1)
      F3 = T**2 + A * T**4
      F4 = BETA/GAMMA * F3 / (1 + A*F3)
      H = GAMMA * PHI**3 * LOG( 1 + F4 )
      FC = ECUNIF + H

C Find correlation energy derivatives
      DRSDD = - THD * RS / DT
      DKFDD =   THD * KF / DT
      DKSDD = HALF * KS * DKFDD / KF
      DZDD(1) =   1 / DT - ZETA / DT
      DZDD(2) = - 1 / DT - ZETA / DT
      DPDZ = HALF * TWOTHD * ( 1/(1+ZETA)**THD - 1/(1-ZETA)**THD )
      DO 40 IS = 1,2
        DECUDD = ( VCUNIF(IS) - ECUNIF ) / DT
        DPDD = DPDZ * DZDD(IS)
        DTDD = - T * ( DPDD/PHI + DKSDD/KS + 1/DT )
        DF1DD = F1 * ( DECUDD/ECUNIF - 3*DPDD/PHI )
        DF2DD = - F2 * DF1DD
        DADD = - A * DF2DD / (F2-1)
        DF3DD = (2*T + 4*A*T**3) * DTDD + DADD * T**4
        DF4DD = F4 * ( DF3DD/F3 - (DADD*F3+A*DF3DD)/(1+A*F3) )
        DHDD = 3 * H * DPDD / PHI
        DHDD = DHDD + GAMMA * PHI**3 * DF4DD / (1+F4)
        DFCDD(IS) = VCUNIF(IS) + H + DT * DHDD

        DO 30 IX = 1,3
          DTDGD = (T / GDMT) * GDT(IX) / GDMT
          DF3DGD = DTDGD * ( 2 * T + 4 * A * T**3 )
          DF4DGD = F4 * DF3DGD * ( 1/F3 - A/(1+A*F3) ) 
          DHDGD = GAMMA * PHI**3 * DF4DGD / (1+F4)
          DFCDGD(IX,IS) = DT * DHDGD
   30   CONTINUE
   40 CONTINUE

C Find exchange energy and potential
      FX = 0
      DO 60 IS = 1,2
        DS   = MAX( DENMIN, 2 * D(IS) )
        GDMS = MAX( GDMIN, 2 * GDM(IS) )
        KFS = (3 * PI**2 * DS)**THD
        S = GDMS / (2 * KFS * DS)
        F1 = 1 + MU * S**2 / KAPPA
        F = 1 + KAPPA - KAPPA / F1
        CALL EXCHNG( IREL, 1, DS, EXUNIF, VXUNIF )
        FX = FX + DS * EXUNIF * F

        DKFDD = THD * KFS / DS
        DSDD = S * ( -DKFDD/KFS - 1/DS )
        DF1DD = 2 * (F1-1) * DSDD / S
        DFDD = KAPPA * DF1DD / F1**2
        DFXDD(IS) = VXUNIF * F + DS * EXUNIF * DFDD

        DO 50 IX = 1,3
          GDS = 2 * GD(IX,IS)
          DSDGD = (S / GDMS) * GDS / GDMS
          DF1DGD = 2 * MU * S * DSDGD / KAPPA
          DFDGD = KAPPA * DF1DGD / F1**2
          DFXDGD(IX,IS) = DS * EXUNIF * DFDGD
   50   CONTINUE
   60 CONTINUE
      FX = HALF * FX / DT

C Set output arguments
      EX = FX
      EC = FC
      DO 90 IS = 1,NSPIN
        DEXDD(IS) = DFXDD(IS)
        DECDD(IS) = DFCDD(IS)
        DO 80 IX = 1,3
          DEXDGD(IX,IS) = DFXDGD(IX,IS)
          DECDGD(IX,IS) = DFCDGD(IX,IS)
   80   CONTINUE
   90 CONTINUE

      END




      SUBROUTINE PW92C( NSPIN, DENS, EC, VC )

C ********************************************************************
C Implements the Perdew-Wang'92 local correlation (beyond RPA).
C Ref: J.P.Perdew & Y.Wang, PRB, 45, 13244 (1992)
C Written by L.C.Balbas and J.M.Soler. Dec'96.  Version 0.5.
C ********* INPUT ****************************************************
C INTEGER NSPIN       : Number of spin polarizations (1 or 2)
C REAL*8  DENS(NSPIN) : Local (spin) density
C ********* OUTPUT ***************************************************
C REAL*8  EC        : Correlation energy density
C REAL*8  VC(NSPIN) : Correlation (spin) potential
C ********* UNITS ****************************************************
C Densities in electrons per Bohr**3
C Energies in Hartrees
C ********* ROUTINES CALLED ******************************************
C None
C ********************************************************************

C Next line is nonstandard but may be supressed
      IMPLICIT          NONE

C Argument types and dimensions
      INTEGER           NSPIN
      DOUBLE PRECISION  DENS(NSPIN), EC, VC(NSPIN)

C Internal variable declarations
      INTEGER           IG
      DOUBLE PRECISION  A(0:2), ALPHA1(0:2), B, BETA(0:2,4), C,
     .                  DBDRS, DECDD(2), DECDRS, DECDZ, DENMIN, DFDZ,
     .                  DGDRS(0:2), DCDRS, DRSDD, DTOT, DZDD(2),
     .                  F, FPP0, FOUTHD, G(0:2), HALF,
     .                  P(0:2), PI, RS, THD, THRHLF, ZETA

C Fix lower bound of density to avoid division by zero
      PARAMETER ( DENMIN = 1.D-12 )

C Fix some numerical constants
      PARAMETER ( FOUTHD=4.D0/3.D0, HALF=0.5D0,
     .            THD=1.D0/3.D0, THRHLF=1.5D0 )

C Parameters from Table I of Perdew & Wang, PRB, 45, 13244 (92)
      DATA P      / 1.00,     1.00,     1.00     /
      DATA A      / 0.031091, 0.015545, 0.016887 /
      DATA ALPHA1 / 0.21370,  0.20548,  0.11125  /
      DATA BETA   / 7.5957,  14.1189,  10.357,
     .              3.5876,   6.1977,   3.6231,
     .              1.6382,   3.3662,   0.88026,
     .              0.49294,  0.62517,  0.49671 /

C Find rs and zeta
      PI = 4 * ATAN(1.D0)
      IF (NSPIN .EQ. 1) THEN
        DTOT = MAX( DENMIN, DENS(1) )
        ZETA = 0
        RS = ( 3 / (4*PI*DTOT) )**THD
C       Find derivatives dRs/dDens and dZeta/dDens
        DRSDD = - RS / DTOT / 3
        DZDD(1) = 0
      ELSE
        DTOT = MAX( DENMIN, DENS(1)+DENS(2) )
        ZETA = ( DENS(1) - DENS(2) ) / DTOT
        RS = ( 3 / (4*PI*DTOT) )**THD
        DRSDD = - RS / DTOT / 3
        DZDD(1) =   1 / DTOT - ZETA / DTOT
        DZDD(2) = - 1 / DTOT - ZETA / DTOT
      ENDIF

C Find eps_c(rs,0)=G(0), eps_c(rs,1)=G(1) and -alpha_c(rs)=G(2)
C using eq.(10) of cited reference (Perdew & Wang, PRB, 45, 13244 (92))
      DO 20 IG = 0,2
        B = BETA(IG,1) * RS**HALF   +
     .      BETA(IG,2) * RS         +
     .      BETA(IG,3) * RS**THRHLF +
     .      BETA(IG,4) * RS**(P(IG)+1)
        DBDRS = BETA(IG,1) * HALF      / RS**HALF +
     .          BETA(IG,2)                         +
     .          BETA(IG,3) * THRHLF    * RS**HALF +
     .          BETA(IG,4) * (P(IG)+1) * RS**P(IG)
        C = 1 + 1 / (2 * A(IG) * B)
        DCDRS = - (C-1) * DBDRS / B
        G(IG) = - 2 * A(IG) * ( 1 + ALPHA1(IG)*RS ) * LOG(C)
        DGDRS(IG) = - 2*A(IG) * ( ALPHA1(IG) * LOG(C) +
     .                            (1+ALPHA1(IG)*RS) * DCDRS / C )
   20 CONTINUE

C Find f''(0) and f(zeta) from eq.(9)
      C = 1 / (2**FOUTHD - 2)
      FPP0 = 8 * C / 9
      F = ( (1+ZETA)**FOUTHD + (1-ZETA)**FOUTHD - 2 ) * C
      DFDZ = FOUTHD * ( (1+ZETA)**THD - (1-ZETA)**THD ) * C

C Find eps_c(rs,zeta) from eq.(8)
      EC = G(0) - G(2) * F / FPP0 * (1-ZETA**4) +
     .    (G(1)-G(0)) * F * ZETA**4
      DECDRS = DGDRS(0) - DGDRS(2) * F / FPP0 * (1-ZETA**4) +
     .        (DGDRS(1)-DGDRS(0)) * F * ZETA**4
      DECDZ = - G(2) / FPP0 * ( DFDZ*(1-ZETA**4) - F*4*ZETA**3 ) +
     .        (G(1)-G(0)) * ( DFDZ*ZETA**4 + F*4*ZETA**3 )
      
C Find correlation potential
      IF (NSPIN .EQ. 1) THEN
        DECDD(1) = DECDRS * DRSDD
        VC(1) = EC + DTOT * DECDD(1)
      ELSE
        DECDD(1) = DECDRS * DRSDD + DECDZ * DZDD(1)
        DECDD(2) = DECDRS * DRSDD + DECDZ * DZDD(2)
        VC(1) = EC + DTOT * DECDD(1)
        VC(2) = EC + DTOT * DECDD(2)
      ENDIF

      END



      SUBROUTINE PW92XC( IREL, NSPIN, DENS, EPSX, EPSC, VX, VC )

C ********************************************************************
C Implements the Perdew-Wang'92 LDA/LSD exchange correlation
C Ref: J.P.Perdew & Y.Wang, PRB, 45, 13244 (1992)
C Written by L.C.Balbas and J.M.Soler. Dec'96. Version 0.5.
C ********* INPUT ****************************************************
C INTEGER IREL        : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER NSPIN       : Number of spin polarizations (1 or 2)
C REAL*8  DENS(NSPIN) : Local (spin) density
C ********* OUTPUT ***************************************************
C REAL*8  EPSX       : Exchange energy density
C REAL*8  EPSC       : Correlation energy density
C REAL*8  VX(NSPIN)  : Exchange (spin) potential
C REAL*8  VC(NSPIN)  : Correlation (spin) potential
C ********* UNITS ****************************************************
C Densities in electrons per Bohr**3
C Energies in Hartrees
C ********* ROUTINES CALLED ******************************************
C EXCHNG, PW92C
C ********************************************************************

      IMPLICIT          NONE
      INTEGER           IREL, NSPIN
      DOUBLE PRECISION  DENS(NSPIN), EPSX, EPSC, VC(NSPIN), VX(NSPIN)

      CALL EXCHNG( IREL, NSPIN, DENS, EPSX, VX )
      CALL PW92C( NSPIN, DENS, EPSC, VC )
      END



      SUBROUTINE PZXC( IREL, NSP, DS, EX, EC, VX, VC )

C *****************************************************************
C  Perdew-Zunger parameterization of Ceperley-Alder exchange and 
C  correlation. Ref: Perdew & Zunger, Phys. Rev. B 23 5075 (1981).
C  Adapted by J.M.Soler from routine velect of Froyen's 
C    pseudopotential generation program. Madrid, Jan'97. Version 0.5.
C **** Input *****************************************************
C INTEGER IREL    : relativistic-exchange switch (0=no, 1=yes)
C INTEGER NSP     : spin-polarizations (1=>unpolarized, 2=>polarized)
C REAL*8  DS(NSP) : total (nsp=1) or spin (nsp=2) electron density
C **** Output *****************************************************
C REAL*8  EX      : exchange energy density
C REAL*8  EC      : correlation energy density
C REAL*8  VX(NSP) : (spin-dependent) exchange potential
C REAL*8  VC(NSP) : (spin-dependent) correlation potential
C **** Units *******************************************************
C Densities in electrons/Bohr**3
C Energies in Hartrees
C *****************************************************************

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION DS(NSP), VX(NSP), VC(NSP)

       PARAMETER (ZERO=0.D0,ONE=1.D0,PFIVE=.5D0,OPF=1.5D0,PNN=.99D0)
       PARAMETER (PTHREE=0.3D0,PSEVF=0.75D0,C0504=0.0504D0) 
       PARAMETER (C0254=0.0254D0,C014=0.014D0,C0406=0.0406D0)
       PARAMETER (C15P9=15.9D0,C0666=0.0666D0,C11P4=11.4D0)
       PARAMETER (C045=0.045D0,C7P8=7.8D0,C88=0.88D0,C20P59=20.592D0)
       PARAMETER (C3P52=3.52D0,C0311=0.0311D0,C0014=0.0014D0)
       PARAMETER (C0538=0.0538D0,C0096=0.0096D0,C096=0.096D0)
       PARAMETER (C0622=0.0622D0,C004=0.004D0,C0232=0.0232D0)
       PARAMETER (C1686=0.1686D0,C1P398=1.3981D0,C2611=0.2611D0)
       PARAMETER (C2846=0.2846D0,C1P053=1.0529D0,C3334=0.3334D0)
Cray       PARAMETER (ZERO=0.0,ONE=1.0,PFIVE=0.5,OPF=1.5,PNN=0.99)
Cray       PARAMETER (PTHREE=0.3,PSEVF=0.75,C0504=0.0504) 
Cray       PARAMETER (C0254=0.0254,C014=0.014,C0406=0.0406)
Cray       PARAMETER (C15P9=15.9,C0666=0.0666,C11P4=11.4)
Cray       PARAMETER (C045=0.045,C7P8=7.8,C88=0.88,C20P59=20.592)
Cray       PARAMETER (C3P52=3.52,C0311=0.0311,C0014=0.0014)
Cray       PARAMETER (C0538=0.0538,C0096=0.0096,C096=0.096)
Cray       PARAMETER (C0622=0.0622,C004=0.004,C0232=0.0232)
Cray       PARAMETER (C1686=0.1686,C1P398=1.3981,C2611=0.2611)
Cray       PARAMETER (C2846=0.2846,C1P053=1.0529,C3334=0.3334)

C    Ceperly-Alder 'ca' constants. Internal energies in Rydbergs.
       PARAMETER (CON1=1.D0/6, CON2=0.008D0/3, CON3=0.3502D0/3) 
       PARAMETER (CON4=0.0504D0/3, CON5=0.0028D0/3, CON6=0.1925D0/3)
       PARAMETER (CON7=0.0206D0/3, CON8=9.7867D0/6, CON9=1.0444D0/3)
       PARAMETER (CON10=7.3703D0/6, CON11=1.3336D0/3)
Cray       PARAMETER (CON1=1.0/6, CON2=0.008/3, CON3=0.3502/3) 
Cray       PARAMETER (CON4=0.0504/3, CON5=0.0028/3, CON6=0.1925/3)
Cray       PARAMETER (CON7=0.0206/3, CON8=9.7867/6, CON9=1.0444/3)
Cray       PARAMETER (CON10=7.3703/6, CON11=1.3336/3) 

C      X-alpha parameter:
       PARAMETER ( ALP = 2.D0 / 3.D0 )

C      Other variables converted into parameters by J.M.Soler
       PARAMETER ( PI   = 3.14159265358979312D0 )
       PARAMETER ( HALF = 0.5D0 ) 
       PARAMETER ( TRD  = 1.D0 / 3.D0 ) 
       PARAMETER ( FTRD = 4.D0 / 3.D0 )
       PARAMETER ( TFTM = 0.51984209978974638D0 )
       PARAMETER ( A0   = 0.52106176119784808D0 )
       PARAMETER ( CRS  = 0.620350490899400087D0 )
       PARAMETER ( CXP  = - 3.D0 * ALP / (PI*A0) )
       PARAMETER ( CXF  = 1.25992104989487319D0 )

C      Find density and polarization
       IF (NSP .EQ. 2) THEN
         D = DS(1) + DS(2)
         IF (D .LE. ZERO) THEN
           EX = ZERO
           EC = ZERO
           VX(1) = ZERO
           VX(2) = ZERO
           VC(1) = ZERO
           VC(2) = ZERO
           RETURN
         ENDIF
         Z = (DS(1) - DS(2)) / D
         FZ = ((1+Z)**FTRD+(1-Z)**FTRD-2)/TFTM
         FZP = FTRD*((1+Z)**TRD-(1-Z)**TRD)/TFTM 
       ELSE
         D = DS(1)
         IF (D .LE. ZERO) THEN
           EX = ZERO
           EC = ZERO
           VX(1) = ZERO
           VC(1) = ZERO
           RETURN
         ENDIF
         Z = ZERO
         FZ = ZERO
         FZP = ZERO
       ENDIF
       RS = CRS / D**TRD

C      Exchange
       VXP = CXP / RS
       EXP = 0.75D0 * VXP
       IF (IREL .EQ. 1) THEN
         BETA = C014/RS
         SB = SQRT(1+BETA*BETA)
         ALB = LOG(BETA+SB)
         VXP = VXP * (-PFIVE + OPF * ALB / (BETA*SB))
         EXP = EXP *(ONE-OPF*((BETA*SB-ALB)/BETA**2)**2) 
       ENDIF
       VXF = CXF * VXP
       EXF = CXF * EXP

C      Correlation 
       IF (RS .GT. ONE) THEN  
         SQRS=SQRT(RS)
         TE = ONE+CON10*SQRS+CON11*RS
         BE = ONE+C1P053*SQRS+C3334*RS
         ECP = -C2846/BE
         VCP = ECP*TE/BE
         TE = ONE+CON8*SQRS+CON9*RS
         BE = ONE+C1P398*SQRS+C2611*RS
         ECF = -C1686/BE
         VCF = ECF*TE/BE
       ELSE
         RSLOG=LOG(RS)
         ECP=(C0622+C004*RS)*RSLOG-C096-C0232*RS
         VCP=(C0622+CON2*RS)*RSLOG-CON3-CON4*RS
         ECF=(C0311+C0014*RS)*RSLOG-C0538-C0096*RS
         VCF=(C0311+CON5*RS)*RSLOG-CON6-CON7*RS
       ENDIF

C      Find up and down potentials
       IF (NSP .EQ. 2) THEN
         EX    = EXP + FZ*(EXF-EXP)
         EC    = ECP + FZ*(ECF-ECP)
         VX(1) = VXP + FZ*(VXF-VXP) + (1-Z)*FZP*(EXF-EXP)
         VX(2) = VXP + FZ*(VXF-VXP) - (1+Z)*FZP*(EXF-EXP)
         VC(1) = VCP + FZ*(VCF-VCP) + (1-Z)*FZP*(ECF-ECP)
         VC(2) = VCP + FZ*(VCF-VCP) - (1+Z)*FZP*(ECF-ECP)
       ELSE
         EX    = EXP
         EC    = ECP
         VX(1) = VXP
         VC(1) = VCP
       ENDIF

C      Change from Rydbergs to Hartrees
       EX = HALF * EX
       EC = HALF * EC
       DO 10 ISP = 1,NSP
         VX(ISP) = HALF * VX(ISP)
         VC(ISP) = HALF * VC(ISP)
   10  CONTINUE
      END

