C $Id: xc.f,v 1.13 1999/03/04 16:21:11 jose Exp $

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

C Internal parameters
C NN     : order of the numerical derivatives: the number of radial 
C          points used is 2*NN+1
C MAXAUX : must be larger than the number of radial mesh points
C MAXSPN : maximum number of spin values (4 for non-collinear spin)
      INTEGER NN, MAXAUX, MAXSPN
      PARAMETER ( NN     =    5 )
      PARAMETER ( MAXAUX = 2000 )
      PARAMETER ( MAXSPN =    4 )

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
     .  AUX(MAXAUX), D(MAXSPN), DECDD(MAXSPN), DECDGD(3,MAXSPN),
     .  DEXDD(MAXSPN), DEXDGD(3,MAXSPN),
     .  DGDM(-NN:NN), DGIDFJ(-NN:NN), DRDM, DVOL, 
     .  EPSC, EPSX, F1, F2, GD(3,MAXSPN), PI
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



      SUBROUTINE CELLXC( FUNCTL, AUTHOR, IREL,
     .                   CELL, NMESH, NSPAN, MAXP, MTYPE, XMESH,
     .                   NSPIN, DENS,
     .                   EX, EC, DX, DC, VXC, STRESS, MAXAUX, AUX )

C *******************************************************************
C Finds total exchange-correlation energy and potential in a
C   periodic cell.
C This version implements the Local (spin) Density Approximation and
C   the Generalized-Gradient-Aproximation with the 'explicit mesh 
C   functional' method of White & Bird, PRB 50, 4954 (1994).
C Gradients are 'defined' by numerical derivatives, using 2*NN+1 mesh
C   points, where NN is a parameter defined below
C Coded by L.C.Balbas and J.M.Soler. December 1996.
C Stress added by J.L.Martins and J.M.Soler. August 1997.
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
C REAL*8  CELL(3,3)    : Unit cell vectors CELL(ixyz,ivector)
C INTEGER NMESH(3)     : Number of mesh divisions of each vector
C INTEGER NSPAN(3)     : Physical dimensions of arrays XMESH, DENS and
C                        VXC (or memory span between array elements)
C                        See usage section for more information
C INTEGER MAXP         : Physical dimension of XMESH, DENS, and VXC
C INTEGER MTYPE        : Mesh type:
C                        0 => Uniform mesh
C                        1 => Adaptive mesh, given in cartesian coord
C                        2 => Adaptive mesh, given in cell-vector coord
C REAL    XMESH(3,MAXP): Mesh point coordinates (not used if MTYPE=0)
C                        When MTYPE=2, cartesian coordinates are
C                        Xcart(ix,im) = Sum_iv(CELL(ix,iv)*XMESH(iv,ip))
C                        Notice single precision in this version
C INTEGER NSPIN        : NSPIN=1 => unpolarized; NSPIN=2 => polarized
C REAL    DENS(MAXP,NSPIN) : Total (NSPIN=1) or spin (NSPIN=2) electron
C                            density at mesh points, ordered as
C                            IP = I1+NSPAN(1)*((I2-1)+NSPAN(2)*(I3-1)),
C                            with I1=1,...,NMESH(1), etc
C                            Notice single precision in this version
C INTEGER MAXAUX       : Physical dimension of array AUX.
C                        If FUNCTL.EQ.'GGA' and MTYPE.NE.0, must be
C                        MAXAUX .GE. NSPAN(1)*NSPAN(2)*NSPAN(3)
C                        Not used otherwise
C ************************* OUTPUT **********************************
C REAL*8  EX              : Total exchange energy
C REAL*8  EC              : Total correlation energy
C REAL*8  DX              : IntegralOf( rho * (eps_x - v_x) )
C REAL*8  DC              : IntegralOf( rho * (eps_c - v_c) )
C REAL    VXC(MAXP,NSPIN) : (Spin) exch-corr potential
C                           Notice single precision in this version
C REAL*8  STRESS(3,3)     : xc contribution to the stress tensor,
C                           assuming constant density (not charge),
C                           i.e. r->r' => rho'(r') = rho(r)
C                           For plane-wave and grid (finite diff) basis
C                           sets, density rescaling gives an extra term
C                           (not included) (DX+DC-EX-EC)/cell_volume for
C                           the diagonal elements of stress. For other
C                           basis sets, the extra term is, in general:
C                           IntegralOf(v_xc * d_rho/d_strain) / cell_vol
C ************************* AUXILIARY *******************************
C REAL    AUX(MAXAUX)     : Not used unless 
C                           FUNCTL.EQ.'GGA' .AND. MTYPE.NE.0
C ************************ UNITS ************************************
C Distances in atomic units (Bohr).
C Densities in atomic units (electrons/Bohr**3)
C Energy unit depending of parameter EUNIT below
C Stress in EUNIT/Bohr**3
C ************************ USAGE ************************************
C Typical calls for different array dimensions:
C     PARAMETER ( MAXP = 1000000 )
C     DIMENSION NMESH(3), DENS(MAXP,2), VXC(MAXP,2)
C     Find CELL vectors
C     Find density at N1*N2*N3 mesh points (less than MAXP) and place 
C       them consecutively in array DENS
C     NMESH(1) = N1
C     NMESH(2) = N2
C     NMESH(3) = N3
C     CALL CELLXC( 'GGA', 'PBE', 0,
C    .              CELL, NMESH, NMESH, MAXP, 0, 0,
C    .              2, DENS,
C    .              EX, EC, DX, DC, VXC, STRESS, 0, 0 )
C Or alternatively:
C     PARAMETER ( M1=100, M2=100, M3=100 )
C     DIMENSION NMESH(3), NSPAN(3), DENS(M1,M2,M3,2), VXC(M1,M2,M3,2)
C     DATA NSPAN / M1, M2, M3 /
C     Find CELL vectors
C     Find DENS at N1*N2*N3 mesh points
C     NMESH(1) = N1
C     NMESH(2) = N2
C     NMESH(3) = N3
C     CALL CELLXC( 'GGA', 'PBE', 0,
C    .              CELL, NMESH, NSPAN, M1*M2*M3, 0, 0,
C    .              2, DENS,
C    .              EX, EC, DX, DC, VXC, STRESS, 0, 0 )
C ********* BEHAVIOUR ***********************************************
C - Notice that XMESH and AUX are not used if MTYPE=0. This means that
C   you do not even need to declare them in the calling program. See
C   usage section for calling examples.
C - If FNCTNL='LDA', DENS and VXC may be the same physical array.
C - Stops and prints a warning if arguments FUNCTL or AUTHOR are not
C   one of the allowed possibilities.
C - Since the exchange and correlation part is usually a small fraction
C   of a typical electronic structure calculation, this routine has
C   been coded with emphasis on simplicity and functionality, not in
C   eficiency. It is not designed to vectorize or parallelize well.
C ********* ROUTINES CALLED *****************************************
C GGAXC, LDAXC, RECLAT, TIMER, VOLCEL
C *******************************************************************

C Next line is nonstandard but may be suppressed
      IMPLICIT NONE

C Argument types and dimensions
      CHARACTER*(*)     FUNCTL, AUTHOR
      INTEGER           IREL, MAXAUX, MAXP, MTYPE,
     .                  NMESH(3), NSPAN(3), NSPIN
      DOUBLE PRECISION  CELL(3,3), DC, DX, EC, EX, STRESS(3,3)
C     If you change next line, please change also the argument
C     explanations in the header comments
*     DOUBLE PRECISION
      REAL
     .                  AUX(MAXAUX), DENS(MAXP,NSPIN),
     .                  VXC(MAXP,NSPIN), XMESH(3,MAXP)

C Fix the order of the numerical derivatives
C NN is the number of points used in each coordinate and direction,
C i.e. a total of 6*NN neighbour points is used to find the gradients
      INTEGER NN
      PARAMETER ( NN = 3 )

C Fix energy unit:  EUNIT=1.0 => Hartrees,
C                   EUNIT=0.5 => Rydbergs,
C                   EUNIT=0.03674903 => eV
      DOUBLE PRECISION EUNIT
      PARAMETER ( EUNIT = 0.5D0 )

C Fix switch to skip points with zero density
      LOGICAL SKIP0
      DOUBLE PRECISION DENMIN
      PARAMETER ( SKIP0  = .TRUE. )
      PARAMETER ( DENMIN = 1.D-15 )

C Local variables and arrays
      INTEGER MAXSPN
      PARAMETER ( MAXSPN = 4 )
      LOGICAL           GGA
      INTEGER           I1, I2, I3, IC, IN, IP, IS, IX,
     .                  J1, J2, J3, JN, JP(3,-NN:NN), JX, NP
      DOUBLE PRECISION  D(MAXSPN), DECDD(MAXSPN), DECDGD(3,MAXSPN),
     .                  DENTOT, DEXDD(MAXSPN), DEXDGD(3,MAXSPN),
     .                  DGDM(-NN:NN), DGIDFJ(3,3,-NN:NN),
     .                  DMDX(3,3), DVOL, DXDM(3,3),
     .                  EPSC, EPSX, F1, F2, GD(3,MAXSPN),
     .                  VOLCEL, VOLUME
      EXTERNAL          GGAXC, LDAXC, RECLAT, VOLCEL


C Start time counter (intended only for debugging and development)
      CALL TIMER( 'CELLXC', 1 )

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

C Find weights of numerical derivation from Lagrange interp. formula
      DO 15 IN = -NN,NN
        F1 = 1
        F2 = 1
        DO 10 JN = -NN,NN
          IF (JN.NE.IN .AND. JN.NE.0) F1 = F1 * (0  - JN)
          IF (JN.NE.IN)               F2 = F2 * (IN - JN)
   10   CONTINUE
        DGDM(IN) = F1 / F2
   15 CONTINUE
      DGDM(0) = 0

C Find total number of mesh points
      NP = NMESH(1) * NMESH(2) * NMESH(3)

C Find Jacobian matrix dx/dmesh for uniform mesh
      IF (MTYPE.EQ.0) THEN

C       Find mesh cell volume 
        DVOL = VOLCEL( CELL ) / NP

        IF (GGA) THEN

C         Find mesh unit vectors and reciprocal mesh vectors
          DO 25 IC = 1,3
            DO 20 IX = 1,3
              DXDM(IX,IC) = CELL(IX,IC) / NMESH(IC)
   20       CONTINUE
   25     CONTINUE
          CALL RECLAT( DXDM, DMDX, 0 )

C         Find the weights for the derivative d(gradF(i))/d(F(j)) of
C         the gradient at point i with respect to the value at point j
          DO 50 IN = -NN,NN
            DO 40 IC = 1,3
              DO 30 IX = 1,3
                DGIDFJ(IX,IC,IN) = DMDX(IX,IC) * DGDM(IN)
   30         CONTINUE
   40       CONTINUE
   50     CONTINUE

        ENDIF

      ENDIF

C Initialize output
      EX = 0
      EC = 0
      DX = 0
      DC = 0
      IF (GGA) THEN
        DO 70 IS = 1,NSPIN
          DO 66 I3 = 0,NMESH(3)-1
          DO 64 I2 = 0,NMESH(2)-1
          DO 62 I1 = 0,NMESH(1)-1
            IP = 1 + I1 + NSPAN(1) * I2 + NSPAN(1) * NSPAN(2) * I3
            VXC(IP,IS) = 0
   62     CONTINUE
   64     CONTINUE
   66     CONTINUE
   70   CONTINUE
      ENDIF
      DO 74 JX = 1,3
        DO 72 IX = 1,3
          STRESS(IX,JX) = 0.D0
   72   CONTINUE
   74 CONTINUE

C Loop on mesh points
      DO 230 I3 = 0,NMESH(3)-1
      DO 220 I2 = 0,NMESH(2)-1
      DO 210 I1 = 0,NMESH(1)-1

C       Find mesh index of this point
        IP = 1 + I1 + NSPAN(1) * I2 + NSPAN(1) * NSPAN(2) * I3

C       Skip point if density=0
        IF (SKIP0) THEN
          DENTOT = 0.D0
          DO 75 IS = 1,NSPIN
            DENTOT = DENTOT + ABS(DENS(IP,IS))
   75     CONTINUE
          IF (DENTOT .LT. DENMIN) THEN
            DO 76 IS = 1,NSPIN
              VXC(IP,IS) = 0
   76       CONTINUE
            GOTO 210
          ENDIF
        ENDIF

C       Find mesh indexes of neighbour points
        IF (GGA .OR. MTYPE.NE.0) THEN
          DO 80 IN = -NN,NN
            J1 = MOD( I1+IN+100*NMESH(1), NMESH(1) )
            J2 = MOD( I2+IN+100*NMESH(2), NMESH(2) )
            J3 = MOD( I3+IN+100*NMESH(3), NMESH(3) )
            JP(1,IN) = 1 + J1 + NSPAN(1) * I2 + NSPAN(1) * NSPAN(2) * I3
            JP(2,IN) = 1 + I1 + NSPAN(1) * J2 + NSPAN(1) * NSPAN(2) * I3
            JP(3,IN) = 1 + I1 + NSPAN(1) * I2 + NSPAN(1) * NSPAN(2) * J3
   80     CONTINUE
        ENDIF

C       Find Jacobian matrix dx/dmesh for adaptative mesh
        IF (MTYPE .NE. 0) THEN

C         Find dx/dmesh
          DO 110 IC = 1,3
            DO 100 IX = 1,3
              DXDM(IX,IC) = 0
              DO 90 IN = -NN,NN
                IF (MTYPE .EQ. 1) THEN
                  DXDM(IX,IC) = DXDM(IX,IC) +
     .                          XMESH(IX,JP(IC,IN)) * DGDM(IN)
                ELSE
                  DXDM(IX,IC) = DXDM(IX,IC) +
     .                   ( CELL(IX,1) * XMESH(1,JP(IC,IN)) +
     .                     CELL(IX,2) * XMESH(2,JP(IC,IN)) +
     .                     CELL(IX,3) * XMESH(3,JP(IC,IN)) ) * DGDM(IN)
                ENDIF
   90         CONTINUE
  100       CONTINUE
  110     CONTINUE

C         Find inverse of matrix dx/dmesh
          CALL RECLAT( DXDM, DMDX, 0 )

C         Find differential of volume = determinant of Jacobian matrix
          DVOL = VOLCEL( DXDM )
          IF (GGA) AUX(IP) = DVOL

C         Find the weights for the derivative d(gradF(i))/d(F(j)), of
C         the gradient at point i with respect to the value at point j
          IF (GGA) THEN
            DO 140 IN = -NN,NN
              DO 130 IC = 1,3
                DO 120 IX = 1,3
                  DGIDFJ(IX,IC,IN) = DMDX(IX,IC) * DGDM(IN)
  120           CONTINUE
  130         CONTINUE
  140       CONTINUE
          ENDIF

        ENDIF

C       Find density and gradient of density at this point
        DO 145 IS = 1,NSPIN
          D(IS) = DENS(IP,IS)
  145   CONTINUE
        IF (GGA) THEN
          DO 165 IS = 1,NSPIN
            DO 160 IX = 1,3
              GD(IX,IS) = 0
              DO 150 IN = -NN,NN
                GD(IX,IS) = GD(IX,IS) +
     .                      DGIDFJ(IX,1,IN) * DENS(JP(1,IN),IS) +
     .                      DGIDFJ(IX,2,IN) * DENS(JP(2,IN),IS) +
     .                      DGIDFJ(IX,3,IN) * DENS(JP(3,IN),IS)
  150         CONTINUE
  160       CONTINUE
  165     CONTINUE
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
        DO 170 IS = 1,MIN(NSPIN,2)
          EX = EX + DVOL * D(IS) * EPSX
          EC = EC + DVOL * D(IS) * EPSC
          DX = DX + DVOL * D(IS) * EPSX
          DC = DC + DVOL * D(IS) * EPSC
  170   CONTINUE
        DO 200 IS = 1,NSPIN
          DX = DX - DVOL * D(IS) * DEXDD(IS)
          DC = DC - DVOL * D(IS) * DECDD(IS)
          IF (GGA) THEN
            VXC(IP,IS) = VXC(IP,IS) + DVOL * ( DEXDD(IS) + DECDD(IS) )
            DO 190 IN = -NN,NN
              DO 180 IC = 1,3
                DO 175 IX = 1,3
                  DX = DX - DVOL * DENS(JP(IC,IN),IS) *
     .                      DEXDGD(IX,IS) * DGIDFJ(IX,IC,IN)
                  DC = DC - DVOL * DENS(JP(IC,IN),IS) *
     .                      DECDGD(IX,IS) * DGIDFJ(IX,IC,IN)
                  VXC(JP(IC,IN),IS) = VXC(JP(IC,IN),IS) + DVOL *
     .               (DEXDGD(IX,IS) + DECDGD(IX,IS)) * DGIDFJ(IX,IC,IN)
  175           CONTINUE
  180         CONTINUE
  190       CONTINUE
          ELSE
            VXC(IP,IS) = DEXDD(IS) + DECDD(IS)
          ENDIF
  200   CONTINUE

C       Add contribution to stress due to change in gradient of density
C       originated by the deformation of the mesh with strain
        IF (GGA) THEN
          DO 206 JX = 1,3
            DO 204 IX = 1,3
              DO 202 IS = 1,NSPIN
                STRESS(IX,JX) = STRESS(IX,JX) - DVOL * GD(IX,IS) *
     .                           ( DEXDGD(JX,IS) + DECDGD(JX,IS) )
  202         CONTINUE
  204       CONTINUE
  206     CONTINUE
        ENDIF

  210 CONTINUE
  220 CONTINUE
  230 CONTINUE

C Divide by volume element to obtain the potential (per electron)
      IF (GGA) THEN
        DO 270 IS = 1,NSPIN
          DO 260 I3 = 0,NMESH(3)-1
          DO 250 I2 = 0,NMESH(2)-1
          DO 240 I1 = 0,NMESH(1)-1
            IP = 1 + I1 + NSPAN(1) * I2 + NSPAN(1) * NSPAN(2) * I3
            IF (MTYPE .NE. 0) DVOL = AUX(IP)
            VXC(IP,IS) = VXC(IP,IS) / DVOL
  240     CONTINUE
  250     CONTINUE
  260     CONTINUE
  270   CONTINUE
      ENDIF

C Add contribution to stress from the change of volume with strain and
C divide by volume to get correct stress definition (dE/dStrain)/Vol
      VOLUME = VOLCEL( CELL )
      DO 274 JX = 1,3
        STRESS(JX,JX) = STRESS(JX,JX) + EX + EC
        DO 272 IX = 1,3
          STRESS(IX,JX) = STRESS(IX,JX) / VOLUME
  272   CONTINUE
  274 CONTINUE
      
C Divide by energy unit
      EX = EX / EUNIT
      EC = EC / EUNIT
      DX = DX / EUNIT
      DC = DC / EUNIT
      DO 310 IS = 1,NSPIN
        DO 300 I3 = 0,NMESH(3)-1
        DO 290 I2 = 0,NMESH(2)-1
        DO 280 I1 = 0,NMESH(1)-1
          IP = 1 + I1 + NSPAN(1) * I2 + NSPAN(1) * NSPAN(2) * I3
          VXC(IP,IS) = VXC(IP,IS) / EUNIT
  280   CONTINUE
  290   CONTINUE
  300   CONTINUE
  310 CONTINUE
      DO 330 JX = 1,3
        DO 320 IX = 1,3
          STRESS(IX,JX) = STRESS(IX,JX) / EUNIT
  320   CONTINUE
  330 CONTINUE

C Stop time counter
      CALL TIMER( 'CELLXC', 2 )

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
         D1 = MAX(DS(1),0.D0)
         D2 = MAX(DS(2),0.D0)
         D = D1 + D2
         IF (D .LE. ZERO) THEN
           EX = ZERO
           VX(1) = ZERO
           VX(2) = ZERO
           RETURN
         ENDIF
         Z = (D1 - D2) / D
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
       VXP = -(3*ALP/(2*PI*A0*RS))
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
      DOUBLE PRECISION  D(NSPIN), DECDD(NSPIN), DECDGD(3,NSPIN),
     .                  DEXDD(NSPIN), DEXDGD(3,NSPIN),
     .                  EPSC, EPSX, GD(3,NSPIN)

      IF (NSPIN .GT. 2) STOP 'GGAXC: Not prepared for this NSPIN'

      IF (AUTHOR.EQ.'PBE' .OR. AUTHOR.EQ.'pbe') THEN
        CALL PBEXC( IREL, NSPIN, D, GD,
     .              EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
      ELSE
        WRITE(6,*) 'GGAXC: Unknown author ', AUTHOR
        STOP
      ENDIF
      END




      SUBROUTINE LDAXC( AUTHOR, IREL, NSPIN, D, EPSX, EPSC, VX, VC )

C ******************************************************************
C Finds the exchange and correlation energies and potentials, in the
C Local (spin) Density Approximation.
C Written by L.C.Balbas and J.M.Soler, Dec'96.
C Non-collinear spin added by J.M.Soler, May'98
C *********** INPUT ************************************************
C CHARACTER*(*) AUTHOR : Parametrization desired:
C     'CA' or 'PZ' => LSD Perdew & Zunger, PRB 23, 5075 (1981)
C           'PW92' => LSD Perdew & Wang, PRB, 45, 13244 (1992)
C                     Uppercase is optional
C INTEGER IREL     : Relativistic exchange? (0=>no, 1=>yes)
C INTEGER NSPIN    : NSPIN=1 => unpolarized; NSPIN=2 => polarized;
C                    NSPIN=4 => non-collinear polarization
C REAL*8  D(NSPIN) : Local (spin) density. For non-collinear
C                    polarization, the density matrix is given by:
C                    D(1)=D11, D(2)=D22, D(3)=Real(D12), D(4)=Im(D12)
C *********** OUTPUT ***********************************************
C REAL*8 EPSX, EPSC : Exchange and correlation energy densities
C REAL*8 VX(NSPIN), VC(NSPIN) : Exchange and correlation potentials,
C                               defined as dExc/dD(ispin)
C *********** UNITS ************************************************
C Lengths in Bohr, energies in Hartrees
C ******************************************************************

      IMPLICIT          NONE
      CHARACTER*(*)     AUTHOR
      INTEGER           IREL, NSPIN
      DOUBLE PRECISION  D(NSPIN), EPSC, EPSX, VX(NSPIN), VC(NSPIN)

      INTEGER           IS, NS
      DOUBLE PRECISION  DD(2), DPOL, DTOT, TINY, VCD(2), VPOL, VXD(2)

      PARAMETER ( TINY = 1.D-12 )

      IF (NSPIN .EQ. 4) THEN
C       Find eigenvalues of density matrix (up and down densities
C       along the spin direction)
C       Note: D(1)=D11, D(2)=D22, D(3)=Real(D12), D(4)=Im(D12)
        NS = 2
        DTOT = D(1) + D(2)
        DPOL = SQRT( (D(1)-D(2))**2 + 4.D0*(D(3)**2+D(4)**2) )
        DD(1) = 0.5D0 * ( DTOT + DPOL )
        DD(2) = 0.5D0 * ( DTOT - DPOL )
      ELSE
        NS = NSPIN
        DO 10 IS = 1,NSPIN
          DD(IS) = D(IS)
   10   CONTINUE
      ENDIF

      IF ( AUTHOR.EQ.'CA' .OR. AUTHOR.EQ.'ca' .OR.
     .     AUTHOR.EQ.'PZ' .OR. AUTHOR.EQ.'pz') THEN
        CALL PZXC( IREL, NS, DD, EPSX, EPSC, VXD, VCD )
      ELSEIF ( AUTHOR.EQ.'PW92' .OR. AUTHOR.EQ.'pw92' ) THEN
        CALL PW92XC( IREL, NS, DD, EPSX, EPSC, VXD, VCD )
      ELSE
        WRITE(6,*) 'LDAXC: Unknown author ', AUTHOR
        STOP
      ENDIF

      IF (NSPIN .EQ. 4) THEN
C       Find dE/dD(ispin) = dE/dDup * dDup/dD(ispin) +
C                           dE/dDdown * dDown/dD(ispin)
        VPOL  = (VXD(1)-VXD(2)) * (D(1)-D(2)) / (DPOL+TINY)
        VX(1) = 0.5D0 * ( VXD(1) + VXD(2) + VPOL )
        VX(2) = 0.5D0 * ( VXD(1) + VXD(2) - VPOL )
        VX(3) = (VXD(1)-VXD(2)) * D(3) / (DPOL+TINY)
        VX(4) = (VXD(1)-VXD(2)) * D(4) / (DPOL+TINY)
        VPOL  = (VCD(1)-VCD(2)) * (D(1)-D(2)) / (DPOL+TINY)
        VC(1) = 0.5D0 * ( VCD(1) + VCD(2) + VPOL )
        VC(2) = 0.5D0 * ( VCD(1) + VCD(2) - VPOL )
        VC(3) = (VCD(1)-VCD(2)) * D(3) / (DPOL+TINY)
        VC(4) = (VCD(1)-VCD(2)) * D(4) / (DPOL+TINY)
      ELSE
        DO 20 IS = 1,NSPIN
          VX(IS) = VXD(IS)
          VC(IS) = VCD(IS)
   20   CONTINUE
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
     .  DS(2), DSDD, DSDGD, DT, DTDD, DTDGD, DZDD(2), 
     .  EC, ECUNIF, EX, EXUNIF,
     .  F, F1, F2, F3, F4, FC, FX, FOUTHD,
     .  GAMMA, GD(3,2), GDM(2), GDMIN, GDMS, GDMT, GDS, GDT(3),
     .  H, HALF, KAPPA, KF, KFS, KS, MU, PHI, PI, RS, S,
     .  T, THD, THRHLF, TWO, TWOTHD, VCUNIF(2), VXUNIF(2), ZETA

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
      DRSDD = - (THD * RS / DT)
      DKFDD =   THD * KF / DT
      DKSDD = HALF * KS * DKFDD / KF
      DZDD(1) =   1 / DT - ZETA / DT
      DZDD(2) = - (1 / DT) - ZETA / DT
      DPDZ = HALF * TWOTHD * ( 1/(1+ZETA)**THD - 1/(1-ZETA)**THD )
      DO 40 IS = 1,2
        DECUDD = ( VCUNIF(IS) - ECUNIF ) / DT
        DPDD = DPDZ * DZDD(IS)
        DTDD = (- T) * ( DPDD/PHI + DKSDD/KS + 1/DT )
        DF1DD = F1 * ( DECUDD/ECUNIF - 3*DPDD/PHI )
        DF2DD = (- F2) * DF1DD
        DADD = (- A) * DF2DD / (F2-1)
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
        DS(IS)   = MAX( DENMIN, 2 * D(IS) )
        GDMS = MAX( GDMIN, 2 * GDM(IS) )
        KFS = (3 * PI**2 * DS(IS))**THD
        S = GDMS / (2 * KFS * DS(IS))
        F1 = 1 + MU * S**2 / KAPPA
        F = 1 + KAPPA - KAPPA / F1
c
c       Note nspin=1 in call to exchng...
c
        CALL EXCHNG( IREL, 1, DS(IS), EXUNIF, VXUNIF(IS) )
        FX = FX + DS(IS) * EXUNIF * F

        DKFDD = THD * KFS / DS(IS)
        DSDD = S * ( -(DKFDD/KFS) - 1/DS(IS) )
        DF1DD = 2 * (F1-1) * DSDD / S
        DFDD = KAPPA * DF1DD / F1**2
        DFXDD(IS) = VXUNIF(IS) * F + DS(IS) * EXUNIF * DFDD

        DO 50 IX = 1,3
          GDS = 2 * GD(IX,IS)
          DSDGD = (S / GDMS) * GDS / GDMS
          DF1DGD = 2 * MU * S * DSDGD / KAPPA
          DFDGD = KAPPA * DF1DGD / F1**2
          DFXDGD(IX,IS) = DS(IS) * EXUNIF * DFDGD
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
      DATA P      / 1.00d0,     1.00d0,     1.00d0     /
      DATA A      / 0.031091d0, 0.015545d0, 0.016887d0 /
      DATA ALPHA1 / 0.21370d0,  0.20548d0,  0.11125d0  /
      DATA BETA   / 7.5957d0,  14.1189d0,  10.357d0,
     .              3.5876d0,   6.1977d0,   3.6231d0,
     .              1.6382d0,   3.3662d0,   0.88026d0,
     .              0.49294d0,  0.62517d0,  0.49671d0 /

C Find rs and zeta
      PI = 4 * ATAN(1.D0)
      IF (NSPIN .EQ. 1) THEN
        DTOT = MAX( DENMIN, DENS(1) )
        ZETA = 0
        RS = ( 3 / (4*PI*DTOT) )**THD
C       Find derivatives dRs/dDens and dZeta/dDens
        DRSDD = (- RS) / DTOT / 3
        DZDD(1) = 0
      ELSE
        DTOT = MAX( DENMIN, DENS(1)+DENS(2) )
        ZETA = ( DENS(1) - DENS(2) ) / DTOT
        RS = ( 3 / (4*PI*DTOT) )**THD
        DRSDD = (- RS) / DTOT / 3
        DZDD(1) =   1 / DTOT - ZETA / DTOT
        DZDD(2) = - (1 / DTOT) - (ZETA / DTOT)
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
        DCDRS = - ( (C-1) * DBDRS / B )
        G(IG) = (- 2) * A(IG) * ( 1 + ALPHA1(IG)*RS ) * LOG(C)
        DGDRS(IG) = (- 2) *A(IG) * ( ALPHA1(IG) * LOG(C) +
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
      DECDZ = (- G(2)) / FPP0 * ( DFDZ*(1-ZETA**4) - F*4*ZETA**3 ) +
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
       PARAMETER ( CXP  = (- 3.D0) * ALP / (PI*A0) )
       PARAMETER ( CXF  = 1.25992104989487319D0 )

C      Find density and polarization
       IF (NSP .EQ. 2) THEN
         D1 = MAX(DS(1),0.D0)
         D2 = MAX(DS(2),0.D0)
         D = D1 + D2
         IF (D .LE. ZERO) THEN
           EX = ZERO
           EC = ZERO
           VX(1) = ZERO
           VX(2) = ZERO
           VC(1) = ZERO
           VC(2) = ZERO
           RETURN
         ENDIF
         Z = (D1 - D2) / D
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
         ECP = -(C2846/BE)
         VCP = ECP*TE/BE
         TE = ONE+CON8*SQRS+CON9*RS
         BE = ONE+C1P398*SQRS+C2611*RS
         ECF = -(C1686/BE)
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






