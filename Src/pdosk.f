      SUBROUTINE PDOSK( NSPIN, NUO, NO, MAXSPN, MAXUO, MAXNH, 
     .                  MAXO, NUMH, LISTHPTR, LISTH, H, S,
     .                  E1, E2, NHIST, SIGMA, 
     .                  XIJ, INDXUO, NK, KPOINT, EO, 
     .                  HAUX, SAUX, PSI, AUX, DTOT, DPR, NUOTOT )

C **********************************************************************
C Find the density of states projected onto the atomic orbitals
C     D_mu(E) = Sum(n,k,nu) C(mu,n,k) C(nu,n,k) S(mu,nu,k) Delta(E-E(n,k))
C where n run over all the bands between two given energies
C Written by J. Junquera and E. Artacho. Nov' 99
C ****  INPUT  *********************************************************
C INTEGER NSPIN             : Number of spin components (1 or 2)
C INTEGER NUO               : Number of atomic orbitals in the unit cell
C INTEGER NO                : Number of atomic orbitals in the supercell
C INTEGER MAXSPN            : Second dimension of eo and qo 
C                             (maximum number of differents spin polarizations)
C INTEGER MAXUO             : Maximum number of atomic orbitals in the unit cell
C INTEGER MAXNH             : Maximum number of orbitals interacting
C                             with any orbital
C INTEGER MAXO              : First dimension of eo
C INTEGER NUMH(NUO)         : Number of nonzero elements of each row
C                             of hamiltonian matrix
C INTEGER LISTHPTR(NUO)     : Pointer to each row (-1) of the
C                             hamiltonian matrix
C INTEGER LISTH(MAXNH)      : Nonzero hamiltonian-matrix element
C                             column indexes for each matrix row
C REAL*8  H(MAXNH,NSPIN)    : Hamiltonian in sparse format
C REAL*8  S(MAXNH)          : Overlap in sparse format
C REAL*8  E1, E2            : Energy range for density-matrix states
C                             (to find local density of states)
C                             Not used if e1 > e2
C INTEGER NHIST             : Number of the subdivisions of the histogram
C REAL*8  SIGMA             : Width of the gaussian to expand the eigenvectors
C REAL*8  XIJ(3,MAXNH)      : Vectors between orbital centers (sparse)
C                             (not used if only gamma point)
C INTEGER INDXUO(NO)        : Index of equivalent orbital in unit cell
C INTEGER NK                : Number of k points
C REAL*8  KPOINT(3,NK)      : k point vectors
C REAL*8  EO(MAXO,MAXSPN,NK): Eigenvalues
C INTEGER NUOTOT            : Total number of orbitals per unit cell
C ****  AUXILIARY  *****************************************************
C REAL*8  HAUX(2,NUO,NUO)   : Auxiliary space for the hamiltonian matrix
C REAL*8  SAUX(2,NUO,NUO)   : Auxiliary space for the overlap matrix
C REAL*8  PSI(2,NUO,NUO)    : Auxiliary space for the eigenvectors
C REAL*8  AUX(2*NUO*5)      : Extra auxiliary space
C ****  OUTPUT  ********************************************************
C REAL*8  DTOT(NHIST,2)   : Total density of states
C REAL*8  DPR(NHIST,NUO,2): Proyected density of states
C **********************************************************************

      IMPLICIT NONE

      INTEGER
     .  NSPIN, NUO, NO, MAXSPN, MAXUO, MAXNH, NK, 
     .  MAXO, NHIST, NUOTOT

      INTEGER
     .  NUMH(NUO), LISTHPTR(NUO), LISTH(MAXNH),
     .  INDXUO(NO)

      DOUBLE PRECISION
     .  H(MAXNH,NSPIN), S(MAXNH), E1, E2, SIGMA, 
     .  XIJ(3,MAXNH), KPOINT(3,NK), EO(MAXO,MAXSPN,NK),
     .  HAUX(2,NUOTOT,NUO), SAUX(2,NUOTOT,NUO), PSI(2,NUOTOT,NUO),
     .  AUX(2*NUOTOT*5), DTOT(NHIST,2), DPR(NHIST,MAXUO,2) 

C Internal variables ---------------------------------------------------
      INTEGER
     .  IK, ISPIN, IUO, JUO, IO, J, JO, IHIST, IBAND, JOTY, IND,
     .  IERROR

      DOUBLE PRECISION
     .  KXIJ, CKXIJ, SKXIJ, DELTA, ENER, DIFF, PI, PIPJ1, PIPJ2, 
     .  PIPJS1, PIPJS2, GAUSS, NORM

      EXTERNAL
     .  CDIAG

C Initialize some variables --------------------------------------------
      DELTA = (E2 - E1)/NHIST
      PI = 4.0D0 * ATAN(1.0D0)

C Solve eigenvalue problem for each k-point ----------------------------
      DO 100 ISPIN = 1, NSPIN

        DO 110 IK = 1, NK

C         Initialize auxiliar variables --------------------------------
          DO IUO = 1, NUO
            DO JUO = 1, NUOTOT
              SAUX(1,JUO,IUO) = 0.D0
              SAUX(2,JUO,IUO) = 0.D0
              HAUX(1,JUO,IUO) = 0.D0
              HAUX(2,JUO,IUO) = 0.D0
            ENDDO
          ENDDO

          DO 120 IO = 1, NUO
            DO 130 J = 1, NUMH(IO)
              IND = J + LISTHPTR(IO)
              JO = LISTH(IND)
              IUO= INDXUO(IO)
              JUO= INDXUO(JO)
C             Calculates the phases k*r_ij -----------------------------
              KXIJ = KPOINT(1,IK) * XIJ(1,IND) +
     .               KPOINT(2,IK) * XIJ(2,IND) +
     .               KPOINT(3,IK) * XIJ(3,IND) 
              CKXIJ = COS(KXIJ)
              SKXIJ = SIN(KXIJ)
C             Calculates the hamiltonian and the overlap in k space ----
C             H(k) = Sum(R) exp(i*k*R) * H(R) --------------------------
              SAUX(1,IUO,JUO) = SAUX(1,IUO,JUO) + S(IND) * CKXIJ
              SAUX(2,IUO,JUO) = SAUX(2,IUO,JUO) + S(IND) * SKXIJ
              HAUX(1,IUO,JUO) = HAUX(1,IUO,JUO) + H(IND,ISPIN) * CKXIJ
              HAUX(2,IUO,JUO) = HAUX(2,IUO,JUO) + H(IND,ISPIN) * SKXIJ
 130        ENDDO
 120      ENDDO

C         Symetrice the hamiltonian and overlap matrices ---------------
          DO 140 IUO = 1, NUO
            DO 150 JUO = 1, IUO-1
              SAUX(1,JUO,IUO) = 0.5D0 * ( SAUX(1,JUO,IUO) +
     .                                    SAUX(1,IUO,JUO) )
              SAUX(1,IUO,JUO) = SAUX(1,JUO,IUO)
              SAUX(2,JUO,IUO) = 0.5D0 * ( SAUX(2,JUO,IUO) -
     .                                    SAUX(2,IUO,JUO) )
              SAUX(2,IUO,JUO) = - SAUX(2,JUO,IUO)
              HAUX(1,JUO,IUO) = 0.5D0 * ( HAUX(1,JUO,IUO) +
     .                                    HAUX(1,IUO,JUO) )
              HAUX(1,IUO,JUO) = HAUX(1,JUO,IUO)
              HAUX(2,JUO,IUO) = 0.5D0 * ( HAUX(2,JUO,IUO) -
     .                                    HAUX(2,IUO,JUO) )
              HAUX(2,IUO,JUO) = - HAUX(2,JUO,IUO)
 150        ENDDO
            SAUX(2,IUO,IUO) = 0.D0
            HAUX(2,IUO,IUO) = 0.D0
 140      ENDDO


C         Diagonalize for each k point ---------------------------------
          CALL CDIAG( HAUX, NUOTOT, SAUX, NUOTOT, NUO,
     .                EO(1,ISPIN,IK), PSI, NUOTOT, NUO, IERROR )

C         Recalculate again the overlap matrix in k-space -------------
          DO IUO = 1, NUO
            DO JUO = 1, NUOTOT
              SAUX(1,JUO,IUO) = 0.D0
              SAUX(2,JUO,IUO) = 0.D0
            ENDDO
          ENDDO

          DO IO = 1, NUO
            DO  J = 1, NUMH(IO)
              IND = LISTHPTR(IO) + J
              JO = LISTH(IND)
              IUO= INDXUO(IO)
              JUO= INDXUO(JO)
C             Calculates the phases k*r_ij -----------------------------
              KXIJ = KPOINT(1,IK) * XIJ(1,IND) +
     .               KPOINT(2,IK) * XIJ(2,IND) +
     .               KPOINT(3,IK) * XIJ(3,IND) 
              CKXIJ = COS(KXIJ)
              SKXIJ = SIN(KXIJ)
C             Calculates the hamiltonian and the overlap in k space ----
C             H(k) = Sum(R) exp(i*k*R) * H(R) --------------------------
              SAUX(1,IUO,JUO) = SAUX(1,IUO,JUO) + S(IND) * CKXIJ
              SAUX(2,IUO,JUO) = SAUX(2,IUO,JUO) + S(IND) * SKXIJ
            ENDDO
          ENDDO

C         Symetrice the hamiltonian and overlap matrices ---------------
          DO IUO = 1, NUO
            DO JUO = 1, IUO-1
              SAUX(1,JUO,IUO) = 0.5D0 * ( SAUX(1,JUO,IUO) +
     .                                    SAUX(1,IUO,JUO) )
              SAUX(1,IUO,JUO) = SAUX(1,JUO,IUO)
              SAUX(2,JUO,IUO) = 0.5D0 * ( SAUX(2,JUO,IUO) -
     .                                    SAUX(2,IUO,JUO) )
              SAUX(2,IUO,JUO) = - SAUX(2,JUO,IUO)
            ENDDO
            SAUX(2,IUO,IUO) = 0.D0
          ENDDO

C         Loop over all the energy range -------------------------------
          DO 160 IHIST = 1, NHIST
            ENER = E1 + (IHIST - 1) * DELTA
            DO 170 IBAND = 1, NUO
              DIFF = (ENER - EO(IBAND,ISPIN,IK))**2 / (SIGMA ** 2)
              IF (DIFF .GT. 15.0D0) THEN
                GOTO 170
              ELSE
                GAUSS = ( EXP(-DIFF) )
                DTOT(IHIST,ISPIN) = DTOT(IHIST,ISPIN) + GAUSS
                DO 180 IUO = 1, NUO
C  Solo para los Juo que satisfagan el criterio del record...
                  DO 190 JUO = 1, NUO
                    PIPJ1 = PSI(1,IUO,IBAND) * PSI(1,JUO,IBAND) +
     .                      PSI(2,IUO,IBAND) * PSI(2,JUO,IBAND)
                    PIPJ2 = PSI(1,IUO,IBAND) * PSI(2,JUO,IBAND) -
     .                      PSI(2,IUO,IBAND) * PSI(1,JUO,IBAND)
                    PIPJS1= PIPJ1*SAUX(1,IUO,JUO)-PIPJ2*SAUX(2,IUO,JUO)
                    PIPJS2= PIPJ1*SAUX(2,IUO,JUO)+PIPJ2*SAUX(1,IUO,JUO)
                    DPR(IHIST,JUO,ISPIN)= DPR(IHIST,JUO,ISPIN) + 
     .                                     PIPJS1*GAUSS
C                   dpr(ihist,record) = dpr(ihist,record) + pipjs1*gauss
 190              ENDDO
 180            ENDDO
              ENDIF
 170        ENDDO

 160      ENDDO


 110    ENDDO

 100  ENDDO

      NORM = SIGMA * SQRT(PI) * DFLOAT(NK)

      DO IHIST = 1, NHIST
        DO ISPIN = 1, NSPIN
          DTOT(IHIST,ISPIN) = DTOT(IHIST,ISPIN) / NORM
          DO IUO = 1, MAXUO
            DPR(IHIST,IUO,ISPIN) = DPR(IHIST,IUO,ISPIN) /NORM
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END

 

