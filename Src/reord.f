      SUBROUTINE REORD( FCLUST, FSEQ, NM, NSM, ITR, MAXAUX, AUX )

C ********************************************************************
C Re-orders a clustered data array into a sequential one and viceversa
C Written by J.M.Soler. May'95.
C **** INPUT *********************************************************
C INTEGER NM(3)  : Number of Mesh divisions in each cell vector
C INTEGER NSM    : Number of Sub-divisions in each Mesh division
C INTEGER ITR    : TRanslation-direction switch
C                  ITR=+1 => From clustered to sequential
C                  ITR=-1 => From sequential to clustered
C INTEGER MAXAUX : MAXimum AUXiliary space
C **** INPUT OR OUTPUT (Depending on ITR ) ***************************
C REAL*4 FCLUST(NSM,NSM,NSM,NM1,NM2,NM3) : CLUSTered data
C REAL*4 FSEQ(NSM*NM1,NSM*NM2,NSM*NM3)   : SEQuential data
C **** AUXILIARY *****************************************************
C REAL*4 AUX(NM1*NM2*NSM**3) : Space that can be used freely outside
C ********************************************************************

      IMPLICIT NONE
      INTEGER
     .  NM(3), NSM, ITR, MAXAUX
      REAL
     .  FCLUST(*), FSEQ(*), AUX(*)
     
      INTEGER MAXNS3
      PARAMETER ( MAXNS3 = 8 )
      INTEGER 
     .  I, I0, I1, I2, I3, IS, IS1, IS2, IS3,
     .  J, J0, JS(MAXNS3), NAUX, NSM2, NSM3, NTM(3)
     
      CALL TIMER('REORD',1)

      NTM(1) = NM(1) * NSM
      NTM(2) = NM(2) * NSM
      NTM(3) = NM(3) * NSM
      NSM2 = NSM**2
      NSM3 = NSM**3
      NAUX = NM(1) * NM(2) * NSM3
      CALL CHKDIM( 'REORD', 'MAXAUX', MAXAUX, NAUX, 1 )
      CALL CHKDIM( 'REORD', 'MAXNS3', MAXNS3, NSM3, 1 )

      IS = 0
      DO IS3 = 0,NSM-1
        DO IS2 = 0,NSM-1
          DO IS1 = 0,NSM-1
            IS = IS + 1
            JS(IS) = 1 + IS1 + NTM(1)*IS2 + NTM(1)*NTM(2)*IS3
          ENDDO
        ENDDO
      ENDDO

      IF (ITR .GT. 0) THEN
        DO I3 = 0, NM(3)-1
          DO I2 = 0, NM(2)-1
            I0 = NSM3 * ( NM(1)*I2 + NM(1)*NM(2)*I3 )
            J0 = NTM(1)*NSM*I2
            DO IS = 1,NSM3
              I = I0 + IS
              J = J0 + JS(IS)
              DO I1 = 1,NM(1)
                AUX(J) = FCLUST(I)
                I = I + NSM3
                J = J + NSM
              ENDDO
            ENDDO
          ENDDO
          I = NM(1) * NM(2) * NSM3 * I3
          DO J = 1,NM(1)*NM(2)*NSM3
            FSEQ(I+J) = AUX(J)
          ENDDO
        ENDDO
      ELSE
        DO I3 = 0, NM(3)-1
          I = NM(1) * NM(2) * NSM3 * I3
          DO J = 1,NM(1)*NM(2)*NSM3
            AUX(J) = FSEQ(I+J)
          ENDDO
          DO I2 = 0, NM(2)-1
            I0 = NSM3 * ( NM(1)*I2 + NM(1)*NM(2)*I3 )
            J0 = NTM(1)*NSM*I2
            DO IS = 1,NSM3
              I = I0 + IS
              J = J0 + JS(IS)
              DO I1 = 1,NM(1)
                FCLUST(I) = AUX(J)
                I = I + NSM3
                J = J + NSM
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      CALL TIMER('REORD',2)
      END

