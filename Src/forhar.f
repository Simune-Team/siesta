
      SUBROUTINE FORHAR( NTPL, NSPIN, NML, NTML, NTM, NPCC,
     .                   XCFUNC, XCAUTH, CELL, RHOATM, RHOPCC,
     .                   VNA, DRHOOUT, 
     .                   VHARRIS1, VHARRIS2, Node, Nodes )

C **********************************************************************
C Build the potentials needed for computing Harris forces:
C (V_NA + V_Hartree(DeltaRho_in) - DV_xc(Rho_in)/Dn * (Rho_out-Rho_in))
C (V_NA + V_Hartree(DeltaRho_in) + V_xc(Rho_in)
C **** BEHAVIOUR *******************************************************
C   In the first SCF step, V_Hartree(DeltaRho_in) is zero, because
C in that case, Rho_SCF(r) = Rho_atm(r) and therefore, DeltaRho(r) = 0
C This calculation will be skipped.
C   If Harris + Spin polarized in the first SCF step, then Vharris2 will
C multiply to D Rho(Harris)/D R inside dfscf, and the change of 
C the harris density respect the displacement
C on one atom will not depend on spin. We add in Vharris2 the
C contributions of both spins.
C Coded by J. Junquera 09/00
C **********************************************************************

      USE MESH

      IMPLICIT NONE

      CHARACTER XCFUNC*10, XCAUTH*10

      INTEGER, INTENT(IN) :: NTPL, NSPIN, 
     .                       NPCC, Node, Nodes
      INTEGER, DIMENSION(3), INTENT(IN) :: NML, NTML, NTM
 
      REAL*8, DIMENSION(3,3), INTENT(IN) :: CELL(3,3) 
      REAL, DIMENSION(NTPL), INTENT(IN) :: VNA, RHOATM, RHOPCC
      REAL, DIMENSION(NTPL,NSPIN), INTENT(INOUT) :: DRHOOUT
      REAL, DIMENSION(NTPL,NSPIN), INTENT(INOUT) :: VHARRIS1
      REAL, DIMENSION(NTPL), INTENT(INOUT) :: VHARRIS2
      
      EXTERNAL REORD, CELLXC
       
C ***** INPUT **********************************************************
C INTEGER NTPL                 : Number of Mesh Total Points in unit cell 
C                                (including subpoints) locally. 
C INTEGER NSPIN                : Spin polarizations
C INTEGER NTM(3)               : Number of mesh divisions of each cell
C                                vector, including subgrid
C INTEGER NPCC                 : Partial core corrections? (0=no,1=yes)
C CHARACTER*10 XCFUNC          : Name of xc functional
C CHARACTER*10 XCAUTH          : Initials of xc version authors
C REAL*8 CELL(3,3)             : Cell vectors
C REAL*4 RHOATM(NTPL)          : Harris density at mesh points
C REAL*4 RHOPCC(NTPL)          : Partial-core-correction density for xc
C REAL*4 VNA(NTPL)             : Sum of neutral atoms potentials
C REAL*4 DRHOOUT(NTPL,NSPIN)   : Charge density at the mesh points
C                                in current step.
C                                The charge density that enters in forhar
C                                is Drho_out-Rhoatm.
C ***** OUTPUT *********************************************************
C REAL*4 VHARRIS1(NTPL,NSPIN)  : Vna + V_Hartree(DeltaRho_in) + V_xc(Rho_in)
C REAL*4 VHARRIS2(NTPL)        : Vna + V_Hartree(Rho_in) +
C                              - DV_xc(Rho_in)/DRho_in * (Rho_out-Rho_in)
C                                If Harris forces are computed in the 
C                                first SCF step, it does not depend on spin
C ***** INTERNAL VARIABLES *********************************************
C REAL*4 DVXDN(NTPL,NSPIN,NSPIN): Derivative of exchange-correlation
C                                potential respect the charge density
C **********************************************************************

C ----------------------------------------------------------------------
C Internal variables and arrays
C ----------------------------------------------------------------------
      INTEGER IP, ISPIN, ISPIN2
      REAL*8 DUSCF, DSTRES(3,3), RHOTOT, EX, EC, DEX, DEC, STRESSL(3,3)
      real aux1(1), aux3(3,1)   !! dummy arrays for cellxc
      REAL, DIMENSION(:,:), ALLOCATABLE :: DRHOIN
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: DVXCDN

      ALLOCATE(DRHOIN(NTPL,NSPIN))
      CALL MEMORY('A','S',ntpl*nspin,'forhar')
      ALLOCATE(DVXCDN(NTPL,NSPIN,NSPIN))
      CALL MEMORY('A','S',ntpl,'forhar')

C ----------------------------------------------------------------------
C Initialize some variables
C ----------------------------------------------------------------------
      VHARRIS1(:,:) = 0.0     
      VHARRIS2(:)   = 0.0     
      DRHOIN(:,:)   = 0.0
      DVXCDN(:,:,:) = 0.0

C ----------------------------------------------------------------------
C Compute exchange-correlation energy and potential and
C their derivatives respect the input charge, that is, Harris charge
C or the sum of atomic charges.
C ----------------------------------------------------------------------

      DO ISPIN = 1, NSPIN
        DRHOIN(1:NTPL,ISPIN) =  RHOATM(1:NTPL)/NSPIN 
        IF( NPCC .EQ. 1) 
     .    DRHOIN(1:NTPL,ISPIN) = DRHOIN(1:NTPL,ISPIN) + 
     .                           RHOPCC(1:NTPL)/NSPIN
      ENDDO

      DO ISPIN = 1, NSPIN
        CALL REORD(DRHOIN(1,ISPIN),DRHOIN(1,ISPIN),NML,NSM,+1)
      ENDDO

      CALL CELLXC( XCFUNC, XCAUTH, 0, 1, CELL, NTML, NTML, NTPL, 0,
     .             AUX3, NSPIN, DRHOIN, EX, EC, DEX, DEC, VHARRIS1,
     .             NTPL, DVXCDN, STRESSL, 
     .             1, AUX1, NTM, NSM, Node, Nodes )

      DO ISPIN = 1, NSPIN
        CALL REORD(DRHOIN(1,ISPIN),DRHOIN(1,ISPIN),NML,NSM,-1)
        CALL REORD(VHARRIS1(1,ISPIN),VHARRIS1(1,ISPIN),NML,NSM,-1)
        DO ISPIN2 = 1, NSPIN
          CALL REORD(DVXCDN(1,ISPIN,ISPIN2),DVXCDN(1,ISPIN,ISPIN2),
     .               NML,NSM,-1)
        ENDDO
      ENDDO

      DO ISPIN = 1, NSPIN
        IF( NPCC .EQ. 1) 
     .    DRHOIN(1:NTPL,ISPIN) = DRHOIN(1:NTPL,ISPIN) - 
     .                           RHOPCC(1:NTPL)/NSPIN
        DO IP = 1, NTPL
          VHARRIS1(IP,ISPIN) = VHARRIS1(IP,ISPIN) + VNA(IP)
        ENDDO
      ENDDO

C ----------------------------------------------------------------------
C Compute the product DV_xc(Rho_in)/DRho_in * (Rho_out - Rho_in).
C Since the charge that enters into forhar is DRHOOUT = Rho_out-Rhoatm
C no extra transformation on the charge density is needed.
C ----------------------------------------------------------------------

      DO ISPIN = 1, NSPIN
        DO ISPIN2 = 1, NSPIN
          DO IP = 1, NTPL
            VHARRIS2(IP) = VHARRIS2(IP) + 
     .                     DVXCDN(IP,ISPIN2,ISPIN) * DRHOOUT(IP,ISPIN2)
          ENDDO
        ENDDO
      ENDDO

C ----------------------------------------------------------------------
C Since V_Hartree(DeltaRho_in) = 0.0, we only add to vharris2 the neutral
C atom potential
C ----------------------------------------------------------------------
      DO IP = 1, NTPL
        VHARRIS2(IP) = VNA(IP) - VHARRIS2(IP)
      ENDDO

      IF (ALLOCATED(DRHOIN)) THEN
        CALL MEMORY('D','S',SIZE(DRHOIN),'forhar')
        DEALLOCATE(DRHOIN)
      ENDIF
      IF (ALLOCATED(DVXCDN)) THEN
        CALL MEMORY('D','S',SIZE(DVXCDN),'forhar')
        DEALLOCATE(DVXCDN)
      ENDIF

      END SUBROUTINE FORHAR



