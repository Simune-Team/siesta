! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module m_forhar
      use precision,    only : dp, grid_p
      use alloc,        only : re_alloc, de_alloc
      use mesh
      use moreMeshSubs, only : setMeshDistr, distMeshData
      implicit none

      public :: forhar
      private

      CONTAINS
      subroutine forhar( NTPL, NSPIN, NML, NTML, NPCC, CELL, RHOATM,
     &                   RHOPCC, VNA, DRHOOUT, VHARRIS1, VHARRIS2 )

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

      INTEGER             :: NTPL, NML(3), NTML(3)
      INTEGER, INTENT(IN) :: NSPIN, NPCC
 
      REAL(dp),                INTENT(IN) :: CELL(3,3)
      REAL(grid_p),            INTENT(IN) :: VNA(NTPL), RHOATM(NTPL),
     &                                       RHOPCC(NTPL)
      REAL(grid_p),         INTENT(INOUT) :: DRHOOUT(NTPL,NSPIN)
      REAL(grid_p), TARGET, INTENT(INOUT) :: VHARRIS1(NTPL,NSPIN)
      REAL(grid_p),         INTENT(INOUT) :: VHARRIS2(NTPL)
      
      EXTERNAL REORD, CELLXC

! AG: Note:  REAL*4 variables are really REAL(kind=grid_p)
!
C ***** INPUT **********************************************************
C INTEGER NTPL                 : Number of Mesh Total Points in unit cell 
C                                (including subpoints) locally. 
C INTEGER NSPIN                : Spin polarizations
C INTEGER NTM(3)               : Number of mesh divisions of each cell
C                                vector, including subgrid
C INTEGER NPCC                 : Partial core corrections? (0=no,1=yes)
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
      INTEGER                :: IP, ISPIN, ISPIN2, NMPL
      REAL(dp)               :: EX, EC, DEX, DEC, STRESSL(3,3)
      real(grid_p)           :: aux3(3,1)   !! dummy arrays for cellxc

      real(grid_p),  pointer :: drhoin(:,:), drhoin_par(:,:),
     &                          dvxcdn(:,:,:), dvxcdn_par(:,:,:),
     &                          vharris1_par(:,:), fsrc(:), fdst(:)
      INTEGER                :: idistr, odistr, ntpl_1, ntpl_3

      nullify( drhoin, dvxcdn )
      call re_alloc( drhoin, 1, ntpl, 1, nspin, 'drhoin', 'forhar' )
      call re_alloc( dvxcdn, 1, ntpl, 1, nspin, 1, nspin,
     &               'dvxcdn', 'forhar' )

C ----------------------------------------------------------------------
C Initialize some variables
C ----------------------------------------------------------------------
      VHARRIS1(:,:) = 0.0_grid_p
      VHARRIS2(:)   = 0.0_grid_p
      DRHOIN(:,:)   = 0.0_grid_p
      DVXCDN(:,:,:) = 0.0_grid_p
      STRESSL(:,:)  = 0.0_dp

C ----------------------------------------------------------------------
C Compute exchange-correlation energy and potential and
C their derivatives respect the input charge, that is, Harris charge
C or the sum of atomic charges.
C ----------------------------------------------------------------------

      DO ISPIN = 1, NSPIN
        DRHOIN(1:NTPL,ISPIN) =  RHOATM(1:NTPL)/NSPIN 
        IF (NPCC .EQ. 1) 
     .    DRHOIN(1:NTPL,ISPIN) = DRHOIN(1:NTPL,ISPIN) + 
     .                           RHOPCC(1:NTPL)/NSPIN
      ENDDO


      ntpl_1 = ntpl
#ifdef MPI
      idistr = 1
      odistr = 3
      call setMeshDistr( odistr, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
      ntpl_3 = ntpl

      nullify( drhoin_par, vharris1_par, dvxcdn_par )
      call re_alloc( drhoin_par, 1, ntpl_3, 1, nspin,
     &               'drhoin_par', 'forhar' )
      call re_alloc( vharris1_par, 1, ntpl_3, 1, nspin,
     &               'vharris1_par', 'forhar' )
      call re_alloc( dvxcdn_par, 1, ntpl_3, 1, nspin, 1, nspin,
     &               'dvxcdn_par', 'forhar' )

      DO ISPIN = 1, NSPIN
        fsrc => drhoin(:,ispin)
        fdst => drhoin_par(:,ispin)
        call distMeshData( idistr, fsrc, odistr, fdst, 0 )
      ENDDO

      CALL CELLXC( 0, 1, CELL, NTML, NTML, NTPL, 0, AUX3, NSPIN,
     &             DRHOIN_PAR, EX, EC, DEX, DEC, VHARRIS1_PAR,
     &             DVXCDN_PAR, STRESSL )

      DO ISPIN = 1, NSPIN
        fsrc => VHARRIS1_PAR(:,ISPIN)
        fdst => VHARRIS1(:,ispin)
        call distMeshData( odistr, fsrc, idistr, fdst, 0 )
        DO ISPIN2 = 1, NSPIN
          fsrc => DVXCDN_PAR(:,ISPIN,ISPIN2)
          fdst => DVXCDN(:,ISPIN,ISPIN2)
          call distMeshData( odistr, fsrc, idistr, fdst, 0 )
        ENDDO
      ENDDO

      call de_alloc( dvxcdn_par,   'dvxcdn_par',   'forhar' )
      call de_alloc( vharris1_par, 'vharris1_par', 'forhar' )
      call de_alloc( drhoin_par,   'drhoin_par',   'forhar' )

#ifdef MPI
      call setMeshDistr( idistr, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
      DO ISPIN = 1, NSPIN
        IF( NPCC .EQ. 1) 
     &    DRHOIN(1:NTPL,ISPIN) = DRHOIN(1:NTPL,ISPIN) - 
     &                           RHOPCC(1:NTPL)/NSPIN
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
     &                     DVXCDN(IP,ISPIN2,ISPIN) * DRHOOUT(IP,ISPIN2)
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

      call de_alloc( dvxcdn, 'dvxcdn', 'forhar')
      call de_alloc( drhoin, 'drhoin', 'forhar')
      end subroutine forhar

      end module m_forhar
