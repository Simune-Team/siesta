! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      SUBROUTINE STM( NA, NO, NO_U, MAXNA, nspin_blocks, non_coll,
     .                ISA, IPHORB, INDXUO, LASTO, XA, CELL, UCELL,
     .                wf_unit, NK, gamma_wfsx,
     .                ZREF, ZMIN, ZMAX, NPX, NPY, NPZ, NSCX, NSCY,
     .                V0, EMAX, EMIN,
     .                ARMUNI, IUNITCD, RMAXO )

C **********************************************************************
C Simulate STM images in the Tersoff-Hamann approximation, by
C extrapolating the wavefunctions into vacuum
C
C Coded by P. Ordejon and N. Lorente,  November 2004
C
C Modified by N. Lorente, August 2005
C **********************************************************************

      use precision, only: dp, sp
      USE ATMFUNCS
      USE FDF
      USE CHEMICAL

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::
     .  NA, NO, NO_U, NPX, NPY, NPZ, IUNITCD,
     .  nspin_blocks, MAXNA, NK,
     .  ISA(NA), IPHORB(NO), INDXUO(NO), LASTO(0:NA),
     .  NSCX, NSCY
      integer, intent(in) :: wf_unit
      logical, intent(in) :: non_coll, gamma_wfsx

      REAL(DP), INTENT(IN) ::
     .  ZMIN, ZMAX, ZREF, 
     .  ARMUNI, RMAXO, V0, EMAX, EMIN

      REAL(DP), INTENT(IN) :: CELL(3,3)
      REAL(DP) :: UCELL(3,3), VOLCEL, XA(3,NA)

      EXTERNAL :: VOLCEL
      
C ****** INPUT *********************************************************
C INTEGER NA               : Total number of atoms in Supercell
C INTEGER NO               : Total number of orbitals in Supercell
C INTEGER NO_U              : Total number of orbitals in Unit Cell
C INTEGER MAXNA            : Maximum number of neighbours of any atom
C INTEGER ISA(NA)          : Species index of each atom
C INTEGER IPHORB(NO)       : Orital index of each orbital in its atom
C INTEGER INDXUO(NO)       : Equivalent orbital in unit cell
C INTEGER LASTO(0:NA)      : Last orbital of each atom in array iphorb
C REAL*8  XA(3,NA)         : Atomic positions in cartesian coordinates
C                            (in bohr), can be modified in the cube format
C REAL*8  CELL(3,3)        : Supercell vectors CELL(IXYZ,IVECT)
C                            (in bohr)
C REAL*8  UCELL(3,3)       : Unit cell vectors CELL(IXYZ,IVECT)
C                            (in bohr)
C INTEGER NK               : Number of k-points
c REAL*8 ZREF              : Position of reference plane for wf. estrapol.
C REAL*8  ZMIN, ZMAX       : Limits of the z-direction for the STM scan
C INTEGER NPX,NPY,NPZ      : Number of points along x and y and z
C INTEGER NSCX, NSCY       : Number of cells in x and y direction to plot
C                            in cube file
C REAL*8  V0               : Value of the potential at the vacuum region in eV
C REAL*8  EMAX             : Maximum value for the energy window for STM in eV
C REAL*8  EMIN             : Minimum value for the energy window for STM in eV
C REAL*8  ARMUNI           : Conversion factor for the charge density
C INTEGER IUNITCD          : Unit of the charge density
C REAL*8  RMAXO            : Maximum range of basis orbitals
C **********************************************************************

      INTEGER, DIMENSION(:), ALLOCATABLE ::  JNA
      REAL(DP), DIMENSION(:), ALLOCATABLE :: R2IJ
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: XIJ

      INTEGER
     .  IA, ISEL, NNA, I, J, IN, IAT1, IO, IUO, IAVEC1, 
     .  IS1, IPHI1, NX, NY, NZ, IWF, IK, ISPIN, UNITRE1,
     .  IX, IY, IZ, NSX, NSY, NAU

      REAL(DP)
     .  DOT, RMAX, XPO(3), RMAX2, XVEC1(3),
     .  PHIMU, GRPHIMU(3),
     .  PHASE, SI, CO, ENER, PMIKR, SIMIKR, COMIKR, USAVE, VC, VU

      real(dp) :: total_weight, k(3)
      
      REAL(DP), ALLOCATABLE :: RHO(:,:,:)
      REAL(SP), ALLOCATABLE :: wf_single(:,:)
      COMPLEX(DP), ALLOCATABLE :: wf(:,:)
      REAL(DP), ALLOCATABLE :: wk(:)

      COMPLEX(DP)  CWAVE, EXPPHI, EXMIKR
      COMPLEX(DP), ALLOCATABLE :: CW(:,:), CWE(:,:,:)
 
      LOGICAL FIRST
      integer :: idummy, number_of_wfns

      CHARACTER  SNAME*40, FNAME*60, stm_label*60

      EXTERNAL ::  NEIGHB, IO_ASSIGN, IO_CLOSE

C **********************************************************************
C INTEGER IA               : Atom whose neighbours are needed.
C                            A routine initialization must be done by
C                            a first call with IA = 0
C                            If IA0=0, point X0 is used as origin instead
C INTEGER ISEL             : Single-counting switch (0=No, 1=Yes). If ISEL=1,
C                            only neighbours with JA.LE.IA are included in JNA
C INTEGER NNA              : Number of non-zero orbitals at a point in 
C                            real space
C INTEGER JNA(MAXNA)       : Atom index of neighbours. The neighbours
C                            atoms might be in the supercell
C REAL*8  XIJ(3,MAXNA)     : Vectors from point in real space to orbitals
C REAL*8  R2IJ(MAXNA)      : Squared distance to atomic orbitals
C REAL*8  XPO(3)           : Coordinates of the point of the plane respect
C                            we are going to calculate the neighbours orbitals
C INTEGER IZA(NA)          : Atomic number of each atom
C **********************************************************************

      
      ! The first dimension is the number of real numbers per orbital                                     
      ! 1 for real wfs, 2 for complex, and four for the two spinor components                             

      if (non_coll) then
        allocate(wf_single(4,1:no_u))
        allocate(wf(1:no_u,2))
      else
        if (gamma_wfsx) then
           allocate(wf_single(1,1:no_u))
           allocate(wf(1:no_u,1))
        else
           allocate(wf_single(2,1:no_u))
           allocate(wf(1:no_u,1))
        endif
      endif

C Initialize neighbour subroutine --------------------------------------
      IA = 0
      ISEL = 0
      RMAX = RMAXO
      NNA  = MAXNA
      IF (ALLOCATED(JNA)) THEN
        CALL MEMORY('D','I',SIZE(JNA),'stm')
        DEALLOCATE(JNA)
      ENDIF
      IF (ALLOCATED(R2IJ)) THEN
        CALL MEMORY('D','D',SIZE(R2IJ),'stm')
        DEALLOCATE(R2IJ)
      ENDIF
      IF (ALLOCATED(XIJ)) THEN
        CALL MEMORY('D','D',SIZE(XIJ),'stm')
        DEALLOCATE(XIJ)
      ENDIF

      ALLOCATE(JNA(MAXNA))
      CALL MEMORY('A','I',MAXNA,'stm')
      ALLOCATE(R2IJ(MAXNA))
      CALL MEMORY('A','D',MAXNA,'stm')
      ALLOCATE(XIJ(3,MAXNA))
      CALL MEMORY('A','D',3*MAXNA,'stm')

      ALLOCATE(CW(0:NPX-1,0:NPY-1))
      CALL MEMORY('A','Z',NPX*NPY,'stm')
      ALLOCATE(CWE(0:NPX-1,0:NPY-1,0:NPZ-1))
      CALL MEMORY('A','Z',NPX*NPY*NPZ,'stm')
      ALLOCATE(RHO(0:NPX-1,0:NPY-1,0:NPZ-1))
      CALL MEMORY('A','D',NPX*NPY*NPZ,'stm')

      FIRST = .TRUE.
      DO I = 1,3
        XPO(I) = 0.D0
      ENDDO

      CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .             NNA, JNA, XIJ, R2IJ, FIRST )
 
      FIRST = .FALSE.
      RMAX2 =  RMAXO**2

      IF (non_coll)  STOP 'stm: Not ready for NC/SOC yet'

      IF (.not. monoclinic(ucell)) then
        WRITE(6,*) 'error: the code only accepts monoclinic cells'
        WRITE(6,*) '       with Z as the vertical axis'
        STOP
      ENDIF

! Initialize density

      RHO = 0
C Loop over k-points and wavefunctions to include in the STM image

      allocate(wk(nk))
      DO IK  = 1, NK
         read(wf_unit) idummy, k(1:3), wk(ik)
         if (idummy /= ik) stop "ik index mismatch in WFS file"
         read(wf_unit) ispin
         read(wf_unit) number_of_wfns

         WRITE(6,*) 'stm:  Processing kpoint ',IK
         WRITE(6,*) 'stm:  nwf: ', number_of_wfns
         WRITE(6,*) '     --------------------------------'
         DO IWF = 1, number_of_wfns
            read(wf_unit) idummy
            if (idummy /= iwf) then
               ! The file holds a subset of wfs, with the original indexes...
               WRITE(6,*) 'Original wf index: ', idummy
            endif
            read(wf_unit) ener

C Check that we have a bound state (E below vacuum level)

            IF (ENER .LT. EMIN .OR. ENER .GT. EMAX) then
               read(wf_unit)  ! skip wfn info
               CYCLE
            ENDIF

            IF (ENER .GT. V0) THEN
               WRITE(6,*) 'ERROR: ENERGY EIGENVALUE ',IWF,
     .              ' FOR K-POINT ', IK, 'FOR SPIN ',ISPIN
               WRITE(6,*) '       IS ABOVE VACUUM LEVEL'
               STOP
            ENDIF

            WRITE(6,"(a,i5,i2)") 'stm: wf (spin) in window: ',iwf,ispin

            read(wf_unit) (wf_single(:,io), io=1,no_u)
            ! Use a double precision form in what follows
            if ( non_coll) then
               wf(:,1) = cmplx(wf_single(1,:), wf_single(2,:), kind=dp)
               wf(:,2) = cmplx(wf_single(3,:), wf_single(4,:), kind=dp)
            else
               if (gamma_wfsx) then
                  wf(:,1) = cmplx(wf_single(1,:), 0.0_sp, kind=dp)
               else
                  wf(:,1) = cmplx(wf_single(1,:),wf_single(2,:),kind=dp)
               endif
            endif
               
! Loop over all points in real space -----------------------------------

             DO NZ = 1,NPZ

                XPO(3) = ZMIN + (NZ-1)*(ZMAX-ZMIN)/NPZ
                if ( XPO(3) < Zref ) then
                  ! Initialize density to unextrapolated density
          
                   WRITE(6,"(a,f10.4)") 'stm: Using plain LDOS for z =',
     $                                  xpo(3)
                   DO NY = 1,NPY
                      DO NX = 1,NPX

                         XPO(1) = (NX-1)*UCELL(1,1)/NPX +
     $                            (NY-1)*UCELL(1,2)/NPY 
                         XPO(2) = (NX-1)*UCELL(2,1)/NPX +
     $                            (NY-1)*UCELL(2,2)/NPY 

                         call get_cwave(wf(:,1))

                         RHO(NX-1,NY-1,NZ-1)  = RHO (NX-1,NY-1,NZ-1)    
     &                    + DREAL(CWAVE*DCONJG(CWAVE))* ARMUNI * WK(IK)

                      ENDDO  
                   ENDDO

                else

                   ! Extrapolate from reference plane
                   ! Compute value of the wfn at this reference plane
                   WRITE(6,"(a,i4)") 'stm: Extrapolating from nz:', nz

                   DO NY = 1,NPY
                      DO NX = 1,NPX

                         XPO(1) = (NX-1)*UCELL(1,1)/NPX +
     $                            (NY-1)*UCELL(1,2)/NPY 
                         XPO(2) = (NX-1)*UCELL(2,1)/NPX +
     $                            (NY-1)*UCELL(2,2)/NPY 
                         XPO(3) = ZREF

                         call get_cwave(wf(:,1))
                         CW(NX-1,NY-1)  = CWAVE 

                      ENDDO  
                   ENDDO  

                   CALL EXTRAPOLATE(NPX,NPY,NPZ,ZREF,ZMIN,ZMAX,UCELL,V0,
     .                     CW,ENER,K,CWE)
                   ! Be careful not to overwrite the z<zref parts...
                   RHO(:,:,NZ-1:) = RHO(:,:,NZ-1:) +
     $                         DREAL(CWE(:,:,NZ-1:)*
     $                               DCONJG(CWE(:,:,NZ-1:)))
     $                         * WK(IK) * ARMUNI

                   ! And we are done with the z planes
                   EXIT  ! loop over NZ

                endif    ! z below or above Zref

             ENDDO       ! NZ
          ! Note that values for different spins are now accumulated
       ENDDO     ! wfn number
      ENDDO      ! k-point

      ! This should not be necessary if a proper BZ-sampled set of wfs is used
      total_weight = sum(wk(1:nk))
      rho = rho / total_weight
      ! Normalize if not spin-polarized
      if ((nspin_blocks == 1) .and. (.not. non_coll))  then
         rho = 2.0_dp * rho
      endif
         

      call io_assign(unitre1)
      SNAME = FDF_STRING('SystemLabel','siesta')
      stm_label = FDF_STRING('stm-label','')
      if (stm_label == '') then
         FNAME = trim(SNAME) // '.STM.cube'
      else
         FNAME = trim(SNAME) // '.' // trim(stm_label) // '.STM.cube'
      endif
      
        WRITE(6,*)
        WRITE(6,*) 'stm: writing cube format file ',FNAME
        WRITE(6,*)
        WRITE(6,*) '     ',NSCX,' x ',NSCY,' cells in cube plot'

C Calculate number of atoms in unit cell
      VC = VOLCEL(CELL)
      VU = VOLCEL(UCELL)
      NAU = NA / IDNINT(VC/VU)

      open(unitre1,file=FNAME,form='formatted',status='unknown')
      WRITE(UNITRE1,*) 'STM'
      WRITE(UNITRE1,*) 'STM'
      !! VOXEL sizes are wrong (should divide by NPX, etc)
      WRITE(UNITRE1,'(i5,4f12.6)') NAU*NSCX*NSCY, 0.0, 0.0, ZMIN
      WRITE(UNITRE1,'(i5,4f12.6)') NPX*NSCX,(UCELL(1,J)/(NPX-1),J=1,3)
      WRITE(UNITRE1,'(i5,4f12.6)') NPY*NSCY,(UCELL(2,J)/(NPY-1),J=1,3)
      WRITE(UNITRE1,'(i5,4f12.6)') NPZ,0.0,0.0,((ZMAX-ZMIN)/(NPZ-1))


      DO NSX = 1,NSCX
        DO NSY = 1,NSCY
           DO IA = 1, NAU
              WRITE(UNITRE1,'(i5,4f12.6)') ATOMIC_NUMBER(ISA(IA)),0.0, 
     .     (XA(IX,IA)+(NSX-1)*UCELL(IX,1)+(NSY-1)*UCELL(IX,2) ,IX=1,3)
           ENDDO
        ENDDO
      ENDDO

      DO NSX = 1,NSCX
         DO NX=0,NPX-1
            DO NSY = 1,NSCY
               DO NY=0,NPY-1
                  WRITE(UNITRE1,'(6e13.5)')
     .                 (RHO(NX,NY,NZ),NZ=0,NPZ-1)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      call io_close(unitre1)

! Write charge density in Siesta format

      call io_assign(unitre1)
      SNAME = FDF_STRING('SystemLabel','siesta')
      stm_label = FDF_STRING('stm-label','')
      if (stm_label == '') then
         FNAME = trim(SNAME) // '.STM.siesta'
      else
         FNAME = trim(SNAME) // '.' // trim(stm_label) // '.STM.siesta'
      endif

      WRITE(6,*)
      WRITE(6,*) 'stm: writing SIESTA format file ', FNAME
      WRITE(6,*)
      open(unitre1,file=FNAME,form='unformatted',
     .         status='unknown')
C write the correct range of the z axis: temporarily override ucell
      USAVE = UCELL(3,3)
      UCELL(3,3) = abs(ZMAX-ZMIN)
      WRITE(unitre1) UCELL
C restore ucell
      UCELL(3,3)=USAVE
      WRITE(unitre1) NPX, NPY, NPZ, 1

      DO IZ=0,NPZ-1
        DO IY=0,NPY-1
          WRITE(unitre1) (REAL(RHO(IX,IY,IZ)),IX=0,NPX-1)
        ENDDO
      ENDDO

      call io_close(unitre1)



C CLOSE ALLOCATABLE ARRAYS


      DEALLOCATE(RHO)
      CALL MEMORY('D','D',NPX*NPY*NPZ,'stm')
      DEALLOCATE(CWE)
      CALL MEMORY('D','Z',NPX*NPY*NPZ,'stm')
      DEALLOCATE(CW)
      CALL MEMORY('D','Z',NPX*NPY,'stm')

      CONTAINS

      subroutine get_cwave(psi)
      complex(dp), intent(in) :: psi(:)
      
      ! Inherits all data by host association

! The periodic part of the Bloch functions is defined by                                                 
! \begin{equation}                                                                                       
!   u_{n \vec{k}} (\vec{r}) =                                                                            
!   \sum_{\vec{R} \mu} c_{n \mu}(\vec{k})                                                                
!        e^{i \vec{k} \cdot ( \vec{r}_{\mu} + \vec{R} - \vec{r} )}                                       
!        \phi_{\mu} (\vec{r} - \vec{r}_{\mu} - \vec{R} ) ,                                               
!\end{equation}                                                                                          
!                                                                                                        
!\noindent where $\phi_{\mu} (\vec{r} - \vec{r}_{\mu} - \vec{R} )$                                       
! is an atomic orbital of the basis set centered on atom $\mu$ in                                        
! the unit cell $\vec{R}$, and $c_{n \mu}(\vec{k})$ are the coefficients                                 
! of the wave function

      CWAVE   = (0.0D0, 0.0D0)

!     CWAVE is meant to be the periodic part of the wavefunction,
!     for both the computation of the charge density directly (exp(ikr) phase
!     is irrelevant) and for propagation of the wave function (?)
      
      ! First step (see below)
      ! Phase to cancel the phase of the wave function: -i.k.r
      
      PMIKR = -(K(1)*XPO(1) + K(2)*XPO(2) + K(3)*XPO(3))
      SIMIKR=DSIN(PMIKR)
      COMIKR=DCOS(PMIKR)
      EXMIKR=DCMPLX(COMIKR,SIMIKR)

C Localize non-zero orbitals at each point in real space ---------------
     
      IA   = 0
      ISEL = 0
      NNA  = MAXNA
      ! Get neighbors of point xpo
      CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .                   NNA, JNA, XIJ, R2IJ, FIRST )

C     Loop over Non-zero orbitals ------------------------------------------
      ! NOTE: If the z-extent of the box is large enough, we might be getting
      ! contributions from orbitals in the next periodic image of the slab.
      ! We should get rid of them
      DO  IAT1 = 1, NNA
         IF( R2IJ(IAT1) .GT. RMAX2 ) CYCLE

         IAVEC1   = JNA(IAT1)
         IS1      = ISA(IAVEC1)
         XVEC1(:) = -XIJ(:,IAT1) ! position of XPO with respect to
                                       ! the atom

         !  XPO + XIJ(IAT1) is just the absolute position of atom IAT1
         !  We could cancel the phase above and keep only k*xij

         PHASE = K(1)*(XPO(1)+XIJ(1,IAT1))+
     .        K(2)*(XPO(2)+XIJ(2,IAT1))+
     .        K(3)*(XPO(3)+XIJ(3,IAT1))

         SI=DSIN(PHASE)
         CO=DCOS(PHASE)
         EXPPHI=DCMPLX(CO,SI)

         DO IO = LASTO(IAVEC1-1) + 1, LASTO(IAVEC1)
            IPHI1 = IPHORB(IO)
            IUO   = INDXUO(IO)
            CALL PHIATM( IS1, IPHI1, XVEC1, PHIMU, GRPHIMU )

            CWAVE = CWAVE + PHIMU * psi(iuo) * EXPPHI * EXMIKR

         ENDDO
      ENDDO

      end subroutine get_cwave
      
      function monoclinic(cell)
      real(dp), intent(in) :: cell(3,3)
      logical monoclinic

      real(dp), parameter :: tol = 1.0e-8_dp

      ! This is too naive. It should be checking that a_3 is orthogonal
      ! to both a_1 and a_2
      
      monoclinic =  (abs(CELL(3,1)) < tol
     $         .and. abs(CELL(3,2)) < tol
     $         .and. abs(CELL(1,3)) < tol
     $         .and. abs(CELL(2,3)) < tol )

      print *, "monoclinic: ", monoclinic
      print *, cell

         end function monoclinic

      END
