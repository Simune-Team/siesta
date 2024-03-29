! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      SUBROUTINE STM( NA, NO, NO_U, MAXNA, nspin,nspin_blocks,non_coll,
     .                ISA, IPHORB, INDXUO, LASTO, XA, CELL, UCELL,
     .                wf_unit, NK, gamma_wfsx,
     .                ZREF, ZMIN, ZMAX, NPX, NPY, NPZ, 
     .                V0, EMAX, EMIN,
     .                ARMUNI, IUNITCD, RMAXO )

C **********************************************************************
C Simulate STM images in the Tersoff-Hamann approximation, by
C extrapolating the wavefunctions into vacuum
C
C Coded by P. Ordejon and N. Lorente,  November 2004
C
C     Modified by N. Lorente, August 2005
      ! Restructured by A. Garcia, March 2019
C **********************************************************************

      use precision, only: dp, sp
      USE ATMFUNCS
      USE FDF
      USE CHEMICAL

      IMPLICIT NONE

      real(dp), parameter :: Ang    = 1.0_dp / 0.529177_dp
      
      INTEGER, INTENT(IN) ::
     .  NA, NO, NO_U, NPX, NPY, NPZ, IUNITCD,
     .  nspin, nspin_blocks, MAXNA, NK,
     .  ISA(NA), IPHORB(NO), INDXUO(NO), LASTO(0:NA)
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
     .  IS1, IPHI1, NX, NY, NZ, IWF, IK, ISPIN, grid_u, str_u,
     .  IX, IY, IZ, NSX, NSY, NAU, iv, is

      REAL(DP)
     .  DOT, RMAX, XPO(3), RMAX2, XVEC1(3),
     .  PHIMU, GRPHIMU(3),
     .  PHASE, SI, CO, ENER, PMIKR, SIMIKR, COMIKR, USAVE, VC, VU

      real(dp) :: total_weight, k(3)
      
      REAL(DP), ALLOCATABLE :: RHO(:,:,:,:)
      REAL(SP), ALLOCATABLE :: wf_single(:,:)
      COMPLEX(DP), ALLOCATABLE :: wf(:,:)
      REAL(DP), ALLOCATABLE :: wk(:)

      COMPLEX(DP)  EXPPHI, EXMIKR, d11, d12, d21, d22

      ! The last dimension of these is the number of spinor components
      ! 1 for collinear, and 2 for NC/SOC
      COMPLEX(DP), ALLOCATABLE :: CW(:,:,:), CWE(:,:,:,:), CWAVE(:)
 
      LOGICAL FIRST
      integer :: idummy, number_of_wfns, spinor_comps
      integer :: nspin_rho  ! Number of components needed for rho array ("spin%Grid")

      CHARACTER ::  SNAME*40, FNAME*256, stm_label*60

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

      
      ! The first dimension of wf_single is the number of real numbers per orbital
      ! to be read from the WFSX file:
      ! 1 for real wfs, 2 for complex, and four for the two spinor components
      ! wf is a complex array which holds either a wfn or a two-component spinor.

      if (non_coll) then
        allocate(wf_single(4,1:no_u))
        allocate(wf(1:no_u,2))
        spinor_comps = 2
      else
        spinor_comps = 1
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
        DEALLOCATE(JNA)
      ENDIF
      IF (ALLOCATED(R2IJ)) THEN
        DEALLOCATE(R2IJ)
      ENDIF
      IF (ALLOCATED(XIJ)) THEN
        DEALLOCATE(XIJ)
      ENDIF

      ALLOCATE(JNA(MAXNA))
      ALLOCATE(R2IJ(MAXNA))
      ALLOCATE(XIJ(3,MAXNA))

      nspin_rho = min(4,nspin)

      allocate(CWAVE(spinor_comps))
      ALLOCATE(CW(0:NPX-1,0:NPY-1,spinor_comps))
      ALLOCATE(CWE(0:NPX-1,0:NPY-1,0:NPZ-1,spinor_comps))
      ALLOCATE(RHO(0:NPX-1,0:NPY-1,0:NPZ-1,nspin_rho))

      FIRST = .TRUE.
      DO I = 1,3
        XPO(I) = 0.D0
      ENDDO
      CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .             NNA, JNA, XIJ, R2IJ, FIRST )
      FIRST = .FALSE.
      RMAX2 =  RMAXO**2

      IF (.not. monoclinic(ucell)) then
        WRITE(6,*) 'error: the code only accepts monoclinic cells'
        WRITE(6,*) '       with Z as the vertical axis'
        STOP
      ENDIF

! Initialize density

      RHO = 0

!     Loop over k-points and wavefunctions to include in the STM image

      allocate(wk(nk))
      DO IK  = 1, NK
        do ispin = 1, nspin_blocks
         read(wf_unit) idummy, k(1:3), wk(ik)
            if (idummy /= ik) then
               write(6,*) "ik index mismatch in WFS file"
               WRITE(6,*) "ik in file, ik: ", idummy, ik
            endif
         read(wf_unit) idummy
            if (idummy /= ispin) then
               write(6,*) "ispin index mismatch in WFS file"
               WRITE(6,*) "ispin in file, ispin: ", idummy, ispin
            endif
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

            ! Check that we have a bound state (E below vacuum level),
            ! in the chosen window

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
            ! Use a double precision complex form in what follows
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
               
             ! Loop over all points in real space
             ! The last point (zmax) is now included by making stepz=(Zmax-Zmin)/(NPZ-1)
             ! This forces the definition of a slightly larger c vector below.
             ! In this way, the last plane recorded in the file will correspond to Z=Zmax
             DO NZ = 0, NPZ-1

                if (npz == 1) then
                   XPO(3) = ZMIN
                else
                   XPO(3) = ZMIN + NZ*(ZMAX-ZMIN)/(NPZ-1)
                endif

                if ( XPO(3) < Zref ) then
                  ! Initialize density to unextrapolated density
          
                   WRITE(6,"(a,f10.4)") 'stm: Using plain LDOS for z =',
     $                                  xpo(3)
                   DO NY = 0,NPY-1
                      DO NX = 0,NPX-1

                         ! Note that the (periodic) X and Y directions
                         ! are treated as usual, with smaller step
                         XPO(1) = NX*UCELL(1,1)/NPX +
     $                            NY*UCELL(1,2)/NPY 
                         XPO(2) = NX*UCELL(2,1)/NPX +
     $                            NY*UCELL(2,2)/NPY 

                         call get_cwave(wf(:,1:spinor_comps))

                         ! Now for the various cases
                         if (nspin <= 2) then
                            RHO(NX,NY,NZ,ispin)  = RHO (NX,NY,NZ,ispin)    
     &                           + REAL(CWAVE(1)*CONJG(CWAVE(1)), dp)
     $                             * ARMUNI * WK(IK)
                         else   ! non-collinear
                            ! CHECK THIS
                            d11 = cwave(1) * conjg(cwave(1))
                            d12 = cwave(1) * conjg(cwave(2))
                            d21 = cwave(2) * conjg(cwave(1))
                            d22 = cwave(2) * conjg(cwave(2))

                            ! Hermitify?
                            D12 = 0.5_dp * (D12 + conjg(D21))
                            
                            ! Recall: dm(:,3) = real(d12);  dm(:,4) = -aimag(d12)
                            rho(nx,ny,nz,1) = rho(nx,ny,nz,1) 
     $                                 + real(d11,dp) * armuni * wk(ik)
                            rho(nx,ny,nz,2) = rho(nx,ny,nz,2) 
     $                                 + real(d22,dp) * armuni * wk(ik)
                            rho(nx,ny,nz,3) = rho(nx,ny,nz,3)
     $                                 + real(d12,dp) * armuni * wk(ik)
                            rho(nx,ny,nz,4) = rho(nx,ny,nz,4)
     $                                 - aimag(d12) * armuni * wk(ik)
                         endif
                      ENDDO  
                   ENDDO

                else

                   ! Extrapolate from reference plane
                   ! Compute value of the wfn at this reference plane
                   WRITE(6,"(a,i4)") 'stm: Extrapolating from nz:', nz

                   DO NY = 0,NPY-1
                      DO NX = 0,NPX-1

                         XPO(1) = NX*UCELL(1,1)/NPX +
     $                            NY*UCELL(1,2)/NPY 
                         XPO(2) = NX*UCELL(2,1)/NPX +
     $                            NY*UCELL(2,2)/NPY 
                         XPO(3) = ZREF

                         call get_cwave(wf(:,1:spinor_comps))
                         CW(NX,NY,1:spinor_comps) =
     $                                      CWAVE(1:spinor_comps) 

                      ENDDO  
                   ENDDO  

                   ! This is mildly wasteful in terms of initialization
                   ! of FFTW, but it will do for now
                   ! We assume that both spinor components are propagated
                   ! in the same way
                   do is = 1, spinor_comps
                      CALL EXTRAPOLATE(NPX,NPY,NPZ,ZREF,ZMIN,ZMAX,
     $                         UCELL,V0,
     .                         CW(0,0,is),ENER,K,CWE(0,0,0,is))
                   enddo
                   
                   ! Now for the various cases
                   ! Be careful not to overwrite the z<zref parts...
                   if (nspin <= 2) then
                      RHO(:,:,NZ:,ispin)  = RHO (:,:,NZ:,ispin)    
     &                     + REAL(CWE(:,:,NZ:,1)*
     $                     CONJG(CWE(:,:,NZ:,1)), dp)
     $                     * ARMUNI * WK(IK)
                   else         ! non-collinear

                      do iz = nz, npz-1
                         do ny=0,npy-1
                            do nx=0,npx-1
                         d11 = cwe(nx,ny,iz,1) * conjg(cwe(nx,ny,iz,1))
                         d12 = cwe(nx,ny,iz,1) * conjg(cwe(nx,ny,iz,2))
                         d21 = cwe(nx,ny,iz,2) * conjg(cwe(nx,ny,iz,1))
                         d22 = cwe(nx,ny,iz,2) * conjg(cwe(nx,ny,iz,2))

                         !     Hermitify?
                         D12 = 0.5_dp * (D12 + conjg(D21))
                            
                      ! Recall: dm(:,3) = real(d12);  dm(:,4) = -aimag(d12)
                         rho(nx,ny,iz,1) = rho(nx,ny,iz,1) 
     $                                 + real(d11,dp) * armuni * wk(ik)
                         rho(nx,ny,iz,2) = rho(nx,ny,iz,2) 
     $                                 + real(d22,dp) * armuni * wk(ik)
                         rho(nx,ny,iz,3) = rho(nx,ny,iz,3)
     $                                 + real(d12,dp) * armuni * wk(ik)
                         rho(nx,ny,iz,4) = rho(nx,ny,iz,4)
     $                                - aimag(d12) * armuni * wk(ik)
                        enddo
                      enddo
                    enddo
                   endif  ! collinear or not

                   ! And we are done with the z planes
                   EXIT  ! loop over NZ

                endif    ! z below or above Zref

             ENDDO       ! NZ
          ENDDO     ! ispin = 1, nspin_blocks
                    
       ENDDO     ! wfn number
      ENDDO      ! k-point

      ! This should not be necessary if a proper BZ-sampled set of wfs is used
      total_weight = sum(wk(1:nk))
      rho = rho / total_weight
      ! Normalize if not spin-polarized
      if ((nspin_blocks == 1) .and. (.not. non_coll))  then
         rho = 2.0_dp * rho
      endif

! Write charge density in Siesta format

      call io_assign(grid_u)
      SNAME = FDF_STRING('SystemLabel','siesta')
      stm_label = FDF_STRING('stm-label','')
      if (stm_label == '') then
         FNAME = trim(SNAME) // '.STM.LDOS'
      else
         FNAME = trim(SNAME) // '.' // trim(stm_label) // '.STM.LDOS'
      endif

      WRITE(6,*)
      WRITE(6,*) 'stm: writing SIESTA format file ', FNAME
      WRITE(6,*)

      open(grid_u,file=FNAME,form='unformatted',
     .         status='unknown')
      ! write the correct range of the z axis: temporarily override ucell
      USAVE = UCELL(3,3)
      ! Make the cell slightly taller, so that the last (NPZ-1) plane corresponds
      ! to Z=Zmax
      if (npz == 1 ) then
         UCELL(3,3) = 1.0_dp   ! We have a single plane. Use a 1.0 bohr-thick height
      else
         UCELL(3,3) = NPZ * abs(ZMAX-ZMIN) / (NPZ-1)
      endif
      
      WRITE(grid_u) UCELL,
     $              [0.0_dp, 0.0_dp, ZMIN], ! Extra info for origin
     $              [.true.,.true.,.false.] ! Periodic ?

      WRITE(grid_u) NPX, NPY, NPZ, nspin_rho 

      do ispin = 1, nspin_rho
         DO IZ=0,NPZ-1
            DO IY=0,NPY-1
               WRITE(grid_u) (REAL(RHO(IX,IY,IZ,ispin),sp),IX=0,NPX-1)
            ENDDO
         ENDDO
      enddo
      call io_close(grid_u)

      ! Write a dummy STRUCT file for use with the 'grid to cube' converter
      call io_assign(str_u)
      open(str_u,file=trim(SNAME)//'.CELL_STRUCT',
     $     form='formatted', status='unknown')
      write(str_u,'(3x,3f18.9)') ((UCELL(ix,iv)/Ang,ix=1,3),iv=1,3)
      write(str_u,*) 0  ! number of atoms
      close(str_u)
      
      ! restore ucell
      UCELL(3,3)=USAVE

      DEALLOCATE(RHO)
      DEALLOCATE(CWE)
      DEALLOCATE(CW)
      deallocate(cwave)

      CONTAINS

      subroutine get_cwave(psi)
      complex(dp), intent(in) :: psi(:,:)  ! Can deal with spinors
      
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
      SIMIKR=SIN(PMIKR)
      COMIKR=COS(PMIKR)
      EXMIKR=CMPLX(COMIKR,SIMIKR,kind=dp)

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

         SI=SIN(PHASE)
         CO=COS(PHASE)
         EXPPHI=CMPLX(CO,SI,kind=dp)

         DO IO = LASTO(IAVEC1-1) + 1, LASTO(IAVEC1)
            IPHI1 = IPHORB(IO)
            IUO   = INDXUO(IO)
            CALL PHIATM( IS1, IPHI1, XVEC1, PHIMU, GRPHIMU )

            ! Note implicit loop over spinor components
            CWAVE(:) = CWAVE(:) + PHIMU * psi(iuo,:) * EXPPHI * EXMIKR

         ENDDO
      ENDDO

      end subroutine get_cwave
      
      function monoclinic(cell)
      real(dp), intent(in) :: cell(3,3)
      logical monoclinic

      real(dp), parameter :: tol = 1.0e-8_dp

      ! This is too naive. It should be checking that a_3 is orthogonal
      ! to both a_1 and a_2, but it is fine for this use case.
      
      monoclinic =  (abs(CELL(3,1)) < tol
     $         .and. abs(CELL(3,2)) < tol
     $         .and. abs(CELL(1,3)) < tol
     $         .and. abs(CELL(2,3)) < tol )

      end function monoclinic

      END
