      module periodic_table

      use precision
      use sys

      implicit none

      CONTAINS


      SUBROUTINE QVLOFZ( Z, QVAL )  
C Returns the atomic LSD ground state configurations, as an array with 
C valence population for each angular momentum L. The valence orbitals 
C for each L are those returned by routine LMXOFZ
C Originally written by A.R.Williams. Modified by J.M.Soler

      integer,  intent(in)  :: Z        ! Atomic number
      real(dp), intent(out) :: QVAL(0:) ! Valence charge for each L

      integer, parameter :: LMAX=3, ZMAX=98, NCHNG=15
      integer :: I, ICHNG, L, LMXATM, LMXCHM, LCHNG(NCHNG),
     .           N(0:LMAX), NVAL(0:5), Q(0:LMAX,0:ZMAX), ZCHNG(NCHNG)

      ! Notice: s valence orbital switched for p occupation = 4
      !           Li, F,Na,Cl, K,Ga,Br,Rb,In, I,Cs,Hf,Tl,At,Fr
      DATA ZCHNG / 3, 9,11,17,19,31,35,37,49,53,55,72,81,85,87/
      DATA LCHNG / 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1/

      ! Notice that s valence charge must be consitent with ZCHNG
      DATA (Q(0,I),I= 0, 2) /0,1,2/
      DATA (Q(0,I),I= 3,10) /1,2,                      2,2,2,2,0,0/
      DATA (Q(0,I),I=11,18) /1,2,                      2,2,2,2,0,0/
      DATA (Q(0,I),I=19,36) /1,2, 2,2,2,1,2,2,2,2,1,2, 2,2,2,2,0,0/
      DATA (Q(0,I),I=37,54) /1,2, 2,2,1,1,2,1,1,0,1,2, 2,2,2,2,0,0/
      DATA (Q(0,I),I=55,71) /1,2, 2,                               14*2/
      DATA (Q(0,I),I=72,86) /       2,2,2,2,2,2,0,1,2, 2,2,2,2,0,0/
      DATA (Q(0,I),I=87,98) /1,2, 10*2/

      DATA (Q(1,I),I= 0, 2) / 3*0/                ! H-He
      DATA (Q(1,I),I= 3,10) / 2*0, 1,2,3,4,5,6/   ! Li-Ne
      DATA (Q(1,I),I=11,18) / 2*0, 1,2,3,4,5,6/   ! Na-Ar
      DATA (Q(1,I),I=19,36) /12*0, 1,2,3,4,5,6/   ! K-Kr
      DATA (Q(1,I),I=37,54) /12*0, 1,2,3,4,5,6/   ! Rb-Xe
      DATA (Q(1,I),I=55,86) /26*0, 1,2,3,4,5,6/   ! Cs-Rn
      DATA (Q(1,I),I=87,98) /12*0/                ! Fr-Cf

      DATA (Q(2,I),I= 0,36) /21*0, 1,2,3,5,5,6,7, 8,10,10, 6*0/
      DATA (Q(2,I),I=37,54) / 2*0, 1,2,4,5,5,7,8,10,10,10, 6*0/
      DATA (Q(2,I),I=55,71) / 2*0, 1,                      6*0,1,6*0,1/
      DATA (Q(2,I),I=72,86) /        2,3,4,5,6,7,10,10,10, 6*0/
      DATA (Q(2,I),I=87,98) / 2*0, 1,2,1,1,3*0,1,2,1/

      DATA (Q(3,I),I= 0,71) /58*0, 2,3,4,5,6,7,7,9,10,11,12,13,14,14/
      DATA (Q(3,I),I=72,98) /18*0, 0,2,3,5,6,7,7,7,9/
 
      IF (Z.GT.ZMAX) call die('QVLOFZ: ERROR: Z out of range')

      ! Find the principal quantum numbers assigned by my own data
      DO L=0,LMAX
         N(L)=L+1
      END DO
      DO ICHNG=1,NCHNG
         IF (ZCHNG(ICHNG).GT.Z) EXIT
         L=LCHNG(ICHNG)
         N(L)=N(L)+1
      END DO

      ! Find the valence principal quantum numbers assigned by GNFIG
      CALL LMXOFZ (Z,LMXCHM,LMXATM)
      CALL CNFIG (Z,NVAL)

      ! Check size of QVAL and initialize it (for L>LMXCHM)
      IF (UBOUND(QVAL,1).LT.LMXCHM) 
     .  call die('QVLOFZ: ERROR: Size of QVAL too small')
      QVAL(:) = 0

      ! Make valence occupation consistent with CNFIG valence assignment
      DO L=0,LMXCHM
        IF (NVAL(L).GT.N(L)) THEN
          QVAL(L)=0
        ELSE
          QVAL(L) = Q(L,Z) + 2*(2*L+1)*(N(L)-NVAL(L))
        ENDIF
      END DO

      END subroutine qvlofz


      SUBROUTINE LMXOFZ( Z, LMXVAL, LMXATM ) 
C Given the atomic number Z, returns the maximum angular mometum L
C which is populated in the atomic ground state configuration:
C For the valence orbitals => LMXVAL. For any orbitals => LMXATM
C Originally written by A.R.Williams. Modified by J.M.Soler

      integer,intent(in) :: Z      ! Atomic number
      integer,intent(out):: LMXVAL ! Max. L for valence states
      integer,intent(out):: LMXATM ! Max. L for valence and core states

      integer, parameter :: NCHNG=17
      integer :: ICHNG, LCHNG(NCHNG), ZCHNG(NCHNG)

      !            B,Na,Al, K,Sc,Ga,Rb, Y,In,Cs,La,Ce,Hf,Tl,Fr,Ac,Th
      DATA ZCHNG / 5,11,13,19,21,31,37,39,49,55,57,58,72,81,87,89,90/
      DATA LCHNG / 1, 0, 1, 0, 2, 1, 0, 2, 1, 0, 2, 3, 2, 1, 0, 2, 3/

      LMXVAL=0
      LMXATM=0
      DO ICHNG=1,NCHNG
         IF (ZCHNG(ICHNG).GT.Z) EXIT
         LMXVAL=LCHNG(ICHNG)
         LMXATM=MAX(LMXATM,LMXVAL)
      ENDDO

      END subroutine lmxofz


      SUBROUTINE CNFIG( Z, NVAL ) 
C Returns the valence configuration for atomic ground state, i.e.
C the principal quantum number NVAL of the valence orbilas for each L
C Originally written by A.R.Williams. Modified by J.M.Soler

      integer,intent(in) :: Z        ! Atomic number
      integer,intent(out):: NVAL(0:) ! Valence electrons for each L

      integer, parameter :: LMAX=3, NCHNG=15
      integer :: ICHNG, L, LCHNG(NCHNG), ZCHNG(NCHNG)

      ! Originally: s valence orbital switched for p occupation = 4
      !           Li, F,Na,Cl, K,Ga,Br,Rb,In, I,Cs,Hf,Tl,At,Fr
*     DATA ZCHNG / 3, 9,11,17,19,31,35,37,49,53,55,72,81,85,87/

      ! Changed to: s valence orbital switched for full p occupation
      !           Li,Na,Na, K, K,Ga,Rb,Rb,In,Cs,Cs,Hf,Tl,Fr,Fr
      DATA ZCHNG / 3,11,11,19,19,31,37,37,49,55,55,72,81,87,87/
      DATA LCHNG / 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1/
      DO L=0,LMAX
         NVAL(L)=L+1
      END DO
      DO ICHNG=1,NCHNG
         IF (ZCHNG(ICHNG).GT.Z) EXIT
         L=LCHNG(ICHNG)
         NVAL(L)=NVAL(L)+1
      END DO

      END subroutine cnfig


      FUNCTION SYMBOL( Z )
C Given the atomic number, returns the atomic symbol (e.g. 'Na')
C Written by J. Soler

      character(len=2)    :: SYMBOL  ! Atomic symbol
      integer, intent(in) :: Z       ! Atomic number

      integer, parameter  :: NZ=103
      character(len=2) :: NAME(NZ)
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

      IF (Z.EQ.0 .OR. Z.EQ.-100) THEN
         SYMBOL = 'BS'
      ELSE IF (ABS(Z).LE.NZ) THEN
         SYMBOL = NAME(ABS(Z))
      ELSE
         WRITE(6,*) 'SYMBOL: ERROR: No data for Z =', Z
         SYMBOL = ' '
      ENDIF

      END function symbol


      FUNCTION ATMASS( Z )
C Returns the average atomic mass from the atomic number Z.
C Dta taken from VCH periodic table.
C Written by J.M.Soler. April'97.

      real(dp)             :: ATMASS ! Average atomic mass, in amu
      integer, intent(in)  :: Z      ! Atomic number

      integer, PARAMETER  :: NZ=94
      character(len=50) message

      DOUBLE PRECISION AMASS(0:NZ)
      DATA AMASS / 0.00,
     .     1.01,  4.00,  6.94,  9.01, 10.81, 12.01, 14.01, 16.00,
     .    19.00, 20.18, 22.99, 24.31, 26.98, 28.09, 30.97, 32.07,
     .    35.45, 39.95, 39.10, 40.08, 44.96, 47.88, 50.94, 52.00,
     .    54.94, 55.85, 58.93, 58.69, 63.55, 65.39, 69.72, 72.61,
     .    74.92, 78.96, 79.90, 83.80, 85.47, 87.62, 88.91, 91.22,
     .    92.91, 95.94, 98.91,101.07,102.91,106.42,107.87,112.41,
     .   114.82,118.71,121.75,127.60,126.90,131.29,132.91,137.33,
     .   138.91,140.12,140.91,144.24,146.92,150.36,151.97,157.25,
     .   158.93,162.50,164.93,167.26,168.93,173.04,174.97,178.49,
     .   180.95,183.85,186.21,190.2 ,192.22,195.08,196.97,200.59,
     .   204.38,207.2 ,208.98,208.98,209.99,222.02,223.02,226.03,
     .   227.03,232.04,231.04,238.03,237.05,244.06/

      IF (Z.LT.0 .OR. Z.GT.NZ) THEN
         write(message,'(a,i4)') 'ATMASS: ERROR: No data for Z =',Z
         call die(message)
      ELSE
         ATMASS=AMASS(Z)
      ENDIF

      END function atmass

      end module periodic_table







