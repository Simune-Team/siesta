      module periodic_table

      use precision
      use sys

      implicit none

      CONTAINS

      SUBROUTINE QVLOFZ(IZ,QVAL)  

C Returns an array with the atomic ground states 
C populations.
C Input: IZ atomic number
C Written by A.R.Williams 

      integer, intent(in)  :: iz
      real(dp), intent(out) :: qval(0:)

      integer lmax, izmax, nchng
      PARAMETER (LMAX=3,IZMAX=98,NCHNG=15)
      real(dp)  :: Q(0:LMAX,0:IZMAX)
      integer   :: NVAL(0:5),N(0:LMAX),IZCHNG(NCHNG),LCHNG(NCHNG)

      integer i, l, ichng, lmxchm, lmxatm

      DATA IZCHNG /3,9,11,17,19,31,35,37,49,53,55,72,81,85,87/
      DATA LCHNG  /0,0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1/

      DATA (Q(0,I),I= 0, 2) /0.,1.,2./
      DATA (Q(0,I),I= 3,10) /1.,2.,  2.,2.,2.,2.,0.,0./
      DATA (Q(0,I),I=11,18) /1.,2.,  2.,2.,2.,2.,0.,0./
      DATA (Q(0,I),I=19,30) /1.,2.,  2.,2.,2.,1.,2.,2.,2.,2.,1.,2./
      DATA (Q(0,I),I=31,36) /        2.,2.,2.,2.,0.,0./
      DATA (Q(0,I),I=37,48) /1.,2.,  2.,2.,1.,1.,2.,1.,1.,0.,1.,2./
      DATA (Q(0,I),I=49,54) /        2.,2.,2.,2.,0.,0./
      DATA (Q(0,I),I=55,71) /1.,2.,  2.,  14*2./
      DATA (Q(0,I),I=72,80) /           2.,2.,2.,2.,2.,2.,0.,1.,2./
      DATA (Q(0,I),I=81,86) /        2.,2.,2.,2.,0.,0./
      DATA (Q(0,I),I=87,98) /1.,2.,  10*2./

      DATA (Q(1,I),I= 0,10) / 5*0.,1.,2.,3.,4.,5.,6./
      DATA (Q(1,I),I=11,18) / 2*0.,1.,2.,3.,4.,5.,6./
      DATA (Q(1,I),I=19,36) /12*0.,1.,2.,3.,4.,5.,6./
      DATA (Q(1,I),I=37,54) /12*0.,1.,2.,3.,4.,5.,6./
      DATA (Q(1,I),I=55,86) /26*0.,1.,2.,3.,4.,5.,6./
      DATA (Q(1,I),I=87,98) /12*0./

      DATA (Q(2,I),I= 0,30) /21*0.,1.,2.,3.,5.,5.,6.,7., 8.,10.,10./
      DATA (Q(2,I),I=31,48) / 8*0.,1.,2.,4.,5.,5.,7.,8.,10.,10.,10./
      DATA (Q(2,I),I=49,71) / 8*0.,1.,  6*0.,1.,6*0.,1./
      DATA (Q(2,I),I=72,80) /         2.,3.,4.,5.,6.,7.,10.,10.,10./
      DATA (Q(2,I),I=81,98) / 8*0.,1.,2.,1.,1.,3*0.,1.,2.,1./

      DATA (Q(3,I),I= 0,64) /58*0.,2.,3.,4.,5.,6.,7.,7./
      DATA (Q(3,I),I=65,71) /      9.,10.,11.,12.,13.,14.,14./
      DATA (Q(3,I),I=72,98) /18*0.,0.,2.,3.,5.,6.,7.,7.,7.,9./
 
      qval(:) = 0.0_dp

      IF (IZ.GT.IZMAX) call die('QVLOFZ: IZ GREATER THAN IZMAX')

      DO 10 L=0,LMAX
         N(L)=L+1
  10  CONTINUE
      DO 20 ICHNG=1,NCHNG
         IF (IZ.LT.IZCHNG(ICHNG)) GOTO 20
         L=LCHNG(ICHNG)
         N(L)=N(L)+1
  20  CONTINUE
      CALL LMXOFZ (IZ,LMXCHM,LMXATM)
      CALL CNFIG (IZ,NVAL)
      DO 30 L=0,LMXCHM
         IF (NVAL(L).LT.N(L)) QVAL(L)=2*(2*L+1)+Q(L,IZ)
         IF (NVAL(L).EQ.N(L)) QVAL(L)=Q(L,IZ)
         IF (NVAL(L).GT.N(L)) QVAL(L)=0
  30  CONTINUE
      END subroutine qvlofz


      SUBROUTINE LMXOFZ (Z,LMXCHM,LMXATM) 

C Returns the maximum angular mometum populated in the 
C atomic ground state configuration LMXATM, and for the 
C valence shell LMXCHM
C Input: Z atomic number
C Written by A.R.Williams 

      integer, intent(in)  ::  z
      integer, intent(out) ::  lmxchm, lmxatm

      integer, PARAMETER   :: NCHNG=17
      INTEGER ZCHNG(NCHNG),LCHNG(NCHNG)

      integer ichng

      DATA ZCHNG /5,11,13,19,21,31,37,39,49,55,57,58,72,81,87,89,90/
      DATA LCHNG /1, 0, 1, 0, 2, 1, 0, 2, 1, 0, 2, 3, 2, 1, 0, 2, 3/

      LMXCHM=0
      LMXATM=0
      DO ICHNG=1,NCHNG
         IF (Z.LT.ZCHNG(ICHNG)) RETURN
         LMXCHM=LCHNG(ICHNG)
         LMXATM=MAX(LMXATM,LMXCHM)
      ENDDO

      END subroutine lmxofz


      SUBROUTINE CNFIG (Z,CONFIG) 
C Returns the valence configuration for atomic ground state
C Input: Z Atomic number
C Written by A.R.Williams 

      integer, intent(in)  :: z
      integer, intent(out) :: config(0:)

      integer lmax, nchng
      PARAMETER (LMAX=3,NCHNG=15)
      INTEGER ZCHNG(NCHNG),LCHNG(NCHNG)
      integer l, ichng

*     DATA ZCHNG /3,9,11,17,19,31,35,37,49,53,55,72,81,85,87/
*     DATA ZCHNG /3,11,11,17,19,31,35,37,49,53,55,72,81,85,87/  

      DATA ZCHNG /3,11,11,19,19,31,37,37,49,55,55,72,81,87,87/
      DATA LCHNG /0,0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1/
      DO 10 L=0,LMAX
         CONFIG(L)=L+1
   10 CONTINUE
      DO 20 ICHNG=1,NCHNG
         IF (Z.LT.ZCHNG(ICHNG)) GOTO 30
         L=LCHNG(ICHNG)
         CONFIG(L)=CONFIG(L)+1
   20 CONTINUE
   30 continue

      end subroutine cnfig

      FUNCTION SYMBOL( IZ )
      character*2 symbol
      integer, intent(in) :: iz
C RETURNS THE SYMBOL OF THE ELEMENT OF ATOMIC NUMBER IZ
C Written by J. Soler

      integer, PARAMETER  :: NZ=103
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
      END function symbol

      FUNCTION ATMASS(IZ)
      real(dp)             :: atmass
      integer, intent(in)  :: iz

C Returns the average atomic mass from the atomic number IZ.
C Dta taken from VCH periodic table.
C Written by J.M.Soler. April'97.

      integer, PARAMETER  :: NA=94
      character(len=100) message

      DOUBLE PRECISION AMASS(0:NA)
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


      IF (IZ.LT.0 .OR. IZ.GT.NA) THEN
         write(message,'(a,i4)') 'ATMASS: NO DATA FOR Z =',IZ
         call die(message)
      ELSE
         ATMASS=AMASS(IZ)
      ENDIF
      END function atmass

      end module periodic_table







