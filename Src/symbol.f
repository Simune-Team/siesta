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

