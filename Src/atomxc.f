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
      subroutine atomxc( IREL, NR, MAXR, RMESH, nspin, Dens, 
     .                   EX, EC, DX, DC, VXC )

C *******************************************************************
C Finds total exchange-correlation energy and potential for a
C spherical electron density distribution.
C This version implements the Local (spin) Density Approximation and
C the Generalized-Gradient-Aproximation with the 'explicit mesh 
C functional' approach of White & Bird, PRB 50, 4954 (1994).
C Gradients are 'defined' by numerical derivatives, using 2*NN+1 mesh
C   points, where NN is a parameter defined below
C Ref: L.C.Balbas et al, PRB 64, 165110 (2001)
C Wrtten by J.M.Soler using algorithms developed by 
C   L.C.Balbas, J.L.Martins and J.M.Soler, Dec.1996
C ************************* INPUT ***********************************
C CHARACTER*(*) FUNCTL : Functional to be used:
C              'LDA' or 'LSD' => Local (spin) Density Approximation
C                       'GGA' => Generalized Gradient Corrections
C                                Uppercase is optional
C CHARACTER*(*) AUTHOR : Parametrization desired:
C     'CA' or 'PZ' => LSD Perdew & Zunger, PRB 23, 5075 (1981)
C           'PW91' => GGA Perdew & Wang, JCP, 100, 1290 (1994) 
C           'PW92' => LSD Perdew & Wang, PRB, 45, 13244 (1992). This is
C                     the local density limit of the next:
C            'PBE' => GGA Perdew, Burke & Ernzerhof, PRL 77, 3865 (1996)
C           'RPBE' => GGA Hammer, Hansen & Norskov, PRB 59, 7413 (1999)
C         'REVPBE' => GGA Zhang & Yang, PRL 80,890(1998)
C            'LYP' => GGA Becke-Lee-Yang-Parr (see subroutine blypxc)
C            'WC'  => GGA Wu-Cohen (see subroutine wcxc)
C                     Uppercase is optional
C INTEGER IREL         : Relativistic exchange? (0=>no, 1=>yes)
C INTEGER NR           : Number of radial mesh points
C INTEGER MAXR         : Physical first dimension of RMESH, Dens and VXC
C REAL*8  RMESH(MAXR)  : Radial mesh points
C INTEGER nspin        : nspin=1 => unpolarized; nspin=2 => polarized
C REAL*8  Dens(MAXR,nspin) : Total (nspin=1) or spin (nspin=2) electron
C                            density at mesh points
C ************************* OUTPUT **********************************
C REAL*8  EX              : Total exchange energy
C REAL*8  EC              : Total correlation energy
C REAL*8  DX              : IntegralOf( rho * (eps_x - v_x) )
C REAL*8  DC              : IntegralOf( rho * (eps_c - v_c) )
C REAL*8  VXC(MAXR,nspin) : (Spin) exch-corr potential
C ************************ UNITS ************************************
C Distances in atomic units (Bohr).
C Densities in atomic units (electrons/Bohr**3)
C Energy unit depending of parameter EUNIT below
C ********* ROUTINES CALLED *****************************************
C GGAXC, LDAXC
C *******************************************************************

      use precision, only : dp
      use xcmod,     only : nXCfunc, XCfunc, XCauth
      use xcmod,     only : XCweightX, XCweightC
      use sys,       only: die
      use alloc,     only: re_alloc, de_alloc

C Next line is nonstandard but may be suppressed
      implicit none

C Argument types and dimensions
      integer,   intent(in)  :: IREL
      integer,   intent(in)  :: MAXR
      integer,   intent(in)  :: NR
      integer,   intent(in)  :: nspin
      real(dp),  intent(in)  :: Dens(MAXR,nspin)
      real(dp),  intent(in)  :: RMESH(MAXR)
      real(dp),  intent(out) :: VXC(MAXR,nspin)
      real(dp),  intent(out) :: DC
      real(dp),  intent(out) :: DX
      real(dp),  intent(out) :: EC
      real(dp),  intent(out) :: EX

C Internal parameters
C NN    : order of the numerical derivatives: the number of radial 
C          points used is 2*NN+1
C mspin : must be equal or larger than nspin (4 for non-collinear spin)
      integer,   parameter   :: mspin = 4
      integer,   parameter   :: NN = 5

C Fix energy unit:  EUNIT=1.0 => Hartrees,
C                   EUNIT=0.5 => Rydbergs,
C                   EUNIT=0.03674903 => eV
      real(dp),  parameter   :: EUNIT = 0.5_dp

C DVMIN is added to differential of volume to avoid division by zero
      real(dp),  parameter   :: DVMIN = 1.0d-12

C Local variables and arrays
      logical
     .  GGA, GGAfunc, VDWfunc
      integer
     .  IN, IN1, IN2, IR, IS, JN, NF
      real(dp)
     .  D(mspin), DECDD(mspin), DECDGD(3,mspin),
     .  DEXDD(mspin), DEXDGD(3,mspin),
     .  DGDM(-NN:NN), DGIDFJ(-NN:NN), DRDM, DVol, 
     .  DVCDN(mspin,mspin), DVXDN(mspin,mspin),
     .  EPSC, EPSX, F1, F2, GD(3,mspin), PI
      real(dp), pointer :: Aux(:)
      external
     .  GGAXC, LDAXC

C Set GGA switch
      GGA = .false.
      do nf = 1,nXCfunc
        if ( XCfunc(nf).eq.'VDW' .or. XCfunc(nf).eq.'vdw') then
          GGA = .true.
        else if ( XCfunc(nf).eq.'GGA' .or. XCfunc(nf).eq.'gga') then
          GGA = .true.
        else
          if ( XCfunc(nf).ne.'LDA' .and. XCfunc(nf).ne.'lda' .and.
     .         XCfunc(nf).ne.'LSD' .and. XCfunc(nf).ne.'lsd' ) then
            call die('ATOMXC: Unknown functional ' // XCfunc(nf))
          endif 
        endif
      enddo

C Initialize output
      EX = 0.0_dp
      EC = 0.0_dp
      DX = 0.0_dp
      DC = 0.0_dp
      do IS = 1,nspin
        do IR = 1,NR
          VXC(IR,IS) = 0.0_dp
        enddo
      enddo

C Set up workspace array
      if (GGA) then
        nullify( Aux )
        call re_alloc( Aux, 1, NR, name='Aux', routine='atomxc' )
      endif

C Get number pi
      PI = 4.0_dp * ATAN(1.0_dp)

C Loop on mesh points
      do IR = 1,NR

C Find interval of neighbour points to calculate derivatives
        IN1 = MAX(  1, IR-NN ) - IR
        IN2 = MIN( NR, IR+NN ) - IR

C Find weights of numerical derivation from Lagrange
C interpolation formula
        do IN = IN1,IN2
          IF (IN .EQ. 0) THEN
            DGDM(IN) = 0
            do JN = IN1,IN2
              IF (JN.NE.0) DGDM(IN) = DGDM(IN) + 1.D0 / (0 - JN)
            enddo
          ELSE
            F1 = 1
            F2 = 1
            do JN = IN1,IN2
              IF (JN.NE.IN .AND. JN.NE.0) F1 = F1 * (0  - JN)
              IF (JN.NE.IN)               F2 = F2 * (IN - JN)
            enddo
            DGDM(IN) = F1 / F2
          ENDIF
        enddo

C Find dr/dmesh
        DRDM = 0.0_dp
        do IN = IN1,IN2
          DRDM = DRDM + RMESH(IR+IN) * DGDM(IN)
        enddo

C Find differential of volume. Use trapezoidal integration rule
        DVol = 4.0_dp * PI * RMESH(IR)**2 * DRDM
C DVMIN is a small number added to avoid a division by zero
        DVol = DVol + DVMIN
        if (IR.eq.1 .or. IR.eq.NR) DVol = 0.5_dp*DVol
        if (GGA) Aux(IR) = DVol

C Find the weights for the derivative d(gradF(i))/d(F(j)), of
C the gradient at point i with respect to the value at point j
        if (GGA) then
          do IN = IN1,IN2
            DGIDFJ(IN) = DGDM(IN) / DRDM
          enddo
        endif

C Find density and gradient of density at this point
        do IS = 1,nspin
          D(IS) = Dens(IR,IS)
        enddo
        if (GGA) then
          do IS = 1,nspin
            GD(1,IS) = 0.0_dp
            GD(2,IS) = 0.0_dp
            GD(3,IS) = 0.0_dp
            do IN = IN1,IN2
              GD(3,IS) = GD(3,IS) + DGIDFJ(IN) * Dens(IR+IN,IS)
            enddo
          enddo
        endif

C Loop over exchange-correlation functions
        do nf = 1,nXCfunc

C Is this a GGA or VDW?
          if (XCfunc(nf).eq.'VDW' .or. XCfunc(nf).eq.'vdw') then
            VDWfunc = .true.
            GGAfunc = .true.
          else if (XCfunc(nf).eq.'GGA' .or. XCfunc(nf).eq.'gga') then
            VDWfunc = .false.
            GGAfunc = .true.
          else
            VDWfunc = .false.
            GGAfunc = .false.
          endif

C Find exchange and correlation energy densities and their 
C derivatives with respect to density and density gradient
          if (VDWfunc) then
            CALL GGAXC( 'revPBE', IREL, nspin, D, GD,
     .                  EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
          else if (GGAfunc) then
            CALL GGAXC( XCauth(nf), IREL, nspin, D, GD,
     .                  EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
          else
            CALL LDAXC( XCauth(nf), IREL, nspin, D, EPSX, EPSC, DEXDD, 
     .                  DECDD, DVXDN, DVCDN )
          endif

C Scale terms by weights
          EPSX = EPSX*XCweightX(nf)
          EPSC = EPSC*XCweightC(nf)
          DEXDD(1:nspin) = DEXDD(1:nspin)*XCweightX(nf)
          DECDD(1:nspin) = DECDD(1:nspin)*XCweightC(nf)
          if (GGAfunc) then
            DEXDGD(1:3,1:nspin) = DEXDGD(1:3,1:nspin)*XCweightX(nf)
            DECDGD(1:3,1:nspin) = DECDGD(1:3,1:nspin)*XCweightC(nf)
          endif

C Add contributions to exchange-correlation energy and its
C derivatives with respect to density at all points
          do IS = 1,nspin
            EX = EX + DVol*D(IS)*EPSX
            EC = EC + DVol*D(IS)*EPSC
            DX = DX + DVol*D(IS)*(EPSX - DEXDD(IS))
            DC = DC + DVol*D(IS)*(EPSC - DECDD(IS))
            if (GGAfunc) then
              VXC(IR,IS) = VXC(IR,IS) + DVol*(DEXDD(IS) + DECDD(IS))
              do IN = IN1,IN2
                DX= DX - DVol*Dens(IR+IN,IS)*DEXDGD(3,IS)*DGIDFJ(IN)
                DC= DC - DVol*Dens(IR+IN,IS)*DECDGD(3,IS)*DGIDFJ(IN)
                VXC(IR+IN,IS) = VXC(IR+IN,IS) + 
     .            DVol*(DEXDGD(3,IS) + DECDGD(3,IS))*DGIDFJ(IN)
              enddo
            else
              if (GGA) then
                VXC(IR,IS) = VXC(IR,IS) + DVol*(DEXDD(IS) + DECDD(IS))
              else
                VXC(IR,IS) = VXC(IR,IS) + DEXDD(IS) + DECDD(IS)
              endif
            endif
          enddo

        enddo

      enddo

C Divide by volume element to obtain the potential (per electron)
      if (GGA) then
        do IS = 1,NSPIN
          do IR = 1,NR
            DVol = AUX(IR)
            VXC(IR,IS) = VXC(IR,IS) / DVol
          enddo
        enddo
        call de_alloc( aux,  name='aux' )
      endif

C Divide by energy unit
      EX = EX / EUNIT
      EC = EC / EUNIT
      DX = DX / EUNIT
      DC = DC / EUNIT
      do IS = 1,nspin
        do IR = 1,NR
          VXC(IR,IS) = VXC(IR,IS) / EUNIT
        enddo
      enddo

      end
