      module m_fermid

      use precision
      use parallel,   only : Node
      use fdf
#ifdef MPI
      use mpi_siesta
#endif

      public :: fermid, fermispin, stepf
      private :: enpy, derfc, derf, whg, hp

      private

      CONTAINS

      subroutine fermid( nspin, MAXSPN, NK, WK, MAXE, NE, E, 
     .                   TEMP, qtot, WKE, EF, entropy )

C *********************************************************************
C Finds the Fermi energy and the occupation weights of states.
C Written by J.M.Soler. August'96.
C Simple single excitation introduced by E. Artacho August 2002.
C Alternative occupation functions introduced by P. Ordejon, May'03.
C ********** INPUT ****************************************************
C INTEGER nspin    : Number of different spin polarizations (1 or 2)
C INTEGER MAXSPN   : Maximum number of different spin polarizations (1 or 2)
C                    for E and WKE matrices dimensions
C INTEGER NK       : Number of K-points
C REAL*8  WK(NK)   : Sampling weights of k-points (must sum 1)
C INTEGER MAXE     : First dimension of E and WKE
C INTEGER NE       : Number of bands
C REAL*8  E(MAXE,MAXSPN,NK) : State eigenvalues
C REAL*8  TEMP     : Temperature (in the same units of E)
C REAL*8  qtot     : Total valence charge (number of electrons)
C ********** OUTPUT ***************************************************
C REAL*8  WKE(MAXE,MAXSPN,NK) : Occupations multiplied by k-point weights
C                               (sum qtot)
C REAL*8  EF                 : Fermi energy
C REAL*8  entropy            : Entropy contribution to the electronic
C                              Free Energy
C *********************************************************************

C
C  Modules
C

CCCCCCCCCC      implicit real(dp)(A-H,O-Z)

C Passed variables
      integer        :: maxe
      integer        :: maxspn
      integer        :: ne
      integer        :: nk
      integer        :: nspin

      real(dp)       :: e(maxe,maxspn,nk)
      real(dp)       :: entropy
      real(dp)       :: qtot
      real(dp)       :: temp
      real(dp)       :: wke(maxe,maxspn,nk)
      real(dp)       :: wk(nk)
      real(dp)       :: ef

C Local variables
      integer        :: ie
      integer        :: ief
      integer        :: ik
      integer        :: ispin
      integer        :: iter
      integer,  save :: nitmax = 150
      logical,  save :: blread = .false.
      logical,  save :: excitd = .false.
      real(dp), save :: tol = 1.0d-10
      real(dp)       :: sumq, emin, emax, t, drange, wkebuf, w, eik

#ifdef MPI
      integer MPIerror
#endif
      
C Reading whether to do an excited state
      if ( .not. blread) then
        if ( Node .eq. 0 ) then
          excitd = fdf_boolean('SingleExcitation', .false.)
          if ( excitd ) write(6,'(/a)') 
     .        'fermid: Calculating for lowest-exciton excited state'
        endif
#ifdef MPI
        call MPI_Bcast(excitd,1,MPI_logical,0,
     .                 MPI_Comm_World,MPIerror)
#endif
        blread = .true.
      endif   

C Determine Fermi level
      sumq = 0.0d0
      emin = e(1,1,1)
      emax = e(1,1,1)
      do ik = 1,nk
        do ispin = 1,nspin
          do ie = 1,ne
            wke(ie,ispin,ik) = wk(ik)*2.0d0/nspin
            sumq = sumq + wke(ie,ispin,ik)
            emin = min(emin,e(ie,ispin,ik))
            emax = max(emax,e(ie,ispin,ik))
          enddo
        enddo
      enddo

      ef = emax

      if (abs(sumq-qtot).lt.tol) then
        if (excitd) then
          if (Node.eq.0) then
            write (6,'(/a)') 
     .            'Fermid: Bands full, no excitation possible'
          endif
          stop
        else
          return
        endif
      endif
      if (sumq.lt.qtot) then
        if (Node.eq.0) then
          write(6,*) 'Fermid: Not enough states'
          write(6,*) 'Fermid: qtot,sumq=',qtot,sumq
        endif
        stop
      endif
      T = max(temp,1.d-6)
      DRange = T*sqrt(-log(tol*0.01d0))
      emin = emin - DRange
      emax = emax + DRange
      DO 50 iter = 1,nitmax
        ef = 0.5d0*(emin + emax)
        sumq = 0.0d0
        do ik = 1,nk
          do ispin = 1,nspin
            do ie = 1,ne
              wke(ie,ispin,ik) = wk(ik)*
     .            stepf((e(ie,ispin,ik)-ef)/T)*2.0d0/nspin
              sumq = sumq + wke(ie,ispin,ik)
            enddo
          enddo
        enddo

C If the Fermi level was found..................... 
        if (abs(sumq-qtot).lt.tol) then

C If excited state is to be calculated, find the level above Ef for
C k=1 and spin=1
          if (excitd) then
            do ie = 1, ne
              if ( e(ie,1,1) .gt. ef ) then
                ief = ie
                go to 314
              endif
            enddo
  314       continue
C and swap populations (meaningful only for T close to 0):
C if nspin =1 populations of homo and lumo are just swapped for is=1
C if nspin =2 populations of homo and lumo become equal.
            wkebuf = ( wke(ief-1,1,1) - wke(ief,1,1) )/(3-nspin)
            wke(ief,1,1) = wke(ief,1,1) + wkebuf
            wke(ief-1,1,1) = wke(ief-1,1,1) - wkebuf
          endif

C Obtain the electronic entropy

100       entropy=0.0d0

          do ik = 1,nk
            do ie = 1,ne
              do ispin = 1,nspin

                W = (nspin / 2.0d0) * wke(ie,ispin,ik) / wk(ik)
                eik = (e(ie,ispin,ik)-ef) / T

                entropy = entropy + ( 2.0d0 * wk(ik) / nspin ) *
     .                    enpy(eik,w)

              enddo
            enddo
          enddo

          return

        endif

        if (sumq.le.qtot) emin = ef
        if (sumq.ge.qtot) emax = ef
  50  continue

      if (Node.eq.0) then
        write(6,*) 'Fermid: Iteration has not converged.'
        write(6,*) 'Fermid: qtot,sumq=',qtot,sumq
        stop 'Fermid: Iteration has not converged.'
      else
        stop
      endif

      end subroutine fermid

!---------------------------------------------------------


      subroutine fermispin( nspin, MAXSPN, NK, WK, MAXE, NE, E, 
     .                      TEMP, qtot, WKE, EF, entropy )

C *********************************************************************
C Finds the Fermi energy and the occupation weights of states,
C for the case where the total spin of the calculation is fixed
C to a given value.
C Written by J.M.Soler. August'96.
C Version modified for fixed spin configurations: P. Ordejon'03-04
C ********** INPUT ****************************************************
C INTEGER nspin    : Number of different spin polarizations (1 or 2)
C INTEGER MAXSPN   : Maximum number of different spin polarizations (1 or 2)
C                    for E and WKE matrices dimensions
C INTEGER NK       : Number of K-points
C REAL*8  WK(NK)   : Sampling weights of k-points (must sum 1)
C INTEGER MAXE     : First dimension of E and WKE
C INTEGER NE       : Number of bands
C REAL*8  E(MAXE,MAXSPN,NK) : State eigenvalues
C REAL*8  TEMP     : Temperature (in the same units of E)
C REAL*8  qtot(MAXSPN) : Total valence charge (number of electrons)
C                         for each spin component
C ********** OUTPUT ***************************************************
C REAL*8  WKE(MAXE,MAXSPN,NK) : Occupations multiplied by k-point weights
C                               (sum qtot)
C REAL*8  EF(nspin)           : Fermi energy (for each spin, if qtot
C                              is different for each spin component.
C REAL*8  entropy            : Entropy contribution to the electronic
C                              Free Energy
C *********************************************************************

      integer :: maxe, maxspn, nk, nspin, nitmax, ispin
      integer :: ik, ie, ne, iter
      real(dp)   :: entropy

      real(dp)         E(MAXE,MAXSPN,NK),emin(4),emax(4),
     .                 EF(nspin),qtot(MAXSPN),sumq(4),TEMP,
     .                 WKE(MAXE,MAXSPN,NK),WK(NK)
      real(dp)   ::    t, drange, w, eik, tol
      logical          conv
      parameter (TOL=1.0D-10,NITMAX=150)

C Local variables
#ifdef MPI
      integer MPIerror
#endif

      conv = .FALSE.

      DO ispin = 1,nspin
       sumq(ispin)=0.0D0
      ENDDO
      DO ispin = 1,nspin
        emin(ispin)=E(1,ispin,1)
        emax(ispin)=E(1,ispin,1)
      ENDDO
      DO IK=1,NK
        DO ispin=1,nspin
          DO IE=1,NE
            WKE(IE,ispin,IK)=WK(IK)*2.D0/nspin
            sumq(ispin)=sumq(ispin)+WKE(IE,ispin,IK)
            emin(ispin)=MIN(emin(ispin),E(IE,ispin,IK))
            emax(ispin)=MAX(emax(ispin),E(IE,ispin,IK))
          enddo
        enddo
      enddo
      DO ispin=1,nspin
        EF(ispin)=emax(ispin)
      ENDDO
      conv = .TRUE.
      DO ispin = 1,nspin
        IF (ABS(sumq(ispin)-qtot(ispin)).GT.TOL) conv = .FALSE.
      ENDDO

      IF (conv) GOTO 100

      DO ispin = 1,nspin
        IF (sumq(ispin).LT.qtot(ispin)) THEN
          if (Node.eq.0) then
            write(6,*) 'Fermid: Not enough states'
            write(6,*) 'Fermid: ispin,qtot,sumq=',
     .                 ispin,qtot(ispin),sumq(ispin)
          endif
          STOP
        ENDIF
      ENDDO
      T=MAX(TEMP,1.D-6)
      DRANGE=T*SQRT(-LOG(TOL*.01D0))
      DO ispin = 1,nspin
        emin(ispin)=emin(ispin)-DRANGE
        emax(ispin)=emax(ispin)+DRANGE
      ENDDO
      DO 50 ITER=1,NITMAX
        DO ispin = 1,nspin
          EF(ispin)=0.5D0*(emin(ispin)+emax(ispin))
          sumq(ispin)=0.D0
        ENDDO
        DO 40 IK=1,NK
          DO 35 ispin=1,nspin
            DO 30 IE=1,NE
              WKE(IE,ispin,IK)=WK(IK)*
     .             STEPF((E(IE,ispin,IK)-EF(ispin))/T)*2.0/nspin
              sumq(ispin)=sumq(ispin)+WKE(IE,ispin,IK)
  30        CONTINUE
  35      CONTINUE
  40    CONTINUE
        conv = .TRUE.
        DO ispin = 1,nspin
          IF (ABS(sumq(ispin)-qtot(ispin)).GT.TOL) conv = .FALSE.
        ENDDO
        IF (conv) GOTO 100
        DO ispin = 1,nspin
          IF (sumq(ispin).LE.qtot(ispin)) emin(ispin)=EF(ispin)
          IF (sumq(ispin).GE.qtot(ispin)) emax(ispin)=EF(ispin)
        ENDDO
  50  CONTINUE

      if (Node.eq.0) then
        WRITE (6,*) 'Fermid: Iteration has not converged.'
        DO ispin = 1,nspin
          WRITE (6,*) 'Fermid: ispin,qtot,sumq=',
     .               ispin,qtot(ispin),sumq(ispin)
        ENDDO 
        STOP 'Fermid: Iteration has not converged.'
      else
        STOP
      endif
100   CONTINUE
      entropy = 0.0d0

      do IK=1,NK
        do IE=1,NE
          do ispin=1,nspin

            W = (nspin / 2.0D0) * WKE(IE,ispin,IK) / WK(IK)
            EIK = (E(IE,ispin,IK)-EF(ispin)) / T

            entropy = entropy + ( 2.0D0 * WK(IK) / nspin ) *
     .                ENPY(EIK,W)

          enddo
        enddo
      enddo

      return

      end subroutine fermispin


      real(dp) function stepf(X)
      real(dp), intent(in) :: x

      integer  :: i, j
      real(dp) :: pi, a, gauss

C Local variables
      character(len=22), save :: ocf = 'FD'
      integer,           save :: nh = 1
      integer,           save :: ocupfnct
      logical                 :: leqi
      logical,           save :: ocfread = .false.
#ifdef MPI
      integer MPIerror
#endif
      
      parameter (PI = 3.14159265358979D0)

C Reading which electronic occupation function to use -------

      if ( .not. ocfread) then

        if ( Node .eq. 0 ) then
          ocf = fdf_string('OccupationFunction','FD')

          if (leqi(ocf,'FD')) then
            ocupfnct=1
            write(6,'(/a)') 
     .     'stepf: Fermi-Dirac step function'
          else if (leqi(ocf,'MP')) then
            ocupfnct=2
            nh = fdf_integer('OccupationMPOrder',1)
            write(6,'(/a,i2)') 
     .     'stepf: Methfessel-Paxton step function'
            write(6,'(a,i2)') 
     .     '       Using Hermite-Gauss polynomials of order ',nh
          else
            stop 'fermid: Error: Allowed values for OccupationFunction 
     . are FD and MP'
          endif
        endif
#ifdef MPI
        call MPI_Bcast(ocupfnct,1,MPI_integer,0,
     .                 MPI_Comm_World,MPIerror)
        call MPI_Bcast(nh,1,MPI_integer,0,
     .                 MPI_Comm_World,MPIerror)
#endif
        ocfread = .true.
      endif   
c--------------------------------------

C Complementary error function. Ref: Fu & Ho, PRB 28, 5480 (1983)
*     STEPF=DERFC(X)  -  not available

      if (ocupfnct .eq. 1) then

C Fermi-Dirac distribution
        if (x.gt.100.D0) then
          stepf = 0.D0
        elseif (x.lt.-100.d0) then
          stepf = 1.d0
        else
          stepf = 1.d0 / ( 1.d0 + exp(x) )
        endif


      else if (ocupfnct .eq. 2) then

C Improved step function. Ref: Methfessel & Paxton PRB40 (15/Aug/89)
C NH is the order of the Hemite polynomial expansion.

        stepf =  0.5d0 * derfc(x)
        a = 1.0d0/sqrt(pi)

        do i = 1,nh

C Get coefficients in Hermite-Gauss expansion
          A = -A / (I * 4.0d0)

C Get contribution to step function at order I
          gauss = dexp(-x*x)
          J = 2*I -1 
          if (gauss .gt. 1.d-20) stepf = stepf + a * hp(x,j) * gauss

        enddO

      else

        stop 'Stepf: Incorrect step function'

      endif

      end function stepf


      real(dp) function enpy(E,W)
      real(dp), intent(in)  :: e, w

C Computes the contribution of a given state with energy E
C (refered to the Fermi energy, in units of the smearing 
C temperature) and occupation W to the electronic entropy
C P. Ordejon, June 2003

C Local variables
      character(len=22), save :: ocf = 'FD'
      integer,           save :: nh = 1
      integer,           save :: ocupfnct
      logical                 :: leqi
      logical,           save :: ocfread = .false.
      real(dp)                :: tiny, wo, we

#ifdef MPI
      integer MPIerror
#endif
      
      data tiny /1.d-15/

c Reading which electronic occupation function to use -------

      if ( .not. ocfread) then

        if ( Node .eq. 0 ) then
          ocf = fdf_string('OccupationFunction', 'FD')

          if (leqi(ocf,'FD')) then
            ocupfnct = 1
          else if (leqi(ocf,'MP')) then
            ocupfnct = 2
            nh = fdf_integer('OccupationMPOrder',1)
            write(6,'(/a,i2)')
     .     'Stepf: Methfessel-Paxton step function'
            write(6,'(a,i2)')
     .     '       Using Hermite-Gauss polynomials of order ',nh
          else
            stop 'Fermid: Error: Allowed values for OccupationFunction 
     . are FD and MP'
          endif
        endif
#ifdef MPI
        call MPI_Bcast(ocupfnct,1,MPI_integer,0,
     .                  MPI_Comm_World,MPIerror)
        call MPI_Bcast(nh,1,MPI_integer,0,
     .                  MPI_Comm_World,MPIerror)
#endif
        ocfread = .true.
      endif   
c--------------------------------------


      if (ocupfnct .eq. 1) then

C Mermin entropy for the Fermi-Dirac distribution
        wo = max( w, tiny )
        we = 1.0d0 - wo
        we = max( we, tiny )

        enpy = - 1.0 * ( wo*log(wo) + we*log(we) )
     .         

      else if (ocupfnct .eq. 2) then

C Entropy for the Improved step function. 
C Ref: Methfessel & Paxton PRB40 (15/Aug/89)

        enpy = whg(e,nh)

      else

        stop 'Stepf: Incorrect step function'

      endif

      END function enpy


      real(dp) function derfc(X)
      real(dp), intent(in) :: x

C  COMPLEMENTARY ERROR FUNCTION FROM "NUMERICAL RECIPES"
C  NOTE: SINGLE PRECISION ACCURACY

      real(dp) :: z, t

      Z = ABS(X)
      T = 1.0d0/(1.0d0 + 0.5d0*Z)
      DERFC=T*EXP(-(Z*Z)-1.26551223D0+T*(1.00002368D0+T*(0.37409196D0+
     .      T*(0.09678418D0+T*(-0.18628806D0+
     .      T*(0.27886807D0+T*(-1.13520398D0+
     .      T*(1.48851587D0+T*(-0.82215223D0+T*.17087277D0)))))))))
      IF (X.LT.0.D0) DERFC=2.D0-DERFC

      end function derfc


      real(dp) function derf(X)
      real(dp), intent(in) :: x

C  ERROR FUNCTION FROM "NUMERICAL RECIPES"
C  NOTE: SINGLE PRECISION ACCURACY

      real(dp) :: z, t

      Z = ABS(X)
      T = 1.0d0/(1.0d0 + 0.5d0*Z)
      DERF= T*EXP(-(Z*Z)-1.26551223D0+T*(1.00002368D0+T*(0.37409196D0+
     .      T*(0.09678418D0+T*(-0.18628806D0+
     .      T*(0.27886807D0+T*(-1.13520398D0+
     .      T*(1.48851587D0+T*(-0.82215223D0+T*.17087277D0)))))))))
      IF (X.LT.0.D0) DERF=2.D0-DERF

      derf = 1.0d0 - derf
      end function derf


      real(dp) function whg(x,n)
C
C  Computes the factors to get the entropy term 
C  for the Methfessel-Paxton smearing with Hermite
C  polynomials of order N
C
C  P. Ordejon, June '03
C

C Passed variables
      integer        :: n
      real(dp)       :: x

C Local variables
      integer        :: i
      real(dp)       :: a
      real(dp)       :: gauss
      real(dp), save :: pi = 3.14159265358979d0
      real(dp)       :: x2

      x2 = x**2.0d0

C Get coefficients

      a = 1.0d0/sqrt(pi)
      do i = 1,n
        a = - a / (dble(I) * 4.0d0)
      enddo

      gauss = dexp(-x2)
      whg = 0.0d0
      if (gauss .gt. 1.0d-20) whg = 0.5d0*a*hp(x,2*n)*gauss

      return
      end function whg


      real(dp) function hp(x,n)
C 
C  Returns the value of the Hermite polynomial of degree N
C  evaluated at X.
C
C  H_0  (x) = 1
C  H_1  (x) = 2x
C  ...
C  H_n+1(x) = 2 x H_n(x) - 2 n H_n-1(x)
C
C  P. Ordejon, June 2003
C

C Passed variables
      integer  :: n
      real(dp) :: x

C Local variables
      integer  :: i
      real(dp) :: hm1
      real(dp) :: hm2

      if (n .gt. 1000) 
     .     stop 'Fermid: Order of Hermite polynomial too large'


      hp = 1.0d0
      hm2 = 0.0d0
      hm1 = 1.0d0

      do i = 1,n

        hp = 2.0_dp * (x * hm1 - dble(i-1) * hm2)
        hm2 = hm1
        hm1 = hp

      enddo
      
      end function hp

      end module m_fermid
