C $Id: vmb.f,v 1.5 1999/02/26 14:39:35 wdpgaara Exp $

      subroutine vmb(nat,ttemp,mass,va)
C *************************************************************************
C This routine creates a velocity distribution according to the
C Maxwell-Boltzman distribution with a given temperature.
C It also imposes the constraint of zero total velocity 
C (to avoid center of mass drifts).
C
C Temperatures in Kelvin
C Mass in atomic mass units
C Velocities in bohr/fs
C
C Writen by J. Junquera  Nov. 96,  P. Ordejon  Nov 97.
C ************************** INPUT ****************************************
C integer nat              : Number of atoms
C real*8  ttemp             : Temperature desired for the distribution (K)
C real*8  mass(nat)        : Atomic masses (in amu)
C ************************* OUTPUT ****************************************
C real*8  va(3,nat)        : Atomic Velocities
C *************************************************************************
      implicit none

      integer 
     .  nat

      double precision
     .  mass(nat), ttemp, va(3,nat)

C Internal variables ..................

      integer 
     .  iseed, i, ix

      double precision
     .  massi, tempe, velo, vtot(3)

      external velo
C .....................

      do ix=1,3
        vtot(ix)=0.0
      enddo

      iseed=-17

C Loop over atoms to assing velocities .................
      do i = 1,nat
        massi = mass(i)
           
        va(1,i) = velo(iseed,ttemp,massi)
        va(2,i) = velo(iseed,ttemp,massi)
        va(3,i) = velo(iseed,ttemp,massi)

        vtot(1) = vtot(1) + va(1,i)
        vtot(2) = vtot(2) + va(2,i)
        vtot(3) = vtot(3) + va(3,i)
      enddo
C ...............
        
C Impose constraint on zero center of mass velocity ....................
      do i = 1,Nat
        do ix=1,3
          va(ix,i) = va(ix,i) - vtot(ix)/nat
        enddo
      enddo
C ...............

      if (nat .le. 1) return

C Correct velocity to exact temperature .....................
      call temp(2,Nat,mass,va,tempe)

      if (abs(tempe-ttemp) .ge. 1.e-4 .and. tempe .ge. 1.e-4) then
        do i = 1,Nat
          do ix=1,3
            va(ix,i) = va(ix,i) * dsqrt(ttemp/tempe)
          enddo
        enddo
      endif

      return
      end



      real*8 function velo(iseed,temp,mass)
C *************************************************************************
C This function assigns velocities to atoms according to the 
C Maxwell-Boltzman distribution.
C It generates random numbers drawn from the normal distribution,
C using as input random numbers distributed uniformly from 0 to 1,  
C which are provided by routine ran3.
C
C Temperatures in Kelvin
C Mass in atomic mass units
C Velocities in bohr/fs
C
C Writen by J. Junquera  Nov. 96,  P. Ordejon  Nov 97.
C ************************** INPUT ****************************************
C integer iseed            : Seed for random number generator
C real*8  temp             : Temperature desired for the distribution (K)
C real*8  mass             : Atomic mass (in amu)
C ************************* OUTPUT ****************************************
C real*8  velo             : Velocity component 
C *************************************************************************

      implicit none

      integer 
     . iseed

      real*8
     . mass,temp

C Internal variables .................

      real*8
     .  arg1, arg2, gauss, med, pi, ran3, var
      
      external
     .  ran3
        
C ...........

C  For other distributions med may be different from cero.
      med = 0.0
      pi = acos(-1.0)
      var = sqrt(temp/mass)
C  conversion factor to bohr/fs
      var = var * 0.00172309
 
      arg1 = sqrt((-2.) * log(ran3(iseed)))
        
      arg2 = 2. * pi * ran3(iseed)
      gauss = arg1 * cos(arg2)

      velo = med + var * gauss

      return
      end



      subroutine temp(iunit,natoms,ma,va,tempe)
C *************************************************************************
C Subroutine to calculate instantaneous temperature
C
C Written by P.Ordejon, November'97
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Atomic masses in atomic mass units
C
C Space units depend on input option:
C
C   Option iunit = 1:
C     Distances are in Angstrom
C   Option iunit = 2:
C     Distances are in Bohr
C ************************* INPUT *********************************************
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer natoms        : Number of atoms in the simulation cell
C real*8 ma(natoms)     : Atomic masses 
C real*8 va(3,natoms)   : Atomic velocities 
C ************************* OUTOPUT *******************************************
C real*8 tempe          : Instantaneous system temperature 
C *****************************************************************************
      implicit none

      integer 
     .   natoms,iunit

      double precision
     .  ma(natoms),tempe,va(3,natoms)

C Internal variables ..........................................................
 
      integer
     .  i,ia

      double precision
     .  Ang,eV,fovermp,kin

C ........................

      if (iunit .ne. 1 .and. iunit .ne. 2) then
        write(6,*) 'temp: Wrong iunit option;  must be 1 or 2'
        stop
      endif

C Define constants and conversion factors .....................................

      Ang = 1.d0 / 0.529177d0
      eV  = 1.d0 / 13.60580d0

      if (iunit .eq. 1) then
C  convert F/m in (eV/Amstrong)/amu  to  Amstrong/fs**2
        fovermp = 0.009579038
      else
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038 *Ang**2 / eV
      endif

C ........................

C Calculate kinetic energy and temperature ...................
C Kinetic energy of atoms
      kin = 0.0d0
      do ia = 1,natoms
        do i = 1,3
          kin = kin + 0.5d0 * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        tempe = 2.0d0*kin/(3.0d0*natoms-3.)/8.617d-5
      else
        tempe = 2.0d0*kin/(3.0d0*natoms-3.)/8.617d-5/eV
      endif

C .....................

      return
      end
    
