      module m_dynamics

      use precision
      use parallel,   only : Node
      use m_ioxv,     only : xv_file_read
      use sys,        only : die, stop_flag
      use atomlist,   only : iza
      use units,      only : Ang, eV
      use m_mpi_utils, only : broadcast

      implicit none

      public :: npr, nose, verlet2, pr, anneal
      private

      CONTAINS

      subroutine npr(istep,iunit,natoms,fa,stress,tp,tt,dt,
     .               ma,mn,mpr,ntcon,va,xa,hdot,h,kin,kn,kpr,vn,vpr,
     .               temp,pressin)
C *************************************************************************
C Subroutine for MD simulations with CONTROLLED PRESSURE AND TEMPERATURE.
C The temperature is controlled with a NOSE thermostat.
C The Pressure is controlled with the PARRINELLO-RAHMAN method.
C (See Allen-Tildesley, Computer Simulations of Liquids, pp 227-238
C (Oxford Science Publications), and references therein).
C It allows for VARIABLE CELL SHAPE to accomodate the external        
C pressure.
C The combined use of Nose and Parrinello-Rahman dynamics provide 
C trajectories which sample the isothermal-isobaric ensamble.
C The equations of motion are integrated with a modified Verlet 
C algorithm, with a selfconsistent set of equations for the 
C Nose and Parrinello-Rahman variables.
C
C Written by P.Ordejon, November'96
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Time in femtoseconds
C Atomic masses in atomic mass units
C
C Other units depend on input option:
C
C   Option iunit = 1:
C     Pressures are in eV/A**3 ( = 160.20506 GPa)
C     Energies are in eV
C     Distances are in Angstrom
C     Nose and PR masses in eV * fs**2
C   Option iunit = 2:
C     Pressures are in Ry/Bohr**3 
C     Energies are in Ry
C     Distances are in Bohr
C     Nose and PR masses in Ry * fs**2
C ************************* INPUT *********************************************
C integer istep         : Number of time step during simulation
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer natoms        : Number of atoms in the simulation cell
C real*8 fa(3,natoms)   : Atomic forces
C real*8 stress(3,3)    : Stress tensor components
C real*8 tp             : Target pressure
C real*8 tt             : Target temperature
C real*8 dt             : Length of the time step
C real*8 ma(natoms)     : Atomic masses
C real*8 mn             : Mass of Nose thermostat
C real*8 mpr            : Mass of Parrinello-Rahman variables
C integer ntcon         : Total number of position constraints imposed
C ******************* INPUT AND OUTPUT ****************************************
C real*8 va(3,natoms)   : Atomic velocities
C real*8 xa(3,natoms)   : Atomic coordinates
C                        (input: current time step; output: next time step)
C real*8 hdot(3,3)      : Matrix of time derivatives of
C                         the vectors defining the unit cell
C real*8 h(3,3)         : Matrix of the vectors defining the unit
C                         cell at the current time step
C                         h(i,j) is the ith component of jth basis vector
C                        (input: current time step; output: next time step)
C ************************* OUTOPUT *******************************************
C real*8 kin            : Kinetic energy of the atomic system
C real*8 kn             : Kinetic energy of Nose variable
C real*8 kpr            : Kinetic energy of Parrinello-Rahman variable
C real*8 vn             : Potential energyy of Nose variable
C real*8 vpr            : Potential energyy of P-R variables
C real*8 temp           : Instantaneous system temperature
C real*8 pressin        : Instantaneous system pressure
C *****************************************************************************
C

      implicit none

      integer 
     .   natoms, ntcon, istep, iunit

      real(dp)
     .  dt,fa(3,natoms),h(3,3),hdot(3,3),kin,kn,kpr,
     .  ma(natoms),mn,mpr,stress(3,3),tp,tt,
     .  va(3,natoms),vn,vpr,xa(3,natoms)

C Internal variables .............................................................

      integer
     .  ct,i,ia,info,j,k

      real(dp)
     .  aux1(3,3),aux2(3,3),diff,dt2,dtby2,
     .  f(3,3),fi(3,3),fovermp,
     .  g(3,3),gdot(3,3),gi(3,3),
     .  hi(3,3),hnew(3,3),hlast(3,3),hold(3,3),hs(3),
     .  pgas,press(3,3),pressin,
     .  tdiff,tekin,temp,tol,twodt,
     .  vol,volcel,
     .  x,xdot,xlast,xnew,xold

      real(dp), dimension(:,:), allocatable, save ::
     .  s,sdot,snew,sold,sunc,suncdot

      save x,xold,hold

      external
     .  volcel, memory
C ......................................................................

      if (iunit .ne. 1 .and. iunit .ne. 2) then
        if (Node.eq.0) then
          write(6,*) 'npr: Wrong iunit option;  must be 1 or 2'
        endif
        stop
      endif

C Allocate local memory
      allocate(s(3,natoms))
      call memory('A','D',3*natoms,'npr')
      allocate(sdot(3,natoms))
      call memory('A','D',3*natoms,'npr')
      allocate(snew(3,natoms))
      call memory('A','D',3*natoms,'npr')
      allocate(sunc(3,natoms))
      call memory('A','D',3*natoms,'npr')
      allocate(suncdot(3,natoms))
      call memory('A','D',3*natoms,'npr')
      if (.not.allocated(sold)) then
        allocate(sold(3,natoms))
        call memory('A','D',3*natoms,'npr')
      endif

      ct = 3 + ntcon
      if (natoms .eq. 1) ct = 0

C Define constants and conversion factors .......................................
      dt2   = dt**2
      dtby2 = dt/2.0d0
      twodt = dt*2.0d0
      tol   = 1.0d-12

      if (iunit .eq. 1) then
C  convert target temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in eV if Temp is in Kelvin)
        tekin = 0.5d0 * (3.d0 * natoms - ct) * 8.617d-5 * tt
C  convert F/m in (eV/Angstrom)/amu  to  Angstrom/fs**2
        fovermp = 0.009579038
      else
C  convert target temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in Ry if Temp is in Kelvin)
        tekin = eV * 0.5d0 * (3.d0 * natoms - ct) * 8.617d-5 * tt
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038 * Ang**2 / eV
      endif

C  calculate cell volume at current time
      vol = volcel( h )
C ........................

C Compute Parrinello-Rahman variables (H and scaled coordinates) ............
C Compute G=HtH at current time 
      do i = 1,3
        do j = 1,3
          g(i,j) = 0.0d0
          do k = 1,3
            g(i,j) = g(i,j) + h(k,i) * h(k,j)
          enddo
        enddo
      enddo

C Compute Inverse of H and G at current time 
      call inver(h,hi,3,3,info)
      if (info .ne. 0) stop 'npr: INVER failed'
      call inver(g,gi,3,3,info)
      if (info .ne. 0) stop 'npr: INVER failed'

C Calculate scaled coordinates (referred to matrix H) at current time
      do ia = 1,natoms
        do i = 1,3
          s(i,ia) = 0.0d0
          do j = 1,3
            s(i,ia) = s(i,ia) + hi(i,j) * xa(j,ia)
          enddo
        enddo
      enddo

C Initialize variables if current time step is the first of the simulation
      if (istep .eq. 1) then
        x = 0.0
        xold = 0.0
        do i = 1,3
          do j = 1,3
            hold(i,j) = h(i,j) - dt*hdot(i,j)
          enddo
        enddo
        do ia = 1,natoms
          do i = 1,3
            sold(i,ia) = 0.0d0
            do j = 1,3
              sold(i,ia) = sold(i,ia) + hi(i,j)*(xa(j,ia)-dt*va(j,ia)
     .                     + (dt2/2.0d0) * fovermp * fa(j,ia) / ma(ia))
            enddo
          enddo
        enddo
      endif
C ..................

C Compute uncorrected next positions .....................................
      do ia = 1,natoms
        do i = 1,3
          sunc(i,ia) = -sold(i,ia) + 2.0d0 * s(i,ia)
          do k = 1,3
            sunc(i,ia) = sunc(i,ia) + 
     .                   dt2 * hi(i,k) * fovermp * fa(k,ia) / ma(ia)
          enddo
        enddo
      enddo
C ...................

C Compute initial guess for Nose and Parrinello-Rahman 
C   variables at next time step ...........................................
      xnew = 2.0d0 * x - xold
      do j = 1,3
        do i = 1,3
          hnew(i,j) = 2.0d0 * h(i,j) - hold(i,j)
        enddo
      enddo
C ...................

C Start selfconsistency loop to calculate Nose and P-R variables ..........
10    continue

      xlast = xnew
      do j = 1,3
        do i = 1,3
          hlast(i,j) = hnew(i,j)
        enddo
      enddo
        
C xdot and hdot (time derivatives at current time), and related stuff
      xdot = (xnew - xold) / twodt
      do j = 1,3
        do i = 1,3
          hdot(i,j) = (hnew(i,j) - hold(i,j)) / twodt
        enddo
      enddo

      do j = 1,3
        do i = 1,3
          gdot(i,j) = 0.0d0
          do k = 1,3
            gdot(i,j) = gdot(i,j) + h(k,i) * hdot(k,j)
     .                            + hdot(k,i) * h(k,j)
          enddo
        enddo
      enddo

      do j = 1,3
        do i = 1,3
          f(i,j) = 0.0d0
          do k = 1,3
             f(i,j) = f(i,j) + gi(i,k) * gdot(k,j)
          enddo
        enddo
      enddo

      do i = 1,3
        f(i,i) = f(i,i) + xdot
      enddo

      do j = 1,3
        do i = 1,3
          f(i,j) = dtby2 * f(i,j)
        enddo
      enddo

      do i = 1,3
        do j = 1,3
          aux1(i,j) = f(i,j)
        enddo
      enddo

      do i = 1,3
        aux1(i,i) = aux1(i,i) + 1.0d0
      enddo

      call inver(aux1,fi,3,3,info)
      if (info .ne. 0) stop 'npr: INVER failed'

      do j = 1,3
        do i = 1,3
          fi(i,j) = fi(i,j) / twodt
        enddo
      enddo

C Calculate corrected velocities at current time
      do ia = 1,natoms
        do i = 1,3
          sdot(i,ia) = 0.0d0
          do j = 1,3
            sdot(i,ia) = sdot(i,ia) + fi(i,j)*(sunc(j,ia) - sold(j,ia))
          enddo
        enddo
      enddo

C Calculate pressure tensor at current time and ideal gas pressure
      do i = 1,3
        do j = 1,3
          press(i,j) = 0.0d0
        enddo
      enddo
      do ia = 1,natoms
        do i = 1,3
          hs(i) = 0.0d0
          do j = 1,3
            hs(i) = hs(i) + h(i,j) * sdot(j,ia)
          enddo
        enddo
        do j = 1,3
          do i = 1,3
            press(i,j) = press(i,j) + ma(ia) * hs(i) * hs(j) / fovermp
          enddo
        enddo
      enddo
      pgas = 0.0d0
      do i = 1,3
        pgas = pgas + press(i,i) / vol
      enddo
      pgas = pgas / 3.0d0
      do j = 1,3
        do i = 1,3
          press(i,j) = press(i,j) / vol - stress(i,j)
        enddo
      enddo

C Compute internal pressure  (pressin = 1/3 Tr (press))   at current time
      pressin = 0.0
      do i = 1,3
        pressin = pressin + press(i,i)
      enddo
      pressin = pressin / 3.0d0

C  Compute Nose and Parrinello-Rahman variables for next time step 
      xnew = 2.0d0 * x - xold 
     .       + (dt2 / mn) * (3.0d0 * vol * pgas - 2.0d0 * tekin)

      do j = 1,3
        do i = 1,3
          aux1(i,j) = 0.0d0
          aux2(i,j) = 0.0d0
        enddo
      enddo
      do i = 1,3
        aux1(i,i) = -tp
      enddo
      do j = 1,3
        do i = 1,3
          aux1(i,j) = aux1(i,j) + press(i,j)
        enddo
      enddo
      do j = 1,3
        do i = 1,3
          aux2(i,j ) = 0.0d0
          do k = 1,3
             aux2(i,j) = aux2(i,j) + aux1(i,k) * hi(j,k)
          enddo
        enddo
      enddo

      do j = 1,3
        do i = 1,3
          hnew(i,j) = (2.0d0 * h(i,j) 
     .                + (dt2 * exp(2 * x) * vol / mpr) * aux2(i,j)
     .                - (1.0d0 + dtby2 * xdot) * hold(i,j))
     .                / (1.0d0 - dtby2 * xdot)
        enddo
      enddo

C Check if selfconsistency has been reached
      diff = abs(xnew - xlast)
      if (xlast .eq. 0.0d0) then
        if (diff .gt. tol)  goto 10
      else
        if (diff/abs(xlast) .gt. tol)  goto 10
      endif

      diff = 0.0d0
      tdiff = 0.0d0
      do j = 1,3
        do i = 1,3
          diff = diff + abs(hnew(i,j) - hlast(i,j))
          tdiff = tdiff + abs(hlast(i,j))
        enddo
      enddo
      if (tdiff .eq. 0.0d0) then
        if (diff .gt. tol) goto 10
      else
        if (diff/tdiff .gt. tol) goto 10
      endif
C ...................

C Calculate corrected atomic coordinates at next time step ................
      do ia = 1,natoms
        do i = 1,3
          suncdot(i,ia) = 0.0d0
          do j = 1,3
            suncdot(i,ia) = suncdot(i,ia) + f(i,j) * sold(j,ia)
          enddo
        enddo
      enddo
      do ia = 1,natoms
        do i = 1,3
          snew(i,ia) = 0.0d0
          do j = 1,3
            snew(i,ia) = snew(i,ia) + twodt * fi(i,j)
     .                                * (sunc(j,ia) + suncdot(j,ia))
          enddo
        enddo
      enddo

C Save current atomic positions as old ones, 
C   and next positions as current ones
      do i = 1,3
        do ia = 1,natoms
          sold(i,ia) = s(i,ia)
          s(i,ia) = snew(i,ia)
        enddo
      enddo

      do j = 1,3
        do i = 1,3
          hold(i,j) = h(i,j) 
          h(i,j) = hnew(i,j)
        enddo
      enddo

      xold = x
      x = xnew

C Transform back to absolute coordinates 
      do ia = 1,natoms
        do i = 1,3
          xa(i,ia) = 0.0d0
          va(i,ia) = 0.0d0
          do j = 1,3
            xa(i,ia) = xa(i,ia) + h(i,j) * s(j,ia)
            va(i,ia) = va(i,ia) + h(i,j) * sdot(j,ia)
          enddo
        enddo
      enddo
C ....................

C Calculate Kinetic and potential energies ................................
C Kinetic energy of atoms
      kin = (3.0d0 / 2.0d0) * pgas * vol

C Kinetic energy of Nose variable
      kn = (1.0d0 / 2.0d0) * mn * xdot**2

C Kinetic energy of Parrinello-Rahman variables
      kpr = 0.0d0
      do i = 1,3
        do j = 1,3
          kpr = kpr + hdot(j,i)**2
        enddo
      enddo
      kpr = (1.0d0 / 2.0d0) * mpr * kpr / exp(2. * x)

C Potential energy of Nose variable
      vn = 2.0d0 * tekin * xold

C Potential energy of Parrinello-Rahman variables
      vpr = tp * vol

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        temp = 3.0d0*vol*pgas/(3.0d0*natoms-ct)/8.617d-5
      else
        temp = 3.0d0*vol*pgas/(3.0d0*natoms-ct)/8.617d-5/eV
      endif

C .....................

C Deallocate local memory
      call memory('D','D',size(s),'npr')
      deallocate(s)
      call memory('D','D',size(sdot),'npr')
      deallocate(sdot)
      call memory('D','D',size(snew),'npr')
      deallocate(snew)
      call memory('D','D',size(sunc),'npr')
      deallocate(sunc)
      call memory('D','D',size(suncdot),'npr')
      deallocate(suncdot)

      end subroutine npr
    
      subroutine pr(istep,iunit,iquench,natoms,fa,stress,tp,dt,
     .               ma,mpr,ntcon,va,xa,hdot,h,kin,kpr,vpr,
     .               temp,pressin)
C *************************************************************************
C Subroutine for MD simulations with CONTROLLED PRESSURE.
C The Pressure is controlled with the PARRINELLO-RAHMAN method.
C (See Allen-Tildesley, Computer Simulations of Liquids, pp 227-238
C (Oxforf Science Publications), and references therein).
C It allows for VARIABLE CELL SHAPE to accomodate the external
C pressure.
C The use of Parrinello-Rahman dynamics provide 
C trajectories which sample the constant NPE ensamble.
C The equations of motion are integrated with a modified Verlet 
C algorithm, with a selfconsistent set of equations for the 
C Parrinello-Rahman variables.
C
C Written by P.Ordejon, November'96
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Time in femtoseconds
C Atomic masses in atomic mass units
C
C Other units depend on input option:
C
C   Option iunit = 1:
C     Pressures are in eV/A**3 ( = 160.20506 GPa)
C     Energies are in eV
C     Distances are in Angstrom
C     PR mass in eV * fs**2
C   Option iunit = 2:
C     Pressures are in Ry/Bohr**3 
C     Energies are in Ry
C     Distances are in Bohr
C     PR mass in Ry * fs**2
C ************************* INPUT *********************************************
C integer istep         : Number of time step during simulation
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer iquench       : Option for quenching:   
C                              0 = no quenching (standard dynamics)
C                              1 = power quenching (set to cero velocity 
C                                 components opposite to force)
C integer natoms        : Number of atoms in the simulation cell
C real*8 fa(3,natoms)   : Atomic forces 
C real*8 stress(3,3)    : Stress tensor components 
C real*8 tp             : Target pressure
C real*8 dt             : Length of the time step 
C real*8 ma(natoms)     : Atomic masses
C real*8 mpr            : Mass of Parrinello-Rahman variables
C integer ntcon         : Total number of position constraints imposed
C ******************* INPUT AND OUTPUT ****************************************
C real*8 va(3,natoms)   : Atomic velocities
C real*8 xa(3,natoms)   : Atomic coordinates
C                        (input: current time step; output: next time step)
C real*8 hdot(3,3)      : Matrix of time derivatives of
C                         the vectors defining the unit cell
C real*8 h(3,3)         : Matrix of the vectors defining the unit
C                         cell at the current time step 
C                         h(i,j) is the ith component of jth basis vector
C                        (input: current time step; output: next time step)
C ************************* OUTOPUT *******************************************
C real*8 kin            : Kinetic energy of the atomic system 
C real*8 kpr            : Kinetic energy of Parrinello-Rahman var 
C real*8 vpr            : Potential energyy of P-R variables 
C real*8 temp           : Instantaneous system temperature 
C real*8 pressin        : Instantaneous system pressure 
C *****************************************************************************

      implicit none

      integer 
     .   natoms,ntcon,istep,iquench,iunit

      real(dp)
     .  dt,fa(3,natoms),h(3,3),hdot(3,3),kin,kpr,
     .  ma(natoms),mpr,stress(3,3),tp,
     .  va(3,natoms),vpr,xa(3,natoms)

C Internal variables .............................................................

      integer
     .  ct,i,info,ia,j,k

      real(dp)
     .  a1,a2,Ang,aux1(3,3),aux2(3,3),diff,dot,dt2,dtby2,
     .  f(3,3),fi(3,3),fovermp,
     .  g(3,3),gdot(3,3),gi(3,3),
     .  hi(3,3),hnew(3,3),hlast(3,3),hold(3,3),hs(3),
     .  pgas,press(3,3),pressin,
     .  tdiff,temp,tol,twodt,
     .  vol,volcel

      real(dp), dimension(:,:), allocatable, save ::
     .  s,sdot,snew,sold,sunc,suncdot

      save hold

      external
     .  volcel, memory
C ...............................................................................

      if (iunit .ne. 1 .and. iunit .ne. 2) then
        if (Node.eq.0) then
          write(6,*) 'pr: Wrong iunit option;  must be 1 or 2'
        endif
        stop
      endif
      ct = 3 + ntcon
      if (natoms .eq. 1) ct = 0

C Allocate local memory
      allocate(s(3,natoms))
      call memory('A','D',3*natoms,'pr')
      allocate(sdot(3,natoms))
      call memory('A','D',3*natoms,'pr')
      allocate(snew(3,natoms))
      call memory('A','D',3*natoms,'pr')
      allocate(sunc(3,natoms))
      call memory('A','D',3*natoms,'pr')
      allocate(suncdot(3,natoms))
      call memory('A','D',3*natoms,'pr')
      if (.not.allocated(sold)) then
        allocate(sold(3,natoms))
        call memory('A','D',3*natoms,'pr')
      endif

C Define constants and conversion factors .......................................
      dt2   = dt**2
      dtby2 = dt/2.0d0
      twodt = dt*2.0d0
      tol   = 1.0d-12

      if (iunit .eq. 1) then
C  convert F/m in (eV/Angstrom)/amu  to  Angstrom/fs**2
        fovermp = 0.009579038
      else
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038 * Ang**2 / eV
      endif
C  calculate cell volume at current time
      vol = volcel( h )
C ........................

C Compute Parrinello-Rahman variables (H and scaled coordinates) ............
C Compute G=HtH at current time 
      do i = 1,3
        do j = 1,3
          g(i,j) = 0.0d0
          do k = 1,3
            g(i,j) = g(i,j) + h(k,i) * h(k,j)
          enddo
        enddo
      enddo

C Compute Inverse of H and G at current time 
      call inver(h,hi,3,3,info)
      if (info .ne. 0) stop 'pr: INVER failed'
      call inver(g,gi,3,3,info)
      if (info .ne. 0) stop 'pr: INVER failed'

C Calculate scaled coordinates (referred to matrix H) at current time
      do ia = 1,natoms
        do i = 1,3
          s(i,ia) = 0.0d0
          do j = 1,3
            s(i,ia) = s(i,ia) + hi(i,j) * xa(j,ia)
          enddo
        enddo
      enddo

C Initialize variables if current time step is the first of the simulation
      if (istep .eq. 1) then
        do i = 1,3
          do j = 1,3
            hold(i,j) = h(i,j) - dt*hdot(i,j)
          enddo
        enddo
        do ia = 1,natoms
          do i = 1,3
            sold(i,ia) = 0.0d0
            do j = 1,3
              sold(i,ia) = sold(i,ia) + hi(i,j)*(xa(j,ia)-dt*va(j,ia)
     .                     + (dt2/2.0d0) * fovermp * fa(j,ia) / ma(ia))
            enddo
          enddo
        enddo
      endif
C ..................

C Compute uncorrected next positions .....................................
      do ia = 1,natoms
        do i = 1,3
          sunc(i,ia) = -sold(i,ia) + 2.0d0 * s(i,ia)
          do k = 1,3
            sunc(i,ia) = sunc(i,ia) + 
     .                   dt2 * hi(i,k) * fovermp * fa(k,ia) / ma(ia)
          enddo
        enddo
      enddo
C ...................

C Compute initial guess for Parrinello-Rahman 
C   variables at next time step ...........................................
      do j = 1,3
        do i = 1,3
          hnew(i,j) = 2.0d0 * h(i,j) - hold(i,j)
        enddo
      enddo
C ...................

C Start selfconsistency loop to calculate P-R variables ..........
10    continue

      do j = 1,3
        do i = 1,3
          hlast(i,j) = hnew(i,j)
        enddo
      enddo
        
C hdot (time derivatives at current time), and related stuff
      do j = 1,3
        do i = 1,3
          hdot(i,j) = (hnew(i,j) - hold(i,j)) / twodt
        enddo
      enddo

      do j = 1,3
        do i = 1,3
          gdot(i,j) = 0.0d0
          do k = 1,3
            gdot(i,j) = gdot(i,j) + h(k,i) * hdot(k,j)
     .                            + hdot(k,i) * h(k,j)
          enddo
        enddo
      enddo

      do j = 1,3
        do i = 1,3
          f(i,j) = 0.0d0
          do k = 1,3
             f(i,j) = f(i,j) + gi(i,k) * gdot(k,j)
          enddo
        enddo
      enddo

      do j = 1,3
        do i = 1,3
          f(i,j) = dtby2 * f(i,j)
        enddo
      enddo

      do i = 1,3
        do j = 1,3
          aux1(i,j) = f(i,j)
        enddo
      enddo

      do i = 1,3
        aux1(i,i) = aux1(i,i) + 1.0d0
      enddo

      call inver(aux1,fi,3,3,info)
      if (info .ne. 0) stop 'pr: INVER failed'

      do j = 1,3
        do i = 1,3
          fi(i,j) = fi(i,j) / twodt
        enddo
      enddo

C Calculate corrected velocities at current time
      do ia = 1,natoms
        do i = 1,3
          sdot(i,ia) = 0.0d0
          do j = 1,3
            sdot(i,ia) = sdot(i,ia) + fi(i,j)*(sunc(j,ia) - sold(j,ia))
          enddo
        enddo
      enddo

C Calculate pressure tensor at current time and ideal gas pressure
      do i = 1,3
        do j = 1,3
          press(i,j) = 0.0d0
        enddo
      enddo
      do ia = 1,natoms
        do i = 1,3
          hs(i) = 0.0d0
          do j = 1,3
            hs(i) = hs(i) + h(i,j) * sdot(j,ia)
          enddo
        enddo
        do j = 1,3
          do i = 1,3
            press(i,j) = press(i,j) + ma(ia) * hs(i) * hs(j) / fovermp
          enddo
        enddo
      enddo
      pgas = 0.0d0
      do i = 1,3
        pgas = pgas + press(i,i) / vol
      enddo
      pgas = pgas / 3.0d0
      do j = 1,3
        do i = 1,3
          press(i,j) = press(i,j) / vol - stress(i,j)
        enddo
      enddo

C Compute internal pressure  (pressin = 1/3 Tr (press))   at current time
      pressin = 0.0
      do i = 1,3
        pressin = pressin + press(i,i)
      enddo
      pressin = pressin / 3.0d0

C  Compute Parrinello-Rahman variables for next time step 
      do j = 1,3
        do i = 1,3
          aux1(i,j) = 0.0d0
          aux2(i,j) = 0.0d0
        enddo
      enddo
      do i = 1,3
        aux1(i,i) = -tp
      enddo
      do j = 1,3
        do i = 1,3
          aux1(i,j) = aux1(i,j) + press(i,j)
        enddo
      enddo
      do j = 1,3
        do i = 1,3
          aux2(i,j ) = 0.0d0
          do k = 1,3
             aux2(i,j) = aux2(i,j) + aux1(i,k) * hi(j,k)
          enddo
        enddo
      enddo

      do j = 1,3
        do i = 1,3
          hnew(i,j) = (2.0d0 * h(i,j) 
     .                + (dt2 * vol / mpr) * aux2(i,j)
     .                - hold(i,j))
        enddo
      enddo

C Check if selfconsistency has been reached
      diff = 0.0d0
      tdiff = 0.0d0
      do j = 1,3
        do i = 1,3
          diff = diff + abs(hnew(i,j) - hlast(i,j))
          tdiff = tdiff + abs(hlast(i,j))
        enddo
      enddo
      if (tdiff .eq. 0.0d0) then
        if (diff .gt. tol) goto 10
      else
        if (diff/tdiff .gt. tol) goto 10
      endif
C ...................

C Calculate corrected atomic coordinates at next time step ................
      do ia = 1,natoms
        do i = 1,3
          suncdot(i,ia) = 0.0d0
          do j = 1,3
            suncdot(i,ia) = suncdot(i,ia) + f(i,j) * sold(j,ia)
          enddo
        enddo
      enddo
      do ia = 1,natoms
        do i = 1,3
          snew(i,ia) = 0.0d0
          do j = 1,3
            snew(i,ia) = snew(i,ia) + twodt * fi(i,j)
     .                                * (sunc(j,ia) + suncdot(j,ia))
          enddo
        enddo
      enddo

C Quench option if iquench = 0 ..............................................

C Quench velocity components going uphill
      if (iquench .eq. 1) then
        do ia = 1,natoms
          do i = 1,3
            a1 = 0.0d0
            a2 = 0.0d0
            do j = 1,3
              a1 = a1 + hi(i,j) * fovermp * fa(j,ia) / ma(ia) 
              a2 = a2 - f(i,j) * sdot(j,ia) / dtby2
            enddo
            dot = a1 * sdot(i,ia)
            if (dot .lt. 0.0) then
              sdot(i,ia) = 0.0
              snew(i,ia) = s(i,ia)
            endif
          enddo
        enddo
  
        do i = 1,3
          do j = 1,3
            dot = hdot(i,j) * aux2(i,j)
            if (dot .le. 0.0) then
              hdot(i,j) = 0.0
              hnew(i,j) = h(i,j)
            endif
          enddo
        enddo
            
C Compute gas pressure again, in case quench has happened
        do i = 1,3
          do j = 1,3
            press(i,j) = 0.0d0
          enddo
        enddo
        do ia = 1,natoms
          do i = 1,3
            hs(i) = 0.0d0
            do j = 1,3
              hs(i) = hs(i) + h(i,j) * sdot(j,ia)
            enddo
          enddo
          do j = 1,3
            do i = 1,3
              press(i,j) = press(i,j) + ma(ia) * hs(i) * hs(j) / fovermp
            enddo
          enddo
        enddo
        pgas = 0.0d0
        do i = 1,3
          pgas = pgas + press(i,i) / vol
        enddo

      endif
C ....................
          
C Save current atomic positions as old ones, 
C   and next positions as current ones
      do i = 1,3
        do ia = 1,natoms
          sold(i,ia) = s(i,ia)
          s(i,ia) = snew(i,ia)
        enddo
      enddo

      do j = 1,3
        do i = 1,3
          hold(i,j) = h(i,j) 
          h(i,j) = hnew(i,j)
        enddo
      enddo

C Transform back to absolute coordinates 
      do ia = 1,natoms
        do i = 1,3
          xa(i,ia) = 0.0d0
          va(i,ia) = 0.0d0
          do j = 1,3
            xa(i,ia) = xa(i,ia) + h(i,j) * s(j,ia)
            va(i,ia) = va(i,ia) + h(i,j) * sdot(j,ia)
          enddo
        enddo
      enddo
C ....................

C Calculate Kinetic and potential energies ................................
C Kinetic energy of atoms
      kin = (3.0d0 / 2.0d0) * pgas * vol

C Kinetic energy of Parrinello-Rahman variables
      kpr = 0.0d0
      do i = 1,3
        do j = 1,3
          kpr = kpr + hdot(j,i)**2
        enddo
      enddo
      kpr = (1.0d0 / 2.0d0) * mpr * kpr 

C Potential energy of Parrinello-Rahman variables
      vpr = tp * vol

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        temp = 3.0d0*vol*pgas/(3.0d0*natoms-ct)/8.617d-5
      else
        temp = 3.0d0*vol*pgas/(3.0d0*natoms-ct)/8.617d-5/eV
      endif
C .....................

C Deallocate local memory
      call memory('D','D',size(s),'pr')
      deallocate(s)
      call memory('D','D',size(sdot),'pr')
      deallocate(sdot)
      call memory('D','D',size(snew),'pr')
      deallocate(snew)
      call memory('D','D',size(sunc),'pr')
      deallocate(sunc)
      call memory('D','D',size(suncdot),'pr')
      deallocate(suncdot)

      end subroutine pr
    
      subroutine nose(istep,iunit,natoms,fa,tt,dt,
     .                ma,mn,ntcon,va,xa,kin,kn,vn,
     .                temp)
C *************************************************************************
C Subroutine for MD simulations with CONTROLLED TEMPERATURE.
C The temperature is controlled with a NOSE thermostat.
C The use of Nose dynamics provides trajectories which sample the 
C isothermal ensamble.
C The equations of motion are integrated with a modified Verlet 
C algorithm, with a selfconsistent set of equations for the 
C Nose variables.
C
C Written by P.Ordejon, November'96
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Time in femtoseconds
C Atomic masses in atomic mass units
C
C Other units depend on input option:
C
C   Option iunit = 1:
C     Energies are in eV
C     Distances are in Angstrom
C     Nose mass in eV * fs**2
C   Option iunit = 2:
C     Energies are in Ry
C     Distances are in Bohr
C     Nose mass in Ry * fs**2
C ************************* INPUT *********************************************
C integer istep         : Number of time step during simulation
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer natoms        : Number of atoms in the simulation cell
C real*8 fa(3,natoms)   : Atomic forces
C real*8 tt             : Target temperature
C real*8 dt             : Length of the time step 
C real*8 ma(natoms)     : Atomic masses
C real*8 mn             : Mass of Nose thermostat
C integer ntcon         : Total number of position constraints imposed
C real*8 va(3,natoms)   : Atomic velocities
C                         (used only if istep = 1)
C ******************* INPUT AND OUTPUT ****************************************
C real*8 xa(3,natoms)   : Atomic coordinates
C                        (input: current time step; output: next time step)
C ************************* OUTOPUT *******************************************
C real*8 kin            : Kinetic energy of the atomic system
C real*8 kn             : Kinetic energy of Nose variable 
C real*8 vn             : Potential energyy of Nose var
C real*8 temp           : Instantaneous system temperature 
C *****************************************************************************
C

      integer 
     .   natoms,ntcon,istep,iunit

      real(dp)
     .  dt,fa(3,natoms),kin,kn,
     .  ma(natoms),mn,tt,
     .  va(3,natoms),vn,xa(3,natoms)

      external
     .  memory
C Internal variables .........................................................

      integer
     .  ct,i,ia

      integer  :: iacc, dummy_iza, old_natoms
      real(dp) :: old_dt

      save x,xold

      real(dp)
     .  diff,dt2,dtby2,fact,fovermp,
     .  tekin,temp,tol,twodt,
     .  x,xdot,xlast,xnew,xold

      real(dp), dimension(:,:), allocatable, save ::
     .  xanew,xaold
C .............................................................................

      if (iunit .ne. 1 .and. iunit .ne. 2) then
        if (Node.eq.0) then
          write(6,*) 'nose: Wrong iunit option;  must be 1 or 2'
        endif
        stop
      endif
      ct = 3 + ntcon
      if (natoms .eq. 1) ct = 0

C Allocate local memory and initialise
      if (.not.allocated(xanew)) then
        allocate(xanew(3,natoms))
        call memory('A','D',3*natoms,'nose')
      endif
      if (.not.allocated(xaold)) then
        allocate(xaold(3,natoms))
        call memory('A','D',3*natoms,'nose')
        do ia = 1,natoms
          do i = 1,3
            xaold(i,ia)=0.0d0
          enddo
        enddo
      endif

C Define constants and conversion factors .....................................
      dt2   = dt**2
      dtby2 = dt/2.0d0
      twodt = dt*2.0d0
      tol   = 1.0d-12

      if (iunit .eq. 1) then
C  convert target ionic temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in eV if Temp is in Kelvin)
         tekin = 0.5d0 * (3.d0 * natoms - ct) * 8.617d-5 * tt
C  convert F/m in (eV/Angstrom)/amu  to  Angstrom/fs**2
        fovermp = 0.009579038
      else
C  convert target temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in Ry if Temp is in Kelvin)
        tekin = eV * 0.5d0 * (3.d0 * natoms - ct) * 8.617d-5 * tt
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038 * Ang**2 / eV
      endif
C Initialize variables if current time step is the first of the simulation
      if (istep .eq. 1) then

         if (.not. xv_file_read) then

C     Compute old positions in terms of current positions and velocities
C     if the time step is the first of the simulation 
!     and we start from x(t), v(t) *at the same time*.
!     (e.g., when the velocities are constructed from
!      the Boltzmann distribution).
!     In this case the algorithm works out well.
!     Nose variables are set to zero, as there is currently no
!     better way to initialize them...

            x = 0.0d0
            xold = 0.0d0
            do ia = 1,natoms
               do i = 1,3
                  xaold(i,ia) = xa(i,ia) - dt * va(i,ia)
     .                 + (dt2/2.0d0) * fovermp * fa(i,ia) / ma(ia)
               enddo
            enddo

         else

!         For restarts, we need information about the old 
!         positions, and the Nose variables
!
           if (Node .eq. 0) then
            call io_assign(iacc)
            open(unit=iacc,file="NOSE_RESTART", form="formatted",
     $           status="old", action="read", position="rewind")
            read(iacc,*) old_natoms, old_dt
            read(iacc,*) x, xold
            if (old_natoms .ne. natoms) then
               write(6,"(a)") "Wrong number of atoms in NOSE_RESTART"
               stop_flag = .true.
            else
               do ia = 1, natoms
                  read(iacc,*) dummy_iza, (xaold(i,ia),i=1,3) ! old positions
                  if (dummy_iza .ne. iza(ia)) then
                     write(6,"(a)")
     $                     "Wrong species number in NOSE_RESTART"
                     stop_flag = .true.
                     exit  ! loop
                  endif
               enddo
            endif
            call io_close(iacc)
            if (.not. stop_flag) then
               write(6,*)
     $         "MD restart: Read old positions and Nose variables",
     $              " from NOSE_RESTART"
               if (abs(old_dt - dt) .gt. 1.0d-8) then
                  write(6,*) "**WARNING: Timestep has changed. Old: ",
     $                 old_dt, " New: ", dt
                  write(6,*) "**WARNING: Approximating old positions."
                  ! First order, using the positions and velocities 
                  ! at t-old_dt (positions from NOSE_RESTART, velocities
                  !              from XV file)
                  xaold(1:3,1:natoms) = xaold(1:3,1:natoms) -
     $                              (dt-old_dt) * va(1:3,1:natoms)
               endif  ! dt /= old_dt
            endif     ! still processing

            endif     ! IONode
            
            call broadcast(stop_flag)
            if (stop_flag) then
               stop_flag = .false.
               call die()  ! Proper way to stop MPI...
            endif

            call broadcast(x)
            call broadcast(xold)
            call broadcast(xaold(1:3,1:natoms))

         endif     ! xv_file_read
      endif        ! istep == 1
C ..................

C Compute uncorrected next positions .....................................
      do ia = 1,natoms
        do i = 1,3
          xanew(i,ia) =  2.0d0 * xa(i,ia) - xaold(i,ia) +
     .                   dt2 * fovermp * fa(i,ia) / ma(ia)
        enddo
      enddo
C ...................

C Compute uncorrected velocities and kinetic energy ......................
      kin = 0.d0
      do ia = 1,natoms
        do i = 1,3
          va(i,ia) = (xanew(i,ia) - xaold(i,ia)) / twodt
          kin = kin + 0.5d0 * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo
C ..................

C Compute initial guess for Nose variables at next time step .............
      xnew = 2.0d0 * x - xold
C ...................

C Start selfconsistency loop to calculate Nose variable ..................
10    continue

      xlast = xnew
        
C xdot and hdot (time derivatives at current time), and related stuff
      xdot = (xnew - xold) / twodt
      fact = (1.0/(1.0+xdot*dtby2))

C  Compute Nose variable for next iteration
      xnew = 2.0d0 * x - xold 
     .       + (dt2/mn) * 2.0 * (fact**2 * kin - tekin)

C Check if selfconsistency has been reached
      diff = abs(xnew - xlast)
      if (xlast .eq. 0.0d0) then
        if (diff .gt. tol)  goto 10
      else
        if (diff/abs(xlast) .gt. tol)  goto 10
      endif
C ...................

C Calculate corrected atomic coordinates at next time step, 
C and corrected velocities and kinetic energy at current time step .........
      do ia = 1,natoms
        do i = 1,3
          xanew(i,ia) = fact * ( xanew (i,ia) +
     .                   dtby2 * xdot * xaold(i,ia))
          va(i,ia) = fact * va(i,ia)
        enddo
      enddo
      kin = kin * fact**2 
C ...................

C Save current atomic positions as old ones, 
C   and next positions as current ones

      xaold(1:3,1:natoms) = xa(1:3,1:natoms)
      xa(1:3,1:natoms) = xanew(1:3,1:natoms)

      xold = x
      x = xnew

C Calculate Kinetic and potential energies ................................
C Kinetic energy of Nose variable
      kn = (1.0d0 / 2.0d0) * mn * xdot**2

C Potential energy of Nose variable (in eV)
      vn = 2.0d0 * tekin * xold

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        temp = kin / (0.5d0 * (3.d0 * natoms - ct) * 8.617d-5)
      else
        temp = kin / (0.5d0 * (3.d0 * natoms - ct) * 8.617d-5 * eV)
      endif

      if (Node .eq. 0) then
!
!       Save (now old) positions and nose variables to NOSE_RESTART
!
         call io_assign(iacc)
         open(unit=iacc,file="NOSE_RESTART", form="formatted",
     $        status="unknown", action= "write", position="rewind")
         write(iacc,*) natoms, dt
         write(iacc,*) x, xold
         do ia = 1, natoms
            write(iacc,*) iza(ia), (xaold(i,ia),i=1,3) ! forces
         enddo
         call io_close(iacc)
      endif

C .....................

      end subroutine nose

      subroutine anneal(istep,iunit,ianneal,taurelax,bulkm,
     .               natoms,fa,stress,tp,tt,dt,
     .               ma,ntcon,va,xa,h,kin,
     .               temp,pressin)
C *************************************************************************
C Subroutine for MD simulations with a TARGET PRESSURE AND TEMPERATURE.
C The system is driven to a desired temperature and pressure in
C a given time, by rescaling the velocities and the cell shape and size.
C It needs an estimate of the bulk modulus of the system, to determine
C the rate of change of cell shape to accomodate to the target
C pressure in the required time. A wrong estimate will simply
C drive the system to the desired pressure, but in a different time
C than the especified in the input. (See Kittel for representative
C values of bulk moduli for materials).
C
C Written by P.Ordejon, November'96
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Time in femtoseconds
C Atomic masses in atomic mass units
C
C Other units depend on input option:
C
C   Option iunit = 1:
C     Pressures are in eV/A**3 ( = 160.20506 GPa)
C     Energies are in eV
C     Distances are in Angstrom
C     Bulk modulus in eV/A**3
C   Option iunit = 2:
C     Pressures are in Ry/Bohr**3 
C     Energies are in Ry
C     Distances are in Bohr
C     Bulk modulus in Ry/Bohr**3
C ************************* INPUT *********************************************
C integer istep       : Number of time step during simulation
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer ianneal     : Work mode option:
C                       1 = reach target temperature only
C                       2 = reach target pressure only
C                       3 = reach target temperature and pressure
C real*8 taurelax     : Relaxation time to reach desired T and P
C real*8 bulkm        : Estimate of the Bulk Modulus of the system
C integer natoms      : Number of atoms in the simulation cell
C real*8 fa(3,natoms) : Atomic forces 
C real*8 stress(3,3)  : Stress tensor components 
C real*8 tp           : Target pressure 
C real*8 tt           : Target temperature 
C real*8 dt           : Length of the time step 
C real*8 ma(natoms)   : Atomic masses 
C integer ntcon         : Total number of position constraints imposed
C real*8 va(3,natoms) : Atomic velocities
C                       (used only if istep = 1)
C ******************* INPUT AND OUTPUT ****************************************
C real*8 xa(3,natoms) : Atomic coordinates 
C                      (input: current time step; output: next time step)
C real*8 h(3,3)       : Matrix of the vectors defining the unit
C                       cell at the current time step 
C                       h(i,j) is the ith component of jth basis vector
C                      (input: current time step; output: next time step)
C ************************* OUTPUT *******************************************
C real*8 kin            : Kinetic energy of the atomic system 
C real*8 temp           : Instantaneous system temperature 
C real*8 pressin        : Instantaneous system pressure 
C *****************************************************************************


      integer 
     .   natoms,ntcon,istep,ianneal,iunit

      real(dp)
     .  bulkm,dt,fa(3,natoms),h(3,3),kin,
     .  ma(natoms),stress(3,3),taurelax,tp,tt,
     .  va(3,natoms),xa(3,natoms)

C Internal variables .............................................................

      integer
     .  ct,i,ia,info,j,k

      real(dp)
     .  dt2,fovermp,hi(3,3),hs(3),
     .  pgas,press(3,3),pressin,rfac,rfac2,
     .  tekin,temp,twodt,vol,volcel

      real(dp), dimension(:,:), allocatable, save ::
     .  s,sdot,snew,sold,sunc

       external
     .  volcel, memory
C ....................................................................

      if (iunit .ne. 1 .and. iunit .ne. 2) then
        if (Node.eq.0) then
          write(6,*) 'anneal: Wrong iunit option;  must be 1 or 2'
        endif
        stop
      endif

!!!      if (taurelax/dt .lt. 0.1) return

      ct = 3 + ntcon
      if (natoms .eq. 1) ct = 0

C Allocate local memory
      allocate(s(3,natoms))
      call memory('A','D',3*natoms,'anneal')
      allocate(sdot(3,natoms))
      call memory('A','D',3*natoms,'anneal')
      allocate(snew(3,natoms))
      call memory('A','D',3*natoms,'anneal')
      allocate(sunc(3,natoms))
      call memory('A','D',3*natoms,'anneal')
      if (.not.allocated(sold)) then
        allocate(sold(3,natoms))
        call memory('A','D',3*natoms,'anneal')
      endif

C Define constants and conversion factors .......................................
      dt2   = dt**2
      twodt = dt*2.0d0

      if (iunit .eq. 1) then
C  convert target ionic temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in eV if Temp is in Kelvin)
        tekin = 0.5d0 * (3.d0 * natoms - ct) * 8.617d-5 * tt
C  convert F/m in (eV/Angstrom)/amu  to  Angstrom/fs**2
        fovermp = 0.009579038
      else
C  convert target temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in Ry if Temp is in Kelvin)
        tekin = eV * 0.5d0 * (3.d0 * natoms - ct) * 8.617d-5 * tt
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038 * Ang**2 / eV
      endif

C  calculate cell volume at current time
      vol = volcel( h )
C ........................

C Compute Parrinello-Rahman variables (H and scaled coordinates) ............
C Compute Inverse of H and G at current time 
      call inver(h,hi,3,3,info)
      if (info .ne. 0) stop 'anneal: INVER failed'

C Calculate scaled coordinates (referred to matrix H) at current time
      do ia = 1,natoms
        do i = 1,3
          s(i,ia) = 0.0d0
          do j = 1,3
            s(i,ia) = s(i,ia) + hi(i,j) * xa(j,ia)
          enddo
        enddo
      enddo

C Initialize variables if current time step is the first of the simulation
      if (istep .eq. 1) then
        do ia = 1,natoms
          do i = 1,3
            sold(i,ia) = 0.0d0
            do j = 1,3
              sold(i,ia) = sold(i,ia) + hi(i,j)*(xa(j,ia)-dt*va(j,ia)
     .                     + (dt2/2.0d0) * fovermp * fa(j,ia) / ma(ia))
            enddo
          enddo
        enddo
      endif
C ..................

C Compute uncorrected next positions .....................................
      do ia = 1,natoms
        do i = 1,3
          sunc(i,ia) = -sold(i,ia) + 2.0d0 * s(i,ia)
          do k = 1,3
            sunc(i,ia) = sunc(i,ia) + 
     .                   dt2 * hi(i,k) * fovermp * fa(k,ia) / ma(ia)
          enddo
        enddo
      enddo
C ...................

C Calculate uncorrected velocities at current time
      do ia = 1,natoms
        do i = 1,3
          sdot(i,ia) = (sunc(i,ia) - sold(i,ia)) / twodt
        enddo
      enddo

C Calculate pressure tensor at current time and ideal gas pressure
      do i = 1,3
        do j = 1,3
          press(i,j) = 0.0d0
        enddo
      enddo
      do ia = 1,natoms
        do i = 1,3
          hs(i) = 0.0d0
          do j = 1,3
            hs(i) = hs(i) + h(i,j) * sdot(j,ia)
          enddo
        enddo
        do j = 1,3
          do i = 1,3
            press(i,j) = press(i,j) + ma(ia) * hs(i) * hs(j) / fovermp
          enddo
        enddo
      enddo
      pgas = 0.0d0
      do i = 1,3
        pgas = pgas + press(i,i) / vol
      enddo
      pgas = pgas / 3.0d0
      do j = 1,3
        do i = 1,3
          press(i,j) = press(i,j) / vol - stress(i,j)
        enddo
      enddo

C Compute internal pressure  (pressin = 1/3 Tr (press))   at current time
      pressin = 0.0d0
      do i = 1,3
        pressin = pressin + press(i,i)
      enddo
      pressin = pressin / 3.0d0

C Compute kinetic energy
      kin = (3.0d0 / 2.0d0) * pgas * vol

      write(*,*) "Anneal: Kinetic Energy= ", kin

      if (ianneal .eq. 1 .or. ianneal .eq. 3) then
C Correct velocities to reach target termperature
      if (kin .eq. 0.0) then
        rfac2 = 1.0d0 + dt/taurelax
      else
         
        rfac2 = (1.0d0 + dt/taurelax * (tekin/kin -1.0d0))
      endif
      if (rfac2 .le. 0.0) call die('Wrong anneal parameter')
      rfac = sqrt(rfac2)
      write(*,*) "Anneal: Velocity scale factor = ", rfac

      sdot(1:3,1:natoms) = rfac * sdot(1:3,1:natoms)

C Compute again pressure, with corrected velocities

      press(1:3,1:3) = 0.0_dp

      do ia = 1,natoms
        do i = 1,3
          hs(i) = 0.0d0
          do j = 1,3
            hs(i) = hs(i) + h(i,j) * sdot(j,ia)
          enddo
        enddo
        do j = 1,3
          do i = 1,3
            press(i,j) = press(i,j) + ma(ia) * hs(i) * hs(j) / fovermp
          enddo
        enddo
      enddo
      pgas = 0.0d0
      do i = 1,3
        pgas = pgas + press(i,i) / vol
      enddo
      pgas = pgas / 3.0d0

      press(1:3,1:3) = press(1:3,1:3)/vol - stress(1:3,1:3)

C Compute internal pressure  (pressin = 1/3 Tr (press))   at current time
      pressin = 0.0
      do i = 1,3
        pressin = pressin + press(i,i)
      enddo
      pressin = pressin / 3.0d0
      endif


C Correct new possitions according to corrected velocities
      do ia = 1,natoms
        do i = 1,3
          snew(i,ia) = sold(i,ia) + twodt * sdot(i,ia)
        enddo
      enddo

      if (ianneal .eq. 2 .or. ianneal .eq. 3) then

C Correct cell shape to reach target pressure
      do i = 1,3
        do j = 1,3
          if (i .ne. j) then
            rfac2 = 1.0 + (dt / taurelax) * press(i,j) 
     .                    / (0.5 * bulkm)
          else
            rfac2 = 1.0 + (dt / taurelax) * (press(i,i) - tp) 
     .                    / (0.5 * bulkm)
          endif
          if (rfac2 .le. 0.0) then
            write(6,*) 'Wrong anneal parameter'
            stop
          endif
          rfac = sqrt(rfac2)
          h(i,j) = rfac * h(i,j)
        enddo
      enddo
      write(*,*) "Anneal: Cell scale factor = ", rfac
      endif

C Save current atomic positions as old ones, 
C   and next positions as current ones
      do i = 1,3
        do ia = 1,natoms
          sold(i,ia) = s(i,ia)
          s(i,ia) = snew(i,ia)
        enddo
      enddo

C Transform back to absolute coordinates 
      do ia = 1,natoms
        do i = 1,3
          xa(i,ia) = 0.0d0
          va(i,ia) = 0.0d0
          do j = 1,3
            xa(i,ia) = xa(i,ia) + h(i,j) * s(j,ia)
            va(i,ia) = va(i,ia) + h(i,j) * sdot(j,ia)
          enddo
        enddo
      enddo
C ....................

C Calculate Kinetic and potential energies ................................
C Kinetic energy of atoms 
      kin = (3.0d0 / 2.0d0) * pgas * vol

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        temp = 3.0d0*vol*pgas/(3.0d0*natoms-ct)/8.617d-5
      else
        temp = 3.0d0*vol*pgas/(3.0d0*natoms-ct)/8.617d-5/eV
      endif

C .....................

C Deallocate local memory
      call memory('D','D',size(s),'anneal')
      deallocate(s)
      call memory('D','D',size(sdot),'anneal')
      deallocate(sdot)
      call memory('D','D',size(snew),'anneal')
      deallocate(snew)
      call memory('D','D',size(sunc),'anneal')
      deallocate(sunc)

!!!      taurelax = taurelax - dt
      return
      end subroutine anneal


      subroutine verlet1(istep,iunit,iquench,natoms,fa,dt,ma,ntcon,va,
     .                   xa,kin,temp)
C *************************************************************************
C Subroutine for MD simulations using the Original Verlet Algrithm.
C (See Allen-Tildesley, Computer Simulations of Liquids, pg. 78)
C
C Written by P.Ordejon, November'96
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Time in femtoseconds
C Atomic masses in atomic mass units
C
C Other units depend on input option:
C
C   Option iunit = 1:
C     Energies are in eV
C     Distances are in Angstrom
C   Option iunit = 2:
C     Energies are in Ry
C     Distances are in Bohr
C ************************* INPUT *********************************************
C integer istep         : Number of time step during the simulation
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer iquench       : Option for quenching:   
C                              0 = no quenching (standard dynamics)
C                              1 = power quenching (set to cero velocity 
C                                 components opposite to force)
C integer natoms        : Number of atoms in the simulation cell
C real*8 fa(3,natoms)   : Atomic forces 
C real*8 dt             : Length of the time step
C real*8 ma(natoms)     : Atomic masses 
C integer ntcon         : Total number of position constraints imposed
C real*8 va(3,natoms)   : Atomic velocities 
C                         (used only if istep = 1)
C ******************* INPUT AND OUTPUT ****************************************
C real*8 xa(3,natoms)   : Atomic coordinates 
C                        (input: current time step; output: next time step)
C ************************* OUTOPUT *******************************************
C real*8 kin            : Kinetic energy at current time step 
C real*8 temp           : Instantaneous system temperature 
C *****************************************************************************
C
      integer 
     .   natoms,ntcon,istep,iquench,iunit

      real(dp)
     .  dt,fa(3,natoms),kin,ma(natoms),
     .  va(3,natoms),xa(3,natoms)

      external
     .  memory

C Internal variables ..........................................................
 
      integer
     .  ct,i,ia

      real(dp)
     .  dot,dt2,fovermp,temp,twodt

      real(dp), dimension(:,:), allocatable, save ::
     .  xanew
      real(dp), dimension(:,:), allocatable, save ::
     .  xaold

C ........................

      if (iunit .ne. 1 .and. iunit .ne. 2) then
        if (Node.eq.0) then
          write(6,*) 'verlet1: Wrong iunit option;  must be 1 or 2'
        endif
        stop
      endif
      ct = 3 + ntcon
      if (natoms .eq. 1) ct = 0

C Allocate local memory
      allocate(xanew(3,natoms))
      call memory('A','D',3*natoms,'verlet1')
      if (.not.allocated(xaold)) then
        allocate(xaold(3,natoms))
        call memory('A','D',3*natoms,'verlet1')
      endif

C Define constants and conversion factors .....................................
      dt2   = dt**2
      twodt = dt*2.0d0

      if (iunit .eq. 1) then
C  convert F/m in (eV/Amstrong)/amu  to  Amstrong/fs**2
        fovermp = 0.009579038
      else
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038  * Ang**2 / eV
      endif

C ........................

C Compute old coordinates if the time step is the first of the simulation .....
      if (istep .eq. 1) then
        do ia = 1,natoms
          do i = 1,3
            xaold(i,ia) = xa(i,ia) - dt * va(i,ia)
     .                + (dt2/2.0d0) * fovermp * fa(i,ia) / ma(ia)
          enddo
        enddo
      endif

C Compute positions at next time step.....................................
      do ia = 1,natoms
        do i = 1,3
          xanew(i,ia) = - xaold(i,ia) + 2.0d0 * xa(i,ia)
     .                  + dt2 * fovermp * fa(i,ia) / ma(ia)
        enddo
      enddo
C ...................

C Calculate velocities at current time .....................................
      do ia = 1,natoms
        do i = 1,3
          va(i,ia) = (xanew(i,ia) - xaold(i,ia)) / twodt
        enddo
      enddo
C ...................

C Quench option if iquench = 0 ..............................................
      if (iquench .eq. 1) then

C Quench velocity components going uphill
        do ia = 1,natoms
          do i = 1,3
            dot = va(i,ia) * fa(i,ia)
            if (dot .lt. 0.0) then
              va(i,ia) = 0.0
              xanew(i,ia) = xa(i,ia)
            endif
          enddo
        enddo

      endif
C......................


C Save current atomic positions as old ones, 
C   and next positions as current ones .....................................
      do i = 1,3
        do ia = 1,natoms
          xaold(i,ia) = xa(i,ia)
          xa(i,ia) = xanew(i,ia)
        enddo
      enddo
C ....................

C Calculate kinetic energy and temperature at current time ...................
C Kinetic energy of atoms
      kin = 0.0d0
      do ia = 1,natoms
        do i = 1,3
          kin = kin + 0.5d0 * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        temp = 2.0d0*kin/(3.0d0*natoms-ct)/8.617d-5
      else
        temp = 2.0d0*kin/(3.0d0*natoms-ct)/8.617d-5/eV
      endif

C .....................

C Deallocate local memory
      call memory('D','D',size(xanew),'verlet1')
      deallocate(xanew)

      end subroutine verlet1
    
      subroutine verlet2(istep,iunit,iquench,natoms,fa,dt,ma,ntcon,va,
     .                   xa,kin,temp)
C *************************************************************************
C Subroutine for MD simulations using the velocity-Verlet Algrithm.
C (See Allen-Tildesley, Computer Simulations of Liquids, pg. 81)
C
C Written by P.Ordejon, November'96
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Time in femtoseconds
C Atomic masses in atomic mass units
C
C Other units depend on input option:
C
C   Option iunit = 1:
C     Energies are in eV
C     Distances are in Angstrom
C   Option iunit = 2:
C     Energies are in Ry
C     Distances are in Bohr
C ************************* INPUT *********************************************
C integer istep         : Number of time step during the simulation
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer iquench       : Option for quenching:
C                              0 = no quenching (standard dynamics)
C                              1 = power quenching (set to cero velocity
C                                 components opposite to force)
C integer natoms        : Number of atoms in the simulation cell
C real*8 fa(3,natoms)   : Atomic forces 
C real*8 dt             : Length of the time step
C real*8 ma(natoms)     : Atomic masses 
C integer ntcon         : Total number of position constraints imposed
C real*8 va(3,natoms)   : Atomic velocities
C                         (used only if istep = 1)
C ******************* INPUT AND OUTPUT ****************************************
C real*8 xa(3,natoms)   : Atomic coordinates 
C                        (input: current time step; output: next time step)
C ************************* OUTOPUT *******************************************
C real*8 kin            : Kinetic energy at current time step 
C real*8 temp           : Instantaneous system temperature
C *****************************************************************************

      integer 
     .   natoms,ntcon,istep,iquench,iunit

      real(dp)
     .  dt,fa(3,natoms),kin,ma(natoms),
     .  va(3,natoms),xa(3,natoms)

      external
     .  memory

C Internal variables ..........................................................
 
      integer
     .  ct,i,ia

      integer :: old_natoms, iacc, dummy_iza
      real(dp) :: old_dt

      real(dp)
     .  dot,dt2,dtby2,fovermp,temp

      real(dp), dimension(:,:), allocatable, save ::
     .  accold,vold
C ........................

      if (iunit .ne. 1 .and. iunit .ne. 2) then
        if (Node.eq.0) then
          write(6,*) 'verlet2: Wrong iunit option;  must be 1 or 2'
        endif
        stop
      endif
      ct = 3 + ntcon
      if (natoms .eq. 1) ct = 0

C Allocate local memory - only done once as data must be saved. As a
C result the memory is not deallocated at the end of the routine.
      if (.not.allocated(accold)) then
        allocate(accold(3,natoms))
        call memory('A','D',3*natoms,'verlet2')
      endif
      if (.not.allocated(vold)) then
        allocate(vold(3,natoms))
        call memory('A','D',3*natoms,'verlet2')
      endif

C Define constants and conversion factors .....................................
      dt2   = dt**2
      dtby2 = dt/2.0d0

      if (iunit .eq. 1) then
C  convert F/m in (eV/Amstrong)/amu  to  Amstrong/fs**2
        fovermp = 0.009579038
      else
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038 * Ang**2 / eV
      endif
C ........................

      
      if (istep .eq. 1) then

         if (.not. xv_file_read) then

C     Compute old accelerations and velocities 
C     if the time step is the first of the simulation ...........................
!     and we start from x(t), v(t) *at the same time*.
!     (e.g., when the velocities are constructed from
!      the Boltzmann distribution).
!     In this case the algorithm works out well.

            do ia = 1,natoms
               do i = 1,3
                  accold(i,ia) = fovermp * fa(i,ia) / ma(ia)
                  vold(i,ia) = va(i,ia) - dt * accold(i,ia)
               enddo
            enddo

        else

!         For restarts, we need information about the old 
!         forces, in order to match the velocities
!         correctly (the velocities in the XV file are
!         one time step behind, so they are already the
!         "old" velocities).
!
          if (Node .eq. 0) then
           call io_assign(iacc)
           open(unit=iacc,file="VERLET_FORCES", form="formatted",
     $          status="old", action="read", position="rewind")
           read(iacc,*) old_natoms, old_dt
           if (old_natoms .ne. natoms) then
               write(6,"(a)") "Wrong number of atoms in VERLET_FORCES"
               stop_flag = .true.
           else
              do ia = 1, natoms
                 read(iacc,*) dummy_iza, (accold(i,ia),i=1,3) ! forces
                 if (dummy_iza .ne. iza(ia)) then
                    write(6,"(a)")
     $                   "Wrong species number in VERLET_FORCES"
                    stop_flag = .true.
                    exit        ! loop
                 endif
                 accold(:,ia) = fovermp * accold(:,ia) / ma(ia)
                 vold(:,ia)  = va(:,ia)
              enddo
           endif
           call io_close(iacc)
           if (.not. stop_flag) then
            write(6,*) "MD restart: Read old forces from VERLET_FORCES"
            if (abs(old_dt - dt) .gt. 1.0d-8) then
               write(6,*) "Timestep has changed. Old: ", old_dt,
     $                     " New: ", dt
            endif
           endif ! still processing

          endif             ! IONode

          call broadcast(stop_flag)
          if (stop_flag) then
             stop_flag = .false.
             call die()         ! Proper way to stop MPI...
          endif

          call broadcast(accold(1:3,1:natoms))
          call broadcast(vold(1:3,1:natoms))

       endif      ! XV file read

      endif    ! first step

C ....................
C Compute velocities at current time step, 
! using the previous step's velocities and the previous and current forces.

      if ((istep .eq. 1) .and. xv_file_read) then
!
!        Use old time step in case it is different, only in 
!        the first step.
!
         do ia = 1,natoms
            do i = 1,3
               va(i,ia) = vold(i,ia) + 0.5d0 * old_dt
     .              * (accold(i,ia) + fovermp * fa(i,ia) / ma(ia))
            enddo
         enddo
         
      else
!
!        Current timestep.
!
         do ia = 1,natoms
            do i = 1,3
               va(i,ia) = vold(i,ia) + dtby2 
     .              * (accold(i,ia) + fovermp * fa(i,ia) / ma(ia))
            enddo
         enddo
      endif    

C Quench option if iquench = 0 ..............................................
      if (iquench .eq. 1) then

C Quench velocity components going uphill
         do ia = 1,natoms
           do i = 1,3
             dot = va(i,ia) * fa(i,ia)
             if (dot .lt. 0.0) va(i,ia) = 0.0
           enddo
         enddo

      endif
C ................

C Compute positions at next time step.....................................
      do ia = 1,natoms
        do i = 1,3
          xa(i,ia) = xa(i,ia) + dt * va(i,ia) 
     .                  + dt2 / 2.0d0 * fovermp * fa(i,ia) / ma(ia)
        enddo
      enddo
C ...................

C Save current velocities and accelerations as old ones .....................
      do i = 1,3
        do ia = 1,natoms
          vold(i,ia) = va(i,ia)
          accold(i,ia) = fovermp * fa(i,ia) / ma(ia)
        enddo
      enddo
C ....................

C Calculate kinetic energy and temperature at current time ...................
C Kinetic energy of atoms 
      kin = 0.0d0
      do ia = 1,natoms
        do i = 1,3
          kin = kin + 0.5d0 * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        temp = 2.0d0*kin/(3.0d0*natoms-ct)/8.617d-5
      else
        temp = 2.0d0*kin/(3.0d0*natoms-ct)/8.617d-5/eV
      endif

!
!       Save (now old) forces to VERLET_FORCES
!
      if (Node .eq. 0) then
         call io_assign(iacc)
         open(unit=iacc,file="VERLET_FORCES", form="formatted",
     $        status="unknown", action= "write", position="rewind")
         write(iacc,*) natoms, dt
         do ia = 1, natoms
            write(iacc,*) iza(ia), (fa(i,ia),i=1,3) ! forces
         enddo
         call io_close(iacc)
      endif
C .....................

      end subroutine verlet2

      end module m_dynamics

    
