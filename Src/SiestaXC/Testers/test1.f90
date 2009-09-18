PROGRAM siestaXCtest1

  ! Checks the consistency of the energies and its derivativas, returned by  
  ! the GGA (and LDA) XC routines contained in m_ggaxc and m_ldaxc modules.
  ! Two things are checked for each functional:
  !   That the xc potentials Vx=dEx/dDens and Vc=dEc/dDens are always negative.
  !   That the (analytical) returned values dE/dDens and dE/dGrad(Dens) coincide
  !   with their numerical counterparts, calculated as finite differences.
  ! J.M.Soler. Sept.2009

  ! Used module procedures and parameters
  USE siestaXC, only: dp => siestaXC_std_p
  USE siestaXC, only: ggaxc
  USE siestaXC, only: ldaxc

  implicit none

  ! Tester parameters
  integer, parameter:: irel    =  1      ! Relativistic? 0=>no, 1=>yes
  integer, parameter:: nCalls  = 20      ! Number of calls to each routine
  integer, parameter:: nfTot   =  9      ! Total number of functionals
  integer, parameter:: nSpin   =  2      ! Number of spin components
  real(dp),parameter:: densMax = 0.1_dp  ! Upper limit to density value
  real(dp),parameter:: gradMax = 0.5_dp  ! Upper limit to density gradient
  real(dp),parameter:: deltaDens  = 1.e-6_dp  ! Finite diff. change
  real(dp),parameter:: deltaGrad  = 1.e-6_dp  ! Finite diff. change
  real(dp),parameter:: dEdDensTol = 1.e-6_dp  ! Tolerance for dE/dDens
  real(dp),parameter:: dEdGradTol = 1.e-6_dp  ! Tolerance for dE/dGadDens

  ! Functionals to be tested
  !                  1,       2,       3,       4,       5,   
  !                  6,       7,       8,       9   
  character(len=3):: &
    func(nfTot) = (/'LDA',   'LDA',   'GGA',   'GGA',   'GGA',    &
                    'GGA',   'GGA',   'GGA',   'GGA'             /)
  character(len=6):: &
    auth(nfTot) = (/'PZ    ','PW92  ','PW91  ','PBE   ','RPBE  ', &
                    'revPBE','LYP   ','WC    ','PBESOL'          /) 

  ! Tester variables and arrays
  integer :: errorC, errorX, iCall, iDelta, iDeriv, iSpin, ix, &
             jf, mSpin, one, two, three, four
  real(dp):: epsC, epsX, dens(nSpin), &
             dEcdDens(nSpin), dEcdGrad(3,nSpin), &
             dExdDens(nSpin), dExdGrad(3,nSpin), densPol, densTot, &
             dVxdDens(nSpin,nSpin), dVcdDens(nSpin,nSpin), &
             Ec, Ex, gradDens(3,nSpin), ran(4)
  real(dp):: epsC0, epsX0, dens0(nSpin), gradDens0(3,nSpin)
  real(dp):: dEcdDens0(nSpin), dEcdGrad0(3,nSpin), &
             dExdDens0(nSpin), dExdGrad0(3,nSpin)
  real(dp):: dEcdDensN(nSpin), dEcdGradN(3,nSpin), &
             dExdDensN(nSpin), dExdGradN(3,nSpin)
  logical :: errorsFound = .false.

  ! Print header
  print'(a4,a7,a3,2a6,3a15,a4)', &
    'func', 'auth  ', 'xyz', 'iSpin', 'Ex|Ec', &
    'deriv(anal)', 'deriv(num)', 'diff  ', 'err'

  ! Loop on calls with different densities
  do iCall = 1,nCalls

    ! Choose random density and gradient
    trial: do
      do iSpin = 1,nSpin
        call random_number( ran(1:4) )
        dens0(iSpin) = densMax * ran(4)
        dens0(iSpin) = dens0(iSpin) + 2*deltaDens  ! Avoid dens0-deltaDens<0
        gradDens0(1:3,iSpin) = gradMax * (2*ran(1:3)-1)
      end do
      if (nSpin==4) then           ! Non-collinear spin
        one   = 1   ! A silly thing to satisfy the compiler when nSpin<4
        two   = 2
        three = 3
        four  = 4
        densTot = dens0(one) + dens0(two)
        densPol = sqrt( (dens0(one)-dens0(two))**2 &
                      + 4*(dens0(three)**2+dens0(four)**2) )
        if (densTot<densPol) then  ! Invalid density matrix => try again
          cycle trial
        else                       ! Valid density matrix => accept it
          exit trial
        end if
      else                         ! Collinear spin => no check needed
        exit trial
      end if
    end do trial

    ! Make density constant, to check that GGAs reproduce LDA
!    gradDens0(:,:) = 0

    ! Make spin collinear to check that nSpin=4 reproduces nSpin=2
!    dens0(3:4) = 0
!    gradDens0(:,3:4) = 0

    ! Loop on functionals
    do jf = 1,nfTot

      ! Find energy density and its derivatives
      if (func(jf)=='LDA') then
        call ldaxc( auth(jf), irel, nSpin, dens0, &
                    epsX0, epsC0, dExdDens0, dEcdDens0, dVxdDens, dVcdDens )
      else if (func(jf)=='GGA') then
        call ggaxc( auth(jf), irel, nSpin, dens0, gradDens0, &
                    epsX0, epsC0, dExdDens0, dEcdDens0, dExdGrad0, dEcdGrad0 )
      else
        stop 'ERROR: unknown functional'
      end if

      ! Loop on numerical derivative and finite-difference delta
      do iDeriv = 0,4*nSpin-1

        ! Select magnitude to change
        iSpin = iDeriv/4 + 1 ! Spin index: iSpin=1,...,nSpin
        ix = mod(iDeriv,4)   ! Cartes. comp. of grad(dens) (0 for dens itself)

        ! For LDA, do not perform derivative with respect to gradient
        if (ix>0 .and. func(jf)=='LDA') cycle

        ! Perform numerical derivative
        dExdDensN = 0
        dEcdDensN = 0
        dExdGradN = 0
        dEcdGradN = 0
        do iDelta = -1,1,2

          ! Change select magnitude
          dens = dens0
          gradDens = gradDens0
          if (ix==0) then
            dens(iSpin) = dens0(iSpin) + iDelta * deltaDens
          else
            gradDens(ix,iSpin) = gradDens0(ix,iSpin) + iDelta * deltaGrad
          end if

          ! Find energy density at modified density or gradient
          if (func(jf)=='LDA') then
            call ldaxc( auth(jf), irel, nSpin, dens,   &
                        epsX, epsC, dExdDens, dEcdDens, dVxdDens, dVcdDens )
          else
            call ggaxc( auth(jf), irel, nSpin, dens, gradDens, &
                        epsX, epsC, dExdDens, dEcdDens, dExdGrad, dEcdGrad )
          end if

          ! Add to finite-difference derivative
          mSpin = min(nSpin,2)
          densTot = sum(dens(1:mSpin))
          Ex = densTot * epsX           ! Exchange energy per electron
          Ec = densTot * epsC           ! Correlation energy per electron
          if (ix==0) then
            dExdDensN(iSpin) = dExdDensN(iSpin) + iDelta*Ex/(2*deltaDens)
            dEcdDensN(iSpin) = dEcdDensN(iSpin) + iDelta*Ec/(2*deltaDens)
          else
            dExdGradN(ix,iSpin) = dExdGradN(ix,iSpin) + iDelta*Ex/(2*deltaGrad)
            dEcdGradN(ix,iSpin) = dEcdGradN(ix,iSpin) + iDelta*Ec/(2*deltaGrad)
          end if

        end do ! iDelta

        ! Divide by 2 the non-diagonal components for non-collinear spin (why?)
        dExdDensN(mSpin+1:nSpin) = dExdDensN(mSpin+1:nSpin) / 2
        dEcdDensN(mSpin+1:nSpin) = dEcdDensN(mSpin+1:nSpin) / 2
        dExdGradN(:,mSpin+1:nSpin) = dExdGradN(:,mSpin+1:nSpin) / 2
        dEcdGradN(:,mSpin+1:nSpin) = dEcdGradN(:,mSpin+1:nSpin) / 2

        ! Check for errors (positive Vxc or inconsistent derivatives)
        errorX = 0
        errorC = 0
        if (ix==0) then
          if (abs(dExdDens0(iSpin)-dExdDensN(iSpin)) > dEdDensTol) errorX=2
          if (abs(dEcdDens0(iSpin)-dEcdDensN(iSpin)) > dEdDensTol) errorC=2
          if (iSpin<=2 .and. dExdDens0(iSpin)>0._dp) errorX = 1
          if (iSpin<=2 .and. dEcdDens0(iSpin)>0._dp) errorC = 1
        else
          if (abs(dExdGrad0(ix,iSpin)-dExdGradN(ix,iSpin))>dEdGradTol) errorX=3
          if (abs(dEcdGrad0(ix,iSpin)-dEcdGradN(ix,iSpin))>dEdGradTol) errorC=3
        end if
        if (errorX/=0 .or. errorC/=0) errorsFound = .true.

        ! Print comparison of analytical and numerical derivatives
        ! Do it only in first call, or if there is an error
        if (ix==0) then
          if (iCall==1 .or. errorX/=0) &
            print'(a4,a7,i3,i6,a6,3f15.9,i4)', &
              func(jf), auth(jf), ix, iSpin, 'Ex', &
              dExdDens0(iSpin), dExdDensN(iSpin), &
              dExdDens0(iSpin)-dExdDensN(iSpin), errorX
          if (iCall==1 .or. errorC/=0) &
            print'(a4,a7,i3,i6,a6,3f15.9,i4)', &
              func(jf), auth(jf), ix, iSpin, 'Ec', &
              dEcdDens0(iSpin), dEcdDensN(iSpin), &
              dEcdDens0(iSpin)-dEcdDensN(iSpin), errorC
        else
          if (iCall==1 .or. errorX/=0) &
            print'(a4,a7,i3,i6,a6,3f15.9,i4)', &
              func(jf), auth(jf), ix, iSpin, 'Ex', &
              dExdGrad0(ix,iSpin), dExdGradN(ix,iSpin), &
              dExdGrad0(ix,iSpin)-dExdGradN(ix,iSpin), errorX
          if (iCall==1 .or. errorC/=0) &
            print'(a4,a7,i3,i6,a6,3f15.9,i4)', &
              func(jf), auth(jf), ix, iSpin, 'Ec', &
              dEcdGrad0(ix,iSpin), dEcdGradN(ix,iSpin), &
              dEcdGrad0(ix,iSpin)-dEcdGradN(ix,iSpin), errorC
        end if

      end do ! iDeriv
      if (iCall==1) print*, ' '  ! Separate functionals for better visualization
    end do ! jf
  end do ! iCall

  ! Print error codes
  if (errorsFound) then
    print'(/,a)', 'Some errors found. Error codes:'
    print'(a)', 'error 1 => V=dE/dDens is positive'
    print'(a)', 'error 2 => analytical and numerical value of dE/dDens differ'
    print'(a)', &
      'error 3 => analytical and numerical value of dE/dGradDens differ'
  else
    print'(/,a)', 'No errors found.'
  end if

END PROGRAM siestaXCtest1
