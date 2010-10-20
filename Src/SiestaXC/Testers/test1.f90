PROGRAM siestaXCtest1

  ! Checks the consistency of the energies and its derivativas, returned by  
  ! the GGA (and LDA) XC routines contained in m_ggaxc and m_ldaxc modules.
  ! Two things are checked for each functional:
  !   That the xc potentials Vx=dEx/dDens and Vc=dEc/dDens are always negative.
  !   That the (analytical) returned values dE/dDens and dE/dGrad(Dens) coincide
  !   with their numerical counterparts, calculated as finite differences.
  ! J.M.Soler. Sept.2009

  ! Used module procedures and parameters
  USE siestaXC, only: dp
  USE siestaXC, only: ggaxc
  USE siestaXC, only: ldaxc

  implicit none

  ! Tester parameters
  integer, parameter:: irel    =  1      ! Relativistic? 0=>no, 1=>yes
  integer, parameter:: nCalls  = 100      ! Number of calls to each routine
  integer, parameter:: nfTot   = 10      ! Total number of functionals
  integer, parameter:: nSpin   =  4      ! Number of spin components (1|2|4)
  real(dp),parameter:: densMax = 0.1_dp  ! Upper limit to density value
  real(dp),parameter:: gradMax = 0.5_dp  ! Upper limit to density gradient
  real(dp),parameter:: deltaDens  = 1.e-6_dp  ! Finite diff. change
  real(dp),parameter:: deltaGrad  = 1.e-6_dp  ! Finite diff. change
  real(dp),parameter:: dEdDensTol = 1.e-6_dp  ! Tolerance for dE/dDens
  real(dp),parameter:: dEdGradTol = 1.e-6_dp  ! Tolerance for dE/dGadDens
  logical, parameter:: collinear  = .false.   ! Make grad(dens) collinear
                                              ! with density?

  ! List of functionals to be tested
!  integer, parameter:: nf = 10        ! Number of tested functionals
!  integer:: indexf(nf) = (/1,2,3,4,5,6,7,8,9,10/)  ! Indexes from list below
  integer, parameter:: nf = 8         ! Skip problematic functionals
  integer:: indexf(nf) = (/1,2,  4,5,6,  8,9,10/)  ! Indexes from list below

  ! All functionals available
  !                  1,       2,       3,       4,       5,   
  !                  6,       7,       8,       9       10
  character(len=3):: &
    func(nfTot) = (/'LDA',   'LDA',   'GGA',   'GGA',   'GGA',    &
                    'GGA',   'GGA',   'GGA',   'GGA',   'GGA'    /)
  character(len=6):: &
    auth(nfTot) = (/'PZ    ','PW92  ','PW91  ','PBE   ','RPBE  ', &
                    'revPBE','LYP   ','WC    ','PBESOL','AM05  ' /) 

  ! Tester variables and arrays
  complex(dp),parameter:: c0 = (0._dp,0._dp)  ! Complex zero
  complex(dp),parameter:: c1 = (1._dp,0._dp)  ! Complex one
  complex(dp),parameter:: ci = (0._dp,1._dp)  ! Complex i
  complex(dp):: densMat(2,2), gradDensMat(3,2,2), Imat(2,2), &
                pauli(3,2,2), polDirMat(2,2)
  real(dp):: epsC, epsX, dens(nSpin), &
             dEcdDens(nSpin), dEcdGrad(3,nSpin), &
             dExdDens(nSpin), dExdGrad(3,nSpin), densPol, densTot, &
             dVxdDens(nSpin,nSpin), dVcdDens(nSpin,nSpin), &
             Ec, Ex, gradDens(3,nSpin), gradDensPol(3), gradDensTot(3)
  real(dp):: cosTheta, phi, pi, ran(6), sinTheta, theta
  real(dp):: epsC0, epsX0, dens0(nSpin), gradDens0(3,nSpin)
  real(dp):: dEcdDens0(nSpin), dEcdGrad0(3,nSpin), &
             dExdDens0(nSpin), dExdGrad0(3,nSpin)
  real(dp):: dEcdDensN(nSpin), dEcdGradN(3,nSpin), &
             dExdDensN(nSpin), dExdGradN(3,nSpin)
  integer :: errorC, errorX, iCall, iDelta, iDeriv, iSpin, ix, &
             jf, kf, mSpin, one, two, three, four
  logical :: errorsFound
  integer :: mySeed(8) = (/39302655,65443109,09887367,99836827, &
                           70376268,32727926,46717826,35271017/)

  ! Initialize random seed
  call random_seed( put=mySeed )

  ! Set Pauli matrices
  Imat(1,:)    = (/c1, c0/)  ! Identity matrix
  Imat(2,:)    = (/c0, c1/)
  pauli(1,1,:) = (/c0, c1/)  ! PauliX
  pauli(1,2,:) = (/c1, c0/)
  pauli(2,1,:) = (/c0,-ci/)  ! PauliY
  pauli(2,2,:) = (/ci, c0/)
  pauli(3,1,:) = (/c1, c0/)  ! PauliZ
  pauli(3,2,:) = (/c0,-c1/)

  ! A silly trick to satisfy the compiler when nSpin<4
  one   = 1
  two   = 2
  three = 3
  four  = 4
  pi = acos(-1._dp)

  ! Print header
  print'(a4,a7,a3,2a6,3a15,a4)', &
    'func', 'auth  ', 'xyz', 'iSpin', 'Ex|Ec', &
    'deriv(anal)', 'deriv(num)', 'diff  ', 'err'

  ! Loop on calls with different densities
  errorsFound = .false.
  do iCall = 1,nCalls

    ! Choose random density, but large enough for safe numerical derivatives
    do ! density selection trial
      call random_number( ran(1:2) )
      densTot = densMax * ran(1)
      densPol = densTot * ran(2)
      if (densTot-densPol>10*deltaDens) exit ! density selection trial
    end do
!    densPol = 0                 ! Unpolarized, to compare with nSpin=1
    if (nSpin==1) then                   ! Non spin polarized
      dens0(one) = densTot
    else if (nSpin==2) then              ! Spin polarized
      dens0(one) = (densTot+densPol)/2
      dens0(two) = (densTot-densPol)/2
    else                                 ! Non collinear spin
      call random_number( ran(1:2) )
      cosTheta = 2*ran(1) - 1            ! Random spin direction
!      cosTheta = 0                      ! Point spin along x-y
!      cosTheta = 1                      ! Point spin along z
      sinTheta = sqrt(1-cosTheta**2)
      phi = 2*pi * ran(2)                ! Random spin direction
!      phi = 0                           ! Point spin along x
!      phi = pi/2                        ! Point spin along y
      polDirMat(:,:) = sinTheta*cos(phi)*pauli(1,:,:) &
                     + sinTheta*sin(phi)*pauli(2,:,:) &
                     + cosTheta*         pauli(3,:,:)
      densMat(:,:) = densTot/2 * Imat(:,:) + densPol/2 * polDirMat(:,:)
      dens0(one)   = real(densMat(1,1))
      dens0(two)   = real(densMat(2,2))
      dens0(three) = real(densMat(1,2))
      dens0(four)  = imag(densMat(1,2))
    end if

    ! Choose random gradient of density
    do iSpin = 1,nSpin
      call random_number( ran(1:3) )
      gradDens0(1:3,iSpin) = gradMax * (2*ran(1:3)-1)
    end do

    ! Choose density gradient collinear with density itself
    if (nSpin==4 .and. collinear) then
      call random_number( ran(1:6) )
      gradDensTot(1:3) = gradMax * (2*ran(1:3)-1)
      gradDensPol(1:3) = gradMax * (2*ran(4:6)-1)
      do ix = 1,3
        gradDensMat(ix,:,:) = gradDensTot(ix)/2 * Imat(:,:) &
                            + gradDensPol(ix)/2 * polDirMat(:,:)
        gradDens0(ix,one)   = real(gradDensMat(ix,1,1))
        gradDens0(ix,two)   = real(gradDensMat(ix,2,2))
        gradDens0(ix,three) = real(gradDensMat(ix,1,2))
        gradDens0(ix,four)  = imag(gradDensMat(ix,1,2))
      end do ! ix
    end if ! (nSpin==4)

    ! Print density (matrix)
!    if (nSpin==4) then
!      print'(a,/,(2f15.9,3x,2f15.9))', 'Density matrix =', &
!        complex(dens0(one),0._dp), complex(dens0(three),dens0(four)), &
!        complex(dens0(three),-dens0(four)), complex(dens0(two),0._dp)
!    else
!      print'(a,/,(2f15.9,3x,2f15.9))', 'Density =', dens0(1:nSpin)
!    end if

    ! Loop on functionals
    do kf = 1,nf
      jf = indexf(kf)

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

      ! Loop on variable with respect to which derivate
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

        ! Divide by 2 the iSpin=3,4 components for non-collinear spin because 
        ! the two nondiagonal elements of density matrix depend on them
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
              abs(dExdDens0(iSpin)-dExdDensN(iSpin)), errorX
          if (iCall==1 .or. errorC/=0) &
            print'(a4,a7,i3,i6,a6,3f15.9,i4)', &
              func(jf), auth(jf), ix, iSpin, 'Ec', &
              dEcdDens0(iSpin), dEcdDensN(iSpin), &
              abs(dEcdDens0(iSpin)-dEcdDensN(iSpin)), errorC
        else
          if (iCall==1 .or. errorX/=0) &
            print'(a4,a7,i3,i6,a6,3f15.9,i4)', &
              func(jf), auth(jf), ix, iSpin, 'Ex', &
              dExdGrad0(ix,iSpin), dExdGradN(ix,iSpin), &
              abs(dExdGrad0(ix,iSpin)-dExdGradN(ix,iSpin)), errorX
          if (iCall==1 .or. errorC/=0) &
            print'(a4,a7,i3,i6,a6,3f15.9,i4)', &
              func(jf), auth(jf), ix, iSpin, 'Ec', &
              dEcdGrad0(ix,iSpin), dEcdGradN(ix,iSpin), &
              abs(dEcdGrad0(ix,iSpin)-dEcdGradN(ix,iSpin)), errorC
        end if

!       Print density and gradient in case of error
!        if (errorX/=0 .or. errorC/=0) then
!          print'(a,4f15.8)','    Dens=', dens0
!          print'(a,4f15.8)','dDens/dx=', gradDens0(1,:)
!          print'(a,4f15.8)','dDens/dy=', gradDens0(2,:)
!          print'(a,4f15.8)','dDens/dz=', gradDens0(3,:)
!        end if

      end do ! iDeriv
      if (iCall==1) print*, ' ' ! Separate functionals for better visualization
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
