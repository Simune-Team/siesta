      program pseudooperator

! Stand alone program to:
! 1. Read the pseudopotential files
! 2. Generate the local part and the Kleynman-Bylander projector
!
! Written by Javier Junquera using code from atom.F in Siesta.
! Re-written by A. Garcia to use the PSML library
!
! Input required:
!
!     The name of the file where the pseudopotential is stored 
!     (here we assume that this name can be constructed from the systemlabel)
!     The angular momentum cutoff for KB non-local pseudopotential
!     (here lmxkb)         
!     The number of KB projectors per angular momentum shell
!     The reference energy for the calculation of the KB projectors
!
! Output:
!
!     the pseudopotential read, of type pseudopotential_t

      use precision,        only: dp
      use atmparams,        only: nrmax, nkbmx
      use m_ncps, only: pseudopotential_t => froyen_ps_t, pseudo_read
      use m_pseudooperator, only: radii_ps, vlocal1, vlocal2, KBgen
      use SiestaXC, only: xc_id_t, get_xc_id_from_atom_id, setXC
      use SiestaXC, only: atomxc

      implicit none

!
!     INPUT VARIABLES REQUIRED TO GENERATE THE LOCAL PART AND KB PROJECTORS
!     THEY MUST BE PROVIDED BY THE USER
!  
      character(len=30)       :: systemlabel ! System label, used to identify
                                             !   the file where the pseudo is
                                             !   stored
      character*10            :: functl      ! Exchange and correlation function
      character*10            :: author      ! Exchange and correlation parametr
      integer                 :: lmxkb       ! Angular momentum cutoff for 
                                             !   Kleinman-Bylander nonlocal 
                                             !   pseudopotential
      integer,  allocatable    :: nkbl(:)    ! Number of KB projectors for
                                             !   each angular momentum
      real(dp), allocatable    :: erefkb(:,:)! Reference energies (in Ry) for
                                             !   the calculation of the KB
                                             !   projectors
      integer                 :: is          ! Species index               

!
!     DERIVED VARIABLE WHERE THE PSEUDO IS STORED
!
      type(pseudopotential_t) :: psr
      type(xc_id_t)           :: xc_id
      integer                 :: status

!
!     VARIABLES RELATED WITH THE RADIAL LOGARITHMIC GRID 
!
      integer                  :: nrval      ! Number of points required to 
                                             !   store the pseudopotential and
                                             !   the wave functions in the
                                             !   logarithmic grid
                                             !   (directly read from the 
                                             !   pseudopotential file)
      real(dp)                 :: a          ! Step parameter of log. grid
                                             !   (directly read from the 
                                             !   pseudopotential file)
      real(dp)                 :: b          ! Scale parameter of log. grid
                                             !   (directly read from the 
                                             !   pseudopotential file)
      real(dp), allocatable    :: rofi(:)    ! Radial points of the 
                                             !   logarithmic grid 
                                             !   rofi(r)=b*[exp(a*(i-1)) - 1]
                                             !   (directly read from the 
                                             !   pseudopotential file)
      real(dp), allocatable    :: drdi(:)    ! Derivative of the radial 
                                             !   distance respect the mesh index
                                             !   Computed after the radial mesh
                                             !    is read
      real(dp), allocatable    :: s(:)       ! Metric array
                                             !   Computed after the radial mesh
                                             !    is read
      real(dp)                 :: rpb, ea    ! Local variables used in the 
                                             !   calculation of the log. grid

!
!     VARIABLE USED TO STORE THE SEMILOCAL COMPONENTS OF THE PSEUDOPOTENTIAL
!
      real(dp), allocatable    :: vps(:,:)   ! Semilocal components of the
                                             !   pseudopotentials 
                                             !   (directly read from the 
                                             !   pseudopotential file)
!
!     Local variables
!  
      real(dp)                 :: Zval       ! Valence charge of the atom   
                                             !   (directly read from the 
                                             !   pseudopotential file)
      integer                  :: nrgauss    ! Number of points in the logarith.
                                             !   grid required to describe the
                                             !   KB projectors
                                             !   Output of radii_ps.
      real(dp)                 :: rgauss     ! Approximately the maximum cut-off
                                             !   radius used in the
                                             !   pseudopotential generation.
                                             !   Output of radii_ps.
      real(dp)                 :: rgauss2    ! Radius where the pseudopotentials
                                             !   reach  the asymptotic behaviour
                                             !   2*Zval/r.
!
!     VARIABLES READ FROM THE PSEUDO FILE
!
      character*4              ::  nicore    ! Flag that determines whether
                                             !   non-linear core corrections
                                             !   are included
      character*3              ::  irel      ! Flag that determines whether
                                             !   the atomic calculation is 
                                             !   relativistic or not
!
!     LOCAL VARIABLES TO DEFINE THE LOCAL PART OF THE PSEUDOPOTENTIAL
!
      integer                  :: nchloc     ! Number of radial points required
                                             !   to describe chlocal
      real(dp)                 :: rchloc     ! Radius where all the semilocal
                                             !   pseudopotentials get the 
                                             !   asymptotic behaviour
      real(dp), allocatable    :: vlocal(:)  ! Local component of the pseudopot.
                                             !   Output of vlocal1 or vlocal2.
      real(dp), allocatable    :: chlocal(:) ! Charge distribution that 
                                             !   generates vlocal
                                             !   Output of vlocal1 or vlocal2.
!
!     LOCAL VARIABLES TO DEFINE THE KLEINMAN-BYLANDER PROJECTORS 
!
      real(dp), allocatable    :: rho(:)     ! Valence charge density 
                                             !   As read from the pseudo file,
                                             !   it is angularly integrated
                                             !   (i.e. multiplied by 4*pi*r^2).
      real(dp), allocatable    :: ve(:)      ! Electrostatic potential
                                             !   generated by the valence charge
                                             !   density, readed from the 
                                             !   pseudo file
      real(dp), allocatable    :: chcore(:)  ! Core charge density 
                                             !   As read from the pseudo file,
                                             !   it is angularly integrated
                                             !   (i.e. multiplied by 4*pi*r^2).
      real(dp), allocatable    :: auxrho(:)  !  Sum of the valence charge and 
                                             !   core charge (if NLCC included)
                                             !   densities to compute the 
                                             !   atomic exchange and correl.
                                             !   potential. 
                                             !   auxrho is NOT angularly integr.
                                             !   (not multiplied by 4*pi*r^2)
      integer                  :: irelt      ! Flag that determines whether the
                                             !   atomic calculation to
                                             !   generate the pseudopotential
                                             !   was relativistic (irelt = 1)
                                             !   or no relativistic (irelt = 0)
      real(dp), allocatable    :: vxc(:) ! Exchange and correlation potential
!     ps%gen_zval                            ! Valence charge of the pseudoion 
                                             !   for which the pseudo was 
                                             !   generated in the ATM code
                                             !   (it might not correspond with 
                                             !   the nominal valence charge
                                             !   of the atom if the pseudo 
                                             !   has been generated for an ionic
                                             !   configuration, for instance 
                                             !   when semicore has been
                                             !   explicitly included in the 
                                             !   valence).
                                             !   For instance, for Ba with 
                                             !   the semicore in valence,
                                             !   (atomic reference configuration
                                             !   5s2 5p6 5d0 4f0),
                                             !   chgvps = 8  (two in the 5s 
                                             !                and six in the 5p)
                                             !   zval   = 10 (two in the 5s, 
                                             !                six in the 5p, 
                                             !                and two in the 6s.
                                             !   These two last not included in
                                             !   reference atomic configuration)
      real(dp)                 :: ex         ! Total exchange energy 
      real(dp)                 :: ec         ! Total correlation energy 
      real(dp)                 :: dx         ! IntegralOf( rho * (eps_x - v_x) )
      real(dp)                 :: dc         ! IntegralOf( rho * (eps_c - v_c) )
!
!     LOCAL VARIABLES
!
      integer                  :: ndown, ir  ! Counters for the do loops
      integer                  :: l          ! Angular momentum of the channel
                                             !   that is read 
      real(dp)                 :: pi         ! pi
      real(dp)                 :: r2         ! Local variables
      integer                  :: nkb        ! Number of KB projectors

      pi = dacos(-1.0_dp)
!
!     DEFINE THE INPUT VARIABLES:
!     This must be done by the user before compiling and running the code
!
      is          = 1
!     Define the system label, used to identify the file where the pseudo is 
!     stored (here systemlabel.psf or systemlabel.vps)
      systemlabel = "Name"

!     Define the angular momentum cutoff for Kleinman-Bylander nonlocal pseudopo
!     In this example, we will expand up to the f-shell (l=3).
!     Therefore, we will include (lmxkb + 1) shells
!     l = 0 (s)
!     l = 1 (p)
!     l = 2 (d)
!     l = 3 (f)
      lmxkb       = 2
!     Define the number of KB projectors for each angular momentum
      allocate( nkbl(0:lmxkb) )
      nkbl(:) = 1
!     Define the reference energies (in Ry) for the calculation of the KB proj.
      allocate( erefkb(nkbmx,0:lmxkb) )
      erefkb(:,:) = huge(1.0_dp)
!
!     READ THE PSEUDOPOTENTIAL FILE
!
      call pseudo_read( systemlabel, psr)
!
!     STORE IN LOCAL VARIABLES SOME OF THE PARAMETERS READ IN THE 
!     PSEUDOPOTENTIAL, AND DEFINITION OF THE RADIAL LOGARITHMIC GRID
!
      nrval         = psr%nrval
      a             = psr%a
      b             = psr%b

      allocate( rofi(nrmax) )
      allocate( drdi(nrmax) )
      allocate( s(nrmax)    )

      rofi(1:nrval) = psr%r(1:nrval)

!     Calculate drdi and s
!     drdi is the derivative of the radial distance respect to the mesh index
!     i.e. rofi(ir)= b*[ exp( a*(i-1) ) - 1 ] and therefore
!     drdi=dr/di =a*b*exp(a*(i-1))= a*[rofi(ir)+b]

      rpb = b
      ea  = dexp(a)
      do ir = 1, nrval
        drdi(ir) = a * rpb
        s(ir)    = dsqrt( a * rpb )
        rpb      = rpb * ea
      enddo

!!     For debugging
!      do ir = 1, nrval
!        write(6,'(i5,3f20.12)') ir, rofi(ir), drdi(ir), s(ir)
!      enddo
!!     End debugging

! 
!     STORE THE IONIC PSEUDOPOTENTIALS IN A LOCAL VARIABLE 
!     Only the 'down' component is used
!
      allocate( vps(nrmax,0:lmxkb) )

      do ndown = 1, lmxkb + 1
         l = psr%ldown(ndown)
         if( l .ne. ndown-1) then
           write(6,'(a)')
     .       'atom: Unexpected angular momentum  for pseudopotential'
           write(6,'(a)')
     .       'atom: Pseudopotential should be ordered by increasing l'
         endif
         vps(1:nrval,l) = psr%vdown(ndown,1:nrval)
!        The semilocal components of the pseudopotential are read in r*V format.
!        Here we compute the pseudopotential part dividing by r.
         do ir = 2, nrval
           vps(ir,l) = vps(ir,l) / rofi(ir)
         enddo
         vps(1,l) = vps(2,l)     ! AG
      enddo

!     Storing locally other variables
      Zval   = psr%zval
      nicore = psr%nicore
      irel   = psr%irel

      call get_xc_id_from_atom_id(psr%icorr,xc_id,status)
      if (status == 0) then
         functl = xc_id%siestaxc_id%family
         author = xc_id%siestaxc_id%authors
      else
         print *, "**** Cannot process XC info ***"
         functl = "LDA"
         author = "PZ"
      endif
      call setxc(1,(/functl/),(/author/),(/1.0_dp/),(/1.0_dp/))
!
!     COMPUTE THE LOCAL PART OF THE PSEUDOPOTENTIAL
!
!     Rgauss is approximately the maximum cut-off radius used in the
!     pseudopotential generation.
!     Rgauss is determined by comparison of the pseudopot.
!     Corresponding to different l (this is not possible if we have
!     just one pseudopotential shell)
!     Rgauss2 is the radius where the pseudopotentials reach  the
!     asymptotic behaviour 2*Zval/r.
!     For just one pseudopotential Rgauss is taken equal to Rgauss2
!
      call radii_ps( vps, rofi, Zval, nrval, lmxkb,
     .               nrgauss, rgauss, rgauss2 )

! 
!     Calculate local pseudopotential
!
      allocate( vlocal(nrmax)  )
      allocate( chlocal(nrmax) )

      if( rgauss2 .gt. 1.3_dp * rgauss ) then
        write(6,'(a)') 'Using large-core scheme for Vlocal'

!       In this case the atom core is so big that we do not have an asymptotic
!       of 2*Zval/r until Rgauss2 (> Rc) . To retain the same asymptotic
!       behaviour as in the pseudopotentials we modify the definition
!       of the local potential, making it join the Vps's smoothly at rgauss.
!
        write(6,'(/,a,f10.5)') 'atom: Estimated core radius ', rgauss2
        if( nicore .eq. 'nc ')
     .    write(6,'(/,2a)')'atom: Including non-local core corrections',
     .                     ' could be a good idea'

!       As all potentials are equal beyond rgauss, we can just use the
!       s-potential here.
        call vlocal2( Zval, nrval, a, rofi, drdi, s, vps(:,0),
     .                nrgauss, vlocal, nchloc, chlocal )

      else
!       In this case the pseudopotential reach to asymptotic behaviour 2*Zval/r
!       for a radius approximately equal to Rc. We build a generalized-gaussian
!       "local charge density" and set Vlocal as the potential generated by
!       it. Note that chlocal is negative.
        call vlocal1( Zval, nrval, a, rofi, drdi, s, rgauss,
     .                vlocal, nchloc, chlocal )
      endif

!
! Save local-pseudopotential charge
! 
      rchloc=rofi(nchloc)
      write(6,'(2a,f10.5)') 'atom: Maximum radius for' ,
     .  ' 4*pi*r*r*local-pseudopot. charge ',rchloc

!!     For debugging
!      do ir = 1, nrval
!        write(6,*) rofi(ir), vlocal(ir), chlocal(ir)
!      enddo
!!     End debugging

!
!     COMPUTE THE NON-LOCAL KLEINMAN-BYLANDER PROJECTORS
!
!     Allocate internal variables
!
      allocate( ve(nrmax)     )
      allocate( vxc(nrmax)    )
      allocate( rho(nrmax)    )
      allocate( auxrho(nrmax) )
      allocate( chcore(nrmax) )

!     Read the valence charge density from the pseudo file
!     and scale it if the ionic charge of the reference configuration
!     is not the same as the nominal valence charge of the atom

!     AG: What are we trying to achieve here?

      do ir = 1, nrval
        rho(ir) = psr%gen_zval * psr%chval(ir)/zval
      enddo

!     Find the Hartree potential created by a radial electron density
!     using the Numerov's method to integrate the radial Poisson equation.
!     The input charge density at this point has to be angularly integrated.
      call vhrtre( rho, ve, rofi, drdi, s, nrval, a )

!     Read the core charge density from the pseudo file
      chcore(1:nrval) = psr%chcore(1:nrval)

!     Compute the exchange and correlation potential in the atom
!     Note that atomxc expects true rho(r), not 4 * pi * r^2 * rho(r)
!     We use auxrho for that
!
      do ir = 2, nrval
        r2 = rofi(ir)**2
        r2 = 4.0_dp * pi * r2
        dc = rho(ir) / r2
        if( nicore .ne. 'nc  ')  dc = dc + chcore(ir) / r2
        auxrho(ir) = dc
      enddo

      r2        = rofi(2) / (rofi(3)-rofi(2))
      auxrho(1) = auxrho(2) - ( auxrho(3) - auxrho(2) ) * r2

!     Determine whether the atomic calculation to generate the pseudopotential
!     is relativistic or not
      if (irel.eq.'rel') irelt=1
      if (irel.ne.'rel') irelt=0

      call atomxc( irelt, nrval, nrmax, rofi,
     .             1, auxrho, ex, ec, dx, dc, vxc )

!     Add the exchange and correlation potential to the Hartree potential
      ve(1:nrval) = ve(1:nrval) + vxc(1:nrval)
      
!!     For debugging
!      do ir = 1, nrval
!        write(6,'(3f20.12)') rofi(ir), auxrho(ir), ve(ir)
!      enddo
!!     End debugging

!
!     Redefine the array s for the Schrodinger equation integration 
!
      s(2:nrval) = drdi(2:nrval) * drdi(2:nrval)
      s(1) = s(2)

!
!     Calculation of the Kleinman-Bylander projector functions
!
      call KBgen( is, a, b, rofi, drdi, s,
     .            vps, vlocal, ve, nrval, Zval, lmxkb,
     .            nkbl, erefkb, nkb )

      deallocate( rofi    )
      deallocate( drdi    )
      deallocate( s       )
      deallocate( vps     )
      deallocate( vlocal  )
      deallocate( chlocal )
      deallocate( ve      )
      deallocate( vxc     )
      deallocate( rho     )
      deallocate( auxrho  )
      deallocate( erefkb  )
      deallocate( nkbl    )

      end program pseudooperator
