module m_ncps_translators

  public :: ncps_xml2froyen_new
  public :: ncps_xml2froyen
  interface ncps_xml2froyen
     module procedure ncps_xml2froyen_new
  end interface

CONTAINS

  subroutine ncps_xml2froyen_new( psxml, p, new_grid, a, b, rmax )

! Translate the more complete xml-adapted data structure 
! into the 'Froyen' ps type used in Atom and Siesta.
! Use accessors, and assume nothing about the grid in 
! the XML file. For this,
!
! p%nr, p%a, and p%b should be set on entry if we want
! a new grid, or if the grid in the file is not logarithmic.
!

        use m_ncps_xml_ps_t  !,     only: xml_ps_t
        use m_ncps_froyen_ps_t,  only: froyen_ps_t

        implicit none 

        integer, parameter  :: dp = selected_real_kind(14)


        type(xml_ps_t), intent(in)           :: psxml
        type(froyen_ps_t), intent(inout)     :: p
        logical, intent(in), optional        :: new_grid 
        real(dp), intent(in), optional       :: a
        real(dp), intent(in), optional       :: b
        real(dp), intent(in), optional       :: rmax


        logical          :: want_new_grid
        integer          :: position, i, il, ir, l, n
        character(len=1) :: ispp, lshell
        logical          :: polarized
        real(dp)         :: zeld(0:4), zelu(0:4)
        real(dp)         :: r2, rc, rmax_grid, znuc
        character(len=64):: xc_string
        character(len=40):: method_string
        character(len=1), dimension(0:4) :: &
                         sym = (/ "s", "p", "d", "f", "g" /)

        ! These are the "current" ATOM parameters
        ! There is another set which is turned on by the UCB_COMPAT flag in ATOM
        real(dp), parameter :: aa_def = 6.0_dp
        real(dp), parameter :: bb_def = 80.0_dp    ! UCB_COMPAT: 40.0
        real(dp), parameter :: rmax_def = 120.0_dp ! UCB_COMPAT: 80.0


        p%name = psxmlAtomicSymbol(psxml)
        p%zval         = psxmlPseudoZval(psxml)
        znuc           = psxmlAtomicNumber(psxml)
        ! This needs to be generalized
        p%gen_zval     = psxmlGenerationZval(psxml)
        
!
!       Need to include "universal" codes, such as those in LibXC
!
        xc_string = psxmlXCFunctional(psxml)

        select case(xc_string)
          case('Ceperley-Alder')
             p%icorr = 'ca'
          case('Wigner')
             p%icorr = 'wi'
          case('Hedin-Lundqvist')
             p%icorr = 'hl'
          case('Gunnarson-Lundqvist')
             p%icorr = 'gl'
          case('von Barth-Hedin')
             p%icorr = 'bh'
          case('Perdew-Burke-Ernzerhof')
             p%icorr = 'pb'
          case('RPBE - Hammer et al')
             p%icorr = 'rp'
          case('revPBE Zhang+Yang')
             p%icorr = 'rv'
          case('Becke-Lee-Yang-Parr')
             p%icorr = 'bl'
          case('Dion-et-al')
             p%icorr = 'vw'
          case('Wu-Cohen')
             p%icorr = 'wc'
          case('Perdew-Burke-Ernzerhof-solid')
             p%icorr = 'ps'
        end select

!
!       Note that most (all?) formats include only "scalar-relativistic" plus
!       maybe spin-orbit components, but never "polarized" pseudos.
!
        if (psxmlIsRelativistic(psxml)) then
           p%irel    = 'rel'
           ispp      = 'r'
           polarized = .false.
        else
           if (psxmlIsSpinPolarized(psxml)) then
              p%irel    = 'isp'
              ispp      = 's'
              polarized = .true.
           else
              p%irel    = 'nrl'
              ispp      = ' '
              polarized = .false.
           end if
        endif

        if (psxmlHasCoreCorrections(psxml)) then
            p%nicore = 'pcec'
         else
            p%nicore = 'nc'
         endif

!
!       Grid handling. Make sure that we cover the case in which
!       the file does not use a logarithmic grid.
!
        want_new_grid = .false.
        if (present(new_grid)) then
           want_new_grid = new_grid
        endif

        if (want_new_grid) then
           if (.not. present(a)) call die("new grid: a not present")
           if (.not. present(b)) call die("new grid: b not present")
           if (present(rmax)) then
              rmax_grid = rmax
              if (rmax == 0.0_dp) rmax_grid = psxmlGridRmax(psxml)
           else
              rmax_grid = psxmlGridRmax(psxml)
           endif
           p%a = a
           p%b = b
           p%nr = nint(log(rmax_grid/b+1.0d0)/a)

        else
           if (.not. psxmlHasGlobalLogGrid(psxml)) then
              print *, "Using ATOM defaults for log grid..."
              ! use the ATOM defaults 
              ! Note that a and b are interchanged...
              p%a = 1.0_dp / bb_def
              p%b = exp(-aa_def)/znuc
              rmax_grid = rmax_def
              p%nr = nint(log(rmax_grid/(p%b)+1.0d0)/(p%a))
              ! call die("Do not have logarithmic grid...")
           else
              p%nr = psxmlGridNpoints(psxml)
              p%a  = psxmlLogGridStep(psxml)
              p%b  = psxmlLogGridScale(psxml)
           endif
           
        endif

        p%method(1)    = psxmlCreator(psxml)
        p%method(2)    = psxmlDate(psxml)
        method_string = psxmlPseudoFlavor(psxml)
        read(method_string,'(4a10)') (p%method(i),i=3,6) 

        p%npotu        = psxmlPotentialsUp(psxml)
        p%npotd        = psxmlPotentialsDown(psxml)

! Allocate the radial variables and semilocal potentials

        ! Note that p's functions start at r=0, whereas
        ! the file's do not include r=0

        p%nrval        = p%nr + 1

        allocate(p%r(1:p%nrval))
        allocate(p%chcore(1:p%nrval))
        allocate(p%chval(1:p%nrval))

        if (p%npotd.gt.0) then
           allocate(p%vdown(1:p%npotd,1:p%nrval))
           allocate(p%ldown(1:p%npotd))
        endif

        if (p%npotu.gt.0) then
           allocate(p%vup(1:p%npotu,1:p%nrval))
           allocate(p%lup(1:p%npotu))
        endif
! ---

! Calculate the points of the logarithmic radial grid 
        ! The first point here (NOT in the files produced by ATOM) is r=0
        do ir = 1, p%nrval
           p%r(ir) = p%b * (exp(p%a*(ir-1))-1)
        enddo
! ---

! Translate the valence charge density and the pseudo-core charge density,
! and define the value at the first point of the logarithmic grid

        if (psxmlHasCoreCorrections(psxml)) then
           do ir = 2, p%nrval
              p%chcore(ir) = psxmlEvaluateCoreCharge(psxml,p%r(ir))
           enddo
        else
           p%chcore(2:p%nrval) = 0.0_dp
        endif

        do ir = 2, p%nrval
           p%chval(ir) = psxmlEvaluateValenceCharge(psxml,p%r(ir))
        enddo
        r2=p%r(2)/(p%r(3)-p%r(2))
        p%chcore(1) = p%chcore(2) - r2*(p%chcore(3)-p%chcore(2))
        p%chval(1) = p%chval(2) - r2*(p%chval(3)-p%chval(2))
! ---

        zeld(:) = 0.0d0
        zelu(:) = 0.0d0

        do il = 1, p%npotd
          p%ldown(il) = psxmlPotAngMomentum(psxml,"d",il)
          zeld(p%ldown(il)) = psxmlOccupation(psxml,"d",il)
          do ir = 2, p%nrval
             p%vdown(il,ir) = p%r(ir) * &
                           psxmlEvaluatePotential(psxml,"d",il,p%r(ir))
          enddo
          p%vdown(il,1) = p%vdown(il,2) - r2*(p%vdown(il,3)-p%vdown(il,2))
        enddo

        do il = 1, p%npotu
           p%lup(il) = psxmlPotAngMomentum(psxml,"u",il)
           zelu(p%lup(il)) = psxmlOccupation(psxml,"u",il)
           do ir = 2, p%nrval
              p%vup(il,ir) = p%r(ir) * &
                          psxmlEvaluatePotential(psxml,"u",il,p%r(ir))
           enddo
          p%vup(il,1) = p%vup(il,2) - r2*(p%vup(il,3)-p%vup(il,2))
        enddo

        ! Encode generation configuration and cutoffs
        p%text = ' '
        position = 1
        do il = 1, p%npotd
           n = psxmlPrincipalN(psxml,"d",il)
           l = psxmlPotAngMomentum(psxml,"d",il)
           rc = psxmlGenerationCutoff(psxml,"d",il)
           if ( .not. polarized) then
              write(p%text(position:),9070)     &
                   n, sym(l), zeld(l)+zelu(l),ispp, rc
 9070         format(i1,a1,f5.2,a1,' r=',f5.2,'/')
              position = position + 17
           else
              write(p%text(position:),9090)    &
                   n, sym(l), zeld(l),zelu(l),ispp, rc
 9090         format(i1,a1,f4.2,',',f4.2,a1,f4.2,'/')
              position = position + 17
           end if
        enddo

      end subroutine ncps_xml2froyen_new
    end module m_ncps_translators
