module m_ncps_translators

  public :: ncps_xml2froyen, ncps_xml2froyen_new

CONTAINS

  subroutine ncps_xml2froyen( psxml, p )

! Translate the more complete xml-adapted data structure 
! into the 'Froyen' ps type used in Atom and Siesta.

        use m_ncps_xml_ps_t,     only: xml_ps_t
        use m_ncps_froyen_ps_t,  only: froyen_ps_t

        implicit none 

        integer, parameter  :: dp = selected_real_kind(14)


        type(xml_ps_t), intent(in)           :: psxml
        type(froyen_ps_t), intent(out)       :: p

        integer          :: position, i, il, ir, lmax, lshellint
        character(len=1) :: ispp, lshell
        logical          :: polarized
        real(dp)         :: zeld(0:4), zelu(0:4)
        real(dp)         :: r2


        p%name         = psxml%header%symbol
        p%nr           = psxml%pot(1)%V%grid%npts
        p%nrval        = p%nr + 1
        p%zval         = psxml%header%zval
        p%gen_zval     = psxml%header%gen_zval
! relativistic and correlation are not needed in Siesta, so they are
! not translated
        select case(psxml%header%xcfunctionalparametrization)
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

        select case(psxml%header%relativistic)
          case(.true.)
            p%irel    = 'rel'
            ispp      = 'r'
            polarized = .false.
          case(.false.)
            select case(psxml%header%polarized)
              case(.true.)
                p%irel    = 'isp'
                ispp      = 's'
                polarized = .true.
              case(.false.)
                p%irel    = 'nrl'
                ispp      = ' '
                polarized = .false.
            end select
        end select

        select case(psxml%header%core_corrections)
          case("yes")
            p%nicore = 'pcec'
          case("no")
            p%nicore = 'nc'
        end select

        p%a            = psxml%pot(1)%V%grid%step
        p%b            = psxml%pot(1)%V%grid%scale

        p%method(1)    = psxml%header%creator
        p%method(2)    = psxml%header%date
        read(psxml%header%flavor,'(4a10)') (p%method(i),i=3,6) 

        p%npotu        = psxml%npots_up
        p%npotd        = psxml%npots_down

! Allocate the radial variables and semilocal potentials
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
        do ir = 1, p%nrval
           p%r(ir) = p%b * (exp(p%a*(ir-1))-1)
        enddo
! ---

! Translate the valence charge density and the pseudo-core charge density,
! and define the value at the first point of the logarithmic grid
        p%chcore(2:p%nrval) = psxml%core_charge%data(1:p%nr)
        p%chval(2:p%nrval)  = psxml%valence_charge%data(1:p%nr)
        r2=p%r(2)/(p%r(3)-p%r(2))
        p%chcore(1) = p%chcore(2) - r2*(p%chcore(3)-p%chcore(2))
        p%chval(1) = p%chval(2) - r2*(p%chval(3)-p%chval(2))
! ---

        zeld(:) = 0.0d0
        zelu(:) = 0.0d0

        do il = 1, p%npotd
          lshell = psxml%pot(il)%l
          select case(lshell)
            case('s') 
              lshellint = 0
            case('p') 
              lshellint = 1
            case('d') 
              lshellint = 2
            case('f') 
              lshellint = 3
          end select
          p%ldown(il) = lshellint
          p%vdown(il,2:p%nrval) = psxml%pot(il)%V%data(1:p%nr)
          zeld(lshellint) = psxml%pot(il)%occupation
        enddo

        do il = 1, p%npotu
          lshell = psxml%pot(p%npotd+il)%l
          select case(lshell)
            case('s') 
              lshellint = 0
            case('p') 
              lshellint = 1
            case('d') 
              lshellint = 2
            case('f') 
              lshellint = 3
          end select
          p%lup(il) = lshellint
          p%vup(il,2:p%nrval) = psxml%pot(p%npotd+il)%V%data(1:p%nr)
          zelu(lshellint) = psxml%pot(p%npotd+il)%occupation
        enddo

        p%text = ' '
        position = 1
        lmax = max(p%npotd, p%npotu)
        do il = 1, lmax
           if ( .not. polarized) then
              write(p%text(position:),9070)     &
                                           psxml%pot(il)%n, &
                                           psxml%pot(il)%l, &
                                           zeld(il-1)+zelu(il-1), &
                                           ispp,  &
                                           psxml%pot(il)%cutoff
 9070         format(i1,a1,f5.2,a1,' r=',f5.2,'/')
              position = position + 17
           else
              write(p%text(position:),9090)    &
                                           psxml%pot(il)%n, &
                                           psxml%pot(il)%l, &
                                           zeld(il-1),      &
                                           zelu(il-1),      &
                                           ispp,            &
                                           psxml%pot(il)%cutoff
 9090         format(i1,a1,f4.2,',',f4.2,a1,f4.2,'/')
              position = position + 17
           end if
        enddo

       
! For debugging -------------------------------------------------------
!        write(6,*)p%name
!        write(6,*)p%nr
!        write(6,*)p%nrval
!        write(6,*)p%zval
!        write(6,*)p%icorr
!        write(6,*)p%irel
!        write(6,*)p%nicore
!        write(6,*)p%a
!        write(6,*)p%b
!        write(6,*)p%method(1)
!        write(6,*)p%method(2)
!        write(6,*)p%method(3)
!        write(6,*)p%method(4)
!        write(6,*)p%method(5)
!        write(6,*)p%method(6)
!        write(6,*)p%text
!        write(6,*)p%npotu
!        write(6,*)p%npotd
!        write(6,'(4f20.12)')p%r(1:4)
!        write(6,'(4f20.12)')p%r(101:104)
!        write(6,'(4f20.12)')p%vdown(1,1:4)
!        write(6,'(4f20.12)')p%vdown(2,1:4)
!        write(6,'(4f20.12)')p%vdown(3,1:4)
!        write(6,'(4f20.12)')p%vdown(4,1:4)
!        write(6,'(4f20.12)')p%vup(1,1:4)
!        write(6,'(4f20.12)')p%vup(2,1:4)
!        write(6,'(4f20.12)')p%vup(3,1:4)
!        write(6,'(4f20.12)')p%chcore(1:4)
!        write(6,'(4f20.12)')p%chval(1:4)
!        write(6,*)p%ldown
!        write(6,*)p%lup
!         do ir = 1, p%nrval
!           write(6,*)p%r(ir), pseudo%pswf(1)%V%data(ir)
!         enddo
! -----
      end subroutine ncps_xml2froyen


  subroutine ncps_xml2froyen_new( psxml, p, new_grid )

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

        logical          :: want_new_grid
        integer          :: position, i, il, ir, l, n
        character(len=1) :: ispp, lshell
        logical          :: polarized
        real(dp)         :: zeld(0:4), zelu(0:4)
        real(dp)         :: r2, rc
        character(len=64):: xc_string
        character(len=40):: method_string
        character(len=1), dimension(0:4) :: &
                         sym = (/ "s", "p", "d", "f", "g" /)


!        p%name         = psxml%header%symbol
        p%name = psxmlAtomicSymbol(psxml)
!        p%zval         = psxml%header%zval
        p%zval         = psxmlPseudoZval(psxml)
!        p%gen_zval     = psxml%header%gen_zval
        p%gen_zval     = psxmlGenerationZval(psxml)
        
! relativistic and correlation are not needed in Siesta, so they are
! not translated
        xc_string = psxmlXCFunctional(psxml)
!        select case(psxml%header%xcfunctionalparametrization)
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

!        select case(psxml%header%relativistic)
!          case(.true.)
        if (psxmlIsRelativistic(psxml)) then
           p%irel    = 'rel'
           ispp      = 'r'
           polarized = .false.
        else
!            select case(psxml%header%polarized)
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

!        select case(psxml%header%core_corrections)
        if (psxmlHasCoreCorrections(psxml)) then
            p%nicore = 'pcec'
         else
            p%nicore = 'nc'
         endif

        want_new_grid = .false.
        if (present(new_grid)) then
           want_new_grid = new_grid
        endif

        if (psxmlHasGlobalLogGrid(psxml)) then
!        These have to come from outside
!!        p%nr           = psxml%pot(1)%V%grid%npts
!        p%a            = psxml%pot(1)%V%grid%step
!        p%b            = psxml%pot(1)%V%grid%scale
           if (want_new_grid) then
              ! do nothing, assume we have entered nr, a, and b
           else
              p%nr = psxmlGridNpoints(psxml)
              p%a  = psxmlLogGridStep(psxml)
              p%b  = psxmlLogGridScale(psxml)
           endif
        else
           if (want_new_grid) then
              ! do nothing
           else
              call die("Do not have logarithmic grid...")
           endif
        endif

!        p%method(1)    = psxml%header%creator
!        p%method(2)    = psxml%header%date
!        read(psxml%header%flavor,'(4a10)') (p%method(i),i=3,6) 
        p%method(1)    = psxmlCreator(psxml)
        p%method(2)    = psxmlDate(psxml)
        method_string = psxmlPseudoFlavor(psxml)
        read(method_string,'(4a10)') (p%method(i),i=3,6) 

!        p%npotu        = psxml%npots_up
!        p%npotd        = psxml%npots_down
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

!        p%chcore(2:p%nrval) = psxml%core_charge%data(1:p%nr)
        do ir = 2, p%nrval
           p%chcore(ir) = psxmlEvaluateCoreCharge(psxml,p%r(ir))
           p%chval(ir) = psxmlEvaluateValenceCharge(psxml,p%r(ir),debug=.true.)
        enddo
        r2=p%r(2)/(p%r(3)-p%r(2))
        p%chcore(1) = p%chcore(2) - r2*(p%chcore(3)-p%chcore(2))
        p%chval(1) = p%chval(2) - r2*(p%chval(3)-p%chval(2))
! ---

        zeld(:) = 0.0d0
        zelu(:) = 0.0d0

        do il = 1, p%npotd
!          p%ldown(il) = lshellint
!          zeld(lshellint) = psxml%pot(il)%occupation
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
