

      subroutine write_psml( ray, npotd, npotu, zion )

      use flib_wxml
      use SiestaXC, only: xc_id_t, get_xc_id_from_atom_id
      use SiestaXC, only: xc_nfuncs_libxc

      implicit none

      include 'param.h'
      include 'radial.h'
      include 'ion.h'
      include 'orbital.h'
      include 'charge.h'
      include 'pseudowave.h'
      include 'corecorr.h'

      type(xmlf_t)  :: xf
      type(xc_id_t) :: xc_id

      integer npotd, npotu
      double precision  :: zion

      character*4      :: polattrib, relattrib, coreattrib
      character*10     :: ray(6)
      character*30 xcfuntype, xcfunparam
      character*30 gridtype, gridunits, gridscale, gridstep

      integer                         :: ivps, ip, i
      double precision, allocatable   :: chval(:), f(:)

      double precision                :: total_valence_charge
!
      integer :: stat
      character(len=132) :: line, msg

      allocate(f(1:nr))

! Digest and dump the information about the exchange and correlation functional

      call get_xc_id_from_atom_id(icorr,xc_id,stat)
      if (stat /= 0) then
         stop "Wrong icorr code!"
      endif

! Digest and dump the information about the pseudopotential flavor
      select case(irel) 

        case('isp') 
          polattrib   = 'yes'
          relattrib   = 'no'

        case('rel') 
          polattrib   = 'no'
          relattrib   = 'yes'

        case('nrl') 
          polattrib   = 'no'
          relattrib   = 'no'

      end select

! Digest and dump the information about the non-linear core corrections
      select case(nicore) 

        case('pcec', 'fcec', 'fche', 'pche') 
          coreattrib  = 'yes'

        case default
          coreattrib  = 'no'

      end select

! Digest and dump the information about the grid
      gridtype    = 'log'
      gridunits   = 'bohr'
      gridscale   = str(a)
      gridstep    = str(b)
      ! Note that in ATOM r(1) = 0.
      

! Allocate and define the valence charge density
      allocate(chval(1:nr))

      ! Note that we do not renormalize the charge 
      ! to make it integrate to zion
      do ip = 2, nr
        chval(ip) = (cdd(ip)+cdu(ip))
      enddo

! ---------------------------------------------------------------------
                                                                               
      call xml_OpenFile("PSML",xf, indent=.false.)

      call xml_AddXMLDeclaration(xf,"UTF-8")

      call xml_NewElement(xf,"psml")
      call my_add_attribute(xf,"version","0.7")
      call my_add_attribute(xf,"energy_unit","hartree")
      call my_add_attribute(xf,"length_unit","bohr")


      call xml_NewElement(xf,"provenance")
      call my_add_attribute(xf,"creator",ray(1))
      call my_add_attribute(xf,"date",ray(2))
      call xml_NewElement(xf,"input-file")
      call my_add_attribute(xf,"name","INP")
!
!     Note that a file already connected to one unit
!     must not be re-opened with another unit...
!     Since INP is still open at this time, we use
!     INP_COPY (generated in atm.f)
!
      open(44,file="INP_COPY",form="formatted",status="old",
     $     position="rewind",action="read")
      do
         read(44,fmt="(a)",iostat=stat) line
         if (stat .ne. 0) exit
         call xml_AddPcData(xf,trim(line),line_feed=.true.)
      enddo
      close(44)
!
      call xml_EndElement(xf,"input-file")
      call xml_EndElement(xf,"provenance")

        call xml_NewElement(xf,"header")
          call my_add_attribute(xf,"atomic-label",nameat)
          call my_add_attribute(xf,"atomic-number",str(znuc))
          call my_add_attribute(xf,"z-pseudo",str(zion))
          call my_add_attribute(xf,"flavor",ray(3)//ray(4))
          call my_add_attribute(xf,"relativistic",relattrib)
          call my_add_attribute(xf,"polarized",polattrib)
          call my_add_attribute(xf,"core-corrections",coreattrib)
 
          call xml_NewElement(xf,"exchange-correlation")
          call xml_NewElement(xf,"annotation")
          call my_add_attribute(xf,"atom-xc-code",icorr)
          call my_add_attribute(xf,"siesta-xc-type",
     $                          trim(xc_id%siestaxc_id%family))
          call my_add_attribute(xf,"siesta-xc-authors",
     $                          trim(xc_id%siestaxc_id%authors))
          call xml_EndElement(xf,"annotation")

          call xml_NewElement(xf,"libxc-info")
          call my_add_attribute(xf,"number-of-functionals",
     $                          str(xc_nfuncs_libxc(xc_id)))
           do i = 1, xc_nfuncs_libxc(xc_id)
              call xml_NewElement(xf,"functional")
               call my_add_attribute(xf,"name",
     $                      trim(xc_id%libxc_id(i)%name))
               call my_add_attribute(xf,"type",
     $                      trim(xc_id%libxc_id(i)%xc_kind%str))
               call my_add_attribute(xf,"id",
     $                      str(xc_id%libxc_id(i)%code))
              call xml_EndElement(xf,"functional")
           enddo
          call xml_EndElement(xf,"libxc-info")
          call xml_EndElement(xf,"exchange-correlation")
          !

        ! This is a child element of <header> for now, as
        ! implied in the psml paper.
        call do_configuration(total_valence_charge)

        call xml_EndElement(xf,"header")


        call xml_NewElement(xf,"grid")
          call my_add_attribute(xf,"npts",str(nr))
          write(msg,"(3(a,':',a,'/'))")  "type","log",
     $                                   "scale",trim(gridscale),
     $                                   "step",trim(gridstep)
         call my_add_attribute(xf,"annotation",trim(msg))

          call xml_NewElement(xf,"grid-data")
            call xml_AddArray(xf,r(1:nr))
          call xml_EndElement(xf,"grid-data")
        call xml_EndElement(xf,"grid")

        call xml_NewElement(xf,"semilocal-potentials")
          call my_add_attribute(xf,"npots-major",str(npotd))
          call my_add_attribute(xf,"npots-minor",str(npotu))
  
! Down pseudopotentials follows
! (Major)

      vpsd: do ivps = 1, lmax
           if (indd(ivps) .eq. 0) cycle
           call xml_NewElement(xf,"vps")
             call my_add_attribute(xf,"set","major")
             call my_add_attribute(xf,"n",str(no(indd(ivps))))
             call my_add_attribute(xf,"l",il(ivps))
             call my_add_attribute(xf,"rc",str(rc(ivps)))
             call my_add_attribute(xf,"flavor",ray(3)//ray(4))

             call xml_NewElement(xf,"radfunc")
               call xml_NewElement(xf,"data")
                 call remove_r(viod(ivps,:),r(:),f(:))
                 call xml_AddArray(xf, 0.5d0 * f(1:nr))
!               call xml_AddArray(xf, 0.5d0 * viod(ivps,2:nr))
               call xml_EndElement(xf,"data")
             call xml_EndElement(xf,"radfunc")
           call xml_EndElement(xf,"vps")
         enddo vpsd

! Up pseudopotentials follows
! (Minor)
         vpsu: do ivps = 1, lmax
           if (indu(ivps) .eq. 0) cycle
           call xml_NewElement(xf,"vps")
             call my_add_attribute(xf,"set","minor")
             call my_add_attribute(xf,"n",str(no(indu(ivps))))
             call my_add_attribute(xf,"l",il(ivps))
             call my_add_attribute(xf,"rc",str(rc(ivps)))
             call my_add_attribute(xf,"flavor",ray(3)//ray(4))

             call xml_NewElement(xf,"radfunc")

               call xml_NewElement(xf,"data")
                 call remove_r(viou(ivps,:),r(:),f(:))
                 call xml_AddArray(xf, 0.5d0 * f(1:nr))
!                 call xml_AddArray(xf, 0.5d0 * viou(ivps,2:nr))
               call xml_EndElement(xf,"data")
             call xml_EndElement(xf,"radfunc")
           call xml_EndElement(xf,"vps")
        enddo vpsu
        call xml_EndElement(xf,"semilocal-potentials")

! Dump of the pseudowave functions
        call xml_NewElement(xf,"pseudo-wave-functions")
          call my_add_attribute(xf,"npswfs",str(nshells_stored))
          call my_add_attribute(xf,"npswfs-major",str(nshells_stored))
          call my_add_attribute(xf,"npswfs-minor",str(0))
  
        ! We only write the "test" wavefunctions ("major" = "average")
        pswfd: do i = 1, nshells_stored
           call xml_NewElement(xf,"pswf")
             call my_add_attribute(xf,"n",str(n_pswf(i)))
             call my_add_attribute(xf,"l",il(l_pswf(i)+1))
             call my_add_attribute(xf,"set","major")
             call xml_NewElement(xf,"radfunc")

               call xml_NewElement(xf,"data")
                 call remove_r(pswf(i,:),r(:),f(:))
                 call xml_AddArray(xf,force_underflow(f(1:nr)))
               call xml_EndElement(xf,"data")
             call xml_EndElement(xf,"radfunc")
           call xml_EndElement(xf,"pswf")
        enddo pswfd

        call xml_EndElement(xf,"pseudo-wave-functions")

        call xml_NewElement(xf,"valence-charge")
          call my_add_attribute(xf,"total-charge",
     $                    str(total_valence_charge))
          call xml_NewElement(xf,"radfunc")

            call xml_NewElement(xf,"data")
            call remove_r2(chval(:),r(:),f(:))
            call xml_AddArray(xf,force_underflow(f(1:nr)))
            call xml_EndElement(xf,"data")
          call xml_EndElement(xf,"radfunc")
        call xml_EndElement(xf,"valence-charge")

        if (coreattrib(1:3) .eq. "yes") then

           call xml_NewElement(xf,"pseudocore-charge")
           call my_add_attribute(xf,"matching-radius",str(rc_core))
           call my_add_attribute(xf,"number-of-continuous-derivatives",
     $                               str(n_of_continuous_derivs))
           if (n_of_continuous_derivs == 1) then
              call my_add_attribute(xf,"annotation",
     $             "original Froyen-Cohen-Louie method")
           else
              call my_add_attribute(xf,"annotation",
     $             "three-parameter Martins method")
           endif
           call xml_NewElement(xf,"radfunc")

           call xml_NewElement(xf,"data")
            call remove_r2(cdc(:),r(:),f(:))
            call xml_AddArray(xf,force_underflow(f(1:nr)))
           call xml_EndElement(xf,"data")
           call xml_EndElement(xf,"radfunc")
           call xml_EndElement(xf,"pseudocore-charge")
        endif

        call xml_EndElement(xf,"psml")
      call xml_Close(xf)

      deallocate(chval)

      CONTAINS

      subroutine do_configuration(total_valence_charge)
      integer, parameter :: dp = selected_real_kind(10,100)

      real(dp), intent(out) :: total_valence_charge
      integer  :: i, lp      
      real(dp) :: occ_down, occ_up, occupation

      call  get_total_valence_charge(total_valence_charge)

      call xml_NewElement(xf,"valence-configuration")

        ! this call is needed here before creating any sub-elements...
        call my_add_attribute(xf,"total-valence-charge",
     $           str(total_valence_charge))

      i = ncore
      do 
         i = i + 1
         if (i > norb) exit
         lp = lo(i) + 1
         occ_down = zo(i)
         occ_up   = 0.0_dp
         if (split_shell(lp)) then
            i = i + 1
            occ_up = zo(i)
         endif
         occupation = occ_down + occ_up
         if (occupation .lt. 1.0e-10_dp) cycle
         call xml_NewElement(xf,"shell")
         call my_add_attribute(xf,"n",str(no(i)))
         call my_add_attribute(xf,"l",il(lp))
         call my_add_attribute(xf,"occupation",str(occupation))
         if (polarized .and. split_shell(lp)) then
            call my_add_attribute(xf,"occupation-down",str(occ_down))
            call my_add_attribute(xf,"occupation-up",str(occ_up))
         endif
         call xml_EndElement(xf,"shell")
      enddo
      call xml_EndElement(xf,"valence-configuration")
      end subroutine do_configuration

      subroutine get_total_valence_charge(tot_occ)
      integer, parameter :: dp = selected_real_kind(10,100)
      real(dp), intent(out) :: tot_occ 

      integer :: i, lp
      real(dp):: occ_down, occ_up

      tot_occ = 0.0_dp

      i = ncore
      do 
         i = i + 1
         if (i > norb) exit
         lp = lo(i) + 1
         occ_down = zo(i)
         occ_up   = 0.0_dp
         if (split_shell(lp)) then
            i = i + 1
            occ_up = zo(i)
         endif
         tot_occ = tot_occ + occ_down + occ_up
      enddo
      end subroutine get_total_valence_charge

      logical function split_shell(lp)
      integer, intent(in) :: lp


      split_shell = .false.
      if (polarized) then
         split_shell = .true.
      else if (relativistic) then
         if (lp /= 1) split_shell = .true.
      endif
      end function split_shell
         
      subroutine my_add_attribute(xf,name,value)
      type(xmlf_t), intent(inout)   :: xf
      character(len=*), intent(in)  :: name
      character(len=*), intent(in)  :: value

       call xml_AddAttribute(xf,name,trim(value))
      end subroutine my_add_attribute

      subroutine remove_r(rf,r,f)
      ! Removes a factor of r from rf to get f
      ! At the same time, it extrapolates to zero (r(1))

      double precision, intent(in)  :: rf(:), r(:)
      double precision, intent(out) :: f(:)
      
      integer i
      double precision :: r2

      do i = 2, nr
         f(i) = rf(i)/r(i)
      enddo

      r2 = r(2)/(r(3)-r(2))
      f(1) = f(2) - (f(3)-f(2))*r2

      end subroutine remove_r

      subroutine remove_r2(rf,r,f)
      ! Removes a factor of r**2 from rf to get f
      ! At the same time, it extrapolates to zero (r(1))

      double precision, intent(in)  :: rf(:), r(:)
      double precision, intent(out) :: f(:)
      
      integer i
      double precision :: r2

      do i = 2, nr
         f(i) = rf(i)/(r(i)*r(i))
      enddo

      r2 = r(2)/(r(3)-r(2))
      f(1) = f(2) - (f(3)-f(2))*r2

      end subroutine remove_r2

      pure elemental function force_underflow(x)
      double precision, intent(in) ::  x
      double precision  force_underflow
C
C     Avoid very small numbers that might need a three-character
C     exponent field in formatted output
C      
      if (abs(x) .lt. 1.0d-99) then
         force_underflow = 0.0d0
      else
         force_underflow = x
      endif

      end function force_underflow

      end subroutine write_psml



