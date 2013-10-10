module m_ncps_parsing_helpers
!
!  This module reads a pseudopotential file written in XML.
!  A full example of the building up of a data structure using
!  the SAX paradigm.
!
use flib_sax
use m_ncps_xml_ps_t        ! Data types

private

!
! It defines the routines that are called from xml_parser in response
! to particular events.
!
public  :: begin_element, end_element, pcdata_chunk
private :: die

logical, private  :: in_vps = .false. , in_radfunc = .false.
logical, private  :: in_semilocal = .false. , in_header = .false.
logical, private  :: in_coreCharge = .false. , in_data = .false.
logical, private  :: in_grid_data = .false. , in_grid = .false.
logical, private  :: in_valenceCharge = .false.
logical, private  :: in_pseudowavefun = .false. , in_pswf = .false.

integer, private, save  :: ndata, ndata_grid

integer, parameter, private    :: dp = selected_real_kind(14)
real(dp), private, save :: zval_generation

type(xml_ps_t), public, target, save :: pseudo

type(grid_t), private, save, pointer  :: grid => null()
!
! Pointers to make it easier to manage the data
!
type(header_t), private, pointer   :: hp => null()
type(vps_t), private, pointer      :: pp => null()
type(pswf_t), private, pointer     :: pw => null()
type(radfunc_t), private, pointer  :: rp => null()



CONTAINS  !===========================================================

!----------------------------------------------------------------------
subroutine begin_element(name,attributes)
character(len=*), intent(in)    :: name
type(dictionary_t), intent(in)  :: attributes

character(len=100)  :: value
integer             :: status


select case(name)

      case ("pseudo")
         pseudo%npots  = 0
         pseudo%npswfs = 0

!         call get_value(attributes,"version",value,status)
!         if (value == "0.5") then
!            print *, "Processing a PSEUDO version 0.5 XML file"
!            pseudo%npots  = 0
!            pseudo%npswfs = 0
!            global_grid%npts = 0
!         else
!            print *, "Can only work with PSEUDO version 0.5 XML files"
!            STOP
!         endif

      case ("header")
         in_header = .true.
         hp => pseudo%header
         
         call get_value(attributes,"symbol",hp%symbol,status)
         if (status /= 0 ) call die("Cannot determine atomic symbol")

         call get_value(attributes,"zval",value,status)
         if (status /= 0 ) call die("Cannot determine zval")
         read(unit=value,fmt=*) hp%zval

         call get_value(attributes,"xc-functional-parametrization", &
                        hp%xcfunctionalparametrization,status)
         if (status /= 0 ) &
            call die("Cannot determine xc-functional-parametrization ")

         call get_value(attributes,"creator",hp%creator,status)
         if (status /= 0 ) hp%creator="unknown"

         call get_value(attributes,"date",hp%date,status)
         if (status /= 0 ) hp%date="unknown"

         call get_value(attributes,"flavor",hp%flavor,status)
         if (status /= 0 ) hp%flavor="unknown"

         call get_value(attributes,"relativistic",value,status)
         if (status /= 0 ) value = "no"
         hp%relativistic = (value == "yes")

         call get_value(attributes,"polarized",value,status)
         if (status /= 0 ) value = "no"
         hp%polarized = (value == "yes")

         call get_value(attributes,"core-corrections", &
                                    hp%core_corrections,status)
         if (status /= 0 ) hp%core_corrections = "nc"

      case ("vps")
         in_vps = .true.

         pseudo%npots = pseudo%npots + 1
         pp => pseudo%pot(pseudo%npots)
         rp => pp%V                       ! Pointer to radial function

         call get_value(attributes,"l",pp%l,status)
         if (status /= 0 ) call die("Cannot determine l for Vps")

         call get_value(attributes,"principal-n",value,status)
         if (status /= 0 ) call die("Cannot determine n for Vps")
         read(unit=value,fmt=*) pp%n

         call get_value(attributes,"cutoff",value,status)
         if (status /= 0 ) call die("Cannot determine cutoff for Vps")
         read(unit=value,fmt=*) pp%cutoff

         call get_value(attributes,"occupation",value,status)
         if (status /= 0 ) call die("Cannot determine occupation for Vps")
         read(unit=value,fmt=*) pp%occupation
         zval_generation = zval_generation + pp%occupation

         call get_value(attributes,"spin",value,status)
         if (status /= 0 ) call die("Cannot determine spin for Vps")
         read(unit=value,fmt=*) pp%spin

      case ("grid")

         in_grid = .true.
         allocate(grid)   ! Will forget about previous allocation

         call get_value(attributes,"type",grid%type,status)
         if (status /= 0 ) call die("Cannot determine grid type")

         call get_value(attributes,"npts",value,status)
         if (status /= 0 ) call die("Cannot determine grid npts")
         read(unit=value,fmt=*) grid%npts

         call get_value(attributes,"scale",value,status)
         if (status /= 0 ) call die("Cannot determine grid scale")
         read(unit=value,fmt=*) grid%scale

         call get_value(attributes,"step",value,status)
         if (status /= 0 ) call die("Cannot determine grid step")
         read(unit=value,fmt=*) grid%step

         !
         ! In this way we allow for a private grid for each radfunc,
         ! or for a global grid specification
         !
         if (in_radfunc) then
            rp%grid => grid
         else
            print *, "Associated global grid"
            pseudo%global_grid => grid
         endif

      case ("data")
         if (.not. in_radfunc) STOP "<data> element outside <rad_func> element"
         in_data = .true.
         if (.not. associated(rp%grid)) STOP "Cannot find grid data for radfunc"
         if (rp%grid%npts == 0) STOP "Grid not specified correctly"
         allocate(rp%data(rp%grid%npts))
         ndata = 0             ! To start the build up

      case ("grid_data")
         if (.not. in_grid) STOP "Grid_data element outside grid element"
         in_grid_data = .true.
         if (grid%npts == 0) STOP "Grid npts attribute not specified correctly"
         allocate(grid%grid_data(grid%npts))
         ndata_grid = 0             ! To start the build up

      case ("radfunc")
         in_radfunc = .true.
         rp%grid => pseudo%global_grid    ! Might be null
                                          ! There should then be a local grid element
                                          ! read later

      case ("pseudocore-charge")
         in_coreCharge = .true.
         rp => pseudo%core_charge

      case ("valence-charge")
         in_valenceCharge = .true.
         rp => pseudo%valence_charge

      case ("semilocal")
         in_semilocal = .true.
         zval_generation = 0.0_dp

         call get_value(attributes,"npots-down",value,status)
         if (status /= 0 ) call die("Cannot determine npots-down")
         read(unit=value,fmt=*) pseudo%npots_down

         call get_value(attributes,"npots-up",value,status)
         if (status /= 0 ) call die("Cannot determine npots-up")
         read(unit=value,fmt=*) pseudo%npots_up

      case ("pseudowave-functions")
         in_pseudowavefun = .true. 

      case ("pswf")
         in_pswf = .true.

         pseudo%npswfs = pseudo%npswfs + 1

         pw => pseudo%pswf(pseudo%npswfs)
         rp => pw%V                       ! Pointer to radial function

         call get_value(attributes,"l",pw%l,status)
         if (status /= 0 ) call die("Cannot determine l for Vps")
                                                                                 
         call get_value(attributes,"principal-n",value,status)
         if (status /= 0 ) call die("Cannot determine n for Vps")
         read(unit=value,fmt=*) pw%n
                                                                                 
         call get_value(attributes,"spin",value,status)
         if (status /= 0 ) call die("Cannot determine spin for Vps")
         read(unit=value,fmt=*) pw%spin

end select

end subroutine begin_element
!----------------------------------------------------------------------

subroutine end_element(name)
character(len=*), intent(in)     :: name

select case(name)

      case ("vps")
         in_vps = .false.

      case ("radfunc")
         in_radfunc = .false.

      case ("grid")
         in_grid = .false.

      case ("data")
      !
      ! We are done filling up the radfunc data
      ! Check that we got the advertised number of items
      !
         in_data = .false.
         if (ndata /= size(rp%data)) STOP "npts mismatch in radfunc data"

      case ("grid_data")
      !
      ! We are done filling up the grid data
      ! Check that we got the advertised number of items
      !
         in_grid_data = .false.
         if (ndata_grid /= size(grid%grid_data)) STOP "npts mismatch in grid"

      case ("pseudocore-charge")
         in_coreCharge = .false.

      case ("valence-charge")
         in_valenceCharge = .false.

      case ("semilocal")
         in_semilocal = .false.
         pseudo%header%gen_zval = zval_generation

      case ("pseudowave-functions")
         in_pseudowavefun = .false. 

      case ("pswf")
         in_pswf = .false.

      case ("pseudo")
!         call dump_pseudo(pseudo)

end select

end subroutine end_element
!----------------------------------------------------------------------

subroutine pcdata_chunk(chunk)
character(len=*), intent(in) :: chunk


if (len_trim(chunk) == 0) RETURN     ! skip empty chunk

if (in_data) then
!
! Note that we know where we need to put it through the pointer rp...
!
      call build_data_array(chunk,rp%data,ndata)

else if (in_grid_data) then
!
!     Fill the explicit grid data
      call build_data_array(chunk,grid%grid_data,ndata_grid)

else if (in_header) then
      !
      ! There should not be any pcdata in header in this version...

      print *, "Header data:"
      print *, trim(chunk)

endif

end subroutine pcdata_chunk
!----------------------------------------------------------------------

      subroutine die(str)
      character(len=*), intent(in), optional   :: str
      if (present(str)) then
         write(unit=0,fmt="(a)") trim(str)
      endif
      write(unit=0,fmt="(a)") "Stopping Program"
      stop
      end subroutine die


end module m_ncps_parsing_helpers











