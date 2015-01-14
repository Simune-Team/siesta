module m_psml_parsing_helpers
!
!  This module reads a pseudopotential file written in XML (PSML format)
!  A full example of the building up of a data structure using
!  the SAX paradigm.
!
 use m_psml_core                   ! For data types
 use external_interfaces, only: die

implicit none

private

!
! It defines the routines that are called from xml_parser in response
! to particular events.
!
public  :: begin_element, end_element, pcdata_chunk

!
! The data will be stored in this public variable
! There are some design issues to decide:
! -- Should this be a pointer, associated by the client
!    program to its own variable? In that case, the
!    client should make sure that the variable is "clean"
!    before calling this routine, as some fields will be
!    allocated here.
! -- Perhaps it should be a pointer allocated here (and
!    then destroyed when done by the client). It should be
!    allocated at the beginning of processing, maybe detected
!    with a (default) "begin_Document" handler, or by 
!    "begin_Element" after checking for association.
!    This is the cleanest option, as the caller might want
!    to keep several instances alive at the same time...
!    (... but this should be handled by the user)
! -- If "pseudo" here is a normal variable, it should also
!    be "cleaned" before the next use. The current usage
!    in Abinit falls in this category: psxml is a pointer
!    associated to "pseudo", and cleaned after use.
!
!    We implement the first option now

type(ps_t), pointer, public, save :: pseudo => null() 
logical, public, save             :: debug_parsing = .false.

logical, private, save  :: in_vps = .false. , in_radfunc = .false.
logical, private, save  :: in_config_val = .false.
logical, private, save  :: in_semilocal = .false. , in_header = .false.
logical, private, save  :: in_coreCharge = .false. , in_data = .false.
logical, private, save  :: in_grid_data = .false. , in_grid = .false.
logical, private, save  :: in_valenceCharge = .false.
logical, private, save  :: in_provenance = .false.
logical, private, save  :: in_valence_config = .false.
logical, private, save  :: in_xc = .false., in_libxc_info = .false.
logical, private, save  :: in_pseudowavefun = .false. , in_pswf = .false.
logical, private, save  :: in_psoperator = .false. , in_projectors = .false.
logical, private, save  :: in_proj = .false. , in_local_potential = .false.
logical, private, save  :: got_explicit_grid_data

integer, private, save  :: ndata, ndata_grid
integer, private, save  :: n_funct

integer, parameter, private    :: dp = selected_real_kind(14)
real(dp), private, save        :: zval_generation

type(grid_t), private, save, pointer  :: grid => null()
!
! Pointers to make it easier to manage the data
!
type(provenance_t), private, pointer      :: pp => null()
type(header_t), private, pointer          :: hp => null()
type(config_val_t), private, pointer      :: cp => null()
type(xc_t), private, pointer              :: xp => null()
type(pswfs_t), private, pointer           :: wfp => null()
type(semilocal_t), private, pointer       :: slp => null()
type(psoperator_t), private, pointer      :: pop => null()
type(valence_charge_t), private, pointer  :: valp => null()
type(core_charge_t), private, pointer     :: corep => null()
type(radfunc_t), private, pointer         :: rp => null()

CONTAINS  !===========================================================

!----------------------------------------------------------------------
#ifndef PSML_USE_FOX
subroutine begin_element(name,attributes)
use xmlf90_sax, only: dictionary_t, get_value
#else
subroutine begin_element(namespaceURI,localName,name,attributes)
use Fox_sax, only: dictionary_t
use fox_extra, only: get_value=>get_value_by_key
character(len=*), intent(in)    :: namespaceURI
character(len=*), intent(in)    :: localName
#endif
character(len=*), intent(in)    :: name
type(dictionary_t), intent(in)  :: attributes

character(len=100)  :: value
integer             :: status

integer             :: i

if (debug_parsing) print *, "Element: ", trim(name)

select case(name)

      case ("psml")

         ! Make sure that pseudo is pointing to something

         if (.not. associated(pseudo)) then
            call die("ps_t object not initialized by client")
         endif

         call get_value(attributes,"version",value,status)
         if (value /= "0.7") then
            call die("Can only work with PSML version 0.7 files")
         endif

      case ("provenance")
         in_provenance = .true.
         if (in_psoperator) then
            pp => pseudo%psoperator%provenance
         else
            pp => pseudo%provenance
         endif

         call get_value(attributes,"creator",pp%creator,status)
         if (status /= 0 ) pp%creator="unknown"
 
         call get_value(attributes,"date",pp%date,status)
         if (status /= 0 ) pp%date="unknown"

      case ("header")
         in_header = .true.
         hp => pseudo%header
         
         call get_value(attributes,"atomic-label",hp%atomic_label,status)
         if (status /= 0 ) call die("Cannot determine atomic-label")

         call get_value(attributes,"z-pseudo",value,status)
         if (status /= 0 ) call die("Cannot determine z-pseudo")
         read(unit=value,fmt=*) hp%zpseudo

         call get_value(attributes,"atomic-number",value,status)
         if (status /= 0 ) call die("Cannot determine atomic number")
         read(unit=value,fmt=*) hp%z

         call get_value(attributes,"flavor",hp%flavor,status)
         if (status /= 0 ) hp%flavor="ATTEND"

         call get_value(attributes,"relativistic",value,status)
         if (status /= 0 ) value = "no"
         hp%relativistic = (value == "yes")

         call get_value(attributes,"polarized",value,status)
         if (status /= 0 ) value = "no"
         hp%polarized = (value == "yes")
         if (hp%polarized .and. hp%relativistic) then
            call die("Cannot be polarized and relativistic at the same time")
         endif

         call get_value(attributes,"core-corrections", &
                                    hp%core_corrections,status)
         if (status /= 0 ) hp%core_corrections = "no"

      case ("exchange-correlation")
         in_xc = .true.
         xp => pseudo%xc_info

      case ("libxc-info")
         if (.not. in_xc) call die("Orphan <libxc-info>")
         in_libxc_info = .true.
         call get_value(attributes,"number-of-functionals", &
                                    value,status)
         if (status /= 0 ) call die("Error reading number of libxc functs")
         read(unit=value,fmt=*)  xp%n_functs_libxc 
         if (xp%n_functs_libxc /= 2 ) then
            call die("Non-conventional number of libxc functionals")
         endif
         n_funct = 0

      case ("functional")
         if (.not. in_libxc_info) call die("Orphan <functional>")
         n_funct = n_funct + 1
         if (n_funct > 2) call die("Too many libxc functionals")

         call get_value(attributes,"name", &
                                    xp%libxc_name(n_funct),status)
         if (status /= 0 ) call die("Error reading libxc name")

         call get_value(attributes,"id", value, status)
         if (status /= 0 ) call die("Error reading libxc id")
         read(unit=value,fmt=*)  xp%libxc_id(n_funct) 

      case ("valence-configuration")
         in_valence_config = .true.

         pseudo%config_val%nshells = 0
         call get_value(attributes,"total-valence-charge",value,status)
         if (status /= 0 ) call die("Cannot determine total-valence-charge")
         read(unit=value,fmt=*) pseudo%config_val%total_charge

      case ("shell")

         if (in_valence_config) then
            cp => pseudo%config_val
!         else if (in_core_config) then
!            cp => pseudo%config_core
         else
            call die("Orphan <shell> element")
         endif

         cp%nshells = cp%nshells + 1
         call get_value(attributes,"l",cp%l(cp%nshells),status)
         if (status /= 0 ) call die("Cannot determine l for shell")

         call get_value(attributes,"n",value,status)
         if (status /= 0 ) call die("Cannot determine n for shell")
         read(unit=value,fmt=*) cp%n(cp%nshells)

         call get_value(attributes,"occupation",value,status)
         if (status /= 0 ) call die("Cannot determine occupation for shell")
         read(unit=value,fmt=*) cp%occ(cp%nshells)

         call get_value(attributes,"occupation-up",value,status)
         if (status == 0 ) then
            read(unit=value,fmt=*) cp%occ_up(cp%nshells)
         endif
         call get_value(attributes,"occupation-down",value,status)
         if (status == 0 ) then
            read(unit=value,fmt=*) cp%occ_down(cp%nshells)
         endif

      case ("vps")
         in_vps = .true.
         if (.not. in_semilocal) call die("Orphan <vps> element")

         slp => pseudo%semilocal
         slp%npots = slp%npots + 1
         i = slp%npots
         rp => slp%V(i)

         call get_value(attributes,"l",slp%l(i),status)
         if (status /= 0 ) call die("Cannot determine l for Vps")

         call get_value(attributes,"n",value,status)
         if (status /= 0 ) call die("Cannot determine n for Vps")
         read(unit=value,fmt=*) slp%n(i)

         call get_value(attributes,"rc",value,status)
         if (status /= 0 ) call die("Cannot determine rc for Vps")
         read(unit=value,fmt=*) slp%rc(i)

         call get_value(attributes,"set",slp%set(i),status)
         if (status /= 0 ) call die("Cannot determine set for Vps")

         call get_value(attributes,"flavor",slp%flavor(i),status)
         if (status /= 0 ) call die("Cannot determine flavor for Vps")

      case ("proj")
         in_proj = .true.
         if (.not. in_projectors) call die("Orphan <proj> element")

         pop => pseudo%psoperator
         pop%nprojs = pop%nprojs + 1
         i = pop%nprojs
         rp => pop%proj(i)

         call get_value(attributes,"l",pop%l(i),status)
         if (status /= 0 ) call die("Cannot determine l for proj")

         call get_value(attributes,"n",value,status)
         if (status /= 0 ) call die("Cannot determine n for proj")
         read(unit=value,fmt=*) pop%n(i)

         call get_value(attributes,"ekb",value,status)
         if (status /= 0 ) call die("Cannot determine Ekb for proj")
         read(unit=value,fmt=*) pop%ekb(i)

         call get_value(attributes,"set",pop%set(i),status)
         if (status /= 0 ) call die("Cannot determine set for proj")

         call get_value(attributes,"type",pop%type(i),status)
         if (status /= 0 ) call die("Cannot determine type of proj")

      case ("pswf")

         if (.not. in_pseudowavefun) call die("Orphan <pswf> element")
         in_pswf = .true.

         wfp => pseudo%pswfs
         wfp%npswfs = wfp%npswfs + 1
         i = wfp%npswfs
         rp => wfp%Phi(i)

         call get_value(attributes,"l",wfp%l(i),status)
         if (status /= 0 ) call die("Cannot determine l for PSwf")
                                                                              
         call get_value(attributes,"n",value,status)
         if (status /= 0 ) call die("Cannot determine n for PSwf")
         read(unit=value,fmt=*) wfp%n(i)

         call get_value(attributes,"set",value,status)
         if (status /= 0 ) call die("Cannot determine set for PSwf")
         read(unit=value,fmt=*) wfp%set(i)

      case ("grid")
         in_grid = .true.
         allocate(grid)   ! Will forget about previous allocation

         got_explicit_grid_data = .false.

         ! This attribute is mandatory
         call get_value(attributes,"npts",value,status)
         if (status /= 0 ) call die("Cannot determine grid npts")
         read(unit=value,fmt=*) grid%npts

         !
         ! In this way we allow for a private grid for each radfunc,
         ! or for a global grid specification
         !
         if (in_radfunc) then
            if (debug_parsing) print *, "Found grid in radfunc"
            if (associated(rp%grid)) then
               call die("psml: Two grids specified for a radfunc")
            endif
            rp%grid => grid

         ! We check whether we are at the top level,
         ! and not at a level that allows a "regional" grids

         else if (in_psoperator) then

            if (associated(pop%grid)) then
               call die("psml: Two psoperator grids specified")
            endif
            if (debug_parsing) print *, "Found psoperator grid"
            pop%grid => grid

         else  ! We are at the top level

            if (debug_parsing) print *, "Found grid at the top level"
            if (associated(pseudo%global_grid)) then
               call die("psml: Two global grids specified")
            endif
            pseudo%global_grid => grid
         endif

      case ("data")
         if (.not. in_radfunc) then
            call die("<data> element outside <rad_func> element")
         endif
         in_data = .true.
         if (.not. associated(rp%grid)) then
            ! Try regional and global grids...
            if (in_psoperator) then
               if (associated(pop%grid)) then
                  rp%grid => pop%grid
               endif
            else ! try global
               if (associated(pseudo%global_grid)) then
                  rp%grid => pseudo%global_grid
               endif
            endif
         endif
         if (.not. associated(rp%grid)) &
               call die("Cannot find grid data for radfunc")
         if (rp%grid%npts == 0) call die("Grid not specified correctly")
         allocate(rp%data(rp%grid%npts))
         ndata = 0             ! To start the build up

      case ("grid-data")
         if (.not. in_grid) call die("Grid_data element outside grid element")
         in_grid_data = .true.
         got_explicit_grid_data = .true.
         if (grid%npts == 0) call die("Grid npts attribute faulty")
         allocate(grid%grid_data(grid%npts))
         ndata_grid = 0             ! To start the build up

      case ("radfunc")
         in_radfunc = .true.

      case ("pseudocore-charge")
         in_coreCharge = .true.
         corep => pseudo%core_charge
         rp => corep%rho_core

         call get_value(attributes,"matching-radius",value,status)
         if (status /= 0 ) call die("Cannot determine radius for pseudocore")
         read(unit=value,fmt=*) corep%rcore
                                                                              
         call get_value(attributes,"number-of-continuous-derivatives", &
                                   value,status)
         if (status /= 0 )  &
               call die("Cannot determine n-cont-derivs for pseudocore")
         read(unit=value,fmt=*) corep%n_cont_derivs

      case ("valence-charge")
         in_valenceCharge = .true.
         valp => pseudo%valence_charge
         rp => valp%rho_val

         call get_value(attributes,"total-charge",value,status)
         if (status /= 0 ) call die("Cannot determine total valence charge")
         read(unit=value,fmt=*) valp%total_charge
                                                                              
      case ("semilocal-potentials")
         in_semilocal = .true.
         slp => pseudo%semilocal
         slp%npots = 0

         call get_value(attributes,"npots-major",value,status)
         if (status /= 0 ) call die("Cannot determine npots-major")
         read(unit=value,fmt=*) slp%npots_major

         call get_value(attributes,"npots-minor",value,status)
         if (status /= 0 ) call die("Cannot determine npots-minor")
         read(unit=value,fmt=*) slp%npots_minor

      case ("pseudopotential-operator")
         in_psoperator = .true.
         pop => pseudo%psoperator

      case ("projectors")
         in_projectors = .true.
         pop => pseudo%psoperator
         pop%nprojs = 0

      case ("local-potential")
         in_local_potential = .true.
         pop => pseudo%psoperator
         rp => pop%vlocal

         call get_value(attributes,"kind",pop%vlocal_kind,status)
         if (status /= 0 ) call die("Cannot determine kind of vlocal")
                                                                              
      case ("pseudo-wave-functions")
         in_pseudowavefun = .true. 

         wfp => pseudo%pswfs
         wfp%npswfs = 0

         call get_value(attributes,"npswfs-major",value,status)
         if (status /= 0 ) call die("Cannot determine npswfs-major")
         read(unit=value,fmt=*) wfp%npswfs_major

         call get_value(attributes,"npswfs-minor",value,status)
         if (status /= 0 ) call die("Cannot determine npswfs-minor")
         read(unit=value,fmt=*) wfp%npswfs_minor

      case ("annotation")
         if (in_grid) then
            call save_annotation(attributes,grid%annotation)
         else if (in_xc) then
            call save_annotation(attributes,xp%annotation)
         else
            call die("unknown <annotation> element")
         endif
                  

end select

end subroutine begin_element
!----------------------------------------------------------------------

#ifndef PSML_USE_FOX
subroutine end_element(name)
#else
subroutine end_element(namespaceURI,localName,name)
character(len=*), intent(in)    :: namespaceURI
character(len=*), intent(in)    :: localName
#endif

character(len=*), intent(in)     :: name

integer :: i, nmajor, nminor

if (debug_parsing) print *, "-- end Element: ", trim(name)

select case(name)

      case ("radfunc")
         in_radfunc = .false.
         if (.not. associated(rp%data)) then
            call die("No data for radfunc!")
         endif

      case ("grid")
         in_grid = .false.
         !
         if (.not. got_explicit_grid_data) then
            call die("Need explicit grid data!")
         endif

      case ("data")
      !
      ! We are done filling up the radfunc data
      ! Check that we got the advertised number of items
      !
         in_data = .false.
         if (ndata /= size(rp%data)) then
            call die("npts mismatch in radfunc data")
         endif

      case ("grid-data")
      !
      ! We are done filling up the grid data
      ! Check that we got the advertised number of items
      !
         in_grid_data = .false.
         if (ndata_grid /= size(grid%grid_data)) then
            call die("npts mismatch in grid")
         endif
         if (debug_parsing) print *, "Got grid data: ", got_explicit_grid_data

      case ("pseudocore-charge")
         in_coreCharge = .false.

      case ("valence-charge")
         in_valenceCharge = .false.

      case ("semilocal-potentials")
         in_semilocal = .false.

         ! Generate indexes

         slp => pseudo%semilocal

         nmajor = 0
         nminor = 0
         do i = 1, slp%npots
            if (slp%set(i) == "major") then
               nmajor = nmajor + 1
               slp%major(nmajor) = i
            else if (slp%set(i) == "minor") then
               nminor = nminor + 1
               slp%minor(nminor) = i
            else
               call die("wrong set in potential")
            endif
         enddo

         if (nmajor /= slp%npots_major) then
            call die("wrong number of major potentials")
         endif
         if (nminor /= slp%npots_minor) then
            call die("wrong number of minor potentials")
         endif
               
      case ("pseudopotential-operator")
         in_psoperator = .false.

      case ("projectors")
         in_projectors = .false.

         ! Generate indexes

         pop => pseudo%psoperator

         nmajor = 0
         nminor = 0
         do i = 1, pop%nprojs
            if (pop%set(i) == "major") then
               nmajor = nmajor + 1
               pop%major(nmajor) = i
            else if (pop%set(i) == "minor") then
               nminor = nminor + 1
               pop%minor(nminor) = i
            else
               call die("wrong set in projector")
            endif
         enddo
         pop%nprojs_major = nmajor
         pop%nprojs_minor = nminor

      case ("vps")
         in_vps = .false.

      case ("proj")
         in_proj = .false.

      case ("local-potential")
         in_local_potential = .false.

      case ("pseudo-wave-functions")
         in_pseudowavefun = .false. 

         ! Generate indexes

         wfp => pseudo%pswfs

         nmajor = 0
         nminor = 0
         do i = 1, wfp%npswfs
            if (wfp%set(i) == "major") then
               nmajor = nmajor + 1
               wfp%major(nmajor) = i
            else if (wfp%set(i) == "minor") then
               nminor = nminor + 1
               wfp%minor(nminor) = i
            else
               call die("wrong set in pseudo-wave-function")
            endif
         enddo

         if (nmajor /= wfp%npswfs_major) then
            call die("wrong number of major pswfs")
         endif
         if (nminor /= wfp%npswfs_minor) then
            call die("wrong number of minor pswfs")
         endif

      case ("pswf")
         in_pswf = .false.

      case ("valence-configuration")
         in_config_val = .false.

      case ("exchange-correlation")
         in_xc = .false.

      case ("libxc-info")
         in_libxc_info = .false.

      case ("provenance")
         in_provenance = .false.

      case ("header")
         in_header = .false.

      case ("psml")
!         call dump_pseudo(pseudo)

end select

end subroutine end_element
!----------------------------------------------------------------------

subroutine pcdata_chunk(chunk)
#ifdef PSML_USE_FOX
  use fox_extra, only: build_data_array
#else
  use xmlf90_sax, only: build_data_array
#endif
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

!      print *, "Header data:"
!      print *, trim(chunk)

endif

end subroutine pcdata_chunk
!----------------------------------------------------------------------

     ! Annotations are encoded as an association list
     ! in a couple of arrays
     ! ( (key "value") (key "value") ...)
     ! 
subroutine save_annotation(atts,annotation)
  use assoc_list, ps_annotation_t => assoc_list_t
#ifdef PSML_USE_FOX
  use Fox_common, only: dictionary_t, len=>getLength
  use fox_extra, only: get_value=>get_value_by_index, get_key
#else
  use xmlf90_sax, only: dictionary_t, get_value, get_key, len
#endif
       type(dictionary_t), intent(in) :: atts
       type(ps_annotation_t), intent(out) :: annotation
       
       integer :: n, i, status
       character(len=50) :: key, value

       n = len(atts)
       call assoc_list_init(annotation,n,status)
       if (status /= 0) call die("Failed to init annotation object")
       do i = 1, n
          call get_key(atts,i,key,status)
          if (status /= 0) call die("cannot get key in atts dict")
          call get_value(atts,i,value,status)
          if (status /= 0) call die("cannot get value in atts dict")
          call assoc_list_insert(annotation,key,value,status)
          if (status /= 0) call die("Failed to insert annotation pair")
       enddo
     end subroutine save_annotation

end module m_psml_parsing_helpers












