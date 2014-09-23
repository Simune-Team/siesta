module m_ncps_parsing_helpers
!
!  This module reads a pseudopotential file written in XML (PSML format)
!  A full example of the building up of a data structure using
!  the SAX paradigm.
!
#ifdef XMLF90
 use flib_sax
#else
!! use FoX_common, only: get_value=>getValue
 use FoX_sax
#endif

 use m_ncps_xml_ps_t        ! Data types

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
!    to keep several instances alive at the same time.
! -- If "pseudo" here is a normal variable, it should also
!    be "cleaned" before the next use. The current usage
!    in Abinit falls in this category: psxml is a pointer
!    associated to "pseudo", and cleaned after use.
!
type(ps_t), pointer, public, save :: pseudo => null()

!**AG**
! Make sure that this is provided by the user
private :: die

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
type(valence_charge_t), private, pointer  :: valp => null()
type(core_charge_t), private, pointer     :: corep => null()
type(radfunc_t), private, pointer         :: rp => null()

#ifndef XMLF90
!
! Helper routines taken from xmlf90

public :: build_data_array

interface build_data_array
      module procedure build_data_array_real_sp,  &
                       build_data_array_real_dp,  &
                       build_data_array_integer
end interface
#endif

CONTAINS  !===========================================================

!----------------------------------------------------------------------
#ifdef XMLF90
subroutine begin_element(name,attributes)
#else
subroutine begin_element(namespaceURI,localName,name,attributes)
character(len=*), intent(in)    :: namespaceURI
character(len=*), intent(in)    :: localName
#endif
character(len=*), intent(in)    :: name
type(dictionary_t), intent(in)  :: attributes

character(len=100)  :: value
integer             :: status

integer             :: i

!print *, "Element: ", trim(name)

select case(name)

      case ("psml")

         ! Allocate unconditionally, as
         ! we must avoid re-using a previous
         ! version (say, when dealing with
         ! two or more elements)

         allocate(pseudo)


         call get_value(attributes,"version",value,status)
         if (value /= "0.7") then
            call die("Can only work with PSML version 0.7 files")
         endif

      case ("provenance")
         in_provenance = .true.
         pp => pseudo%provenance

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

!!         ! This attribute is optional
!!         call get_value(attributes,"annotation",grid%annotation,status)
!!         if (status /= 0 ) grid%annotation=""

         !
         ! In this way we allow for a private grid for each radfunc,
         ! or for a global grid specification
         !
         if (in_radfunc) then
            if (associated(rp%grid)) then
               call die("psml: Two grids specified for a radfunc")
            endif
            rp%grid => grid
         else
            ! We should really check that we are at the top level,
            ! and not, say, at the semilocal or pswf level (although
            ! it could be useful to allow these "regional" grids)

            if (associated(pseudo%global_grid)) then
               call die("psml: Two global grids specified")
            endif
            !print *, "Found global grid"
            pseudo%global_grid => grid
         endif

      case ("data")
         if (.not. in_radfunc) then
            call die("<data> element outside <rad_func> element")
         endif
         in_data = .true.
         if (.not. associated(rp%grid)) then
            if (associated(pseudo%global_grid)) then
               rp%grid => pseudo%global_grid
            else
               call die("Cannot find grid data for radfunc")
            endif
         endif
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

#ifdef XMLF90
subroutine end_element(name)
#else
subroutine end_element(namespaceURI,localName,name)
character(len=*), intent(in)    :: namespaceURI
character(len=*), intent(in)    :: localName
#endif

character(len=*), intent(in)     :: name

integer :: i, nmajor, nminor

!print *, "-- end Element: ", trim(name)

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
!         print *, "Got grid data: ", got_explicit_grid_data

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
               
      case ("vps")
         in_vps = .false.

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

      ! To be bypassed by the global "die" of a client program
      ! Use interface mechanism as in newer versions of Siesta
      !
      subroutine die(str)
      character(len=*), intent(in), optional   :: str
      if (present(str)) then
         write(unit=0,fmt="(a)") trim(str)
      endif
      write(unit=0,fmt="(a)") "Stopping Program"
      stop
      end subroutine die

#ifndef XMLF90
      subroutine get_value(dict,attr,value,status)
        use FoX_common, only: getValue
        type(dictionary_t), intent(in) :: dict
        character(len=*), intent(in)   :: attr
        character(len=*), intent(out)  :: value
        integer, intent(out)           :: status

        status = 0
        value = getValue(dict,attr)
        if (value == "") status = 1
      end subroutine get_value

! (Copy of "m_converters.f90" in the xmlf90 distribution...)
!
! Takes a string and turns it into useful data structures,
! such as numerical arrays.
!
! NOTE: The string must contain *homogeneous* data, i.e.: all real numbers,
! all integers, etc.
!

!---------------------------------------------------------------
subroutine build_data_array_real_dp(str,x,n)
integer, parameter  :: dp = selected_real_kind(14)
!
character(len=*), intent(in)                ::  str
real(kind=dp), dimension(:), intent(inout)  ::    x
integer, intent(inout)                      ::    n

integer                            :: ntokens, status, last_pos
character(len=len(str))  :: s

s = str
call token_analysis(s,ntokens,last_pos)
!if (debug) print *, "ntokens, last_pos ", ntokens, last_pos
!if (debug) print *, s
if ((n + ntokens) > size(x)) STOP "data array full"
read(unit=s(1:last_pos),fmt=*,iostat=status) x(n+1:n+ntokens)
if (status /= 0) STOP "real conversion error"
n = n + ntokens

end subroutine build_data_array_real_dp
!---------------------------------------------------------------

subroutine build_data_array_real_sp(str,x,n)
integer, parameter  :: sp = selected_real_kind(6)
!
character(len=*), intent(in)                :: str
real(kind=sp), dimension(:), intent(inout)  ::    x
integer, intent(inout)                      ::    n

integer                            :: ntokens, status, last_pos
character(len=len(str))  :: s

s = str
call token_analysis(s,ntokens,last_pos)
!if (debug) print *, "ntokens, last_pos ", ntokens, last_pos
!if (debug) print *, s
if ((n + ntokens) > size(x)) STOP "data array full"
read(unit=s(1:last_pos),fmt=*,iostat=status) x(n+1:n+ntokens)
if (status /= 0) STOP "real conversion error"
n = n + ntokens

end subroutine build_data_array_real_sp

!---------------------------------------------------------------
subroutine build_data_array_integer(str,x,n)
integer, parameter  :: sp = selected_real_kind(14)
!
character(len=*), intent(in)                :: str
integer, dimension(:), intent(inout)        ::    x
integer, intent(inout)                      ::    n

integer                            :: ntokens, status, last_pos
character(len=len(str))  :: s

s = str
call token_analysis(s,ntokens,last_pos)
!if (debug) print *, "ntokens, last_pos ", ntokens, last_pos
!if (debug) print *, s
if ((n + ntokens) > size(x)) STOP "data array full"
read(unit=s(1:last_pos),fmt=*,iostat=status) x(n+1:n+ntokens)
if (status /= 0) STOP "integer conversion error"
n = n + ntokens

end subroutine build_data_array_integer


!==================================================================

function is_separator(c) result(sep)
character(len=1), intent(in)          :: c
logical                               :: sep

 sep = ((c == char(32)) .or. (c == char(10))             &
         .or. (c == char(9)) .or. (c == char(13)))

end function is_separator
!----------------------------------------------------------------
function is_CR_or_LF(c) result(res)
character(len=1), intent(in)          :: c
logical                               :: res

 res = ((c == char(10)) .or. (c == char(13)))

end function is_CR_or_LF

!==================================================================

subroutine token_analysis(str,ntokens,last_pos)
!
character(len=*), intent(inout)          :: str
integer, intent(out)                     :: ntokens, last_pos
!
!
! Checks the contents of a string and finds the number of tokens it contains
! The standard separator is generalized whitespace (space, tab, CR, or LF)
! It also returns the last useful position in the string (excluding
! separator characters which are not blanks, and thus not caught by the
! (len_)trim fortran intrinsic). This is necessary to perform list-directed
! I/O in the string as an internal file.
! 
! Also, replace on the fly CR and LF by blanks. This is necessary if
! str spans more than one record. In that case, internal reads only 
! look at the first record. 
! -- ** Compiler limits on size of internal record??
!
integer           :: i, str_length
logical           :: in_token
character(len=1)  :: c

in_token = .false.
ntokens = 0
last_pos = 0

str_length = len_trim(str)
!print *, "string length: ", str_length

do i = 1, str_length
      c = str(i:i)

      if (in_token) then
         if (is_separator(c)) then
            in_token = .false.
            if (is_CR_or_LF(c)) str(i:i) = " "
         else
            last_pos = i
         endif

      else   ! not in token
         
         if (is_separator(c)) then
            if (is_CR_or_LF(c)) str(i:i) = " "
            ! do nothing
         else
            in_token = .true.
            last_pos = i
            ntokens = ntokens + 1
         endif
      endif
enddo
!print *, "ntokens, last_pos: ", ntokens, last_pos

end subroutine token_analysis
#endif

     ! Annotations are encoded as an association list
     ! in a couple of arrays
     ! ( (key "value") (key "value") ...)
     ! 
     subroutine save_annotation(atts,annotation)
       use assoc_list, ps_annotation_t => assoc_list_t

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

end module m_ncps_parsing_helpers












