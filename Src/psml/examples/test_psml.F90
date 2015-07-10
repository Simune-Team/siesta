program test_psml

!
! Example of extraction of information from PSML file
! This is quite a "low-level" example, in which the
! program queries all the contents and prints the information.
!
! Examples closer to the use cases in electronic-structure
! codes are in preparation
!
! No examples of extraction of "j" or "s" values yet.
!
use m_psml
use m_getopts

integer, parameter :: dp = selected_real_kind(10,100)
character(len=1), dimension(0:4) :: sym = (/ "s", "p", "d", "f", "g" /)

type(ps_t)   :: ps

      character(len=200) :: filename
      logical            :: debug, plot
      character(len=200) :: opt_arg, mflnm, ref_line
      character(len=10)  :: opt_name 
      character(len=10)  :: interp_type
      integer            :: interp_quality
      integer :: nargs, iostat, n_opts, nlabels, iorb, ikb
      
      integer :: i, j, l, n, num, nfun, set, seq
      real(dp) :: ekb, rc, jval

      real(dp), allocatable, dimension(:) :: r
      real(dp) :: Rmax, delta, val
      integer  :: npts, ir
      character(len=30) :: outname,  target_set, set_str
      integer, allocatable   :: setidx(:), idxl(:)

      real(dp), allocatable   :: raw_r(:), raw_core(:)

      ! There is a copy of dpnint1 inside the library
      ! This external declaration is to support calls
      ! to ps_SetInterpolator
      external :: interpolate_drh

      ! Uncomment the following when a new interpolator
      ! is added to the source tree. The one used previously
      ! was removed due to license issues.
!!      external :: interpolate_other  

!
!     Process options
!
      n_opts = 0
      debug = .false.
      plot = .false.

      Rmax = 4.0_dp   ! default maximum radius for plots
      npts = 200      ! default number of points

      interp_type = "drh"
      interp_quality = 7  ! this is npoint=2 in OTHER interpolator
      do
         call getopts('dR:n:pi:q:',opt_name,opt_arg,n_opts,iostat)
         if (iostat /= 0) exit
         select case(opt_name)
           case ('d')
              debug = .true.
           case ('p')
              plot = .true.
           case ('R')
              read(opt_arg,*) Rmax
           case ('n')
              read(opt_arg,*) npts
           case ('i')
              read(opt_arg,*) interp_type
           case ('q')
              read(opt_arg,*) interp_quality
           case ('?',':')
             write(0,*) "Invalid option: ", opt_arg(1:1)
             write(0,*) "Usage: test_psml [ -d -R Rmax -n npts ] FILE"
             write(0,*) " -d for debug output"
             write(0,*) " -R Rmax : maximum range of output grid (def: 4 bohr)"
             write(0,*) " -n npts : number of points of output grid (def: 200)"
             write(0,*) " -i type : interpolator type (drh) ('other' not implemented)"
             write(0,*) " -q nq   : interpolator quality (npoly for drh)"
             STOP
          end select
       enddo

       nargs = command_argument_count()
       nlabels = nargs - n_opts + 1
       if (nlabels /= 1)  then
          write(0,*) "Usage: test_psml [ -d -R Rmax -n npts ] FILE"
          write(0,*) " -d for debug output"
          write(0,*) " -R Rmax : maximum range of output grid (def: 4 bohr)"
          write(0,*) " -n npts : number of points of output grid (def: 200)"
          STOP
       endif

       call get_command_argument(n_opts,value=filename,status=iostat)
       if (iostat /= 0) then
          STOP "Cannot get filename"
       endif
!
call ps_destroy(ps)
print "(a)", "Processing: " // trim(filename)
call psml_reader(filename,ps,debug=debug)

#ifdef __NO_PROC_POINTERS__
if (trim(interp_type)=="other") then
   print "(a)", "No support for other interpolators if no proc pointers"
   STOP 
endif
call ps_SetInterpolatorQuality(interp_quality)
print "(a,i3)", "Using DRH (default) interpolator with npoly:", interp_quality

#else

if (trim(interp_type)=="other") then
!   call ps_SetInterpolator(interpolate_other,interp_quality)
!   print "(a,i3)", "Using OTHER interpolator with quality index:",interp_quality
   print "(a)", "Please implement 'other' interpolator first"
   STOP 
else if (trim(interp_type)=="drh") then
   call ps_SetInterpolator(interpolate_drh,interp_quality)
   print "(a,i3)", "Using DRH interpolator with npoly:",interp_quality
else
   STOP "unknown interpolator"
endif

#endif

! Set up our grid
allocate(r(npts))
delta = Rmax / (npts - 1)
do ir = 1, npts
   r(ir) = (ir-1)*delta
enddo

print *, "Relativity: ", trim(ps_Relativity(ps))

print *, "Exchange-correlation info:"
nfun = ps_NLibxcFunctionals(ps)
print *, "Number of functionals: ", nfun
do i = 1, nfun
  print "(a,1x,a,i4,f10.6)", "name, id, weight:", &
    trim(ps_LibxcName(ps,i)), ps_LibxcId(ps,i), ps_LibxcWeight(ps,i)
enddo
print *
!
npots = ps_Number_Of_Potentials(ps,SET_ALL)
print *, "There are ", npots, " semi-local pots"
!
setidx = ps_Potential_Indexes(ps,SET_ALL)
do i = 1, npots
   l = ps_Potential_L(ps,setidx(i))
   n = ps_Potential_N(ps,setidx(i))
   rc = ps_Potential_Rc(ps,setidx(i))
   set = ps_Potential_Set(ps,setidx(i))
   set_str = str_of_set(set)
   if (set == SET_LJ) then
      jval = ps_Potential_J(ps,setidx(i))
      print "(a,1x,i1,a1,1x,f3.1,f10.3)", trim(set_str), n, sym(l), jval, rc
   else
      print "(a,1x,i1,a1,f10.3)", trim(set_str), n, sym(l), rc
   endif

   if (plot) then
   write(outname,"(a,i0,a1)") "Vsl." //trim(set_str) // ".", n, sym(l)
   open(4,file=outname,form="formatted",status="unknown",&
        action="write",position="rewind")
   do ir = 1, npts
      val = ps_Potential_Value(ps,setidx(i),r(ir))
      write(4,*) r(ir), val
   enddo
   close(4)
   endif
enddo

if (.not. ps_HasPsOperator(ps)) then
   print *, "This file does not have a psoperator section"
else
   print *, "Vlocal type: " // trim(ps_LocalPotential_Type(ps))

   if (plot) then
   write(outname,"(a)") "Vlocal"
   open(4,file=outname,form="formatted",status="unknown",&
        action="write",position="rewind")
   do ir = 1, npts
      val = ps_LocalPotential_Value(ps,r(ir))
      write(4,*) r(ir), val
   enddo
   close(4)
   endif

nprojs = ps_Number_Of_Projectors(ps,SET_ALL)
print *, "There are ", nprojs, " projectors"

setidx = ps_Projector_Indexes(ps,SET_ALL)
print "(a20,a8,a4,a12)", "set", "l (j)", "seq", "Ekb"
do l = 0, 4
   idxl = ps_Projector_Indexes_byL(ps,l,setidx)
   num = size(idxl)
   do j = 1, num
      ekb = ps_Projector_Ekb(ps,idxl(j))
      set = ps_Projector_Set(ps,idxl(j))
      seq = ps_Projector_Seq(ps,idxl(j))
      if (set == SET_LJ) then
         jval = ps_Projector_J(ps,idxl(j))
         print "(a20,i4,1x,f3.1,i4,f12.6)", set_str, l, jval, seq, ekb
         write(outname,"(a,i0,a1,f3.1,a1,i0)") "Proj." //trim(set_str) // ".", &
                                                l, ".", jval,".", seq
      else
         print "(a20,i4,i4,f12.6)", set_str, l, seq, ekb
         write(outname,"(a,i0,a1,i0)") "Proj." //trim(set_str) // ".", l, ".", seq
      endif

      if (plot) then

      open(4,file=outname,form="formatted",status="unknown",&
           action="write",position="rewind")
      do ir = 1, npts
         val = ps_Projector_Value(ps,idxl(j),r(ir))
         write(4,*) r(ir), val
      enddo
      close(4)
   endif

   enddo
enddo
endif

! Valence charge density

   if (plot) then
      write(outname,"(a)") "Valence.charge"
      open(4,file=outname,form="formatted",status="unknown",&
           action="write",position="rewind")
      do ir = 1, npts
         val = ps_ValenceCharge_Value(ps,r(ir))
         write(4,*) r(ir), val
      enddo
      close(4)
   endif

! Pseudo-Core charge

   if (ps_HasCoreCorrections(ps)) then
    if (plot) then
      write(outname,"(a)") "Core.charge"
      open(4,file=outname,form="formatted",status="unknown",&
           action="write",position="rewind")
      do ir = 1, npts
         val = ps_CoreCharge_Value(ps,r(ir))
         write(4,"(1p,2e21.13)") r(ir), val
      enddo
      close(4)
    endif
    if (plot) then
      call ps_CoreCharge_GetRawData(ps,raw_r,raw_core)
      write(outname,"(a)") "Raw.Core.charge"
      open(4,file=outname,form="formatted",status="unknown",&
           action="write",position="rewind")
      do ir = 1, size(raw_r)
         write(4,*) raw_r(ir), raw_core(ir)
      enddo
      close(4)
    endif
   endif


end program test_psml

