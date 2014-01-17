program semicore_info
!
! Reads a psf or vps file and attempts to produce an ATOM inp file from it
!
  use pseudopotential, only: pseudopotential_t
  use pseudopotential, only: pseudo_read, get_ps_conf
  use m_ground_state
  use periodic_table, only: atnumber

  implicit none

      integer, parameter :: dp = selected_real_kind(16,100)

  character(len=20) :: label
  type(pseudopotential_t) :: p

  type(ground_state_t)    :: gs

  integer   :: n_val_first, ncore, lmax, inp_lun, l, n, z

  character(len=1) :: code_val_first
  character(len=1) :: ispp
  character(len=2) :: job
  character(len=3) :: flavor = "tm2"   ! for now

  character(len=*), parameter ::   &
       ruler = "#23456789012345678901234567890123456789012345678901234567890      Ruler"

  character(len=2), allocatable :: orb_arr(:)
  real(dp), allocatable         :: zdown_arr(:)
  real(dp), allocatable         :: zup_arr(:)
  real(dp), allocatable         :: rc_arr(:)
  integer,  allocatable         :: gen_n(:)

  real(dp) :: aa, bb, rmax

  read(5,*) label
  call pseudo_read(label,p)

  z = atnumber(p%name)
  call ground_state(z,gs)

  lmax = p%npotd-1
  allocate(orb_arr(0:lmax))
  allocate(zdown_arr(0:lmax))
  allocate(zup_arr(0:lmax))
  allocate(rc_arr(0:lmax))
  allocate(gen_n(0:lmax))

  call get_ps_conf(p%irel,lmax,p%text,p%gen_zval, &
                   orb_arr,zdown_arr,zup_arr,rc_arr)

  print *, " -- Semicore analysis"
  do l = 0, lmax
     read(orb_arr(l),"(i1)") gen_n(l)
     print "(a,3i3)", "l, gen_n, val_n:", l, gen_n(l), gs%n(l)
     if (gen_n(l) < gs%n(l)) then
        print "(a,i2,1x,a)", "* There are ", gs%n(l) - gen_n(l), " semicore shells"
     endif
  enddo

  call core_below_orb(lmax,orb_arr,ncore)
  print *, "There are ", ncore, " core orbitals"


CONTAINS

subroutine core_below_orb(lmax,orb_arr,ncore)
integer, intent(in)  :: lmax
character(len=2), intent(in)  :: orb_arr(0:lmax)
integer, intent(out) :: ncore

integer :: l, i
character(len=1) :: symbol
logical  :: found


!  data for orbitals:                                                                                
!                                                                                                    
!              1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p,4f,5d,6s,6p,                                         
!                     7s, 5f, 6d, 7p                                                                 
!                                                                                                    
integer, parameter :: nshells = 19
integer, dimension(nshells) ::  nc = (/ 1, 2, 2, 3, 3, 3, 4, 4, 4,  &
                                        5, 5, 4, 5, 6, 6, 7, 5, 6, 7 /)
integer, dimension(nshells) ::  lc = (/ 0, 0, 1, 0, 1, 2, 0, 1, 2,  &
                                        0, 1, 3, 2, 0, 1, 0, 3, 2, 1 /)
character(len=1), dimension(4) :: sym = (/ "s", "p", "d", "f" /)


ncore = 0
do i = 1, nshells

! Search for orbital among the valence ones
  found = .false.
  do l = 0, lmax
     read(orb_arr(l),"(i1)") n
     read(orb_arr(l),"(1x,a1)") symbol
     if ((nc(i) == n .AND. sym(lc(i)+1) == symbol)) then
        found = .true.
        exit
     endif
  enddo

  if (found) then
     exit
  else
!     print *, "-- ", nc(i), lc(i)
     ncore = ncore + 1
  endif

enddo

if (ncore == nshells) call die("Major snafu in ncore_below_orb")
!print *, "Number of core orbs: ", ncore

end subroutine core_below_orb

end program semicore_info
