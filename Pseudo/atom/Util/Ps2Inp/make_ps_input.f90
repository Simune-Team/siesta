program make_ps_input
!
! Reads a psf or vps file and attempts to produce an ATOM inp file from it
!
  use pseudopotential, only: pseudopotential_t
  use pseudopotential, only: pseudo_read, get_ps_conf

  implicit none

      integer, parameter :: dp = selected_real_kind(16,100)

  character(len=20) :: label
  type(pseudopotential_t) :: p

  integer   :: n_val_first, ncore, lmax, inp_lun, l, n

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

  real(dp) :: aa, bb, rmax

  read(5,*) label
  call pseudo_read(label,p)

  lmax = p%npotd-1
  allocate(orb_arr(0:lmax))
  allocate(zdown_arr(0:lmax))
  allocate(zup_arr(0:lmax))
  allocate(rc_arr(0:lmax))

  call get_ps_conf(p%irel,lmax,p%text,p%gen_zval, &
                   orb_arr,zdown_arr,zup_arr,rc_arr)

  read(orb_arr(0),"(i1)") n_val_first
  read(orb_arr(0),"(1x,a1)") code_val_first
  print *, "First valence orbital: ", orb_arr(0)
  print *, "Pieces: ", n_val_first, code_val_first
  if (code_val_first /= "s") call die("First valence pseudo is not s!")

  call core_below_orb(lmax,orb_arr,ncore)
  print *, "There are ", ncore, " core orbitals"

!
!  call io_assign(inp_lun)
  inp_lun = 2
  open(inp_lun, file=trim(label)//".inp", form="formatted",  &
      status="replace",action="write",position="rewind")

  if (p%nicore(1:2) == "pc" ) then
     job = "pe"
  else
     job = "pg"
  endif
  write(inp_lun,"(3x,a2,a50)") job , " -- file generated from " // trim(label) // " ps file"
  write(inp_lun,"(8x,a3)") flavor
  select case (p%irel)
     case ("isp") 
        ispp = "s"
     case ("rel")
        ispp = "r"
     case default
        ispp = " "
  end select
  write(inp_lun,"(3x,a2,3x,a2,a1)") p%name, p%icorr, ispp
  aa = 0.0
  bb = 0.0
  rmax = 0.0
  write(inp_lun,"(6f10.3)") 0.0, 0.0, 0.0, aa, bb, rmax
  write(inp_lun,"(2i5)") ncore, lmax+1

  do l = 0, lmax
     read(orb_arr(l),"(i1)") n
     write(inp_lun,"(2i5,2f10.3,2x,a5)") n, l, zdown_arr(l), zup_arr(l), "  #" // orb_arr(l)
  end do
  do l = 0, lmax
     write(inp_lun,"(f10.5)",advance="no") rc_arr(l)
  enddo
  do l = lmax + 1, 3
     write(inp_lun,"(f10.5)",advance="no") 0.0
  end do

  if (job == "pe") then
     write(inp_lun,"(2f10.5,a)") 0.0, -1.0, " FIX CORE RADIUS"
  else
     write(inp_lun,"(2f10.5,a)") 0.0, 0.0
  endif

!
! Final blank line and ruler
  write(inp_lun,"(/,a)") ruler

!  call io_close(inp_lun)
  close(inp_lun)

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
  
end program make_ps_input
