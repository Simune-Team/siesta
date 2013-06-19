! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!

program Eig2DOS

  use m_getopts
  use f2kcli

  implicit none

! Utility for obtaining the total density of states out of the .EIG 
! file generated by SIESTA. Energies in eV.
!
! Written by Emilio Artacho, April 1999.
! Rewritten by Alberto Garcia, April 2012.

! Read from file "systemlabel.EIG" generated by SIESTA

!  -s     - peak width (in eV) for broadening (gaussian or lorentzian).
!  -n     - number of points in the energy window
!  -m, -M - energy window: Emin and Emax 
!  -f     - shift E_F to zero
!
! Density of states in standard output, for both spins.
! If nspin = 1, the first spin component is multiplied by 2.
!
!   NOT printed in this version:
!   the number of electrons (states) in the energy window, 
!   with and without broadening (number of eigenvalues and 
!   DOS integral, with simple sum)
!
! Equal weight for all k-points is in principle assumed. This condition can be 
! removed by using the "-k file" option and reading the weights from "systemlabel.KP"

  integer, parameter :: dp = selected_real_kind(10,100)

  character(len=200) :: opt_arg
  character(len=10)  :: opt_name 
  integer :: nargs, iostat, n_opts, nlabels

  integer    nk, nspin, nband, ne, ie, ik, ika, is, ib, nel
  integer    nk_kpoints, ik_read

  real(dp)   e, eincr, ef, pi, x, sta, norm
  real(dp)   emin_file, emax_file
  real(dp)   integral(2), weight, k(3)

  real(dp), allocatable    :: eig(:,:), dos(:,:)

  logical ::  debug    = .false.
  character(len=256) :: eig_file, kpoint_file

  integer  :: npts_energy = 200
  real(dp) :: emin        = huge(1.0_dp)
  real(dp) :: emax        = -huge(1.0_dp)
  real(dp) :: smear       = 0.2_dp
  logical  :: loren       = .false.
  logical  :: emin_given  = .false.
  logical  :: emax_given  = .false.
  logical  :: using_weights = .false.
  logical  :: energies_only = .false.
  logical  :: shift_efermi = .false.
  integer  :: min_band = 0
  integer  :: max_band = 0
!
!     Process options
!
  n_opts = 0
  do
     call getopts('dhefls:n:m:M:b:B:k:',opt_name,opt_arg,n_opts,iostat)
     if (iostat /= 0) exit
     select case(opt_name)
     case ('d')
        debug = .true.
     case ('e')
        energies_only = .true.
     case ('f')
        shift_efermi = .true.
     case ('s')
        read(opt_arg,*) smear
     case ('n')
        read(opt_arg,*) npts_energy
     case ('m')
        emin_given = .true.
        read(opt_arg,*) emin
     case ('M')
        emax_given = .true.
        read(opt_arg,*) emax
     case ('b')
        read(opt_arg,*) min_band
     case ('B')
        read(opt_arg,*) max_band
     case ('l')
        loren = .true.
     case ('k')
        using_weights = .true.
        read(opt_arg,*) kpoint_file
     case ('h')
        call manual()
        STOP
     case ('?',':')
        write(0,*) "Invalid option: ", opt_arg(1:1)
        write(0,*) "Usage: Eig2DOS [-d] [-s] [-l] [-m] [-M] [-b] [-B] [-h] EIGfile"
        write(0,*) "Use -h option for manual"
        STOP
     end select
  enddo

  nargs = command_argument_count()
  nlabels = nargs - n_opts + 1
  if (nlabels /= 1)  then
     write(0,*) "Usage: Eig2DOS [-d] [-s] [-l] [-m] [-M] [-b] [-B] [-h] EIGfile"
     write(0,*) "Use -h option for manual"
     STOP
  endif

  call get_command_argument(n_opts,value=eig_file,status=iostat)
  if (iostat /= 0) then
     STOP "Cannot get EIGfile"
  endif

!==================================================

! reading and initializing --------------------------------------------

  pi = acos(-1.0d0)

  write(*,"(2a)") '# EIG2DOS: Utility for SIESTA to obtain the ',  &
       'electronic density of states'
  write(*,"(2a)") '#  E. Artacho, Apr 1999, A. Garcia, Apr 2012'
  write(*,"(2a)") '# ------------------------------------------'

  open(unit=1,file=trim(eig_file),form="formatted",status="old", &
       action="read")
  read(1,*) ef
  read(1,*) nband, nspin, nk
  write(*,"(a)") "# Eigenvalues read from " // trim(eig_file)
  if (debug) print *, "Ef, nband, nspin, nk:", Ef, nband, nspin, nk

  if (using_weights) then
     open(unit=2,file=trim(kpoint_file),form="formatted",status="old", &
          action="read")
     read(2,*) nk_kpoints
     if (nk_kpoints /= nk) then
        write(*,*) "****: nk is different in EIG and KP files:", nk, nk_kpoints
        STOP
     endif
     write(*,"(a)") "# Kpoint weights read from " // trim(kpoint_file)
  endif

  allocate(eig(nband,nspin))
  allocate(dos(npts_energy,nspin))

!     Sanity checks

  if (min_band == 0) min_band = 1
  if (max_band == 0) max_band = nband

  if (min_band > max_band) STOP "min_band > max_band"
  if (max_band > nband) then
     write(0,*) "max_band reset to maximum in file: ", nband
  endif
  if (min_band < 1) then
     write(0,*) "min_band reset to 1. "
  endif
  if (min_band > nband) then
     write(0,*) "min_band reset to maximum in file: ", nband
  endif


  write(*,"(a,f7.3)") "# Using smearing parameter: ", smear
  write(*,"(a,i6,a)") "# Using ", npts_energy, " points in energy range"
  write(*,"(2(a,i0))") "# Selected bands: ", min_band ," to: ", max_band

  emin_file = huge(1.0_dp)
  emax_file = -huge(1.0_dp)
  do ik = 1, nk
     read(1,*) ika, ((eig(ib,is), ib = 1, nband), is = 1, nspin)
     emin_file = min(emin_file,minval(eig(min_band:max_band,1:nspin)))
     emax_file = max(emax_file,maxval(eig(min_band:max_band,1:nspin)))
  enddo

  write(*,"(a,2f15.7)") "# Emin, emax in file for selected band(s):", emin_file, emax_file
  if (energies_only) then
     deallocate(eig,dos)
     close(1)
     if (using_weights) then
        close(2)
     end if
     stop
  end if

! rewind and place file handle at the right point
  rewind(1)
  read(1,*)
  read(1,*)

  if (.not. emin_given) emin = emin_file - 6._dp*smear
  if (.not. emax_given) emax = emax_file + 6._dp*smear

  if (npts_energy .lt. 2) npts_energy = 2
  eincr = (emax - emin) / dfloat(npts_energy - 1)

  nel = 0
  dos(:,:) = 0.0_dp

! For each eigenvalue a smearing is applied

  do ik = 1, nk

     read(1,*) ika, ((eig(ib,is), ib = 1, nband), is = 1, nspin)
     if (using_weights) then
        read(2,*) ik_read, k(:), weight
        if (ik_read /= ik) STOP "ik mismatch"
     else
        weight = 1.0_dp / nk
     endif

     do is = 1, nspin
        do ib = min_band, max_band
           e = eig(ib,is) 
           if (shift_efermi)  e = e - ef
           if ( (e.ge.emin) .and. (e.le.emax) ) nel = nel + 1
           do ie = 1, npts_energy
              x = emin + (ie-1)*eincr - e
              if (loren) then
                 dos(ie,is) = dos(ie,is) + weight * smear / (smear*smear + x*x) 
              else
                 dos(ie,is) = dos(ie,is) + weight * exp( - x*x/(smear*smear) )
              endif
           enddo
        enddo
     enddo
  enddo

  if (loren) then
     norm = pi 
  else
     norm = sqrt(pi) * smear
  endif

  do is = 1, nspin
     do ie = 1, npts_energy
        dos(ie,is) = dos(ie,is)/norm
     enddo
  end do

! integral, extremely sophisticated -----------------------------------

  do is = 1, nspin
     integral(is) = 0.0_dp
     do ie = 1, npts_energy
        integral(is) = integral(is) + dos(ie,is)*eincr
     enddo
  enddo

! number of electrons -------------------------------------------------

  if (nspin .eq. 1) then
     nel = nel * 2
  endif
  sta = nel / dfloat(nk)

! output, prepared for gnuplot ----------------------------------------

  write(*,"(a,3i6)") '# Nbands, Nspin, Nk   = ', nband, nspin, nk
  if (shift_efermi) then
     write(*,"(a,f10.4,a)") '# E_F                 = ', ef , ' eV --> (shifted to ZERO)'
  else
     write(*,"(a,f10.4,a)") '# E_F                 = ', ef , ' eV (NOT shifted)'
  endif
  write(*,"(a,f10.4,a)") '# Broadening          = ', smear, ' eV'
!!      write(*,"(a,2f10.4,a)") '# Number of electrons = ', sta, sum(integral(:)),  &
!!                              ' per cell'
  write(*,"(a)") '#'

  if (nspin == 1) then
     write(*,"(a)") '#        E            N(up)  (=)   N(down)         Ntot'
     do ie = 1, npts_energy
        write(*,"(4f14.6)") emin + (ie-1)*eincr, dos(ie,1), dos(ie,1), 2._dp*dos(ie,1)
     end do
  else
     write(*,"(a)") '#        E            N(up)        N(down)         Ntot'
     do ie = 1, npts_energy
        write(*,"(4f14.6)") emin + (ie-1)*eincr, (dos(ie,is),is=1,nspin), sum(dos(ie,1:2))
     end do
  endif

  deallocate(eig,dos)
  close(1)
  if (using_weights) then
     close(2)
  end if

contains
  subroutine manual()
    write(0,"(a)") " -------------------"
    write(0,"(a)") " Usage: Eig2DOS [options] eigfile"
    write(0,"(a)") "  "
    write(0,"(a)") "        eigfile : SIESTA .EIG file"
    write(0,"(a)") " "
    write(0,"(a)") " OPTIONS: "
    write(0,"(a)") " "
    write(0,"(a)") " -h         Print this help"
    write(0,"(a)") " -d         Print debugging info"
    write(0,"(a)") " -e         Stop after printing Emin,Emax in file for selected bands"
    write(0,"(a)") " -f         Shift energy axis so that Efermi is at 0"
    write(0,"(a)") " -l         Use Lorentzian instead of Gaussian broadening"
    write(0,"(a)") " -s arg     Broadening parameter in eV"
    write(0,"(a)") " -n arg     Number of energy points at which to compute the DOS"
    write(0,"(a)") " -m emin    Minimum energy in range"
    write(0,"(a)") " -M emax    Maximum energy in range"
    write(0,"(a)") " -b arg     Index of first band to consider"
    write(0,"(a)") " -B arg     Index of last band to consider"
    write(0,"(a)") " -k kfile   Use kfile for k-point weight information (.KP format)"
    write(0,"(a)") " -------------------"

  end subroutine manual
end program Eig2DOS
