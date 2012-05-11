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

      program fatbands

        use m_getopts
        use f2kcli

      implicit none

      integer           nk, nspin, nband, ik, is, ib

      double precision, allocatable ::  e(:,:,:), k(:), fat(:,:,:)
      double precision  ef, kmin, kmax, emin, emax


  integer  :: min_band = -huge(1)
  integer  :: max_band =  huge(1)
  logical  :: min_band_set = .false.
  logical  :: max_band_set = .false.
  logical  :: band_interval_set = .false.


  integer  :: bands_u = 1, fat_u = 2
  integer  :: n_opts, iostat, nargs, nlabels
  character(len=132) :: opt_name, opt_arg
  character(len=132) :: bandfile, fatfile
  logical  :: want_fat = .false.
  integer  :: ika, nband_fat, nk_fat, nspin_fat

  !
  !     Process options
  !
  n_opts = 0
  do
     call getopts('hF:b:B:',opt_name,opt_arg,n_opts,iostat)
     if (iostat /= 0) exit
     select case(opt_name)
     case ('F')
        read(opt_arg,*) fatfile
        want_fat = .true.
     case ('b')
        read(opt_arg,*) min_band
        min_band_set = .true.
     case ('B')
        read(opt_arg,*) max_band
        max_band_set = .true.
     case ('h')
      !        call manual()
     case ('?',':')
        write(0,*) "Invalid option: ", opt_arg(1:1)
        write(0,*) "Usage: fatbands  [ options ] bandsfile"
        write(0,*) "Use -h option for manual"
        STOP
     end select
  enddo

  nargs = command_argument_count()
  nlabels = nargs - n_opts + 1
  if (nlabels /= 1)  then
     write(0,*) "Usage: mprop [ -d ] [ -h ] MPROP_FILE_ROOT"
     write(0,*) "Use -h option for manual"
     STOP
  endif

  call get_command_argument(n_opts,value=bandfile,status=iostat)
  if (iostat /= 0) then
     STOP "Cannot get bands file"
  endif

  open(bands_u,file=trim(bandfile), status='old',action="read",position="rewind")

      read(bands_u,*) ef
      read(bands_u,*) kmin, kmax
      read(bands_u,*) emin, emax
      read(bands_u,*) nband, nspin, nk

      allocate(k(nk))
      allocate(e(nband,nspin,nk))

      write(6,"(2a)") '# GNUBANDS: Utility for SIESTA to transform ',  &
                     'bands output into Gnuplot format'
      write(6,"(a)") '#'
      write(6,"(2a)") '#                                           ',  &
                     '       Emilio Artacho, Feb. 1999'
      write(6,"(2a)") '# ------------------------------------------', &
                     '--------------------------------'
      write(6,"(a,f10.4)")  '# E_F               = ', ef
      write(6,"(a,2f10.4)") '# k_min, k_max      = ', kmin, kmax
      write(6,"(a,2f10.4)") '# E_min, E_max      = ', emin, emax
      write(6,"(a,3i6)")    '# Nbands, Nspin, Nk = ', nband, nspin, nk
      write(6,"(a)") '#'
      write(6,"(a)") '#        k            E'
      write(6,"(2a)") '# ------------------------------------------',   &
                     '--------------------------------'

      read(bands_u,*) (k(ik),((e(ib,is,ik),ib=1,nband), is=1,nspin), ik=1,nk)

      close(bands_u)

     if (min_band_set .and. (min_band < 1)) then
        print "(a)", " ** Min_band implicitly reset to 1..."
        min_band = 1
     endif
     if (min_band_set .and. (min_band > nband)) then
        print "(a,2i5)", " ** Min_band is too large for some k-points: (min_band, nband):", min_band, nband
        STOP
     endif
     if (max_band_set .and. (max_band > nband)) then
        print "(a,2i5)", " ** Max_band is too large for some k-points: (max_band, nband):", max_band, nband
        print "(a)", " ** Max_band will be effectively reset to its maximum allowed value"
        max_band = nband
     endif
     if (max_band_set .and. (max_band < max(1,min_band))) then
        print "(a,2i5)", " ** Max_band is less than the effective min_band: (max_band, eff min_band):", max_band, max(1,min_band)
        STOP
     endif

     if (want_fat) then
        open(fat_u,file=trim(fatfile), status='old',action="read",position="rewind")
        read(fat_u,*)
        read(fat_u,*) nband_fat, nspin_fat, nk_fat
        if (nk_fat /= nk) STOP "wrong number of k-points"
        if (nspin_fat /= nspin) STOP "wrong number of spins"
        allocate(fat(nband_fat,nspin,nk))
        do ik = 1, nk
           read(fat_u,*) ika, ((fat(ib,is,ik), ib = 1, nband_fat), is = 1, nspin)
        enddo
        close(fat_u)
     endif

      
      do is = 1, nspin
        do ib = min_band, max_band
           if (want_fat) then
              write(6,"(3f14.6)") ( k(ik), e(ib,is,ik), fat(ib,is,ik), ik = 1, nk)
           else
              write(6,"(3f14.6)") ( k(ik), e(ib,is,ik),1.0, ik = 1, nk)
           endif
           write(6,"(/)")
        enddo
      enddo

      if (want_fat) then
         deallocate(fat)
      endif

    end program fatbands
