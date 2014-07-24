program fhi2psml

  ! Reads fort.X files from FHIPP98 and generates
  ! a PSML file

  ! Alberto Garcia, May 23, 2014
  ! Alberto Garcia, July 24, 2014

  use flib_wxml

  implicit none

  type(xmlf_t) :: xf

  integer, parameter :: dp = selected_real_kind(10,100)

  integer   :: i, j, npts, i_dummy, stat
  real(dp)  :: rho, rhop, rhopp, dummy
  real(dp)  :: znuc, zion, rcore, v
  real(dp)  :: total_valence_charge
  real(dp), allocatable :: r(:), chval(:), chcore(:)
  real(dp), allocatable :: vps(:,:), pswfs(:,:)
  real(dp), allocatable :: vps_raw(:)

  character(len=100)    :: line, fname
  character(len=100)    :: msg
  character*4      :: polattrib, relattrib, coreattrib
  character*1      :: pscode, char_dummy
  character(len=2) :: nameat
  character(len=40):: psflavor

  character(len=1), dimension(0:4) :: lsymb = (/'s','p','d','f','g'/)

  character*30 xcfuntype, xcfunparam
  integer          :: ncore, nval, ncp, norbs, lmax, npots
  integer          :: xc_code
  integer, allocatable  :: n(:), l(:)
  integer, allocatable  :: nn(:), ll(:)
  real(dp), allocatable :: f(:), ff(:), rc(:)
  real(dp), allocatable :: fdown(:), fup(:)

  logical :: tdopsp, nonrel, polarized, there_is_core, found
  !
  !     Read grid info from valence-charge file 
  !
  open(unit=1,file="fort.25",form="formatted")
  npts = 0
  do
     read(1,fmt=*,iostat=stat)
     if (stat /= 0) exit
     npts = npts + 1
  enddo
  close(1)
  !
  !     Valence charge
  !
  allocate(r(npts),chval(npts))
  open(unit=20,file="fort.25",form="formatted")
  do i = 1, npts
     read(20,*) r(i), rho, rhop, rhopp
     chval(i) = rho
  enddo


  call xml_OpenFile("FHIPSML",xf, indent=.false.)

  call xml_AddXMLDeclaration(xf,"UTF-8")

  call xml_NewElement(xf,"psml")
  call my_add_attribute(xf,"version","0.7")
  call my_add_attribute(xf,"energy_unit","hartree")
  call my_add_attribute(xf,"length_unit","bohr")

  call xml_NewElement(xf,"provenance")
  call my_add_attribute(xf,"creator","FHIPP98-2003")
  call my_add_attribute(xf,"translator","fhi2psml v0.2")
  call my_add_attribute(xf,"date","01-01-01")
  call xml_NewElement(xf,"input-file")
  call my_add_attribute(xf,"name","fort.20")
  open(1,file="fort.20",form="formatted",status="old", &
       position="rewind",action="read")
  do
     read(1,fmt="(a)",iostat=stat) line
     if (stat .ne. 0) exit
     call xml_AddPcData(xf,trim(line),line_feed=.true.)
  enddo
  close(1)
  call xml_EndElement(xf,"input-file")
  !
  call xml_NewElement(xf,"input-file")
  call my_add_attribute(xf,"name","fort.22")
  open(1,file="fort.22",form="formatted",status="old", &
       position="rewind",action="read")
  do
     read(1,fmt="(a)",iostat=stat) line
     if (stat .ne. 0) exit
     call xml_AddPcData(xf,trim(line),line_feed=.true.)
  enddo
  close(1)
  call xml_EndElement(xf,"input-file")
  call xml_EndElement(xf,"provenance")

  open(1,file="fort.20",form="formatted",status="old", &
       position="rewind",action="read")
  read(1,*) tdopsp, nonrel
  read(1,*) polarized
  close(1)

  print *, "Done reading fort.20"
  print *, "Polarized: ", polarized

  open(1,file="fort.22",form="formatted",status="old", &
       position="rewind",action="read")
  read(1,*) znuc, ncore, nval, xc_code, rcore

  nameat = symbol(nint(znuc))

  norbs = ncore + nval
  allocate (n(norbs), l(norbs), f(norbs))
  if (polarized) then
     allocate (fup(norbs),fdown(norbs))
  endif

  total_valence_charge = 0.0_dp
  ncp = ncore + 1
  do i = 1, norbs
     if (polarized) then
        read(1,*) n(i), l(i), fup(i), fdown(i)
        f(i) = fup(i) + fdown(i)
     else
        read(1,*) n(i), l(i), f(i)
     endif
     if (i > ncore) then
        total_valence_charge =   total_valence_charge + f(i)
     endif
  enddo
  read(1,*) lmax, pscode
  close(1)

  print *, "Done reading fort.20"

  zion = 0.0_dp
  do i = ncp, norbs
     zion = zion + f(i)
  enddo

  call xml_NewElement(xf,"valence-configuration")
  call my_add_attribute(xf,"total-valence-charge", str(total_valence_charge))
  do i = ncp, norbs
     if (f(i) .lt. 1.0e-10_dp) cycle
     call xml_NewElement(xf,"shell")
     call my_add_attribute(xf,"n",str(n(i)))
     call my_add_attribute(xf,"l",lsymb(l(i)))
     call my_add_attribute(xf,"occupation",str(f(i)))
     if (polarized) then
        call my_add_attribute(xf,"occupation-down",str(fdown(i)))
        call my_add_attribute(xf,"occupation-up",str(fup(i)))
     endif
     call xml_EndElement(xf,"shell")
  enddo
  call xml_EndElement(xf,"valence-configuration")
          
  npots = lmax + 1
  allocate (rc(npots), ll(npots), nn(npots), ff(npots))
  do i = 1, npots
     ll(i) = i - 1
     found = .false.
     ! look for the appropriate shell in the valence
     do j = ncp, norbs
        if (l(j) == ll(i)) then
           found = .true.
           nn(i) = n(j)
           ff(i) = f(j)
           exit
        endif
     enddo
     if (.not. found) then
        ! generate the appropriate effective n
        nn(i) = ll(i) + 1
        do j = 1, ncore
           if (l(j) == ll(i)) then
              nn(i) = nn(i) + 1
           endif
           ff(i) = 0.0_dp
        enddo
     endif
  enddo


  select case (pscode)
  case ("h","H")
     psflavor ="Hamann GNCPP"
  case ("t","T")
     psflavor ="Troullier-Martins"
  end select
  !
  !
  if (nonrel) then
     relattrib = "no"
  else
     relattrib = "yes"
  endif
  if (polarized) then
     polattrib = "yes"
  else
     polattrib = "no"
  endif

  inquire(file="fort.27",exist=there_is_core)
  if (there_is_core) then
     coreattrib = "yes"
  else
     coreattrib = "no"
  endif

  select case(xc_code)

  case(1) 
     xcfuntype    = 'LDA'
     xcfunparam   = 'Wigner'

  case(2) 
     xcfuntype    = 'LDA'
     xcfunparam   = 'Hedin-Lundqvist'

  case(3) 
     xcfuntype    = 'LDA'
     xcfunparam   = 'Ceperley-Alder PZ'

  case(4) 
     xcfuntype    = 'GGA'
     xcfunparam   = 'Perdew-Wang 91'

  case(5) 
     xcfuntype    = 'GGA'
     xcfunparam   = 'Becke X, Perdew C'

  case(6) 
     xcfuntype    = 'GGA'
     xcfunparam   = 'Perdew-Burke-Ernzerhof'

  case(7) 
     xcfuntype    = 'LDA'
     xcfunparam   = 'Zhao/Parr'

  case(8) 
     xcfuntype    = 'LDA'
     xcfunparam   = 'Ceperley-Alder PW'

  case(9) 
     xcfuntype    = 'GGA'
     xcfunparam   = 'Becke-Lee-Yang-Parr'

  case(10) 
     xcfuntype    = 'GGA'
     xcfunparam   = 'Perdew-Wang X Lee-Yang-Parr'

  case(11) 
     xcfuntype    = 'LDA'
     xcfunparam   = 'Exchange only'

  case(12,13) 
     xcfuntype    = 'KLI'
     xcfunparam   = '-----'

  case(14) 
     xcfuntype    = 'GGA'
     xcfunparam   = 'RPBE - Hammer et al'

  case(15) 
     xcfuntype    = 'GGA'
     xcfunparam   = 'revPBE Zhang+Yang'

  case(16) 
     xcfuntype    = 'MGGA'
     xcfunparam   = 'Perdew/Kurth/Zupan/Blaha'

  case(17) 
     xcfuntype    = 'GGA'
     xcfunparam   = 'PBE X, LDA C'

  case(18,19) 
     xcfuntype    = 'KLI'
     xcfunparam   = '-----'

  end select

  !
  !
  call xml_NewElement(xf,"header")
  call my_add_attribute(xf,"atomic-label",nameat)
  call my_add_attribute(xf,"atomic-number",str(znuc))
  call my_add_attribute(xf,"z-pseudo",str(zion))
  call my_add_attribute(xf,"flavor",psflavor)
  call my_add_attribute(xf,"relativistic",relattrib)
  call my_add_attribute(xf,"polarized",polattrib)
  call my_add_attribute(xf,"core-corrections",coreattrib)
  call my_add_attribute(xf,"xc-functional",trim(xcfuntype)//trim(xcfunparam))
  call my_add_attribute(xf,"xc-libxc-exchange","TBA")
  call my_add_attribute(xf,"xc-libxc-correlation","TBA")
  call xml_EndElement(xf,"header")

  call xml_NewElement(xf,"grid")
  call my_add_attribute(xf,"npts",str(npts))
  write(msg,"(3(a,':',a,'--'))")  "type","fhi",          &
                                 "r_1",trim(str(r(1))),       &
                                 "factor",trim(str(r(2)/r(1)))
  call my_add_attribute(xf,"annotation",trim(msg))
  call xml_NewElement(xf,"grid-data")
  call xml_AddArray(xf,r(1:npts))
  call xml_EndElement(xf,"grid-data")
  call xml_EndElement(xf,"grid")

  allocate (vps(npts,npots), pswfs(npts,npots))
  allocate (vps_raw(npts))
  do i = 1, npots

     write(fname,"(a,i2)") "fort.", 40 + i - 1
     open(unit=1,file=fname,form="formatted")
     read(1,*) char_dummy, i_dummy, dummy, ll(i), rc(i)
     do j = 1, npts
        read(1,fmt=*,iostat=stat) i_dummy, dummy, pswfs(j,i), v
        ! keep it as V and in hartree units
        !*** Need to cut off the non-coulomb tail!!!
        vps_raw(j) =  v 
        ! remove factor of r
        pswfs(j,i) =  pswfs(j,i) / r(j)
        if (stat /=0) stop "reading ps, vps"
     enddo
     close(1)
     vps(:,i) = vps_raw(:)
!     call cutoff_tail(r(:),vps_raw(:),zion,vps(:,i))
  enddo

  call xml_NewElement(xf,"semilocal-potentials")
  call my_add_attribute(xf,"npots-major",str(npots))
  call my_add_attribute(xf,"npots-minor","0")
  !         

  vpsd: do i = 1, npots
     call xml_NewElement(xf,"vps")
     call my_add_attribute(xf,"set","major")
     call my_add_attribute(xf,"n",str(nn(i)))
     call my_add_attribute(xf,"l",lsymb(ll(i)))
     call my_add_attribute(xf,"rc",str(rc(i)))
     call my_add_attribute(xf,"flavor",psflavor)

     call xml_NewElement(xf,"radfunc")
     call xml_NewElement(xf,"data")
     call xml_AddArray(xf,vps(1:npts,i))
     call xml_EndElement(xf,"data")
     call xml_EndElement(xf,"radfunc")
     call xml_EndElement(xf,"vps")
  enddo vpsd
  call xml_EndElement(xf,"semilocal-potentials")

  ! Dump of the pseudowave functions
  call xml_NewElement(xf,"pseudo-wave-functions")
  call my_add_attribute(xf,"npswfs",str(npots))
  call my_add_attribute(xf,"npswfs-major",str(npots))
  call my_add_attribute(xf,"npswfs-minor","0")

  ! Down pseudowave function follows

  pswfd: do i = 1, npots
     call xml_NewElement(xf,"pswf")
     call my_add_attribute(xf,"n",str(nn(i)))
     call my_add_attribute(xf,"l",lsymb(ll(i)))
     call my_add_attribute(xf,"set","major")

     call xml_NewElement(xf,"radfunc")

     call xml_NewElement(xf,"data")
     call xml_AddArray(xf,pswfs(1:npts,i))
     call xml_EndElement(xf,"data")
     call xml_EndElement(xf,"radfunc")
     call xml_EndElement(xf,"pswf")
  enddo pswfd
  call xml_EndElement(xf,"pseudo-wave-functions")

  call xml_NewElement(xf,"valence-charge")
  call my_add_attribute(xf,"total-charge",  &
                      str(total_valence_charge))
  call xml_NewElement(xf,"radfunc")

  call xml_NewElement(xf,"data")
  call xml_AddArray(xf,chval(1:npts))
  call xml_EndElement(xf,"data")
  call xml_EndElement(xf,"radfunc")
  call xml_EndElement(xf,"valence-charge")

  if (there_is_core) then
     allocate(chcore(npts))
     open(unit=1,file="fort.27",form="formatted")
     do i = 1, npts
        read(1,*) dummy, rho, rhop, rhopp
        chcore(i) = rho
     enddo
     close(1)

     call xml_NewElement(xf,"pseudocore-charge")
     call my_add_attribute(xf,"matching-radius",str(rcore))  
     call my_add_attribute(xf,"number-of-continuous-derivatives", &
                                    str(2)) ! ****
     call xml_NewElement(xf,"radfunc")

     call xml_NewElement(xf,"data")
     call xml_AddArray(xf,chcore(1:npts))
     call xml_EndElement(xf,"data")
     call xml_EndElement(xf,"radfunc")
     call xml_EndElement(xf,"pseudocore-charge")
     deallocate(chcore)
  endif

  call xml_EndElement(xf,"psml")


  call xml_Close(xf)

  deallocate(chval,r,vps,pswfs)

   CONTAINS

     subroutine my_add_attribute(xf,name,value)
       type(xmlf_t), intent(inout)   :: xf
       character(len=*), intent(in)  :: name
       character(len=*), intent(in)  :: value

       call xml_AddAttribute(xf,name,trim(value))
     end subroutine my_add_attribute

      FUNCTION SYMBOL( Z )

!! ** This function should not be called from within an I/O statement

! Given the atomic number, returns the atomic symbol (e.g. 'Na')
! Written by J. Soler

      character(len=2)    :: SYMBOL  ! Atomic symbol
      integer, intent(in) :: Z       ! Atomic number

      integer, parameter  :: NZ=103
      character(len=2), parameter :: NAME(NZ) =   &
              (/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',  &
                'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',  &
                'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',  &
                'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',  &
                'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',  &
                'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',  &
                'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',  &
                'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',  &
                'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',  &
                'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',  &
                'Md','No','Lr'/)

      IF (ABS(Z).LE.NZ) THEN
         SYMBOL = NAME(ABS(Z))
      ELSE
         WRITE(6,*) 'SYMBOL: ERROR: No data for Z =', Z
         SYMBOL = '??'
      ENDIF

      END function symbol

      subroutine add_zero_r(f,r,f0)
      ! Adds an r=0 element to a grid function f0, extrapolating

      double precision, intent(in)  :: f(:), r(:)
      double precision, intent(out) :: f0(:)
      
      integer i
      double precision :: r2

      do i = 2, npts
         f0(i) = f(i-1)
      enddo
      r2 = r(1)/(r(2)-r(1))
      f0(1) = f(2) - (f(3)-f(2))*r2

    end subroutine add_zero_r

    subroutine cutoff_tail(r,f,Z,fp)
      ! cuts off the non-coulombic tail of a pseudopotential
      ! one has to be careful and determine a safe area, avoiding
      ! early crossings and large-r mis-behavior

      real(dp), intent(in) :: r(:)
      real(dp), intent(in) :: f(:)
      real(dp), intent(in) :: Z
      real(dp), intent(out) :: fp(:)

      real(dp), parameter :: tolerance = 1.0e-4_dp

      integer :: n, j, jcut, nf, nflat
      real(dp) :: fcut, vp2z
      logical :: in_flat_region(size(r))

      n = size(r)

      fp(:) = f(:)

      in_flat_region(:) = .false.
      do j = 1, n
         vp2z = r(j)*f(j) + Z
         in_flat_region(j) = (abs(vp2z) .lt. tolerance)
      enddo
      nflat = 0
      do j=1,n
         if (in_flat_region(j)) nflat = nflat + 1
      enddo
      print *, "nflat: ", nflat

      if (nflat < 20) then
         print *, "Only ", nflat, " points below tolerance..."
         stop
      endif
      ! Choose a point safely into the flat region
      nf = 0
      do j = 1, n
         if (in_flat_region(j)) then
            nf = nf + 1
         endif
         if (nf > 5) then
            jcut = j
            exit
         endif
      enddo
      print *, "jcut: ", jcut
         
!      ...
!           Default cutoff function: f(r)=exp(-5*(r-r_cut)). It damps
!           down the residual of rV+2*Zion.
!           Should be made smoother... Vps ends up with a kink at rcut.
!           Maybe use one of the Vanderbilt generalized gaussians.

            do  j = jcut, n
               fcut = cutoff_function(r(j)-r(jcut))
               fp(j) = (-Z + fcut*(r(j)*f(j)+Z))/r(j)
            enddo
       end subroutine cutoff_tail

      function cutoff_function(r) result(x)

        real(dp), intent(in) :: r
        real(dp)             :: x

        x = exp(-5.0_dp*r)

      end function cutoff_function

    end program fhi2psml