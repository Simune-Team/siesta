program fhi2xml

  ! Reads fort.X files from FHIPP98 and generates
  ! a PSXML file

  ! Alberto Garcia, May 23, 2014

  use flib_wxml

  implicit none

  type(xmlf_t) :: xf

  integer, parameter :: dp = selected_real_kind(10,100)

  integer   :: i, j, npts, i_dummy, stat
  real(dp)  :: rho, rhop, rhopp, dummy
  real(dp)  :: znuc, zion, rcore, v
  real(dp), allocatable :: r(:), chval(:), chcore(:)
  real(dp), allocatable :: rvps(:,:), pswfs(:,:)

  character(len=100)    :: line, fname
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
     chval(i) = rho*r(i)*r(i)
  enddo


  call xml_OpenFile("FHIPSML",xf, indent=.false.)

  call xml_AddXMLDeclaration(xf,"UTF-8")

  call xml_NewElement(xf,"pseudo")

  call xml_NewElement(xf,"provenance")
  call my_add_attribute(xf,"creator","FHIPP98-2003")
  call my_add_attribute(xf,"translator","fhi2xml v0.1")
  call xml_NewElement(xf,"fort.20")
  open(1,file="fort.20",form="formatted",status="old", &
       position="rewind",action="read")
  do
     read(1,fmt="(a)",iostat=stat) line
     if (stat .ne. 0) exit
     call xml_AddPcData(xf,trim(line),line_feed=.true.)
  enddo
  close(1)
  call xml_EndElement(xf,"fort.20")
  !
  call xml_NewElement(xf,"fort.22")
  open(1,file="fort.22",form="formatted",status="old", &
       position="rewind",action="read")
  do
     read(1,fmt="(a)",iostat=stat) line
     if (stat .ne. 0) exit
     call xml_AddPcData(xf,trim(line),line_feed=.true.)
  enddo
  close(1)
  call xml_EndElement(xf,"fort.22")

  call xml_EndElement(xf,"provenance")

  print *, "Done provenance"

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

  ncp = ncore + 1
  do i = 1, norbs
     if (polarized) then
        read(1,*) n(i), l(i), fup(i), fdown(i)
        f(i) = fup(i) + fdown(i)
     else
        read(1,*) n(i), l(i), f(i)
     endif
  enddo
  read(1,*) lmax, pscode
  close(1)

  print *, "Done reading fort.20"

  zion = 0.0_dp
  do i = ncp, norbs
     zion = zion + f(i)
  enddo

  call xml_NewElement(xf,"ps-generation-configuration")
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
  call xml_EndElement(xf,"ps-generation-configuration")
          
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
  call my_add_attribute(xf,"symbol",nameat)
  call my_add_attribute(xf,"atomic-number",str(znuc))
  call my_add_attribute(xf,"zval",str(zion))
  call my_add_attribute(xf,"creator","FHIPP98-2003")
  call my_add_attribute(xf,"date","01-01-01")
  call my_add_attribute(xf,"flavor",psflavor)
  call my_add_attribute(xf,"relativistic",relattrib)
  call my_add_attribute(xf,"polarized",polattrib)
  call my_add_attribute(xf,"core-corrections",coreattrib)
  call my_add_attribute(xf,"xc-functional-type",xcfuntype)
  call my_add_attribute(xf,"xc-functional-parametrization", &
       xcfunparam)
  call xml_EndElement(xf,"header")

  call xml_NewElement(xf,"grid")
  call my_add_attribute(xf,"npts",str(npts))
  call xml_NewElement(xf,"grid_data")
  call xml_AddArray(xf,r(1:npts))
  call xml_EndElement(xf,"grid_data")
  call xml_EndElement(xf,"grid")

  allocate (rvps(npts,npots), pswfs(npts,npots))
  do i = 1, npots

     write(fname,"(a,i2)") "fort.", 40 + i - 1
     open(unit=1,file=fname,form="formatted")
     read(1,*) char_dummy, i_dummy, dummy, ll(i), rc(i)
     do j = 1, npts
        read(1,fmt=*,iostat=stat) i_dummy, dummy, pswfs(j,i), v
        ! make it rV and in rydberg units
        rvps(j,i) = 2.0_dp* v * r(j)
        if (stat /=0) stop "reading ps, vps"
     enddo
     close(1)
  enddo

  call xml_NewElement(xf,"semilocal")
  call my_add_attribute(xf,"units","Rydberg")
  call my_add_attribute(xf,"format","r*V")
  call my_add_attribute(xf,"npots-down",str(npots))
  call my_add_attribute(xf,"npots-up","0")

  !         

  vpsd: do i = 1, npots
     call xml_NewElement(xf,"vps")
     call my_add_attribute(xf,"principal-n",str(nn(i)))
     call my_add_attribute(xf,"l",lsymb(ll(i)))
     call my_add_attribute(xf,"cutoff",str(rc(i)))
     call my_add_attribute(xf,"occupation",str(ff(i)))
     call my_add_attribute(xf,"spin","-1")


     call xml_NewElement(xf,"radfunc")
     call xml_NewElement(xf,"data")
     call xml_AddArray(xf,rvps(1:npts,i))
     call xml_EndElement(xf,"data")
     call xml_EndElement(xf,"radfunc")
     call xml_EndElement(xf,"vps")
  enddo vpsd
  call xml_EndElement(xf,"semilocal")

  ! Dump of the pseudowave functions
  call xml_NewElement(xf,"pseudowave-functions")
  call my_add_attribute(xf,"units","none")
  call my_add_attribute(xf,"format","?")
  call my_add_attribute(xf,"n-pseudowave-functions-down", &
       str(npots))
  call my_add_attribute(xf,"n-pseudowave-functions-up","0")

  ! Down pseudowave function follows

  pswfd: do i = 1, npots
     call xml_NewElement(xf,"pswf")
     call my_add_attribute(xf,"principal-n",str(nn(i)))
     call my_add_attribute(xf,"l",lsymb(ll(i)))
     call my_add_attribute(xf,"spin","-1")

     call xml_NewElement(xf,"radfunc")

     call xml_NewElement(xf,"data")
     call xml_AddArray(xf,pswfs(1:npts,i))
     call xml_EndElement(xf,"data")
     call xml_EndElement(xf,"radfunc")
     call xml_EndElement(xf,"pswf")
  enddo pswfd
  call xml_EndElement(xf,"pseudowave-functions")

  call xml_NewElement(xf,"valence-charge")
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
        chcore(i) = rho*r(i)*r(i)
     enddo
     close(1)

     call xml_NewElement(xf,"pseudocore-charge")
     call xml_NewElement(xf,"radfunc")

     call xml_NewElement(xf,"data")
     call xml_AddArray(xf,chcore(1:npts))
     call xml_EndElement(xf,"data")
     call xml_EndElement(xf,"radfunc")
     call xml_EndElement(xf,"pseudocore-charge")
     deallocate(chcore)
  endif

  call xml_EndElement(xf,"pseudo")


  call xml_Close(xf)

  deallocate(chval,r,rvps,pswfs)

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

   end program fhi2xml
