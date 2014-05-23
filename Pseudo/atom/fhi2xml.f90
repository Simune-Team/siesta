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
  real(dp)  :: znuc, zion, rcore
  real(dp), allocatable :: r(:), chval(:), chcore(:)
  real(dp), allocatable :: vps(:,:), pswfs(:,:)

  character(len=100)    :: line, fname
  character*4      :: polattrib, relattrib, coreattrib
  character*1      :: pscode, char_dummy
  character(len=40):: psflavor

  character*30 xcfuntype, xcfunparam
  integer          :: ncore, nval, ncp, norbs, lmax, npots
  integer          :: xc_code
  integer, allocatable  :: n(:), l(:)
  integer, allocatable  :: nn(:), ll(:)
  real(dp), allocatable :: f(:), ff(:), rc(:)

  logical :: tdopsp, nonrel, polarized, there_is_core
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


  call xml_OpenFile("VPSXML",xf, indent=.false.)

  call xml_AddXMLDeclaration(xf,"UTF-8")

  call xml_NewElement(xf,"pseudo")

  call xml_NewElement(xf,"provenance")
  call my_add_attribute(xf,"creator","FHIPP98-2003")
  call xml_NewElement(xf,"fort.20")
  open(1,file="fort.20",form="formatted",status="old", &
       position="rewind",action="read")
  do
     read(1,fmt="(a)",iostat=stat) line
     if (stat .ne. 0) exit
     call xml_AddPcData(xf,line,line_feed=.true.)
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
     call xml_AddPcData(xf,line,line_feed=.true.)
  enddo
  close(1)
  call xml_EndElement(xf,"fort.22")
  call xml_EndElement(xf,"provenance")

  open(1,file="fort.22",form="formatted",status="old", &
       position="rewind",action="read")
  read(1,*) znuc, ncore, nval, xc_code, rcore
  norbs = ncore + nval
  allocate (n(norbs), l(norbs), f(norbs))
  ncp = ncore + 1
  do i = 1, norbs
     read(1,*) n(i), l(i), f(i)
  enddo
  read(1,*) lmax, psflavor

  npots = lmax + 1
  allocate (rc(npots), ll(npots), nn(npots), ff(npots))
  do i = 1, nval
     nn(i) = n(ncore + i)
     ff(i) = f(ncore + i)
  enddo
  ! horrible kludge for now
  do i = nval + 1, npots
     nn(i) = n(ncore+nval) + 1
     ff(i) = 0.0_dp
  enddo
          

  select case (pscode)
  case ("h","H")
     psflavor ="Hamann GNCPP"
  case ("t","T")
     psflavor ="Troullier-Martins"
  end select

  zion = 0.0_dp
  do i = ncp, norbs
     zion = zion + f(i)
  enddo
  close(1)
  !
  !
  !
  open(1,file="fort.20",form="formatted",status="old", &
       position="rewind",action="read")
  read(1,*) tdopsp, nonrel
  read(1,*) polarized
  close(1)


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
  call my_add_attribute(xf,"symbol","Symbol")
  call my_add_attribute(xf,"atomic-number",str(znuc))
  call my_add_attribute(xf,"zval",str(zion))
  call my_add_attribute(xf,"creator","FHIPP98-2003")
  call my_add_attribute(xf,"date","1-1-1")
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

  allocate (vps(npts,npots), pswfs(npts,npots))
  do i = 1, npots

     write(fname,"(a,i2)") "fort.", 40 + i - 1
     open(unit=1,file=fname,form="formatted")
     read(1,*) char_dummy, i_dummy, dummy, ll(i), rc(i)
     do j = 1, npts
        read(1,fmt=*,iostat=stat) i_dummy, dummy, pswfs(j,i), vps(j,i)
        if (stat /=0) stop "reading ps, vps"
     enddo
     close(1)
  enddo

  call xml_NewElement(xf,"semilocal")
  call my_add_attribute(xf,"units","Rydberg?")
  call my_add_attribute(xf,"format","V")
  call my_add_attribute(xf,"npots-down",str(npots))
  call my_add_attribute(xf,"npots-up","0")

  !         

  vpsd: do i = 1, npots
     call xml_NewElement(xf,"vps")
     call my_add_attribute(xf,"principal-n",str(nn(i)))
     call my_add_attribute(xf,"l",str(ll(i)))
     call my_add_attribute(xf,"cutoff",str(0.0))
     call my_add_attribute(xf,"occupation",str(ff(i)))
     call my_add_attribute(xf,"spin","-1")


     call xml_NewElement(xf,"radfunc")
     call xml_NewElement(xf,"data")
     call xml_AddArray(xf,vps(1:npts,i))
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
     call my_add_attribute(xf,"l",str(ll(i)))
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
        read(20,*) dummy, rho, rhop, rhopp
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

  deallocate(chval,r,vps,pswfs)

   CONTAINS

     subroutine my_add_attribute(xf,name,value)
       type(xmlf_t), intent(inout)   :: xf
       character(len=*), intent(in)  :: name
       character(len=*), intent(in)  :: value

       call xml_AddAttribute(xf,name,trim(value))
     end subroutine my_add_attribute

   end program fhi2xml
