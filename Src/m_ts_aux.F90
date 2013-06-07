! Auxilary functions for transiesta

module m_ts_aux
      
  use precision, only : dp

  implicit none

contains
  
  subroutine csolve(N,A,B,ipvt,ierr)

! INPUT
    integer, intent(in) :: N ! problem size
    complex(dp), intent(inout) :: A(N,N)

! INPUT/OUTPUT
    complex(dp), intent(inout) :: B(N,N)

! OUTPUT
    integer :: ipvt(N)           !pivoting vector
    integer :: ierr              ! error in inversion
! helpers, tempos ...
    call  csolveg(N,N,A,B,ipvt,ierr)

  end subroutine csolve


  subroutine csolveg(N,NRHS,A,B,ipvt,ierr)
      
! Finds the solution to a general complex N x N matrix equation by calling
! lapack or ESSL.
!
!     A.X = B;   returns X = Inv(A).B
!

! INPUT
    integer, intent(in) :: N, NRHS ! problem size
    complex(dp), intent(inout) :: A(N,N)

! INPUT/OUTPUT
    complex(dp), intent(inout) ::  B(N,NRHS)

! OUTPUT
    integer, intent(out) :: ipvt(N)           !pivoting vector
    integer, intent(out) :: ierr

! ESSL: 
!      CALL zgef(A,N,N,ipvt) ! factorize A for inv.
!      CALL zgesm('N',A,N,N,ipvt,B,N,N) ! inversion
! -----------
! LAPACK:
    CALL zgesv(N,NRHS,A,N,ipvt,B,N,ierr)
    
  end subroutine csolveg

!
!     Fermi Function
!
  function nf(z)
    complex(dp), intent(in) :: z
    complex(dp) :: nf
    nf = dcmplx(1._dp,0._dp)/(dcmplx(1._dp,0._dp) + exp(z))
  end function nf

  function nf1(z)
    real(dp), intent(in) :: z
    real(dp) :: nf1
    if ( z < -20._dp ) then
       nf1 = 1.0_dp
    else if ( z > 20._dp ) then
       nf1 = 0.0_dp
    else
       nf1 = 1._dp/(1._dp + exp(z))
    end if
  end function nf1

! ==================================================================
!
!  This subroutine returns points and weights for
!  Gaussian quadrature with the Fermi function as weight-function
!  for the integral starting at  0 kT and continuing to infinity
!  The integration variable must be (x-EF)/kT.
!
!  Roots and weights were calculated with Mathematica 
!  and only np <= 8 are possible.
!
! ==================================================================

  subroutine gaufermi0(np,r,w)
  
! INPUT
    integer, intent(in) :: np
    
! OUTPUT
    real(dp), intent(out) :: r(np), w(np)
    
! BEGIN

    if(np.EQ.1) then
       r(1)=1.083907926912384d0
       w(1)=0.6749972526421359d0
    else if(np.EQ.2) then
       r(1)=0.5382565051423951d0
       r(2)=2.542596344437450d0
       w(1)=0.4912393871473518d0
       w(2)=0.1837578654947843d0
    else if(np.EQ.3) then
       r(1)=0.3190572096226809d0
       r(2)=1.543368974102370d0
       r(3)=3.265580302380743d0
       w(1)=0.3364877459882430d0
       w(2)=0.2793814775930595d0
       w(3)=0.05912802906083348d0
    else if(np.EQ.4) then
       r(1)=0.2102419194790040d0
       r(2)=1.039620984654129d0
       r(3)=2.309891141094663d0
       r(4)=3.582671706563394d0
       w(1)=0.2375767017286547d0
       w(2)=0.2856670842991763d0
       w(3)=0.1249092420167874d0
       w(4)=0.02684422459751746d0
    else if(np.EQ.5) then
       r(1)=0.1487138474089773d0
       r(2)=0.7471353545323619d0
       r(3)=1.699933176888322d0
       r(4)=2.816341626065425d0
       r(5)=3.735921432431975d0
       w(1)=0.1745761309829550d0
       w(2)=0.2571062173241764d0
       w(3)=0.1662379825812521d0
       w(4)=0.06194888745327785d0
       w(5)=0.01512803430047453d0
    else if(np.EQ.6) then
       r(1)=0.1106442935476619d0
       r(2)=0.5620427996140370d0
       r(3)=1.299759391298092d0
       r(4)=2.216324366704579d0
       r(5)=3.142990449796749d0
       r(6)=3.819140236582460d0
       w(1)=0.1329364100888368d0
       w(2)=0.2212759675407208d0
       w(3)=0.1817594254611726d0
       w(4)=0.0941024437279011d0
       w(5)=0.03521803221252896d0
       w(6)=0.00970497361097562d0
    else if(np.EQ.7) then
       r(1)=0.08565069296593567d0
       r(2)=0.4383923908099421d0
       r(3)=1.025656129209570d0
       r(4)=1.779052498009138d0
       r(5)=2.603212676577210d0
       r(6)=3.358817473382426d0
       r(7)=3.869034064351875d0
       w(1)=0.1044673695976933d0
       w(2)=0.1880946611580349d0
       w(3)=0.1803883253731436d0
       w(4)=0.1166653019921832d0
       w(5)=0.05633616872844913d0
       w(6)=0.02229436891401237d0
       w(7)=0.006751056878619468d0
    else if(np.EQ.8) then
       r(1)=0.06505281409533842d0
       r(2)=0.3365059666811694d0
       r(3)=0.8007379593863171d0
       r(4)=1.419515309473741d0
       r(5)=2.136448210698268d0
       r(6)=2.867246805128238d0
       r(7)=3.496476499451854d0
       r(8)=3.899485946812052d0
       w(1)=0.08040619745384284d0
       w(2)=0.1553971715533348d0
       w(3)=0.1701652217622570d0
       w(4)=0.1322663627104325d0
       w(5)=0.07822398307939787d0
       w(6)=0.03772561389875811d0
       w(7)=0.01575324256952425d0
       w(8)=0.005059459614588813d0
    else
       call die('ERROR in contour: No. points in Gauss-Fermi not valid')
    end if

  end subroutine gaufermi0

! ==================================================================
!
!  This subroutine returns points and weights for
!  Gaussian quadrature with the Fermi function as weight-function
!  for the integral starting at - 2 kT and continuing to infinity
!  The integration variable must be (x-EF)/kT.
!
!  Roots and weights were calculated with Mathematica 
!  and only np <= 10 are possible.
!
! ==================================================================

  subroutine gaufermi2(np,r,w)

! INPUT
    integer, intent(in) :: np

! OUTPUT
    real(dp), intent(out) :: r(np), w(np)

! BEGIN
    
    if(np.EQ.1) then
       r(1)=-0.3939563898360445d0
       w(1)= 2.108778083125163d0
    else if(np.EQ.2) then
       r(1)=-1.170526120464507d0
       r(2)=1.538491358782026d0
       w(1)=1.504273593762035d0
       w(2)=0.6045044893631276d0
    else if(np.EQ.3) then
       r(1)=-1.504503786523456d0
       r(2)=0.2468248611220023d0
       r(3)=2.697778990181126d0
       w(1)=0.997660057171304d0
       w(2)=0.949562501130609d0
       w(3)=0.1615555248232495d0
    else if(np.EQ.4) then
       r(1)=-1.673943849304428d0
       r(2)=-0.4474674370058036d0
       r(3)=1.330534184865143d0
       r(4)=3.266373721011666d0
       w(1)=0.6867074911761266d0
       w(2)=0.951474287837990d0
       w(3)=0.4094587617749555d0
       w(4)=0.06113754233609050d0
    else if(np.EQ.5) then
       r(1)=-1.770062653937712d0
       r(2)=-0.8709555383368669d0
       r(3)=0.4901712495177842d0
       r(4)=2.093895368503766d0
       r(5)=3.545170037285158d0
       w(1)=0.4955853946090569d0
       w(2)=0.8273713699688027d0
       w(3)=0.5754603330610385d0
       w(4)=0.1796939640008372d0
       w(5)=0.03066702148542732d0
    else if(np.EQ.6) then
       r(1)=-1.829439518979821d0
       r(2)=-1.146204397890313d0
       r(3)=-0.07131730547282354d0
       r(4)=1.235731668828154d0
       r(5)=2.610998494506139d0
       r(6)=3.694342270990330d0
       w(1)=0.3725590149625833d0
       w(2)=0.6892484328427096d0
       w(3)=0.6285588378016490d0
       w(4)=0.3102264543081084d0
       w(5)=0.08988471674032117d0
       w(6)=0.01830062646979192d0
    else if(np.EQ.7) then
       r(1)=-1.868566706265513d0
       r(2)=-1.333750753492947d0
       r(3)=-0.4671819299710864d0
       r(4)=0.6155596633627794d0
       r(5)=1.810132047792538d0
       r(6)=2.959567717565726d0
       r(7)=3.781636349808305d0
       w(1)=0.2895226487127946d0
       w(2)=0.5697417025838409d0
       w(3)=0.6103950906250394d0
       w(4)=0.4056872955348393d0
       w(5)=0.1699103200674708d0
       w(6)=0.05135034388595542d0
       w(7)=0.01217068171522311d0
    else if(np.EQ.8) then
       r(1)=-1.895672249027231d0
       r(2)=-1.466614717531430d0
       r(3)=-0.7561290787790193d0
       r(4)=0.1533431699819214d0
       r(5)=1.182947486839513d0
       r(6)=2.246075346436202d0
       r(7)=3.198653145685052d0
       r(8)=3.836613235441707d0
       w(1)=0.2311174008433614d0
       w(2)=0.4732806927315964d0
       w(3)=0.5608286729022941d0
       w(4)=0.4533435994189919d0
       w(5)=0.2495564466955304d0
       w(6)=0.0993480484945936d0
       w(7)=0.03260622012705790d0
       w(8)=0.008697001911737815d0
    else if(np.EQ.9) then
       r(1)=-1.914824523014234d0
       r(2)=-1.562038368765970d0
       r(3)=-0.968910863577704d0
       r(4)=-0.1949288400949761d0
       r(5)=0.6979256033432659d0
       r(6)=1.650096838706900d0
       r(7)=2.580856447048328d0
       r(8)=3.369436690752455d0
       r(9)=3.873842852170074d0
       w(1)=0.1894170780936077d0
       w(2)=0.3981890665581406d0
       w(3)=0.5032913972644507d0
       w(4)=0.4625541099882755d0
       w(5)=0.3100463445522872d0
       w(6)=0.1543598572529525d0
       w(7)=0.06213664525098128d0
       w(8)=0.02227715398894096d0
       w(9)=0.006506430175526240d0
    else if(np.EQ.10) then
       r(1)=-1.935027961292951d0
       r(2)=-1.664572001282213d0
       r(3)=-1.201668446818815d0
       r(4)=-0.5760373502793187d0
       r(5)=0.1783245856273426d0
       r(6)=1.022542228584406d0
       r(7)=1.906103776470704d0
       r(8)=2.751594285665453d0
       r(9)=3.451723464834261d0
       r(10)=3.891156181888537d0
       w(1)=0.1449623786381726d0
       w(2)=0.3121591556563253d0
       w(3)=0.4225212785341255d0
       w(4)=0.4455687626558362d0
       w(5)=0.3672963596094461d0
       w(6)=0.2312162310669771d0
       w(7)=0.1139204247840664d0
       w(8)=0.04754537984898786d0
       w(9)=0.01805810935560718d0
       w(10)=0.005530002975619004d0
    else
       call die('ERROR in contour: No. points in Gauss-Fermi not valid')
    end if

  end subroutine gaufermi2


!==================================================================
!
!  This subroutine returns points and weights for
!  Gaussian quadrature with the Fermi function as weight-function
!  for the integral starting at -10 kT and continuing to infinity
!  The integration variable must be (x-EF)/kT.
!
!  Roots and weights were calculated with Mathematica 
!  and only np <= 11 are possible.
!
!==================================================================
  subroutine gaufermi10(np,r,w)
  
! INPUT
    integer, intent(in) :: np

! OUTPUT
    real(dp), intent(out) :: r(np), w(np)

! BEGIN

    if(np.EQ.1) then
       r(1)=0.0d0
       w(1)= 10.00004539889922d0
    else if(np.EQ.2) then
       r(1)=-7.589291604880055d0
       r(2)=-1.221802100721404d0
       w(1)=5.675311882692906d0
       w(2)=4.324733516206313d0
    else if(np.EQ.3) then
       r(1)=-8.585496555093303d0
       r(2)=-3.820530664777789d0
       r(3)=0.986066695581377d0
       w(1)=3.477526165098058d0
       w(2)=5.186811694320817d0
       w(3)=1.335707539480341d0
    else if(np.EQ.4) then
       r(1)=-9.04515460499372d0
       r(2)=-5.50546361929916d0
       r(3)=-1.157816323315441d0
       r(4)=3.295961913561931d0
       w(1)=2.387438931678946d0
       w(2)=4.34278682840718d0
       w(3)=3.060088862366506d0
       w(4)=0.2097307764465852d0
    else if(np.EQ.5) then
       r(1)=-9.29715099738137d0
       r(2)=-6.568033857374533d0
       r(3)=-2.735520320810778d0
       r(4)=0.934176762578089d0
       r(5)=6.006574317458234d0
       w(1)=1.772164747565569d0
       w(2)=3.511187681918489d0
       w(3)=3.627022219007695d0
       w(4)=1.07272306589775d0
       w(5)=0.01694768450971697d0
    else if(np.EQ.6) then
       r(1)=-9.45209964409028d0
       r(2)=-7.268242767388545d0
       r(3)=-3.963938304412221d0
       r(4)=-0.531068137612203d0
       r(5)=3.132688689824056d0
       r(6)=8.993027010392847d0
       w(1)=1.388124608118496d0
       w(2)=2.87983050576294d0
       w(3)=3.481266741849019d0
       w(4)=2.063661054137707d0
       w(5)=0.1862097379558001d0
       w(6)=0.000952751075253517d0
    else if(np.EQ.7) then
       r(1)=-9.55532055243929d0
       r(2)=-7.753725776689058d0
       r(3)=-4.904511670069002d0
       r(4)=-1.685779131853237d0
       r(5)=1.44864595543743d0
       r(6)=5.647385498519065d0
       r(7)=12.14658196402703d0
       w(1)=1.130015410503081d0
       w(2)=2.4110492495816d0
       w(3)=3.149315366918614d0
       w(4)=2.644314743151264d0
       w(5)=0.6472142083156056d0
       w(6)=0.01809243818046256d0
       w(7)=0.00004398224858996766d0
    else if(np.EQ.8) then
       r(1)=-9.62822843268492d0
       r(2)=-8.105484817611283d0
       r(3)=-5.627717747784063d0
       r(4)=-2.656520346093119d0
       r(5)=0.2633552229335414d0
       r(6)=3.616143170792033d0
       r(7)=8.393930661380015d0
       r(8)=15.41279614015821d0
       w(1)=0.94666040430034d0
       w(2)=2.057582109167864d0
       w(3)=2.811270590144065d0
       w(4)=2.809366434050096d0
       w(5)=1.269983201747899d0
       w(6)=0.1038955516793378d0
       w(7)=0.001285321320416654d0
       w(8)=1.786489199868805d-6
    else if(np.EQ.9) then
       r(1)=-9.6820652654789d0
       r(2)=-8.369781210774285d0
       r(3)=-6.192345665037222d0
       r(4)=-3.473805652661154d0
       r(5)=-0.6831205856841883d0
       r(6)=2.205969578015823d0
       r(7)=6.046663086537123d0
       r(8)=11.29297056119824d0
       r(9)=18.76137286350282d0
       w(1)=0.8107341271298983d0
       w(2)=1.785079818170019d0
       w(3)=2.51103258451205d0
       w(4)=2.750687100048243d0
       w(5)=1.809690141620464d0
       w(6)=0.3221677848019469d0
       w(7)=0.01057747756805342d0
       w(8)=0.00007629887742744467d0
       w(9)=6.617111537452023d-8
    else if(np.EQ.10) then
       r(1)=-9.72322591294335d0
       r(2)=-8.57437551293854d0
       r(3)=-6.641129380914527d0
       r(4)=-4.159028746042904d0
       r(5)=-1.494255554367855d0
       r(6)=1.147456409041345d0
       r(7)=4.375411858099028d0
       r(8)=8.664596261501208d0
       r(9)=14.30216830430639d0
       r(10)=22.173270847707d0
       w(1)=0.7065181875587684d0
       w(2)=1.570348197982865d0
       w(3)=2.254158866270095d0
       w(4)=2.604370153468784d0
       w(5)=2.143421464243213d0
       w(6)=0.6740721035836645d0
       w(7)=0.04630860182648593d0
       w(8)=0.0008438204906516487d0
       w(9)=4.001190273759472d-6
       w(10)=2.284414108877665d-9
    else if(np.EQ.11) then
       r(1)=-9.75558463809598d0
       r(2)=-8.736712816014405d0
       r(3)=-7.004092507085979d0
       r(4)=-4.734305854683416d0
       r(5)=-2.207092188781701d0
       r(6)=0.2974975611383293d0
       r(7)=3.109612898457216d0
       r(8)=6.758483175691933d0
       r(9)=11.41448003156238d0
       r(10)=17.39664226438993d0
       r(11)=25.63558657929966d0
       w(1)=0.6244149508070781d0
       w(2)=1.39775477661742d0
       w(3)=2.036301124424707d0
       w(4)=2.433863384193428d0
       w(5)=2.285934600421204d0
       w(6)=1.080706625117594d0
       w(7)=0.1362296513422132d0
       w(8)=0.004782303192092462d0
       w(9)=0.00005779176576227073d0
       w(10)=1.90943134104802d-7
       w(11)=7.458786235767771d-11
    else
       call die('ERROR in contour: No. points in Gauss-Fermi not valid')
    end if

  end subroutine gaufermi10


!==================================================================
!
!  This subroutine returns points and weights for
!  Gaussian quadrature with the Fermi function as weight-function
!  for the integral starting at -20 kT and continuing to infinity
!  The integration variable must be (x-EF)/kT.
!
!  Roots and weights were calculated with Mathematica 
!  and only np <= 12 are possible.
!
!==================================================================

  subroutine gaufermi20(np,r,w)

! INPUT
    integer, intent(in) :: np
      
! OUTPUT
    real(dp), intent(out) :: r(np), w(np)
      
! BEGIN
    
    if(np.EQ.1) then
       r(1)=0d0
       w(1)=20.00000000206115d0
    else if(np.EQ.2) then
       r(1)=-15.59867300181492d0
       r(2)=-3.761795334794256d0
       w(1)=10.4013205793134d0
       w(2)=9.59867942274775d0
    else if(np.EQ.3) then
       r(1)=-17.56069617714904d0
       r(2)=-9.2132959955567d0
       r(3)=-1.155162436732018d0
       w(1)=6.008953676528447d0
       w(2)=9.51482429458606d0
       w(3)=4.476222030946645d0
    else if(np.EQ.4) then
       r(1)=-18.43109871253151d0
       r(2)=-12.55785501966099d0
       r(3)=-5.008576968354179d0
       r(4)=0.6703031973794247d0
       w(1)=3.928604574772423d0
       w(2)=7.327554038101563d0
       w(3)=7.006527638623332d0
       w(4)=1.737313750563833d0
    else if(np.EQ.5) then
       r(1)=-18.89086253068388d0
       r(2)=-14.55265566977275d0
       r(3)=-8.251795056131431d0
       r(4)=-2.226669302745686d0
       r(5)=2.575330427221925d0
       w(1)=2.800001804214125d0
       w(2)=5.634869133395311d0
       w(3)=6.590038222067974d0
       w(4)=4.558760058277085d0
       w(5)=0.4163307841066546d0
    else if(np.EQ.6) then
       r(1)=-19.16408385964757d0
       r(2)=-15.81251227255939d0
       r(3)=-10.62333214685802d0
       r(4)=-4.894452357027714d0
       r(5)=-0.2096250897412564d0
       r(6)=4.860746544780049d0
       w(1)=2.120046040392251d0
       w(2)=4.449478155957893d0
       w(3)=5.711100387817238d0
       w(4)=5.411409622995354d0
       w(5)=2.254377505926917d0
       w(6)=0.05358828897149676d0
    else if(np.EQ.7) then
       r(1)=-19.34059078585103d0
       r(2)=-16.65569013503401d0
       r(3)=-12.33651621649692d0
       r(4)=-7.191313942517983d0
       r(5)=-2.328646905154673d0
       r(6)=1.650831568040055d0
       r(7)=7.48601526930274d0
       w(1)=1.677225432192864d0
       w(2)=3.612032868283468d0
       w(3)=4.890449861748301d0
       w(4)=5.214321072289275d0
       w(5)=3.909057363121214d0
       w(6)=0.6925216813788421d0
       w(7)=0.004391723047198199d0
    else if(np.EQ.8) then
       r(1)=-19.46189127711878d0
       r(2)=-17.24833759177651d0
       r(3)=-13.59760191684534d0
       r(4)=-9.04368442685297d0
       r(5)=-4.321229543192369d0
       r(6)=-0.3263829332688342d0
       r(7)=3.765173397835948d0
       r(8)=10.34189376685321d0
       w(1)=1.371330711552274d0
       w(2)=3.003944350247329d0
       w(3)=4.207819513446733d0
       w(4)=4.778581486266408d0
       w(5)=4.437164268712301d0
       w(6)=2.086403937563884d0
       w(7)=0.1144814526588105d0
       w(8)=0.0002742816134181913d0
    else if(np.EQ.9) then
       r(1)=-19.54929920213494d0
       r(2)=-17.68197584463833d0
       r(3)=-14.54887398428311d0
       r(4)=-10.51940501772798d0
       r(5)=-6.103558181297253d0
       r(6)=-1.992655192194112d0
       r(7)=1.525571123999186d0
       r(8)=6.193010336092172d0
       r(9)=13.358832440893d0
       w(1)=1.150127890629439d0
       w(2)=2.549412354432156d0
       w(3)=3.654799257745899d0
       w(4)=4.319442139790358d0
       w(5)=4.389145910496666d0
       w(6)=3.255216039160415d0
       w(7)=0.670167198771067d0
       w(8)=0.01167490049403081d0
       w(9)=0.00001431054112259171d0
    else if(np.EQ.10) then
       r(1)=-19.61468272035192d0
       r(2)=-18.00990625972198d0
       r(3)=-15.28374170887942d0
       r(4)=-11.70144577307719d0
       r(5)=-7.631842448220589d0
       r(6)=-3.567593257971375d0
       r(7)=-0.07719172772979343d0
       r(8)=3.590736265289387d0
       r(9)=8.845658657625209d0
       r(10)=16.49542210281781d0
       w(1)=0.98424471631702d0
       w(2)=2.200538114598689d0
       w(3)=3.207319415278673d0
       w(4)=3.895860252445527d0
       w(5)=4.161376340361231d0
       w(6)=3.741730275965698d0
       w(7)=1.69116571629326d0
       w(8)=0.116866651654976d0
       w(9)=0.0008978640577459201d0
       w(10)=6.55088330179433d-7
    else if(np.EQ.11) then
       r(1)=-19.66508766518369d0
       r(2)=-18.26476095764296d0
       r(3)=-15.86380413793839d0
       r(4)=-12.65866266515433d0
       r(5)=-8.924531074202514d0
       r(6)=-5.022612993647351d0
       r(7)=-1.450026876652112d0
       r(8)=1.766952686027662d0
       r(9)=5.932823052968517d0
       r(10)=11.65525938581522d0
       r(11)=19.72466968095453d0
       w(1)=0.8561241966042454d0
       w(2)=1.92646055850799d0
       w(3)=2.842542307761985d0
       w(4)=3.521870657374105d0
       w(5)=3.886694053012067d0
       w(6)=3.807609024119133d0
       w(7)=2.636620120451863d0
       w(8)=0.5090785166986918d0
       w(9)=0.0129428062879975d0
       w(10)=0.00005773411646617993d0
       w(11)=2.712661548517718d-8
    else if(np.EQ.12) then
       r(1)=-19.70491970277043d0
       r(2)=-18.46740137216554d0
       r(3)=-16.33047296145924d0
       r(4)=-13.44345073197975d0
       r(5)=-10.01690683045319d0
       r(6)=-6.325019889856077d0
       r(7)=-2.744101111786821d0
       r(8)=0.3777065401860518d0
       r(9)=3.827502559994769d0
       r(10)=8.477132233008886d0
       r(11)=14.58255327804397d0
       r(12)=23.02779435446008d0
       w(1)=0.7547335534779741d0
       w(2)=1.70675470412558d0
       w(3)=2.542100551909292d0
       w(4)=3.196832534647743d0
       w(5)=3.610860082637905d0
       w(6)=3.705337391439684d0
       w(7)=3.159915864329894d0
       w(8)=1.236955758816305d0
       w(9)=0.08540068678220045d0
       w(10)=0.00110561158564446d0
       w(11)=3.261271441086333d-6
       w(12)=1.037491566971821d-9
    else
       call die('ERROR in contour: No. points in Gauss-Fermi not valid')
    end if

  end subroutine gaufermi20


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  gauss.f: Points and weights for Gaussian quadrature                 c
!								       c
!  taken from: "Projects in Computational Physics" by Landau and Paez  c 
!	       copyrighted by John Wiley and Sons, New York            c      
!                                                                      c
!  written by: Oregon State University Nuclear Theory Group            c
!	       Guangliang He & Rubin H. Landau                         c
!  supported by: US National Science Foundation, Northwest Alliance    c
!                for Computational Science and Engineering (NACSE),    c
!                US Department of Energy 	                       c
!								       c
!  comment: error message occurs if subroutine called without a main   c
!  comment: this file has to reside in the same directory as integ.c   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     rescale rescales the gauss-legendre grid points and weights
!
!     npts     number of points
!     job = 0  rescalling uniformly between (a,b)
!           1  for integral (0,b) with 50% points inside (0, ab/(a+b))
!           2  for integral (a,inf) with 50% inside (a,b+2a)
!     x, w     output grid points and weights.
!
  subroutine gauss(npts,job,a,b,x,w) 
    integer, intent(in) :: npts, job
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: x(npts),w(npts)
    real(dp) :: xi, t,t1,pp,p1,p2,p3,aj

    integer :: m,i,j 
    real(dp), parameter :: eps = 3.0d-14
    real(dp), parameter :: pi = 3.14159265358979323846264338328_dp
    real(dp), parameter :: zero = 0._dp, one = 1._dp, two = 2._dp
    real(dp), parameter :: half = 0.5_dp, quarter = .25_dp

!
! *** FIRST EXECTUABLE *************************************************
!
    m = (npts+1)/2
    do i = 1 , m
       t  = cos(pi*(i-quarter)/(npts+half))
       t1 = huge(1._dp)
       do while ( abs(t-t1) > eps )
          p1 = one
          p2 = zero
          aj = zero
          do j = 1 , npts
             p3 = p2
             p2 = p1
             aj = aj + one
             p1 = ((two*aj-one)*t*p2-(aj-one)*p3)/aj
          end do
          pp = npts*(t*p1-p2)/(t*t-one)
          t1 = t
          t  = t1 - p1/pp
       end do

       x(i)        = -t
       x(npts+1-i) = t
       w(i)        = two/((one-t*t)*pp*pp)
       w(npts+1-i) = w(i)
    end do

!
! rescale the grid points 
    if (job.eq.0) then
!     scale to (a,b) uniformly
       do i = 1 , npts
          x(i) = x(i)*(b-a)/two+(b+a)/two
          w(i) = w(i)*(b-a)/two
       end do
    else if (job.eq.1) then
! scale to (0,b) with 50% points inside (0,ab/(a+b))
       do i = 1 , npts
          xi   = x(i)
          x(i) = a*b*(one+xi)/(b+a-(b-a)*xi)
          w(i) = w(i)*two*a*b*b/((b+a-(b-a)*xi)*(b+a-(b-a)*xi))
       end do
    else if (job.eq.2) then
! scale to (a,inf) with 50% points inside (a,b+2a)
       do i = 1 , npts
          xi   = x(i)
          x(i) = (b*xi+b+a+a)/(one-xi)
          w(i) = w(i)*two*(a+b)/((one-xi)*(one-xi))
       end do
    else
       call die('gauss: Wrong value of job for Gaussian quadrature')
    endif
    
  end subroutine gauss

end module m_ts_aux
