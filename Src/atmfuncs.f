      module atmfuncs

C     This file contains a set of routines which provide all the information
C     about the basis set, pseudopotential, atomic mass, etc... of all the
C     chemical species present in the calculation.

C     Written by D. Sanchez-Portal, Sept. 1998
C     Modified by  DSP, July 1999

C     The routines contained in this file can only be called after they
C     are initialized by calling the subroutine 'atom' for all the 
C     different chemical species in the calculation:

      use precision

      implicit none 
!
!----------------------------------------------------------------
!
!    Hard-wired parameters to be eliminated in the future
!
C INTEGER NSMAX     : Maximum number of species.
C INTEGER NTBMAX    : Maximum number of points in the tables defining
C                      orbitals,projectors and local neutral-atom 
C                      pseudopotential.
C INTEGER LMAXD     : Maximum angular momentum for both orbitals and 
C                      projectors.
C INTEGER  NZETMX   : Maximum number of PAOs or polarization orbitals
C                      with the same angular  momentum and 
C                      for the same specie.       
C INTEGER  NKBMX    : Maximum number of Kleinman-Bylander projectors
C                      for each angular momentum
C INTEGER  NSMX    : Maximum number of semicore shells for each angular
C                    momentum present in the atom ( for normal atom nsmx=0)
C INTEGER  NRMAX    : Maximum number of points in the functions read
C                      from file '.vps' or '.psatom.data' (this number is
C                      determined by the parameter nrmax in the
C                      program atm, which generates the files with
C                      the pseudopotential information).

       integer  ntbmax, lmaxd, nzetmx, nkbmx, nrmax, lmx2, nsemx, nsmx
         parameter ( ntbmax =  500 )
         parameter ( lmaxd  =    4 )
         parameter ( nzetmx =    3 )
         parameter ( nkbmx  =    2 )
         parameter ( nsmx  =    1 )
         parameter ( nsemx = 1 + nsmx) 
         parameter ( nrmax  = 2000 )
         parameter ( lmx2   = (lmaxd+1)*(lmaxd+1) )
!----------------------------------------------------------------

      integer, parameter      :: maxos=2*nzetmx*lmx2*nsemx

      integer, save      :: Node     ! To be initialized by calling atom

      integer, save                   ::  ismax
      integer, save                   ::  nsmax
      integer, save, allocatable      ::  izsave(:)
      integer, save, allocatable      ::  nomax(:)
      integer, save, allocatable      ::  nkbmax(:)

      integer, save, allocatable      ::  lmxosave(:)
      integer, save, allocatable      ::  nkblsave(:,:)
      integer, save, allocatable      ::  npolorbsave(:,:,:)
      integer, save, allocatable      ::  lsemicsave(:,:)
      integer, save, allocatable      ::  nzetasave(:,:,:)
      integer, save, allocatable      ::  cnfigtb(:,:,:)

      logical, save, allocatable      ::  semicsave(:)
     
      integer, save, allocatable           :: izvaltb(:)
      integer, save, allocatable           :: lmxkbsave(:)
      real*8, save, allocatable            :: smasstb(:)
      real*8, save, allocatable            :: chargesave(:)
      real*8, save, allocatable            :: slfe(:)
      real*8, save, allocatable            :: lambdatb(:,:,:,:)

      double precision, save, allocatable  :: qtb(:,:)
      double precision, save, allocatable  :: qltb(:,:,:)

      real*8, save, allocatable  :: table(:,:,:)
      real*8, save, allocatable  :: tabpol(:,:,:)
      real*8, save, allocatable  :: tab2(:,:,:)
      real*8, save, allocatable  :: tab2pol(:,:,:)


      real*8, save, allocatable  ::  coretab(:,:,:)
      real*8, save, allocatable  ::  chloctab(:,:,:)
      real*8, save, allocatable  ::  corrtab(:,:,:)
      real*8, save, allocatable  ::  rctb(:,:,:)
      real*8, save, allocatable  ::  rcotb(:,:,:,:)
      real*8, save, allocatable  ::  rcpoltb(:,:,:,:)
!
!     PGI compiler cannot allocate these!!!!
!
      character(len=20), save    :: label_save(40)
      character(len=10), save    :: basistype_save(40)  


      contains
!
!
      subroutine check_is(name,is)
      character(len=*), intent(in) :: name
      integer, intent(in) :: is

      if ((is.lt.1).or.(is.gt.ismax)) then 
         if (Node.eq.0) then
            write(6,*) trim(name),': THERE ARE NO DATA FOR IS=',IS
            write(6,*) trim(name),': ISMIN= 1, ISMAX= ',ismax
         endif
         call die
      endif
      end subroutine check_is
!
!
      FUNCTION IZOFIS( IS )
      integer :: izofis ! Atomic number
      integer, intent(in) :: is ! Species index

      call check_is('izofis',is)
      izofis=izsave(is)

      end function izofis
!
      FUNCTION IZVALFIS( IS )
      integer :: izvalfis          ! Valence charge
      integer, intent(in) :: is            ! Species index

      call check_is('izvalfis',is)
 
      izvalfis=izvaltb(is)
      end function izvalfis
!
      FUNCTION LABELFIS (IS)
      character(len=2) ::  labelfis  ! Atomic label
      integer, intent(in) :: is            ! Species index

      call check_is('labelfis',is)
      labelfis=label_save(is)
      end function labelfis
!
      FUNCTION LMXKBFIS (IS)
      integer :: lmxkbfis    ! Maximum ang mom of the KB projectors
      integer, intent(in) :: is            ! Species index

      call check_is('lmxkbfis',is)
      lmxkbfis=lmxkbsave(is)
      end function lmxkbfis
!
      FUNCTION LOMAXFIS (IS)
      integer :: lomaxfis  ! Maximum ang mom of the Basis Functions
      integer, intent(in) :: is            ! Species index

      integer lmx, nsm

      call check_is('lomaxfis',is)

      lomaxfis=0           
      lmx=lmxosave(is)
      do nsm=1,lsemicsave(lmx,is)+1
         if(npolorbsave(lmx,nsm,is).gt.0)   lomaxfis=lmx+1
      enddo     
      
      lomaxfis=max(lomaxfis,lmx)
      end function lomaxfis
!

      FUNCTION MASSFIS(IS)
      real*8 :: massfis            ! Mass
      integer, intent(in) :: is            ! Species index

      call check_is('massfis',is)
      massfis=smasstb(is)
      end function massfis
!
      FUNCTION NKBFIS(IS)
      integer :: nkbfis    ! Total number of KB projectors
      integer, intent(in) :: is            ! Species index

      call check_is('nkbfis',is)
      nkbfis=nkbmax(is)
      end function nkbfis
!

      FUNCTION NOFIS(IS)
      integer :: nofis    ! Total number of Basis functions
      integer, intent(in) :: is            ! Species index

      call check_is('nofis',is)
      nofis=nomax(is)
      end function nofis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! AMENoFIS

      FUNCTION ATMPOPFIO (IS,IO)
      real*8 atmpopfio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns the population of the atomic basis orbitals in the atomic 
C ground state configuration.

      call check_is('atmpopfio',is)
      if((io.gt.nomax(is)).or.(io.lt.1)) then
         if (Node.eq.0) then
            write(6,*) 'ATMPOPFIO: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'ATMPOPFIO: IOMIN= 1', ' IOMAX= ',nomax(is)
         endif
         call die
      endif
 
      atmpopfio=qtb(io,is) 
      end function atmpopfio

!
!
      FUNCTION CNFIGFIO(IS,IO)
      integer cnfigfio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns the valence-shell configuration in the atomic ground state
C (i.e. the principal quatum number for orbitals of angular momentum l)

C   INTEGER CNFIGFIO: Principal quantum number of the shell to what 
C                     the orbital belongs ( for polarization orbitals
C                     the quantum number corresponds to the shell which
C                     is polarized by the orbital io) 

      integer l, norb, lorb, izeta, ipol,nsm
      integer  indx, nsmorb

C
      call check_is('cnfigfio',is)
      if ((io.gt.nomax(is)).or.(io.lt.1)) then
         if (Node.eq.0) then
            write(6,*) 'CNFIGFIO: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'CNFIGFIO: IOMIN= 1',
     .           ' IOMAX= ',nomax(is)
         endif
         call die
      endif

        norb=0
        indx=0
        do 10 l=0,lmxosave(is)
         do 8 nsm=1,lsemicsave(l,is)+1
          do 5 izeta=1,nzetasave(l,nsm,is)
            norb=norb+(2*l+1)
            indx=indx+1
            if(norb.ge.io) goto 30
 5        continue
 8       continue
10      continue

        indx=0
        do  20 l=0,lmxosave(is)
          do 18 nsm=1,lsemicsave(l,is)+1
            do 15 ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              indx=indx+1
              if(norb.ge.io) goto 40
15          continue
18        continue
20      continue
        if (Node.eq.0) then
           write(6,*) 'CNFIGFIO: ERROR: ORBITAL INDEX IO=',IO
           write(6,*) 'CNFIGFIO: ERROR: NOT FOUND'
        endif
        call die

30      lorb=l
        nsmorb=nsm
        cnfigfio=cnfigfl(is,lorb,nsmorb)
        return

40      lorb=l
        nsmorb=nsm
        cnfigfio=cnfigfl(is,lorb,nsmorb)  
        return

      end function cnfigfio
!
!
      FUNCTION LOFIO (IS,IO)
      integer lofio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns total angular momentum quantum number of a given atomic basis
C   basis orbital or Kleynman-Bylander projector.

C    INTEGER  IO   : Orbital index (within atom)
C                    IO > 0 => Basis orbitals
C                    IO < 0 => Kleynman-Bylander projectors
C************************OUTPUT*****************************************
C   INTEGER LOFIO  : Quantum number L of orbital or KB projector

      integer l, norb, izeta, ipol, nkb, nsm

      call check_is('lofio',is)
      if ((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then 
         if (Node.eq.0) then
            write(6,*) 'LOFIO: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'LOFIO: IOMIN= ',-nkbmax(is),
     .           ' IOMAX= ',nomax(is)
         endif
         CALL DIE
      endif
 
       if (io.gt.0) then

        norb=0
        do 10 l=0,lmxosave(is)
          do 8 nsm=1,lsemicsave(l,is)+1
            do 5 izeta=1,nzetasave(l,nsm,is)
              norb=norb+(2*l+1)
              if(norb.ge.io) goto 30
 5          continue
 8        continue
10      continue

        do  20 l=0,lmxosave(is)
          do 18 nsm=1,lsemicsave(l,is)+1
            do 15 ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              if(norb.ge.io) goto 40
15          continue
18        continue
20      continue
        if (Node.eq.0) then
          write(6,*) 'LOFIO: ERROR: ORBITAL INDEX IO=',IO
          write(6,*) 'LOFIO: ERROR: NOT FOUND'
        endif
        call die

30      lofio=l
        return

40      lofio=l+1
        return

       elseif(io.lt.0) then

        nkb=0
        do 50 l=0,lmxkbsave(is)
          do 45 izeta=1,nkblsave(l,is)
             nkb=nkb-(2*l+1)
             if(nkb.le.io) goto 60
45        continue
50      continue 

60      lofio=l       

c       elseif (io.eq.0) then
        else

        lofio=0

        endif
      end  function lofio
!
!
      FUNCTION MOFIO (IS,IO)
      integer mofio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns magnetic quantum number of a given atomic basis
C   basis orbital or Kleynman-Bylander projector.

C    INTEGER  IO   : Orbital index (within atom)
C                    IO > 0 => Basis orbitals
C                    IO < 0 => Kleynman-Bylander projectors
C************************OUTPUT*****************************************
C   INTEGER MOFIO  : Quantum number M of orbital or KB projector

      integer l, norb, izeta, ipol, nkb, lorb, lkb, nsm

      call check_is('mofio',is)
      if((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
         if (Node.eq.0) then
            write(6,*) 'MOFIO: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'MOFIO: IOMIN= ',-nkbmax(is),
     .           ' IOMAX= ',nomax(is)
         endif
         CALL DIE
      endif

      if (io.gt.0) then

        norb=0
        do 10 l=0,lmxosave(is)
          do 8 nsm=1,lsemicsave(l,is)+1
            do 5 izeta=1,nzetasave(l,nsm,is)
              norb=norb+(2*l+1)
              if(norb.ge.io) goto 30
 5          continue
 8        continue
10      continue 

        do  20 l=0,lmxosave(is)
          do 18 nsm=1,lsemicsave(l,is)+1
            do 15 ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              if(norb.ge.io) goto 40
15          continue
18        continue
20      continue
        if (Node.eq.0) then
          write(6,*) 'MOFIO: ERROR: ORBITAL INDEX IO=',IO
          write(6,*) 'MOFIO: ERROR: NOT FOUND'
        endif
        call die

30      lorb=l 
        mofio=io-norb+lorb
        return

40      lorb=l+1 
        mofio=io-norb+lorb
        return


       elseif(io.lt.0) then


        nkb=0
        do 50 l=0,lmxkbsave(is)
          do 45 izeta=1,nkblsave(l,is)
             nkb=nkb-(2*l+1)
             if(nkb.le.io) goto 60
45        continue
50      continue

60      lkb=l
        mofio=-io+nkb+lkb
c       elseif (io.eq.0) then
        else

        mofio=0

        endif
        
      end function mofio
!
!
      FUNCTION SYMFIO (IS,IO)
      character(len=20) symfio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns a label describing the symmetry of the
C   basis orbital or Kleynman-Bylander projector.
C    INTEGER  IO   : Orbital index (within atom)
C                    IO > 0 => Basis orbitals
C                    IO < 0 => Kleynman-Bylander projectors

C   INTEGER SYMFIO  : Symmetry of the orbital or KB projector
C  2) Returns 's' for IO = 0


      integer ilm, i, lorb, morb
      integer, parameter  :: lmax_sym=3
   
      character(len=6)  sym_label((lmax_sym+1)*(lmax_sym+1)) 
      character(len=7)  paste
         
      external paste
C
        data  sym_label(1)          / 's' /
        data (sym_label(i),i=2,4)   / 'py', 'pz', 'px' /
        data (sym_label(i),i=5,9)   / 'dxy', 'dyz', 'dz2',
     .                                'dxz', 'dx2-y2' / 
        data (sym_label(i),i=10,16) / 'f', 'f', 'f', 'f', 
     .                                'f', 'f', 'f' /

        call check_is('symfio',is)
        if((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
          if (Node.eq.0) then
          write(6,*) 'SYMFIO: THERE ARE NO DATA FOR IO=',IO
          write(6,*) 'SYMFIO: IOMIN= ',-nkbmax(is),
     .     ' IOMAX= ',nomax(is)
          endif
          CALL DIE
        endif

        lorb=lofio(is,io)
        morb=mofio(is,io)

        if(lorb.gt.lmax_sym ) then 
            symfio=' '
        else
            ilm=lorb*lorb+lorb+morb+1  
            if(pol(is,io)) then 
              symfio=paste('P',sym_label(ilm))
            else
              symfio=sym_label(ilm) 
            endif 
        endif         

      end function symfio
!
!  End of FIOs
!


      FUNCTION ATMPOPFL(IS,L,NSM)
      real*8 atmpopfl
      integer, intent(in)  :: is   ! Species index
      integer, intent(in)  :: l    ! Angular momentum of the shell
      integer, intent(in)  :: nsm  ! Index of the shell (if there are
                                   !  semicore states)
      
C    REAL*8 ATMPOPFL: Population of the valence l-shell in the
C                     ground state 

       call check_is('atmpopfl',is)
        if ((nsm.lt.1).or.(nsm.gt.lsemicsave(l,is)+1)) then
          if (Node.eq.0) then
            write(6,*) 'ATMPOPFL: THERE ARE NO DATA FOR NSM=',NSM
            write(6,*) 'ATMPOPFL: NSMMIN= 1, NSMMAX= ',
     .           lsemicsave(l,is)+1
          endif
          CALL DIE
        endif

        if(l.gt.3) then  
            write(6,*) '*** Test for l .gt. 3 ***'
            atmpopfl=0.0d0
        else
            atmpopfl=qltb(l,nsm,is)
        endif

      end function atmpopfl
!
!

      FUNCTION CNFIGFL(IS,L,NSM)
      integer cnfigfl
      integer, intent(in)  :: is   ! Species index
      integer, intent(in)  :: l    ! Angular momentum of the shell
      integer, intent(in)  :: nsm  ! Index of the shell (if there are
                                   !  semicore states)

C Returns the valence-shell configuration in the atomic ground state
C (i.e. the principal quatum number for orbitals of angular momentum l)

C   INTEGER CNFIGFL: Valence-shell configuration in the atomic 
C                     ground state (Principal quantum number).

      call check_is('cnfigfl',is)
        if ((nsm.lt.1).or.(nsm.gt.lsemicsave(l,is)+1)) then
          if (Node.eq.0) then
            write(6,*) 'CNFIGFL: THERE ARE NO DATA FOR NSM=',NSM
            write(6,*) 'CNFIGFL: NSMMIN= 1, NSMMAX= ',
     .           lsemicsave(l,is)+1
          endif
          CALL DIE
        endif

          if(l.gt.3) then 
             write(6,*) '*** Test for l .gt. 3 ***'
             cnfigfl=l+1
          else 
            cnfigfl=cnfigtb(l,nsm,is) 
          endif 

      end function cnfigfl
!
!

      FUNCTION NZTFL (IS,L)
      integer nztfl
      integer, intent(in)  :: is   ! Species index
      integer, intent(in)  :: l    ! Angular momentum of the basis funcs
C Returns the number of different basis functions
C with the same angular momentum and for a given species

      integer nsm

      call check_is('nztfl',is)
         
          if(l.gt.lmxosave(is)+1) then 
            nztfl=0
          elseif(l.eq.lmxosave(is)+1) then 
            nztfl=0
            do nsm=1,lsemicsave(l-1,is)+1
                nztfl=nztfl+npolorbsave(l-1,nsm,is)  
            enddo 
          elseif(l.eq.0) then 
            nztfl=0
            do nsm=1,lsemicsave(l,is)+1
                  nztfl= nztfl+nzetasave(l,nsm,is)
            enddo 
          else
            nztfl=0
            do nsm=1,lsemicsave(l,is)+1
                nztfl=nztfl+nzetasave(l,nsm,is)
            enddo 
            do nsm=1,lsemicsave(l-1,is)+1
                nztfl=nztfl+npolorbsave(l-1,nsm,is)
            enddo 
          endif 

      end function nztfl
!
!

      FUNCTION EPSKB (IS,IO)
      real*8 epskb
      integer, intent(in)   ::  is   ! Species index
      integer, intent(in)   ::  io   ! KB proyector index (within atom)
                                     ! May be positive or negative 
                                     ! (only ABS(IO) is used).

C  Returns the energies epsKB_l of the Kleynman-Bylander projectors:
C       <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                 Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C  where Phi_lm is returned by subroutine PHIATM.


C  REAL*8 EPSKB  : Kleynman-Bylander projector energy
C  Energy in Rydbergs.

      integer  ionew, nkb, indx, l, ikb
C
C******************************************************************
C*****************Variables in the common blocks*******************
C
      call check_is('epskb',is)
 
         ionew=-abs(io)
         if (ionew.eq.0) then 
           if (Node.eq.0) then
             write(6,*) 'EPSKB: FUNCTION CANNOT BE CALLED WITH'
     .         ,' ARGUMENT EQUAL TO ZERO' 
           endif
           CALL DIE 
         endif

         if(ionew.lt.-nkbmax(is)) then
           if (Node.eq.0) then
             write(6,*) 'EPSKB: THERE ARE NO DATA FOR IO=',IONEW
             write(6,*) 'EPSKB: IOMIN= ',-nkbmax(is)
           endif
           CALL DIE
         endif

         nkb=0
         indx=0
         do 10  l=0,lmxkbsave(is)
             do 5 ikb=1,nkblsave(l,is)
                indx=indx+1
                nkb=nkb-(2*l+1)
                if(nkb.le.ionew) goto 20
5            continue
10       continue 

20        indx=-indx
 
        epskb=table(2,indx,is)

      end function epskb
!
!
      FUNCTION NKBL_FUNC (IS,L)
      integer nkbl_func
      integer, intent(in)  :: is   ! Species index
      integer, intent(in)  :: l    ! Angular momentum of the basis funcs

C Returns the number of different KB projectors
C with the same angular momentum and for a given specie

      call check_is('nkbl_func',is)

          if(l.gt.lmxkbsave(is)) then 
            nkbl_func=0
          else
            nkbl_func=nkblsave(l,is)
          endif  
 
      end function nkbl_func
!
!

      FUNCTION NSEMICORE (IS,L)
      integer nsemicore
      integer, intent(in)  :: is   ! Species index
      integer, intent(in)  :: l    ! Angular momentum of the basis funcs

C Returns the number of semicore shells (usually zero)
C with a given angular momentum and for a given specie

      call check_is('nsemicore',is)

      if(l.gt.lmxosave(is)) then 
         nsemicore=0
      else
         nsemicore=lsemicsave(l,is)
      endif

      end function nsemicore
!
!

      FUNCTION POL (IS,IO)
      logical pol
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)
                                   ! io>0 => basis orbitals

C If true, the orbital IO is a perturbative polarization orbital

C  2) Returns zero for IO = 0   ???????

      integer l, norb, izeta, ipol, nsm

      call check_is('pol',is)
        if(io.gt.nomax(is)) then 
          if (Node.eq.0) then
            write(6,*) 'POL: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'POL: IOMIN= 0',
     .       ' IOMAX= ',nomax(is)
          endif
          CALL DIE
        endif
 
        norb=0
        do 10 l=0,lmxosave(is)
         do 8 nsm=1,lsemicsave(l,is)+1
          do 5 izeta=1,nzetasave(l,nsm,is)
            norb=norb+(2*l+1)
            if(norb.ge.io) goto 30
 5        continue
 8       continue
10      continue  

        do  20 l=0,lmxosave(is)
          do 18 nsm=1,lsemicsave(l,is)+1
            do 15 ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1) 
              if(norb.ge.io) goto 40
15          continue
18        continue
20      continue    

        if (Node.eq.0) then
          write(6,*) 'POL: ERROR: ORBITAL INDEX IO=',IO
          write(6,*) 'POL: ERROR: NOT FOUND'
        endif
        call die

30      pol=.false.
        return

40      pol=.true.

      end function pol
!
!

      FUNCTION UION ( IS )
      real*8 uion
      integer, intent(in) :: is    ! Species index

C  Returns electrostatic self-energy of the 'ions', assigned by routine PSEUDO
C   Energy in Rydbergs

      call check_is('uion',is)
      uion=slfe(is)

      end function uion
!
!
      function rcore(is)
      real*8 rcore
      integer, intent(in) :: is    ! Species index

C  Returns cutoff radius of the pseudo-core charge density for the non-linear
C   core corrections for xc potential.
C  Distances in Bohr

      call check_is('rcore',is)
      rcore=coretab(1,1,is)*dble(ntbmax-1)

      end function rcore
!
!
      function rcut(is,io)
      real*8 rcut
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)
                                   ! io> => basis orbitals
                                   ! io<0  => KB projectors
                                   ! io=0 : Local screened pseudopotential

C  Returns cutoff radius of Kleynman-Bylander projectors and
C  atomic basis orbitals.
C  Distances in Bohr

      integer l, norb, lorb, nzetorb, izeta, ipol, nkb,lkb,nsm
      integer  indx, nsmorb        
C
      call check_is('rcut',is)

      if ((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
        if (Node.eq.0) then
          write(6,*) 'RCUT: THERE ARE NO DATA FOR IO=',IO
          write(6,*) 'RCUT: IOMIN= ',-nkbmax(is),
     .      ' IOMAX= ',nomax(is)
        endif
        CALL DIE
      endif


       if (io.gt.0) then

        norb=0 
        indx=0
        do 10 l=0,lmxosave(is)
         do 8 nsm=1,lsemicsave(l,is)+1
          do 5 izeta=1,nzetasave(l,nsm,is)
            norb=norb+(2*l+1)
            indx=indx+1 
            if(norb.ge.io) goto 30 
 5        continue      
 8       continue 
10      continue   

        indx=0
        do  20 l=0,lmxosave(is)
          do 18 nsm=1,lsemicsave(l,is)+1
            do 15 ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)  
              indx=indx+1  
              if(norb.ge.io) goto 40
15          continue
18        continue
20      continue     
        if (Node.eq.0) then
          write(6,*) 'RCUT: ERROR: ORBITAL INDEX IO=',IO
          write(6,*) 'RCUT: ERROR: NOT FOUND'
        endif
        call die

30      lorb=l
        nzetorb=izeta
        nsmorb=nsm
        rcut=rcotb(nzetorb,lorb,nsmorb,is)  
c       rcut=table(1,indx,is)*(ntbmax-1)
        return  

40      lorb=l
        nzetorb=ipol
        nsmorb=nsm
        rcut=rcpoltb(nzetorb,lorb,nsmorb,is) 
c       rcut=tabpol(1,indx,is)*(ntbmax-1)
        return 


       elseif(io.lt.0) then 


        nkb=0
        do 50 l=0,lmxkbsave(is)
          do 45 izeta=1,nkblsave(l,is)
            nkb=nkb-(2*l+1)
            if(nkb.le.io) goto 60
45        continue
50      continue 

60      lkb=l 
        indx=-(lkb+1)
        rcut=rctb(izeta,lkb,is) 
c       rcut=table(1,indx,is)*(ntbmax-1)
        return 

c       elseif (io.eq.0) then
        else

        rcut=table(2,0,is)

        endif

      end function rcut
!
      subroutine vlocal_sub(is,r,v,grv)
      integer, intent(in) :: is      ! Species index
      real*8, intent(in)  :: r(3)    ! Point vector, relative to atom
      real*8, intent(out) :: v(3)    ! Value of local pseudopotential
      real*8, intent(out) :: grv(3)  ! Gradient of local pseudopotential

C Returns local part of neutral-atom Kleynman-Bylander pseudopotential.
C Written by D.Sanchez-Portal. Oct, 1996
C Distances in Bohr
C Energies in Rydbergs
C  2) Returns exactly zero when |R| > RCUT(IS,0)

         integer i
         double precision  rvmx, rmod, dvdr, delt

         call check_is('vlocal_sub',is)
 
          delt=table(1,0,is)
          rvmx=table(2,0,is)
         
          rmod=0.0d0
          do i=1,3
            rmod=rmod+r(i)*r(i)
          enddo
          rmod=dsqrt(rmod)
           
          if(rmod.gt.rvmx) then
             v=0.0d0
             grv(1:3)=0.0d0
          else
            call splint(delt,table(3,0,is),tab2(1,0,is),ntbmax,
     .        rmod,v,dvdr)
             rmod=rmod+1.0d-20
             grv(1:3) = dvdr * r(1:3)/rmod
          endif
 
      end subroutine vlocal_sub
!
      subroutine psover(is1,is2,r,energ,dedr)
      integer, intent(in) :: is1, is2     ! Species indexes
      real*8, intent(in)  :: r       ! Distance between atoms
      real*8, intent(out) :: energ   ! Value of the correction
                                     !  interaction energy
      real*8, intent(out) :: dedr    ! Radial derivative of the correction

C Returns electrostatic correction to the ions interaction energy
C due to the overlap of the two 'local pseudopotential charge densities'
C Written by D.Sanchez-Portal. March, 1997
C Distances in Bohr
C Energies in Rydbergs
C  2) Returns exactly zero when |R| > Rchloc

      integer ismx, ismn, indx
      real*8 dloc, delt, rcmx, r_local

      call check_is('psover',is1)
      call check_is('psover',is2)

          ismx=max(is1,is2)
          ismn=min(is1,is2)
    
          indx=((ismx-1)*ismx)/2 + ismn

          dloc=corrtab(1,2,indx)
          if(dabs(dloc).lt.1.0d-8) then 
            energ=0.0d0 
            dedr=0.0d0 
            return
          endif

          delt=corrtab(1,1,indx) 
          rcmx=delt*(ntbmax-1)
          if(r.gt.rcmx-1.d-12) then

             energ=0.0d0
             dedr=0.0d0
        
          else

            call splint(delt,corrtab(2,1,indx),corrtab(2,2,indx),ntbmax,
     .        r,energ,dedr)
         
             r_local = r+1.0d-20
             energ=2.0d0*energ/r_local
             dedr=(-energ + 2.0d0*dedr)/r_local
          endif

      end subroutine psover

!
      subroutine chcore_sub(is,r,ch,grch)
      integer, intent(in) :: is      ! Species index
      real*8, intent(in)  :: r(3)    ! Point vector, relative to atom
      real*8, intent(out) :: ch      ! Value of pseudo-core charge dens.
      real*8, intent(out) :: grch(3) ! Gradient of pseudo-core ch. dens.

C Returns returns pseudo-core charge density for non-linear core correction
C in the xc potential.
C Written by D.Sanchez-Portal. May 1996

C Distances in Bohr
C Energies in Rydbergs
C Density in electrons/Bohr**3
C  2) Returns exactly zero when |R| > Rcore

      integer i
      real*8 core, delt, rcmx, rmod, dchdr
          
      call check_is('chcore_sub',is)

          core=coretab(1,2,is)
          if(core.eq.0) then 
            ch=0.0d0 
            grch(1:3) = 0.0d0 
            return
          endif

          delt=coretab(1,1,is) 
          rcmx=delt*(ntbmax-1)
          
          rmod=0.0d0
          do i=1,3
            rmod=rmod+r(i)*r(i)
          enddo
          rmod=dsqrt(rmod)
           
          if(rmod.gt.rcmx-1.d-12) then
             ch=0.0d0
             grch(1:3)=0.0d0
          else
            call splint(delt,coretab(2,1,is),coretab(2,2,is),ntbmax,
     .        rmod,ch,dchdr)
             rmod=rmod+1.0d-20
             grch(1:3)=dchdr*r(1:3)/rmod
          endif

      end subroutine chcore_sub
!
      subroutine psch(is,r,ch,grch)
      integer, intent(in) :: is      ! Species index
      real*8, intent(in)  :: r(3)    ! Point vector, relative to atom
      real*8, intent(out) :: ch      ! Local pseudopot. charge dens.
      real*8, intent(out) :: grch(3) ! Gradient of local ps. ch. dens.

C Returns 'local-pseudotential charge density'.
C Written by D.Sanchez-Portal. March, 1997
C Distances in Bohr
C Energies in Rydbergs
C Density in electrons/Bohr**3
C  2) Returns exactly zero when |R| > Rchloc

      integer i
      double precision rmod, dchdr, dloc, delt, rcmx

      call check_is('psch',is)

          dloc=chloctab(1,2,is)
          if(dabs(dloc).lt.1.0d-8) then 
            ch=0.0d0 
            grch(1:3)=0.0d0 
            return
          endif

          delt=chloctab(1,1,is) 
          rcmx=delt*(ntbmax-1)
          rmod=0.0d0
          do i=1,3
            rmod=rmod+r(i)*r(i)
          enddo
          rmod=dsqrt(rmod)

          if(rmod.gt.rcmx-1.d-12) then
             ch=0.0d0
             grch(1:3)=0.0d0
          else
            call splint(delt,chloctab(2,1,is),chloctab(2,2,is),ntbmax,
     .        rmod,ch,dchdr)

             rmod=rmod+1.0d-20
           
             grch(1:3)=dchdr*r(1:3)/rmod
          endif

      end subroutine psch
!
!
!
      subroutine phiatm(is,io,r,phi,grphi)
      integer, intent(in) :: is      ! Species index
      integer, intent(in) :: io      ! Orbital index (within atom)
!              IO > 0 =>  Basis orbitals
!              IO = 0 =>  Local screened pseudopotential
!              IO < 0 =>  Kleynman-Bylander projectors
      real*8, intent(in)  :: r(3)    ! Point vector, relative to atom
      real*8, intent(out) :: phi     ! Basis orbital, KB projector, or
                                     !  local pseudopotential
      real*8, intent(out) :: grphi(3)! Gradient of BO, KB proj, or Loc ps

C  Returns Kleynman-Bylander local pseudopotential, nonlocal projectors,
C  and atomic basis orbitals (and their gradients).
C  Written by D.Sanchez-Portal. May, 1996. 
C  Modified August, 1998.

C Distances in Bohr
C 1) Each projector and basis function has a well defined total
C    angular momentum (quantum number l).
C 2) Basis functions are normalized and mutually orthogonal
C 3) Projection functions are normalized and mutually orthogonal
C 4) Normalization of KB projectors |Phi_lm> is such that 
C     <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                   Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C    where epsKB_l is returned by function EPSKB
C 5) Prints a message and stops when no data exits for IS and/or IO
C 6) Returns exactly zero when |R| > RCUT(IS,IO)
C 7) PHIATM with IO = 0 is strictly equivalent to VLOCAL_SUB

         integer l, norb, lorb, izeta, ipol, nkb,
     .    indx, morb, ilm, i, nsm
         logical pol

         double precision  rly(lmx2), grly(3,lmx2), rmax, rmod,
     .      phir, dphidr, delt
         
        call check_is('phiatm',is)
        if ((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
          if (Node.eq.0) then
            write(6,*) 'PHIATM: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'PHIATM: IOMIN= ',-nkbmax(is),
     .         ' IOMAX= ',nomax(is)
          endif
          CALL DIE
        endif

       pol=.false.
       if (io.gt.0) then

          norb=0 
          indx=0
          do  l=0,lmxosave(is)
           do nsm=1,lsemicsave(l,is)+1
            do izeta=1,nzetasave(l,nsm,is)
               norb=norb+(2*l+1)
               indx=indx+1
               if(norb.ge.io) then 
                   lorb=l
                   morb=io-norb+lorb
                   goto 20 
               endif 
            enddo 
           enddo 
          enddo 

          indx=0
          do  l=0,lmxosave(is)
           do nsm=1,lsemicsave(l,is)+1
            do ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              indx=indx+1
              if(norb.ge.io) then  
                    lorb=l+1
                    morb=io-norb+lorb
                    pol=.true.
                    goto 20 
             endif   
            enddo 
           enddo 
          enddo  

20       continue

         elseif(io.lt.0) then
         nkb=0
         indx=0
         do l=0,lmxkbsave(is)
            do izeta=1,nkblsave(l,is)
              indx=indx+1
              nkb=nkb-(2*l+1)
              if(nkb.le.io) goto 30
            enddo 
         enddo 

 30      lorb=l
         morb=-io+nkb+lorb 
         indx=-indx
      
         else

             indx=0 
C            Next two lines introduced by J.Soler. 2/7/96.
             lorb=0
             morb=0

         endif 


        if(.not.pol) then 
           delt=table(1,indx,is)  
        else
           delt=tabpol(1,indx,is)
        endif 

        rmax=delt*(ntbmax-1)
  
        rmod=0.0d0
        do i=1,3
           rmod=rmod+r(i)*r(i)
        enddo
        rmod=dsqrt(rmod)+1.0d-20

        if(rmod.gt.rmax-1.d-12) then

           phi=0.0d0
           grphi(1:3)=0.0d0
        else
        
           if(.not.pol) then
               call splint(delt,table(3,indx,is),
     .               tab2(1,indx,is),ntbmax,
     .               rmod,phir,dphidr)
           else
               call splint(delt,tabpol(3,indx,is),
     .               tab2pol(1,indx,is),ntbmax,
     .               rmod,phir,dphidr)
           endif

           if (indx.eq.0) then
              phi=phir
              grphi(1:3)=dphidr*r(1:3)/rmod
           else

              ilm = lorb*lorb + lorb + morb + 1
              call rlylm( lorb, r, rly, grly )

              phi = phir * rly(ilm)
              
              do i = 1,3
                grphi(i)=dphidr*rly(ilm)*r(i)/rmod+phir*grly(i,ilm)
              enddo

*             write(6,'(a,i4,2f12.6)')
*    .         'phiatm: ilm,phi/rl,rl*ylm=', ilm, phi, rly(ilm)

           endif
                   
        endif

      end subroutine phiatm


      subroutine all_phi(is,io,r,nphi,phi,grphi)
      integer, intent(in) :: is     ! Species index
      integer, intent(in) :: io     ! Orbital-type switch:
                                    ! IO > 0 => Basis orbitals
                                    ! IO = 0 => Local screened pseudop
                                    ! IO < 0 => KB projectors
      real*8, intent(in)  :: r(3)   ! Point vector, relative to atom
      integer, intent(out):: nphi   ! Number of phi's
      real*8, intent(out) :: phi(:) ! Basis orbital, KB projector, or
                                    !  local pseudopotential
      real*8, optional, intent(out) :: grphi(:,:) ! Gradient of phi

C  Returns Kleynman-Bylander local pseudopotential, nonlocal projectors,
C  and atomic basis orbitals (and their gradients).
C  Same as phiatm but returns all orbitals or KB projectors of the atom
C  Written by D.Sanchez-Portal and J.M.Soler. Jan. 2000 

C Distances in Bohr
C 1) Each projector and basis function has a well defined total
C    angular momentum (quantum number l).
C 2) Basis functions are normalized and mutually orthogonal
C 3) Projection functions are normalized and mutually orthogonal
C 4) Normalization of KB projectors |Phi_lm> is such that 
C     <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                   Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C    where epsKB_l is returned by function EPSKB
C 5) Prints a message and stops when no data exits for IS
C 6) Returns exactly zero when |R| > RCUT(IS,IO)
C 7) PHIATM with IO = 0 is strictly equivalent to VLOCAL_SUB
C 8) If arrays phi or grphi are too small, returns with the required
C    value of nphi

      integer i, index, ipol, it, izeta, jlm,
     .        l, lmax, m, maxlm, maxo, n, nsm
      logical polar
      double precision  rly(lmx2), grly(3,lmx2), rmod,
     .                  phir, dphidr, delt

      integer, parameter :: maxphi=100, maxs=10
      integer, save, dimension(maxphi,maxs,-1:1) :: ilm, ind=0
      logical, save, dimension(maxphi,maxs,-1:1) :: pol
      double precision, save, dimension(maxphi,maxs,-1:1) :: rmax
      logical, dimension(maxphi) :: within

!     Check that species index is valid
      call check_is('all_phi',is)

!     Find number of orbitals
      if (io.gt.0) then
        it=+1
        nphi=nomax(is)
      elseif (io.lt.0) then
        it=-1
        nphi=nkbmax(is)
      else
        it=0
        nphi=1
      endif

!     Find internal indexes of required orbitals
      if (ind(1,is,it).eq.0) then

!       Check size of tables
        if (nphi.gt.maxphi) stop 'all_phi: maxphi too small'
        if (is.gt.maxs)     stop 'all_phi: maxs too small'

        if (io.gt.0) then
          n=0 
          index=0
          do l=0,lmxosave(is)
            do nsm=1,lsemicsave(l,is)+1
              do izeta=1,nzetasave(l,nsm,is)
                index=index+1
                do m=-l,l
                  n=n+1
                  ind(n,is,it)=index
                  ilm(n,is,it)=l*(l+1)+m+1
                  pol(n,is,it)=.false.
                  rmax(n,is,it)=table(1,index,is)*(ntbmax-1)
                  rmax(n,is,it)=rmax(n,is,it)-1.d-12
                enddo 
              enddo 
            enddo 
          enddo 
          index=0
          do l=0,lmxosave(is)
            do nsm=1,lsemicsave(l,is)+1
              do ipol=1, npolorbsave(l,nsm,is)
                index=index+1
                do m=-(l+1),(l+1)
                  n=n+1
                  ind(n,is,it)=index
                  ilm(n,is,it)=(l+1)*(l+2)+m+1
                  pol(n,is,it)=.true.
                  rmax(n,is,it)=tabpol(1,index,is)*(ntbmax-1)
                  rmax(n,is,it)=rmax(n,is,it)-1.d-12
                enddo 
              enddo 
            enddo 
          enddo  

        elseif (io.lt.0) then
          n=0
          index=0
          do l=0,lmxkbsave(is)
            do izeta=1,nkblsave(l,is)
              index=index-1
              do m=-l,l
                n=n+1
                ind(n,is,it)=index
                ilm(n,is,it)=l*(l+1)+m+1
                pol(n,is,it)=.false.
                rmax(n,is,it)=table(1,index,is)*(ntbmax-1)
                rmax(n,is,it)=rmax(n,is,it)-1.d-12
              enddo 
            enddo 
          enddo 

        else
          n=1
          index=0
          ind(n,is,it)=0 
          ilm(n,is,it)=1
          pol(n,is,it)=.false.
          rmax(n,is,it)=table(1,index,is)*(ntbmax-1)
          rmax(n,is,it)=rmax(n,is,it)-1.d-12
        endif 

!       A safety check
        if (n.ne.nphi) stop 'all_phi: n.ne.nphi'

      endif

!     Check size of output arrays
      if (present(grphi)) then
        if (size(grphi,1).ne.3)
     .    stop 'all_phi: incorrect first dimension of grphi'
        n = min( size(phi), size(grphi,2) )
      else
        n = size(phi)
      endif
      if (n.lt.nphi) return

!     Initialize orbital values
      phi(1:nphi) = 0.d0
      if (present(grphi)) grphi(:,1:nphi) = 0.d0

!     Find for which orbitals rmod < rmax
      rmod = sqrt(sum(r*r)) + 1.0d-20
      within(1:nphi) = rmod .lt. rmax(1:nphi,is,it)
      if (.not.any(within(1:nphi))) return

!     Find spherical harmonics
      maxlm = maxval( ilm(1:nphi,is,it), mask=within )
      lmax=nint(sqrt(dble(maxlm)))-1
      call rlylm(lmax,r,rly,grly)

!     Find orbital values
      index=0
      polar=.false.
      i_loop: do i=1,nphi

!       Check if rmod > rmax
        if (.not.within(i)) cycle i_loop
          
!       Find radial part of orbital
        if ( (ind(i,is,it).ne.index) .or.
     .       (pol(i,is,it).neqv.polar) ) then
          index=ind(i,is,it)
          polar=pol(i,is,it)
          if (polar) then
            delt=tabpol(1,index,is)
            call splint(delt,tabpol(3,index,is),
     .                  tab2pol(1,index,is),ntbmax,
     .                  rmod,phir,dphidr)
          else
            delt=table(1,index,is)  
            call splint(delt,table(3,index,is),
     .                  tab2(1,index,is),ntbmax,
     .                  rmod,phir,dphidr)
          endif

        endif

!       Multiply radial and angular parts
        jlm = ilm(i,is,it)
        phi(i) = phir * rly(jlm)
        if (present(grphi))
     .    grphi(:,i) = dphidr * rly(jlm) * r(:) / rmod + 
     .                 phir * grly(:,jlm)

      enddo i_loop

      end subroutine all_phi
!
        subroutine rphiatm(is,io,r,phi,dphidr)
      integer, intent(in) :: is      ! Species index
      integer, intent(in) :: io      ! Orbital index (within atom)
!              IO > 0 =>  Basis orbitals
!              IO = 0 =>  Local screened pseudopotential
!              IO < 0 =>  Kleynman-Bylander projectors
      real*8, intent(in)  :: r       ! Radial distance, relative to atom
      real*8, intent(out) :: phi     ! Basis orbital, KB projector, or
                                     !  local pseudopotential
      real*8, intent(out) :: dphidr  ! Radial derivative of BO, 
                                     !  KB proj, or Loc pseudopot.

C  Returns the radial component of 
C  Kleynman-Bylander local pseudopotential, nonlocal projectors,
C  and atomic basis orbitals (and their radial drivatives)

C Distances in Bohr

C 1) Each projector and basis function has a well defined total
C    angular momentum (quantum number l).
C 2) Basis functions are normalized and mutually orthogonal
C 3) Projection functions are normalized and mutually orthogonal
C 4) Normalization of KB projectors |Phi_lm> is such that 
C     <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                   Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C    where epsKB_l is returned by function EPSKB
C 6) Returns exactly zero when |R| > RCUT(IS,IO)
C 7) RPHIATM with ITYPE = 0 is strictly equivalent to VLOCAL_SUB

         integer l, norb, lorb, izeta, ipol, nkb,
     .    indx, nsm

         double precision  rmax, rmod, phir, delt
         logical pol

        call check_is('rphiatm',is)
        if((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
          if (Node.eq.0) then
            write(6,*) 'RPHIATM: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'RPHIATM: IOMIN= ',-nkbmax(is),
     .       ' IOMAX= ',nomax(is)
          endif
          CALL DIE
        endif

       pol=.false.
       if (io.gt.0) then

          norb=0 
          indx=0
          do  l=0,lmxosave(is)
           do nsm=1,lsemicsave(l,is)+1
            do izeta=1,nzetasave(l,nsm,is)
               norb=norb+(2*l+1)
               indx=indx+1
               if(norb.ge.io) then
                   lorb=l
                   goto 20
               endif
            enddo
           enddo
          enddo

          indx=0
          do  l=0,lmxosave(is)
           do nsm=1,lsemicsave(l,is)+1
            do ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              indx=indx+1
              if(norb.ge.io) then
                    lorb=l+1
                    pol=.true.
                    goto 20
             endif 
            enddo
           enddo
          enddo 

20       continue

         elseif(io.lt.0) then
         nkb=0
         indx=0
         do l=0,lmxkbsave(is)
            do izeta=1,nkblsave(l,is)
               indx=indx+1
               nkb=nkb-(2*l+1)
               if(nkb.le.io) goto 30
            enddo 
         enddo 

 30      lorb=l
         indx=-indx
      
c        elseif (io.eq.0) then
         else

             indx=0 
C            Next two lines introduced by J.Soler. 2/7/96.
             lorb=0
         endif 

        if(.not.pol) then 
           delt=table(1,indx,is)  
        else
           delt=tabpol(1,indx,is)
        endif 

        rmax=delt*(ntbmax-1)
  
        rmod=r+1.0d-20

        if(rmod.gt.rmax-1.d-12) then

           phi=0.0d0
           dphidr=0.0d0
        else
        
           if(.not.pol) then
               call splint(delt,table(3,indx,is),
     .               tab2(1,indx,is),ntbmax,
     .               rmod,phir,dphidr)
           else
               call splint(delt,tabpol(3,indx,is),
     .               tab2pol(1,indx,is),ntbmax,
     .               rmod,phir,dphidr)
           endif

            if (lorb.eq.0) then
               phi=phir
            elseif (lorb.eq.1) then
               phi=phir*r
               dphidr=dphidr*r
               dphidr=dphidr+phir 
            else
               phi=phir*r**lorb 
               dphidr=dphidr*r**lorb
               dphidr=dphidr+lorb*r**(lorb-1)*phir 
            endif
 
         endif
      end subroutine rphiatm 
!
!
      subroutine xphiatm(is,io,r,phi,grphi)
      integer, intent(in) :: is      ! Species index
      integer, intent(in) :: io      ! Orbital index (within atom)
!              IO > 0 =>  Basis orbitals
!              IO = 0 =>  Local screened pseudopotential
!              IO < 0 =>  Kleynman-Bylander projectors
      real*8, intent(in)  :: r(3)    ! Point vector, relative to atom
      real*8, intent(out) :: phi     ! Basis orbital, KB projector, or
                                     !  local pseudopotential, 
                                     !  multiplied by x
      real*8, intent(out) :: grphi(3)! Gradient of BO, KB proj, or Loc ps
                                     !  multiplied by x

C  Returns Kleynman-Bylander local pseudopotential, nonlocal projectors,
C  and atomic basis orbitals multiply by x (and their gradients)
C  Written by D.Sanchez-Portal. May, 1996.
C  Modified August, 1998.

C Distances in Bohr

C 1) Each projector and basis function has a well defined total
C    angular momentum (quantum number l).
C 2) Basis functions are normalized and mutually orthogonal
C 3) Projection functions are normalized and mutually orthogonal
C 4) Normalization of KB projectors |Phi_lm> is such that
C     <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                   Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C    where epsKB_l is returned by function EPSKB
C 6) Returns exactly zero when |R| > RCUT(IS,IO)
C 7) XPHIATM with IO = 0 is strictly equivalent to VLOCAL_SUB

         integer l, norb, lorb, izeta, ipol, nkb,
     .    indx, morb, ilm, i, nsm
         double precision  rly(lmx2), grly(3,lmx2), rmax, rmod,
     .      phir, dphidr, delt
         logical pol

        call check_is('xphiatm',is)
        if((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
          if (Node.eq.0) then
           write(6,*) 'XPHIATM: THERE ARE NO DATA FOR IO=',IO
           write(6,*) 'XPHIATM: IOMIN= ',-nkbmax(is),
     .       ' IOMAX= ',nomax(is)
          endif
          CALL DIE
        endif
 
       pol=.false.
       if (io.gt.0) then

          norb=0
          indx=0
          do  l=0,lmxosave(is)
           do nsm=1,lsemicsave(l,is)+1
            do izeta=1,nzetasave(l,nsm,is)
               norb=norb+(2*l+1)
               indx=indx+1
               if(norb.ge.io) then
                   lorb=l
                   morb=io-norb+lorb
                   goto 20
               endif
            enddo
           enddo
          enddo

          indx=0
          do  l=0,lmxosave(is)
           do nsm=1,lsemicsave(l,is)+1
            do ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              indx=indx+1
              if(norb.ge.io) then
                    lorb=l+1
                    morb=io-norb+lorb
                    pol=.true.
                    goto 20
             endif 
            enddo
           enddo
          enddo 


20       continue

         elseif(io.lt.0) then
         nkb=0
         indx=0
         do l=0,lmxkbsave(is)
            do izeta=1,nkblsave(l,is)
               indx=indx+1
               nkb=nkb-(2*l+1)
               if(nkb.le.io) goto 30
            enddo 
         enddo

 30      lorb=l
         morb=-io+nkb+lorb
         indx=-indx

         else

             indx=0
C            Next two lines introduced by J.Soler. 2/7/96.
             lorb=0
             morb=0

         endif


        if(.not.pol) then
           delt=table(1,indx,is)
        else
           delt=tabpol(1,indx,is)
        endif

        rmax=delt*(ntbmax-1)
 
        rmod=0.0d0
        do i=1,3
           rmod=rmod+r(i)*r(i)
        enddo
        rmod=dsqrt(rmod)+1.0d-20

        if(rmod.gt.rmax-1.d-12) then

           phi=0.0d0
           grphi(1)=0.0d0
           grphi(2)=0.0d0
           grphi(3)=0.0d0

        else

           if(.not.pol) then
               call splint(delt,table(3,indx,is),
     .               tab2(1,indx,is),ntbmax,
     .               rmod,phir,dphidr)
           else
               call splint(delt,tabpol(3,indx,is),
     .               tab2pol(1,indx,is),ntbmax,
     .               rmod,phir,dphidr)
           endif

           if (indx.eq.0) then
              phi=phir
              do i=1,3
                 grphi(i)=dphidr*r(i)/rmod
              enddo
           else

              ilm = lorb*lorb + lorb + morb + 1
              call rlylm( lorb, r, rly, grly )

              phi = phir * rly(ilm)

              do i = 1,3
                grphi(i)=dphidr*rly(ilm)*r(i)/rmod+phir*grly(i,ilm) 
                grphi(i)= r(i)*grphi(i)
              enddo
               grphi(1)=phi +grphi(1)
             
               phi = r(1) * phi 
 
*             write(6,'(a,i4,2f12.6)')
*    .         'xphiatm: ilm,phi/rl,rl*ylm=', ilm, phi, rly(ilm)

           endif
                  
        endif

      end subroutine xphiatm
!
      subroutine yphiatm(is,io,r,phi,grphi)
      integer, intent(in) :: is      ! Species index
      integer, intent(in) :: io      ! Orbital index (within atom)
!              IO > 0 =>  Basis orbitals
!              IO = 0 =>  Local screened pseudopotential
!              IO < 0 =>  Kleynman-Bylander projectors
      real*8, intent(in)  :: r(3)    ! Point vector, relative to atom
      real*8, intent(out) :: phi     ! Basis orbital, KB projector, or
                                     !  local pseudopotential, 
                                     !  multiplied by y
      real*8, intent(out) :: grphi(3)! Gradient of BO, KB proj, or Loc ps
                                     !  multiplied by y

C  Returns Kleynman-Bylander local pseudopotential, nonlocal projectors,
C  and atomic basis orbitals multiply by y (and their gradients)
C  Written by D.Sanchez-Portal. May, 1996.
C  Modified August, 1998.

C Distances in Bohr

C 1) Each projector and basis function has a well defined total
C    angular momentum (quantum number l).
C 2) Basis functions are normalized and mutually orthogonal
C 3) Projection functions are normalized and mutually orthogonal
C 4) Normalization of KB projectors |Phi_lm> is such that
C     <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                   Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C    where epsKB_l is returned by function EPSKB
C 6) Returns exactly zero when |R| > RCUT(IS,IO)
C 7) YPHIATM with IO = 0 is strictly equivalent to VLOCAL_SUB

         integer l, norb, lorb, izeta, ipol, nkb,
     .    indx, morb, ilm, i, nsm
         double precision  rly(lmx2), grly(3,lmx2), rmax, rmod,
     .      phir, dphidr, delt
         logical pol

        call check_is('yphiatm',is)
        if((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
          if (Node.eq.0) then
           write(6,*) 'YPHIATM: THERE ARE NO DATA FOR IO=',IO
           write(6,*) 'YPHIATM: IOMIN= ',-nkbmax(is),
     .       ' IOMAX= ',nomax(is)
          endif
          CALL DIE
        endif
 
       pol=.false.
       if (io.gt.0) then

          norb=0
          indx=0
          do  l=0,lmxosave(is)
           do nsm=1,lsemicsave(l,is)+1
            do izeta=1,nzetasave(l,nsm,is)
               norb=norb+(2*l+1)
               indx=indx+1
               if(norb.ge.io) then
                   lorb=l
                   morb=io-norb+lorb
                   goto 20
               endif
            enddo
           enddo
          enddo

          indx=0
          do  l=0,lmxosave(is)
           do nsm=1,lsemicsave(l,is)+1
            do ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              indx=indx+1
              if(norb.ge.io) then
                    lorb=l+1
                    morb=io-norb+lorb
                    pol=.true.
                    goto 20
             endif 
            enddo
           enddo
          enddo 


20       continue

         elseif(io.lt.0) then
         nkb=0
         indx=0
         do l=0,lmxkbsave(is)
            do izeta=1,nkblsave(l,is)
               indx=indx+1
               nkb=nkb-(2*l+1)
               if(nkb.le.io) goto 30
            enddo 
         enddo

 30      lorb=l
         morb=-io+nkb+lorb
         indx=-indx

         else

             indx=0
C            Next two lines introduced by J.Soler. 2/7/96.
             lorb=0
             morb=0

         endif


        if(.not.pol) then
           delt=table(1,indx,is)
        else
           delt=tabpol(1,indx,is)
        endif

        rmax=delt*(ntbmax-1)
 
        rmod=0.0d0
        do i=1,3
           rmod=rmod+r(i)*r(i)
        enddo
        rmod=dsqrt(rmod)+1.0d-20

        if(rmod.gt.rmax-1.d-12) then

           phi=0.0d0
           grphi(1)=0.0d0
           grphi(2)=0.0d0
           grphi(3)=0.0d0

        else

           if(.not.pol) then
               call splint(delt,table(3,indx,is),
     .               tab2(1,indx,is),ntbmax,
     .               rmod,phir,dphidr)
           else
               call splint(delt,tabpol(3,indx,is),
     .               tab2pol(1,indx,is),ntbmax,
     .               rmod,phir,dphidr)
           endif

           if (indx.eq.0) then
              phi=phir
              do i=1,3
                 grphi(i)=dphidr*r(i)/rmod
              enddo
           else

              ilm = lorb*lorb + lorb + morb + 1
              call rlylm( lorb, r, rly, grly )

              phi = phir * rly(ilm)

              do i = 1,3
                grphi(i)=dphidr*rly(ilm)*r(i)/rmod+phir*grly(i,ilm) 
                grphi(i)= r(i)*grphi(i)
              enddo
               grphi(2)=phi +grphi(2)
             
               phi = r(2) * phi 
 
*             write(6,'(a,i4,2f12.6)')
*    .         'yphiatm: ilm,phi/rl,rl*ylm=', ilm, phi, rly(ilm)

           endif
                  
        endif

      end subroutine yphiatm

!
      subroutine zphiatm(is,io,r,phi,grphi)
      integer, intent(in) :: is      ! Species index
      integer, intent(in) :: io      ! Orbital index (within atom)
!              IO > 0 =>  Basis orbitals
!              IO = 0 =>  Local screened pseudopotential
!              IO < 0 =>  Kleynman-Bylander projectors
      real*8, intent(in)  :: r(3)    ! Point vector, relative to atom
      real*8, intent(out) :: phi     ! Basis orbital, KB projector, or
                                     !  local pseudopotential, 
                                     !  multiplied by z
      real*8, intent(out) :: grphi(3)! Gradient of BO, KB proj, or Loc ps
                                     !  multiplied by z

C  Returns Kleynman-Bylander local pseudopotential, nonlocal projectors,
C  and atomic basis orbitals multiply by y (and their gradients)
C  Written by D.Sanchez-Portal. May, 1996.
C  Modified August, 1998.

C Distances in Bohr

C 1) Each projector and basis function has a well defined total
C    angular momentum (quantum number l).
C 2) Basis functions are normalized and mutually orthogonal
C 3) Projection functions are normalized and mutually orthogonal
C 4) Normalization of KB projectors |Phi_lm> is such that
C     <Psi|V_KB|Psi'> = <Psi|V_local|Psi'> +
C                   Sum_lm( epsKB_l * <Psi|Phi_lm> * <Phi_lm|Psi'> )
C    where epsKB_l is returned by function EPSKB
C 6) Returns exactly zero when |R| > RCUT(IS,IO)
C 7) ZPHIATM with IO = 0 is strictly equivalent to VLOCAL_SUB

         integer l, norb, lorb, izeta, ipol, nkb,
     .    indx, morb, ilm, i, nsm
         double precision  rly(lmx2), grly(3,lmx2), rmax, rmod,
     .      phir, dphidr, delt
         logical pol

        call check_is('zphiatm',is)
        if((io.gt.nomax(is)).or.(io.lt.-nkbmax(is))) then
          if (Node.eq.0) then
           write(6,*) 'ZPHIATM: THERE ARE NO DATA FOR IO=',IO
           write(6,*) 'ZPHIATM: IOMIN= ',-nkbmax(is),
     .       ' IOMAX= ',nomax(is)
          endif
          CALL DIE
        endif
 
       pol=.false.
       if (io.gt.0) then

          norb=0
          indx=0
          do  l=0,lmxosave(is)
           do nsm=1,lsemicsave(l,is)+1
            do izeta=1,nzetasave(l,nsm,is)
               norb=norb+(2*l+1)
               indx=indx+1
               if(norb.ge.io) then
                   lorb=l
                   morb=io-norb+lorb
                   goto 20
               endif
            enddo
           enddo
          enddo

          indx=0
          do  l=0,lmxosave(is)
           do nsm=1,lsemicsave(l,is)+1
            do ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              indx=indx+1
              if(norb.ge.io) then
                    lorb=l+1
                    morb=io-norb+lorb
                    pol=.true.
                    goto 20
             endif 
            enddo
           enddo
          enddo 


20       continue

         elseif(io.lt.0) then
         nkb=0
         indx=0
         do l=0,lmxkbsave(is)
            do izeta=1,nkblsave(l,is)
               indx=indx+1
               nkb=nkb-(2*l+1)
               if(nkb.le.io) goto 30
            enddo 
         enddo

 30      lorb=l
         morb=-io+nkb+lorb
         indx=-indx

         else

             indx=0
C            Next two lines introduced by J.Soler. 2/7/96.
             lorb=0
             morb=0

         endif


        if(.not.pol) then
           delt=table(1,indx,is)
        else
           delt=tabpol(1,indx,is)
        endif

        rmax=delt*(ntbmax-1)
 
        rmod=0.0d0
        do i=1,3
           rmod=rmod+r(i)*r(i)
        enddo
        rmod=dsqrt(rmod)+1.0d-20

        if(rmod.gt.rmax-1.d-12) then

           phi=0.0d0
           grphi(1)=0.0d0
           grphi(2)=0.0d0
           grphi(3)=0.0d0

        else

           if(.not.pol) then
               call splint(delt,table(3,indx,is),
     .               tab2(1,indx,is),ntbmax,
     .               rmod,phir,dphidr)
           else
               call splint(delt,tabpol(3,indx,is),
     .               tab2pol(1,indx,is),ntbmax,
     .               rmod,phir,dphidr)
           endif

           if (indx.eq.0) then
              phi=phir
              do i=1,3
                 grphi(i)=dphidr*r(i)/rmod
              enddo
           else

              ilm = lorb*lorb + lorb + morb + 1
              call rlylm( lorb, r, rly, grly )

              phi = phir * rly(ilm)

              do i = 1,3
                grphi(i)=dphidr*rly(ilm)*r(i)/rmod+phir*grly(i,ilm) 
                grphi(i)= r(i)*grphi(i)
              enddo
               grphi(3)=phi +grphi(3)
             
               phi = r(3) * phi 
 
*             write(6,'(a,i4,2f12.6)')
*    .         'zphiatm: ilm,phi/rl,rl*ylm=', ilm, phi, rly(ilm)

           endif
                  
        endif

      end subroutine zphiatm

      end module atmfuncs
