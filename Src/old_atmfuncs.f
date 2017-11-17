! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module old_atmfuncs

C     This module contains a set of routines which provide all the information
C     about the basis set, pseudopotential, atomic mass, etc... of all the
C     chemical species present in the calculation.

!     AG: It is now a "legacy" module to interface the old and new data
!     structures.

      use precision, only: dp
      use sys
      use atmparams, only: nzetmx, lmaxd, nsemx
      use atmparams, only: maxos, nkbmx, ntbmax
      use alloc,     only: re_alloc, de_alloc
      implicit none 

      integer,  save, public     ::  nsmax            

      integer,  public, pointer  ::  izsave(:)          
      integer,  public, pointer  ::  nomax(:)           

      integer,  public, pointer  ::  lmxosave(:)        
      integer,  public, pointer  ::  npolorbsave(:,:,:) 
      integer,  public, pointer  ::  nsemicsave(:,:)    
      integer,  public, pointer  ::  nzetasave(:,:,:)   
      integer,  public, pointer  ::  cnfigtb(:,:,:)     

      logical,  public, pointer  ::  semicsave(:)       
     
      real(dp), public, pointer  :: zvaltb(:)
      real(dp), public, pointer  :: smasstb(:)
      real(dp), public, pointer  :: chargesave(:)
      real(dp), public, pointer  :: slfe(:)
      real(dp), public, pointer  :: lambdatb(:,:,:,:)
      real(dp), public, pointer  :: filtercuttb(:,:,:)

      real(dp), public, pointer  :: qtb(:,:)

      real(dp), public, pointer  ::  rcotb(:,:,:,:)
      real(dp), public, pointer  ::  rcpoltb(:,:,:,:)

      character(len=20), save, public, pointer :: label_save(:)
      character(len=10), save, public, pointer :: basistype_save(:)  

!     Public routines
      public :: labelfis, izofis, zvalfis
      public :: massfis, lomaxfis, nofis
      public :: cnfigfio, lofio, mofio
      public :: clear_tables, allocate_old_arrays
      public :: deallocate_old_arrays


      PRIVATE

      CONTAINS !================================================
!
      subroutine allocate_old_arrays()

      !allocate(rcotb(nzetmx,0:lmaxd,nsemx,nsmax))
      nullify( rcotb )
      call re_alloc( rcotb, 1, nzetmx, 0, lmaxd, 1, nsemx, 1, nsmax,
     &               'rcotb', 'old_atmfuncs' )

      !allocate(rcpoltb(nzetmx,0:lmaxd,nsemx,nsmax))
      nullify( rcpoltb )
      call re_alloc( rcpoltb, 1, nzetmx, 0, lmaxd, 1, nsemx, 1, nsmax,
     &               'rcpoltb', 'old_atmfuncs' )
      !allocate(lambdatb(nzetmx,0:lmaxd,nsemx,nsmax))
      nullify( lambdatb )
      call re_alloc( lambdatb, 1, nzetmx, 0, lmaxd, 1, nsemx,
     &               1, nsmax, 'lambdatb', 'old_atmfuncs' )
      !allocate(filtercuttb(0:lmaxd,nsemx,nsmax))
      nullify( filtercuttb )
      call re_alloc( filtercuttb, 0, lmaxd, 1, nsemx,
     &               1, nsmax, 'filtercuttb', 'old_atmfuncs' )
      !allocate(qtb(maxos,nsmax))
      nullify( qtb )
      call re_alloc( qtb, 1, maxos, 1, nsmax,
     &               'qtb', 'old_atmfuncs' )
      !allocate(slfe(nsmax))
      nullify( slfe )
      call re_alloc( slfe, 1, nsmax, 'slfe', 'old_atmfuncs' )

      !allocate(smasstb(nsmax))
      nullify( smasstb )
      call re_alloc( smasstb, 1, nsmax, 'smasstb', 'old_atmfuncs' )
      !allocate(chargesave(nsmax))
      nullify( chargesave )
      call re_alloc( chargesave, 1, nsmax,
     &               'chargesave', 'old_atmfuncs' )
!

      !allocate(izsave(nsmax))
      nullify( izsave )
      call re_alloc( izsave, 1, nsmax, 'izsave', 'old_atmfuncs' )
      !allocate(lmxosave(nsmax))
      nullify( lmxosave )
      call re_alloc( lmxosave, 1, nsmax, 'lmxosave', 'old_atmfuncs' )
      !allocate(npolorbsave(0:lmaxd,nsemx,nsmax))
      nullify( npolorbsave )
      call re_alloc( npolorbsave, 0, lmaxd, 1, nsemx, 1, nsmax,
     &               'npolorbsave', 'old_atmfuncs' )
      !allocate(nsemicsave(0:lmaxd,nsmax))
      nullify( nsemicsave )
      call re_alloc( nsemicsave, 0, lmaxd, 1, nsmax,
     &               'nsemicsave', 'old_atmfuncs' )
      !allocate(nzetasave(0:lmaxd,nsemx,nsmax))
      nullify( nzetasave )
      call re_alloc( nzetasave, 0, lmaxd, 1, nsemx, 1, nsmax,
     &               'nzetasave', 'old_atmfuncs' )
      !allocate(nomax(nsmax))
      nullify( nomax )
      call re_alloc( nomax, 1, nsmax, 'nomax', 'old_atmfuncs' )

      !allocate(zvaltb(nsmax))
      nullify( zvaltb )
      call re_alloc( zvaltb, 1, nsmax, 'zvaltb', 'old_atmfuncs' )
      !allocate(cnfigtb(0:lmaxd,nsemx,nsmax))
      nullify( cnfigtb )
      call re_alloc( cnfigtb, 0, lmaxd, 1, nsemx, 1, nsmax,
     &               'cnfigtb', 'old_atmfuncs' )
!
      nullify (label_save)
      allocate(label_save(nsmax))
!      call re_alloc(label_save,1,nsmax,"label_save",
!     $                routine= "allocate_old_arrays")
      nullify (basistype_save)
      allocate(basistype_save(nsmax))
!      call re_alloc(basistype_save,1,nsmax,"basistype_save",
!     $                routine= "allocate_old_arrays")
      nullify (semicsave)
      call re_alloc( semicsave, 1, nsmax,
     &               'semicsave', 'old_atmfuncs' )

      end subroutine allocate_old_arrays

      subroutine deallocate_old_arrays()

      call de_alloc( rcotb,       'rcotb',       'old_atmfuncs' )
      call de_alloc( rcpoltb,     'rcpoltb',     'old_atmfuncs' )
      call de_alloc( lambdatb,    'lambdatb',    'old_atmfuncs' )
      call de_alloc( filtercuttb, 'filtercuttb', 'old_atmfuncs' )
      call de_alloc( qtb,         'qtb',         'old_atmfuncs' )
      call de_alloc( slfe,        'slfe',        'old_atmfuncs' )
      call de_alloc( smasstb,     'smasstb',     'old_atmfuncs' )
      call de_alloc( chargesave,  'chargesave',  'old_atmfuncs' )
      call de_alloc( izsave,      'izsave',      'old_atmfuncs' )
      call de_alloc( lmxosave,    'lmxosave',    'old_atmfuncs' )
      call de_alloc( npolorbsave, 'npolorbsave', 'old_atmfuncs' )
      call de_alloc( nsemicsave,  'nsemicsave',  'old_atmfuncs' )
      call de_alloc( nzetasave,   'nzetasave',   'old_atmfuncs' )
      call de_alloc( nomax,       'nomax',       'old_atmfuncs' )
      call de_alloc( zvaltb,      'zvaltb',      'old_atmfuncs' )
      call de_alloc( cnfigtb,     'cnfigtb',     'old_atmfuncs' )
      call de_alloc( semicsave,   'semicsave',   'old_atmfuncs' )
      deallocate( label_save )
!      call de_alloc( label_save, 'label_save', 'old_atmfuncs' )
      deallocate( basistype_save )
!      call de_alloc( basistype_save, 'basistype_save', 'old_atmfuncs' )

      end subroutine deallocate_old_arrays

      subroutine clear_tables()

      integer is

      do is=1,nsmax
        izsave(is)=0
        lmxosave(is)=0
        label_save(is)='  '
        nomax(is)=0  
        semicsave(is)=.false.
              
        nsemicsave(:,is) = 0
        nzetasave(:,:,is) = 0
        rcotb(:,:,:,is) = 0.0_dp
        lambdatb(:,:,:,is) = 0.0_dp
        filtercuttb(:,:,is) = 0.0_dp
        rcpoltb(:,:,:,is) = 0.0_dp

        qtb(1:maxos,is)=0.00_dp

      enddo 
      end subroutine clear_tables

      subroutine check_is(name,is)
      character(len=*), intent(in) :: name
      integer, intent(in) :: is

      if ((is.lt.1).or.(is.gt.nsmax)) then 
            write(6,*) trim(name),': THERE ARE NO DATA FOR IS=',IS
            write(6,*) trim(name),': ISMIN= 1, NSMAX= ',nsmax
         call die()
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
      FUNCTION ZVALFIS( IS )
      real(dp) :: zvalfis          ! Valence charge
      integer, intent(in) :: is            ! Species index

      call check_is('zvalfis',is)
 
      zvalfis=zvaltb(is)
      end function zvalfis
!
      FUNCTION LABELFIS (IS)
      character(len=20) ::  labelfis  ! Atomic label
      integer, intent(in) :: is            ! Species index

      call check_is('labelfis',is)
      labelfis=label_save(is)
      end function labelfis
!
      FUNCTION LOMAXFIS (IS)
      integer :: lomaxfis  ! Maximum ang mom of the Basis Functions
      integer, intent(in) :: is            ! Species index

      integer lmx, nsm

      call check_is('lomaxfis',is)

      lomaxfis=0           
      lmx=lmxosave(is)
      do nsm=1,nsemicsave(lmx,is)+1
         if(npolorbsave(lmx,nsm,is).gt.0)   lomaxfis=lmx+1
      enddo     
      
      lomaxfis=max(lomaxfis,lmx)
      end function lomaxfis
!

      FUNCTION MASSFIS(IS)
      real(dp) :: massfis            ! Mass
      integer, intent(in) :: is            ! Species index

      call check_is('massfis',is)
      massfis=smasstb(is)
      end function massfis
!
      FUNCTION NOFIS(IS)
      integer :: nofis    ! Total number of Basis functions
      integer, intent(in) :: is            ! Species index

      call check_is('nofis',is)
      nofis=nomax(is)
      end function nofis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! AMENoFIS
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
            write(6,*) 'CNFIGFIO: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'CNFIGFIO: IOMIN= 1',
     .           ' IOMAX= ',nomax(is)
         call die()
      endif

        norb=0
        indx=0
        do 10 l=0,lmxosave(is)
         do 8 nsm=1,nsemicsave(l,is)+1
          do 5 izeta=1,nzetasave(l,nsm,is)
            norb=norb+(2*l+1)
            indx=indx+1
            if(norb.ge.io) goto 30
 5        continue
 8       continue
10      continue

        indx=0
        do  20 l=0,lmxosave(is)
          do 18 nsm=1,nsemicsave(l,is)+1
            do 15 ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              indx=indx+1
              if(norb.ge.io) goto 40
15          continue
18        continue
20      continue
           write(6,*) 'CNFIGFIO: ERROR: ORBITAL INDEX IO=',IO
           write(6,*) 'CNFIGFIO: ERROR: NOT FOUND'
        call die()

30      lorb=l
        nsmorb=nsm
        cnfigfio=cnfigtb(lorb,nsmorb,is)
        return

40      lorb=l
        nsmorb=nsm
        cnfigfio=cnfigtb(lorb,nsmorb,is)  
        return

      end function cnfigfio
!
!
      FUNCTION LOFIO (IS,IO)
      integer lofio
      integer, intent(in) :: is    ! Species index
      integer, intent(in) :: io    ! Orbital index (within atom)

C Returns total angular momentum quantum number of a given atomic basis
C orbital

C    INTEGER  IO   : Orbital index (within atom)
C                    IO > 0 => Basis orbitals
C************************OUTPUT*****************************************
C   INTEGER LOFIO  : Quantum number L of orbital or KB projector

      integer l, norb, izeta, ipol, nkb, nsm

      call check_is('lofio',is)
      if ((io.gt.nomax(is)).or.(io.lt.1)) then 
            write(6,*) 'LOFIO: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'LOFIO: IOMIN= ',1,
     .           ' IOMAX= ',nomax(is)
         CALL DIE
      endif
 
        norb=0
        do 10 l=0,lmxosave(is)
          do 8 nsm=1,nsemicsave(l,is)+1
            do 5 izeta=1,nzetasave(l,nsm,is)
              norb=norb+(2*l+1)
              if(norb.ge.io) goto 30
 5          continue
 8        continue
10      continue

        do  20 l=0,lmxosave(is)
          do 18 nsm=1,nsemicsave(l,is)+1
            do 15 ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              if(norb.ge.io) goto 40
15          continue
18        continue
20      continue
          write(6,*) 'LOFIO: ERROR: ORBITAL INDEX IO=',IO
          write(6,*) 'LOFIO: ERROR: NOT FOUND'
        call die

30      lofio=l
        return

40      lofio=l+1
        return

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
C************************OUTPUT*****************************************
C   INTEGER MOFIO  : Quantum number M of orbital

      integer l, norb, izeta, ipol, nkb, lorb, lkb, nsm

      call check_is('mofio',is)
      if((io.gt.nomax(is)).or.(io.lt.1)) then
            write(6,*) 'MOFIO: THERE ARE NO DATA FOR IO=',IO
            write(6,*) 'MOFIO: IOMIN= ',1,
     .           ' IOMAX= ',nomax(is)
         CALL DIE
      endif

        norb=0
        do 10 l=0,lmxosave(is)
          do 8 nsm=1,nsemicsave(l,is)+1
            do 5 izeta=1,nzetasave(l,nsm,is)
              norb=norb+(2*l+1)
              if(norb.ge.io) goto 30
 5          continue
 8        continue
10      continue 

        do  20 l=0,lmxosave(is)
          do 18 nsm=1,nsemicsave(l,is)+1
            do 15 ipol=1, npolorbsave(l,nsm,is)
              norb=norb+(2*(l+1)+1)
              if(norb.ge.io) goto 40
15          continue
18        continue
20      continue
        write(6,*) 'MOFIO: ERROR: ORBITAL INDEX IO=',IO
        write(6,*) 'MOFIO: ERROR: NOT FOUND'
        call die()

30      lorb=l 
        mofio=io-norb+lorb
        return

40      lorb=l+1 
        mofio=io-norb+lorb
        return
        
      end function mofio
!

!  End of FIOs 
!

      end module old_atmfuncs
