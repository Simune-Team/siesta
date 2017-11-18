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
      use atmparams, only: nzetmx, lmaxd, nsemx
      use atmparams, only: maxos
      use alloc,     only: re_alloc, de_alloc
      implicit none 

      integer,  save, public     ::  nsmax            

      integer,  public, pointer  ::  izsave(:)          

      integer,  public, pointer  ::  lmxosave(:)        
      integer,  public, pointer  ::  npolorbsave(:,:,:) 
      integer,  public, pointer  ::  nsemicsave(:,:)    
      integer,  public, pointer  ::  nzetasave(:,:,:)   
      integer,  public, pointer  ::  cnfigtb(:,:,:)     

      logical,  public, pointer  ::  semicsave(:)       
     
      real(dp), public, pointer  :: zvaltb(:)
      real(dp), public, pointer  :: chargesave(:)
      real(dp), public, pointer  :: lambdatb(:,:,:,:)
      real(dp), public, pointer  :: filtercuttb(:,:,:)

      real(dp), public, pointer  :: qtb(:,:)

      real(dp), public, pointer  ::  rcotb(:,:,:,:)
      real(dp), public, pointer  ::  rcpoltb(:,:,:,:)

      character(len=20), save, public, pointer :: label_save(:)
      character(len=10), save, public, pointer :: basistype_save(:)  

!     Public routines
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
      call de_alloc( chargesave,  'chargesave',  'old_atmfuncs' )
      call de_alloc( izsave,      'izsave',      'old_atmfuncs' )
      call de_alloc( lmxosave,    'lmxosave',    'old_atmfuncs' )
      call de_alloc( npolorbsave, 'npolorbsave', 'old_atmfuncs' )
      call de_alloc( nsemicsave,  'nsemicsave',  'old_atmfuncs' )
      call de_alloc( nzetasave,   'nzetasave',   'old_atmfuncs' )
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

      end module old_atmfuncs
