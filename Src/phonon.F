! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module phonon
!
!     Implements the interface to K. Parlinski's PHONON program
!
!

#ifdef MPI
      use mpi_siesta
#endif
      use parallel,     only : Node
      use precision,    only : dp
      use alloc,        only : re_alloc, de_alloc
      use sys

      implicit none

      private

      public phonon_setup, phonon_set_coords, phonon_restore_coords
      public phonon_write_forces
      public phonon_num_disps

      integer, parameter  :: max_phonon = 100
      integer phonon_num_disps
      integer phonon_atnum(max_phonon)
      character(len=1) phonon_label(max_phonon)
      real(dp) phonon_disp(max_phonon,3)


      CONTAINS

!     Read atomic displacements
!
      subroutine phonon_setup()

      use fdf
      implicit none

      integer ip, idum
      real(dp) xdum(3)
      character(len=1) lab

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      if (fdf_block('MD.ATforPhonon',bfdf)) then
        if (Node.eq.0) then
          write(6,'(/,a)') 'Displacements for PHONON calculation:'
        endif

        ip = 0
        do while (fdf_bline(bfdf,pline))
          xdum(1) = fdf_breals(pline,1)
          xdum(2) = fdf_breals(pline,2)
          xdum(3) = fdf_breals(pline,3)
          lab = fdf_bnames(pline,1)
          idum = fdf_bintegers(pline,1)

          ip = ip + 1
          phonon_disp(ip,1:3) = xdum(1:3)
          phonon_label(ip)    = lab
          phonon_atnum(ip)    = idum
          if (Node.eq.0) then
            write(6,'(i3,1x,a1,1x,3f10.6)') idum, lab, xdum(1:3)
          endif
        enddo

        phonon_num_disps = ip
      else
        call die('Need block MD.ATforPhonon for Phonon MD option')
      endif

      end subroutine phonon_setup

      subroutine phonon_set_coords(istep,xa,ucell)
      integer, intent(in)      :: istep
      real(dp), intent(inout)  :: xa(:,:)
      real(dp), intent(in)     :: ucell(3,3)

      integer iadispl, ix

      if (node.eq.0)
     $    write(6,'(a,i6)') 'siesta: Begin PHONON step = ',istep
      iadispl = phonon_atnum(istep)

      if (node.eq.0) then
         write(6,'(a,i6)') 'siesta:        displace atom = ',
     .        iadispl
         write(6,'(a,3f10.6,a)') 'siesta:    by  = ',
     $        phonon_disp(istep,1:3), ' in fractional coords'
      endif

      do ix = 1,3
         xa(ix,iadispl) = xa(ix,iadispl) +
     $        ucell(ix,1) * phonon_disp(istep,1) +
     .        ucell(ix,2) * phonon_disp(istep,2) +
     .        ucell(ix,3) * phonon_disp(istep,3)
      enddo

      if (node.eq.0)
     $      write(6,'(a,3f10.6,a)') 'siesta:    New coords = ',
     $             xa(1:3,iadispl), ' (cartesian Bohr)'


      end subroutine phonon_set_coords

      subroutine phonon_restore_coords(istep,xa,ucell)
      integer, intent(in)      :: istep
      real(dp), intent(inout)  :: xa(:,:)
      real(dp), intent(in)     :: ucell(3,3)

      integer iadispl, ix

      iadispl = phonon_atnum(istep)
      do ix = 1,3
         xa(ix,iadispl) = xa(ix,iadispl) -
     $        ucell(ix,1) * phonon_disp(istep,1) -
     .        ucell(ix,2) * phonon_disp(istep,2) -
     .        ucell(ix,3) * phonon_disp(istep,3)
      enddo

      end subroutine phonon_restore_coords

      subroutine phonon_write_forces(fa, na, ucell, istep)

C Writes forces in PHONON format to file
C Input forces are in Ry/Bohr
C Forces written to file are in  eV/Ang

      use fdf
      use files, only : slabel, label_length

      implicit          none

      integer           na            ! Number of atoms
      real(dp)          fa(3,na)      ! Forces in Ry/Bohr
      real(dp)          ucell(3,3)    ! Unit cell vectors in Bohrs
      integer           istep         ! PHONON step


C integer atnum               : number in list of displaced atom 
C real*8 disp(3)              : fractional displacement

      integer           atnum, ns
      real(dp)          disp(3)


      real(dp),         pointer, save ::   xa(:,:)
      integer,          pointer, save ::   isa(:)
      character(len=1), pointer, save ::   phonon_label(:)

      external          io_assign, io_close, paste, timer,
     .                  memory

C Internal variables and arrays
      character(len=label_length+7)       :: paste
      character(len=label_length+7), save :: fname
      character(len=40)                   :: sname
      character(len=2)                    :: atname
      logical,                       save :: frstme = .true.
      integer                             :: i, j, ix
      integer                             :: unit1, ia, isp
      real(dp)                            :: Ang, eV
      character(len=50)                   :: acf

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

C Define conversion factors
      Ang = 1.d0 / 0.529177d0
      eV  = 1.d0 / 13.60580d0

      atnum = phonon_atnum(istep)
      disp(1:3) = phonon_disp(istep,1:3)
      ns = maxval(isa(1:na))          ! Compute number of species

C Find file name
      if (frstme) then
         nullify(xa,isa,phonon_label)
         call re_alloc( xa, 1, 3, 1, na, 'xa', 'phonon_write_forces' )
         call re_alloc( isa, 1, na, 'isa', 'phonon_write_forces' )
         call re_alloc( phonon_label, 1, ns, 'phonon_label',
     &                  'phonon_write_forces' )

        acf = fdf_string('AtomicCoordinatesFormat',' ')
        if (.not. (leqi(acf,'Fractional') .or.
     .             leqi(acf,'ScaledByLatticeVectors'))) then
          call die("phonon_write_forces:" //
     .             "Cannot deal with non-fractional coords")
        endif

        if ( fdf_block('AtomicCoordinatesAndAtomicSpecies',bfdf) ) then
          do ia = 1, na
            if (.not. fdf_bline(bfdf,pline))
     .        call die('phonon_write_forces: ERROR in ' //
     .                 'AtomicCoordinatesAndAtomicSpecies block')
            xa(1,ia) = fdf_breals(pline,1)
            xa(2,ia) = fdf_breals(pline,2)
            xa(3,ia) = fdf_breals(pline,3)
            isa(ia) = fdf_bintegers(pline,1)
          enddo
        else
          call die("phonon_write_forces:" //
     .             "Cannot read fractional atomic coords")
        endif
        if ( fdf_block('PhononLabels',bfdf) )  then
          do ia = 1, ns
            if (.not. fdf_bline(bfdf,pline))
     .        call die('phonon_write_forces: ERROR in ' //
     .                 'PhononLabels block')
            isp = fdf_bintegers(pline,1)
            phonon_label(isp) = fdf_bnames(pline,1)
            atname = fdf_bnames(pline,2)
          enddo
        else
          call die("phonon_write_forces: Cannot read PhononLabels")
        endif
        fname = paste(slabel,'.PHONON')
        sname = fdf_string('SystemName','siesta')
        frstme = .false.
        call io_assign(unit1)
        open( unit1, file=fname, position='rewind' )
        write(unit1,'(2a)') '# PHONON file from SIESTA: ', sname
        write(unit1,'(3f14.8)') ((ucell(i,j)/ang,j=1,3),i=1,3)
        call io_close(unit1)
      endif

      call io_assign(unit1)
      open( unit1, file=fname, position='append' )
      write(unit1,'(3f10.6,4x,3f10.6)') (xa(ix,atnum),ix=1,3),
     .     (disp(ix), ix=1,3)
      do i=1,na
        write(unit1,'(i3,1x,a2,3f10.6,3x,3f12.6)')
     .        i, phonon_label(isa(i)), (xa(ix,i),ix=1,3),
     .        (fa(ix,i)* Ang/eV, ix=1,3)
      enddo

      call io_close(unit1)

      end subroutine phonon_write_forces

      end module phonon
