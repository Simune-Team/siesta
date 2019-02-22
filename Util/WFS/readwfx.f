! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

        program readwfx

c****************************************************************************
c READWF  version 0.0.1
c
c This program READWF reads a the coefficients of wavefunctions
c on the expansion of the atomic orbitals basis set, as written
c by SIESTA (unformatted) and writes them in ascii, user-friendly
c form.
c
c Written by P. Ordejon, June 2003
c****************************************************************************
c USAGE:
c
c This program reads files generated from SIESTA, with information of
c the wavefunctions coefficients (NEW STYLE, SystemLabel.WFSX) and writes
c this information in user-friendly form to an output file.
c
!     ** TO BE UPDATED AND MODERNIZED
        
c The program needs two input files:
c
c 1) Main input file, read by standard input. A sample of input file is:
c
c    --- begin input file ---
c        h2o.WFSX
c        h2o.ascii
c        0.00001
c    --- end input file ---
c
c    where:
c    - The first line is the name of the wavefunctions file genterated
c      by Siesta.
c    - The second line is the name of the output file
c    - The third line is a real number, specifying a thresshold to plot
c      the coefficients of the wavefunctions. If the norm of the weight
c      of a wavefunction on a given orbital is smaller than this thresshold
c      then it is not printed. This is useful for large systems, with
c      a very large number of basis orbitals.
c
c 2) The file with the information genetared by SIESTA, containing the
c    information on the wavefunctions. In example above, h2o.WFSX
c****************************************************************************



        implicit none

        integer, parameter :: dp = selected_real_kind(10,100)
        integer, parameter :: sp = selected_real_kind(5,10)

        integer io,iu, nk, nspin_blocks, ik, iik, ispin, iispin,
     .          nwflist, iw, indwf, j, nuotot, jj, nspin_flag

        character fname*33, oname*33
        integer, allocatable, dimension(:) :: iaorb,iphorb,cnfigfio
        character(len=20), allocatable, dimension(:) :: symfio,labelfis

        real(sp), allocatable, dimension(:,:) :: psi
        logical gamma, non_coll

        real(dp) k(3), energy, thress, norm
 
!        read(5,*) fname
!        read(5,*) oname
!        read(5,*) thress

        fname = "WFSX"
        oname = "FMT_WFSX"
        thress = 0.0001
        thress = thress**2  ! We use the square later on

        iu = 10
        io = 11

        open(iu, file=fname, form='unformatted', status='old' )
        open(io, file=oname, form='formatted', status='unknown' )

        rewind (iu)

        read(iu) nk, gamma
        read(iu) nspin_flag
        read(iu) nuotot

        if (nspin_flag == 4) then
           non_coll = .true.
           nspin_blocks = 1
        else
           non_coll = .false.
           nspin_blocks = nspin_flag
        endif
        
        if (gamma) then
           if (non_coll) then
              allocate(psi(4,nuotot))
           else
              allocate(psi(1,nuotot))
           endif
        else
           if (non_coll) then
              allocate(psi(4,nuotot))
           else
              allocate(psi(2,nuotot))
           endif
        endif

        allocate(iaorb(nuotot), labelfis(nuotot),
     $           iphorb(nuotot), cnfigfio(nuotot),
     $           symfio(nuotot))

         read(iu) (iaorb(j),labelfis(j),
     .            iphorb(j), cnfigfio(j),
     .            symfio(j), j=1,nuotot)

        write(io,*)
        write(io,'(a22,2x,i6)') 'Nr of k-points = ',nk
        if (non_coll) then
           write(io,'(a22,2x,i6)') 'Nr of Spins blocks ' //
     $                             '(non-collinear)  = ', nspin_blocks
        else
           write(io,'(a22,2x,i6)') 'Nr of Spins blocks = ',nspin_blocks
        endif
        write(io,'(a22,2x,i6)') 'Nr of basis orbs = ',nuotot
        write(io,*)


        do iik = 1,nk
          do iispin = 1,nspin_blocks

          read(iu) ik,k(1),k(2),k(3)
          if (ik .ne. iik) stop 'error in index of k-point'
          read(iu) ispin
          if (.not. non_coll) then
             if (ispin .ne. iispin) stop 'error in index of spin'
          endif
          read(iu) nwflist

          write(io,*)
          write(io,'(a72)')    ' ***************************************
     .********************************'
          write(io,'(a22,2x,i6,2x,3f10.6)') 'k-point = ',ik,
     .         k(1),k(2),k(3)
          if (.not. non_coll) then
             write(io,'(a22,2x,i6)') 'Spin component = ',ispin
          endif
          write(io,'(a22,2x,i6)') 'Num. wavefunctions = ',nwflist



C Loop over wavefunctions 

          do iw = 1,nwflist

            read(iu) indwf
            read(iu) energy

            write(io,*)
            write(io,'(a22,2x,i6)') 'Wavefunction = ', indwf
            write(io,'(a22,2x,f10.6)') 'Energy (eV) = ', energy
            write(io,'(a72)')  ' ---------------------------------------
     .--------------------------------'
            if (non_coll) then
               write(io,'(a72)')
     $                      ' Atom  Species Orb-global Orb-in-atom '//
     $                      '  .Orb-type    [Re(psi)  Im(psi)]Up ' //
     $                      ' [Re(psi)  Im(psi)]Down '
            else
               write(io,'(a72)')
     $              ' Atom  Species Orb-global  Orb-in-atom ' //
     $              ' .Orb-type      Re(psi)   Im(psi)'
            endif
            
            read(iu) (psi(1:,j), j=1,nuotot)
            do jj = 1,nuotot
               if (non_coll) then
                  norm = dot_product(psi(1:,jj),psi(1:,jj))  ! same
               else
                  norm = dot_product(psi(1:,jj),psi(1:,jj))
               endif
               
              if (norm .lt. thress) CYCLE
              write(io,
     $            '(i6,5x,a10,1x,i10,2x,i3,2x,i1,a20,2(2(f10.6),2x))') 
     .          iaorb(jj),labelfis(jj),jj,
     .          iphorb(jj), cnfigfio(jj),
     .          symfio(jj), psi(1:,jj)

            enddo

            write(io,'(a72)')  ' ---------------------------------------
     .--------------------------------'

          enddo
        enddo
      enddo


      close (iu)
      end
