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
      subroutine coor(na,cell)

C Reads atomic coordinates and format in which they are given, to be
C transformed into Bohr cartesian for internal handling. 
C It also shifts them all according to AtomicCoordinatesOrigin.
C Written by E. Artacho, December 1997, on the original piece of the
C redata subroutine written by P. Ordejon in December 1996.
C Modified by J.M.Soler. August 1998.

!! Modified by Alberto Garcia, May 16, 2000

C integer na  : Number of atoms in unit cell
!               Note: the user might specify a supercell in the
!                     input file, and the coordinates are generated
!                     accordingly in this routine, but this does
!                     not have anything to do with the 'virtual'
!                     supercell used to deal with k-point sampling.

      use precision, only : dp
      use sys
      use fdf
      use alloc
      use parallel,  only : Node
      use zmatrix
      use m_mpi_utils, only : broadcast
      use siesta_geom
      implicit none

      integer,  intent(out) :: na
      real(dp), intent(out) :: cell(3,3)

C Internal variables and arrays
      character             :: acf*22, acf_default*22
      logical               :: isdiag, leqi
      integer               :: i, ia, ic, iscale, iua, iunit, ix,
     .                         mscell(3,3), ncells, nsc(3), nua
      real(dp)              :: alat, dcell(3,3), dscell(3,3), 
     .                         origin(3),
     .                         volcel, xac(3)
      external          digcel, redcel, volcel

      logical               :: lOrigin
      data origin /3*0.0_dp/


      if (Node.eq.0) then
C Read unit cell and supercell
        call redcel( alat, ucell, scell, mscell )

C Set cell argument to "supercell" for use by the program
        cell = scell

C Format of atomic coordinates
        acf_default = 'NotScaledCartesianBohr'
        acf = fdf_string('AtomicCoordinatesFormat',acf_default)
        if (leqi(acf,'NotScaledCartesianBohr') .or.
     .      leqi(acf,'Bohr') ) then
          iscale = 0
          write(6,'(a,a)')
     .     'coor:   Atomic-coordinates input format  = ',
     .     '    Cartesian coordinates'
          write(6,'(a,a)')
     .     'coor:                                      ',
     .     '      (in Bohr units)'
        else if (leqi(acf,'NotScaledCartesianAng') .or.
     .           leqi(acf,'Ang') ) then
          iscale = 1
          write(6,'(a,a)')
     .     'coor:   Atomic-coordinates input format  = ',
     .     '    Cartesian coordinates'
          write(6,'(a,a)')
     .     'coor:                                    ',
     .     '      (in Angstroms)'
        else if (leqi(acf,'ScaledCartesian')) then
          if (alat.eq.0.0_dp) then
             write(6,"(/,2a)") 'coor: ERROR: Explicit lattice ',
     .         'constant is needed for ScaledCartesian format'
             call die
          endif
          iscale = 2
          write(6,'(a,a)')
     .     'coor:   Atomic-coordinates input format  = ',
     .     '    Cartesian coordinates'
          write(6,'(a,a)')
     .     'coor:                                    ',
     .     '      (in units of alat)'
        else if (leqi(acf,'ScaledByLatticeVectors') .or. 
     .           leqi(acf,'Fractional') ) then
          if (alat.eq.0.0_dp) then
             write(6,"(/,2a)") 'coor: ERROR: Explicit lattice ',
     .         'constant is needed for Fractional format'
             call die
          endif
          iscale = 3
          write(6,'(a,a)')
     .     'coor:   Atomic-coordinates input format  = ',
     .     '    Fractional'
        else
          write(6,"(/,'coor: ',72(1h*))")
          write(6,"('coor:                  INPUT ERROR')")
          write(6,'(a)') 'coor: '
          write(6,'(2a)') 'coor: You must use one of the following',
     .                              ' coordinate scaling options:'
          write(6,'(a)') 'coor:     - NotScaledCartesianBohr (or Bohr)'
          write(6,'(a)') 'coor:     - NotScaledCartesianAng (or Ang) '
          write(6,'(a)') 'coor:     - ScaledCartesian                '
          write(6,'(2a)') 'coor:     - ScaledByLatticeVectors ',
     .                                                 '(or Fractional)'
          write(6,"('coor: ',72(1h*))")

          call die()

        endif

C Read atomic coordinates and species
        na = fdf_integer('NumberOfAtoms',0)
        if (na.eq.0) call die("Must specify number of atoms!")
!
!     Check if we need more space to accommodate supercell
!     (still a "real" supercell, specified for convenience!)
!
        if (volcel(ucell) .lt. 1.d-8) then
          ncells = 1
        else
          ncells = nint( volcel(scell) / volcel(ucell) )
        endif

        nua = na
        na  = na * ncells

C Find origin with which to translate all coordinates
        lOrigin = fdf_block('AtomicCoordinatesOrigin',iunit)
        if (lOrigin) then 
          read(iunit,*) (origin(i),i=1,3)
        endif

      endif
C  end of Node.eq.0

C Distribute data

      call broadcast(iscale)
      call broadcast(mscell(1:3,1:3))
      call broadcast(na)
      call broadcast(nua)
      call broadcast(ncells)
      call broadcast(lOrigin)
      call broadcast(ucell(1:3,1:3))
      call broadcast(cell(1:3,1:3))
      call broadcast(origin(1:3))
      call broadcast(alat)

C Set array dimensions
      nullify(isa,xa)
      call re_alloc(isa,1,na,name='isa',routine='coor')
      call re_alloc(xa,1,3,1,na,name='xa',routine='coor')

C Attempt to read a Z matrix
      call read_Zmatrix(nua,isa,alat,ucell,lOrigin,origin)

      if (lUseZmatrix) then
C Generate Cartesian coordinates from Z-matrix
        call Zmat_to_Cartesian(xa)

      else
        if (Node.eq.0) then
C If Z matrix hasn't been found, read regular atomic coordinates
          if ( fdf_block('AtomicCoordinatesAndAtomicSpecies',iunit) )
     .        then
            do ia = 1,nua
              read(iunit,*) (xa(i,ia), i=1,3), isa(ia)
            enddo
          else
            call die("coo: You must specify the atomic coordinates")
          endif

C Scale atomic coordinates
C   Coord. option = 0 => Do nothing
C   Coord. option = 1 => Multiply by 1./0.529177 (Ang --> Bohr)
C   Coord. option = 2 => Multiply by lattice constant
C   Coord. option = 3 => Multiply by lattice vectors

          if (iscale .eq. 1) then
            xa = xa / 0.529177_dp 
          elseif (iscale .eq. 2) then
            xa = xa * alat
          elseif (iscale .eq. 3) then
            do ia = 1,nua
              do ix = 1,3
                xac(ix) = xa(ix,ia)
              enddo
              do ix = 1,3
                xa(ix,ia) = ucell(ix,1) * xac(1) +
     .                      ucell(ix,2) * xac(2) +
     .                      ucell(ix,3) * xac(3)
              enddo
            enddo
          endif

C Translate atomic coordinates by origin
          if (lOrigin) then
            do ia = 1,nua
              do i = 1,3
                xa(i,ia) = xa(i,ia) + origin(i)
              enddo
            enddo
          endif

        endif
C end of Node.eq.0
        
C distribute data
        call broadcast(isa(1:na))
        call broadcast(xa(1:3,1:na))


C Expand the coordinates to whole supercell

        if (ncells.gt.1) then

C Find equivalent diagonal combination of cell/supercell
          call digcel( ucell, mscell, dcell, dscell, nsc, isdiag )

C Expand coordinates
          call superx_coor( dcell, nsc, nua, na, xa, dscell )

C Expand index isa
          ia = 0
          do ic = 1,ncells
            do iua = 1,nua
              ia = ia + 1
              isa(ia) = isa(iua)
            enddo
          enddo

        endif


      endif
C  end of not lUseZmatrix

        ! Construct references: cannot use isa alone, since 
        ! refs may not start with a number.
        allocate(cisa(na))
        do ia=1, nua
           write(cisa(ia), '("siesta:e",i3.3)') isa(ia)
        enddo

      contains

      subroutine superx_coor(ucell, NSC, na, maxa, XA, scell )
C **********************************************************************
C Generates supercell vectors and atomic positions.
C Written by J.M.Soler, August 1998
C *************** Input ************************************************
C Real*8  ucell(3,3)  : Unit cell vectors ucell(Ixyz,Ivector)
C Integer NSC(3)      : Number of cells in each supercell direction:
C                         scell(ix,i) = ucell(ix,i) * NSC(i)
C Integer NA          : Number of atoms in unit cell
C Integer MAXA        : Second dimension of XA
C *************** Input and output *************************************
C Real*8  XA(3,MAXA)  : Atomic positions in unit cell (input) and
C                       in supercell (output), in cartesian coord.
C Real*8  scell(3,3)  : Supercell vectors
C *********** Units ****************************************************
C Units of CELL and XA are arbitrary but must be the same
C *********** Behavior *************************************************
C - If NA*NCELLS > MAXA (where NCELLS is the total number of cells),
C   the supercell atomic coordinates are not generated.
C - The first supercell atoms are those of the initial unit cell, i.e.
C   the positions XA(i,ia) for (ia.le.NA) are not modified.
C - The remaining atoms are ordered by unit cells, i.e. the atom ia
C   is equivalent to the unit-cell atom ja=MOD(ia-1,NA)+1
C **********************************************************************

      use precision, only : dp

      implicit          none

      integer           maxa, na, NSC(3)
      real(dp)          scell(3,3), ucell(3,3), XA(3,maxa)

C Internal variables
      integer           I, I1, I2, I3, IA, IX, JA, ncells
      real(dp)          XC(3)

C Find supercell vectors
      do I = 1,3
        do IX = 1,3
          scell(IX,I) = ucell(IX,I) * NSC(I)
        enddo
      enddo

C Expand atomic positions to supercell
      ncells = NSC(1) * NSC(2) * NSC(3)
      if (na*ncells .le. maxa) then
        IA = 0
        do I3 = 0,NSC(3)-1
        do I2 = 0,NSC(2)-1
        do I1 = 0,NSC(1)-1
          do IX = 1,3
            XC(IX) = ucell(IX,1)*I1 + ucell(IX,2)*I2 + ucell(IX,3)*I3
          enddo
          do JA = 1,na
            IA = IA + 1
            do IX = 1,3
              XA(IX,IA) = XA(IX,JA) + XC(IX)
            enddo
          enddo
        enddo
        enddo
        enddo
      endif

      end subroutine superx_coor

      end subroutine coor
