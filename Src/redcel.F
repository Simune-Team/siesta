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
      subroutine redcel( alat, ucell, scell, mscell )

c *********************************************************************
c Reads the unit cell (and supercell, if present) vectors from the
c data file.
c Written by J.M.Soler
c ********* Output ****************************************************
c real*8 alat       : Lattice constant
c real*8 ucell(3,3) : Unit cell vectors ucell(ixyz,ivec)
c real*8 scell(3,3) : Supercell vectors (equal to ucell if no supercell)
c integer mscell(3,3): Supercell vectors in units of ucell, i.e.
c                      scell(ix,i) = Sum_j( ucell(ix,j)*mscell(j,i) )
c *********************************************************************
C
C  Modules
C
      use precision, only : dp
      use parallel,  only : Node
      use sys,       only : die
      use fdf
      implicit          none

      integer, intent(out)  ::      mscell(3,3)
      real(dp), intent(out) ::      alat, scell(3,3), ucell(3,3)

C Internal variables
      integer           i, ix, j
      real(dp)          alp, alplp, betlp, blp, cell(3,3), clp,
     .                  gamlp, pi, xxx

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

C Lattice constant
      alat = fdf_physical('LatticeConstant',0.d0,'Bohr')

C Lattice vectors

      if ( fdf_block('LatticeParameters',bfdf) .and.
     .     fdf_block('LatticeVectors',bfdf) ) then
         call die('redcel: ERROR: Lattice vectors doubly ' //
     $       'specified: by LatticeVectors and by LatticeParameters.')
      endif

      if ( fdf_block('LatticeParameters',bfdf) ) then
        if (.not. fdf_bline(bfdf, pline))
     .    call die('redcel: ERROR in LatticeParameters block')
        alp = fdf_breals(pline,1)
        blp = fdf_breals(pline,2)
        clp = fdf_breals(pline,3)
        alplp = fdf_breals(pline,4)
        betlp = fdf_breals(pline,5)
        gamlp = fdf_breals(pline,6)
*       write(6,'(a)')
*    .   'redcel: Lattice Parameters (units of Lattice Constant) ='
*       write(6,'(a,3f10.5,3f9.3)')
*    .   'redcel: ',alp,blp,clp,alplp,betlp,gamlp
        pi = acos(-1.d0)
        alplp = alplp * pi/180.d0
        betlp = betlp * pi/180.d0
        gamlp = gamlp * pi/180.d0
        cell(1,1) = alp
        cell(2,1) = 0.d0
        cell(3,1) = 0.d0
        cell(1,2) = blp * cos(gamlp)
        cell(2,2) = blp * sin(gamlp)
        cell(3,2) = 0.d0
        cell(1,3) = clp * cos(betlp)
        xxx = (cos(alplp) - cos(betlp)*cos(gamlp))/sin(gamlp)
        cell(2,3) = clp * xxx
        cell(3,3) = clp * sqrt(sin(betlp)*sin(betlp) - xxx*xxx)
      elseif ( fdf_block('LatticeVectors',bfdf) ) then
        do i = 1,3
          if (.not. fdf_bline(bfdf, pline))
     .      call die('redcel: ERROR in LatticeVectors block')
          cell(1,i) = fdf_breals(pline,1)
          cell(2,i) = fdf_breals(pline,2)
          cell(3,i) = fdf_breals(pline,3)
        enddo
      else
        do i = 1,3
          do j  = 1,3
            cell(i,j) = 0.d0
          enddo
          cell(i,i) = 1.d0
        enddo
      endif

      cell = alat * cell

C     Super cell 

      if ( fdf_block('SuperCell',bfdf) ) then
        if ( alat .eq. 0.d0 ) then
          call die('redcel: ERROR: LatticeConstant required' //
     .             ' to define SuperCell')
        endif
        do i = 1,3
          if (.not. fdf_bline(bfdf, pline))
     .      call die('redcel: ERROR in SuperCell block')
          mscell(1,i) = fdf_breals(pline,1)
          mscell(2,i) = fdf_breals(pline,2)
          mscell(3,i) = fdf_breals(pline,3)
        enddo
      else
        do i = 1,3
          do j = 1,3
            mscell(j,i) = 0
          enddo
          mscell(i,i) = 1
        enddo
      endif

      ucell = cell

      do i = 1,3
        do ix = 1,3
          scell(ix,i) = ucell(ix,1) * mscell(1,i) +
     .                  ucell(ix,2) * mscell(2,i) +
     .                  ucell(ix,3) * mscell(3,i)
        enddo
      enddo

      end
