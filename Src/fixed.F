! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      MODULE m_fixed
      use precision
      implicit none

      PUBLIC                     :: fixed, resetFixedPointers
      PRIVATE
      ! maxc  : maximum number of constraints
      integer                    :: maxc = 1000
      integer                    :: nc
      real(dp)                   :: xs(6)
      logical                    :: frstme =.true.
      character(len=10), pointer :: ctype(:)
      integer,           pointer :: iac(:)
      real(dp),          pointer :: xc(:,:)

      CONTAINS

      subroutine resetFixedPointers( )
      use alloc, only : de_alloc
      implicit none
      CHARACTER(LEN=*), PARAMETER :: MYNAME =  'fixed'

      if (associated(ctype))
     &  deallocate( ctype )

      if (associated(iac))
     &  call de_alloc( iac, name='iac', routine=MYNAME )

      if (associated(xc))
     &  call de_alloc( xc, name='xc', routine=MYNAME )

      end subroutine resetFixedPointers

c **********************************************************************
c Reads and imposes required constraints to atomic displacements by
c making zero the forces in those directions. Constraints are specified
c by the FDF data block GeometryConstraints (see example below).
c Only types position and routine implemented in this version.
c Written by J.M.Soler. Feb., 1998
c Modified by P. Ordejon to output the total number of constraints
C    imposed.  June, 2003
c *********** INPUT ****************************************************
c real*8  cell(3,3)    : Lattice vectors
c real*8  stress( 3,3) : Stress tensor
c integer na           : Number of atoms
c integer isa(na)      : Species indexes
c real*8  amass(na)    : Atomic masses
c real*8  xa(3,na)     : Atomic cartesian coordinates
c real*8  fa(3,na)     : Atomic forces
c *********** OUTPUT ***************************************************
c real*8  cstress( 3,3) : Constrained stress tensor
c real*8  cfa(3,na)     : Constrained atomic forces
c integer ntcon         : Total number of position constraints imposed
c *********** UNITS ****************************************************
c Units are arbitrary but cell and xa must be in the same units
c *********** BEHAVIOUR ************************************************
c cstress may be the same physical array as stress, and cfa the same 
c as fa, i.e. it is allowed:
c     call fixed( cell, stress, na, isa, amass, xa, fa, stress, fa, ntcon )
C NOTE: This is NOT good practice, and might be flagged by some compilers
c *********** USAGE ****************************************************
c Example: consider a diatomic molecule (atoms 1 and 2) above a surface, 
c represented by a slab of 5 atomic layers, with 10 atoms per layer.
c To fix the cell height, the slab's botom layer (last 10 atoms),
c the molecule's interatomic distance, its height above the surface
c (but not its azimutal orientation and lateral position), and the
c relative height of the two atoms:
c
c   %block GeometryConstraints
c   cellside   c 
c   cellangle  alpha  beta  gamma
c   position  from -1 to -10
c   rigid  1  2
c   center 1  2   0.0  0.0  1.0
c   routine constr
c   stress 1  2  3
c   %endblock GeometryConstraints
c
c where constr is the following user-written subroutine:
c
c      subroutine constr( cell, na, isa, amass, xa, stress, fa, ntcon )
cc real*8  cell(3,3)    : input lattice vectors (Bohr)
cc integer na           : input number of atoms
cc integer isa(na)      : input species indexes
cc real*8  amass(na)    : input atomic masses
cc real*8  xa(3,na)     : input atomic cartesian coordinates (Bohr)
cc real*8  stress( 3,3) : input/output stress tensor (Ry/Bohr**3)
cc real*8  fa(3,na)     : input/output atomic forces (Ry/Bohr)
cc integer ntcon        : total number of position constraints, accumulative
c      integer na, isa(na), ntcon
c      double precision amass(na), cell(3,3), fa(3,na),
c     .                 stress(3,3), xa(3,na), fz
c      fz =  fa(3,1) + fa(3,2) 
c      fa(3,1) = fz*amass(1)/(amass(1)+amass(2))
c      fa(3,2) = fz*amass(2)/(amass(1)+amass(2))
c      ntcon = ntcon+1
c      end
c **********************************************************************

      subroutine fixed( cell, stress, na, isa, amass, xa, fa,
     .                  cstress, cfa, ntcon )
      use precision
      use parallel,   only : Node
      use fdf
      use alloc

      implicit none

      integer,  intent(in)  :: na, isa(na)
      real(dp), intent(in)  :: amass(na), cell(3,3), 
     .                         fa(3,na), stress(3,3), xa(3,na)
      integer, intent(out)  :: ntcon
      real(dp), intent(out) :: cfa(3,na), cstress(3,3)

      ! Internal parameters
      ! maxl  : maximum number of input constraint lines
      ! maxw  : maximum number of items in an input constraint line
      integer, parameter :: maxl= 10000,  maxw=  1000


      ! Internal variables
      logical found
      character name1*10, name2*10, name3*10, name4*10


      integer i, ia, ia1, ia2, ia3, ic, il, ix,
     &        jx, ni, nn, nr, nv, maxcnew



      integer   integs(maxw)
      real(dp)  ddot, fxc, xnorm
      character(len=10), pointer :: ctype_tmp(:)

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      external ddot, constr
#ifdef DEBUG
      call write_debug( '    PRE fixed' )
#endif

      ! Read constraint data only the first time
      if (frstme) then
!       Allocate memory for constraints - Preserve contents
        nullify( ctype )
        allocate( ctype(maxc) )

        nullify( iac )
        call re_alloc( iac, 1, maxc, 'iac', 'fixed' )

        nullify( xc )
        call re_alloc( xc, 1, 3, 1, maxc, 'xc', 'fixed' )

        ! Initialise stress constraints to unconstrained state
        xs(1:6) = 0.0d0

        nc = 0

        ! Look for constraints data block
        found = fdf_block('GeometryConstraints',bfdf)
        if (.not. found) goto 30

        ! Loop on data lines
        il = 1
        do while(fdf_bline(bfdf,pline) .and. (il .le. maxl))

          ! Check if constraints are finished
          nn = fdf_bnnames(pline)
          nv = fdf_bnvalues(pline)
          ni = fdf_bnintegers(pline)
          nr = fdf_bnreals(pline)
          name1 = fdf_bnames(pline,1)

          ! Select type of constraint
          if (name1 .eq. 'routine') then
            if ((nn.eq.1) .or. (fdf_bnames(pline,2).eq.'constr')) then
              nc = nc + 1
              if (nc.eq.maxc) then
                maxcnew = nc + nint(0.1*nc)

                allocate(ctype_tmp(maxcnew))
                ctype_tmp(1:maxc) = ctype
                deallocate(ctype)
                ctype => ctype_tmp
!                call re_alloc( ctype, 1, maxcnew, name='ctype',            
!     &                         routine='fixed', copy=.true. )
                call re_alloc( iac, 1, maxcnew, name='iac',            
     &                          routine='fixed', copy=.true. )
                call re_alloc( xc, 1, 3, 1, maxcnew, name='xc',            
     &                         routine='fixed', copy=.true. )
                maxc = maxcnew
              endif
              ctype(nc) = 'routine'
            else
              if (Node.eq.0) then
                write(6,*) 'fixed: ERROR: user-constraints routine',
     &                   ' must be called constr'
              endif
            endif

          elseif (name1 .eq. 'position') then

            ! Check syntax
            if ((nr .ne. 0) .and. (nr .ne. 3)) then
              if (Node.eq.0) then
                write(6,'(a,/,a)')
     &            'fixed: syntax ERROR in %block GeometryConstraints:',
     &            trim(fdf_getline(bfdf%mark))
              endif
              cycle
            endif

            ! Find constrained atoms
            if (nn .gt. 1) then
              ! Atoms specified by range. Make list of them.
              name2 = fdf_bnames(pline,2)
              name3 = fdf_bnames(pline,3)
              if (nn .eq. 4) name4 = fdf_bnames(pline,4)
              if ((name2 .eq. 'from') .and. (name3 .eq. 'to')) then
                ia1 = fdf_bintegers(pline,1)
                ia2 = fdf_bintegers(pline,2)
                if (ia1 .lt. 0) ia1 = na + ia1 + 1
                if (ia2 .lt. 0) ia2 = na + ia2 + 1
                if ( min(ia1,ia2) < 1 .or. na < max(ia1,ia2) ) then
                 call die('Atom constraints must be in range: [-na;na]')
                end if
                if ((nn .eq. 4) .and. (name4 .eq. 'step')) then
                  ia3 = abs(fdf_bintegers(pline,3))
                else
                  ia3 = 1
                endif
                ni = 0
                do ia= min(ia1,ia2), max(ia1,ia2), ia3
                  ni = ni + 1
                  integs(ni) = ia
                enddo
              else
                if (Node.eq.0) then
                  write(6,'(a,/,a)')
     &             'fixed: syntax ERROR in %block GeometryConstraints:',
     &             trim(fdf_getline(bfdf%mark))
                endif
                cycle
              endif
            else
              do ia = 1, ni
                integs(ia) = fdf_bintegers(pline,ia)
              enddo
            endif

            ! Store position constraints
            do i= 1, ni
              ia = integs(i)
              if (ia .lt. 0) ia = na + ia + 1
              if (nr .eq. 0) then
                ! Make one constraint for each cartesian coordinate
                do ix= 1, 3
                  nc = nc + 1
                  if (nc .eq. maxc) then
                    maxcnew = nc + max(nint(0.1*nc),1)
                    call re_alloc( ctype, 1, maxcnew, name='ctype',            
     &                             routine='fixed', copy=.true. )
                    call re_alloc( iac, 1, maxcnew, name='iac',            
     &                             routine='fixed', copy=.true. )
                    call re_alloc( xc, 1, 3, 1, maxcnew, name='xc',            
     &                             routine='fixed', copy=.true. )
                    maxc = maxcnew
                  endif
                  ctype(nc) = 'position'
                  iac(nc) = ia
                  do jx= 1, 3
                    xc(jx,nc) = 0.0d0
                  enddo
                  xc(ix,nc) = 1.0d0
                enddo
              elseif (nr .eq. 3) then
                ! Make only one constraint in the specified direction
                nc = nc + 1
                if (nc .eq. maxc) then
                  maxcnew = nc + max(nint(0.1*nc),1)
                  call re_alloc( ctype, 1, maxcnew, name='ctype',            
     &                           routine='fixed', copy=.true. )
                  call re_alloc( iac, 1, maxcnew, name='iac',            
     &                           routine='fixed', copy=.true. )
                  call re_alloc( xc, 1, 3, 1, maxcnew, name='xc',            
     &                         routine='fixed', copy=.true. )
                  maxc = maxcnew
                endif
                ctype(nc) = 'position'
                iac(nc) = ia
                xnorm = sqrt(fdf_breals(pline,1)**2 +
     &                       fdf_breals(pline,2)**2 +
     &                       fdf_breals(pline,3)**2)
                do ix = 1,3
                  xc(ix,nc) = fdf_breals(pline,ix) / xnorm
                enddo
              endif
            enddo

          elseif (name1 .eq. 'stress') then

            do ic= 1, nv
              nc = nc + 1
              ctype(nc) = 'stress'
              ix = nint(fdf_bvalues(pline,ic))
              if ((ix .ge. 1) .and. (ix .le. 6)) xs(ix) = 1.0d0
            enddo

          else
            if (Node.eq.0) then
              write(6,*) 'fixed: ERROR: sorry, constraint type ',
     &                  name1, ' not implemented yet'
            endif
          endif
          il = il + 1
        enddo
   30   continue
        frstme = .false.
      endif

      ! Copy stress and forces to output arrays 
      do ix = 1,3
        do jx = 1,3
          cstress(jx,ix) = stress(jx,ix)
        enddo
      enddo
      do ia = 1,na
        do ix = 1,3
          cfa(ix,ia) = fa(ix,ia)
        enddo
      enddo

      ! Apply constraints
      ntcon = 0
      do ic = 1,nc
        if (ctype(ic) .eq. 'routine') then
          call constr( cell, na, isa, amass, xa, cstress, cfa, ntcon )
        elseif (ctype(ic) .eq. 'position') then
          ntcon = ntcon + 1
          ia = iac(ic)
          fxc = ddot(3,cfa(1,ia),1,xc(1,ic),1)
          do ix = 1,3
            cfa(ix,ia) = cfa(ix,ia) - fxc * xc(ix,ic)
          enddo
        endif
      enddo

      ! Apply stress constraints
      cstress(1,1) = cstress(1,1) - xs(1)*cstress(1,1)
      cstress(2,2) = cstress(2,2) - xs(2)*cstress(2,2)
      cstress(3,3) = cstress(3,3) - xs(3)*cstress(3,3)
      cstress(3,2) = cstress(3,2) - xs(4)*cstress(3,2)
      cstress(2,3) = cstress(2,3) - xs(4)*cstress(2,3)
      cstress(3,1) = cstress(3,1) - xs(5)*cstress(3,1)
      cstress(1,3) = cstress(1,3) - xs(5)*cstress(1,3)
      cstress(1,2) = cstress(1,2) - xs(6)*cstress(1,2)
      cstress(2,1) = cstress(2,1) - xs(6)*cstress(2,1)

!
!     We could treat here extra "constraints" affecting
!     the stress which are implicitly applied by all
!     "relaxation" routines:
!
!     1. A possible "target stress" (this can be thought of
!        as an external generalized pressure).
!     2. Symmetrization
!
!     3. (maybe) The "constant-volume" constraint.
!
!     If carried out, this would mean that 'cstress' on exit
!     from this routine is the "stress_dif" of the relaxation routines.
!
!     It should be up to the user to make sure that the options are consistent.
!
!     - what happens with a user-defined routine? (precedence?)
!     - should we do this only for the "relaxation" MD routines?
!
!       In pr and npr, there is currently only a "target pressure", but its effect
!       on the 'aux' variable could be generalized to a "target stress", and the
!       'potential energy of the PR variables' generalized from pt*vol to sigma_target*epsilon (?)
!       The problem here is that the equations of motion for the cell variables might
!       get confused.
!
!       In anneal there is an explicit "target stress"
!     
!     To avoid changing the current defaults (which, for example, do not allow a "target_stress"
!     in pr and npr), it might be better to apply the (current) "target stress" option conditions in 
!     each block of "siesta_move". For idyn=0, this is straightforward, and it has the added benefit
!     that a lot of repeated code can be removed from the routines. 
!
!     It can be also added as a new feature to 'pr' if iquench>0. (and set 'tp=0')
!



#ifdef DEBUG
      call write_debug( '    POS fixed' )
#endif
      end subroutine fixed
      
      END MODULE m_fixed
