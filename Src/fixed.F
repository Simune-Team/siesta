      subroutine fixed( cell, stress, na, isa, amass, xa, fa,
     .                  cstress, cfa )

c **********************************************************************
c Reads and imposes required constraints to atomic displacements by
c making zero the forces in those directions. Constraints are specified
c by the FDF data block GeometryConstraints (see example below).
c Only types position and routine implemented in this version.
c Written by J.M.Soler. Feb., 1998
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
c *********** UNITS ****************************************************
c Units are arbitrary but cell and xa must be in the same units
c *********** BEHAVIOUR ************************************************
c cstress may be the same physical array as stress, and cfa the same 
c as fa, i.e. it is allowed:
c     call fixed( cell, stress, na, isa, amass, xa, fa, stress, fa )
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
c   %endblock GeometryConstraints
c
c where constr is the following user-written subroutine:
c
c      subroutine constr( cell, na, isa, amass, xa, stress, fa )
cc real*8  cell(3,3)    : input lattice vectors (Bohr)
cc integer na           : input number of atoms
cc integer isa(na)      : input species indexes
cc real*8  amass(na)    : input atomic masses
cc real*8  xa(3,na)     : input atomic cartesian coordinates (Bohr)
cc real*8  stress( 3,3) : input/output stress tensor (Ry/Bohr**3)
cc real*8  fa(3,na)     : input/output atomic forces (Ry/Bohr)
c      integer na, isa(na)
c      double precision amass(na), cell(3,3), fa(3,na),
c     .                 stress(3,3), xa(3,na), fz
c      fz = ( fa(3,1) + fa(3,2) ) / 2.d0
c      fa(3,1) = fz
c      fa(3,2) = fz
c      end
c **********************************************************************

C
C  Modules
C
      use precision
      use fdf
#ifdef MPI
      use mpi
#endif

      implicit          none

      integer           na, isa(na)
      double precision  amass(na), cell(3,3), cfa(3,na), cstress(3,3),
     .                  fa(3,na), stress(3,3), xa(3,na)

c Internal parameters
c maxc  : maximum number of constraints
c maxl  : maximum number of input constraint lines
c maxw  : maximum number of items in an input constraint line
      integer maxl, maxw
      parameter ( maxl  = 10000 )
      parameter ( maxw  =  1000 )
      integer, save :: maxc

c Internal variables
      logical
     .  found, frstme
      character
     .  line*130,
     .  name1*10, name2*10, name3*10, name4*10, names*130
      character(len=10), dimension(:), allocatable, save ::
     .  ctype, cbuffer
      integer
     .  i, ia, ia1, ia2, ia3, ic, il, integs(maxw), iu, ix,
     .  jx, lastch, lch(0:maxw), nc, ni, nn, nr, nv, Node,
     .  maxcnew
      integer, dimension(:), allocatable, save ::
     .  iac, ibuffer
#ifdef MPI
      integer MPIerror
#endif
      double precision
     .  dot, fxc, reals(maxw), values(maxw), xnorm
      double precision, dimension(:,:), allocatable, save ::
     .  xc, dpbuffer2
      external
     .  dot, constr, parse, memory
      save
     .  frstme, nc

      data
     .  frstme /.true./,
     .  maxc   / 1000 /

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

C Read constraint data only the first time
      if (frstme) then

C Allocate memory for constraints - do only once and preserve contents
        allocate(ctype(maxc))
        allocate(iac(maxc))
        call memory('A','I',maxc,'fixed')
        allocate(xc(3,maxc))
        call memory('A','D',3*maxc,'fixed')

        nc = 0

C       Look for constraints data block
        if (Node.eq.0) then
          found = fdf_block('GeometryConstraints',iu)
        endif
#ifdef MPI
        call MPI_Bcast(found,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
        if (.not.found) goto 30

C       Loop on data lines
        do 20 il = 1,maxl

C         Read and parse data line
          read(iu,'(a)',end=30) line
          lastch = index(line,'#') - 1
          if (lastch .le. 0) lastch = len(line)
          call parse( line(1:lastch), nn, lch, names, nv, values,
     .                ni, integs, nr, reals )

c         Check if constraints are finished
          name1 = names(lch(0)+1:lch(1))
          if (name1 .eq. '%end' .or.
     .        name1 .eq. '%endblock') then
            goto 30

c         Select type of constraint
          elseif (name1 .eq. 'routine') then
            if (nn.eq.1 .or. names(lch(1)+1:lch(2)) .eq. 'constr') then
              nc = nc + 1
              if (nc.eq.maxc) then
                maxcnew = nc + nint(0.1*nc)
C
                allocate(cbuffer(maxc))
                cbuffer(1:maxc) = ctype(1:maxc)
                deallocate(ctype)
                allocate(ctype(maxcnew))
                ctype(1:maxc) = cbuffer(1:maxc)
                deallocate(cbuffer)
C
                allocate(ibuffer(maxc))
                call memory('A','I',maxc,'fixed')
                ibuffer(1:maxc) = iac(1:maxc)
                call memory('D','I',size(iac),'fixed')
                deallocate(iac)
                allocate(iac(maxcnew))
                call memory('A','I',maxcnew,'fixed')
                iac(1:maxc) = ibuffer(1:maxc)
                call memory('D','I',size(ibuffer),'fixed')
                deallocate(ibuffer)
C
                allocate(dpbuffer2(3,maxc))
                call memory('A','D',maxc,'fixed')
                dpbuffer2(1:3,1:maxc) = xc(1:3,1:maxc)
                call memory('D','D',size(xc),'fixed')
                deallocate(xc)
                allocate(xc(3,maxcnew))
                call memory('A','D',3*maxcnew,'fixed')
                xc(1:3,1:maxc) = dpbuffer2(1:3,1:maxc)
                call memory('D','D',size(dpbuffer2),'fixed')
                deallocate(dpbuffer2)
C
                maxc = maxcnew
              endif
              ctype(nc) = 'routine'
            else
              if (Node.eq.0) then
                write(6,*) 'fixed: ERROR: user-constraints routine',
     .                   ' must be called constr'
              endif
            endif

          elseif (name1 .eq. 'position') then

c           Check syntax
            if (nr.ne.0 .and. nr.ne.3) then
              if (Node.eq.0) then
                write(6,'(a,/,a)')
     .            'fixed: syntax ERROR in %block GeometryConstraints:',
     .            line(1:lastch)
              endif
              goto 20
            endif

c           Find constrained atoms
            if (nn .gt. 1) then
c             Atoms specified by range. Make list of them.
              name2 = names(lch(1)+1:lch(2))
              name3 = names(lch(2)+1:lch(3))
              if (nn.eq.4) name4 = names(lch(3)+1:lch(4))
              if (name2.eq.'from' .and. name3.eq.'to') then
                ia1 = integs(1)
                ia2 = integs(2)
                if (ia1 .lt. 0) ia1 = na + ia1 + 1
                if (ia2 .lt. 0) ia2 = na + ia2 + 1
                if (nn.eq.4 .and. name4.eq.'step') then
                  ia3 = abs(integs(3))
                else
                  ia3 = 1
                endif
                ni = 0
                do ia = min(ia1,ia2),max(ia1,ia2),ia3
                  ni = ni + 1
                  integs(ni) = ia
                enddo
              else
                if (Node.eq.0) then
                  write(6,'(a,/,a)')
     .            'fixed: syntax ERROR in %block GeometryConstraints:',
     .            line(1:lastch)
                endif
                goto 20
              endif
            endif

c           Store position constraints
            do i = 1,ni
              ia = integs(i)
              if (ia .lt. 0) ia = na + ia + 1
              if (nr .eq. 0) then
c               Make one constraint for each cartesian coordinate
                do ix = 1,3
                  nc = nc + 1
                  if (nc.eq.maxc) then
                    maxcnew = nc + nint(0.1*nc)
C
                    allocate(cbuffer(maxc))
                    cbuffer(1:maxc) = ctype(1:maxc)
                    deallocate(ctype)
                    allocate(ctype(maxcnew))
                    ctype(1:maxc) = cbuffer(1:maxc)
                    deallocate(cbuffer)
C
                    allocate(ibuffer(maxc))
                    call memory('A','I',maxc,'fixed')
                    ibuffer(1:maxc) = iac(1:maxc)
                    call memory('D','I',size(iac),'fixed')
                    deallocate(iac)
                    allocate(iac(maxcnew))
                    call memory('A','I',maxcnew,'fixed')
                    iac(1:maxc) = ibuffer(1:maxc)
                    call memory('D','I',size(ibuffer),'fixed')
                    deallocate(ibuffer)
C
                    allocate(dpbuffer2(3,maxc))
                    call memory('A','D',maxc,'fixed')
                    dpbuffer2(1:3,1:maxc) = xc(1:3,1:maxc)
                    call memory('D','D',size(xc),'fixed')
                    deallocate(xc)
                    allocate(xc(3,maxcnew))
                    call memory('A','D',3*maxcnew,'fixed')
                    xc(1:3,1:maxc) = dpbuffer2(1:3,1:maxc)
                    call memory('D','D',size(dpbuffer2),'fixed')
                    deallocate(dpbuffer2)
C
                    maxc = maxcnew
                  endif
                  ctype(nc) = 'position'
                  iac(nc) = ia
                  do jx = 1,3
                    xc(jx,nc) = 0.0d0
                  enddo
                  xc(ix,nc) = 1.0d0
                enddo
              elseif (nr .eq. 3) then
c               Make only one constraint in the specified direction
                nc = nc + 1
                if (nc.eq.maxc) then
                  maxcnew = nc + nint(0.1*nc)
C
                  allocate(cbuffer(maxc))
                  cbuffer(1:maxc) = ctype(1:maxc)
                  deallocate(ctype)
                  allocate(ctype(maxcnew))
                  ctype(1:maxc) = cbuffer(1:maxc)
                  deallocate(cbuffer)
C
                  allocate(ibuffer(maxc))
                  call memory('A','I',maxc,'fixed')
                  ibuffer(1:maxc) = iac(1:maxc)
                  call memory('D','I',size(iac),'fixed')
                  deallocate(iac)
                  allocate(iac(maxcnew))
                  call memory('A','I',maxcnew,'fixed')
                  iac(1:maxc) = ibuffer(1:maxc)
                  call memory('D','I',size(ibuffer),'fixed')
                  deallocate(ibuffer)
C
                  allocate(dpbuffer2(3,maxc))
                  call memory('A','D',maxc,'fixed')
                  dpbuffer2(1:3,1:maxc) = xc(1:3,1:maxc)
                  call memory('D','D',size(xc),'fixed')
                  deallocate(xc)
                  allocate(xc(3,maxcnew))
                  call memory('A','D',3*maxcnew,'fixed')
                  xc(1:3,1:maxc) = dpbuffer2(1:3,1:maxc)
                  call memory('D','D',size(dpbuffer2),'fixed')
                  deallocate(dpbuffer2)
C
                  maxc = maxcnew
                endif
                ctype(nc) = 'position'
                iac(nc) = ia
                xnorm = sqrt(reals(1)**2 + reals(2)**2 + reals(3)**2)
                do ix = 1,3
                  xc(ix,nc) = reals(ix) / xnorm
                enddo
              endif
            enddo

          else
            if (Node.eq.0) then
              write(6,*) 'fixed: ERROR: sorry, constraint type ',
     .                  name1, ' not implemented yet'
            endif
          endif
   20   continue
   30   continue
        frstme = .false.
      endif

c Copy stress and forces to output arrays
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

c Apply constraints
      do ic = 1,nc
        if (ctype(ic) .eq. 'routine') then
          call constr( cell, na, isa, amass, xa, cstress, cfa )
        elseif (ctype(ic) .eq. 'position') then
          ia = iac(ic)
          fxc = dot( cfa(1,ia), xc(1,ic), 3 )
          do ix = 1,3
            cfa(ix,ia) = cfa(ix,ia) - fxc * xc(ix,ic)
          enddo
        endif
      enddo

      end


