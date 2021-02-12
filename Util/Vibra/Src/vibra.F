! 
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.

!#define SANKEY_DIAG     ! Do this in the makefile
!#define VIBRA_DEBUG
      program vibra

c *********************************************************************
c Calculation of vibrational modes for clusters, linear chains, slabs 
c and 3D xtals.
c 
c Uses the FDF (Flexible Data Format) package (version 0.66.1.5)
c of J.M.Soler and A. Garcia, 
c
c Written by P.Ordejon, August'98
c
c Modified to add Infra-red intensity calculation by J.D. Gale, July'06
c
c **********************************************************************

      use fdf
      use units, only: Pi
      use precision, only: dp

      implicit none

c Internal variables ...

      character(len=150) ::
     .  filein, fileout, fname, 
     .  slabel, sname

      integer 
     .  i, i1, i2, iatom, iunit,
     .  ik, iunit2, iunit3, ind, ix, j,
     .  ij, ii, iq, jj, jx, icall, 
     .  lx, ly, lz, lxmax, lymax, lzmax, 
     .  llx, lly, llz,
     .  natoms, ncells, nnat, ndof

      ! imass(maxa), imasssc(maxasc)
      integer, allocatable :: imass(:), imasssc(:)
#ifndef SANKEY_DIAG
      character(len=1) :: heev_mode
#endif

      integer maxlin, maxk, maxnq
      parameter (maxlin =  1000)
      parameter (maxk   =  5000)
!     maxnq is the maximum number of different q-points for
!     the Wiegner-Sietz cells. Due to the loop it will never be more
!     than 27 (3 x -1:1)
      parameter (maxnq = 27)
     
      integer 
     .  in, neq, nk, nlines, lastk(maxlin)

      character label(maxlin)*8

      real(dp) :: kpoint(3,maxk)

      real(dp) :: dx, alat, alp, blp, clp, alplp, betlp, gamlp, xxx

      real(dp) :: zpe, planck
      logical :: eigen, intensity

      ! ek(maxd,maxk)
      real(dp), allocatable :: ek(:,:)
      ! phi0(3,maxa,3,maxa,-maxx:maxx,-maxy:maxy,-maxz:maxz),
      real(dp), allocatable :: phi0(:,:,:,:, :,:,:)
      ! phi(3,maxa,3,maxa,-maxx:maxx,-maxy:maxy,-maxz:maxz)
      real(dp), allocatable :: phi(:,:,:,:, :,:,:)
      ! IRinten(3*maxa), BornQ(3,3,maxa)
      real(dp), allocatable :: IRinten(:), BornQ(:,:,:)

      ! b(3,maxa), xa(3,maxasc), xmass(maxa)
      real(dp), allocatable :: b(:,:), xa(:,:), xmass(:)
      real(dp) :: cell(3,3), r(3), scell(3,3)

      real(dp) :: correct, dmin, q(3), qr(maxnq), r2, rl(3), rmass

c Correction terms to satisfy translational modes
      real(dp), allocatable :: zero(:,:,:)
      real(dp) :: zeroo(3,3)

c Work space for diagonalization.
      complex(dp) :: phase, IRtrm, vecmw(3)
      ! dc(maxd,maxd)
      complex(dp), allocatable :: dc(:,:)
#ifdef SANKEY_DIAG
      ! work(maxd),work2(2,maxd)
      real(dp), allocatable :: work(:),work2(:,:)
      ! dd(maxd,maxd),zr(maxd,maxd),zi(maxd,maxd)
      real(dp), allocatable :: dd(:,:), zr(:,:), zi(:,:)
#else
      ! work(3*maxd),work2(4,maxd)
      real(dp), allocatable :: work(:),work2(:,:)
#endif
      real(dp), allocatable :: omega(:)

c Conversion factor from sqrt(K/M) in eV and Ang to cm**-1 is 519.6
      real(dp), parameter :: xmagic = 519.6_dp

      data nk / 0 /

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

c ZPE contains the zero point energy of the system
      zpe = 0.0d0

c Planck contains planck's constant in eV.cm
      planck = 6.62618d-5*2.997925d0/1.602189d0
c ...
     
C ****************** READ DATA FROM FDF FILE *********************

c Set up fdf ...
      filein = 'stdin'
      fileout = 'out.fdf'
      call fdf_init(filein,fileout)
c ...

c Defile Name of the system ...
      sname = fdf_get('SystemName'," ")
      write(6,'(a,a)')
     . 'redata: System Name                      = ',sname
c ...

c Defile System Label (short name to label files) ...
      slabel = fdf_get('SystemLabel', "vibra")
      write(6,'(a,a)')
     . 'redata: System Label                     = ',slabel
c ...

c Read Number of Atoms in Unit cell ...
      natoms = fdf_get('NumberOfAtoms', 0)
      if (natoms .le. 0) then
        write(6,'(a)') 
     . 'ERROR: Number of atoms must be larger than zero.'
        write(6,'(a)') 
     . '       You MUST specify NumberOfatoms in fdf file.'
        stop
      endif
      write(6,'(a,i5)') 
     . 'Number of Atoms                  = ',natoms
c ...

      ! Number of degrees of freedom
      ndof = natoms * 3

c Lattice constant of unit cell...
      alat = fdf_get('LatticeConstant',0.d0,'Bohr')
      if (alat .eq. 0.d0) then
        write(6,'(a)') 
     . 'ERROR: No valid lattice constant specified.'
        write(6,'(a)') 
     . '       You MUST specify LatticeConstant in fdf file.'
        stop
      endif
      write(6,'(a,f10.5,a)') 'Lattice Constant    = ',alat,'  Bohr'
c ...

c Lattice vectors of unit cell...
      if ( fdf_block('LatticeParameters',bfdf) .and.
     .     fdf_block('LatticeVectors',bfdf) ) then
         write(6,'(2a)')'ERROR: Lattice vectors doubly ',
     .     'specified: by LatticeVectors and by LatticeParameters.' 
         stop 
      endif

      if ( fdf_block('LatticeParameters',bfdf) ) then
         if (.not. fdf_bline(bfdf, pline))
     $        call die("No LatticeParameters")
         if (.not. (fdf_bmatch(pline, 'vvvvvv') )) then
            call die ('LatticeParameters: Error in syntax')
         endif
         alp = fdf_bvalues(pline,1)
         blp = fdf_bvalues(pline,2)
         clp = fdf_bvalues(pline,3)
         alplp = fdf_bvalues(pline,4)
         betlp = fdf_bvalues(pline,5)
         gamlp = fdf_bvalues(pline,6)
         write(6,'(a)')
     .    'Lattice Parameters (units of Lattice Constant) ='
         write(6,'(a,3f10.5,3f9.3)')
     .    '    ',alp,blp,clp,alplp,betlp,gamlp
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
          cell(1,i) = fdf_bvalues(pline,1)
          cell(2,i) = fdf_bvalues(pline,2)
          cell(3,i) = fdf_bvalues(pline,3)
        enddo
      else
        do i = 1,3
          do j  = 1,3
            cell(i,j) = 0.d0
          enddo
          cell(i,i) = 1.d0
        enddo
      endif
      write(6,'(a)') 
     .   'Lattice vectors (in units of Lattice Constant) ='
      do i = 1,3
        write(6,'(a,3f10.5)')
     .   '        ',(cell(j,i), j=1,3)
      enddo
c ...

c Multiply cell vectors by lattice constant ...
      do i = 1,3
        do ix = 1,3
          cell(ix,i) = alat * cell(ix,i)
        enddo
      enddo
      write(6,'(a)') 
     .   'Lattice vectors (in Bohr) ='
      do i = 1,3
        write(6,'(a,3f10.5)') '        ',(cell(j,i), j=1,3)
      enddo
c ...

c     Define number of unit cells in the supercell ...
      lxmax = fdf_get('SuperCell_1', 0)
      lymax = fdf_get('SuperCell_2', 0)
      lzmax = fdf_get('SuperCell_3', 0)
      ncells = (2*lxmax+1)*(2*lymax+1)*(2*lzmax+1)
      write(6,'(a,i5)') 'lxmax    = ',lxmax
      write(6,'(a,i5)') 'lymax    = ',lymax
      write(6,'(a,i5)') 'lzmax    = ',lzmax
      write(6,'(a,i5)') 'Number of unit cells in Supercell  = ',ncells
c ...

!     Now we are ready for allocation of matrices
      allocate(phi0(3,natoms,3,natoms,
     &    -lxmax:lxmax,-lymax:lymax,-lzmax:lzmax))
      allocate(phi(3,natoms,3,natoms,
     &    -lxmax:lxmax,-lymax:lymax,-lzmax:lzmax))
      allocate(IRinten(ndof),BornQ(3,3,natoms))
      allocate(zero(3,3,natoms))
      allocate(b(3,natoms), xmass(natoms))
      allocate(xa(3,ncells*natoms))
      allocate(imass(natoms),imasssc(ncells*natoms))

c Read atomic coordinates and species of unit cell...
      call recoor(cell,alat,b,imass,xmass,natoms)
c ...

c Determine q points to compute the phonon dispersion relations ...
      call klines(maxk, nk, nlines, lastk, label, kpoint)
c ...

c Determine whether IR intensities are to be computed
      intensity = fdf_get('Intensities', .false.)

c If intensities are requested then eigenvectors must be computed
      if ( intensity ) then
        eigen = .true.
      else
        eigen = fdf_get('Eigenvectors', .false.)
      end if

c Determine whether eigenvectors are computed, or only eigenvalues ...
      if (eigen) then
         icall = 3
#ifndef SANKEY_DIAG
         heev_mode = 'V'
#endif
         write(6,'(a,i5)') 'Eigenvectors =   True'
         write(6,'(a,i5)') 'Computing Eigenvalues and Eigenvectors'
      else
         icall = 1
#ifndef SANKEY_DIAG
         heev_mode = 'N'
#endif
         write(6,'(a,i5)') 'Eigenvectors =   False'
         write(6,'(a,i5)') 'Computing Eigenvalues only'
      endif

C *************** END READ DATA FROM FDF FILE ********************

c Build lattice vector of the supercell ...
      do ix=1,3
        scell(ix,1) = (2*lxmax+1)*cell(ix,1)
        scell(ix,2) = (2*lymax+1)*cell(ix,2)
        scell(ix,3) = (2*lzmax+1)*cell(ix,3)
      enddo
c ...

C Build atomic coordinates in the supercell ...
c loop over unit cells within supercell
      iatom=0
      do lx=-lxmax,lxmax
      do ly=-lymax,lymax
      do lz=-lzmax,lzmax
        r(1) = lx*cell(1,1) + ly*cell(1,2) + lz*cell(1,3)
        r(2) = lx*cell(2,1) + ly*cell(2,2) + lz*cell(2,3)
        r(3) = lx*cell(3,1) + ly*cell(3,2) + lz*cell(3,3)
c loop over atoms in unit cell
        do i=1,natoms
          iatom=iatom+1
          imasssc(iatom)=imass(i)
          do ix=1,3
            xa(ix,iatom) = b(ix,i) + r(ix)
          enddo
        enddo
      enddo
      enddo
      enddo

      nnat = natoms * (2*lxmax+1) * (2*lymax+1) * (2*lzmax+1)
      if (iatom .ne. nnat) stop 'Error computing number of atoms'
c ...

c  Determine the indices of the atoms in the central (lx=ly=lz=0) cell
c  (those that need to be displaced to calculate the force constants)  ...
      i1 = natoms * (4*lxmax*lymax*lzmax + 
     .               2*(lxmax*lymax + lxmax*lzmax + lymax*lzmax) +
     .               lxmax + lymax + lzmax) + 1
      i2 = i1 + natoms - 1
c ...

c Read Force Constants Matrix ...
      call io_assign(iunit2)
      fname = fdf_get('Vibra.FC', trim(slabel)//'.FC')
      open(iunit2,file=fname,status='old')
      read(iunit2,*)
! Loop displaced atoms
      do j=1,natoms
!  displacement directions x, y, z
        do ij=1,3
c Negative displacements
          do lx=-lxmax,lxmax
          do ly=-lymax,lymax
          do lz=-lzmax,lzmax
            do i=1,natoms
              !read(iunit2,*) (pn(ii,i,lx,ly,lz),ii=1,3)
              read(iunit2,*) (phi0(ii,i,ij,j,lx,ly,lz),ii=1,3)
            enddo
          enddo
          enddo
          enddo
c Positive displacements
          do lx=-lxmax,lxmax
          do ly=-lymax,lymax
          do lz=-lzmax,lzmax
            do i=1,natoms
              !read(iunit2,*) (pp(ii,i,lx,ly,lz),ii=1,3)
              read(iunit2,*) q(:)
              do ii = 1, 3
!               Average displacements
!               Dimensions are:
! [xyz, atoms, displacement direction, displacement atoms, (supercell)]
                phi0(ii,i,ij,j,lx,ly,lz) =
     &              (phi0(ii,i,ij,j,lx,ly,lz) + q(ii))*0.5
              end do
            enddo
          enddo
          enddo
          enddo
        enddo
      enddo
      call io_close(iunit2)
c ...

c If intensities are required then read Born effective charges
      if (intensity) then
        call io_assign(iunit2)
        fname = trim(slabel)//'.BC'
        open(iunit2,file=fname,status='old')
        read(iunit2,*)
        do j = 1,natoms
          do ii = 1,3
            read(iunit2,*) (BornQ(ij,ii,j),ij=1,3)
          enddo
        enddo
        call io_close(iunit2)
      endif

c =================================================================
c Now form phibar(1). This force constant matrix is hermitian.

      do i=1,natoms
        do j=1,natoms
          do ii=1,3
            do ij=1,3
              do lx=-lxmax,lxmax
              do ly=-lymax,lymax
              do lz=-lzmax,lzmax
c                phi(ii,i,ij,j,lx,ly,lz) = phi0(ii,i,ij,j,lx,ly,lz)
                phi(ii,i,ij,j,lx,ly,lz) = 
     .            0.5d0 * (phi0(ii,i,ij,j,lx,ly,lz) 
     .            + phi0(ij,j,ii,i,-lx,-ly,-lz))
              enddo
              enddo
              enddo
            enddo
          enddo
        enddo
      enddo

!     Deallocate phi0 array
      deallocate(phi0)

c     goto 100
c =================================================================
c Now form zero(ix,jx,j), defined as zijbeta in the notes.

      do j=1,natoms
        do ii=1,3
          do ij=1,3
            zero(ii,ij,j) = 0.0d0
            do i=1,natoms
              do lx=-lxmax,lxmax
              do ly=-lymax,lymax
              do lz=-lzmax,lzmax
                zero(ii,ij,j) = zero(ii,ij,j) + phi(ii,i,ij,j,lx,ly,lz)
              enddo
              enddo
              enddo
            enddo
            zero(ii,ij,j) = zero(ii,ij,j) / (natoms*ncells)
          enddo
        enddo
      enddo
         

c =================================================================
c     Now form zeroo(ix,jx), the sum ofer beta of z(i,j,beta)

      do ii=1,3
        do ij=1,3
          zeroo(ii,ij) = 0.0d0
          do j=1,natoms
            zeroo(ii,ij) = zeroo (ii,ij) + zero(ii,ij,j)
          enddo
          zeroo(ii,ij) = zeroo(ii,ij) / natoms
        enddo
      enddo

#ifdef VIBRA_DEBUG
      print '(a)', "DEBUG"
      print '(a)', "zero sum of atoms due to displacements"
      do j = 1, natoms
        print '(a,i0)', "displaced atom ", j
        print '(a)', "   displ dir x"
        print '(tr1,3(tr1,es16.7))', zero(:,1,j)
        print '(a)', "   displ dir y"
        print '(tr1,3(tr1,es16.7))', zero(:,2,j)
        print '(a)', "   displ dir z"
        print '(tr1,3(tr1,es16.7))', zero(:,3,j)
      end do
      print '(a)', "zeroo sum zero"
      print '(a)', "   displ dir x"
      print '(tr1,3(tr1,es16.7))', zeroo(:,1)
      print '(a)', "   displ dir y"
      print '(tr1,3(tr1,es16.7))', zeroo(:,2)
      print '(a)', "   displ dir z"
      print '(tr1,3(tr1,es16.7))', zeroo(:,3)

      print '(a)', "before symmetrization"
      do i=1,natoms
        print '(a,i0)', "displaced atom ", i
        do ii = 1 , 3
          select case ( ii )
          case ( 1 )
            print '(a)', "   displ dir x"
          case ( 2 )
            print '(a)', "   displ dir y"
          case ( 3 )
            print '(a)', "   displ dir z"
          end select
          do lx=-lxmax,lxmax
            do ly=-lymax,lymax
              do lz=-lzmax,lzmax
                do j=1,natoms
                  print '(tr1,3(tr1,es16.7))', phi(:,j,ii,i,lx,ly,lz)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

#endif

c =================================================================
c Now form phibar(2)

      do i=1,natoms
        do j=1,natoms
          do ii=1,3
            do ij=1,3
              correct = (zeroo(ii,ij)+zeroo(ij,ii))/2.0d0 -
     .                  (zero(ii,ij,j)+zero(ij,ii,i))
              do lx=-lxmax,lxmax
              do ly=-lymax,lymax
              do lz=-lzmax,lzmax
                phi(ii,i,ij,j,lx,ly,lz) = phi(ii,i,ij,j,lx,ly,lz) +
     .                                    correct
              enddo
              enddo
              enddo
            enddo
          enddo
        enddo
      enddo

!     Deallocate zero
      deallocate(zero)

#ifdef VIBRA_DEBUG
      print '(a)', "DEBUG"
      print '(a)', "after symmetrization"
      do i=1,natoms
        print '(a,i0)', "displaced atom ", i
        do ii = 1 , 3
          select case ( ii )
          case ( 1 )
            print '(a)', "   displ dir x"
          case ( 2 )
            print '(a)', "   displ dir y"
          case ( 3 )
            print '(a)', "   displ dir z"
          end select
          do lx=-lxmax,lxmax
            do ly=-lymax,lymax
              do lz=-lzmax,lzmax
                do j=1,natoms
                  print '(tr1,3(tr1,es16.7))', phi(:,j,ii,i,lx,ly,lz)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      print '(a)', "atomic masses"
      do j = 1, natoms
        print '(tr1,i3,tr1,es16.7)', j,xmass(j)
      enddo

#endif

100   continue

c =================================================================
c Now finally compute the dynamical matrix. We store in real-hermitian 
c form because hermdp.f wants input this way.

      allocate(dc(ndof,ndof), omega(ndof))
#ifdef SANKEY_DIAG
      allocate(dd(ndof,ndof), zr(ndof,ndof), zi(ndof,ndof))
      allocate(work(ndof), work2(2,ndof))
#else
      ! work2 should optimally be (NB + 1) * ndof, but minimally
      ! 2*ndof (in complex form)
      allocate(work(3*ndof), work2(2,16*ndof))
#endif
      allocate(ek(ndof,maxk))
c Loop over k points 

      do ik=1,nk
#ifdef VIBRA_DEBUG
        print "(a,i0)", "running ik == ", ik
#endif

        q(1) = kpoint(1,ik)
        q(2) = kpoint(2,ik)
        q(3) = kpoint(3,ik)

        do jx=1,ndof
        do ix=1,ndof
            dc(ix,jx)= cmplx(0.0_dp,0.0_dp, dp)
        enddo
        enddo

        do i=1,natoms
        do j=1,natoms
          do lx=-lxmax,lxmax
          do ly=-lymax,lymax
          do lz=-lzmax,lzmax
            do jj=1,3
              r(jj) = b(jj,i)-b(jj,j)+
     .                 lx*cell(jj,1) + ly*cell(jj,2) + lz*cell(jj,3)
            enddo
! I *guess* this is the Wigner-Sietz reduction due
! to supercell interactions. But I do not have the notes.
            dmin=1.e8
            do llx=-1,1
            do lly=-1,1
            do llz=-1,1
              R2=0.
              do jj=1,3
                rl(jj)=llx*scell(jj,1)+lly*scell(jj,2)+llz*scell(jj,3)
                R2=R2+(rl(jj)+r(jj))**2
              enddo
              if (abs(R2-dmin) .gt. 0.0001) then
                if (R2 .lt. dmin) then
                  neq = 1
                  dmin = R2
                  qr(1)=q(1)*(rl(1)+r(1))+
     .                  q(2)*(rl(2)+r(2))+
     .                  q(3)*(rl(3)+r(3))
                  goto 10
                endif
              endif
              if (abs(R2-dmin) .lt. 0.0001) then
                neq = neq+1
                qr(neq)=q(1)*(rl(1)+r(1))+
     .                  q(2)*(rl(2)+r(2))+
     .                  q(3)*(rl(3)+r(3))
              endif
10            continue
            enddo
            enddo
            enddo
#ifdef VIBRA_DEBUG
      print "(a,'(',i2,',',i2,')',tr2,'[',2(i0,','),i0,']',tr2,i4)",
     &          "(i,j) [x,y,z] neq     ",i, j, lx, ly, lz, neq
#endif
c           qr = q(1)*r(1) + q(2)*r(2) + q(3)*r(3)
            do in = 1,neq
              phase = cmplx(cos(qr(in)), sin(qr(in)), dp)/neq
#ifdef VIBRA_DEBUG
              print "(tr2,2(tr1,es16.7))", phase
#endif
              do ij=1,3
              jx = (j-1)*3+ij
              do ii=1,3
                ix = (i-1)*3+ii
                dc(ix,jx)=dc(ix,jx)+phi(ii,i,ij,j,lx,ly,lz)*phase
              enddo
              enddo
            enddo
          enddo
          enddo
          enddo
        enddo
        enddo

        do j = 1, natoms
          do ij = 1, 3
            jx = (j-1)*3+ij
            do i = 1, natoms
              do ii = 1, 3
                ix = (i-1)*3+ii
                dc(ix,jx) = dc(ix,jx) / sqrt(xmass(i)*xmass(j))
              end do
            end do
          end do
        end do

#ifdef VIBRA_DEBUG
        print '(a)', "DEBUG"
        print '(a,tr1,3(tr1,f10.3))', "dc @ ", q(:)
        do ix=1,ndof
          print '(tr1,2000(tr1,es16.8))', real(dc(1:ndof,ix), dp)
        enddo
        print '(a)', ""
#endif

#ifdef SANKEY_DIAG
        do jx=ix,ndof
          do ix=1,ndof
            dd(ix,jx)=aimag(dc(jx,ix))
            dd(jx,ix)=real(dc(jx,ix), dp)
          enddo
        enddo
#endif

#ifdef SANKEY_DIAG
        call hermdp(dd,zr,zi,omega,ndof,ndof,work,work2,icall)
#else
        call zheev(heev_mode,'L',ndof,dc,ndof,omega,work2,2*ndof,work,i)
#endif

c the eigenvalues omega are actually omega**2.
c zr(i,j)=i'th component of j'th eigenvector. zr is the real part, and
c zi is the imaginary part. Our matrix is real, so we should get zero inmaginary
c part.!!
c =================================================================
c Compute the IR intensities if requested
      if (intensity) then

c Loop over modes
        do i = 1,ndof
          IRinten(i) = 0._dp
          do ij = 1,3
            IRtrm = cmplx(0._dp,0._dp, dp)
            ind = - 3
            do j = 1,natoms
              ind = ind + 3

c Mass weight eigenvectors
              rmass = 1._dp/sqrt(xmass(j))
              do ii = 1,3
#ifdef SANKEY_DIAG
                vecmw(ii) = cmplx(zr(ind+ii,i),zi(ind+ii,i), dp)
                vecmw(ii) = vecmw(ii)*rmass
#else
                vecmw(ii) = dc(ind+ii,i)*rmass
#endif
              enddo

c Multiply eigenvectors by Born effective charges
              do ii = 1,3
                IRtrm = IRtrm + BornQ(ii,ij,j)*vecmw(ii)
              enddo

            enddo
            IRinten(i) = IRinten(i) + IRtrm*IRtrm
          enddo
        enddo

      endif

c =================================================================
c write out eigenvalues.
        do i=1,ndof
c convert to cm**-1. the conversion factor is xmagic (511**2)
          omega(i)=xmagic*xmagic*omega(i)
          if(omega(i).gt.0.0d0)then
            omega(i)=sqrt(omega(i))
          else
            omega(i)=-omega(i)
            omega(i)=-sqrt(omega(i))
c        write(*,*)' Caution: omega**2 .lt.0.0 .....',
c     1 ' sqrt(abs(omega2))=',omega(i)
          end if
          write(*,*)' eigenvalue #',i,' omega=',omega(i)

c Add zero point energy contribution for positive modes
          if (omega(i).gt.0.0d0) then
            zpe = zpe + omega(i)
          endif

          ek(i,ik) = omega(i)
        end do

Cc =============================================================
c write the eigenvalues and eigenvectors to output data file.
        if (icall .eq. 3) then
          call io_assign(iunit3)
          fname = trim(slabel)//'.vectors'
          open(iunit3,file=fname,status='unknown',position='append')
          if (ik.eq.1) then
            rewind(iunit3)
            write(6,'(/,a)')' Writing eigenvalues and eigenvectors'
            write(6,'(2a,/)')' to output file ', fname
          endif
          write(iunit3,'(/,a,3f12.6)') 'k            = ',
     .                                  q(1),q(2),q(3)
          do i=1,ndof
            write(iunit3,'(a,i5)')     'Eigenvector  = ',i
            write(iunit3,'(a,f12.6)')  'Frequency    = ',omega(i)
            if (intensity) then
              write(iunit3,'(a,f12.6)')  'IR Intensity = ',IRinten(i)
            endif
            write(iunit3,'(a)')        'Eigenmode (real part)'
#ifdef SANKEY_DIAG
            write(iunit3,'(3e12.4)') (zr(j,i),j=1,ndof)
#else
            write(iunit3,'(3e12.4)') (real(dc(j,i),dp),j=1,ndof)
#endif
            write(iunit3,'(a)')        'Eigenmode (imaginary part)'
#ifdef SANKEY_DIAG
            write(iunit3,'(3e12.4)') (zi(j,i),j=1,ndof)
#else
            write(iunit3,'(3e12.4)') (aimag(dc(j,i)),j=1,ndof)
#endif
          end do
          call io_close(iunit3)
        endif
Cc =================================================================

      enddo

c Finish computation of zero point energy (averaged across k points) and output
      if (nk.gt.0) then
        zpe = 0.5d0*planck*zpe/dble(nk)
        write(*,'('' Zero point energy = '',f20.6,'' eV'')') zpe
      endif

c Write eigenvalues ...
      if (nk .gt. 0) then
        call outbands(0, 1, ndof, ndof, maxk, nk, nlines, lastk,
     .                label, kpoint, ek, 0.0d0)
      endif
C ...


   33 format(7(f10.4,2x))
Cc
      stop

      CONTAINS

      subroutine die(str)
      character(len=*), intent(in), optional:: str
      if (present(str)) then
         write(6,"(a)") str
         write(0,"(a)") str
      endif
      STOP
      end subroutine die


      end
