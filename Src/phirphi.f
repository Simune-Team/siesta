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
      subroutine phirphi(nua, na, nuo, no, scell, xa, rmaxo,
     .                   maxna, maxnh, lasto, iphorb, isa, 
     .                   numh, listhptr, listh, dk, 
     .                   jna, xij, r2ij, S)
C *********************************************************************
C Finds the matrix elements of the dk*r between basis
C orbitals, where dk is a given vector and r is the position operator.
C
C     S_a,b(R_a,R_b) = 0.5*( <phi_a(r-R_a)|dk*(r-R_b)|phi_b(r-R_b)>
C                      +      <phi_b(r-R_b)|dk*(r-R_a)|phi_a(r-R_a)> )
C
C Energies in Ry. Lengths in Bohr.
C Writen by DSP March 1999 ( Based in routine overfsm )
C **************************** INPUT **********************************
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C integer nuo              : Number of basis orbitals in unit cell
C integer no               : Number of basis orbitals 
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
C real*8  rmaxo            : Maximum cutoff for atomic orbitals
C integer maxna            : Maximum number of neighbours of any atom
C integer maxnh            : First dimension of listh
C integer lasto(0:na)      : Last orbital index of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom
C integer isa(na)          : Species index of each atom
C integer numh(nuo)        : Number of nonzero elements of each row
C                            of the overlap matrix
C integer listhptr(nuo)    : Pointer to the start of each row 
C                            of the overlap matrix
C integer listh(maxnh)     : Column indexes of the nonzero elements  
C                            of each row of the overlap matrix
C real*8  dk(3)            : Vector in k-space
C integer jna(maxna)       : Aux. space to find neighbours (indexes)
C real*8  xij(3,maxna)     : Aux. space to find neighbours (vectors)
C real*8  r2ij(maxna)      : Aux. space to find neighbours (distances)
C **************************** OUTPUT *********************************
C real*8  S(maxnh)         : Sparse overlap matrix
C *********************************************************************

      use precision
      use atmfuncs,     only : rcut
      use parallel,     only : Node, Nodes
      use parallelsubs, only : GlobalToLocalOrb

      implicit none

C Passed variables 
      integer
     .  maxna, maxnh, na, no, nua, nuo

      integer
     .  iphorb(no), isa(na), jna(maxna), lasto(0:na), 
     .  listh(maxnh), listhptr(nuo), numh(nuo)
      
      real(dp)
     .  scell(3,3), r2ij(maxna), rmaxo, dk(3),
     .  S(maxnh), xa(3,na), xij(3,maxna)

C Internal variables 
      integer
     .  ia, ind, iio, io, ioa, is, ix, 
     .  j, ja, jn, jo, joa, js, nnia

      real(dp)
     .  grSij(3), rij, Sij, xinv(3)

      real(dp), dimension(:), allocatable, save :: 
     .  Si

      real(dp), save :: tiny = 1.0d-9

      external
     .  neighb, memory
C ......................

C Initialize neighb subroutine 
      nnia = maxna
      call neighb( scell, 2.d0*rmaxo, na, xa, 0, 0,
     .             nnia, jna, xij, r2ij )

C Allocate local memory
      allocate(Si(no))
      call memory('A','D',no,'phirphi')

      do jo = 1,no
        Si(jo) = 0.0d0
      enddo
      do ia = 1,nua 
        nnia = maxna 
        is = isa(ia)
        call neighb( scell, 2.0d0*rmaxo, na, xa, ia, 0,
     .               nnia, jna, xij, r2ij )  
           
        do iio = lasto(ia-1)+1,lasto(ia)  
          call GlobalToLocalOrb(iio,Node,Nodes,io)
          if (io .gt. 0) then
            ioa = iphorb(iio)
            do jn = 1,nnia 
              do ix = 1,3
                xinv(ix) = - xij(ix,jn)
              enddo 
              ja = jna(jn)
              rij = sqrt( r2ij(jn) )
              do jo = lasto(ja-1)+1,lasto(ja)
                joa = iphorb(jo)
                js = isa(ja) 
                if (rcut(is,ioa)+rcut(js,joa) .gt. rij) then  

                  if (abs(dk(1)).gt.tiny) then
                    call matel('X', is, js, ioa, joa, xij(1,jn),
     .                          Sij, grSij ) 
                    Si(jo) = Si(jo) + 0.5d0*Sij*dk(1)  
 
                    call matel('X', js, is, joa, ioa, xinv,
     .                          Sij, grSij )
                    Si(jo) = Si(jo) + 0.5d0*Sij*dk(1)  
                  endif
                     
                  if (abs(dk(2)).gt.tiny) then
                    call matel('Y', is, js, ioa, joa, xij(1,jn),
     .                          Sij, grSij )
                    Si(jo) = Si(jo) + 0.5d0*Sij*dk(2) 
                
                    call matel('Y', js, is, joa, ioa, xinv,
     .                          Sij, grSij )
                    Si(jo) = Si(jo) + 0.5d0*Sij*dk(2)  
                  endif
 
                  if (abs(dk(3)).gt.tiny) then
                    call matel('Z', is, js, ioa, joa, xij(1,jn),
     .                          Sij, grSij )
                    Si(jo) = Si(jo) + 0.5d0*Sij*dk(3) 
 
                    call matel('Z', js, is, joa, ioa, xinv,
     .                          Sij, grSij )
                    Si(jo) = Si(jo) + 0.5d0*Sij*dk(3) 
                  endif
                endif
              enddo
            enddo
            do j = 1,numh(io)
              ind = listhptr(io)+j
              jo = listh(ind)
              S(ind) = Si(jo) 
              Si(jo) = 0.0d0 
            enddo
          endif
        enddo
      enddo

C Deallocate local memory
      call memory('D','D',size(Si),'phirphi')
      deallocate(Si)

      return
      end
