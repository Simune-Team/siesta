C $Id: diag2g.f,v 1.3 1999/02/23 12:05:23 wdpgaara Exp $

      subroutine diag2g( no, maxo, maxno, maxnd, 
     .                   numh, listh, numd, listd, H, S,
     .                   getD, qtot, temp, e1, e2,
     .                   eo, qo, Dnew, Enew, ef,
     .                   Haux, Saux, psi, aux )
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices.
C This version is for non-colinear spin at gamma point.
C Writen by J.Soler, May and August 1998.
C **************************** INPUT **********************************
C integer no                  : Number of basis orbitals
C integer maxo                : Maximum number of basis  orbitals
C integer maxno               : Maximum number of orbitals interacting  
C                               with any orbital
C integer maxnd               : Maximum number of nonzero elements of 
C                               each row of density matrix
C integer numh(no)            : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listh(maxno,no)     : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer numd(no)            : Number of nonzero elements of each row 
C                               ofdensity matrix
C integer listd(maxnd,no)     : Nonzero density-matrix element column 
C                               indexes for each matrix row
C real*8  H(maxno,maxo,4) : Hamiltonian in sparse form
C real*8  S(maxno,maxo)       : Overlap in sparse form
C logical getD                : Find occupations and density matrices?
C real*8  qtot                : Number of electrons in unit cell
C real*8  temp                : Electronic temperature 
C real*8  e1, e2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
C *************************** OUTPUT **********************************
C real*8 eo(maxo*2)         : Eigenvalues
C real*8 qo(maxo*2)         : Occupations of eigenstates
C real*8 Dnew(maxnd,maxo,4) : Output Density Matrix
C real*8 Enew(maxnd,maxo,4) : Output Energy-Density Matrix
C real*8 ef                     : Fermi energy
C *************************** AUXILIARY *******************************
C real*8 Haux(2,2,no,2,no): Auxiliary space for the hamiltonian matrix
C real*8 Saux(2,2,no,2,no): Auxiliary space for the overlap matrix
C real*8 psi(2,2,no,2*no) : Auxiliary space for the eigenvectors
C real*8 aux(2*no*5)      : Extra auxiliary space
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C temp and H must be in the same energy units.
C eo, Enew and ef returned in the units of H.
C *********************************************************************
      implicit none

      integer
     .  maxo, maxnd, maxno, no

      integer 
     .  listh(maxno,no), numh(no),
     .  listd(maxnd,no), numd(no)

      double precision
     .  Dnew(maxnd,maxo,4),
     .  e1, e2, ef, Enew(maxnd,maxo,4), eo(maxo*2),
     .  H(maxno,maxo,4), qo(maxo*2), qtot,
     .  S(maxno,maxo), stepf, temp
     
      double precision
     .  aux(2*no*5), Haux(2,2,no,2,no),
     .  psi(2,2,no,2*no), Saux(2,2,no,2,no)

      logical
     .  getD

      external
     .  cdiag, fermid, stepf

C  Internal variables .............................................
      integer           i, ie, io, is, ispin, j, jo, js
      double precision  ee, pipj, qe, t
C  ....................

C Initialize Hamiltonian and overlap matrices in full format .....
C Index i is for real/imag parts
C Indices is and js are for spin components
C Indices iuo and juo are for orbital components:
C Haux(i,js,juo,is,iuo) = <js,juo|H|is,iuo>
      do io = 1,no
        do is = 1,2
          do jo = 1,no
            do js = 1,2
              do i = 1,2
                Saux(i,js,jo,is,io) = 0.d0
                Haux(i,js,jo,is,io) = 0.d0
              enddo
            enddo
          enddo
        enddo
      enddo
C  ....................

C Copy S,H matrices from sparse to full format ....................
C Convention: ispin=1 => H11, ispin=2 => H22, ispin=3 => Real(H12),
C             ispin=4 => Imag(H12)
      do io = 1,no
        do j = 1,numh(io)
          jo = listh(j,io)
          Saux(1,1,jo,1,io) =  S(j,io)
          Saux(1,2,jo,2,io) =  S(j,io)
          Haux(1,1,jo,1,io) =  H(j,io,1)
          Haux(1,2,jo,2,io) =  H(j,io,2)
          Haux(1,1,jo,2,io) =  H(j,io,3)
          Haux(1,2,jo,1,io) =  H(j,io,3)
          Haux(2,1,jo,2,io) = -H(j,io,4)
          Haux(2,2,jo,1,io) =  H(j,io,4)
        enddo
      enddo
C ....................

C Symmetrize S and H matrices .....................................
*     do io = 1,no
*       do is = 1,2
*         do jo = 1,io-1
*           do js = 1,2
*             do i = 1,2
*               Saux(i,js,jo,is,io) = ( Saux(i,js,jo,is,io) +
*    .                                  Saux(i,js,io,is,jo) ) / 2
*               Saux(i,js,io,is,jo) =   Saux(i,js,jo,is,io)
*               Haux(i,js,jo,is,io) = ( Haux(i,js,jo,is,io) +
*    .                                  Haux(i,js,io,is,jo) ) / 2
*               Haux(i,js,io,is,jo) =   Haux(i,js,jo,is,io)
*             enddo
*           enddo
*         enddo
*       enddo
*     enddo
C ....................

C Solve the eigenvalue problem .......................................
      call cdiag( Haux, 2*no, Saux, 2*no, 2*no,
     .            eo, psi, 2*no, aux )
C ....................

C Check if we are done ................................................
      if (.not.getD) return
C ....................

C Find new Fermi energy and occupation weights ........................
      call fermid( 2, 2, 1, 1.d0, no, no, eo, 
     .             temp, qtot, qo, ef )
C ....................

*     write(6,'(/,a,/,(10f7.2))') 'diag2g: eo =', eo
*     write(6,'(/,a,/,(10f7.2))') 'diag2g: qo =', qo

C Find weights for local density of states ............................
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.d-6 )
        do io = 1,no*2
          qo(io) =  ( stepf( (eo(io)-e2)/t ) -
     .                stepf( (eo(io)-e1)/t ) ) / 2
        enddo
      endif
C ....................
      
c New density and energy-density matrices of unit-cell orbitals .......
      do ispin = 1,4
        do io = 1,no
          do j = 1,numd(io)
            Dnew(j,io,ispin) = 0.d0
            Enew(j,io,ispin) = 0.d0
          enddo
        enddo
      enddo

      do ie = 1,no*2
        qe = qo(ie)
        ee = qo(ie) * eo(ie)
        do io = 1,no
          do j = 1,numd(io)
            jo = listd(j,io)
            pipj = psi(1,1,io,ie) * psi(1,1,jo,ie) +
     .             psi(2,1,io,ie) * psi(2,1,jo,ie)
            Dnew(j,io,1) = Dnew(j,io,1) + qe * pipj
            Enew(j,io,1) = Enew(j,io,1) + ee * pipj
            pipj = psi(1,2,io,ie) * psi(1,2,jo,ie) +
     .             psi(2,2,io,ie) * psi(2,2,jo,ie)
            Dnew(j,io,2) = Dnew(j,io,2) + qe * pipj
            Enew(j,io,2) = Enew(j,io,2) + ee * pipj
            pipj = psi(1,1,io,ie) * psi(1,2,jo,ie) +
     .             psi(2,1,io,ie) * psi(2,2,jo,ie)
            Dnew(j,io,3) = Dnew(j,io,3) + qe * pipj
            Enew(j,io,3) = Enew(j,io,3) + ee * pipj
            pipj = psi(1,1,io,ie) * psi(2,2,jo,ie) -
     .             psi(2,1,io,ie) * psi(1,2,jo,ie)
            Dnew(j,io,4) = Dnew(j,io,4) + qe * pipj
            Enew(j,io,4) = Enew(j,io,4) + ee * pipj
          enddo
        enddo
      enddo
C ....................

      end


