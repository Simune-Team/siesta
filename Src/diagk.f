C $Id: diagk.f,v 1.4 1999/05/18 17:17:32 ordejon Exp $

      subroutine diagk( nspin, nuo, no, maxspn, maxuo, maxno, maxnd, 
     .                  numh, listh, numd, listd, H, S,
     .                  getD, fixspin, qtot, qs, temp, e1, e2,
     .                  xij, indxuo, nk, kpoint, wk,
     .                  eo, qo, Dnew, Enew, ef, efs,
     .                  Haux, Saux, psi, Dk, Ek, aux )
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization). K-sampling version.
C Writen by J.Soler, August 1998.
C **************************** INPUT **********************************
C integer nspin               : Number of spin components (1 or 2)
C integer nuo                 : Number of basis orbitals in unit cell
C integer no                  : Number of basis orbitals in supercell
C integer maxspn              : Second dimension of eo and qo
C integer maxuo               : First dimension of eo, qo, last of xij
C                               Must be at least max(indxuo)
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
C real*8  H(maxno,maxuo,nspin): Hamiltonian in sparse form
C real*8  S(maxno,maxuo)      : Overlap in sparse form
C logical getD                : Find occupations and density matrices?
C logical fixspin             : Fix the spin of the system?
C real*8  qtot                : Number of electrons in unit cell
C real*8  qs(nspin)           : Number of electrons in unit cell for each 
C                               spin component (if fixed spin option is used)
C real*8  temp                : Electronic temperature 
C real*8  e1, e2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
C real*8  xij(3,maxno,maxuo)  : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C integer nk                  : Number of k points
C real*8  kpoint(3,nk)        : k point vectors
C real*8  wk(nk)              : k point weights (must sum one)
C *************************** OUTPUT **********************************
C real*8 eo(maxuo,maxspn,nk)     : Eigenvalues
C ******************** OUTPUT (only if getD=.true.) *******************
C real*8 qo(maxuo,maxspn,nk)     : Occupations of eigenstates
C real*8 Dnew(maxno,maxuo,nspin) : Output Density Matrix
C real*8 Enew(maxno,maxuo,nspin) : Output Energy-Density Matrix
C real*8 ef                      : Fermi energy
C real*8 efs(nspin)              : Fermi energy for each spin
C                                  (for fixed spin calculations)
C *************************** AUXILIARY *******************************
C real*8 Haux(2,nuo,nuo) : Auxiliary space for the hamiltonian matrix
C real*8 Saux(2,nuo,nuo) : Auxiliary space for the overlap matrix
C real*8 psi(2,nuo,nuo)  : Auxiliary space for the eigenvectors
C real*8 aux(2*nuo*5)    : Extra auxiliary space
C real*8 Dk(2,nuo,nuo)   : Aux. space that may be the same as Haux
C real*8 Ek(2,nuo,nuo)   : Aux. space that may be the same as Saux
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C temp and H must be in the same energy units.
C eo, Enew and ef returned in the units of H.
C *********************************************************************
      implicit          none
      integer           maxnd, maxno, maxspn, maxuo, nk, no,
     .                  nspin, nuo
      integer           indxuo(no), listh(maxno,no), numh(no),
     .                  listd(maxnd,no), numd(no)
      double precision  Dnew(maxnd,maxuo,nspin),
     .                  e1, e2, ef, efs(nspin), Enew(maxnd,maxuo,nspin),
     .                  eo(maxuo,maxspn,nk), H(maxno,maxuo,nspin),
     .                  kpoint(3,nk), qo(maxuo,maxspn,nk), qs(nspin), 
     .                  qtot, S(maxno,maxuo), stepf, temp, wk(nk),
     .                  xij(3,maxno,maxuo)
      double precision  Dk(2,nuo,nuo), Ek(2,nuo,nuo),
     .                  Haux(2,nuo,nuo), Saux(2,nuo,nuo),
     .                  psi(2,nuo,nuo), aux(2*nuo*5)
      logical           fixspin, getD
      external          cdiag, fermid, fermispin, stepf

C  Internal variables .............................................
      integer
     .  ie, ik, io, ispin, iuo, j, jo, juo
      double precision
     .  ckxij, ee, kxij, pipj1, pipj2, qe, skxij, t
C  ....................

C Solve eigenvalue problem .........................................
      do ik = 1,nk
        do ispin = 1,nspin
          do iuo = 1,nuo
            do juo = 1,nuo
              Saux(1,juo,iuo) = 0.d0
              Saux(2,juo,iuo) = 0.d0
              Haux(1,juo,iuo) = 0.d0
              Haux(2,juo,iuo) = 0.d0
            enddo
          enddo
          do io = 1,nuo
            do j = 1,numh(io)
              jo = listh(j,io)
              iuo = indxuo(io)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xij(1,j,io) +
     .               kpoint(2,ik) * xij(2,j,io) +
     .               kpoint(3,ik) * xij(3,j,io)
              ckxij = cos(kxij)
              skxij = sin(kxij)
              Saux(1,iuo,juo)=Saux(1,iuo,juo)+S(j,io)*ckxij
              Saux(2,iuo,juo)=Saux(2,iuo,juo)+S(j,io)*skxij
              Haux(1,iuo,juo)=Haux(1,iuo,juo)+H(j,io,ispin)*ckxij
              Haux(2,iuo,juo)=Haux(2,iuo,juo)+H(j,io,ispin)*skxij
            enddo
          enddo
          do iuo = 1,nuo
            do juo = 1,iuo-1
              Saux(1,juo,iuo) = 0.5d0 * ( Saux(1,juo,iuo) +
     .                                    Saux(1,iuo,juo) )
              Saux(1,iuo,juo) = Saux(1,juo,iuo)
              Saux(2,juo,iuo) = 0.5d0 * ( Saux(2,juo,iuo) -
     .                                    Saux(2,iuo,juo) )
              Saux(2,iuo,juo) = -Saux(2,juo,iuo)
              Haux(1,juo,iuo) = 0.5d0 * ( Haux(1,juo,iuo) +
     .                                    Haux(1,iuo,juo) )
              Haux(1,iuo,juo) = Haux(1,juo,iuo)
              Haux(2,juo,iuo) = 0.5d0 * ( Haux(2,juo,iuo) -
     .                                    Haux(2,iuo,juo) )
              Haux(2,iuo,juo) = -Haux(2,juo,iuo)
            enddo
            Saux(2,iuo,iuo) = 0.d0
            Haux(2,iuo,iuo) = 0.d0
          enddo
          call cdiag( Haux, nuo, Saux, nuo, nuo,
     .                eo(1,ispin,ik), psi, nuo, aux )
        enddo
      enddo
C ....................

C Check if we are done ................................................
      if (.not.getD) return
C ....................

C Find new Fermi energy and occupation weights ........................
      if (fixspin) then
        call fermispin( nspin, maxspn, nk, wk, maxuo, nuo, eo,
     .               temp, qs, qo, efs )
      else
        call fermid( nspin, maxspn, nk, wk, maxuo, nuo, eo, 
     .               temp, qtot, qo, ef )
      endif
C ....................

C Find weights for local density of states ............................
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.d-6 )
        do ik = 1,nk
          do ispin = 1,nspin
            do io = 1,nuo
              qo(io,ispin,ik) = wk(ik) * 
     .             ( stepf( (eo(io,ispin,ik)-e2)/t ) -
     .               stepf( (eo(io,ispin,ik)-e1)/t ) ) / nspin
            enddo
          enddo
        enddo
      endif
C ....................
      
c New density and energy-density matrices of unit-cell orbitals .......

      do ispin = 1,nspin
        do io = 1,nuo
          do j = 1,numd(io)
            Dnew(j,io,ispin) = 0.d0
            Enew(j,io,ispin) = 0.d0
          enddo
        enddo
      enddo

      do ik = 1,nk
        do ispin = 1,nspin

C         Find eigenvectors again (were stored only for one k point)
          do iuo = 1,nuo
            do juo = 1,nuo
              Saux(1,juo,iuo) = 0.d0
              Saux(2,juo,iuo) = 0.d0
              Haux(1,juo,iuo) = 0.d0
              Haux(2,juo,iuo) = 0.d0
            enddo
          enddo
          do io = 1,nuo
            do j = 1,numh(io)
              jo = listh(j,io)
              iuo = indxuo(io)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xij(1,j,io) +
     .               kpoint(2,ik) * xij(2,j,io) +
     .               kpoint(3,ik) * xij(3,j,io)
              ckxij = cos(kxij)
              skxij = sin(kxij)
              Saux(1,iuo,juo)=Saux(1,iuo,juo)+S(j,io)*ckxij
              Saux(2,iuo,juo)=Saux(2,iuo,juo)+S(j,io)*skxij
              Haux(1,iuo,juo)=Haux(1,iuo,juo)+H(j,io,ispin)*ckxij
              Haux(2,iuo,juo)=Haux(2,iuo,juo)+H(j,io,ispin)*skxij
            enddo
          enddo
          do iuo = 1,nuo
            do juo = 1,iuo-1
              Saux(1,juo,iuo) = 0.5d0 * ( Saux(1,juo,iuo) +
     .                                    Saux(1,iuo,juo) )
              Saux(1,iuo,juo) = Saux(1,juo,iuo)
              Saux(2,juo,iuo) = 0.5d0 * ( Saux(2,juo,iuo) -
     .                                    Saux(2,iuo,juo) )
              Saux(2,iuo,juo) = -Saux(2,juo,iuo)
              Haux(1,juo,iuo) = 0.5d0 * ( Haux(1,juo,iuo) +
     .                                    Haux(1,iuo,juo) )
              Haux(1,iuo,juo) = Haux(1,juo,iuo)
              Haux(2,juo,iuo) = 0.5d0 * ( Haux(2,juo,iuo) -
     .                                    Haux(2,iuo,juo) )
              Haux(2,iuo,juo) = -Haux(2,juo,iuo)
            enddo
            Saux(2,iuo,iuo) = 0.d0
            Haux(2,iuo,iuo) = 0.d0
          enddo
          call cdiag( Haux, nuo, Saux, nuo, nuo,
     .                eo(1,ispin,ik), psi, nuo, aux )

C         Add contribution to density matrices of unit-cell orbitals
C         WARNING: Dk and Ek may be EQUIVALENCE'd to Haux and Saux
          do iuo = 1,nuo
            do juo = 1,nuo
              Dk(1,juo,iuo) = 0.d0
              DK(2,juo,iuo) = 0.d0
              Ek(1,juo,iuo) = 0.d0
              Ek(2,juo,iuo) = 0.d0
            enddo
          enddo
          do ie = 1,nuo
            qe = qo(ie,ispin,ik)
            ee = qo(ie,ispin,ik) * eo(ie,ispin,ik)
            do iuo = 1,nuo
              do juo = 1,nuo
                pipj1 = psi(1,iuo,ie) * psi(1,juo,ie) +
     .                  psi(2,iuo,ie) * psi(2,juo,ie)
                pipj2 = psi(1,iuo,ie) * psi(2,juo,ie) -
     .                  psi(2,iuo,ie) * psi(1,juo,ie)
                Dk(1,juo,iuo) = Dk(1,juo,iuo) + qe * pipj1
                Dk(2,juo,iuo) = Dk(2,juo,iuo) + qe * pipj2
                Ek(1,juo,iuo) = Ek(1,juo,iuo) + ee * pipj1
                Ek(2,juo,iuo) = Ek(2,juo,iuo) + ee * pipj2
              enddo
            enddo
          enddo
          do io = 1,nuo
            do j = 1,numd(io)
              jo = listd(j,io)
              iuo = indxuo(io)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xij(1,j,io) +
     .               kpoint(2,ik) * xij(2,j,io) +
     .               kpoint(3,ik) * xij(3,j,io)
              ckxij = cos(kxij)
              skxij = sin(kxij)
              Dnew(j,io,ispin)=Dnew(j,io,ispin)+ Dk(1,juo,iuo)*ckxij -
     .                                           Dk(2,juo,iuo)*skxij
              Enew(j,io,ispin)=Enew(j,io,ispin)+ Ek(1,juo,iuo)*ckxij -
     .                                           Ek(2,juo,iuo)*skxij
            enddo
          enddo
          
        enddo
      enddo
C ....................

      end


