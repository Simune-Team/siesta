      subroutine ener3(c,grad,lam,eta,enum,h,s,
     .                 nbasis,nbands,ncmax,nctmax,
     .                 nfmax,nhmax,nhijmax,
     .                 numc,listc,numct,listct,cttoc,numf,listf,
     .                 numh,listh,numhij,listhij,
     .                 ener)

C ************************************************************************
C Finds the energy at three points of the line passing thru C in the
C direction of GRAD. LAM is the distance (in units of GRAD) between 
C points.
C Uses the functional of Kim et al (PRB 52, 1640 (95))
C Works only with spin-unpolarized systems
C Written by P.Ordejon. October'96
C ****************************** INPUT ***********************************
C real*8 c(ncmax,nbasis)       : Current point (wave function coeff.
C                                  in sparse form)
C real*8 grad(ncmax,nbasis)    : Direction of search (sparse)
C real*8 lam                   : Length of step
C real*8 eta                   : Fermi level parameter of Kim et al.
C real*8 enum                  : Total number of electrons
C real*8 h(nhmax,nbasis)       : Hamiltonian matrix (sparse)
C real*8 s(nhmax,nbasis)       : Overlap matrix (sparse)
C integer nbasis               : Number of basis orbitals
C integer nbands               : Number of LWF's
C integer ncmax                : Max num of <>0 elements of each row of C
C integer nctmax               : Max num of <>0 elements of each col of C
C integer nfmax                : Max num of <>0 elements of each row of 
C                                   F = Ct x H
C integer nhmax                : Max num of <>0 elements of each row of H
C integer nhijmax              : Max num of <>0 elements of each row of 
C                                   Hij=Ct x H x C
C integer numc(nbasis)         : Control vector of C matrix
C                                (number of <>0  elements of each row of C)
C integer listc(ncmax,nbasis)  : Control vector of C matrix 
C                               (list of <>0  elements of each row of C)
C integer numct(nbands)        : Control vector of C transpose matrix
C                               (number of <>0  elements of each col of C)
C integer listct(ncmax,nbands) : Control vector of C transpose matrix
C                               (list of <>0  elements of each col of C)
C integer cttoc(ncmax,nbands)  : Map from Ct to C indexing
C integer numf(nbands)         : Control vector of F matrix
C                                (number of <>0  elements of each row of F)
C integer listf(nfmax,nbands)  : Control vector of F matrix 
C                                (list of <>0  elements of each row of F)
C integer numh(nbasis)         : Control vector of H matrix
C                                (number of <>0  elements of each row of H)
C integer listh(nhmax,nbasis)  : Control vector of H matrix 
C                               (list of <>0  elements of each row of H)
C integer numhij(nbasis)       : Control vector of Hij matrix
C                                (number of <>0  elements of each row of Hij)
C integer listhij(nhijmax,nbands): Control vector of Hij matrix 
C                                (list of <>0  elements of each row of Hij)
C ***************************** OUTPUT ***********************************
C real*8 ener(3)               : Energy at the three points:
C                                     C +     lam * GRAD
C                                     C + 2 * lam * GRAD
C                                     C + 3 * lam * GRAD
C ************************************************************************

      implicit none

      integer
     .  nbasis,nbands,ncmax,nctmax,nfmax,nhmax,nhijmax

      integer
     .  cttoc(nctmax,nbands),listc(ncmax,nbasis),
     .  listct(nctmax,nbands),listf(nfmax,nbands),
     .  listh(nhmax,nbasis),listhij(nhijmax,nbands),
     .  numc(nbasis),numct(nbands),numf(nbands),
     .  numh(nbasis),numhij(nbands)

      double precision
     .  c(ncmax,nbasis),ener(3),eta,enum,grad(ncmax,nbasis),
     .  h(nhmax,nbasis),lam,s(nhmax,nbasis)

C Internal variables ......................................................
C maxo   = maximum number of basis orbitals
C maxlwf = maximum number of LWF's
C maxnct = maximum number of non-zero elements of each column of C
      integer 
     .  maxlwf, maxo,  maxnct

      parameter ( maxo   = 10000 )
      parameter ( maxlwf =  5000 )
      parameter ( maxnct =  5000 )
     

      integer
     .  i,in,jn,k,kn

      double precision
     .  a1,a2,a3,b1,b2,b3,c1,c2,c3,func1(3),func2(3),
     .  lam1,lam2,lam3,
     .  p1(maxnct),p2(maxnct),p3(maxnct),pp1,pp2,pp3,
     .  aux1(maxo),aux2(maxo),aux3(maxo),
     .  aux4(maxo),aux5(maxo),aux6(maxo),
     .  bux1(maxlwf),bux2(maxlwf),bux3(maxlwf),
     .  bux4(maxlwf),bux5(maxlwf),bux6(maxlwf)
C..................

      call timer('ener3',1)

C Check dimensions .........................................................
      call chkdim('ener3','maxo',maxo,nbasis,1)
      call chkdim('ener3','maxlwf',maxlwf,nbands,1)
      call chkdim('ener3','maxnct',maxnct,nctmax,1)
C .................


C Define points to compute energy ..........................................
      lam1 = lam
      lam2 = lam*2.
      lam3 = lam*3.
C..................

C Initialize output and auxiliary varialbles ...............................
      do i = 1,3
        ener(i) = 0.0d0
        func1(i) = 0.0d0
        func2(i) = 0.0d0
      enddo

      do i = 1,nbasis
        aux1(i) = 0.0d0
        aux2(i) = 0.0d0
        aux3(i) = 0.0d0
        aux4(i) = 0.0d0
        aux5(i) = 0.0d0
        aux6(i) = 0.0d0
      enddo

      do i = 1,nbands
        bux1(i) = 0.0d0
        bux2(i) = 0.0d0
        bux3(i) = 0.0d0
        bux4(i) = 0.0d0
        bux5(i) = 0.0d0
        bux6(i) = 0.0d0
      enddo
C..................


C Calculate Functional .....................................................
C F=CtH
C Fs=CtS

      do i = 1,nbands
        do in = 1,numct(i)
          p1(in) = c(cttoc(in,i),listct(in,i)) +
     .             lam1 * grad(cttoc(in,i),listct(in,i))
          p2(in) = c(cttoc(in,i),listct(in,i)) +
     .             lam2 * grad(cttoc(in,i),listct(in,i))
          p3(in) = c(cttoc(in,i),listct(in,i)) +
     .             lam3 * grad(cttoc(in,i),listct(in,i))
        enddo

        do in = 1,numct(i)
          k = listct(in,i)
          pp1 = p1(in)
          pp2 = p2(in)
          pp3 = p3(in)

          do kn = 1,numh(k)
            aux1(listh(kn,k)) = aux1(listh(kn,k)) + pp1 * h(kn,k)
            aux2(listh(kn,k)) = aux2(listh(kn,k)) + pp2 * h(kn,k)
            aux3(listh(kn,k)) = aux3(listh(kn,k)) + pp3 * h(kn,k)
            aux4(listh(kn,k)) = aux4(listh(kn,k)) + pp1 * s(kn,k)
            aux5(listh(kn,k)) = aux5(listh(kn,k)) + pp2 * s(kn,k)
            aux6(listh(kn,k)) = aux6(listh(kn,k)) + pp3 * s(kn,k)
          enddo
        enddo
        do in = 1,numf(i)
          k = listf(in,i)
          a1 = aux1(k)
          a2 = aux2(k)
          a3 = aux3(k)
          b1 = aux4(k)
          b2 = aux5(k)
          b3 = aux6(k)
          aux1(k) = 0.0
          aux2(k) = 0.0
          aux3(k) = 0.0
          aux4(k) = 0.0
          aux5(k) = 0.0
          aux6(k) = 0.0

C Hij=CtHC
C Sij=CtSC
C multiply FxC and FsxC
          do kn = 1,numc(k)
            c1 = c(kn,k) + lam1 * grad(kn,k)
            c2 = c(kn,k) + lam2 * grad(kn,k)
            c3 = c(kn,k) + lam3 * grad(kn,k)
            bux1(listc(kn,k)) = bux1(listc(kn,k)) + a1 * c1
            bux2(listc(kn,k)) = bux2(listc(kn,k)) + b1 * c1
            bux3(listc(kn,k)) = bux3(listc(kn,k)) + a2 * c2
            bux4(listc(kn,k)) = bux4(listc(kn,k)) + b2 * c2
            bux5(listc(kn,k)) = bux5(listc(kn,k)) + a3 * c3
            bux6(listc(kn,k)) = bux6(listc(kn,k)) + b3 * c3
          enddo
        enddo
        func1(1) = func1(1) + bux1(i) - eta * bux2(i)
        func1(2) = func1(2) + bux3(i) - eta * bux4(i)
        func1(3) = func1(3) + bux5(i) - eta * bux6(i)
        do jn = 1,numhij(i)
          func2(1) = func2(1) + (bux1(listhij(jn,i)) -
     .             eta * bux2(listhij(jn,i))) * bux2(listhij(jn,i))
          func2(2) = func2(2) + (bux3(listhij(jn,i)) -
     .             eta * bux4(listhij(jn,i))) * bux4(listhij(jn,i))
          func2(3) = func2(3) + (bux5(listhij(jn,i)) -
     .             eta * bux6(listhij(jn,i))) * bux6(listhij(jn,i))
          bux1(listhij(jn,i)) = 0.0
          bux2(listhij(jn,i)) = 0.0
          bux3(listhij(jn,i)) = 0.0
          bux4(listhij(jn,i)) = 0.0
          bux5(listhij(jn,i)) = 0.0
          bux6(listhij(jn,i)) = 0.0
        enddo
      enddo

C This is valid for an spin-unpolarized sytem
      do i=1,3
        ener(i) = 2.0d0 * func1(i) - func2(i) + eta * enum / 2.0d0
      enddo
C ...................

      call timer('ener3',2)
      return
      end

