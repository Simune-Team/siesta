      subroutine gradient(c,eta,enum,h,s,
     .                 nbasis,nbands,ncmax,nctmax,
     .                 nfmax,nftmax,nhmax,nhijmax,
     .                 numc,listc,numct,listct,cttoc,numf,listf,
     .                 numft,listft,fttof,
     .                 numh,listh,numhij,listhij,indgf,f,fs,
     .                 grad,ener)

C ************************************************************************
C Finds the energy and gradient at point C.
C Uses the functional of Kim et al (PRB 52, 1640 (95))
C Works only with spin-unpolarized systems
C Written by P.Ordejon. October'96
C Last modified: J.M.Soler. 30/04/97
C ****************************** INPUT ***********************************
C real*8 c(ncmax,nbasis)       : Current point (wave function coeff.
C                                  in sparse form)
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
C integer nftmax               : Max num of <>0 elements of each col of F
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
C integer numft(nbasis)        : Control vector of F transpose matrix
C                               (number of <>0  elements of each col of F)
C integer listft(nfmax,nbasis) : Control vector of F transpose matrix
C                               (list of <>0  elements of each col of F)
C integer fttof(nfmax,nbasis)  : Map from Ft to F indexing
C integer numh(nbasis)         : Control vector of H matrix
C                                (number of <>0  elements of each row of H)
C integer listh(nhmax,nbasis)  : Control vector of H matrix 
C                               (list of <>0  elements of each row of H)
C integer numhij(nbasis)       : Control vector of Hij matrix
C                                (number of <>0  elements of each row of Hij)
C integer listhij(nhijmax,nbands): Control vector of Hij matrix 
C                                (list of <>0  elements of each row of Hij)
C integer indgf(ncmax,nbasis)  : Maps elements of GRAD and F
C real*8 f(nfmax,nbands)       : Auxiliary space
C real*8 fs(nfmax,nbands)      : Auxiliary space
C ***************************** OUTPUT ***********************************
C real*8 ener                  : Energy at point C
C real*8 grad(ncmax,nbasis)    : Gradient of functional (sparse)
C ************************************************************************
      implicit none

      integer
     .  nbasis,nbands,ncmax,nctmax,nfmax,nftmax,nhmax,nhijmax

      integer
     .  cttoc(nctmax,nbands),fttof(nftmax,nbasis),
     .  indgf(ncmax,nbasis),
     .  listc(ncmax,nbasis),listct(nctmax,nbands),
     .  listf(nfmax,nbands),listft(nftmax,nbasis),
     .  listh(nhmax,nbasis),listhij(nhijmax,nbands),
     .  numc(nbasis),numct(nbands),numf(nbands),numft(nbasis),
     .  numh(nbasis),numhij(nbands)

      double precision
     .  c(ncmax,nbasis),ener,eta,enum,grad(ncmax,nbasis),
     .  h(nhmax,nbasis),s(nhmax,nbasis),
     .  f(nfmax,nbands),fs(nfmax,nbands)

C Internal variables ......................................................
C maxo   = maximum number of basis orbitals
C maxlwf = maximum number of LWF's
*     integer
*    .  maxlwf, maxo
*     parameter ( maxo   = 10000 )
*     parameter ( maxlwf =  5000 )
      include 'ordern.h'

      integer
     .  i,ik,in,j,jk,jn,k,kn,mu,muk

      double precision
     .  a0,b0,p0,func1,func2,
     .  aux1(maxo),aux2(maxo),
     .  bux1(maxlwf),bux2(maxlwf)
     
C-JMS New local arrays
      double precision ft(maxnf,maxo), fst(maxnf,maxo)
      logical          frstme
      data frstme /.true./     
C..................

      call timer('gradient',1)
      
C-JMS Print array sizes ....................................................
      if (frstme) then
        call prmem( 0, 'gradient', 'ft',  'd', maxnf*maxo )
        call prmem( 0, 'gradient', 'fst', 'd', maxnf*maxo )
        call prmem( 0, 'gradient', ' ',   ' ', 0          )
        frstme = .false.
      endif
C ................

C Check dimensions .........................................................
*     call chkdim('ener3','maxo',maxo,nbasis,1)
*     call chkdim('ener3','maxlwf',maxlwf,nbands,1)
      call chkdim('gradient','maxo',maxo,nbasis,1)
      call chkdim('gradient','maxlwf',maxlwf,nbands,1)
C-JMS Added check
      call chkdim('gradient','maxnf',maxnf,nfmax,1)
C .................

C Initialize output and auxiliary varialbles ...............................

      ener = 0.0d0
      func1 = 0.0d0
      func2 = 0.0d0

      do i = 1,nbasis
        aux1(i) = 0.0d0
        aux2(i) = 0.0d0
      enddo

      do i = 1,nbands
        bux1(i) = 0.0d0
        bux2(i) = 0.0d0
      enddo
C..................


C Calculate Functional .....................................................

C F=CtH  ---> JMS: F=Ct*(H-eta*S)
C Fs=CtS
      do i = 1,nbands
        do in = 1,numct(i)
          k = listct(in,i)
          ik = cttoc(in,i)
          p0 = c(ik,k)
          do kn = 1,numh(k)
*           aux1(listh(kn,k)) = aux1(listh(kn,k)) + p0 * h(kn,k)
            aux1(listh(kn,k)) = aux1(listh(kn,k)) + 
     .                          p0 * ( h(kn,k) - eta*s(kn,k) )
            aux2(listh(kn,k)) = aux2(listh(kn,k)) + p0 * s(kn,k)
          enddo
        enddo
        do in = 1,numf(i)
          k = listf(in,i)
          f(in,i) = aux1(k)
          fs(in,i) = aux2(k)
          aux1(k) = 0.0
          aux2(k) = 0.0
        enddo
      enddo
      
C-JMS Find transpose of F and Fs
      do mu = 1,nbasis
        do muk = 1,numft(mu)
          j = listft(muk,mu)
          jk = fttof(muk,mu)
          ft(muk,mu) = f(jk,j)
          fst(muk,mu) = fs(jk,j)
        enddo
      enddo

C Hij=CtHC
C Sij=CtSC
C multiply FxC and FsxC row by row
      do i = 1,nbands
        do in = 1,numf(i)
          k = listf(in,i)
          a0 = f(in,i)
          b0 = fs(in,i)
          do kn = 1,numc(k)
            bux1(listc(kn,k)) = bux1(listc(kn,k)) + a0 * c(kn,k)
            bux2(listc(kn,k)) = bux2(listc(kn,k)) + b0 * c(kn,k)
          enddo
        enddo
*       func1 = func1 + bux1(i) - eta * bux2(i)
        func1 = func1 + bux1(i)

C multiply Hij x Fs and Sij x F row by row
C (only products of neccesary elements)
        do ik = 1,numct(i)
          mu = listct(ik,i)
          a0 = 0.d0
          do muk = 1,numft(mu)
            j = listft(muk,mu)
*           jk = fttof(muk,mu)
*           a0 = a0 +
*    .           bux1(j) * fs(jk,j) +
*    .           bux2(j) * (f(jk,j) - 2.d0 * eta * fs(jk,j))
            a0 = a0 + bux1(j) * fst(muk,mu) + bux2(j) * ft(muk,mu)
          enddo
          grad(cttoc(ik,i),mu) = (-2.d0) * a0
        enddo

        do jn = 1,numhij(i)
          j = listhij(jn,i)
*         func2 = func2 + (bux1(j) - eta * bux2(j)) * bux2(j)
          func2 = func2 + bux1(j) * bux2(j)
          bux1(j) = 0.0
          bux2(j) = 0.0
        enddo
      enddo

      do k = 1,nbasis
        do ik = 1,numc(k)
	  i = indgf(ik,k)
          j = listc(ik,k)
*         grad(ik,k) = 4.d0 * (f(i,j) - eta * fs(i,j)) + grad(ik,k)
          grad(ik,k) = 4.d0 * f(i,j) + grad(ik,k)
        enddo
      enddo

      ener = 2.d0 * func1 - func2 + eta * enum / 2.d0
C ...................

      call timer('gradient',2)
      return
      end

