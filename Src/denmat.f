C $Id: denmat.f,v 1.12 1999/02/23 12:05:22 wdpgaara Exp $

      subroutine denmat(c,eta,h,s,enum,
     .                 nbasis,nbands,ncmax,nctmax,
     .                 nfmax,nftmax,nhmax,nhijmax,
     .                 numc,listc,numct,listct,cttoc,
     .                 numf,listf,numft,listft,fttof,
     .                 numh,listh,numhij,listhij,chcc,cscc,
     .                 dm,edm)
C *******************************************************************
C Subroutine to compute the Density and Energy Density matrices
C for the Order-N functional of Kim et al. (PRB 52, 1640 (95))
C (generalization of that proposed by Mauri et al, and Ordejon et al)
C
C Density Matrix:
C  D_mu,nu = 2 * C_i,mu * ( 2 * delta_i,j - S_i,j) * C_j,nu
C
C Energy Density Matrix:
C  E_mu,nu = 2 * C_i,mu * ( H_i,j + 2 * eta * (delta_i,j - S_i,j) ) * C_j,nu
C
C (The factor 2 is for spin)
C
C (See Ordejon et al, PRB 51, 1456 (95))
C
C The DM is normalized to the exact number of electrons!!!
C
C Written by P.Ordejon, Noviembre'96
C Modified by J.M.Soler, May'97
C ************************** INPUT **********************************
C real*8 c(ncmax,nbasis)      : Localized Wave Functions (sparse)
C real*8 eta                  : Fermi level parameter of Kim et al.
C real*8 h(nhmax,nbasis)      : Hamiltonian matrix (sparse)
C real*8 s(nhmax,nbasis)      : Overlap matrix (sparse)
C real*8 enum                 : Total number of electrons
C integer nbasis              : Number of atomic orbitals
C integer nbands              : Number of Localized Wave Functions
C integer ncmax               : First dimension of listc and C, and maximum
C                               number of nonzero elements of each row of C
C integer nctmax              : Max num of <>0 elements of each col of C
C integer nfmax               : Max num of <>0 elements of each row of 
C                               F = Ct x H
C integer nftmax              : Max num of <>0 elements of each col of F
C integer nhmax               : First dimension of listh and H, and maximum
C                               number of nonzero elements of each row of H
C integer nhijmax             : Maximum number of non-zero elements of each
C                               row of Hij
C integer numc(nbasis)        : Control vector of C matrix
C                               (number of nonzero elements of each row of C)
C integer listc(ncmax,nbasis) : Control vector of C matrix
C                              (list of nonzero elements of each row of C)
C integer numct(nbands)       : Control vector of C transpose matrix
C                              (number of <>0  elements of each col of C)
C integer listct(ncmax,nbands): Control vector of C transpose matrix
C                              (list of <>0  elements of each col of C)
C integer cttoc(ncmax,nbands) : Map from Ct to C indexing
C integer numf(nbands)        : Control vector of F matrix
C                               (number of <>0  elements of each row of F)
C integer listf(nfmax,nbands) : Control vector of F matrix
C                               (list of <>0  elements of each row of F)
C integer numft(nbasis)        : Control vector of F transpose matrix
C                               (number of <>0  elements of each col of F)
C integer listft(nfmax,nbasis) : Control vector of F transpose matrix
C                               (list of <>0  elements of each col of F)
C integer fttof(nfmax,nbasis)  : Map from Ft to F indexing
C integer numh(nbasis)        : Control vector of H matrix
C                               (number of nonzero elements of each row of H)
C integer listh(nhmax,nbasis) : Control vector of H matrix
C                               (list of nonzero elements of each row of H)
C integer numhij(nbasis)       : Control vector of Hij matrix
C                                (number of <>0  elements of each row of Hij)
C integer listhij(nhijmax,nbands): Control vector of Hij matrix 
C                                (list of <>0  elements of each row of Hij)
C real*8 chcc(nfmax,nbands)   : Auxiliary space
C real*8 cscc(nfmax,nbands)   : Auxiliary space
C ************************* OUTPUT **********************************
C real*8 dm(nhmax,nbasis)     : Density Matrix
C real*8 edm(nhmax,nbasis)    : Energy density matrix
C *******************************************************************
      implicit none

      integer
     .  nbasis,nbands,ncmax,nctmax,nfmax,nhmax,nhijmax

      integer
     .  cttoc(nctmax,nbands),listc(ncmax,nbasis),
     .  listct(nctmax,nbands),listf(nfmax,nbands),
     .  listh(nhmax,nbasis),listhij(nhijmax,nbands),
     .  numc(nbasis),numct(nbands),numf(nbands),
     .  numh(nbasis),numhij(nbands)
     
      integer
     .  nftmax,fttof(nftmax,nbasis),listft(nftmax,nbasis),numft(nbasis)

      double precision
     .  c(ncmax,nbasis),dm(nhmax,nbasis),edm(nhmax,nbasis),enum,eta,
     .  h(nhmax,nbasis),s(nhmax,nbasis),
     .  chcc(nfmax,nbands),cscc(nfmax,nbands)
     
      external
     .  chkdim, iodm, timer

C Internal variales ..................................................
C   Notation hints:
C     m,n : basis orbital inexes (mu,nu)
C     i,j : band (and LWF) indexes
C     im  : index for LWF's of basis orbital m
C     mi  : index for basis orbitals of LWF i
C     nm  : index for basis orbitals connected to basis orbital m

      include 'ordern.h'

      integer 
     .  i, in, im, j, ji, jm, jn, m, mi, mn, n, ni, nm

      double precision
     .  cHrow(maxo), cSrow(maxo),
     .  chcrow(maxlwf), cscrow(maxlwf),
     .  chccCol(maxlwf), csccCol(maxlwf),
     .  cim, cnj, chin, csin, chccim, csccim, cchccmn, ccsccmn,
     .  Hmn, Smn, qout, fact

C ........................

C Start time counter .....................................................
      call timer('denmat',1)
C .......................

C Check dimensions .........................................................
      call chkdim( 'denmat', 'maxo',   maxo,   nbasis, 1 )
      call chkdim( 'denmat', 'maxlwf', maxlwf, nbands, 1 )
C .................

C Initialize temporary arrays ..............................................
      do m = 1,nbasis
        cHrow(m) = 0.d0
        cSrow(m) = 0.d0
      enddo
      do i = 1,nbands
        cscrow(i) = 0.d0
        chcrow(i) = 0.d0
        csccCol(i) = 0.d0
        chccCol(i) = 0.d0
      enddo
C ........................

C Find cscc=(2-ct*S*c)*ct and chcc=(ct*H*c+2eta(1-ct*S*c))*ct.
      do i = 1,nbands
      
C       Find row i of cS=ct*S and cH=ct*H
        do mi = 1,numct(i)
          m = listct(mi,i)
          im = cttoc(mi,i)
          cim = c(im,m)
          do nm = 1,numh(m)
            n = listh(nm,m)
            Smn = S(nm,m)
            Hmn = H(nm,m)
            cSrow(n) = cSrow(n) + cim * Smn
            cHrow(n) = cHrow(n) + cim * Hmn
          enddo
        enddo

C       Find row i of csc=2-ct*S*c and chc=ct*H*c+2eta(1-ct*S*c)
C       First set diagonal terms 2 and 2eta
        cscrow(i) = 2.d0
        chcrow(i) = 2.d0 * eta
C       Now use the list of nonzero elements of f=ct*H
        do ni = 1,numf(i)
          n = listf(ni,i)
          csin = - cSrow(n)
          chin = cHrow(n) - 2.d0*eta*cSrow(n)
          do jn = 1,numc(n)
            j = listc(jn,n)
            cnj = c(jn,n)
            cscrow(j) = cscrow(j) + csin * cnj
            chcrow(j) = chcrow(j) + chin * cnj
          enddo
C         Restore cSrow and cHrow for next row
          cSrow(n) = 0.d0
          cHrow(n) = 0.d0
        enddo
        
C       Find row i of cscc=csc*ct and chcc=chc*ct. 
C       Only the nonzero elements of f=cH will be required.
        do mi = 1,numf(i)
          m = listf(mi,i)
          csccim = 0.d0
          chccim = 0.d0
          do jm = 1,numc(m)
            j = listc(jm,m)
            csccim = csccim + cscrow(j) * c(jm,m)
            chccim = chccim + chcrow(j) * c(jm,m)
          enddo
          cscc(mi,i) = csccim
          chcc(mi,i) = chccim
        enddo
        
C       Restore cscrow and chcrow for next row
        do ji = 1,numhij(i)
          j = listhij(ji,i)
          cscrow(j) = 0.d0
          chcrow(j) = 0.d0
        enddo
      enddo
C ........................

C Find dm=c*cscc and edm=c*chcc. Only the nonzero elements of H.
      do n = 1,nbasis
      
C       Use listft to expand a column of cscc
        do in = 1,numft(n)
          i = listft(in,n)
          ni = fttof(in,n)
          csccCol(i) = cscc(ni,i)
          chccCol(i) = chcc(ni,i)
        enddo
        
C       Find column n of c*cscc and c*chcc
C       Use that H is symmetric to determine required elements
        do mn = 1,numh(n)
          m = listh(mn,n)
C         Find element (m,n) of c*cscc and c*chcc
          ccsccmn = 0.d0
          cchccmn = 0.d0
          do im = 1,numc(m)
            i = listc(im,m)
            ccsccmn = ccsccmn + c(im,m) * csccCol(i)
            cchccmn = cchccmn + c(im,m) * chccCol(i)
          enddo
C         Use that dm and edm are symmetric
          dm(mn,n)  = 2.d0 * ccsccmn
          edm(mn,n) = 2.d0 * cchccmn
        enddo
        
C       Restore csccCol and chccCol for next column
        do in = 1,numft(n)
          i = listft(in,n)
          csccCol(i) = 0.d0
          chccCol(i) = 0.d0
        enddo
      enddo
C ........................

C Normalize DM to exact charge .........................
C Calculate total output charge ...
      qout = 0.0d0
      do i = 1,nbasis
        do in = 1,numh(i)
          qout = qout + dm(in,i) * s(in,i)
        enddo
      enddo
      write(6,"(/a,f12.4)") 'denmat: qtot (before DM normalization) = ',
     .              qout
C ...

      if (dabs(enum-qout) .gt. 0.05d0) then
        fact = enum / qout
      
C Normalize ...
        do i = 1,nbasis
          do in = 1,numh(i)
            dm(in,i) = dm(in,i) * fact
            edm(in,i) = edm(in,i) * fact
          enddo
        enddo
C ...
      endif
C ........................

C Save DM to disk: now it is done at top level ............................
C     call iodm('write',nhmax,nbasis,nbasis,1,numh,listh,dm,found )
C ......................

C Stop time counter and return ..................
      call timer('denmat',2)
      end

