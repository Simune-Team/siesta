C $Id: on_subs.f,v 1.2 1999/01/31 11:20:10 emilio Exp $

      subroutine ctrans(nr,nc,nmax,ntmax,
     .              num,list,numt,listt,cttoc)
C ********************************************************************
C Finds the C transpose matrix control vectors numt and listt,
C and the index vector cttoc.
C Written by P.Ordejon. October'96
C ***************************** INPUT *********************************
C integer nr             : Number of rows of C 
C integer nc             : Number of columns of full C
C integer nmax           : First dimension of list and C, and maximum
C                           number of nonzero elements of each row of C
C integer ntmax          : maximum number of nonzero elements of each
C                           column of C
C integer num(nr)        : Control vector of C matrix
C                           (number of nonzero elements of each row of C)
C integer list(nmax,nr)  : Control vector of C matrix
C                          (list of nonzero elements of each row of C)
C **************************** OUTPUT *********************************
C integer numt(nc)       : Control vector of C matrix
C                         (number of nonzero elements of each column of C)
C integer listt(ntmax,nc) : Control vector of C transpose matrix
C                           (list of nonzero elements of each column of C)
C integer cttoc(ntmax,nc) : Map from C transpose to C indexing
C *********************************************************************
      implicit none

      integer
     .  nc,nmax,nr,ntmax,
     .  cttoc(ntmax,nc),list(nmax,nr),listt(ntmax,nc),
     .  num(nr),numt(nc)

C Internal variables ..................................................
      integer
     .  i,imu,mu,n
C ..........................
C  Initialize numt list ...............................................
      do i = 1,nc
        numt(i) = 0
      enddo
C ..........................
C  Construct information for transpose of C ............................
      do mu = 1,nr
        do imu = 1,num(mu)
          i = list(imu,mu)
          numt(i) = numt(i)+1
          n = numt(i)
          if (n .gt. ntmax) goto 10
          listt(n,i) = mu
          cttoc(n,i) = imu
10        continue
        enddo
      enddo
C ..........................
C Check dimensions .....................................................
      do i = 1,nc
        call chkdim('ctrans','ntmax',ntmax,numt(i),1)
      enddo
C ..........................
      return
      end
      subroutine axb_build(nra,nca,namax,numa,lista,
     .               nrb,ncb,nbmax,numb,listb,
     .               ind,nindv,
     .               ncmax,numc,listc)
C ********************************************************************
C Constructs control indexes of a C matrix in sparse form,
C C being the product of A and B (also in sparse form)
C
C              C = A x B
C
C In full form: A is rectangular, and has dimension:  nra x nca
C               B is rectangular, and has dimension:  nrb x ncb
C and, as a result:
C               C is rectangular, and has dimension:  nra x ncb
C (Of course, nca must be equal to nrb)
C
C Written by P.Ordejon. October'96
C ***************************** INPUT *********************************
C integer nra               : Number of rows of A 
C integer nca               : Number of columns of A
C integer namax             : First dimension of A matrix in sparse form,
C                             as declared in calling routine
C                             (max. number of <>0 elements of each row of A)
C integer numa(nra)         : Control vector of A matrix
C                            (number of nonzero elements of each row of A)
C integer lista(namax,nra)  : Control vector of A matrix
C                           (list of nonzero elements of each row of A)
C integer nrb               : Number of rows of B
C integer ncb               : Number of columns of B
C integer nbmax             : First dimension of B matrix in sparse form,
C                             as declared in calling routine
C                             (max. number of <>0 elements of each row of B)
C integer numb(nrb)          : Control vector of B matrix
C                            (number of nonzero elements of each row of B)
C integer listb(nbmax,nrb)   : Control vector of B matrix
C                            (list of nonzero elements of each row of B)
C integer ind(ncb)          : Auxiliary array to build C in sparse form
C integer nindv(ncmax)      : Auxiliary array to store indexes of nonzero
C                             matrix elements of each row of C
C integer ncmax             : First dimension of C matrix in sparse form,
C                             as declared in calling routine
C                             (max. number of <>0 elements of each row of C)
C **************************** OUTPUT *********************************
C integer numc(nra)          : Control vector of C matrix
C                            (number of nonzero elements of each row of C)
C integer listc(ncmax,nra)   : Control vector of C matrix
C                            (list of nonzero elements of each row of C)
C *********************************************************************
      implicit none

      integer
     .  nca,ncb,namax,nbmax,ncmax,nra,nrb,
     .  ind(ncb),
     .  lista(namax,nra),listb(nbmax,nrb),listc(ncmax,nra),
     .  nindv(ncmax),
     .  numa(nra),numb(nrb),numc(nra)

C Internal variables..................................................
      integer
     .  i,in,j,k,kn,nind
C............................

C Check dimensions ....................................................
      call chkdim('axb_build','nca',nca,nrb,0)
C ...........................
C Initialize internal variables .......................................
      nind=0
      do i = 1,ncb
        ind(i) = 0
      enddo
      do i = 1,ncmax
        nindv(i)=0
      enddo
C ...........................
C Find out C control vectors ..........................................
      do i = 1,nra
        do in = 1,numa(i)
          k = lista(in,i)
          do kn = 1,numb(k)
            j = listb(kn,k)
            if (ind(j) .eq. 0) then
              ind(j) = 1
              nind = nind+1
              nindv(nind) = j
            endif
          enddo
        enddo
        numc(i) = nind
        call chkdim ('axb_build','ncmax',ncmax,nind,1)
        do in = 1,nind
          j = nindv(in)
          nindv(in) = 0
          ind(j) = 0
          listc(in,i) = j
        enddo
        nind = 0
      enddo
C ...........................
      return
      end
      subroutine ind_gf(nr,nc,ncmax,nfmax,
     .              numc,listc,numf,listf,indgf)
C ********************************************************************
C Maps the F matrix into C matrix
C Written by P.Ordejon. October'96
C ***************************** INPUT *********************************
C integer nr             : Number of rows of C (columns of full F)
C integer nc             : Number of columns of full C (rows of F)
C integer ncmax          : First dimension of listc and C, and maximum
C                           number of nonzero elements of each row of C
C integer nfmax          : First dimension of listf and F, and maximum
C                           number of nonzero elements of each row of F
C integer numc(nr)       : Control vector of C matrix
C                           (number of nonzero elements of each row of C)
C integer listc(ncmax,nr): Control vector of C matrix
C                          (list of nonzero elements of each row of C)
C integer numf(nc)       : Control vector of F matrix
C                           (number of nonzero elements of each row of F)
C integer listf(nfmax,nc): Control vector of F matrix
C                          (list of nonzero elements of each row of F)
C **************************** OUTPUT *********************************
C integer indgf(ncmax,nr) : Map from F to C 
C                    indgf(i,j) is the column index (in sparse notation)
C                    of the element of F corresponding to the element
C                    (i,j) of C (in sparse notation)
C *********************************************************************
      implicit none

      integer
     .  nc,ncmax,nfmax,nr,
     .  indgf(ncmax,nr),listc(ncmax,nr),listf(nfmax,nc),
     .  numc(nr),numf(nc)

      integer
     .  i,imu,jk,mu,mmu

      do mu=1,nr
        do imu=1,numc(mu)
          i=listc(imu,mu)
          do jk=1,numf(i)
            mmu=listf(jk,i)
            if (mmu .eq. mu) indgf(imu,mu)=jk
          enddo
        enddo
      enddo

      return
      end
