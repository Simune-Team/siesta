C $Id: pulayx.f,v 1.9 1999/04/13 17:03:24 emilio Exp $

       subroutine pulayx(pulfile,iscf,mix,nbasis,maxo,maxno,numd,
     1                   nspin,maxsav,alpha,nkick,alphakick,
     2                   savedm,savere,dimaux,dmnew,dmold,dmax)
C ***********************************************************************
C Pulay mixing implemented to accelerate the self-consistency
C Mixes MAXSAV previous steps.  Linear mixing if MAXAV =< 0
C 
C Written by In-Ho Lee, Beckman Inst., Univ. of Illinois, Mar. 25 '97
C Modified and partly re-written by P. Ordejon, July'97
C Modified and optimized by P. Ordejon, November'97
C Algorithm modified by P. Ordejon, February'99
C ************************** INPUT **************************************
C logical pulfile            : Use file to store Pulay info
C                                 .true.  = use file
C                                 .false. = use memory
C integer iscf               : Current SCF iteration
C logical mix                : Mix first SCF step (T or F)
C integer nbasis             : Number of atomic orbitals 
C integer maxo               : Maximum number of atomic orbitals (dimensio)
C integer maxno              : First dimension of D.M., and 
C                              maximum number of nonzero elements of 
C                              each row of the D.M.
C integer numd(maxo)         : Control vector of D.M.
C                              (number of nonzero elements of each row)
C integer nspin              : Spin polarization (1=unpolarized, 2=polarized)
C integer maxsav             : Pulay mixing is done every maxsav iterations.
C                              If maxav=<0, linear mixing is done
C real*8 alpha               : Mixing parameter (for linear mixing)
C integer nkick              : Do a linear mixing (kick) each nkick cycles
C real*8 alphakick           : Mixing parameter for kick cycles
C real*4 savedm(dimaux)      : Auxiliary storage (DM in former iterations)
C real*4 savere(dimaux)      : Auxiliary storage (resuduals in former iterations)
C integer dimaux             : Auxiliart matrices size
C ********************* INPUT AND OUTPUT*********************************
C real*8 dmnew(maxno,maxo)   : Density Matrix
C                              Input:  d.m. output in current SCF iter
C                              Output: d.m. input for next SCF iteration
C real*8 dmold(maxno,maxo)   : Density matrix
C                              Input:  d.m. input in current SCF step
C                              Output: d.m. input for next SCF iteration
C ************************** OUTPUT *************************************
C real*8 dmax                : Maximum change of a DM element between 
C                              input and output
C ************************ BEHAVIOUR ************************************
C All data are saved in tape with direct access & unformatted form
C
C Algorithm changed!! (2/19/99)
C Pulay mixing is now done EVERY iteration, except the first one
C (which is done with linear mixing or with no mixing, according 
C to logical variable mix) 
C Mixes the last maxsav iterations (normal linear mixing if maxsav =< 1)
C 
C The mixing is done the following way (Anderson, Pulay):
C
C Modified input and output matrices from the former step are obtained
C by mixing maxsav prefious steps:
C
C   D'_in (n) = Sum_i=1,maxsav  beta_i D_in (n-maxsav+i)
C   D'_out(n) = Sum_i=1,maxsav  beta_i D_out(n-maxsav+i)
C
C The beta coefficients are obtained by minimizing the norm between
C D'_in and D'_out.
C The input charge for step (n+1) is done by simple mixing the D's
C 
C   D_in (n+1) = alpha D'_out(n) + (1-alpha) D_in(n)
C              = Sum_i=1,maxsav D_in (n-maxsav+i) + 
C                alpha Sum_i=1,maxsav deltaD(n-maxsav+i)
C
C where deltaD(n) is the residual of step n.
C
C The density matrices of BOTH spins are mixed at the same
C time, with the same coefficients (to ensure conservation of
C total number of electrons).
C ***********************************************************************

       implicit none

       integer 
     .  dimaux,iscf,maxsav,maxo,maxno,nbasis,nkick,nspin

       integer  
     .  numd(maxo)

       real*8 
     .  alpha,alphakick,dmax,
     .  dmnew(maxno,maxo,nspin),dmold(maxno,maxo,nspin)

       real*8
     .  savedm(dimaux), savere(dimaux)

       logical
     .  pulfile, mix

       character
     .  paste*33

       external
     .  io_assign, io_close, paste

       include 'fdf/fdfdefs.h'

C Internal variables ....................................................
       integer
     .  maxmix,mxddim

       parameter(mxddim=6)

       integer
     .  i0,i,ii,in,is,isite,j,jj,jtape,jtap1,numel

       real*8 
     .  b(mxddim,mxddim),bi(mxddim,mxddim),
     .  coeff(mxddim), sum

       character
     .  fname1*33, fname2*33, sname*30

       logical
     .  frstme

       save coeff, sname, frstme

       data  frstme /.true./
C ........................

       if (frstme) then
         sname = fdf_string('SystemLabel','siesta')
         frstme = .false.
       endif

C Check some input and dimensions .......................................
       if (iscf .eq. 1) then
         if (maxsav+1 .gt. mxddim) then
           write(6,*) 'pulayx: maxsav+1 > mxddim'
           stop
         endif
       endif

       numel = 0
       do is = 1,nspin
       do i = 1,nbasis
       do j = 1,numd(i)
         numel = numel + 1
       enddo
       enddo
       enddo

       if (.not. pulfile) then
         if (dimaux .lt. numel*maxsav) then
           write(6,*) 'pulayx: dimaux too small'
           stop
         endif
       endif
         
C ........................

       if (maxsav .gt. 1) then
C Open direct access files ..............................................
         fname1 = paste(sname,'.P1')
         fname2 = paste(sname,'.P2')

         if (pulfile) then
           call io_assign(jtape)
           call io_assign(jtap1)
           open(unit=jtape,file=fname1,form='unformatted',
     1     access='direct',recl=8*numel,status='unknown')
           open(unit=jtap1,file=fname2,form='unformatted',
     1     access='direct',recl=8*numel,status='unknown')
         endif
C ........................

C Write current D_in and Residual on tape ................................
         isite = mod(iscf,maxsav)
         if (isite .eq. 0) isite = maxsav
         if (pulfile) then
           write(jtape,rec=isite) 
     .      (((dmold(j,i,is),j=1,numd(i)),i=1,nbasis),is=1,nspin)
           write(jtap1,rec=isite) 
     .      ((((dmnew(j,i,is)-dmold(j,i,is)),j=1,numd(i)),
     .                                     i=1,nbasis),is=1,nspin)
         else
           i0 = (isite-1) * numel
           do is = 1,nspin
           do i = 1,nbasis
           do j = 1,numd(i)
             i0 = i0 + 1
             savedm(i0) = dmold(j,i,is)
             savere(i0) = dmnew(j,i,is) - dmold(j,i,is)
           enddo
           enddo
           enddo
         endif
       endif
         
C ........................


C Perform linear mixing if maxsav =< 1
       if (maxsav.le.1) then
         dmax = 0.0d0
         do is = 1,nspin
           do i = 1,nbasis
             do in = 1,numd(i)
               dmax = max(dmax, abs(dmnew(in,i,is) - dmold(in,i,is)))
               if (iscf .gt. 1 .or. mix) then
                 dmnew(in,i,is) = 
     .           (1.0d0-alpha)*dmold(in,i,is) + alpha*dmnew(in,i,is)
               endif
               dmold(in,i,is) = dmnew(in,i,is)
             enddo
           enddo
         enddo
         return
       endif
C .......................

C Pulay mixing otherwise... (except if iscf=1, or if iscf is multiple of nkick)

C Perform linear mixing if iscf = 1
       if (iscf.eq.1) then
         dmax = 0.0d0
         do is = 1,nspin
           do i = 1,nbasis
             do in = 1,numd(i)
               dmax = max(dmax, abs(dmnew(in,i,is) - dmold(in,i,is)))
               if (mix) then
                 dmnew(in,i,is) = 
     .           (1.0d0-alpha)*dmold(in,i,is) + alpha*dmnew(in,i,is)
               endif
               dmold(in,i,is) = dmnew(in,i,is)
             enddo
           enddo
         enddo
         if (pulfile) then
           call io_close(jtape)
           call io_close(jtap1)
         endif
         return
       endif
C .......................

C Perform linear mixing if iscf = N x nkick
       if (nkick .le. 0) goto 100
       if (mod(iscf,nkick).eq.0) then
         dmax = 0.0d0
         do is = 1,nspin
           do i = 1,nbasis
             do in = 1,numd(i)
               dmax = max(dmax, abs(dmnew(in,i,is) - dmold(in,i,is)))
               if (mix) then
                 dmnew(in,i,is) = 
     .      (1.0d0-alphakick)*dmold(in,i,is) + alphakick *dmnew(in,i,is)
               endif
               dmold(in,i,is) = dmnew(in,i,is)
             enddo
           enddo
         enddo
         if (pulfile) then
           call io_close(jtape)
           call io_close(jtap1)
         endif
         return
       endif
100    continue
C .......................

C Perform Pulay mixing if iscf .gt. 2

C  Compute current maximum deviation ...........
       dmax = 0.0d0
       do is = 1,nspin
         do i = 1,nbasis
           do in = 1,numd(i)
             dmax = max(dmax, abs(dmnew(in,i,is) - dmold(in,i,is)))
           enddo
         enddo
       enddo
C .......

C  calculate mixing coefficients, only if mixing the Density Matrix ........

       maxmix=maxsav
       if (iscf.lt.maxsav) maxmix=iscf

       do i=1,maxmix
         if (pulfile) then
           read(jtap1,rec=i) 
     .     (((dmnew(jj,ii,is),jj=1,numd(ii)),ii=1,nbasis),is=1,nspin)
         else
           i0 = (i-1) * numel
           do is = 1,nspin
           do ii = 1,nbasis
           do jj = 1,numd(ii)
             i0 = i0 + 1
             dmnew(jj,ii,is) = savere(i0)
           enddo
           enddo
           enddo
         endif

         b(i,i) = 0.0d0
         sum=0.0d0
         do is=1,nspin
         do ii=1,nbasis
         do jj=1,numd(ii)
           sum=sum+dmnew(jj,ii,is)*dmnew(jj,ii,is)
         enddo
         enddo
         enddo
         b(i,i)=sum

         do j=1,i-1
           if (pulfile) then
             read(jtap1,rec=j) 
     .     (((dmold(jj,ii,is),jj=1,numd(ii)),ii=1,nbasis),is=1,nspin)
           else
             i0 = (j-1) * numel
             do is = 1,nspin
             do ii = 1,nbasis
             do jj = 1,numd(ii)
               i0 = i0 + 1
               dmold(jj,ii,is) = savere(i0)
             enddo
             enddo
             enddo
           endif

           b(i,j)=0.0d0
           sum=0.0d0
           do is=1,nspin
           do ii=1,nbasis
           do jj=1,numd(ii)
             sum=sum+dmold(jj,ii,is)*dmnew(jj,ii,is)
           enddo
           enddo
           enddo
           b(i,j)=sum
           b(j,i)=sum
         enddo
         b(i,maxmix+1)=1.0
         b(maxmix+1,i)=1.0
       enddo

       b(maxmix+1,maxmix+1)=0.0

       call inver2(b,bi,maxmix+1,mxddim)

       do i=1,maxmix
         coeff(i)=bi(i,maxmix+1)
       enddo
C ........
 
C Read former matrices for mixing .........
       do is=1,nspin
       do ii=1,nbasis
       do j=1,numd(ii)
         dmnew(j,ii,is)=0.0d0
       enddo
       enddo
       enddo
       do i=1,maxmix
         if (pulfile) then
           read(jtape,rec=i) 
     .     (((dmold(j,ii,is),j=1,numd(ii)),ii=1,nbasis),is=1,nspin)
         else
           i0 = (i-1) * numel
           do is = 1,nspin
           do ii = 1,nbasis
           do j = 1,numd(ii)
             i0 = i0 + 1
             dmold(j,ii,is) = savedm(i0)
           enddo
           enddo
           enddo
         endif

         do is=1,nspin
         do ii=1,nbasis
         do j=1,numd(ii)
           dmnew(j,ii,is)=dmnew(j,ii,is)+coeff(i)*dmold(j,ii,is)
         enddo
         enddo
         enddo
       enddo

       do i=1,maxmix
         if (pulfile) then
           read(jtap1,rec=i) 
     .     (((dmold(j,ii,is),j=1,numd(ii)),ii=1,nbasis),is=1,nspin)
         else
           i0 = (i-1) * numel
           do is = 1,nspin
           do ii = 1,nbasis
           do j = 1,numd(ii)
             i0 = i0 + 1
             dmold(j,ii,is) = savere(i0)
           enddo
           enddo
           enddo
         endif

         do is=1,nspin
         do ii=1,nbasis
         do j=1,numd(ii)
           dmnew(j,ii,is)=dmnew(j,ii,is)+
     .           alpha*coeff(i)*dmold(j,ii,is)
         enddo
         enddo
         enddo
       enddo

       do is=1,nspin
       do ii=1,nbasis
       do j=1,numd(ii)
         dmold(j,ii,is)=dmnew(j,ii,is)
       enddo
       enddo
       enddo
C ........
       if (pulfile) then
         call io_close(jtape)
         call io_close(jtap1)
       endif

       return
       end


        SUBROUTINE INVER2(A,B,N,mxddim)
        implicit real*8  (a-h,o-z)
        real*8 A(mxddim,mxddim),B(mxddim,mxddim),X
        DO 20 I=1,N
        DO 20 J=1,N
20      B(I,J)=A(I,J)
        DO 4 I=1,N
        X=B(I,I)
        B(I,I)=1.
        DO 1 J=1,N
1       B(J,I)=B(J,I)/X
        DO 4 K=1,N
        IF(K-I) 2,4,2
2       X=B(I,K)
        B(I,K)=0.0d0
        DO 3 J=1,N
3       B(J,K)=B(J,K)-B(J,I)*X
4       CONTINUE
        RETURN
        END
