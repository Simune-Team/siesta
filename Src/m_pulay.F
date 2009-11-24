! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module m_pulay
      use precision, only: dp

      implicit none

      private

      integer            :: nauxpul = 1  ! Max. size of auxiliary 
                                         ! sparse Pulay matrices
      real(dp), pointer  :: auxpul(:,:)  ! Auxilialy matrices 
                                         ! for Pulay mixing

      public :: pulayx, init_pulay_arrays, resetPulayArrays

      CONTAINS

      subroutine init_pulay_arrays()
      use siesta_options, only: pulfile, maxsav
      use alloc
      use atomlist, only: no_l
      use m_spin, only: nspin
      use m_mpi_utils, only: globalize_max
      use sparse_matrices, only: numh

      implicit none

      integer :: ntmp
      logical, save  :: first_time = .true.

      if (first_time) then
         nullify(auxpul)
         first_time = .false.
      endif

      if (pulfile .or. maxsav .le. 0) then
        ! No need for auxiliary arrays
        nauxpul = 1
        call re_alloc( auxpul, 1, nauxpul, 1, 2, 'auxpul', 'pulay' )
      else
        nauxpul = sum(numh(1:no_l)) * nspin * maxsav
#ifdef MPI
        call globalize_max(nauxpul,ntmp)
        nauxpul = ntmp
#endif
! Increase nauxpul with safety margin when atoms move  ! ???
        nauxpul = 1.1 * nauxpul + 10
        call re_alloc( auxpul, 1, nauxpul, 1, 2, 'auxpul', 'pulay' )
      endif
      end subroutine init_pulay_arrays
!
      subroutine pulayx(pulfile,iscf,mix,nbasis,maxnd,numd,
     .                   listdptr,nspin,maxsav,alpha,nkick,alphakick,
     .                   dmnew,dmold,dmax)
C ***********************************************************************
C Pulay mixing implemented to accelerate the self-consistency
C Mixes MAXSAV previous steps.  Linear mixing if MAXAV =< 0
C 
C Written by In-Ho Lee, Beckman Inst., Univ. of Illinois, Mar. 25 '97
C Modified and partly re-written by P. Ordejon, July'97
C Modified and optimized by P. Ordejon, November'97
C
C NOTE : cannot use pulfile=.true. yet in parallel!!!
C
C ************************** INPUT **************************************
C logical pulfile            : Use file to store Pulay info
C                                 .true.  = use file
C                                 .false. = use memory
C integer iscf               : Current SCF iteration
C logical mix                : Mix first SCF step (T or F)
C integer nbasis             : Number of atomic orbitals stored locally
C integer maxnd              : First dimension of D.M., and 
C                              maximum number of nonzero elements of D.M.
C integer numd(nbasis)       : Control vector of D.M.
C                              (number of nonzero elements of each row)
C integer listdptr(nbasis)   : Pointer to start of rows in listd
C integer nspin              : Spin polarization (1=unpolarized, 2=polarized)
C integer maxsav             : Pulay mixing is done every maxsav iterations.
C                              Remaining iterations are done by linear mixing.
C real*8 alpha               : Mixing parameter (for linear mixing)
C integer nkick              : Do a linear mixing (kick) each nkick cycles
C real*8 alphakick           : Mixing parameter for kick cycles

C ********************* INPUT AND OUTPUT*********************************
C real*8 dmnew(maxnd)        : Density Matrix
C                              Input: d.m. output in current SCF step
C                              Output: d.m. input for next SCF iteration
C real*8 dmold(maxnd)        : Density matrix
C                              Input: d.m. input in current SCF step
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

C
C  Modules
C
      use precision
      use parallel,   only : Node, Nodes
      use files,      only : slabel, label_length
      use sys,        only : die
      use alloc
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      integer , intent(in) ::  iscf,maxsav,maxnd,
     $                         nbasis,nkick,nspin

      integer , intent(in) :: numd(*),listdptr(*)
      logical, intent(in) ::  pulfile, mix

      real(dp), intent(in)    ::  alpha,alphakick
      real(dp), intent(inout) ::  dmnew(maxnd,nspin),dmold(maxnd,nspin)
      real(dp), intent(out)   :: dmax


      real(dp), pointer :: savedm(:), savere(:)

      external  io_assign, io_close, memory

C Internal variables ....................................................
      integer
     .  i0,i,ii,in,is,isite,j,jj,jtape,jtap1,numel,ind,info,maxmix

#ifdef MPI
      integer
     .  MPIerror
#endif

      integer  :: real_length        
      real(dp) ::   ssum
      real(dp), dimension(:,:), pointer ::  b, bi
      real(dp), dimension(:), pointer   ::  buffer
      real(dp), dimension(:), pointer   ::  coeff

      character(len=label_length+3) :: fname1
      character(len=label_length+3) :: fname2

      savedm => auxpul(:,1)         ! DM in former iterations
      savere => auxpul(:,2)         ! Residuals in former iterations

#ifdef MPI
C Check that this is not a parallel run with pulfile=.true.
       if (pulfile.and.Nodes.gt.1)
     $     call die('pulayx: cannot use pulay on file in parallel')
#endif

C Check some input and dimensions 
       numel = 0
       do i = 1,nbasis
         numel = numel + numd(i)
       enddo
       numel = numel*nspin

       if (.not.pulfile) then
         if (size(savere) .lt. numel*maxsav)
     $         call die('pulayx: dimensions are too small')
       endif
         

       if (maxsav .gt. 1) then
C Open direct access files
         fname1 = trim(slabel) // '.P1'
         fname2 = trim(slabel) // '.P2'

         if (pulfile) then
           call io_assign(jtape)
           call io_assign(jtap1)
           inquire(iolength=real_length) savere(1)
           open(unit=jtape,file=fname1,form='unformatted',
     .       access='direct',recl=real_length*numel,status='unknown')
           open(unit=jtap1,file=fname2,form='unformatted',
     .       access='direct',recl=real_length*numel,status='unknown')
         endif
C ........................

C Write current D_in and Residual on tape ................................
         isite = mod(iscf,maxsav)
         if (isite .eq. 0) isite = maxsav
         if (pulfile) then
           write(jtape,rec=isite) 
     .      (((dmold(listdptr(i)+j,is),j=1,numd(i)),i=1,nbasis),
     .      is=1,nspin)
           write(jtap1,rec=isite) 
     .      ((((dmnew(listdptr(i)+j,is)-dmold(listdptr(i)+j,is)),
     .      j=1,numd(i)),i=1,nbasis),is=1,nspin)
         else
           i0 = (isite-1) * numel
           do is = 1,nspin
             do i = 1,nbasis
               do j = 1,numd(i)
                 i0 = i0 + 1
                 savedm(i0) = dmold(listdptr(i)+j,is)
                 savere(i0) = dmnew(listdptr(i)+j,is) - 
     .                        dmold(listdptr(i)+j,is)
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
               ind = listdptr(i) + in
               dmax = max(dmax, abs(dmnew(ind,is) - dmold(ind,is)))
               if (iscf .gt. 1 .or. mix) then
                 dmnew(ind,is) =
     .           (1.0d0-alpha)*dmold(ind,is) + alpha*dmnew(ind,is)
               endif
               dmold(ind,is) = dmnew(ind,is)
             enddo
           enddo
         enddo
         return
       endif

C Pulay mixing otherwise... (except if iscf=1, or if iscf is multiple of nkick)

C Perform linear mixing if iscf = 1
       if (iscf.eq.1) then
         dmax = 0.0d0
         do is = 1,nspin
           do i = 1,nbasis
             do in = 1,numd(i)
               ind = listdptr(i) + in
               dmax = max(dmax, abs(dmnew(ind,is) - dmold(ind,is)))
               if (mix) then
                 dmnew(ind,is) =
     .           (1.0d0-alpha)*dmold(ind,is) + alpha*dmnew(ind,is)
               endif
               dmold(ind,is) = dmnew(ind,is)
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

       if (nkick> 0 .AND. mod(iscf,nkick).eq.0) then
         dmax = 0.0d0
         do is = 1,nspin
           do i = 1,nbasis
             do in = 1,numd(i)
               ind = listdptr(i) + in
               dmax = max(dmax, abs(dmnew(ind,is) - dmold(ind,is)))
               if (mix) then
                 dmnew(ind,is) =
     .      (1.0d0-alphakick)*dmold(ind,is) + alphakick *dmnew(ind,is)
               endif
               dmold(ind,is) = dmnew(ind,is)
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

C Perform Pulay mixing if iscf .gt. 2

C Allocate local arrays

        nullify( b )
        call re_alloc( b, 1, maxsav+1, 1, maxsav+1, 'b', 'pulay' )
        nullify( bi )
        call re_alloc( bi, 1, maxsav+1, 1, maxsav+1, 'bi', 'pulay' )
        nullify( coeff )
        call re_alloc( coeff, 1, maxsav+1, 'coeff', 'pulay' )
        nullify( buffer )
        call re_alloc( buffer, 1, maxsav, 'buffer', 'pulay' )

C  Compute current maximum deviation ...........
       dmax = 0.0d0
       do is = 1,nspin
         do i = 1,nbasis
           do in = 1,numd(i)
             ind = listdptr(i) + in
             dmax = max(dmax, abs(dmnew(ind,is) - dmold(ind,is)))
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
     .     (((dmnew(listdptr(ii)+jj,is),jj=1,numd(ii)),ii=1,nbasis),
     .     is=1,nspin)
         else
           i0 = (i-1) * numel
           do is = 1,nspin
             do ii = 1,nbasis
               do jj = 1,numd(ii)
                 ind = listdptr(ii) + jj
                 i0 = i0 + 1
                 dmnew(ind,is) = savere(i0)
               enddo
             enddo
           enddo
         endif

         b(i,i) = 0.0d0
         ssum=0.0d0
         do is=1,nspin
           do ii=1,nbasis
             do jj=1,numd(ii)
               ind = listdptr(ii) + jj
               ssum=ssum+dmnew(ind,is)*dmnew(ind,is)
             enddo
           enddo
         enddo
         b(i,i)=ssum

         do j=1,i-1
           if (pulfile) then
             read(jtap1,rec=j) 
     .       (((dmold(listdptr(ii)+jj,is),jj=1,numd(ii)),ii=1,nbasis),
     .         is=1,nspin)
           else
             i0 = (j-1) * numel
             do is = 1,nspin
               do ii = 1,nbasis
                 do jj = 1,numd(ii)
                   ind = listdptr(ii) + jj
                   i0 = i0 + 1
                   dmold(ind,is) = savere(i0)
                 enddo
               enddo
             enddo
           endif

           b(i,j)=0.0d0
           ssum=0.0d0
           do is=1,nspin
             do ii=1,nbasis
               do jj=1,numd(ii)
                 ind = listdptr(ii) + jj
                 ssum=ssum+dmold(ind,is)*dmnew(ind,is)
               enddo
             enddo
           enddo
           b(i,j)=ssum
           b(j,i)=ssum
         enddo
         b(i,maxmix+1)=1.0d0
         b(maxmix+1,i)=1.0d0
       enddo

       b(maxmix+1,maxmix+1)=0.0d0

#ifdef MPI
C Global operations
       do i=1,maxmix
         call MPI_AllReduce(b(1,i),buffer,maxmix,MPI_double_precision,
     .     MPI_sum,MPI_Comm_World,MPIerror)
         do j=1,maxmix
           b(j,i)=buffer(j)
         enddo
       enddo
#endif

       call inver(b,bi,maxmix+1,maxsav+1,info)

C If inver was successful, get coefficients for Pulay mixing
       if (info .eq. 0) then
         do i=1,maxmix
           coeff(i)=bi(i,maxmix+1)
         enddo
       else
C Otherwise, use only last step
         do i=1,maxmix
           coeff(i)=0.0d0
         enddo
         j=mod(iscf,maxmix)
         if(j.eq.0) j=maxmix
         coeff(j) = 1.0d0
       endif

C ........
 
C Read former matrices for mixing .........
       dmnew(1:maxnd,1:nspin)=0.0d0
       do i=1,maxmix
         if (pulfile) then
           read(jtape,rec=i) 
     .     (((dmold(listdptr(ii)+j,is),j=1,numd(ii)),ii=1,nbasis),
     .       is=1,nspin)
         else
           i0 = (i-1) * numel
           do is = 1,nspin
             do ii = 1,nbasis
               do j = 1,numd(ii)
                 ind = listdptr(ii) + j
                 i0 = i0 + 1
                 dmold(ind,is) = savedm(i0)
               enddo
             enddo
           enddo
         endif

         do is=1,nspin
           do ii=1,nbasis
             do j=1,numd(ii)
               ind = listdptr(ii) + j
               dmnew(ind,is)=dmnew(ind,is)+dmold(ind,is)*coeff(i)
             enddo
           enddo
         enddo
       enddo

       do i=1,maxmix
         if (pulfile) then
           read(jtap1,rec=i)
     .     (((dmold(listdptr(ii)+j,is),j=1,numd(ii)),ii=1,nbasis),
     .       is=1,nspin)
         else
           i0 = (i-1) * numel
           do is = 1,nspin
             do ii = 1,nbasis
               do j = 1,numd(ii)
                 ind = listdptr(ii) + j
                 i0 = i0 + 1
                 dmold(ind,is) = savere(i0)
               enddo
             enddo
           enddo
         endif

         do is=1,nspin
           do ii=1,nbasis
             do j=1,numd(ii)
               ind = listdptr(ii) + j
               dmnew(ind,is)=dmnew(ind,is)+
     .           alpha*coeff(i)*dmold(ind,is)
             enddo
           enddo
         enddo
       enddo

       ! Test in case processor has no orbitals.
       ! 
       if (nbasis>0) then
         do is=1,nspin
           do ii=1,listdptr(nbasis)+numd(nbasis)
             dmold(ii,is)=dmnew(ii,is)
           enddo
         enddo
       endif
C ........
       if (pulfile) then
         call io_close(jtape)
         call io_close(jtap1)
       endif

C Deallocate local arrays

       call de_alloc( b, 'b', 'pulay' )
       call de_alloc( bi, 'bi', 'pulay' )
       call de_alloc( coeff, 'coeff', 'pulay' )
       call de_alloc( buffer, 'buffer', 'pulay' )

       end subroutine pulayx

      subroutine resetPulayArrays( )
      use alloc
      implicit none
      if (associated(auxpul))
     &  call de_alloc( auxpul, 'auxpul', 'pulay' )
      end subroutine resetPulayArrays
      
      end module m_pulay

