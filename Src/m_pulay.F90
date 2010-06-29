module m_pulay
  use precision, only: dp
  !
  implicit none
  !
  private
  !
  real(dp), pointer  :: auxpul(:,:)  => null()
  integer            :: n_records_saved = 0
  !
  public :: pulayx, init_pulay_arrays
  !
CONTAINS
  !
  subroutine init_pulay_arrays()
    use siesta_options, only: maxsav
    use alloc
    use atomlist, only: no_l
    use m_spin, only: nspin
    use sparse_matrices, only: numh
    !
    implicit none
    !
    integer :: ntmp, nauxpul
    !
    n_records_saved = 0

    if (maxsav .le. 0) then
       ! No need for auxiliary arrays
    else
       nauxpul = sum(numh(1:no_l)) * nspin * maxsav
       !
       call re_alloc(auxpul,1,nauxpul,1,2,name="auxpul",        &
            routine="init_pulay_arrrays")
    endif
    !
  end subroutine init_pulay_arrays
  !---------------------------------------------------------------------

  subroutine pulayx( iscf,mix1,no_l,maxnd,numd,             &
       listdptr,nspin,maxsav,alpha,nkick,alphakick,   &
       dmnew,dmold,dmax)

    ! Pulay mixing implemented to accelerate the self-consistency
    ! Mixes MAXSAV previous steps.  Linear mixing if MAXAV =< 0
    ! 
    ! Written by In-Ho Lee, Beckman Inst., Univ. of Illinois, Mar. 25 '97
    ! Modified and partly re-written by P. Ordejon, July'97
    ! Modified and optimized by P. Ordejon, November'97
    ! Modified by A. Garcia, April-June 2010
    !
    ! ************************** INPUT **************************************
    ! integer iscf               : Current SCF iteration
    ! logical mix1               : Mix first SCF step (T or F)
    ! integer no_l               : Number of atomic orbitals stored locally
    ! integer maxnd              : First dimension of D.M., and 
    !                              maximum number of nonzero elements of D.M.
    ! integer numd(no_l)         : Control vector of D.M.
    !                              (number of nonzero elements of each row)
    ! integer listdptr(no_l)     : Pointer to start of rows in listd
    ! integer nspin              : Spin polarization (1= no spin, 2=polarized)
    ! integer maxsav             : Number of history records for Pulay mixing.
    ! real*8 alpha               : Mixing parameter (for linear mixing)
    ! integer nkick              : Do a linear mixing (kick) each nkick cycles
    ! real*8 alphakick           : Mixing parameter for kick cycles
    !
    ! ********************* INPUT AND OUTPUT*********************************
    ! real*8 dmnew(maxnd)        : Density Matrix
    !                              Input: d.m. output in current SCF step
    !                              Output: d.m. input for next SCF iteration
    ! real*8 dmold(maxnd)        : Density matrix
    !                              Input: d.m. input in current SCF step
    !                              Output: d.m. input for next SCF iteration
    ! ************************** OUTPUT *************************************
    ! real*8 dmax                : Maximum change of a DM element between 
    !                              input and output
    ! ************************ BEHAVIOUR ************************************
    !
    ! Algorithm changed!! (2/19/99)
    ! Pulay mixing is now done EVERY iteration, except the first one
    ! (which is done with linear mixing or with no mixing, according
    ! to logical variable mix1)
    ! Mixes the last maxsav iterations (normal linear mixing if maxsav =< 1)
    !
    ! The mixing is done the following way (Anderson, Pulay):
    !
    ! Modified input and output matrices from the former step are obtained
    ! by mixing maxsav prefious steps:
    !
    !   D'_in (n) = Sum_i=1,maxsav  beta_i D_in (n-maxsav+i)
    !   D'_out(n) = Sum_i=1,maxsav  beta_i D_out(n-maxsav+i)
    !
    ! The beta coefficients are obtained by minimizing the norm between
    ! D'_in and D'_out.
    ! The input charge for step (n+1) is done by simple mixing the D's
    !
    !   D_in (n+1) = alpha D'_out(n) + (1-alpha) D_in(n)
    !              = Sum_i=1,maxsav D_in (n-maxsav+i) +
    !                alpha Sum_i=1,maxsav deltaD(n-maxsav+i)
    !
    ! where deltaD(n) is the residual of step n.
    !
    ! The density matrices of BOTH spins are mixed at the same
    ! time, with the same coefficients (to ensure conservation of
    ! total number of electrons).
    ! (and for NONCOL?)
    ! ***********************************************************************
    !
    !
    !  Modules
    !
    use precision,  only : dp
    use parallel,   only : Node
    use sys,        only : die
    use alloc
    use siesta_options, only: avoid_first_after_kick
#ifdef MPI
    use mpi_siesta
#endif
    !
    implicit none
    !
    integer , intent(in) ::  iscf,maxsav,maxnd, no_l,nkick,nspin
    !
    integer , intent(in) :: numd(*),listdptr(*)
    logical, intent(in) ::  mix1
    !
    real(dp), intent(in)    ::  alpha,alphakick
    real(dp), intent(inout) ::  dmnew(maxnd,nspin),dmold(maxnd,nspin)
    real(dp), intent(out)   :: dmax
    !
    !
    real(dp), pointer :: savedm(:), saveres(:)
    !
    ! Internal variables ....................................................
    integer :: i0,i,ii,in,is,isite,j,jj,numel,ind,info,maxmix
    logical :: after_kick
    !
#ifdef MPI
    integer  MPIerror
#endif
    !
    real(dp) ::   ssum
    real(dp), dimension(:,:), pointer ::  b, bi
    real(dp), dimension(:), pointer   ::  buffer
    real(dp), dimension(:), pointer   ::  coeff
    !

    if (maxsav > 1) then
       savedm => auxpul(:,1)          ! DM in former iterations
       saveres => auxpul(:,2)         ! Residuals in former iterations
       !
       ! Check some input and dimensions 
       numel = sum(numd(1:no_l))
       call assert(numel == maxnd, "sum(numd) // maxnd")
       numel = numel*nspin
       call assert(numel*maxsav == size(savedm), "maxnd*spin*maxsav // savedm")

       !
       if (size(saveres) .lt. numel*maxsav)                            &
            call die('pulayx: dimensions are too small')
    endif

    ! For the purposes of history keeping, the very first iteration
    ! is treated as the first step after a kick: a possibly large
    ! jump in the DM

    after_kick = (iscf == 1 .OR.     &
                 (nkick > 0 .AND. mod(iscf,nkick) == 1) ) 

    if (after_kick .and. avoid_first_after_kick) then

       ! Do not keep the residual of the very first iteration, or of 
       ! the first iteration after a kick, if so instructed

    else
       call WriteCurrentDinAndResidualInStore()
       ! updates n_records_saved
    endif

    !  Compute current maximum deviation ...........
    dmax = 0.0_dp
    do is = 1,nspin
       do i = 1,no_l
          do in = 1,numd(i)
             ind = listdptr(i) + in
             dmax = max(dmax, abs(dmnew(ind,is) - dmold(ind,is)))
          enddo
       enddo
    enddo
    ! .......
    !
    ! Perform linear mixing if we do not have enough history

    if (n_records_saved <= 1) then

       ! Could do simply:  
       ! if (iscf > 1 .or. mix1) DMnew = (1.0_dp-alpha)*DMold + alpha*DMnew
       ! DMold = DMnew

       do is = 1,nspin
          do i = 1,no_l
             do in = 1,numd(i)
                ind = listdptr(i) + in
                if (iscf .gt. 1 .or. mix1) then
                   dmnew(ind,is) =                                        &
                        (1.0d0-alpha)*dmold(ind,is) + alpha*dmnew(ind,is)
                endif
                dmold(ind,is) = dmnew(ind,is)
             enddo
          enddo
       enddo
       RETURN

    endif
    !
    ! Perform linear mixing if iscf = N x nkick
    !
    if (nkick > 0 .AND. mod(iscf,nkick).eq.0) then
       ! Reset the history information
       if (Node == 0) then
          write(6,"(a,i4)") "Applying kick at iscf: ", iscf
          write(6,"(a,i4)") "Re-setting Pulay history"
       endif
       n_records_saved = 0

       do is = 1,nspin
          do i = 1,no_l
             do in = 1,numd(i)
                ind = listdptr(i) + in
                dmnew(ind,is) =                                        &
                     (1.0d0-alphakick)*dmold(ind,is) + alphakick *dmnew(ind,is)
                dmold(ind,is) = dmnew(ind,is)
             enddo
          enddo
       enddo
       RETURN
    endif
    !
    ! .......................
    !
    ! Perform Pulay mixing if n_records_saved > 1
    !
    call assert(n_records_saved > 1, "N_SCF_records <= 1")
    ! Use only as many records as we have in the history store
    maxmix= min(n_records_saved, maxsav)

    ! Allocate local arrays
    !
    nullify( b )
    call re_alloc( b, 1, maxmix+1, 1, maxmix+1, name='b',           &
         routine='pulayx' )
    nullify( bi )
    call re_alloc( bi, 1, maxmix+1, 1, maxmix+1, name='bi',         &
         routine='pulayx' )
    nullify( coeff )
    call re_alloc( coeff, 1, maxmix+1, name='coeff',                &
         routine='pulayx' )
    nullify( buffer )
    call re_alloc( buffer, 1, maxmix, name='buffer',                &
         routine='pulayx' )
    !
    !
    !  calculate mixing coefficients
    !
    !
    do i=1,maxmix
       i0 = (i-1) * numel
       do is = 1,nspin
          do ii = 1,no_l
             do jj = 1,numd(ii)
                ind = listdptr(ii) + jj
                i0 = i0 + 1
                dmnew(ind,is) = saveres(i0)
             enddo
          enddo
       enddo
       !
       ! B(i,i) = dot_product(Res(i)*Res(i))
       b(i,i) = 0.0d0
       ssum=0.0d0
       do is=1,nspin
          do ii=1,no_l
             do jj=1,numd(ii)
                ind = listdptr(ii) + jj
                ssum=ssum+dmnew(ind,is)*dmnew(ind,is)
             enddo
          enddo
       enddo
       b(i,i)=ssum
       !
       do j=1,i-1

          i0 = (j-1) * numel
          do is = 1,nspin
             do ii = 1,no_l
                do jj = 1,numd(ii)
                   ind = listdptr(ii) + jj
                   i0 = i0 + 1
                   dmold(ind,is) = saveres(i0)
                enddo
             enddo
          enddo

          !          ! B(i,j) = B(j,i) = dot_product(Res(i)*Res(j))

          b(i,j)=0.0d0
          ssum=0.0d0
          do is=1,nspin
             do ii=1,no_l
                do jj=1,numd(ii)
                   ind = listdptr(ii) + jj
                   ssum=ssum+dmold(ind,is)*dmnew(ind,is)
                enddo
             enddo
          enddo
          b(i,j)=ssum
          b(j,i)=ssum
       enddo

       ! Now extend the matrix with ones in an extra colum
       ! and row ...
       b(i,maxmix+1)=1.0d0
       b(maxmix+1,i)=1.0d0
    enddo
    !      ! ... except in the extra diagonal entry
    b(maxmix+1,maxmix+1)=0.0d0
    !
#ifdef MPI
    ! Global operations, but only for the non-extended entries
    do i=1,maxmix
       call MPI_AllReduce(b(1:maxmix,i),buffer,maxmix,     &
            MPI_double_precision,            &
            MPI_sum,MPI_Comm_World,MPIerror)
       do j=1,maxmix
          b(j,i)=buffer(j)
       enddo
    enddo
#endif
    !
    call inver(b,bi,maxmix+1,maxmix+1,info)
    !
    ! If inver was successful, get coefficients for Pulay mixing
    if (info .eq. 0) then
       do i=1,maxmix
          coeff(i)=bi(i,maxmix+1)
       enddo
    else
       ! Otherwise, use only last step
       if (Node == 0) then
          write(6,"(a)") "FAILED inversion in Pulayx"
!!          write(0,"(a)") "FAILED inversion in Pulayx"
       endif
       do i=1,maxmix
          coeff(i)=0.0d0
       enddo
       j=mod(iscf,maxmix)
       if(j.eq.0) j=maxmix
       coeff(j) = 1.0d0
    endif
    !
    ! ........
    !
    ! Read former matrices for mixing .........
    dmnew(1:maxnd,1:nspin)=0.0d0
    do i=1,maxmix
       i0 = (i-1) * numel
       do is = 1,nspin
          do ii = 1,no_l
             do j = 1,numd(ii)
                ind = listdptr(ii) + j
                i0 = i0 + 1
                dmold(ind,is) = savedm(i0)
             enddo
          enddo
       enddo
       !
       do is=1,nspin
          do ii=1,no_l
             do j=1,numd(ii)
                ind = listdptr(ii) + j
                dmnew(ind,is)=dmnew(ind,is)+dmold(ind,is)*coeff(i)
             enddo
          enddo
       enddo

    enddo
    !
    do i=1,maxmix

       i0 = (i-1) * numel
       do is = 1,nspin
          do ii = 1,no_l
             do j = 1,numd(ii)
                ind = listdptr(ii) + j
                i0 = i0 + 1
                dmold(ind,is) = saveres(i0)
             enddo
          enddo
       enddo

       !
       do is=1,nspin
          do ii=1,no_l
             do j=1,numd(ii)
                ind = listdptr(ii) + j
                dmnew(ind,is) = dmnew(ind,is)  +     &
                     alpha*coeff(i)*dmold(ind,is)
             enddo
          enddo
       enddo

    enddo
    !
    ! Test in case processor has no orbitals.  ?????

    if (no_l>0) then
       do is=1,nspin
          do ii=1,listdptr(no_l)+numd(no_l)   ! 1, maxnd
             dmold(ii,is)=dmnew(ii,is)
          enddo
       enddo
    endif

    ! Deallocate local arrays
    !
    call de_alloc( b, name='b' )
    call de_alloc( bi, name="bi" )
    call de_alloc( coeff, name="coeff" )
    call de_alloc( buffer, name="buffer" )
    !
  CONTAINS

    subroutine WriteCurrentDinAndResidualInStore()

      ! The store is really a circular array of size maxsav. It is
      ! implemented inline here, but it could be done more cleanly
      ! with a finite circular stack.

      if (maxsav > 1) then
         n_records_saved = n_records_saved + 1
         ! isite marks the point to write 
         isite = mod(n_records_saved,maxsav)
         if (isite .eq. 0) isite = maxsav
         i0 = (isite-1) * numel
         if (Node == 0) then
            !write(6,"(a,i4)") "Saving scf history record. iscf: ", iscf
         endif

         ! Both spins are written, one block after the other,
         ! but the position on the array will not be conmmensurate
         ! for the different processors...
         ! ... that is why it is difficult to parallelize the on-file version
   
         do is = 1,nspin
            do i = 1,no_l
               do j = 1,numd(i)
                  i0 = i0 + 1
                  savedm(i0) = dmold(listdptr(i)+j,is)
                  saveres(i0) = dmnew(listdptr(i)+j,is) -                 &
                       dmold(listdptr(i)+j,is)
               enddo
            enddo
         enddo
      endif

    end subroutine WriteCurrentDinAndResidualInStore

    subroutine assert(condition, message)
    ! A simple assertion checker
      logical, intent(in)          :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         if (Node == 0) then
            write(6,"(2a)") "** Assertion failed **: ", message
            write(0,"(2a)") "** Assertion failed **: ", message
         endif
      endif
    end subroutine assert

  end subroutine pulayx
  !
end module m_pulay
!
