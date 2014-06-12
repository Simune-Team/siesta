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
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.
! * It has been heavily inspired by the original authors of the 
!   Transiesta code (hence the references here are still remaining) *

module m_ts_mumps_init

  use precision, only : dp

  implicit none

#ifdef MUMPS  

  integer, public, save :: MUMPS_mem = 20
  integer, public, save :: MUMPS_ordering = 7
  integer, public, save :: MUMPS_block = -8 ! blocking factor

  public :: init_MUMPS
  public :: analyze_MUMPS
  public :: prep_LHS
  public :: prep_RHS_Eq
  public :: prep_RHS_nEq
  public :: insert_Self_Energies

  private
  
contains

  subroutine init_MUMPS(mum,ID)
#ifdef MPI
    use mpi_siesta
#endif
    include 'zmumps_struc.h'
    type(zMUMPS_STRUC), intent(inout) :: mum
    integer, intent(in) :: ID
    character(len=25) :: file

    ! define processors and MATRIX-type
    mum%SYM =  0 ! non-symmetric
    mum%PAR =  1 ! sequential
    mum%COMM = MPI_COMM_SELF ! only communicate with it-self

    mum%JOB = -1 ! initialise
    call zMUMPS(mum)
    if ( mum%INFO(1) < 0 .or. mum%INFOG(1) < 0 ) then
       call die('MUMPS initialization had an error.')
    end if

    ! Setup the control parameters
    mum%ICNTL(5) = 0 ! assembled format
    ! The ordering of the matrix
    mum%ICNTL(7) = MUMPS_ordering

    ! We request a sparse right-hand side: 
    !   20 == 1, we need not initialize RHS_SPARSE!
    mum%ICNTL(20) = 1
    mum%ICNTL(21) = 0
    
    ! Allow memory increase handled by the user
    mum%ICNTL(14) = MUMPS_Mem

    ! Sets the blocking factor
    mum%ICNTL(27) = MUMPS_block

    ! Request specific elements of the inverse matrix: 
    !   30 == 1 we MUST allocate it, but need not initialize it
    mum%ICNTL(30) = 1
    ! this however, posses some other problems.

    ! For each processor, add their own output files.
    ! As they become quite big we make them rewritten 
    ! in every SCF
    write(file,'(a,i0,a)') 'TS_MUMPS_',ID,'.dat'
    call io_assign(mum%ICNTL(1))
    mum%ICNTL(2) = mum%ICNTL(1)
    mum%ICNTL(3) = mum%ICNTL(1)
    open(mum%ICNTL(1),file=trim(file),status='replace', &
         action='write',form='formatted')

  end subroutine init_MUMPS

  subroutine analyze_MUMPS(mum)
    include 'zmumps_struc.h'
    type(zMUMPS_STRUC), intent(inout) :: mum
    integer :: iu

    ! analyse the MUMPS solver, this will determine
    ! factorization strategy
    mum%JOB = 1
    call zMUMPS(mum)
    if ( mum%INFO(1) < 0 .or. mum%INFOG(1) < 0 ) then
       call die('MUMPS analysis step had an error.')
    end if

    ! Write out estimated memory requirements
    iu = mum%ICNTL(1)
    write(iu,'(/,a)')'### Memory estimation from ANALYSIS step...'
    write(iu,'(a,i0,a)')'### Minimum memory required for calculation: ', &
         mum%INFO(15),' MB'
    write(iu,'(a,i0,a)')'### MUMPS is allocating: ', &
         nint(mum%INFO(15)*real(100+mum%ICNTL(14),dp)/100._dp),' MB'
    write(iu,'(a,/)')"### MUMPS memory can be altered using TS.MUMPS.Mem."

  end subroutine analyze_MUMPS

  subroutine prep_LHS(mum,N_Elec,Elecs)
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_Self
#endif
    use class_OrbitalDistribution
    use class_Sparsity
    use m_ts_sparse, only : ts_sp_uc
    use m_ts_electype
    use m_ts_method, only : orb_offset
    use create_Sparsity_Union, only: crtSparsity_Union
    include 'zmumps_struc.h'

    type(zMUMPS_STRUC), intent(inout) :: mum
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)

    type(OrbitalDistribution) :: dit
    type(Sparsity) :: tmpSp1, tmpSp2
    integer :: iEl, idx, no, io, jo, ind, nr, ioff
    integer, pointer :: l_ptr(:),l_ncol(:), l_col(:)

#ifdef MPI
    call newDistribution(nrows_g(ts_sp_uc),MPI_Comm_Self,dit, &
         name='MUMPS UC distribution')
#else    
    call newDistribution(nrows_g(ts_sp_uc),-1,dit, &
         name='MUMPS UC distribution')
#endif

    ! This works as creating a new sparsity deletes the previous
    ! and as it is referenced several times it will not be actually
    ! deleted...
    tmpSp2 = ts_sp_uc
    do iEl = 1 , N_Elec

       idx = Elecs(iEl)%idx_no
       no = TotUsedOrbs(Elecs(iEl))

       ! we first create the super-set sparsity
       tmpSp1 = tmpSp2
       call crtSparsity_Union(dit,tmpSp1, &
            idx,idx,no,no, tmpSp2)
    end do
    call delete(tmpSp1)
    call delete(dit)

    ! Create the index 
    call attach(tmpSp2,list_ptr=l_ptr, &
         n_col=l_ncol,list_col=l_col,nrows_g=nr, &
         nnzs=mum%NZ)

    ! Allocate LHS, first we need to
    ! calculate the actual size of the matrix
    ! We know that it must be the Hamiltonian sparsity
    ! pattern + all Bulk-electrodes!
    allocate( mum%IRN ( mum%NZ ) )
    allocate( mum%JCN ( mum%NZ ) )
    allocate( mum%A   ( mum%NZ ) )

    do io = 1, nr

       if ( l_ncol(io) == 0 ) cycle
       
       ioff = orb_offset(io)
       
       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io) 
             
          mum%JCN(ind) = io - ioff
          mum%IRN(ind) = l_col(ind) - orb_offset(l_col(ind))

       end do
    end do

    call delete(tmpSp2)

  end subroutine prep_LHS

  subroutine allocate_mum(mum,size,N_Elec,Elecs,GF)
    use m_ts_electype
    include 'zmumps_struc.h'
    type(zMUMPS_STRUC), intent(inout) :: mum
    integer, intent(in) :: size, N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    complex(dp), pointer :: Gf(:)

    integer :: no, iEl, io

    ! We allocate for the equilibrium GF
    mum%NZ_RHS = size ! number of non-zero RHS elements
    if ( associated(mum%IRHS_PTR) ) then
       deallocate(mum%IRHS_PTR)
       deallocate(mum%IRHS_SPARSE)
       nullify(mum%IRHS_PTR,mum%IRHS_SPARSE,mum%RHS_SPARSE)
    end if
    allocate( mum%IRHS_PTR(mum%NRHS+1) )
    allocate( mum%IRHS_SPARSE(mum%NZ_RHS) )

    no = 0
    do iEl = 1 , N_Elec
       io = TotUsedOrbs(Elecs(iEl))
       no = no + io ** 2
    end do
    ! Allocate maximum space available
    if ( associated(Gf) ) then
       deallocate(Gf)
       nullify(Gf)
    end if
    ! for large electrodes and not so large device
    ! we need to allocate more space
    allocate( Gf(max(no,mum%NZ_RHS)) )
    mum%RHS_SPARSE => Gf(1:mum%NZ_RHS)

    no = 0
    do iEl = 1 , N_Elec

       ! This seems stupid, however, we never use the Sigma and
       ! GF at the same time. Hence it will be safe
       ! to have them point to the same array.
       ! When the UC_expansion_Sigma_GammaT is called
       ! first the Sigma is assigned and then 
       ! it is required that prepare_GF_inv is called
       ! immediately (which it is)
       ! Hence the GF must NOT be used in between these two calls!
       io = TotUsedOrbs(Elecs(iEl))
       Elecs(iEl)%Sigma => Gf(no+1:no+io**2)
       no = no + io ** 2

    end do

  end subroutine allocate_mum

  subroutine prep_RHS_Eq(mum,no_u_TS,nzs,N_Elec,Elecs,GF)
    use m_ts_electype
    use m_ts_sparse, only : tsup_sp_uc
    use m_ts_method, only : orb_offset, orb_type, TYP_BUFFER
    use class_Sparsity
    include 'zmumps_struc.h'
    type(zMUMPS_STRUC), intent(inout) :: mum
    integer, intent(in) :: no_u_TS, nzs, N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    complex(dp), pointer :: Gf(:)
    integer :: iEl, no, io, j, nr, ind
    integer, pointer :: l_ptr(:), l_ncol(:), l_col(:)

    ! Create the index 
    call attach(tsup_sp_uc,list_ptr=l_ptr, &
         n_col=l_ncol,list_col=l_col,nrows_g=nr)

    call allocate_mum(mum,nzs,N_Elec,Elecs,GF)

    ! TODO, this requires that the sparsity pattern is symmetric
    ! Which it always is!
    mum%IRHS_PTR(:) = 1 

    do io = 1 , nr
       if ( orb_type(io) == TYP_BUFFER ) cycle
       mum%IRHS_PTR(io-orb_offset(io)) = l_ptr(io) + 1
       if ( l_ncol(io) == 0 ) cycle ! no entries
       ! Create the row-index
       do j = 1 , l_ncol(io)
          ind = l_ptr(io) + j
          mum%IRHS_SPARSE(ind) = l_col(ind) - orb_offset(l_col(ind))
       end do

    end do
    mum%IRHS_PTR(no_u_TS+1) = nzs + 1

    if ( ind /= mum%NZ_RHS ) then
       write(*,*)ind,mum%NZ_RHS
       call die('Error in sparsity pattern of equilibrium')
    end if

  end subroutine prep_RHS_Eq

  subroutine prep_RHS_nEq(mum,no_u_TS,N_Elec, Elecs,Gf)
    use m_ts_electype
    use m_ts_sparse, only : tsup_sp_uc
    use m_ts_method, only : ts2s_orb
    use class_Sparsity
    include 'zmumps_struc.h'
    type(zMUMPS_STRUC), intent(inout) :: mum

    integer, intent(in) :: no_u_TS, N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    complex(dp), pointer :: Gf(:)
    integer :: iEl, no, i, io, jo, nr
    integer, pointer :: l_ptr(:), l_ncol(:), l_col(:)
    integer :: ind, j

    ! We only need a partial size of the Green's function
    no = sum(TotUsedOrbs(Elecs))

    call allocate_mum(mum,no*no_u_TS,N_Elec,Elecs,GF)

    ! Create the index 
    call attach(tsup_sp_uc,list_ptr=l_ptr, &
         n_col=l_ncol,list_col=l_col,nrows_g=nr)

    ! TODO, this requires that the sparsity pattern is symmetric
    ! Which it always is!
    mum%IRHS_PTR(1) = 1

    ind = 0
    do i = 1 , no_u_TS
       if ( i > 1 ) then
          ! initialize
          mum%IRHS_PTR(i) = mum%IRHS_PTR(i-1)
       end if
       ! get correct siesta-orbital
       io = ts2s_orb(i)
       iElec: do iEl = 1 , N_Elec
          if ( .not. OrbInElec(Elecs(iEl),io) ) cycle
          ! update pointer
          mum%IRHS_PTR(i) = mum%IRHS_PTR(i-1) + no_u_TS - 1
          ! Create the row-index
          do j = 1 , no_u_TS
             ind = ind + 1
             mum%IRHS_SPARSE(ind) = j
          end do
          exit iElec
       end do iElec
    end do
    mum%IRHS_PTR(no_u_TS+1) = no*no_u_TS + 1

    if ( ind /= mum%NZ_RHS ) then
       call die('Error in sparsity pattern of non-equilibrium')
    end if

  end subroutine prep_RHS_nEq

  subroutine insert_Self_Energies(mum, El)
    use m_ts_electype
    use m_ts_method, only : orb_offset, ts2s_orb
    include 'zmumps_struc.h'
    type(zMUMPS_STRUC), intent(inout) :: mum
    type(Elec), intent(in) :: El
    complex(dp), pointer :: iG(:)

    integer :: off, ii, no, ind, iso, jso

    iG => mum%A(:)

    no = TotUsedOrbs(El)
    off = El%idx_no - 1

    if ( El%Bulk ) then
       do ind = 1 , mum%NZ
          jso = ts2s_orb(mum%JCN(ind))
          if ( .not. OrbInElec(El,jso) ) cycle
          iso = ts2s_orb(mum%IRN(ind))
          if ( .not. OrbInElec(El,iso) ) cycle
          ii = (jso - El%idx_no ) * no + iso - off
          iG(ind) = El%Sigma(ii)
       end do
    else
       do ind = 1 , mum%NZ
          jso = ts2s_orb(mum%JCN(ind))
          if ( .not. OrbInElec(El,jso) ) cycle
          iso = ts2s_orb(mum%IRN(ind))
          if ( .not. OrbInElec(El,iso) ) cycle
          ii = (jso - El%idx_no ) * no + iso - off
          iG(ind) = iG(ind) - El%Sigma(ii)
       end do
    end if

  end subroutine insert_Self_Energies

#endif
end module m_ts_mumps_init
