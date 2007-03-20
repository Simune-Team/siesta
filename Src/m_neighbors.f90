module m_neighbors

  use precision, only: dp

  implicit none

  private

! Holds the arrays needed for calls to neighb
! Users of this module should always initialize neighb
! by calling it with ia=0
! (This is necessary since neighb's internal data structures
! depend on whether the unit cell or the auxiliary supercell
! are considered, even though the actual neighbors do not).

  ! Max. number of neighbor atoms
  integer, public                  :: maxna=200 

  integer,  pointer, save, public  :: jna(:)
  real(dp), pointer, save, public  :: r2ij(:)
  real(dp), pointer, save, public  :: xij(:,:)

  public :: init_neighbor_arrays

CONTAINS

  subroutine init_neighbor_arrays(rmax)

    ! Resizes neighbor arrays so that a range "rmax" is covered

    use siesta_geom,    only: ucell, xa, na_u
    use alloc,          only: re_alloc

    implicit none

    real(dp), intent(in) :: rmax

    integer :: isel, ia, nnia, nnamax

    logical, save :: first_time = .true.

    external :: neighb

!-------------------------------------------------------------- BEGIN
    if (first_time) then
       nullify(jna,xij,r2ij)
       call re_alloc(jna,1,maxna,name="jna",routine="init_neighbor_arrays")
       call re_alloc(xij,1,3,1,maxna,name="xij",routine="init_neighbor_arrays")
       call re_alloc(r2ij,1,maxna,name="r2ij",routine="init_neighbor_arrays")
       first_time = .false.
    endif

    isel = 0  

    ! It is enough to deal with the unit cell, as images are
    ! handled also.

    ! Initialize neighb subroutine  (ia=0)
    nnia = maxna
    call neighb( ucell, rmax, na_u, xa, 0, isel, nnia, jna, xij, r2ij )

    ! Find maximum number of neighbors
    nnamax = 0
    do ia = 1,na_u
       nnia = 0  ! Pass without filling arrays,
                 ! just to get maximum size needed
       call neighb( ucell, rmax, na_u, xa, ia, isel, nnia, jna, xij, r2ij )
       nnamax = max( nnamax, nnia )
    enddo
    if (nnamax .gt. maxna) then
       maxna = nnamax      ! Update maximum size and reallocate
       call re_alloc(jna,1,maxna,name="jna",routine="init_neighbor_arrays")
       call re_alloc(xij,1,3,1,maxna,name="xij",routine="init_neighbor_arrays")
       call re_alloc(r2ij,1,maxna,name="r2ij",routine="init_neighbor_arrays")
    endif
    
  end subroutine init_neighbor_arrays

end module m_neighbors
