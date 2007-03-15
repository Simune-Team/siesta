module m_neighbors

  use precision, only: dp

  implicit none

  private

! First stage: Hold the variables

  ! Max. number of neighbor atoms
  integer, public                  :: maxna=200 

  integer,  pointer, save, public  :: jna(:)
  real(dp), pointer, save, public  :: r2ij(:)
  real(dp), pointer, save, public  :: xij(:,:)

! Second stage: Create a "neighboring atoms" sparse data structure

  public :: init_neighbor_arrays

CONTAINS

  subroutine init_neighbor_arrays(rmax)
    ! Resizes neighbor arrays so that a range "rmax" is covered

    use siesta_geom,    only: scell, xa, na_u, na_s
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

    ! Initialize neighb subroutine  (ia=0)
    nnia = maxna
    call neighb( scell, rmax, na_s, xa, 0, isel, nnia, jna, xij, r2ij )
    nnamax = 0
    do ia = 1,na_u
       nnia = 0  ! Pass without filling arrays,
                 ! just to get maximum size needed
       call neighb( scell, rmax, na_s, xa, ia, isel, nnia, jna, xij, r2ij )
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
