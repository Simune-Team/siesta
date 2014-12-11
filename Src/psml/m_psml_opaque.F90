module m_psml_opaque

  use m_psml_core
  use iso_c_binding

   type, public, bind(c) :: psml_t
      private
      integer(c_int)   :: i = 0
   end type psml_t

   type, private :: ps_store_t
      integer                         :: nitems = 0
      integer                         :: nslots = 0
      type(ps_t),         allocatable :: p(:)
   end type ps_store_t

   public :: psml_reader
   public :: psml_AtomicLabel
   public :: psml_NValenceShells

   type(ps_store_t), target :: ps_store

   private

   CONTAINS

     function psml_AtomicLabel(psml) result(name)

       type(psml_t), intent(in) :: psml
       character(len=20) :: name

       name = trim(ps_AtomicLabel(psml2ps(psml)))

     end function psml_AtomicLabel
!
function psml_NValenceShells(psml) result(nshells)
type(psml_t), intent(in) :: psml
integer                :: nshells

nshells = ps_NValenceShells(psml2ps(psml))

end function psml_NValenceShells

function psml2ps(psml) result(ps)
type(psml_t), intent(in) :: psml
type(ps_t), pointer :: ps

       ps => ps_store%p(psml%i)
end function psml2ps

     subroutine psml_reader(fname,psml)
       use m_psml_reader, reader=>psml_reader

       character(len=*), intent(in) :: fname
       type(psml_t), intent(inout) :: psml

       type(ps_t), pointer :: ps

       call get_free_slot(psml,ps)

       call reader(fname,ps)

     end subroutine psml_reader

       subroutine get_free_slot(psml,ps)
        use m_psml_core, only: ps_destroy
        
         type(psml_t), intent(inout) :: psml
         type(ps_t), pointer :: ps

         if (psml%i /= 0) then
            ! handle has been used. Clean
            ps => ps_store%p(psml%i)
            call ps_destroy(ps)
            RETURN
         endif
         
         ! Get fresh storage
         if (.not. allocated(ps_store%p)) then
            allocate(ps_store%p(4))
            ps_store%nslots = 4
            ps_store%nitems = 0
         endif
       
         nitems = ps_store%nitems
         if (nitems == ps_store%nslots) then
            call die("ps_store is full...")
         endif

         ps => ps_store%p(nitems+1)
         psml%i = nitems + 1
         ps_store%nitems = nitems + 1
         

       end subroutine get_free_slot
end module m_psml_opaque

#ifdef __TEST__
program test_psml
use m_psml_opaque

implicit none

type(psml_t) :: psml(5)
character(len=8) :: label

integer :: i, nshells

do i = 1, 5
   print *, "iteration: ", i
   call psml_reader("PSML",psml(i))
   label = psml_AtomicLabel(psml(i))
   nshells = psml_NValenceShells(psml(i))
   print *, label, nshells
enddo
end program test_psml

subroutine die(str)
character(len=*), optional :: str
if (present(str)) then
   print "(a)", trim(str)
endif
stop
end subroutine die

#endif

     

   
   

