module m_psml_user

  use m_psml_core, only: ps_t
  use iso_c_binding

   type, public, bind(c) :: psh_t
      private
      integer(c_int)   :: i = 0
   end type psh_t

   type, private :: ps_store_t
      integer                         :: nitems = 0
      integer                         :: nslots = 0
      type(ps_t),         allocatable :: p(:)
   end type ps_store_t

   public :: psml_reader
   public :: ps_AtomicLabel

   type(ps_store_t), target :: ps_store

   private

   CONTAINS

     function ps_AtomicLabel(psh) result(name)
!       use m_psml_core, only: atomic_label=>ps_AtomicLabel

       type(psh_t), intent(in) :: psh
       character(len=20) :: name

       type(ps_t), pointer :: ps

       ps => ps_store%p(psh%i)
       name = trim(ps%header%atomic_label)
     end function ps_AtomicLabel
!

     subroutine psml_reader(fname,psh)
       use m_psml_reader, reader=>psml_reader

       character(len=*), intent(in) :: fname
       type(psh_t), intent(inout) :: psh

       type(ps_t), pointer :: ps

       call get_free_slot(psh,ps)

       call reader(fname,ps)

     end subroutine psml_reader


       subroutine get_free_slot(psh,ps)
        use m_psml_core, only: ps_destroy
        
         type(psh_t), intent(inout) :: psh
         type(ps_t), pointer :: ps

         if (psh%i /= 0) then
            ! handle has been used. Clean
            ps => ps_store%p(psh%i)
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
         psh%i = nitems + 1
         ps_store%nitems = nitems + 1
         

       end subroutine get_free_slot
end module m_psml_user

#ifdef __TEST__
program test_psh
use m_psml_user

type(psh_t) :: psh(5)
character(len=8) :: label

do i = 1, 5
   print *, "iteration: ", i
   call psml_reader("PSML",psh(i))
   label = ps_AtomicLabel(psh(i))
   print *, label
enddo
end program test_psh

subroutine die(str)
character(len=*), optional :: str
if (present(str)) then
   print "(a)", trim(str)
endif
stop
end subroutine die

#endif

     

   
   

