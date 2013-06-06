module m_comm
!
! Simple structure to keep track of source, target, starting points, and
! number of items to communicate
!
! This is a stand-alone module to simplify dependencies. Hopefully it
! can be made abstract enough to serve for multiple patterns of communication.
!
type, public :: comm_t
   integer :: src, dst, i1, i2, nitems
end type comm_t

end module m_comm
