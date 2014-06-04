program test_psml

use m_psml, only: psml_t, psml_destroy, psml_reader   ! Clarify this
use m_semicore_info, only: get_n_semicore_shells

type(psml_t), pointer :: psxml
integer      :: nsemic(0:3)

call psml_reader("PSML",psxml)
call get_n_semicore_shells(psxml,nsemic)

print "(a,4i4)", "Number of semicore shells:", nsemic(0:3)
end program test_psml

