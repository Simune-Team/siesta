module m_psml

!
  use m_psml_reader, only: psml_reader 

  use m_psml_core, only: ps_t, ps_destroy

  use m_psml_core, only: str_of_set
  use m_psml_core, only:     SET_SREL     , &
                             SET_NONREL   , &
                             SET_SO       , &
                             SET_LJ       , &
                             SET_UP       , &
                             SET_DOWN     , &
                             SET_SPINAVE  , &
                             SET_SPINDIFF , &
                             SET_USER1    , &
                             SET_USER2    , &
                             SET_ALL

  use m_psml_api

  use assoc_list, only: ps_annotation_t => assoc_list_t
  use assoc_list, only: nitems_annotation => assoc_list_nitems
  use assoc_list, only: get_annotation_key => assoc_list_get_key
  use assoc_list, only: get_annotation_value => assoc_list_get_value

  public

end module m_psml
