module m_psml

!
  use m_psml_reader, only: psml_reader 

  use m_psml_core

  use assoc_list, only: nitems_annotation => assoc_list_nitems
  use assoc_list, only: get_annotation_key => assoc_list_get_key
  use assoc_list, only: get_annotation_value => assoc_list_get_value

  public

end module m_psml
