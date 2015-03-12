module local_xml

  ! Simple module to share the
  ! XML file handle

  use flib_wxml, only: xmlf_t

  type(xmlf_t), public            :: xf
end module local_xml
