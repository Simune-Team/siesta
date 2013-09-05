      module m_ncps_reader

      use m_ncps_froyen_ps_t,    only: pseudopotential_t => froyen_ps_t

      public :: pseudo_read

      CONTAINS

        subroutine pseudo_read(label,p)
        use m_ncps_froyen_reader,  only: pseudo_read_formatted
        use m_ncps_froyen_reader,  only: pseudo_read_unformatted

        character(len=*), intent(in)   :: label
        type(pseudopotential_t)                    :: p

!       PS information can be in a .vps file (unformatted)
!       or in a .psf file (formatted)
!       or in a .xml file 

        character(len=40) fname
        logical found

        fname  = trim(label) // '.vps'
        inquire(file=fname, exist=found)
        if (found) then
           call pseudo_read_unformatted(fname,p)
        else
           fname = trim(label) // '.psf'
           inquire(file=fname, exist=found)
           if (found) then
              call pseudo_read_formatted(fname,p)
           else
              fname = trim(label) // '.xml'
              inquire(file=fname, exist=found)
              if (found) then
                 call pseudo_read_xml(fname,p)
              else
                 write(6,'(/,2a,a20,/)') 'read_pseudo: ERROR: ',
     .                'Pseudopotential file not found: ', fname
                 call die
              endif
           endif
        endif
!        if (write_ion_plot_files)
!     $       call pseudo_dump(trim(label) // ".psdump",p)
        end subroutine pseudo_read
!
        subroutine pseudo_read_xml(fname,p)

        use m_ncps_xml_ps_t, only: xml_ps_t
        use m_ncps_xmlreader, only: ncps_xmlreader
        use m_ncps_translators, only: ncps_xml2froyen

        character(len=*), intent(in)              :: fname
        type(pseudopotential_t), intent(out)      :: p

        type(xml_ps_t)      :: psxml

        call ncps_xmlreader(fname,psxml)
        call ncps_xml2froyen(psxml,p)

        end subroutine pseudo_read_xml
!----
      end module m_ncps_reader
