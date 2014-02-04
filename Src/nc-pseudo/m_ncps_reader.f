      module m_ncps_reader

      use m_ncps_froyen_ps_t,    only: pseudopotential_t => froyen_ps_t

      integer, parameter, private :: dp = selected_real_kind(14,100)

      public :: pseudo_read

      CONTAINS

        subroutine pseudo_read(label,p,new_grid,a,b,rmax)
        use m_ncps_froyen_reader,  only: pseudo_read_formatted
        use m_ncps_froyen_reader,  only: pseudo_read_unformatted
        use m_ncps_froyen_reader,  only: pseudo_reparametrize

        character(len=*), intent(in)   :: label
        type(pseudopotential_t)        :: p
        logical, intent(in), optional  :: new_grid
        real(dp), intent(in), optional :: a
        real(dp), intent(in), optional :: b
        real(dp), intent(in), optional :: rmax

!       PS information can be in a .vps file (unformatted)
!       or in a .psf file (formatted)
!       or in a .xml file 

        character(len=40) fname
        logical found, reparametrize

        reparametrize = .false.
        if (present(new_grid)) then
           reparametrize = new_grid
        endif
        if (reparametrize) then
           if (.not. present(a)) call die("New a not present")
           if (.not. present(b)) call die("New b not present")
        endif

        fname  = trim(label) // '.vps'
        inquire(file=fname, exist=found)
        if (found) then
           call pseudo_read_unformatted(fname,p)
           if (reparametrize) then
              call pseudo_reparametrize(p,a,b,label,rmax)
           endif
        else
           fname = trim(label) // '.psf'
           inquire(file=fname, exist=found)
           if (found) then
              call pseudo_read_formatted(fname,p)
              if (reparametrize) then
                 call pseudo_reparametrize(p,a,b,label,rmax)
              endif
           else
              fname = trim(label) // '.xml'
              inquire(file=fname, exist=found)
              if (found) then
                 call pseudo_read_xml(fname,p,reparametrize,a,b,rmax)
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
        subroutine pseudo_read_xml(fname,p,reparametrize,a,b,rmax)

        use m_ncps_xml_ps_t, only: xml_ps_t, xml_ps_destroy
        use m_ncps_xmlreader, only: ncps_xmlreader
        use m_ncps_translators, only: ncps_xml2froyen_new

        character(len=*), intent(in)              :: fname
        type(pseudopotential_t), intent(out)      :: p
        logical, intent(in), optional  :: reparametrize
        real(dp), intent(in), optional :: a
        real(dp), intent(in), optional :: b
        real(dp), intent(in), optional :: rmax

        type(xml_ps_t), pointer        :: psxml=>null()

        call ncps_xmlreader(fname,psxml)
        call ncps_xml2froyen_new(psxml,p,reparametrize,a,b,rmax)
        call xml_ps_destroy(psxml)

        end subroutine pseudo_read_xml
!----
      end module m_ncps_reader
