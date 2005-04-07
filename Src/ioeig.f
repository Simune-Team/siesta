      subroutine ioeig(eo, ef, no, ns, nk, maxo, maxs, maxk,
     .                 kpoints, kweights)

c *******************************************************************
c Writes eigenvalues of Hamiltonian in k-points of sampling
c Emilio Artacho, Feb. 1999

      use fdf
      use precision, only : dp
      use siesta_cml
      use units, only : eV


      implicit          none

      integer,  intent(in) :: maxo
      integer,  intent(in) :: maxs
      integer,  intent(in) :: maxk
      real(dp), intent(in) :: eo(maxo, maxs, maxk)
      real(dp), intent(in) :: ef
      integer,  intent(in) :: no
      integer,  intent(in) :: ns
      integer,  intent(in) :: nk
      real(dp), intent(in) :: kpoints(3,nk)
      real(dp), intent(in) :: kweights(nk)
      
      external          io_assign, io_close, paste

c Internal 
      integer           ik, iu, io, is, nspin

      character(len=33), save :: fname
      logical, save           :: frstme = .true.
c -------------------------------------------------------------------

      if (frstme) then
        fname = fdf_string( 'SystemLabel', 'siesta' )
        fname = trim(fname) // '.EIG'
        frstme = .false.
      endif
      
      nspin=min(ns,2)

      call io_assign( iu )
      open( iu, file=fname, form='formatted', status='unknown' )      

      write(iu,"(f14.4)") ef/eV
      write(iu,"(3i6)")   no, min(ns,2), nk
      do ik = 1,nk
        write(iu,"(i5,10f12.5,/,(5x,10f12.5))")
     .          ik, ((eo(io,is,ik)/eV,io=1,no),is=1,nspin)
      enddo

      call io_close( iu )

      if (cml_p) then
        call cmlStartPropertyList(xf=mainXML, title="Eigenvalues")
        call cmlAddProperty(xf=mainXML, property=ef/eV, 
     .       title='Fermi Energy', dictref='siesta:E_Fermi', 
     .       fmt='(f12.5)', units='eV')
        call cmlAddProperty(xf=mainXML, property=nk, 
     .       title='Number of k-points', dictRef='siesta:nkpoints')
        call cmlEndPropertyList(mainXML)
        do is = 1,nspin
          call cmlStartBandList(mainXML)
          if (nspin.eq.2) then
            if (is.eq.1) then
              call xml_addAttribute(xf=mainXML, 
     .             name='Spin', value='Up')
            else
              call xml_addAttribute(xf=mainXML, 
     .             name='Spin', value='Down')
            endif
          endif
          do ik = 1, nk
            call cmlAddBand(xf=mainXML, 
     .           kpoint=kpoints(:, ik), kweight=kweights(ik), 
     .           bands=eo(1:no,is,ik))
          enddo
          call cmlEndBandList(mainXML)
        enddo
      endif

      end subroutine ioeig
