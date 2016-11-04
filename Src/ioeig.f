! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine ioeig(eo, ef, no, ns, nk, maxo, nspinS, maxk,
     .                 kpoints, kweights)

c *******************************************************************
c Writes eigenvalues of Hamiltonian in k-points of sampling
c Emilio Artacho, Feb. 1999

      use precision, only : dp
      use siesta_cml
      use units, only : eV
      use files, only : slabel, label_length
      use m_spin, only : SPpol

      implicit          none

      integer,  intent(in) :: maxo
      integer,  intent(in) :: nspinS
      integer,  intent(in) :: maxk
      real(dp), intent(in) :: eo(maxo, nspinS, maxk)
      real(dp), intent(in) :: ef
      integer,  intent(in) :: no
      integer,  intent(in) :: ns
      integer,  intent(in) :: nk
      real(dp), intent(in) :: kpoints(3,nk)
      real(dp), intent(in) :: kweights(nk)
      
      external          io_assign, io_close

c Internal 
      integer           ik, iu, io, is

      character(len=label_length+4) :: fname
c -------------------------------------------------------------------

      fname = slabel
      fname = trim(fname) // '.EIG'
      
      call io_assign( iu )
      open( iu, file=fname, form='formatted', status='unknown' )      

      write(iu,"(f14.4)") ef/eV
      write(iu,"(3i6)")   no, nspinS, nk
      do ik = 1,nk
        write(iu,"(i5,10f12.5,/,(5x,10f12.5))")
     .          ik, ((eo(io,is,ik)/eV,io=1,no),is=1,nspinS)
      enddo

      call io_close( iu )

      if (cml_p) then
         call cmlStartPropertyList(xf=mainXML, title="Eigenvalues")
         call cmlAddProperty(xf=mainXML, value=ef/eV, 
     .        title='Fermi Energy', dictref='siesta:E_Fermi', 
     .        fmt='r5', units='siestaUnits:ev')
         call cmlAddProperty(xf=mainXML, value=nk, 
     .        title='Number of k-points', dictRef='siesta:nkpoints',
     .        units='cmlUnits:countable')
         do is = 1 , nspinS
            call cmlStartPropertyList(mainXML,
     .           dictRef='siesta:kpt_band')
            if ( is == 1 .and. nspinS > 1 ) then
               call cmlAddProperty(xf=mainXML, value="up", 
     .              dictRef="siesta:spin")
            else if ( nspinS > 1 ) then
               call cmlAddProperty(xf=mainXML, value="down", 
     .              dictRef="siesta:spin")
            end if
            do ik = 1, nk
               call cmlAddKPoint(xf=mainXML, coords=kpoints(:, ik), 
     .              weight=kweights(ik))
               call cmlAddProperty(xf=mainXML, value=eo(1:no,is,ik)/eV, 
     .              dictRef='siesta:eigenenergies',
     .              units='siestaUnits:ev')
!     call cmlAddBand(xf=mainXML, 
!     .           kpoint=kpoints(:, ik), kweight=kweights(ik), 
!     .           bands=eo(1:no,is,ik))
            enddo
            call cmlEndPropertyList(mainXML)
         enddo
         call cmlEndPropertyList(mainXML)
      endif
      
      end subroutine ioeig
