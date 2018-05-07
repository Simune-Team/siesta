! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine ioeig(eo, ef, no, ns, nk, maxo, maxs, maxk,
     .                 kpoints, kweights)

c *******************************************************************
c Writes eigenvalues of Hamiltonian in k-points of sampling
c Emilio Artacho, Feb. 1999

      use precision, only : dp
      use siesta_cml
      use units, only : eV
      use files, only : slabel, label_length

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
      
      external          io_assign, io_close

c Internal 
      integer           ik, iu, io, is, nspin

      character(len=label_length+4) :: fname
c -------------------------------------------------------------------

      fname = trim(slabel) // '.EIG'
      
      nspin = min(ns,2)

      call io_assign( iu )
      open( iu, file=fname, form='formatted', status='unknown' )      

      write(iu,"(e17.9)") ef/eV
      if ( ns > 2 ) then
        write(iu,"(tr1,i10,i2,tr1,i10)")   no*2, 1, nk
      else
        write(iu,"(tr1,i10,i2,tr1,i10)")   no, min(ns,2), nk
      end if
      do ik = 1,nk
        write(iu,"(i10,10(tr1,e17.9),/,(tr10,10(tr1,e17.9)))")
     .          ik, ((eo(io,is,ik)/eV,io=1,no),is=1,nspin)
      enddo

      call io_close( iu )

      if (cml_p) then
        call cmlStartPropertyList(xf=mainXML, title="Eigenvalues")
        call cmlAddProperty(xf=mainXML, value=ef/eV, 
     .       title='Fermi Energy', dictref='siesta:E_Fermi', 
     .       fmt='r5', units='siestaUnits:ev')
        call cmlAddProperty(xf=mainXML, value=nk, 
     .       title='Number of k-points', dictRef='siesta:nkpoints',
     .       units='cmlUnits:countable')
        if ( ns > 2 ) then
           call cmlStartPropertyList(mainXML, dictRef='siesta:kpt_band')
           do ik = 1, nk
              call cmlAddKPoint(xf=mainXML, coords=kpoints(:, ik), 
     .             weight=kweights(ik))
              call cmlAddProperty(xf=mainXML,
     .             value=reshape(eo(1:no,1:2,ik)/eV, (/no*2/)),
     .             dictRef='siesta:eigenenergies',
     .             units='siestaUnits:ev')
           end do
           call cmlEndPropertyList(mainXML)
        else
         do is = 1,nspin
           call cmlStartPropertyList(mainXML, dictRef='siesta:kpt_band')
           if (nspin.eq.2) then
            if (is.eq.1) then
              call cmlAddProperty(xf=mainXML, value="up", 
     .                            dictRef="siesta:spin")
            else
              call cmlAddProperty(xf=mainXML, value="down", 
     .                            dictRef="siesta:spin")
            endif
           endif
           do ik = 1, nk
              call cmlAddKPoint(xf=mainXML, coords=kpoints(:, ik), 
     .                          weight=kweights(ik))
              call cmlAddProperty(xf=mainXML, value=eo(1:no,is,ik)/eV, 
     .                          dictRef='siesta:eigenenergies',
     .                          units='siestaUnits:ev')
           end do
           call cmlEndPropertyList(mainXML)
         end do
        end if
        call cmlEndPropertyList(mainXML)
      endif

      end subroutine ioeig
