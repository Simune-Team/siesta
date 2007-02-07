      MODULE m_local_DOS
      private
      public :: local_DOS

      CONTAINS

      subroutine local_DOS( )
      USE SIESTA_TODO

      implicit none

      integer :: dummy_iscf = 1

! Find local density of states
      if (IOnode) then
        genlogic = fdf_block('LocalDensityOfStates',iu)
      endif
      call broadcast(genlogic)

      if ( genlogic ) then

! Find the desired energy range
        if (IOnode) then
          read(iu,'(a)') line
          p=>digest(line)
          if (.not. match(p,"vvn"))
     .       call die("Wrong format in LocalDensityofStates")
          factor = fdf_convfac( names(p,1), 'Ry' )
          e1 = values(p,1)*factor
          e2 = values(p,2)*factor
          call destroy(p)
        endif
        call broadcast(e1)
        call broadcast(e2)

!       Find the density matrix for states between e1 and e2
        if (isolve .eq. 0) then
          call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     .                numh, listhptr, listh, numh, listhptr, listh, 
     .                H, S, qtot, fixspin, qs, temp, e1, e2,
     .                gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     .                eo, qo, Dscf, Escf, ef, efs, Entrop, no_u,
     .                occtol, dummy_iscf, neigwanted)

!       Find the LDOS in the real space mesh
          filrho = paste( slabel, '.LDOS' )
          g2max = g2cut
          call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .              no_u, na_u, na_s, isa, xa, indxua, 
     .              ucell, mscell, g2max, ntm,
     .              0, 0, 0, filrho, ' ', ' ', ' ', ' ', ' ',
     .              maxnh, numh, listhptr, listh, Dscf, Datm,
     .              maxnh, numh, listhptr, listh, H,
     .              Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, Exc, Dxc,
     .              dipol, fa, stress, fal, stressl )
        else
          if (IOnode)  write(6,*)
     .       'siesta: ERROR: LDOS implemented only with diagon'
        endif

      endif ! genlogic


      END subroutine local_DOS

      END module m_local_DOS
