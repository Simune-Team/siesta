      Module write_subs_positions

      use m_steps

      private
      public :: siesta_write_positions
      CONTAINS
      
      subroutine siesta_write_positions()
      use siesta_geom
      use atomlist, only: elem, iza
      use siesta_cml
      use units, only: Ang
      use zmatrix, only: lUseZmatrix, write_canonical_ucell_and_Zmatrix
      use m_iostruct, only: write_struct
      use siesta_options, only: idyn

      implicit none

!
!      Write out structural information in "crystallography" format,
!      and in "canonical" (see Zmatrix module) zmatrix format.
!      By now the structure might have changed in response to forces
!      and stresses.

        call write_struct( ucell, na_u, isa, iza, xa, moved=.true.)
        if (lUseZmatrix .and. (idyn .eq. 0)) then
           call write_canonical_ucell_and_Zmatrix
     $                         (filename="NEXT_ITER.UCELL.ZMATRIX")
        endif

        if (cml_p) then
          call cmlAddMolecule(xf=mainXML, natoms=na_u, elements=elem,
     .           atomRefs=cisa, coords=xa/Ang)
          call cmlAddLattice(xf=mainXML, cell=ucell/Ang, 
     .           units='siestaUnits:Ang', dictref='siesta:ucell')
        endif
      end subroutine siesta_write_positions

      END MODULE write_subs_positions
