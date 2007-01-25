      subroutine siesta_end()
      USE SIESTA_TODO
      implicit none

      external :: automatic_cell,
     .  bonds, cgvc, cgvc_zmatrix, fixed,
     .  dhscf, diagon, dnaefs, extrapol, initatom,
     .  iodm, iozm,
     .  kgrid, kgridinit, kinefsm, mulliken, naefs, neighb,
     .  pulayx, 
     .  reinit, shaper, spnvec, 
     .  timer, xijorb, memory,
     .  ioeig, iofa, iokp, iomd, prversion, eggbox


!------------------------------------------------------------------------- BEGIN
! Output memory use up to the end of the program
      call printmemory( 6, 1 )

! Print allocation report
      call alloc_report( printNow=.true. )

! Stop time counter
      call timer( 'siesta', 2 )
      call timer( 'all', 3 )

! Print final date and time
      if (IOnode) then
        call timestamp('End of run')
        call wallclock('End of run')
      endif

! Finalize MPI
#ifdef MPI
      call MPI_Finalize( MPIerror )
#endif

      if (cml_p) then
        call cmlEndModule(mainXML)
        call siesta_cml_exit()
      endif      
!--------------------------------------------------------------------------- END
      END subroutine siesta_end
