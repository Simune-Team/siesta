      Program SIESTA

      use m_siesta_init
      use m_siesta_analysis
      use m_siesta_move
      use m_siesta_end
      use m_siesta_forces
      USE m_steps, only: inicoor, fincoor
      use alloc

      implicit none

      integer :: istep
      logical :: relaxd

!------------------------------------------------------------------------- BEGIN
      call siesta_init()
#ifdef DEBUG
        call write_debug( 'POS siesta_init' )
#endif

C     Begin of coordinate relaxation iteration
      relaxd = .false.
      istep  = inicoor
      DO WHILE ((istep.le.fincoor) .AND. (.not. relaxd))
#ifdef DEBUG
        call write_debug( 'PRE siesta_forces' )
#endif
        call siesta_forces( istep )

#ifdef DEBUG
        call write_debug( 'PRE siesta_move' )
#endif
        call siesta_move( istep, relaxd )
        if (.not. relaxd) then
          istep = istep + 1
        endif

      ENDDO

C     End of coordinate-relaxation loop 
#ifdef DEBUG
      call write_debug( 'PRE siesta_analysis' )
#endif
      call siesta_analysis( relaxd )

#ifdef DEBUG
      call write_debug( 'SIESTA END' )
#endif
      call siesta_end()
!--------------------------------------------------------------------------- END
      END program siesta
