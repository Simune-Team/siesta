      Program SIESTA

      use m_siesta_init
      use m_siesta_analysis
      use m_siesta_move
      use m_siesta_end
      use m_siesta_forces
      USE m_steps, only: inicoor, fincoor
      use alloc
      use siesta_options, only: td_elec_dyn
      use m_siesta_tddft
      implicit none

      integer :: istep
      logical :: relaxd
      logical :: tdstart

!------------------------------------------------------------------------- BEGIN
      call siesta_init()

C     Begin of coordinate relaxation iteration
      relaxd = .false.
      istep  = inicoor
      tdstart = .true.
      DO WHILE ((istep.le.fincoor) .AND. (.not. relaxd))
        
        if(td_elec_dyn) then
        call siesta_tddft( istep, tdstart )
        else

        call siesta_forces( istep )
        end if 

        call siesta_move( istep, relaxd )

        if (.not. relaxd) then
          istep = istep + 1
        endif

      ENDDO

C     End of coordinate-relaxation loop 
      call siesta_analysis( relaxd )

      call siesta_end()
!--------------------------------------------------------------------------- END
      END program siesta
