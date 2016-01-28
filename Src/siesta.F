! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
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

C     Begin of coordinate relaxation iteration
      relaxd = .false.
      istep  = inicoor
      DO WHILE ((istep.le.fincoor) .AND. (.not. relaxd))
        call siesta_forces( istep )

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
