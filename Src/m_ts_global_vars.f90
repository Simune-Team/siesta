! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!==========================================================================*
!                                                                          *
!  TRANSIESTA MODULE m_ts_global_vars : Declaration of the TS variables    *
!  that are accessed in different parts of the code by using a:            *
!      use m_ts_global_vars                                                *
!  declaration, instead of passing as dummy arguments                      *  
!                                                                          *
!  Written by F.D.Novaes, Apr'10                                           *
!==========================================================================*                                         


module m_ts_global_vars

USE precision, only : dp

integer :: TSiscf=0

logical :: TSinit=.false.,TSrun=.false.

integer :: ts_istep      ! FC step in phonon calculation

end module m_ts_global_vars
