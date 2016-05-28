! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module class_Fstack_Pair_Vectors

use class_Pair_Vectors

character(len=*), parameter :: mod_name="Fstack_Pair_Vectors"

#define _T_ Pair_Vectors
#define FSTACK_NAME Fstack_Pair_Vectors

#include "Fstack.T90"

end module class_Fstack_Pair_Vectors
