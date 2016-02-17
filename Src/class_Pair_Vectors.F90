! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module class_Pair_Vectors

use class_Vector

character(len=*), parameter :: mod_name="Pair_Vectors"

#define PAIR_NAME Pair_Vectors
#define _T1_ Vector
#define _T2_ Vector

#include "Pair.T90"

end module class_Pair_Vectors


