      MODULE MPI_INTERFACES
!
!
!       Michael Hennecke
!       A Fortran 90 interface to MPI version 1.1
!       RZ Uni Karlsruhe, Internal Report 63/96
!
!         03-Oct-1996, hennecke@rz.uni-karlsruhe.de (v0.9c beta)
!
!       Permission is granted to copy and distribute this file
!       or modified versions of this file for no fee, provided the 
!       copyright notice and this permission notice are preserved 
!       on all copies.
!
!       (C) 1996  Michael Hennecke, RZ Universitaet Karlsruhe
!
!
!       Modified by Alberto Garcia, wdpgaara@lg.ehu.es                       
!
! In case we are using this name as a compiler preprocessor symbol...
#ifdef MPI
#define __USING__MPI
#undef MPI
#endif
        USE MPI

!       ... generic overloads for <choice> argument routines ...

        USE MPI__integer_V    
        USE MPI__integer_S
        USE MPI__integer8_V    
        USE MPI__integer8_S
        USE MPI__logical_V   
        USE MPI__logical_S
        USE MPI__character_V  
        USE MPI__character_S

        include "V_S.uses"

!       ... these are for two different-rank <choice> arguments ...

        USE MPI__integer_VS  
        USE MPI__integer_SV
        USE MPI__integer8_VS  
        USE MPI__integer8_SV
        USE MPI__logical_VS   
        USE MPI__logical_SV
        USE MPI__character_VS 
        USE MPI__character_SV

        include "VS.uses"

!       ... this is for *ALL* combinations of type/rank of SENDRECV ...
!       USE MPI__sendrecv

        IMPLICIT NONE
!
        public

      END MODULE MPI_INTERFACES
