dnl @synopsis TW_CHECK_FC_TR15581([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler is fully
dnl compliant with the Fortran 95 Enhanced Datatype Facilities
dnl Extension, ISO TR15581.
dnl
dnl @version 1.0
dnl @author <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_TR15581],[
dnl
AC_MSG_CHECKING([$FC for compliance to the Enhanced Datatype Facilities Extension])
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
      Program test_tr15581

      Type test_type
        Integer, Allocatable :: array(:)
      End Type test_type

      End Program test_tr15581

      Function test_function
        Integer, Allocatable :: test_function(:)

        Allocate(test_function(5))
      End Function test_function

      Subroutine test_subroutine(array)
        Integer, Allocatable :: array(:)

        Allocate(array(5))
      End Subroutine test_subroutine
   ]]),   
   [AC_MSG_RESULT([yes])
    m4_ifval([$1],[$1],[])
   ],
   [AC_MSG_RESULT([no])
    m4_ifval([$2],[$2],
                  [AC_MSG_ERROR([A fully TR15581-compliant compiler is required.])])
   ]
)
AC_LANG_POP(Fortran)
dnl
])
