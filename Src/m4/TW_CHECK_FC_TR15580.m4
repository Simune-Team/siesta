dnl @synopsis TW_CHECK_FC_TR15580([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler is fully
dnl compliant with the Fortran 95 Floating Point Exception Handling
dnl Extension, ISO TR15580.
dnl
dnl @version 1.0
dnl @author <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_TR15580],[
dnl
AC_MSG_CHECKING([$FC for compliance to the Floating Point Exception Handling Extension])
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
      Program test_tr15580

      Use, Intrinsic :: IEEE_Arithmetic
      Use, Intrinsic :: IEEE_Exceptions
      Use, Intrinsic :: IEEE_Features

      End Program test_tr15580
   ]]),   
   [AC_MSG_RESULT([yes])
    m4_ifval([$1],[$1],[])
   ],
   [AC_MSG_RESULT([no])
    m4_ifval([$2],[$2],
                  [AC_MSG_ERROR([A fully TR15580-compliant compiler is required.])])
   ]
)
AC_LANG_POP(Fortran)
dnl
])
