dnl @synopsis TW_CHECK_FC_90([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler is fully
dnl compliant with Fortran 90 (ISO/IEC-1539:1991)
dnl If so, ACTION_IF_TRUE is performed; if not, ACTION_IF_FALSE
dnl 
dnl It currently tests for:
dnl
dnl Modules
dnl Private 
dnl New-style variable declarations.
dnl
dnl @version 1.0
dnl @author Toby White <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_90],[
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_MSG_CHECKING([for Fortran 90 compliance])
dnl
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
      Module test_module

      Implicit None
      Private

      Contains

      Function test_function() Result(out)
      Integer :: out
      out = 0
      End Function test_function

      End Module test_module
   ]]),
   [AC_MSG_RESULT([yes])
    m4_default([$1],[])
   ],
   [AC_MSG_RESULT([no])
    m4_default([$2],
               [AC_MSG_ERROR([ A fully Fortran-90-compliant compiler is required.])])
   ]
)
AC_LANG_POP(Fortran)
dnl
])
