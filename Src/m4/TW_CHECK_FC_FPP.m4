dnl @synopsis TW_CHECK_FC_FPP([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler supports
dnl cpp-like functionality when called on a suitable fixed-format file.
dnl If so, ACTION_IF_TRUE is performed; if not, ACTION_IF_FALSE
dnl 
dnl @version 1.0
dnl @author Toby White <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_FPP],[
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_FC_SRCEXT(F)
dnl
AC_MSG_CHECKING([whether $FC has an integrated Fortran cpp-style preprocessor for fixed-form source])
dnl
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
      Program test_cpp
#if 1
      Integer i
#else
      Integer j
#endif
      End Program
   ]]),
   [AC_MSG_RESULT([yes])
    m4_ifval([$1],[$1],[])
   ],
   [AC_MSG_RESULT([no])
    m4_ifval([$2],[$2],
                  [AC_MSG_ERROR([A Fortran compiler with integrated cpp-style preprocessor for fixed-form source is required.])])
   ]
)
AC_LANG_POP(Fortran)
dnl
])
