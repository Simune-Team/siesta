dnl @synopsis TW_CHECK_FC_95([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler is fully
dnl compliant with Fortran 95 (ISO-IEC 1539-1:1997)
dnl 
dnl It currently tests for:
dnl
dnl Named End Interface
dnl Derived type initialization
dnl The Null() intrinsic
dnl The Forall statement 
dnl The Cpu_Time intrinsic
dnl Pure functions
dnl Elemental functions
dnl 
dnl @version 1.0
dnl @author <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_95],[
dnl
AC_MSG_CHECKING([for Fortran 95 compliance])
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
      Program test_f95

!      Interface test_interface
!      End Interface test_interface

      Type test_type
        Integer :: i = 1
      End Type test_type

      Integer, Pointer :: j => Null()

      Integer :: i
      Real :: a

      Forall (i=1:50)
      End Forall

      Call CPU_TIME(a)

      Contains

      Pure Integer Function test_pure()
        test_pure = 0
      End Function test_pure

      Elemental Integer Function test_elemental(in)
        Integer, Intent(In) :: in
        test_elemental = 0
      End Function test_elemental

      End Program test_f95
   ]]),
   [AC_MSG_RESULT([yes])
    m4_default([$1],[:])
   ],
   [AC_MSG_RESULT([no])
    m4_default([$2], 
               [AC_MSG_ERROR([A fully Fortran-95-compliant compiler is required.])])
   ]
)
AC_LANG_POP(Fortran)
dnl
])
