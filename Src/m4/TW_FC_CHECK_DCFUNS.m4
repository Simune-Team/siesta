dnl Macro to check how to use functions that NAG requires -dcfuns for.
dnl DIMAG used.

AC_DEFUN([_TW_DIMAG_PROGRAM],
[
      PROGRAM TESTDIMAG
      DOUBLE PRECISION A
      DOUBLE COMPLEX B
      A=DIMAG(B)
      END PROGRAM
])

AC_DEFUN([TW_FC_CHECK_DCFUNS], [
AC_REQUIRE([AC_PROG_FC])
AC_LANG_ASSERT(Fortran)
tw_dcfuns_ok=no


AC_MSG_CHECKING([how to compile DIMAG])

AC_LINK_IFELSE([_TW_DIMAG_PROGRAM],
               [tw_dcfuns_ok=yes;
                AC_MSG_RESULT([default])],[])

if test $tw_dcfuns_ok = no; then
   FCFLAGS_save=$FCFLAGS
   FCFLAGS="$FCFLAGS -dcfuns"
   AC_LINK_IFELSE([_TW_DIMAG_PROGRAM],
                  [tw_dcfuns_ok=yes;
                   AC_MSG_RESULT([with -dcfuns])],[]) 
   if test $tw_dcfuns_ok = no; then
      FCFLAGS=$FCFLAGS_save
   fi
fi

AS_IF([test $tw_dcfuns_ok = yes],
      [$1],
      [m4_default([$2],[AC_MSG_ERROR([Cannot compile DIMAG function])])])
])# TW_FC_CHECK_DCFUNS
