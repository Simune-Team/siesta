dnl Macro to check that correct flags have been chosen for 
dnl compilation with BLACS. Only works with Fortran at the moment.

AC_DEFUN([_TW_TRY_BLACS], [
ac_ext=f
AC_LINK_IFELSE(
   [AC_LANG_SOURCE([[
      PROGRAM BLACS                   
      INTEGER MYPNUM, NPROCS        
      CALL BLACS_PINFO(MYPNUM, NPROCS)
      END PROGRAM
   ]])],
[m4_default([$1],[:])],
[m4_default([$2],[:])]
)
])


AC_DEFUN([TW_CHECK_BLACS], [
#AC_REQUIRE(ACX_MPI)
tw_blacs_ok=no

# All tests must be run with the MPI fortran compiler.
save_FC=$FC
FC=$MPIFC

AC_ARG_WITH(blacs,
        [AC_HELP_STRING([--with-blacs=<lib>], [use BLACS library <lib>])])
if test x"$with_blacs" != x; then
   BLACS_LIBS="$with_blacs"
fi

save_LIBS=$LIBS
LIBS="$BLACS_LIBS $LIBS"
AC_MSG_CHECKING([if we can compile a BLACS program])
_TW_TRY_BLACS([tw_blacs_ok=yes],[])
AC_MSG_RESULT([$tw_blacs_ok])
LIBS=$save_LIBS

if test $tw_blacs_ok != yes; then
   AC_MSG_CHECKING([for BLACS in -lblacs])
   save_LIBS=$LIBS
   LIBS="-lblacs $LIBS"
   _TW_TRY_BLACS([tw_blacs_ok=yes;BLACS_LIBS=-lblacs],[])
   AC_MSG_RESULT([$tw_blacs_ok])
   LIBS=$save_LIBS
fi

if test $tw_blacs_ok != yes; then
   AC_MSG_CHECKING([for BLACS in -lblacsF77init -lblacs -lblacsF77init])
   save_LIBS=$LIBS
   LIBS="-lblacsF77init -lblacs -lblacsF77init $LIBS"
   _TW_TRY_BLACS([tw_blacs_ok=yes;BLACS_LIBS="-lblacsF77init -lblacs -lblacsF77init"],[])
   AC_MSG_RESULT([$tw_blacs_ok])
   LIBS=$save_LIBS
fi

# Try in Sun performance library
if test $tw_blacs_ok != yes; then
   AC_MSG_CHECKING([for BLACS in -ls3l])
   save_LIBS=$LIBS
   LIBS="-ls3l $LIBS"
   _TW_TRY_BLACS([tw_blacs_ok=yes;BLACS_LIBS="-ls3l"],[])
   AC_MSG_RESULT([$tw_blacs_ok])
   LIBS=$save_LIBS
fi

FC=$save_FC

AC_SUBST(BLACS_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
AS_IF([test $tw_blacs_ok = yes],
      [$1],
      [m4_default([$2],[AC_MSG_ERROR([Cannot compile correctly with BLACS])])]
      )
])dnl TW_CHECK_BLACS

