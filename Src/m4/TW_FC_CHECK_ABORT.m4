dnl Check how to get at the abort intrinsic.

AC_DEFUN([_TW_TRY_ABORT_BARE],
[
      PROGRAM TESTABORT
      CALL ABORT
      END PROGRAM TESTABORT
])
AC_DEFUN([_TW_TRY_ABORT_NAG],
[
      PROGRAM TESTABORT
      USE F90_UNIX_PROC, ONLY:ABORT
      CALL ABORT
      END PROGRAM TESTABORT
])
AC_DEFUN([_TW_TRY_ABORT_XLF],
[
      PROGRAM TESTABORT
      CALL ABORT_
      END PROGRAM TESTABORT
])

AC_DEFUN([TW_FC_CHECK_ABORT], [
AC_REQUIRE([AC_PROG_FC])dnl
dnl
AC_MSG_CHECKING([how to compile a call to ABORT])
dnl
dnl Try first with nothing
dnl
tw_abort_ok=no
dnl
AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_ABORT_BARE])],
    [tw_abort_ok=yes; TW_ABORT=bare;tw_method=default;DEFS="$DEFS FC_HAVE_ABORT"],
    [])
if test $tw_abort_ok = no; then
   save_LDFLAGS=$LDFLAGS
   LDFLAGS="$LDFLAGS -Vaxlib"
   AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_ABORT_BARE])],
    [tw_abort_ok=yes; TW_ABORT=INTEL;tw_method="with -Vaxlib";DEFS="$DEFS FC_HAVE_ABORT"],
    [])
   if test $tw_abort_ok = no; then
      LDFLAGS=$save_LDFLAGS
   fi
fi
if test $tw_abort_ok = no; then
  AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_ABORT_NAG])],
    [tw_abort_ok=yes; TW_ABORT=NAG;tw_method="with f90_unix_proc";DEFS="$DEFS FC_HAVE_ABORT"],
    [])
fi
if test $tw_abort_ok = no; then
  AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_ABORT_XLF])],
    [tw_abort_ok=yes; TW_ABORT=XLF;tw_method="with underscore"],
    [])
fi
AC_MSG_RESULT([$tw_method])
dnl
AS_IF([test $tw_abort_ok = yes],
      [$1],
      [m4_default([$2],[AC_MSG_ERROR([Cannot compile call to ABORT ])])]
     )
dnl
])# TW_FC_CHECK_ABORT
