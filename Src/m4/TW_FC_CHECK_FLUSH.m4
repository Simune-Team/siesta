dnl Check how to get at the flush intrinsic.

AC_DEFUN([_TW_TRY_FLUSH_BARE],
[
      PROGRAM TESTFLUSH
      PRINT*
      CALL FLUSH(5)
      END PROGRAM TESTFLUSH
])
AC_DEFUN([_TW_TRY_FLUSH_NAG],
[
      PROGRAM TESTFLUSH
      USE F90_UNIX_IO, ONLY:FLUSH
      PRINT*
      CALL FLUSH(5)
      END PROGRAM TESTFLUSH
])
AC_DEFUN([_TW_TRY_FLUSH_XLF],
[
      PROGRAM TESTFLUSH
      PRINT*
      CALL FLUSH_(5)
      END PROGRAM TESTFLUSH
])

AC_DEFUN([TW_FC_CHECK_FLUSH], [
AC_REQUIRE([AC_PROG_FC])dnl
dnl
AC_MSG_CHECKING([how to compile a FLUSH subroutine])
dnl
dnl Try first with nothing
dnl
tw_flush_ok=no
dnl
AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_FLUSH_BARE])],
    [tw_flush_ok=yes; TW_FLUSH=bare;tw_method=default],
    [])
if test $tw_flush_ok = no; then
   save_LDFLAGS=$LDFLAGS
   LDFLAGS="$LDFLAGS -Vaxlib"
   AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_FLUSH_BARE])],
    [tw_flush_ok=yes; TW_FLUSH=INTEL;tw_method="with -Vaxlib"],
    [])
   if test $tw_flush_ok = no; then
      LDFLAGS=$save_LDFLAGS
   fi
fi
if test $tw_flush_ok = no; then
  AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_FLUSH_NAG])],
    [tw_flush_ok=yes; TW_FLUSH=NAG;tw_method="with f90_unix_io"],
    [])
fi
if test $tw_flush_ok = no; then
  AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_FLUSH_XLF])],
    [tw_flush_ok=yes; TW_FLUSH=XLF;tw_method="with underscore"],
    [])
fi
AC_MSG_RESULT([$tw_method])
dnl
AS_IF([test $tw_flush_ok = yes],
      [$1],
      [m4_default([$2],[AC_MSG_ERROR([Cannot compile FLUSH statement])])]
     )
dnl
])# TW_FC_CHECK_FLUSH
