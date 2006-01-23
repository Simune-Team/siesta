dnl autoconf macros for detecting NetCDF (fortan implementation only)
dnl
AC_DEFUN([_TW_TRY_NETCDF], [
ac_ext=f
AC_LINK_IFELSE(
   [AC_LANG_SOURCE([[
      PROGRAM NETCDF
      CALL NF_CLOSE()
      END PROGRAM
   ]])],
   [m4_ifval([$1],[$1],[])],
   [m4_ifval([$2],[$2],[])]
)
])
dnl
AC_DEFUN([TW_PATH_NETCDF],[
tw_netcdf_ok=no
dnl
AC_LANG_PUSH([Fortran])
dnl
case $with_netcdf in
  yes | "") ;;
  no) tw_netcdf_ok=disable ;;
  -* | */* | *.a | *.so | *.so.* | *.o) NETCDF_LIBS="$with_netcdf" ;;
   *) NETCDF_LIBS="-l$with_netcdf" ;;
esac
dnl
if test $tw_netcdf_ok != disable; then
  if test "x$NETCDF_LIBS" = x; then
     NETCDF_LIBS="-lnetcdf"
  fi
dnl
  AC_MSG_CHECKING([for netcdf])
  save_LIBS="$LIBS"
  LIBS="$LIBS $NETCDF_LIBS"
  _TW_TRY_NETCDF([tw_netcdf_ok=yes],[])
  AC_MSG_RESULT([$tw_netcdf_ok])
  LIBS="$save_LIBS"
fi
dnl
dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
AS_IF([test $tw_netcdf_ok = yes],
      [ifelse([$1],,AC_DEFINE(HAVE_NETCDF,1,[Define if you have NetCDF library.]),[$1])],
      [NETCDF_LIBS="";tw_netcdf_ok=no;$2])
AC_LANG_POP
])
