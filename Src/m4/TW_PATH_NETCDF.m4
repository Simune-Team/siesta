# autoconf macros for detecting NetCDF (fortan implementation only)
#

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

AC_DEFUN([TW_PATH_NETCDF],[
tw_netcdf_ok=no

AC_LANG([Fortran])

AC_ARG_WITH(netCDF,
        [AC_HELP_STRING([--with-netCDF=<lib>], [use netCDF library <lib>])])
if test x"$with_netCDF" != x; then
    NETCDF_LIBS="$with_netCDF"
fi
dnl case $with_netCDF in
dnl     -* | */* | *.a | *.so | *.so.* | *.o) NETCDF_LIBS="$with_netCDF" ;;
dnl     *) NETCDF_LIBS="-l$with_netCDF " ;;
dnl esac

if test "x$NETCDF_LIBS" = x; then
   NETCDF_LIBS="-lnetcdf"
fi

AC_MSG_CHECKING([for netcdf])
save_LIBS="$LIBS"
LIBS="$LIBS $NETCDF_LIBS"
#LIBS="-lnetcdf"
_TW_TRY_NETCDF([tw_netcdf_ok=yes],[])
AC_MSG_RESULT([$tw_netcdf_ok])
LIBS="$save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$tw_netcdf_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_NETCDF,1,[Define if you have NetCDF library.]),[$1])
        :
else
        tw_netcdf_ok=no
        $2
fi

])
