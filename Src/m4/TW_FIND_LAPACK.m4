dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_lapack.html
dnl

AC_DEFUN([_TW_TRY_LAPACK], [
ac_ext=f
AC_LINK_IFELSE(
    [AC_LANG_SOURCE([[
      PROGRAM LAPACK
      CHARACTER        JOBZ, UPLO
      INTEGER          INFO, ITYPE, LDA, LDB, LWORK, N
      DOUBLE PRECISION A(1,1), B(1,1), W(1), WORK(1)
      CALL DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,                 &
     &            LWORK, INFO )
      END PROGRAM
   ]])],
   [m4_ifval([$1],[$1],[])],
   [m4_ifval([$2],[$2],[])]
)
])

AC_DEFUN([_TW_TRY_LAPACK_OTHER], [
ac_ext=f
AC_LINK_IFELSE(
    [AC_LANG_SOURCE([[
      PROGRAM LAPACK
      CHARACTER        JOBZ, UPLO
      INTEGER          INFO, ITYPE, LDA, LDB, N
      DOUBLE PRECISION A(1,1), B(1,1)
      CALL DSYGST( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, INFO )
      END PROGRAM
   ]])],
   [m4_ifval([$1],[$1],[])],
   [m4_ifval([$2],[$2],[])]
)
])


AC_DEFUN([TW_FIND_LAPACK], [
AC_REQUIRE([TW_FIND_FC_BLAS])
acx_lapack_ok=no
dnl
AC_ARG_WITH(lapack,
        [AC_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])
case $with_lapack in
        yes | "") ;;
        no) acx_lapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$with_lapack" ;;
        *) LAPACK_LIBS="-l$with_lapack" ;;
esac
dnl
# Get fortran linker name of LAPACK function to check for.
#AC_FC_FUNC(dsygv)
dsygv=dsygv
dnl
# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
        AC_MSG_CHECKING([for $dsygv in $LAPACK_LIBS])
        _TW_TRY_LAPACK([acx_lapack_ok=yes],[:])
        AC_MSG_RESULT($acx_lapack_ok)
        LIBS="$save_LIBS"
        if test acx_lapack_ok = no; then
                LAPACK_LIBS=""
        fi
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $acx_lapack_ok = no; then
        save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
        AC_MSG_CHECKING([LAPACK already linked])
        _TW_TRY_LAPACK([acx_lapack_ok=yes], [:])
	AC_MSG_RESULT([$acx_lapack_ok])
        LIBS="$save_LIBS"
fi

# Generic LAPACK library?
for lapack in lapackmt lapack lapack_rs6k; do
        if test $acx_lapack_ok = no; then
                save_LIBS="$LIBS"; LIBS="-l$lapack $BLAS_LIBS $LIBS"
                AC_MSG_CHECKING([for LAPACK in -l$lapack])
                _TW_TRY_LAPACK([acx_lapack_ok=yes; LAPACK_LIBS="-l$lapack"], [:])
                AC_MSG_RESULT([$acx_lapack_ok])
                LIBS="$save_LIBS"
        fi
done

AC_SUBST(LAPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_lapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
        :
else
        acx_lapack_ok=no
        $2
fi
])dnl ACX_LAPACK
