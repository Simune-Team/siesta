AC_DEFUN([_TW_TRY_DC_LAPACK], [
ac_ext=f
AC_LINK_IFELSE(
    [AC_LANG_SOURCE([[
      PROGRAM DC_LAPACK
      CHARACTER        JOBZ, UPLO
      INTEGER          INFO, ITYPE, LDA, LDB, LIWORK, LRWORK, LWORK, N
      DOUBLE PRECISION A(1,1), B(1,1), W(1), WORK(1), RWORK(1)
      CALL ZHEGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,                &
     &            LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
      END PROGRAM
   ]])],
   [m4_ifval([$1],[$1],[])],
   [m4_ifval([$2],[$2],[])]
)
])

AC_DEFUN([TW_CHECK_DC_LAPACK], [
dnl AC_REQUIRE([TW_FIND_LAPACK])
acx_dc_lapack_ok=no
dnl

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
save_LIBS="$LIBS"; LIBS="$LIBS $LAPACK_LIBS $BLAS_LIBS $FLIBS"
AC_MSG_CHECKING([LAPACK includes divide-and-conquer routines])
_TW_TRY_DC_LAPACK([acx_dc_lapack_ok=yes], [:])
AC_MSG_RESULT([$acx_dc_lapack_ok])
LIBS="$save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
AS_IF([test $acx_dc_lapack_ok = yes],
      [m4_default([$1],[:])],
      [m4_default([$2],[AC_MSG_ERROR([Need a more complete Lapack library])])]
     )
])dnl TW_CHECK_DC_LAPACK
