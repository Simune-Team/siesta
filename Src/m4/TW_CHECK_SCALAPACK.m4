dnl Macro to find what paramaters must be passed to the compiler
dnl and linker to compile with scalapack. Currently only works
dnl with fortran.
dnl Macro to check that correct flags have been chosen for 
dnl ScaLAPACK compilation. Only works with Fortran at the moment.

dnl Must use fixed+free format-compatible line continuation, to avoid
dnl problems with picky compilers. & in column 6 & 81

AC_DEFUN([_TW_TRY_SCALAPACK], [
ac_ext=f
AC_LINK_IFELSE(
   [AC_LANG_SOURCE([[
      PROGRAM SCALAPACK                                                    
      CHARACTER        JOBZ, RANGE, UPLO                                   
      INTEGER          IA, IB, IBTYPE, IL, INFO, IU, IZ, JA, JB, JZ,             & 
     &                 LIWORK, LWORK, M, N, NZ                                   
      DOUBLE PRECISION ABSTOL, ORFAC, VL, VU                                     
      INTEGER          DESCA(1), DESCB(1), DESCZ(1), ICLUSTR(1),                 &
     &                 IFAIL(1), IWORK(1)                              
      DOUBLE PRECISION A(1), B(1), GAP(1), W(1), WORK(1), Z(1) 
      CALL PDSYGVX( IBTYPE, JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA,  B,          &
     &              IB,  JB,  DESCB,  VL, VU, IL, IU, ABSTOL, M, NZ, W,          &
     &              ORFAC,  Z,  IZ,  JZ,  DESCZ,  WORK,  LWORK,  IWORK,          &
     &              LIWORK, IFAIL, ICLUSTR, GAP, INFO )                    
      END PROGRAM
   ]])],
   [m4_ifval([$1],[$1],[])],
   [m4_ifval([$2],[$2],[])]
)
])

AC_DEFUN([_TW_TRY_SCALAPACK_VN], [
ac_ext=f
AC_LINK_IFELSE(
   [AC_LANG_SOURCE([[
      PROGRAM SCALAPACK                                                    
      CHARACTER        NAME, OPTS                                   
      INTEGER          ICTXT, ISPEC, N1, N2, N3, N4
      CALL PJLAENV( ICTXT, ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      END PROGRAM
   ]])],
   [m4_ifval([$1],[$1],[])],
   [m4_ifval([$2],[$2],[])]
)
])


AC_DEFUN([TW_CHECK_SCALAPACK], [
#AC_REQUIRE([TW_FIND_BLAS])
#AC_REQUIRE([TW_CHECK_BLACS])
tw_scalapack_ok=no

# All tests must be run with the MPI fortran compiler.
save_FC=$FC
FC=$MPIFC

AC_ARG_WITH(scalapack,
        [AC_HELP_STRING([--with-scalapack=<lib>], [use ScaLAPACK library <lib>])])
case $with_scalapack in
        yes | "") ;;
        no) tw_scalapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) SCALAPACK_LIBS="$with_scalapack" ;;
        *) SCALAPACK_LIBS="-l$with_scalapack" ;;
esac

save_LIBS=$LIBS
LIBS="$LIBS $SCALAPACK_LIBS $BLACS_LIBS $BLAS_LIBS"
AC_MSG_CHECKING([if we can compile a ScaLAPACK program])
_TW_TRY_SCALAPACK([tw_scalapack_ok=yes], [])
LIBS="$save_LIBS"

if test $tw_scalapack_ok = no; then
  LIBS="$LIBS -lscalapack $BLACS_LIBS $BLAS_LIBS"
  _TW_TRY_SCALAPACK([tw_scalapack_ok=yes; SCALAPACK_LIBS=-lscalapack], [])
  LIBS="$save_LIBS"
fi

AC_MSG_RESULT([$tw_scalapack_ok])
if test $tw_scalapack_ok = yes; then
   AC_MSG_CHECKING([if ScaLAPACK version is sufficiently recent])
   LIBS="$LIBS $SCALAPACK_LIBS $BLACS_LIBS $BLAS_LIBS"
   _TW_TRY_SCALAPACK_VN([AC_MSG_RESULT([yes])],
                        [COMP_LIBS="scalapack_extra.o $COMP_LIBS";
                         AC_MSG_RESULT([no - using additional SIESTA routines])])
   LIBS="$save_LIBS"
fi
FC=$save_FC

AC_SUBST(SCALAPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
AS_IF([test $tw_scalapack_ok = yes],
      [$1],
      [m4_default([$2],[AC_MSG_ERROR([Cannot compile correctly with ScaLAPACK])])]
     )
])# TW_CHECK_SCALAPACK
