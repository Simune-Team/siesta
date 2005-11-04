# AC_CONFIG_AUX_DIR(DIR)
# ----------------------
# Find install-sh, config.sub, config.guess, and Cygnus configure
# in directory DIR.  These are auxiliary files used in configuration.
# DIR can be either absolute or relative to $srcdir.
AC_DEFUN([AC_CONFIG_AUX_DIR],
[AC_CONFIG_AUX_DIRS($1 $srcdir/$1)])

# AC_CONFIG_AUX_DIRS(DIR ...)
# ---------------------------
# Internal subroutine.
# Search for the configuration auxiliary files in directory list $1.
# We look only for install-sh, so users of AC_PROG_INSTALL
# do not automatically need to distribute the other auxiliary files.
# Edited by tow21@cam.ac.uk -> check for config.sub in absence of install.sh
AC_DEFUN([AC_CONFIG_AUX_DIRS],
[ac_aux_dir=
for ac_dir in $1; do
  if test -f $ac_dir/install-sh; then
    ac_aux_dir=$ac_dir
    ac_install_sh="$ac_aux_dir/install-sh -c"
    break
  elif test -f $ac_dir/install.sh; then
    ac_aux_dir=$ac_dir
    ac_install_sh="$ac_aux_dir/install.sh -c"
    break
  elif test -f $ac_dir/shtool; then
    ac_aux_dir=$ac_dir
    ac_install_sh="$ac_aux_dir/shtool install -c"
    break
  fi
  if test -f $ac_dir/config.sub; then
    ac_aux_dir=$ac_dir
    break
  fi
done
if test -z "$ac_aux_dir"; then
  AC_MSG_ERROR([cannot find install-sh or install.sh in $1])
fi
ac_config_guess="$SHELL $ac_aux_dir/config.guess"
ac_config_sub="$SHELL $ac_aux_dir/config.sub"
ac_configure="$SHELL $ac_aux_dir/configure" # This should be Cygnus configure.
AC_PROVIDE([AC_CONFIG_AUX_DIR_DEFAULT])dnl
])# AC_CONFIG_AUX_DIRS
dnl @synopsis TW_CHECK_FC_FPP([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler supports
dnl cpp-like functionality when called on a suitable fixed-format file.
dnl If so, ACTION_IF_TRUE is performed; if not, ACTION_IF_FALSE
dnl 
dnl @version 1.0
dnl @author Toby White <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_FPP],[
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_FC_SRCEXT(F)
dnl
AC_MSG_CHECKING([whether $FC has an integrated Fortran cpp-style preprocessor for fixed-form source])
dnl
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
      Program test_cpp
#if 1
      Integer i
#else
      Integer j
#endif
      End Program
   ]]),
   [AC_MSG_RESULT([yes])
    m4_ifval([$1],[$1],[])
   ],
   [AC_MSG_RESULT([no])
    m4_ifval([$2],[$2],
                  [AC_MSG_ERROR([A Fortran compiler with integrated cpp-style preprocessor for fixed-form source is required.])])
   ]
)
AC_LANG_POP(Fortran)
dnl
])
dnl A macro to set various compiler-dependent things that can't be sensibly
dnl deduced.

AC_DEFUN([TW_FC_ID_FLAGS], [
AC_REQUIRE([TW_FC_ID])

case $FC_ID in

  Absoft)
     FFLAGS_DEBUG="-et -g -Rb -Rc -Rp -Rs"
     ;;

  Digital)
     FFLAGS_FAST=-O2
     FFLAGS_DEBUG="-g -Rabc -ei"
     ;;

  G77)
     ;;

  Gfortran)
     ;;
 
  Intel)
     FFLAGS_DEBUG="-C -g -inline_debug_info"
     ;;

  Lahey)
     FFLAGS_DEBUG="--chk aesux --chkglobal -g --trace"
     FFLAGS_FAST="-O --warn --quiet --tpp --ntrace"
     ;;

  Nag)
     # This is a hack - we should test for these next two
     FCFLAGS="$FCFLAGS -mismatch -kind=byte"
     FFLAGS_MPI="-kind=byte -mismatch"
     FFLAGS_DEBUG="-C=all -g -gline -nan"
     DEFS="$DEFS __NAG__"
     SYS=nag
     ;;
  
  Portland)
     FFLAGS_DEBUG="-g -Mbounds"
     FFLAGS_FAST="-fast"
     ;;

  SGI)
     FFLAGS_DEBUG="-g -O0"
     FFLAGS_FAST="-O3 -OPT:Olimit=0"
     ;;

  Sun)
     FFLAGS_DEBUG="-C -g"
     FFLAGS_FAST="-fast"
     ;;

  Xlf)
     FFLAGS_DEBUG="-g -C -qinitauto -qsave -qmaxmem=16000 -qnolm"
     FFLAGS_FAST="-O3 -qarch=auto -qtune=auto -qcache=auto -qnolm"
     SYS=xlf
     ;;

esac

AC_SUBST(FFLAGS_MPI)

])
AC_DEFUN([ACX_MPI], [
AC_PREREQ(2.50) dnl for AC_LANG_CASE

AC_LANG_CASE([C], [
        AC_REQUIRE([AC_PROG_CC])
        AC_ARG_VAR(MPICC,[MPI C compiler command])
        AC_CHECK_PROGS(MPICC, mpicc hcc mpcc mpcc_r mpxlc, $CC)
        acx_mpi_save_CC="$CC"
        CC="$MPICC"
        AC_SUBST(MPICC)
],
[C++], [
        AC_REQUIRE([AC_PROG_CXX])
        AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
        AC_CHECK_PROGS(MPICXX, mpiCC mpCC, $CXX)
        acx_mpi_save_CXX="$CXX"
        CXX="$MPICXX"
        AC_SUBST(MPICXX)
],
[Fortran 77], [
        AC_REQUIRE([AC_PROG_F77])
        AC_ARG_VAR(MPIF77,[MPI Fortran 77 compiler command])
        AC_CHECK_PROGS(MPIF77, mpif77 hf77 mpxlf mpf77 mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r, $F77)
        acx_mpi_save_F77="$F77"
        F77="$MPIF77"
        AC_SUBST(MPIF77)
],
[Fortran], [
        AC_REQUIRE([AC_PROG_FC])
        AC_ARG_VAR(MPIFC,[MPI Fortran compiler command])
        AC_CHECK_PROGS(MPIFC, mpifc mpxlf mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r, $FC)
        acx_mpi_save_FC="$FC"
        FC="$MPIFC"
        AC_SUBST(MPIFC)
])

if test x = x"$MPILIBS"; then
        AC_LANG_CASE([C], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
                [C++], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
                [Fortran 77], [AC_MSG_CHECKING([for MPI_Init])
                        AC_TRY_LINK([],[      call MPI_Init], [MPILIBS=" "
                                AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])],
		[Fortran], [AC_MSG_CHECKING([for MPI_Init])
                        AC_LINK_IFELSE([      call MPI_Init], [MPILIBS=" "
                                AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])]
		)
fi
if test x = x"$MPILIBS"; then
        AC_CHECK_LIB(mpi, MPI_Init, [MPILIBS="-lmpi"])
fi
if test x = x"$MPILIBS"; then
        AC_CHECK_LIB(mpich, MPI_Init, [MPILIBS="-lmpich"])
fi

dnl We have to use AC_TRY_COMPILE and not AC_CHECK_HEADER because the
dnl latter uses $CPP, not $CC (which may be mpicc).
AC_LANG_CASE([C], [if test x != x"$MPILIBS"; then
        AC_MSG_CHECKING([for mpi.h])
        AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
                AC_MSG_RESULT(no)])
fi],
[C++], [if test x != x"$MPILIBS"; then
        AC_MSG_CHECKING([for mpi.h])
        AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
                AC_MSG_RESULT(no)])
fi])

AC_LANG_CASE([C], [CC="$acx_mpi_save_CC"],
        [C++], [CXX="$acx_mpi_save_CXX"],
        [Fortran 77], [F77="$acx_mpi_save_F77"],
	[Fortran], [FC="$acx_mpi_save_FC"])

AC_SUBST(MPILIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x = x"$MPILIBS"; then
        $2
        :
else
        ifelse([$1],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$1])
        :
fi
])dnl ACX_MPI
dnl @synopsis TW_CHECK_FC_TR15580([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler is fully
dnl compliant with the Fortran 95 Floating Point Exception Handling
dnl Extension, ISO TR15580.
dnl
dnl @version 1.0
dnl @author <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_TR15580],[
dnl
AC_MSG_CHECKING([$FC for compliance to the Floating Point Exception Handling Extension])
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
      Program test_tr15580

      Use, Intrinsic :: IEEE_Arithmetic
      Use, Intrinsic :: IEEE_Exceptions
      Use, Intrinsic :: IEEE_Features

      End Program test_tr15580
   ]]),   
   [AC_MSG_RESULT([yes])
    m4_ifval([$1],[$1],[])
   ],
   [AC_MSG_RESULT([no])
    m4_ifval([$2],[$2],
                  [AC_MSG_ERROR([A fully TR15580-compliant compiler is required.])])
   ]
)
AC_LANG_POP(Fortran)
dnl
])
dnl A macro to determine which compiler is being used, in order that
dnl different flags can be set

AC_DEFUN([TW_FC_ID], [
AC_REQUIRE([AC_PROG_FC])

FC_ID=

dnl Firstly go by compiler name.

case $FC in 
   
   g77*)
      FC_ID=G77
      ;;

   gfortran*)
      FC_ID=Gfortran
      ;;

   if*)
      FC_ID=Intel
      ;;

   lf9*)
      FC_ID=Lahey
      ;;
   
   pgf*)
      FC_ID=Portland
      ;;

   xlf*)
      FC_ID=Xlf 

esac

dnl then try and disambiguate all f77, f90, and f95 types.
dnl We should have a choice between
dnl nag. absoft. sun. sgi. digital. hp. cray. ...?

if test x$FC_ID = x; then
   tw_fc_v_output=$($FC -V 2>&1 )
   if test $?; then
      case $tw_fc_v_output in
         *NAG*)
            FC_ID=Nag
            ;;
         *Sun*)
            FC_ID=Sun # there's more than one compiler here ...
            ;;
      esac
   fi
fi
 if test x$FC_ID = x; then
   tw_fc_v_output=$($FC -version 2>&1)
   if test $?; then
      case $tw_fc_v_output in
         *Compaq*)
            FC_ID=Digital
            ;;
         *Digital*)
            FC_ID=Digital
            ;;
         *SGI*)
            FC_ID=SGI
            ;;
      esac
   fi
fi   
   
AS_IF([test x$FC_ID != x],
      [AC_MSG_NOTICE([$FC seems to be a $FC_ID compiler])],
      [FC_ID=unknown; AC_MSG_NOTICE([Could not determine type of compiler])])

dnl for more fun, try and get the version number now ...


])# TW_FC_ID
dnl @synopsis TW_CHECK_FC_TR15581([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler is fully
dnl compliant with the Fortran 95 Enhanced Datatype Facilities
dnl Extension, ISO TR15581.
dnl
dnl @version 1.0
dnl @author <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_TR15581],[
dnl
AC_MSG_CHECKING([$FC for compliance to the Enhanced Datatype Facilities Extension])
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
      Program test_tr15581

      Type test_type
        Integer, Allocatable :: array(:)
      End Type test_type

      End Program test_tr15581

      Function test_function
        Integer, Allocatable :: test_function(:)

        Allocate(test_function(5))
      End Function test_function

      Subroutine test_subroutine(array)
        Integer, Allocatable :: array(:)

        Allocate(array(5))
      End Subroutine test_subroutine
   ]]),   
   [AC_MSG_RESULT([yes])
    m4_ifval([$1],[$1],[])
   ],
   [AC_MSG_RESULT([no])
    m4_ifval([$2],[$2],
                  [AC_MSG_ERROR([A fully TR15581-compliant compiler is required.])])
   ]
)
AC_LANG_POP(Fortran)
dnl
])
dnl
AC_DEFUN([TW_FIND_FC_BLAS], [
AC_PREREQ(2.50)
acx_blas_ok=no

AC_LANG_PUSH([Fortran])

AC_ARG_WITH(blas,
	[AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) acx_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac

acx_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

sgemm=sgemm
dgemm=dgemm 

# First, check BLAS_LIBS environment variable
if test "$acx_blas_ok" = no; then
  if test "x$BLAS_LIBS" != x; then
    save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
    AC_MSG_CHECKING([for $sgemm in $BLAS_LIBS])
    AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes], [BLAS_LIBS=""])
    AC_MSG_RESULT($acx_blas_ok)
    LIBS="$save_LIBS"
  fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test "$acx_blas_ok" = no; then
  AC_MSG_CHECKING([is BLAS linked by default])
  save_LIBS="$LIBS"; LIBS="$LIBS"
  AC_LINK_IFELSE(
           AC_LANG_SOURCE([
      Program test_blas
      Call sgemm
      End Program
  ]), [acx_blas_ok=yes])
  AC_MSG_RESULT($acx_blas_ok)
  LIBS="$save_LIBS"
fi

# BLAS in libblasmt.a? (shipped with Lahey Fortran)
if test "$acx_blas_ok" = no; then
  AC_MSG_CHECKING([for BLAS in -lblasmt])
  save_LIBS="$LIBS"; LIBS="$LIBS -lblasmt"
  AC_LINK_IFELSE(
           AC_LANG_SOURCE([
      Program test_blas
      Call sgemm
      End Program
  ]),
      [BLAS_LIBS="-lblasmt"
       acx_blas_ok=yes])
  AC_MSG_RESULT($acx_blas_ok)
  LIBS="$save_LIBS"
fi


# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test "$acx_blas_ok" = no; then
  AC_MSG_CHECKING([for BLAS in ATLAS])
  save_LIBS="$LIBS"; LIBS="$LIBS -lcblas -lf77blas -latlas"
  AC_LINK_IFELSE(
           AC_LANG_SOURCE([
      Program test_blas
      Call sgemm
      End Program
  ]),
      [BLAS_LIBS="-lcblas -lf77blas -latlas"
       acx_blas_ok=yes])
  if test "$acx_blas_ok" = no; then
  LIBS="$save_LIBS"
  save_LIBS="$LIBS"; LIBS="$LIBS -lcblas -lf77blas -latlas -lg2c"
  AC_LINK_IFELSE(
           AC_LANG_SOURCE([
      Program test_blas
      Call sgemm
      End Program
  ]),
      [BLAS_LIBS="-lcblas -lf77blas -latlas -lg2c"
       acx_blas_ok=yes])
  fi
  AC_MSG_RESULT($acx_blas_ok)
  LIBS="$save_LIBS"
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(dgemm, $dgemm,
		[AC_CHECK_LIB(sgemm, $sgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
fi

# BLAS in Alpha CXML library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(cxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(dxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $acx_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, $sgemm,
        			[BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 acx_blas_ok=yes],[],[-lsunmath])])
	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(scs, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $sgemm,
		     [acx_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(essl, $sgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
fi

# BLAS in vecLib Framework (MacOSX 10.2 onwards)
if test $acx_blas_ok = no; then
  LDFLAGS_save="$LDFLAGS"
  for ac_flag in "-framework vecLib" "-Wl,-framework -Wl,vecLib"
  do
    LDFLAGS="$LDFLAGS $ac_flag"
    AC_LINK_IFELSE(
           AC_LANG_SOURCE([
      Program test_blas
      Call sgemm
      End Program
    ]),
    [acx_blas_ok=yes;BLAS_LIBS="$ac_flag"],[])
  LDFLAGS="$LDFLAGS_save"
  done
fi

# Generic BLAS library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lblas"])
fi
dnl
AC_SUBST(BLAS_LIBS)
dnl
LIBS="$acx_blas_save_LIBS"
dnl
# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
AS_IF([test x"$acx_blas_ok" = xyes],
      [$1],
      [m4_default([$2],[AC_MSG_ERROR([Could not find BLAS library])])]
      )
AC_LANG_POP
])
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
LIBS="$LIBS $SCALAPACK_LIBS  $BLACS_LIBS $LAPACK_LIBS $BLAS_LIBS"
AC_MSG_CHECKING([if we can compile a ScaLAPACK program])
_TW_TRY_SCALAPACK([tw_scalapack_ok=yes], [])
LIBS="$save_LIBS"

if test $tw_scalapack_ok = no; then
  LIBS="$LIBS -lscalapack $BLACS_LIBS $LAPACK_LIBS $BLAS_LIBS"
  _TW_TRY_SCALAPACK([tw_scalapack_ok=yes; SCALAPACK_LIBS=-lscalapack], [])
  LIBS="$save_LIBS"
fi

AC_MSG_RESULT([$tw_scalapack_ok])
if test $tw_scalapack_ok = yes; then
   AC_MSG_CHECKING([if ScaLAPACK version is sufficiently recent])
   LIBS="$LIBS $SCALAPACK_LIBS  $BLACS_LIBS $LAPACK_LIBS $BLAS_LIBS"
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
dnl @synopsis TW_CHECK_FC_90([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler is fully
dnl compliant with Fortran 90 (ISO/IEC-1539:1991)
dnl If so, ACTION_IF_TRUE is performed; if not, ACTION_IF_FALSE
dnl 
dnl It currently tests for:
dnl
dnl Modules
dnl Private 
dnl New-style variable declarations.
dnl
dnl @version 1.0
dnl @author Toby White <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_90],[
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_MSG_CHECKING([for Fortran 90 compliance])
dnl
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
      Module test_module

      Implicit None
      Private

      Contains

      Function test_function() Result(out)
      Integer :: out
      out = 0
      End Function test_function

      End Module test_module
   ]]),
   [AC_MSG_RESULT([yes])
    m4_default([$1],[])
   ],
   [AC_MSG_RESULT([no])
    m4_default([$2],
               [AC_MSG_ERROR([ A fully Fortran-90-compliant compiler is required.])])
   ]
)
AC_LANG_POP(Fortran)
dnl
])
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
dnl @synopsis TW_CHECK_FC_FPP_90([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler supports
dnl cpp-like functionality when called on a suitable fixed-format file.
dnl If so, ACTION_IF_TRUE is performed; if not, ACTION_IF_FALSE
dnl 
dnl @version 1.0
dnl @author Toby White <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_FPP_90],[
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_FPP_SRCEXT(F90)
ac_ext=F90
dnl
AC_MSG_CHECKING([whether $FC has an integrated Fortran cpp-style preprocessor for free-form source])
dnl
FCFLAGS_save=$FCFLAGS
FCFLAGS="$FCFLAGS $FPPFLAGS_F90 $FCFLAGS_free"
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
Program test_cpp
#if 1
  Integer i
#else
  Integer j
#endif
End Program
   ]]),
   [AC_MSG_RESULT([yes])
    m4_ifval([$1],[$1],[])
   ],
   [AC_MSG_RESULT([no])
    m4_ifval([$2],[$2],
                  [AC_MSG_ERROR([A Fortran compiler with integrated cpp-style preprocessor for free-from source is required.])])
   ]
)
AC_LANG_POP(Fortran)

FCFLAGS=$FCFLAGS_save
dnl
])
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
AC_MSG_CHECKING([how to compile a call to FLUSH])
dnl
dnl Try first with nothing
dnl
tw_flush_ok=no
dnl
AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_FLUSH_BARE])],
    [tw_flush_ok=yes; TW_FLUSH=bare;tw_method=default;DEFS="$DEFS FC_HAVE_FLUSH"],
    [])
if test $tw_flush_ok = no; then
   save_LDFLAGS=$LDFLAGS
   LDFLAGS="$LDFLAGS -Vaxlib"
   AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_FLUSH_BARE])],
    [tw_flush_ok=yes; TW_FLUSH=INTEL;tw_method="with -Vaxlib";DEFS="$DEFS FC_HAVE_FLUSH"],
    [])
   if test $tw_flush_ok = no; then
      LDFLAGS=$save_LDFLAGS
   fi
fi
if test $tw_flush_ok = no; then
  AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_FLUSH_NAG])],
    [tw_flush_ok=yes; TW_FLUSH=NAG;tw_method="with f90_unix_io";DEFS="$DEFS FC_HAVE_FLUSH"],
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
#FIXME: output of real & int kind query looks ugly,
# because output of FREEFORM is interspersed. How
# to avoid?

# AC_FC_REAL_KIND([KIND_DECLARATION], [VARIABLE_SUFFIX], [ACTION_IF_SUCCESS], [ACTION_IF_FAIL])
# -------------------
#
# This macro checks what integer is produced by the kind 
# declaration KIND_DECLARATION. This integer is placed in 
# AC_FC_KIND_<VARIABLE_SUFFIX>. If we successfully find a
# kind integer, ACTION_IF_SUCCESS is performed; otherwise
# ACTION_IF_FAIL.

AC_DEFUN([AC_FC_REAL_KIND], [dnl
AC_REQUIRE([AC_PROG_FC])
AC_CACHE_CHECK([for kind number produced by $1], 
                ac_cv_fc_real_kind_[]$2[],
[dnl
AC_LANG_PUSH([Fortran])
ac_ext=f90
AC_FC_FREEFORM([f90])
FCFLAGS_save="$FCFLAGS"
FCFLAGS="$FCFLAGS $FCFLAGS_free_f90"
AC_MSG_CHECKING([for kind number produced by $1])
ac_fc_kind_test=1
ac_fc_kind_found=no
while test $ac_fc_kind_test -lt 100
do
  cat > conftest.f90 << ACEOF
dnl The program below will fail to compile if 
dnl sp != mysp; ie if the kind produced by the 
dnl supplied kind declaration ($1) is not the same
dnl same as $ac_fc_kind_test. This is because Fortran
dnl pointers & targets must be of the same kind. 
dnl All conforming compilers must fail to compile the
dnl subroutine otherwise.
dnl
dnl This approach is taken since it enables us to
dnl test for kind numbers at compile time rather
dnl than run time, which means the macro will support
dnl crosss-compilation.
dnl
dnl However, note that kind numbers can theoretically
dnl be anything from 1 to the largest default integer
dnl supported by the compiler. Here we only test up to
dnl 99, which is more than enough on all compilers tried
dnl so far
dnl
subroutine kind_explorer
  integer, parameter :: sp = $1
  integer, parameter :: mysp = $ac_fc_kind_test
  real(kind=sp), target :: x
  real(kind=mysp), pointer :: y
  y=>x
end subroutine kind_explorer
ACEOF
dnl
  eval echo $ac_compile
  if eval $ac_compile 2>&5
  then
    ac_fc_kind_found=yes
    break
  fi
  ac_fc_kind_test=`expr $ac_fc_kind_test + 1`
done
if test "$ac_fc_kind_found" = yes; then
  ac_cv_fc_real_kind_[]$2[]=$ac_fc_kind_test
else 
  ac_cv_fc_real_kind_[]$2[]=none
fi

FCFLAGS="$FCFLAGS_save"
AC_LANG_POP([Fortran])
])
AS_IF([test $ac_cv_fc_real_kind_[]$2[] != no],
      [ac_fc_real_kind_[]$2[]=$ac_cv_fc_real_kind_[]$2[]; m4_default([$3],[:])],
      [m4_default([$4],[AC_MSG_ERROR([Could not find Fortran real kind number for $1])])]
     )
]) # AC_FC_REAL_KIND

# AC_FC_INT_KIND([KIND_DECLARATION], [VARIABLE_SUFFIX], [ACTION_IF_SUCCESS], [ACTION_IF_FAIL])
# -------------------
#
# This macro checks what integer is produced by the kind 
# declaration KIND_DECLARATION. This integer is placed in 
# AC_FC_KIND_<VARIABLE_SUFFIX>. If we successfully find a
# kind integer, ACTION_IF_SUCCESS is performed; otherwise
# ACTION_IF_FAIL.

AC_DEFUN([AC_FC_INT_KIND], [dnl
AC_REQUIRE([AC_PROG_FC])
AC_CACHE_CHECK([for kind number produced by $1], 
                ac_cv_fc_int_kind_[]$2[],
[dnl
AC_LANG_PUSH([Fortran])
ac_ext=f90
AC_FC_FREEFORM([f90])
FCFLAGS_save="$FCFLAGS"
FCFLAGS="$FCFLAGS $FCFLAGS_free_f90"
AC_MSG_CHECKING([for kind number produced by $1])
ac_fc_kind_test=1
ac_fc_kind_found=no
while test $ac_fc_kind_test -lt 100
do
  cat > conftest.f90 << ACEOF
dnl The program below will fail to compile if 
dnl sp != mysp; ie if the kind produced by the 
dnl supplied kind declaration ($1) is not the same
dnl same as $ac_fc_kind_test. This is because Fortran
dnl pointers & targets must be of the same kind. 
dnl All conforming compilers must fail to compile the
dnl subroutine otherwise.
dnl
dnl This approach is taken since it enables us to
dnl test for kind numbers at compile time rather
dnl than run time, which means the macro will support
dnl crosss-compilation.
dnl
dnl However, note that kind numbers can theoretically
dnl be anything from 1 to the largest default integer
dnl supported by the compiler. Here we only test up to
dnl 99, which is more than enough on all compilers tried
dnl so far
dnl
subroutine kind_explorer
  integer, parameter :: sp = $1
  integer, parameter :: mysp = $ac_fc_kind_test
  integer(kind=sp), target :: x
  integer(kind=mysp), pointer :: y
  y=>x
end subroutine kind_explorer
ACEOF
dnl
  eval echo $ac_compile
  if eval $ac_compile 2>&5
  then
    ac_fc_kind_found=yes
    break
  fi
  ac_fc_kind_test=`expr $ac_fc_kind_test + 1`
done
if test "$ac_fc_kind_found" = yes; then
  ac_cv_fc_int_kind_[]$2[]=$ac_fc_kind_test
else 
  ac_cv_fc_int_kind_[]$2[]=none
fi

FCFLAGS="$FCFLAGS_save"
AC_LANG_POP([Fortran])
])
AS_IF([test $ac_cv_fc_int_kind_[]$2[] != no],
      [ac_fc_int_kind_[]$2[]=$ac_cv_fc_int_kind_[]$2[]; m4_default([$3],[:])],
      [m4_default([$4],[AC_MSG_ERROR([Could not find Fortran integer kind number for $1])])]
     )
]) # AC_FC_INT_KIND


# AC_FC_GET_REAL_KINDS([ACTION_IF_SUCCESS], [ACTION_IF_FAIL])
# -------------------
#
# This macro attempts to find the Fortran compiler's kinds
# for the following four types of real number:
#  Compiler default (single) precision
#  Compiler double precision
#  IEEE single precision
#  IEEE double precision
# The first two are guaranteed to exist; the second two may
# or may not.
# If all 4 are succesfully found,. ACTION_IF_SUCCESS is
# performed.
# Otherwise, ACTION_IF_FAIL is performed
#
AC_DEFUN([AC_FC_GET_REAL_KINDS], [dnl
AC_REQUIRE([AC_PROG_FC])

ac_fc_got_kinds=yes

AC_FC_REAL_KIND([[kind(1.0)]], [sp], 
                [], [ac_got_kinds=no])
AC_FC_REAL_KIND([[kind(1.0d0)]], [dp], 
                [], [ac_got_kinds=no])
AC_FC_REAL_KIND([[selected_real_kind(6,34)]], [ieee_sp], 
                [], [ac_got_kinds=no])
AC_FC_REAL_KIND([[selected_real_kind(15,300)]], [ieee_dp], 
                [], [ac_got_kinds=no])

AS_IF([test $ac_fc_got_kinds != no],
      [m4_default([$1],[:])],
      [m4_default([$2],[AC_MSG_ERROR([Could not find all Fortran real kinds])])]
      )
]) dnl AC_FC_GET_REAL_KINDS

