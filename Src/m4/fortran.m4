# This file is part of Autoconf.                       -*- Autoconf -*-
# Fortran languages support.
# Copyright (C) 2001, 2003, 2004
# Free Software Foundation, Inc.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.
#
# As a special exception, the Free Software Foundation gives unlimited
# permission to copy, distribute and modify the configure scripts that
# are the output of Autoconf.  You need not follow the terms of the GNU
# General Public License when using or distributing such scripts, even
# though portions of the text of Autoconf appear in them.  The GNU
# General Public License (GPL) does govern all other use of the material
# that constitutes the Autoconf program.
#
# Certain portions of the Autoconf source text are designed to be copied
# (in certain cases, depending on the input) into the output of
# Autoconf.  We call these the "data" portions.  The rest of the Autoconf
# source text consists of comments plus executable code that decides which
# of the data portions to output in any given case.  We call these
# comments and executable code the "non-data" portions.  Autoconf never
# copies any of the non-data portions into its output.
#
# This special exception to the GPL applies to versions of Autoconf
# released by the Free Software Foundation.  When you make and
# distribute a modified version of Autoconf, you may extend this special
# exception to the GPL to apply to your modified version as well, *unless*
# your modified version has the potential to copy into its output some
# of the text that was the non-data portion of the version that you started
# with.  (In other words, unless your change moves or copies text from
# the non-data portions to the data portions.)  If your modification has
# such potential, you must delete any notice of this special exception
# to the GPL from your modified version.
#
# Written by David MacKenzie, with help from
# Franc,ois Pinard, Karl Berry, Richard Pixley, Ian Lance Taylor,
# Roland McGrath, Noah Friedman, david d zuhn, and many others.
#
# Extended by Martin Wilck (2000), Norman Gray and Toby White (2004),
# to add support for preprocessable Fortran (fpp).


# XXX Temporary notes on the fpp work, which should turn into a ChangeLog
# entry on submission to the mainline.
#
# Original fpp work done by Martin Wilck: see
# <http://sources.redhat.com/ml/autoconf-patches/2000-07/msg00287.html>
# and the following thread.
#
# The modifications to Martin's work, done by Norman Gray and Toby White,
# are as follows:
#
# - Change define to m4_define, to match newer (post 2.50) versions of autoconf
# - Change .F to .fpp (on case-insensitive filesystems such as HFS+,
#   file.f and file.F are the same file)
# - Added second, optional, FPP-SRC-EXT parameter to AC_PROG_FPP, to set the
#   extension expected, in case you (persistent!) want to use .F
#   instead of .fpp
# - Slight changes to tests, following comments by Akim Demaille
#   on the mailing list (including Martin's cunning
#   _AC_SHELL_POSITIONAL).
# - renamed $ac_gnu_compiler to $ac_compiler_gnu
# - Added _AC_LANG_ABBREV(Preprocessed Fortran), and a few other
#   functions which match the support apparatus seems to want
# - Added AC_LANG_PREPROC and AC_LANG_COMPILER
# - Added AC_REQUIRE of _AC_PROG_FC_CPP to _AC_PROG_FPP
# - In test of how to run the preprocessor alone, check that our tests
#   don't produce a a.out file -- if they do, that counts as a
#   failure, even with a zero exit status
# - Removed _AC_ECHO calls -- were always internal, now disappeared
# - Changed old >&AC_FD_LOG to >&AS_MESSAGE_LOG_FD
# - Based on "Fortran" language, rather than "Fortran 77"
# - Addition of separate yes/no check for whether FC supports preprocessing.
# - Remove hardcoding of file extensions in preprocessing checks.
# - Improvement of AC_FC_SRCEXT to allow a variable to be passed to it.
# - Add utility macros AC_FC_FIXEDFORM, AC_FC_OPEN_SPECIFIERS,
#   AC_FC_CHECK_INTRINSICS, AC_FC_RECL_UNIT
#
# Further modifications by Toby White
# - rework AC_FC_FREEFORM & AC_FC_FIXEDFORM to take a (mandatory)
#   file extension argument. 
# - add AC_FPP_FREEFORM & AC_FPP_FIXEDFORM along similar lines
# - deprecate AC_FC_SRCEXT


# Fortran vs. Fortran 77:
#   This file contains macros for both "Fortran 77" and "Fortran", where
# the former is the "classic" autoconf Fortran interface and is intended
# for legacy F77 codes, while the latter is intended to support newer Fortran
# dialects.  Fortran 77 uses environment variables F77, FFLAGS, and FLIBS,
# while Fortran uses FC, FCFLAGS, and FCLIBS.  For each user-callable AC_*
# macro, there is generally both an F77 and an FC version, where both versions
# share the same _AC_*_FC_* backend.  This backend macro requires that
# the appropriate language be AC_LANG_PUSH'ed, and uses _AC_LANG_ABBREV and
# _AC_LANG_PREFIX in order to name cache and environment variables, etc.


# _AC_LIST_MEMBER_IF(ELEMENT, LIST, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ---------------------------------------------------------------------------
#
# Processing the elements of a list is tedious in shell programming,
# as lists tend to be implemented as space delimited strings.
#
# This macro searches LIST for ELEMENT, and executes ACTION-IF-FOUND
# if ELEMENT is a member of LIST, otherwise it executes
# ACTION-IF-NOT-FOUND.
AC_DEFUN([_AC_LIST_MEMBER_IF],
[dnl Do some sanity checking of the arguments.
m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl
  ac_exists=false
  for ac_i in $2; do
    if test x"$1" = x"$ac_i"; then
      ac_exists=true
      break
    fi
  done

  AS_IF([test x"$ac_exists" = xtrue], [$3], [$4])[]dnl
])# _AC_LIST_MEMBER_IF


# _AC_LINKER_OPTION(LINKER-OPTIONS, SHELL-VARIABLE)
# -------------------------------------------------
#
# Specifying options to the compiler (whether it be the C, C++ or
# Fortran 77 compiler) that are meant for the linker is compiler
# dependent.  This macro lets you give options to the compiler that
# are meant for the linker in a portable, compiler-independent way.
#
# This macro take two arguments, a list of linker options that the
# compiler should pass to the linker (LINKER-OPTIONS) and the name of
# a shell variable (SHELL-VARIABLE).  The list of linker options are
# appended to the shell variable in a compiler-dependent way.
#
# For example, if the selected language is C, then this:
#
#   _AC_LINKER_OPTION([-R /usr/local/lib/foo], foo_LDFLAGS)
#
# will expand into this if the selected C compiler is gcc:
#
#   foo_LDFLAGS="-Xlinker -R -Xlinker /usr/local/lib/foo"
#
# otherwise, it will expand into this:
#
#   foo_LDFLAGS"-R /usr/local/lib/foo"
#
# You are encouraged to add support for compilers that this macro
# doesn't currently support.
# FIXME: Get rid of this macro.
AC_DEFUN([_AC_LINKER_OPTION],
[if test "$ac_compiler_gnu" = yes; then
  for ac_link_opt in $1; do
    $2="[$]$2 -Xlinker $ac_link_opt"
  done
else
  $2="[$]$2 $1"
fi[]dnl
])# _AC_LINKER_OPTION



## ----------------------- ##
## 1. Language selection.  ##
## ----------------------- ##


# -------------------------- #
# 1d. The Fortran language.  #
# -------------------------- #


# AC_LANG(Fortran 77)
# -------------------
m4_define([AC_LANG(Fortran 77)],
[ac_ext=f
ac_compile='$F77 -c $FFLAGS conftest.$ac_ext >&AS_MESSAGE_LOG_FD'
ac_link='$F77 -o conftest$ac_exeext $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&AS_MESSAGE_LOG_FD'
ac_compiler_gnu=$ac_cv_f77_compiler_gnu
])


# AC_LANG(Fortran)
# ----------------
m4_define([AC_LANG(Fortran)],
[ac_ext=${FC_SRCEXT-f}
ac_compile='$FC -c $FCFLAGS $FCFLAGS_SRCEXT conftest.$ac_ext >&AS_MESSAGE_LOG_FD'
ac_link='$FC -o conftest$ac_exeext $FCFLAGS $LDFLAGS $FCFLAGS_SRCEXT conftest.$ac_ext $LIBS >&AS_MESSAGE_LOG_FD'
ac_compiler_gnu=$ac_cv_fc_compiler_gnu
])

# AC_LANG(Preprocessed Fortran)
# --------------------------------
# We need a separate `preprocessed' language, because not all Fortran 
# compilers have a preprocessor built in.  Therefore we may need to
# resort to an `indirect' compilation, .fpp->.f->.o, including the
# generation of a suitable extra build rule.  The language extension
# is set in macro AC_PROG_FPP, to $FPP_SRC_EXT.
m4_define([AC_LANG(Preprocessed Fortran)],
[ac_ext=$FPP_SRC_EXT
# We need to use variables because compilation depends on whether 
# $F77 supports direct compilation of source with cpp directives
ac_compile=$ac_fpp_compile
ac_link=$ac_fpp_link
ac_compiler_gnu=$ac_cv_fc_compiler_gnu
])


# AC_LANG_FORTRAN77
# -----------------
AU_DEFUN([AC_LANG_FORTRAN77], [AC_LANG(Fortran 77)])


# _AC_FORTRAN_ASSERT
# ------------------
# Current language must be Fortran, Fortran 77, or preprocessed Fortran.
# FIXME: is there any reason why this can't be AC_LANG_CASE?
m4_defun([_AC_FORTRAN_ASSERT],
[m4_if(_AC_LANG, [Fortran], [],
       [m4_if(_AC_LANG, [Fortran 77], [],
              [m4_if(_AC_LANG, [Preprocessed Fortran], []
                     [m4_fatal([$0: current language is not Fortran: ] _AC_LANG)])])])])


# _AC_LANG_ABBREV(Fortran 77)
# ---------------------------
m4_define([_AC_LANG_ABBREV(Fortran 77)], [f77])

# _AC_LANG_ABBREV(Fortran)
# ------------------------
m4_define([_AC_LANG_ABBREV(Fortran)], [fc])

# _AC_LANG_ABBREV(Preprocessed Fortran)
# -------------------------------------
m4_define([_AC_LANG_ABBREV(Preprocessed Fortran)], [fpp])


# _AC_LANG_PREFIX(Fortran 77)
# ---------------------------
m4_define([_AC_LANG_PREFIX(Fortran 77)], [F])

# _AC_LANG_PREFIX(Fortran)
# ------------------------
m4_define([_AC_LANG_PREFIX(Fortran)], [FC])

# _AC_LANG_PREFIX(Preprocessed Fortran)
# -------------------------------------
m4_define([_AC_LANG_PREFIX(Preprocessed Fortran)], [FPP])


# _AC_FC
# ------
# Return F77, FC or PPFC, depending upon the language.
AC_DEFUN([_AC_FC],
[_AC_FORTRAN_ASSERT()dnl
AC_LANG_CASE([Fortran 77],           [F77],
             [Fortran],              [FC],
             [Preprocessed Fortran], [PPFC])])


## ---------------------- ##
## 2.Producing programs.  ##
## ---------------------- ##


# --------------------- #
# 2d. Fortran sources.  #
# --------------------- #

# AC_LANG_SOURCE(Fortran 77)(BODY)
# AC_LANG_SOURCE(Fortran)(BODY)
# --------------------------------
# FIXME: Apparently, according to former AC_TRY_COMPILER, the CPP
# directives must not be included.  But AC_TRY_RUN_NATIVE was not
# avoiding them, so?
m4_define([AC_LANG_SOURCE(Fortran 77)],
[$1])
m4_define([AC_LANG_SOURCE(Fortran)],
[$1])
m4_define([AC_LANG_SOURCE(Preprocessed Fortran)],
[$1])


# AC_LANG_PROGRAM(Fortran 77)([PROLOGUE], [BODY])
# -----------------------------------------------
# Yes, we discard the PROLOGUE.
m4_define([AC_LANG_PROGRAM(Fortran 77)],
[m4_ifval([$1],
       [m4_warn([syntax], [$0: ignoring PROLOGUE: $1])])dnl
      program main
$2
      end])


# AC_LANG_PROGRAM(Fortran)([PROLOGUE], [BODY])
# -----------------------------------------------
# FIXME: can the PROLOGUE be used?
m4_define([AC_LANG_PROGRAM(Fortran)],
[m4_ifval([$1],
       [m4_warn([syntax], [$0: ignoring PROLOGUE: $1])])dnl
      program main
$2
      end])


# AC_LANG_PROGRAM(Preprocessed Fortran)([PROLOGUE], [BODY])
# ---------------------------------------------------------
# FIXME: can the PROLOGUE be used?
m4_define([AC_LANG_PROGRAM(Preprocessed Fortran)],
[$1
      program main
$2
      end])


# AC_LANG_CALL(Fortran 77)(PROLOGUE, FUNCTION)
# --------------------------------------------
# FIXME: This is a guess, help!
# FIXME: ...but it's a good guesss -- what's the problem?
m4_define([AC_LANG_CALL(Fortran 77)],
[AC_LANG_PROGRAM([$1],
[      call $2])])


# AC_LANG_CALL(Fortran)(PROLOGUE, FUNCTION)
# --------------------------------------------
# FIXME: This is a guess, help!
m4_define([AC_LANG_CALL(Fortran)],
[AC_LANG_PROGRAM([$1],
[      call $2])])


# AC_LANG_CALL(Preprocessed Fortran)(PROLOGUE, FUNCTION)
# ------------------------------------------------------
# FIXME: This is a guess, help!
m4_define([AC_LANG_CALL(Preprocessed Fortran)],
[AC_LANG_PROGRAM([$1],
[      call $2])])


# AC_LANG_FUNC_LINK_TRY(Fortran)(FUNCTION)
# ----------------------------------
# FIXME: This is a severely broken implementation.
# It does not distinguish between functions and subroutines.
# It ignores any arguments.
# It is, however, necessary, to make AC_CHECK_FUNC work.
m4_define([AC_LANG_FUNC_LINK_TRY(Fortran)],
[AC_LANG_PROGRAM([],[
c     What to do if $1 is something weird
c     Either already declared as fortran keyword
c     Or something needing quoting?
c#define $1 innocuous_$1

      Program test
      External $1
      Call $1
      End Program

])])


## -------------------------------------------- ##
## 3. Looking for Compilers and Preprocessors.  ##
## -------------------------------------------- ##


# -------------------------- #
# 3d. The Fortran compiler.  #
# -------------------------- #


# AC_LANG_PREPROC(Fortran 77)
# ---------------------------
# Find the Fortran 77 preprocessor.  Must be AC_DEFUN'd to be AC_REQUIRE'able.
AC_DEFUN([AC_LANG_PREPROC(Fortran 77)],
[m4_warn([syntax],
	 [$0: No preprocessor defined for ]_AC_LANG)])

# AC_LANG_PREPROC(Fortran)
# ---------------------------
# Find the Fortran preprocessor.  Must be AC_DEFUN'd to be AC_REQUIRE'able.
AC_DEFUN([AC_LANG_PREPROC(Fortran)],
[m4_warn([syntax],
         [$0: No preprocessor defined for ]_AC_LANG)])

# AC_LANG_PREPROC(Preprocessed Fortran)
# -------------------------------------
# Find the Fortran preprocessor.  Must be AC_DEFUN'd to be AC_REQUIRE'able.
AC_DEFUN([AC_LANG_PREPROC(Preprocessed Fortran)],
[AC_REQUIRE([AC_PROG_FPP])])


# AC_LANG_COMPILER(Fortran 77)
# ----------------------------
# Find the Fortran 77 compiler.  Must be AC_DEFUN'd to be
# AC_REQUIRE'able.
AC_DEFUN([AC_LANG_COMPILER(Fortran 77)],
[AC_REQUIRE([AC_PROG_F77])])

# AC_LANG_COMPILER(Fortran)
# -------------------------
# Find the Fortran compiler.  Must be AC_DEFUN'd to be
# AC_REQUIRE'able.
AC_DEFUN([AC_LANG_COMPILER(Fortran)],
[AC_REQUIRE([AC_PROG_FC])])

# AC_LANG_COMPILER(Preprocessed Fortran)
# --------------------------------------
# Find the Fortran compiler.  Must be AC_DEFUN'd to be
# AC_REQUIRE'able.
AC_DEFUN([AC_LANG_COMPILER(Preprocessed Fortran)],
[AC_REQUIRE([AC_PROG_FC])])


# ac_cv_prog_g77
# --------------
# We used to name the cache variable this way.
AU_DEFUN([ac_cv_prog_g77],
[ac_cv_f77_compiler_gnu])


# _AC_FC_DIALECT_YEAR([DIALECT])
# ------------------------------
# Given a Fortran DIALECT, which is Fortran [YY]YY or simply [YY]YY,
# convert to a 4-digit year.  The dialect must be one of Fortran 77,
# 90, 95, or 2000, currently.  If DIALECT is simply Fortran or the
# empty string, returns the empty string.
AC_DEFUN([_AC_FC_DIALECT_YEAR],
[m4_case(m4_bpatsubsts(m4_tolower([$1]), [fortran],[], [ *],[]),
	 [77],[1977], [1977],[1977],
	 [90],[1990], [1990],[1990],
	 [95],[1995], [1995],[1995],
	 [2000],[2000],
         [],[],
         [m4_fatal([unknown Fortran dialect])])])


# _AC_PROG_FC([DIALECT], [COMPILERS...])
# --------------------------------------
# DIALECT is a Fortran dialect, given by Fortran [YY]YY or simply [YY]YY,
# and must be one of those supported by _AC_FC_DIALECT_YEAR
#
# If DIALECT is supplied, then we search for compilers of that dialect
# first, and then later dialects.  Otherwise, we search for compilers
# of the newest dialect first, and then earlier dialects in increasing age.
# This search order is necessarily imperfect because the dialect cannot
# always be inferred from the compiler name.
#
# Known compilers:
#  f77/f90/f95: generic compiler names
#  g77: GNU Fortran 77 compiler
#  gfortran: putative GNU Fortran 95+ compiler (in progress)
#  fort77: native F77 compiler under HP-UX (and some older Crays)
#  frt: Fujitsu F77 compiler
#  pgf77/pgf90/pgf95: Portland Group F77/F90/F95 compilers
#  xlf/xlf90/xlf95: IBM (AIX) F77/F90/F95 compilers
#  lf95: Lahey-Fujitsu F95 compiler
#  fl32: Microsoft Fortran 77 "PowerStation" compiler
#  af77: Apogee F77 compiler for Intergraph hardware running CLIX
#  epcf90: "Edinburgh Portable Compiler" F90
#  fort: Compaq (now HP) Fortran 90/95 compiler for Tru64 and Linux/Alpha
#  ifc: Intel Fortran 95 compiler for Linux/x86
#  efc: Intel Fortran 95 compiler for IA64
#
# Must check for lf95 before f95 - some Lahey versions ship an f95 binary 
# in the default path that must be avoided.
# Similarly for xlf95.
# (A proper fix would involve being able to go back and try another compiler
#  if the first one fails, but that requires a major reworking of much of
#  autoconf. The same problem arises (with no easy solution) on some Digital
#  compilers: f95 fails on .F files, f90 succeeds.)
#
m4_define([_AC_F95_FC], [xlf95 lf95 gfortran f95 fort ifc efc pgf95])
m4_define([_AC_F90_FC], [f90 xlf90 pgf90 epcf90])
m4_define([_AC_F77_FC], [g77 f77 xlf frt pgf77 fort77 fl32 af77])
AC_DEFUN([_AC_PROG_FC],
[_AC_FORTRAN_ASSERT()dnl
AC_CHECK_TOOLS([]_AC_FC[],
      m4_default([$2],
	m4_case(_AC_FC_DIALECT_YEAR([$1]),
		[1995], [_AC_F95_FC],
		[1990], [_AC_F90_FC _AC_F95_FC],
		[1977], [_AC_F77_FC _AC_F90_FC _AC_F95_FC],
		[_AC_F95_FC _AC_F90_FC _AC_F77_FC])))

# Provide some information about the compiler.
echo "$as_me:__oline__:" \
     "checking for _AC_LANG compiler version" >&AS_MESSAGE_LOG_FD
ac_compiler=`set X $ac_compile; echo $[2]`
_AC_EVAL([$ac_compiler --version </dev/null >&AS_MESSAGE_LOG_FD])
_AC_EVAL([$ac_compiler -v </dev/null >&AS_MESSAGE_LOG_FD])
_AC_EVAL([$ac_compiler -V </dev/null >&AS_MESSAGE_LOG_FD])
rm -f a.out

m4_expand_once([_AC_COMPILER_EXEEXT])[]dnl
m4_expand_once([_AC_COMPILER_OBJEXT])[]dnl
_AC_FC_MOD_SUFFIX dnl This will fail if F77-only compiler is used.
dnl FIXME how should we fix that? Currently warns & continues.
# If we don't use `.F' as extension, the preprocessor is not run on the
# input file.  (Note that this only needs to work for GNU compilers.)
ac_save_ext=$ac_ext
ac_ext=F
_AC_LANG_COMPILER_GNU
ac_ext=$ac_save_ext
_AC_PROG_FC_G
])# _AC_PROG_FC


# _AC_FC_MOD_SUFFIX
# -----------------
# Determines the form of the filename of modules produced
# by the Fortran compiler.
# Tests for all forms of file extension I've (TW) found in the
# wild. Note that at least one compiler (PGI??) changes the
# case of the basename as well. Whether this happens is
# encoded in the variable ac_fc_mod_uppercase.
#
AC_DEFUN([_AC_FC_MOD_SUFFIX],
[
ac_mod_ext=
ac_fc_mod_uppercase=no
cat > conftest.$ac_ext << \_ACEOF
      module conftest
       implicit none
       integer :: i
      end module conftest
_ACEOF
_AC_EVAL_STDERR($ac_compile)
AC_MSG_CHECKING([for suffix of module files])
for ac_mod_file in conftest.mod conftest.MOD conftest.M CONFTEST.MOD none	 
do
  if test -f $ac_mod_file; then
    break;
  fi
done
rm -f conftest.$ac_ext conftest.$ac_exe_ext conftest.mod conftest.MOD conftest.M CONFTEST.MOD
#
case $ac_mod_file in
  conftest.mod)
    ac_mod_ext=mod
    ;;
  conftest.MOD)
    ac_mod_ext=MOD
    ;;
  conftest.M)
    ac_mod_ext=M
    ;;
  CONFTEST.MOD)
    ac_mod_ext=MOD
    ac_fc_mod_uppercase=yes
    ;;
  none)
    AC_MSG_WARN([Could not find Fortran module file extension.])
    ;;
esac

if test $ac_mod_file != none; then
  AC_MSG_RESULT([$ac_mod_ext])
fi
if test $ac_fc_mod_uppercase = yes; then
  AC_MSG_NOTICE([Fortran module filenames are uppercase.])
fi
])


# AC_PROG_F77([COMPILERS...])
# ---------------------------
# COMPILERS is a space separated list of Fortran 77 compilers to search
# for.  See also _AC_PROG_FC.
AC_DEFUN([AC_PROG_F77],
[AC_LANG_PUSH(Fortran 77)dnl
AC_ARG_VAR([F77],    [Fortran 77 compiler command])dnl
AC_ARG_VAR([FFLAGS], [Fortran 77 compiler flags])dnl
_AC_ARG_VAR_LDFLAGS()dnl
_AC_PROG_FC([Fortran 77], [$1])
G77=`test $ac_compiler_gnu = yes && echo yes`
AC_LANG_POP(Fortran 77)dnl
])# AC_PROG_F77


# AC_PROG_FC([COMPILERS...], [DIALECT])
# -------------------------------------
# COMPILERS is a space separated list of Fortran 77 compilers to search
# for, and [DIALECT] is an optional dialect.  See also _AC_PROG_FC.
AC_DEFUN([AC_PROG_FC],
[AC_BEFORE([$0], [AC_PROG_FPP])dnl
AC_LANG_PUSH(Fortran)dnl
AC_ARG_VAR([FC],    [Fortran compiler command])dnl
AC_ARG_VAR([FCFLAGS], [Fortran compiler flags])dnl
_AC_ARG_VAR_LDFLAGS()dnl
_AC_PROG_FC([$2], [$1])
AC_LANG_POP(Fortran)dnl
])# AC_PROG_FC


# _AC_PROG_FC_G
# -------------
# Check whether -g works, even if F[C]FLAGS is set, in case the package
# plays around with F[C]FLAGS (such as to build both debugging and normal
# versions of a library), tasteless as that idea is.
m4_define([_AC_PROG_FC_G],
[_AC_FORTRAN_ASSERT()dnl
ac_test_FFLAGS=${[]_AC_LANG_PREFIX[]FLAGS+set}
ac_save_FFLAGS=$[]_AC_LANG_PREFIX[]FLAGS
_AC_LANG_PREFIX[]FLAGS=
AC_CACHE_CHECK(whether $[]_AC_FC[] accepts -g, ac_cv_prog_[]_AC_LANG_ABBREV[]_g,
[_AC_LANG_PREFIX[]FLAGS=-g
_AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
[ac_cv_prog_[]_AC_LANG_ABBREV[]_g=yes],
[ac_cv_prog_[]_AC_LANG_ABBREV[]_g=no])
])
if test "$ac_test_FFLAGS" = set; then
  _AC_LANG_PREFIX[]FLAGS=$ac_save_FFLAGS
elif test $ac_cv_prog_[]_AC_LANG_ABBREV[]_g = yes; then
  if test "x$ac_cv_[]_AC_LANG_ABBREV[]_compiler_gnu" = xyes; then
    _AC_LANG_PREFIX[]FLAGS="-g -O2"
  else
    _AC_LANG_PREFIX[]FLAGS="-g"
  fi
else
  if test "x$ac_cv_[]_AC_LANG_ABBREV[]_compiler_gnu" = xyes; then
    _AC_LANG_PREFIX[]FLAGS="-O2"
  else
    _AC_LANG_PREFIX[]FLAGS=
  fi
fi[]dnl
])# _AC_PROG_FC_G


# _AC_PROG_FC_C_O
# ---------------
# Test if the Fortran compiler accepts the options `-c' and `-o'
# simultaneously, and define `[F77/FC]_NO_MINUS_C_MINUS_O' if it does not.
#
# The usefulness of this macro is questionable, as I can't really see
# why anyone would use it.  The only reason I include it is for
# completeness, since a similar test exists for the C compiler.
#
# FIXME: it seems like we could merge the C/Fortran versions of this.
AC_DEFUN([_AC_PROG_FC_C_O],
[_AC_FORTRAN_ASSERT()dnl
AC_CACHE_CHECK([whether $[]_AC_FC[] understands -c and -o together],
               [ac_cv_prog_[]_AC_LANG_ABBREV[]_c_o],
[AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])
# We test twice because some compilers refuse to overwrite an existing
# `.o' file with `-o', although they will create one.
ac_try='$[]_AC_FC[] $[]_AC_LANG_PREFIX[]FLAGS -c conftest.$ac_ext -o conftest.$ac_objext >&AS_MESSAGE_LOG_FD'
if AC_TRY_EVAL(ac_try) &&
     test -f conftest.$ac_objext &&
     AC_TRY_EVAL(ac_try); then
  ac_cv_prog_[]_AC_LANG_ABBREV[]_c_o=yes
else
  ac_cv_prog_[]_AC_LANG_ABBREV[]_c_o=no
fi
rm -f conftest*])
if test $ac_cv_prog_[]_AC_LANG_ABBREV[]_c_o = no; then
  AC_DEFINE([]_AC_FC[]_NO_MINUS_C_MINUS_O, 1,
            [Define to 1 if your Fortran compiler doesn't accept
             -c and -o together.])
fi
])# _AC_PROG_FC_C_O


# AC_PROG_F77_C_O
# ---------------
AC_DEFUN([AC_PROG_F77_C_O],
[AC_REQUIRE([AC_PROG_F77])dnl
AC_LANG_PUSH(Fortran 77)dnl
_AC_PROG_FC_C_O
AC_LANG_POP(Fortran 77)dnl
])# AC_PROG_F77_C_O


# AC_PROG_FC_C_O
# ---------------
AC_DEFUN([AC_PROG_FC_C_O],
[AC_REQUIRE([AC_PROG_FC])dnl
AC_LANG_PUSH(Fortran)dnl
_AC_PROG_FC_C_O
AC_LANG_POP(Fortran)dnl
])# AC_PROG_FC_C_O


## ------------------------------- ##
## 4. Compilers' characteristics.  ##
## ------------------------------- ##


# ---------------------------------------- #
# 4d. Fortran 77 compiler characteristics. #
# ---------------------------------------- #


# _AC_PROG_FC_V_OUTPUT([FLAG = $ac_cv_prog_{f77/fc}_v])
# -------------------------------------------------
# Link a trivial Fortran program, compiling with a verbose output FLAG
# (whose default value, $ac_cv_prog_{f77/fc}_v, is computed by
# _AC_PROG_FC_V), and return the output in $ac_{f77/fc}_v_output.  This
# output is processed in the way expected by _AC_FC_LIBRARY_LDFLAGS,
# so that any link flags that are echoed by the compiler appear as
# space-separated items.
AC_DEFUN([_AC_PROG_FC_V_OUTPUT],
[_AC_FORTRAN_ASSERT()dnl
AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])

# Compile and link our simple test program by passing a flag (argument
# 1 to this macro) to the Fortran compiler in order to get
# "verbose" output that we can then parse for the Fortran linker
# flags.
ac_save_FFLAGS=$[]_AC_LANG_PREFIX[]FLAGS
_AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS m4_default([$1], [$ac_cv_prog_[]_AC_LANG_ABBREV[]_v])"
(eval echo $as_me:__oline__: \"$ac_link\") >&AS_MESSAGE_LOG_FD
ac_[]_AC_LANG_ABBREV[]_v_output=`eval $ac_link AS_MESSAGE_LOG_FD>&1 2>&1 | grep -v 'Driving:'`
echo "$ac_[]_AC_LANG_ABBREV[]_v_output" >&AS_MESSAGE_LOG_FD
_AC_LANG_PREFIX[]FLAGS=$ac_save_FFLAGS

rm -f conftest*

# On HP/UX there is a line like: "LPATH is: /foo:/bar:/baz" where
# /foo, /bar, and /baz are search directories for the Fortran linker.
# Here, we change these into -L/foo -L/bar -L/baz (and put it first):
ac_[]_AC_LANG_ABBREV[]_v_output="`echo $ac_[]_AC_LANG_ABBREV[]_v_output |
	grep 'LPATH is:' |
	sed 's,.*LPATH is\(: *[[^ ]]*\).*,\1,;s,: */, -L/,g'` $ac_[]_AC_LANG_ABBREV[]_v_output"

case $ac_[]_AC_LANG_ABBREV[]_v_output in
  # If we are using xlf then replace all the commas with spaces.
  *xlfentry*)
    ac_[]_AC_LANG_ABBREV[]_v_output=`echo $ac_[]_AC_LANG_ABBREV[]_v_output | sed 's/,/ /g'` ;;

  # With Intel ifc, ignore the quoted -mGLOB_options_string stuff (quoted
  # $LIBS confuse us, and the libraries appear later in the output anyway).
  *mGLOB_options_string*)
    ac_[]_AC_LANG_ABBREV[]_v_output=`echo $ac_[]_AC_LANG_ABBREV[]_v_output | sed 's/\"-mGLOB[[^\"]]*\"/ /g'` ;;

  # If we are using Cray Fortran then delete quotes.
  # Use "\"" instead of '"' for font-lock-mode.
  # FIXME: a more general fix for quoted arguments with spaces?
  *cft90*)
    ac_[]_AC_LANG_ABBREV[]_v_output=`echo $ac_[]_AC_LANG_ABBREV[]_v_output | sed "s/\"//g"` ;;
esac

])# _AC_PROG_FC_V_OUTPUT


# _AC_PROG_FC_V
# --------------
#
# Determine the flag that causes the Fortran compiler to print
# information of library and object files (normally -v)
# Needed for _AC_FC_LIBRARY_FLAGS
# Some compilers don't accept -v (Lahey: -verbose, xlf: -V, Fujitsu: -###)
AC_DEFUN([_AC_PROG_FC_V],
[_AC_FORTRAN_ASSERT()dnl
AC_CACHE_CHECK([how to get verbose linking output from $[]_AC_FC[]],
                [ac_cv_prog_[]_AC_LANG_ABBREV[]_v],
[AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
[ac_cv_prog_[]_AC_LANG_ABBREV[]_v=
# Try some options frequently used verbose output
for ac_verb in -v -verbose --verbose -V -\#\#\#; do
  _AC_PROG_FC_V_OUTPUT($ac_verb)
  # look for -l* and *.a constructs in the output
  for ac_arg in $ac_[]_AC_LANG_ABBREV[]_v_output; do
     case $ac_arg in
        [[\\/]]*.a | ?:[[\\/]]*.a | -[[lLRu]]*)
          ac_cv_prog_[]_AC_LANG_ABBREV[]_v=$ac_verb
          break 2 ;;
     esac
  done
done
if test -z "$ac_cv_prog_[]_AC_LANG_ABBREV[]_v"; then
   AC_MSG_WARN([cannot determine how to obtain linking information from $[]_AC_FC[]])
fi],
                  [AC_MSG_WARN([compilation failed])])
])])# _AC_PROG_FC_V


# _AC_FC_LIBRARY_LDFLAGS
# ----------------------
#
# Determine the linker flags (e.g. "-L" and "-l") for the Fortran
# intrinsic and run-time libraries that are required to successfully
# link a Fortran program or shared library.  The output variable
# FLIBS/FCLIBS is set to these flags.
#
# This macro is intended to be used in those situations when it is
# necessary to mix, e.g. C++ and Fortran, source code into a single
# program or shared library.
#
# For example, if object files from a C++ and Fortran compiler must
# be linked together, then the C++ compiler/linker must be used for
# linking (since special C++-ish things need to happen at link time
# like calling global constructors, instantiating templates, enabling
# exception support, etc.).
#
# However, the Fortran intrinsic and run-time libraries must be
# linked in as well, but the C++ compiler/linker doesn't know how to
# add these Fortran libraries.  Hence, the macro
# "AC_F77_LIBRARY_LDFLAGS" was created to determine these Fortran
# libraries.
#
# This macro was packaged in its current form by Matthew D. Langston.
# However, nearly all of this macro came from the "OCTAVE_FLIBS" macro
# in "octave-2.0.13/aclocal.m4", and full credit should go to John
# W. Eaton for writing this extremely useful macro.  Thank you John.
AC_DEFUN([_AC_FC_LIBRARY_LDFLAGS],
[_AC_FORTRAN_ASSERT()dnl
_AC_PROG_FC_V
AC_CACHE_CHECK([for Fortran libraries of $[]_AC_FC[]], ac_cv_[]_AC_LANG_ABBREV[]_libs,
[if test "x$[]_AC_LANG_PREFIX[]LIBS" != "x"; then
  ac_cv_[]_AC_LANG_ABBREV[]_libs="$[]_AC_LANG_PREFIX[]LIBS" # Let the user override the test.
else

_AC_PROG_FC_V_OUTPUT

ac_cv_[]_AC_LANG_ABBREV[]_libs=

# Save positional arguments (if any)
ac_save_positional="$[@]"

set X $ac_[]_AC_LANG_ABBREV[]_v_output
while test $[@%:@] != 1; do
  shift
  ac_arg=$[1]
  case $ac_arg in
        [[\\/]]*.a | ?:[[\\/]]*.a)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
              ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg")
          ;;
        -bI:*)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
             [_AC_LINKER_OPTION([$ac_arg], ac_cv_[]_AC_LANG_ABBREV[]_libs)])
          ;;
          # Ignore these flags.
        -lang* | -lcrt*.o | -lc | -lgcc | -libmil | -LANG:=*)
          ;;
        -lkernel32)
          test x"$CYGWIN" != xyes && ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg"
          ;;
        -[[LRuY]])
          # These flags, when seen by themselves, take an argument.
          # We remove the space between option and argument and re-iterate
          # unless we find an empty arg or a new option (starting with -)
	  case $[2] in
	     "" | -*);;
	     *)
		ac_arg="$ac_arg$[2]"
		shift; shift
		set X $ac_arg "$[@]"
		;;
	  esac
          ;;
        -YP,*)
          for ac_j in `echo $ac_arg | sed -e 's/-YP,/-L/;s/:/ -L/g'`; do
            _AC_LIST_MEMBER_IF($ac_j, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
                               [ac_arg="$ac_arg $ac_j"
                               ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_j"])
          done
          ;;
        -[[lLR]]*)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
                             ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg")
          ;;
          # Ignore everything else.
  esac
done
# restore positional arguments
set X $ac_save_positional; shift

# We only consider "LD_RUN_PATH" on Solaris systems.  If this is seen,
# then we insist that the "run path" must be an absolute path (i.e. it
# must begin with a "/").
case `(uname -sr) 2>/dev/null` in
   "SunOS 5"*)
      ac_ld_run_path=`echo $ac_[]_AC_LANG_ABBREV[]_v_output |
                        sed -n 's,^.*LD_RUN_PATH *= *\(/[[^ ]]*\).*$,-R\1,p'`
      test "x$ac_ld_run_path" != x &&
        _AC_LINKER_OPTION([$ac_ld_run_path], ac_cv_[]_AC_LANG_ABBREV[]_libs)
      ;;
esac
fi # test "x$[]_AC_LANG_PREFIX[]LIBS" = "x"
])
[]_AC_LANG_PREFIX[]LIBS="$ac_cv_[]_AC_LANG_ABBREV[]_libs"
AC_SUBST([]_AC_LANG_PREFIX[]LIBS)
])# _AC_FC_LIBRARY_LDFLAGS


# AC_F77_LIBRARY_LDFLAGS
# ----------------------
AC_DEFUN([AC_F77_LIBRARY_LDFLAGS],
[AC_REQUIRE([AC_PROG_F77])dnl
AC_LANG_PUSH(Fortran 77)dnl
_AC_FC_LIBRARY_LDFLAGS
AC_LANG_POP(Fortran 77)dnl
])# AC_F77_LIBRARY_LDFLAGS


# AC_FC_LIBRARY_LDFLAGS
# ----------------------
AC_DEFUN([AC_FC_LIBRARY_LDFLAGS],
[AC_REQUIRE([AC_PROG_FC])dnl
AC_LANG_PUSH(Fortran)dnl
_AC_FC_LIBRARY_LDFLAGS
AC_LANG_POP(Fortran)dnl
])# AC_FC_LIBRARY_LDFLAGS


# _AC_FC_DUMMY_MAIN([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------------
#
# Detect name of dummy main routine required by the Fortran libraries,
# (if any) and define {F77,FC}_DUMMY_MAIN to this name (which should be
# used for a dummy declaration, if it is defined).  On some systems,
# linking a C program to the Fortran library does not work unless you
# supply a dummy function called something like MAIN__.
#
# Execute ACTION-IF-NOT-FOUND if no way of successfully linking a C
# program with the {F77,FC} libs is found; default to exiting with an error
# message.  Execute ACTION-IF-FOUND if a dummy routine name is needed
# and found or if it is not needed (default to defining {F77,FC}_DUMMY_MAIN
# when needed).
#
# What is technically happening is that the Fortran libraries provide
# their own main() function, which usually initializes Fortran I/O and
# similar stuff, and then calls MAIN__, which is the entry point of
# your program.  Usually, a C program will override this with its own
# main() routine, but the linker sometimes complain if you don't
# provide a dummy (never-called) MAIN__ routine anyway.
#
# Of course, programs that want to allow Fortran subroutines to do
# I/O, etcetera, should call their main routine MAIN__() (or whatever)
# instead of main().  A separate autoconf test (_AC_FC_MAIN) checks
# for the routine to use in this case (since the semantics of the test
# are slightly different).  To link to e.g. purely numerical
# libraries, this is normally not necessary, however, and most C/C++
# programs are reluctant to turn over so much control to Fortran.  =)
#
# The name variants we check for are (in order):
#   MAIN__ (g77, MAIN__ required on some systems; IRIX, MAIN__ optional)
#   MAIN_, __main (SunOS)
#   MAIN _MAIN __MAIN main_ main__ _main (we follow DDD and try these too)
AC_DEFUN([_AC_FC_DUMMY_MAIN],
[_AC_FORTRAN_ASSERT()dnl
m4_define(_AC_LANG_PROGRAM_C_[]_AC_FC[]_HOOKS,
[#ifdef ]_AC_FC[_DUMMY_MAIN
]AC_LANG_CASE([Fortran], [#ifndef FC_DUMMY_MAIN_EQ_F77])
[#  ifdef __cplusplus
     extern "C"
#  endif
   int ]_AC_FC[_DUMMY_MAIN() { return 1; }
]AC_LANG_CASE([Fortran], [#endif])
[#endif
])
AC_CACHE_CHECK([for dummy main to link with Fortran libraries],
               ac_cv_[]_AC_LANG_ABBREV[]_dummy_main,
[ac_[]_AC_LANG_ABBREV[]_dm_save_LIBS=$LIBS
 LIBS="$LIBS $[]_AC_LANG_PREFIX[]LIBS"
 ac_fortran_dm_var=[]_AC_FC[]_DUMMY_MAIN
 AC_LANG_PUSH(C)dnl

 # First, try linking without a dummy main:
 AC_LINK_IFELSE([AC_LANG_PROGRAM([], [])],
                [ac_cv_fortran_dummy_main=none],
                [ac_cv_fortran_dummy_main=unknown])

 if test $ac_cv_fortran_dummy_main = unknown; then
   for ac_func in MAIN__ MAIN_ __main MAIN _MAIN __MAIN main_ main__ _main; do
     AC_LINK_IFELSE([AC_LANG_PROGRAM([[@%:@define $ac_fortran_dm_var $ac_func]])],
                    [ac_cv_fortran_dummy_main=$ac_func; break])
   done
 fi
 AC_LANG_POP(C)dnl
 ac_cv_[]_AC_LANG_ABBREV[]_dummy_main=$ac_cv_fortran_dummy_main
 rm -f conftest*
 LIBS=$ac_[]_AC_LANG_ABBREV[]_dm_save_LIBS
])
[]_AC_FC[]_DUMMY_MAIN=$ac_cv_[]_AC_LANG_ABBREV[]_dummy_main
AS_IF([test "$[]_AC_FC[]_DUMMY_MAIN" != unknown],
      [m4_default([$1],
[if test $[]_AC_FC[]_DUMMY_MAIN != none; then
  AC_DEFINE_UNQUOTED([]_AC_FC[]_DUMMY_MAIN, $[]_AC_FC[]_DUMMY_MAIN,
                     [Define to dummy `main' function (if any) required to
                      link to the Fortran libraries.])
  if test "x$ac_cv_fc_dummy_main" = "x$ac_cv_f77_dummy_main"; then
	AC_DEFINE([FC_DUMMY_MAIN_EQ_F77], 1,
                  [Define if F77 and FC dummy `main' functions are identical.])
  fi
fi])],
      [m4_default([$2],
            [AC_MSG_FAILURE([linking to Fortran libraries from C fails])])])
])# _AC_FC_DUMMY_MAIN


# AC_F77_DUMMY_MAIN
# ----------------------
AC_DEFUN([AC_F77_DUMMY_MAIN],
[AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])dnl
AC_LANG_PUSH(Fortran 77)dnl
_AC_FC_DUMMY_MAIN
AC_LANG_POP(Fortran 77)dnl
])# AC_F77_DUMMY_MAIN


# AC_FC_DUMMY_MAIN
# ----------------------
AC_DEFUN([AC_FC_DUMMY_MAIN],
[AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])dnl
AC_LANG_PUSH(Fortran)dnl
_AC_FC_DUMMY_MAIN
AC_LANG_POP(Fortran)dnl
])# AC_FC_DUMMY_MAIN


# _AC_FC_MAIN
# -----------
# Define {F77,FC}_MAIN to name of alternate main() function for use with
# the Fortran libraries.  (Typically, the libraries may define their
# own main() to initialize I/O, etcetera, that then call your own
# routine called MAIN__ or whatever.)  See _AC_FC_DUMMY_MAIN, above.
# If no such alternate name is found, just define {F77,FC}_MAIN to main.
#
AC_DEFUN([_AC_FC_MAIN],
[_AC_FORTRAN_ASSERT()dnl
AC_CACHE_CHECK([for alternate main to link with Fortran libraries],
               ac_cv_[]_AC_LANG_ABBREV[]_main,
[ac_[]_AC_LANG_ABBREV[]_m_save_LIBS=$LIBS
 LIBS="$LIBS $[]_AC_LANG_PREFIX[]LIBS"
 ac_fortran_dm_var=[]_AC_FC[]_DUMMY_MAIN
 AC_LANG_PUSH(C)dnl
 ac_cv_fortran_main="main" # default entry point name
 for ac_func in MAIN__ MAIN_ __main MAIN _MAIN __MAIN main_ main__ _main; do
   AC_LINK_IFELSE([AC_LANG_PROGRAM([@%:@ifdef FC_DUMMY_MAIN_EQ_F77
@%:@  undef F77_DUMMY_MAIN
@%:@  undef FC_DUMMY_MAIN
@%:@else
@%:@  undef $ac_fortran_dm_var
@%:@endif
@%:@define main $ac_func])],
                  [ac_cv_fortran_main=$ac_func; break])
 done
 AC_LANG_POP(C)dnl
 ac_cv_[]_AC_LANG_ABBREV[]_main=$ac_cv_fortran_main
 rm -f conftest*
 LIBS=$ac_[]_AC_LANG_ABBREV[]_m_save_LIBS
])
AC_DEFINE_UNQUOTED([]_AC_FC[]_MAIN, $ac_cv_[]_AC_LANG_ABBREV[]_main,
                   [Define to alternate name for `main' routine that is
                    called from a `main' in the Fortran libraries.])
])# _AC_FC_MAIN


# AC_F77_MAIN
# -----------
AC_DEFUN([AC_F77_MAIN],
[AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])dnl
AC_LANG_PUSH(Fortran 77)dnl
_AC_FC_MAIN
AC_LANG_POP(Fortran 77)dnl
])# AC_F77_MAIN


# AC_FC_MAIN
# ----------
AC_DEFUN([AC_FC_MAIN],
[AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])dnl
AC_LANG_PUSH(Fortran)dnl
_AC_FC_MAIN
AC_LANG_POP(Fortran)dnl
])# AC_FC_MAIN


# __AC_FC_NAME_MANGLING
# ---------------------
# Test for the name mangling scheme used by the Fortran compiler.
#
# Sets ac_cv_{f77,fc}_mangling. The value contains three fields, separated
# by commas:
#
# lower case / upper case:
#    case translation of the Fortran symbols
# underscore / no underscore:
#    whether the compiler appends "_" to symbol names
# extra underscore / no extra underscore:
#    whether the compiler appends an extra "_" to symbol names already
#    containing at least one underscore
#
AC_DEFUN([__AC_FC_NAME_MANGLING],
[_AC_FORTRAN_ASSERT()dnl
AC_CACHE_CHECK([for Fortran name-mangling scheme],
               ac_cv_[]_AC_LANG_ABBREV[]_mangling,
[AC_COMPILE_IFELSE(
[      subroutine foobar()
      return
      end
      subroutine foo_bar()
      return
      end],
[mv conftest.$ac_objext cfortran_test.$ac_objext

  ac_save_LIBS=$LIBS
  LIBS="cfortran_test.$ac_objext $LIBS $[]_AC_LANG_PREFIX[]LIBS"

  AC_LANG_PUSH(C)dnl
  ac_success=no
  for ac_foobar in foobar FOOBAR; do
    for ac_underscore in "" "_"; do
      ac_func="$ac_foobar$ac_underscore"
      AC_LINK_IFELSE([AC_LANG_CALL([], [$ac_func])],
		     [ac_success=yes; break 2])
    done
  done
  AC_LANG_POP(C)dnl

  if test "$ac_success" = "yes"; then
     case $ac_foobar in
	foobar)
	   ac_case=lower
	   ac_foo_bar=foo_bar
	   ;;
	FOOBAR)
	   ac_case=upper
	   ac_foo_bar=FOO_BAR
	   ;;
     esac

     AC_LANG_PUSH(C)dnl
     ac_success_extra=no
     for ac_extra in "" "_"; do
	ac_func="$ac_foo_bar$ac_underscore$ac_extra"
	AC_LINK_IFELSE([AC_LANG_CALL([], [$ac_func])],
		       [ac_success_extra=yes; break])
     done
     AC_LANG_POP(C)dnl

     if test "$ac_success_extra" = "yes"; then
	ac_cv_[]_AC_LANG_ABBREV[]_mangling="$ac_case case"
        if test -z "$ac_underscore"; then
           ac_cv_[]_AC_LANG_ABBREV[]_mangling="$ac_cv_[]_AC_LANG_ABBREV[]_mangling, no underscore"
	else
           ac_cv_[]_AC_LANG_ABBREV[]_mangling="$ac_cv_[]_AC_LANG_ABBREV[]_mangling, underscore"
        fi
        if test -z "$ac_extra"; then
           ac_cv_[]_AC_LANG_ABBREV[]_mangling="$ac_cv_[]_AC_LANG_ABBREV[]_mangling, no extra underscore"
	else
           ac_cv_[]_AC_LANG_ABBREV[]_mangling="$ac_cv_[]_AC_LANG_ABBREV[]_mangling, extra underscore"
        fi
      else
	ac_cv_[]_AC_LANG_ABBREV[]_mangling="unknown"
      fi
  else
     ac_cv_[]_AC_LANG_ABBREV[]_mangling="unknown"
  fi

  LIBS=$ac_save_LIBS
  rm -f cfortran_test* conftest*],
  [AC_MSG_FAILURE([cannot compile a simple Fortran program])])
])
])# __AC_FC_NAME_MANGLING

# The replacement is empty.
AU_DEFUN([AC_F77_NAME_MANGLING], [])


# _AC_F77_NAME_MANGLING
# ----------------------
AC_DEFUN([_AC_F77_NAME_MANGLING],
[AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])dnl
AC_REQUIRE([AC_F77_DUMMY_MAIN])dnl
AC_LANG_PUSH(Fortran 77)dnl
__AC_FC_NAME_MANGLING
AC_LANG_POP(Fortran 77)dnl
])# _AC_F77_NAME_MANGLING


# _AC_FC_NAME_MANGLING
# ----------------------
AC_DEFUN([_AC_FC_NAME_MANGLING],
[AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])dnl
AC_REQUIRE([AC_FC_DUMMY_MAIN])dnl
AC_LANG_PUSH(Fortran)dnl
__AC_FC_NAME_MANGLING
AC_LANG_POP(Fortran)dnl
])# _AC_FC_NAME_MANGLING


# _AC_FC_WRAPPERS
# ---------------
# Defines C macros {F77,FC}_FUNC(name,NAME) and {F77,FC}_FUNC_(name,NAME) to
# properly mangle the names of C identifiers, and C identifiers with
# underscores, respectively, so that they match the name mangling
# scheme used by the Fortran compiler.
AC_DEFUN([_AC_FC_WRAPPERS],
[_AC_FORTRAN_ASSERT()dnl
AH_TEMPLATE(_AC_FC[_FUNC],
    [Define to a macro mangling the given C identifier (in lower and upper
     case), which must not contain underscores, for linking with Fortran.])dnl
AH_TEMPLATE(_AC_FC[_FUNC_],
    [As ]_AC_FC[_FUNC, but for C identifiers containing underscores.])dnl
case $ac_cv_[]_AC_LANG_ABBREV[]_mangling in
  "lower case, no underscore, no extra underscore")
          AC_DEFINE(_AC_FC[_FUNC(name,NAME)],  [name])
          AC_DEFINE(_AC_FC[_FUNC_(name,NAME)], [name]) ;;
  "lower case, no underscore, extra underscore")
          AC_DEFINE(_AC_FC[_FUNC(name,NAME)],  [name])
          AC_DEFINE(_AC_FC[_FUNC_(name,NAME)], [name ## _]) ;;
  "lower case, underscore, no extra underscore")
          AC_DEFINE(_AC_FC[_FUNC(name,NAME)],  [name ## _])
          AC_DEFINE(_AC_FC[_FUNC_(name,NAME)], [name ## _]) ;;
  "lower case, underscore, extra underscore")
          AC_DEFINE(_AC_FC[_FUNC(name,NAME)],  [name ## _])
          AC_DEFINE(_AC_FC[_FUNC_(name,NAME)], [name ## __]) ;;
  "upper case, no underscore, no extra underscore")
          AC_DEFINE(_AC_FC[_FUNC(name,NAME)],  [NAME])
          AC_DEFINE(_AC_FC[_FUNC_(name,NAME)], [NAME]) ;;
  "upper case, no underscore, extra underscore")
          AC_DEFINE(_AC_FC[_FUNC(name,NAME)],  [NAME])
          AC_DEFINE(_AC_FC[_FUNC_(name,NAME)], [NAME ## _]) ;;
  "upper case, underscore, no extra underscore")
          AC_DEFINE(_AC_FC[_FUNC(name,NAME)],  [NAME ## _])
          AC_DEFINE(_AC_FC[_FUNC_(name,NAME)], [NAME ## _]) ;;
  "upper case, underscore, extra underscore")
          AC_DEFINE(_AC_FC[_FUNC(name,NAME)],  [NAME ## _])
          AC_DEFINE(_AC_FC[_FUNC_(name,NAME)], [NAME ## __]) ;;
  *)
          AC_MSG_WARN([unknown Fortran name-mangling scheme])
          ;;
esac
])# _AC_FC_WRAPPERS


# AC_F77_WRAPPERS
# ---------------
AC_DEFUN([AC_F77_WRAPPERS],
[AC_REQUIRE([_AC_F77_NAME_MANGLING])dnl
AC_LANG_PUSH(Fortran 77)dnl
_AC_FC_WRAPPERS
AC_LANG_POP(Fortran 77)dnl
])# AC_F77_WRAPPERS


# AC_FC_WRAPPERS
# --------------
AC_DEFUN([AC_FC_WRAPPERS],
[AC_REQUIRE([_AC_FC_NAME_MANGLING])dnl
AC_LANG_PUSH(Fortran)dnl
_AC_FC_WRAPPERS
AC_LANG_POP(Fortran)dnl
])# AC_FC_WRAPPERS


# _AC_FC_FUNC(NAME, [SHELLVAR = NAME])
# ------------------------------------
# For a Fortran subroutine of given NAME, define a shell variable
# $SHELLVAR to the Fortran-mangled name.  If the SHELLVAR
# argument is not supplied, it defaults to NAME.
AC_DEFUN([_AC_FC_FUNC],
[_AC_FORTRAN_ASSERT()dnl
case $ac_cv_[]_AC_LANG_ABBREV[]_mangling in
  upper*) ac_val="m4_toupper([$1])" ;;
  lower*) ac_val="m4_tolower([$1])" ;;
  *)      ac_val="unknown" ;;
esac
case $ac_cv_[]_AC_LANG_ABBREV[]_mangling in *," underscore"*) ac_val="$ac_val"_ ;; esac
m4_if(m4_index([$1],[_]),-1,[],
[case $ac_cv_[]_AC_LANG_ABBREV[]_mangling in *," extra underscore"*) ac_val="$ac_val"_ ;; esac
])
m4_default([$2],[$1])="$ac_val"
])# _AC_FC_FUNC


# AC_F77_FUNC(NAME, [SHELLVAR = NAME])
# ------------------------------------
AC_DEFUN([AC_F77_FUNC],
[AC_REQUIRE([_AC_F77_NAME_MANGLING])dnl
AC_LANG_PUSH(Fortran 77)dnl
_AC_FC_FUNC([$1],[$2])
AC_LANG_POP(Fortran 77)dnl
])# AC_F77_FUNC


# AC_FC_FUNC(NAME, [SHELLVAR = NAME])
# -----------------------------------
AC_DEFUN([AC_FC_FUNC],
[AC_REQUIRE([_AC_FC_NAME_MANGLING])dnl
AC_LANG_PUSH(Fortran)dnl
_AC_FC_FUNC([$1],[$2])
AC_LANG_POP(Fortran)dnl
])# AC_FC_FUNC


# AC_FC_SRCEXT(EXT, [ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
# -----------------------------------------------------------
# Set the source-code extension used in Fortran (FC) tests to EXT (which
# defaults to f).  Also, look for any necessary additional FCFLAGS needed
# to allow this extension, and store them in the output variable
# FCFLAGS_<EXT> (e.g. FCFLAGS_f90 for EXT=f90).  If successful,
# call ACTION-IF-SUCCESS.  If unable to compile source code with EXT,
# call ACTION-IF-FAILURE, which defaults to failing with an error
# message.
#
# (The flags for the current source-code extension, if any, are stored
# in the FCFLAGS_SRCEXT variable and are automatically used in subsequent
# autoconf tests.)
#
# For ordinary extensions like f90, etcetera, the modified FCFLAGS
# are currently needed for IBM's xlf* and Intel's ifc (grrr).  Unfortunately,
# xlf* will only take flags to recognize one extension at a time, so if the
# user wants to compile multiple extensions (.f90 and .f95, say), she
# will need to use the FCFLAGS_F90 and FCFLAGS_F95 individually rather
# than just adding them all to FCFLAGS, for example.
#
# Also, for Intel's ifc compiler (which does not accept .f95 by default in
# some versions), the $FCFLAGS_<EXT> variable *must* go immediately before
# the source file on the command line, unlike other $FCFLAGS.  Ugh.
AC_DEFUN([AC_FC_SRCEXT],
[AC_LANG_ASSERT(Fortran)dnl
AC_MSG_WARN([AC@&t@_FC@&t@_SRCEXT is deprecated. Use AC@&t@_FC_FIXEDFORM([srcext]) or AC@&t@_FC_FREEFORM([srcext]) as appropriate.])
AC_CACHE_CHECK([for Fortran flag to compile .$1 files],
                ac_cv_fc_srcext_$1,
[ac_ext=$1
ac_fc_srcext_FCFLAGS_SRCEXT_save=$FCFLAGS_SRCEXT
FCFLAGS_SRCEXT=""
ac_cv_fc_srcext_$1=unknown
for ac_flag in none -qsuffix=f=$1 -Tf; do
  test "x$ac_flag" != xnone && FCFLAGS_SRCEXT="$ac_flag"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM()], [ac_cv_fc_srcext_$1=$ac_flag; break])
done
rm -f conftest.$ac_objext conftest.$1
FCFLAGS_SRCEXT=$ac_fc_srcext_FCFLAGS_SRCEXT_save
])
if test "x$ac_cv_fc_srcext_$1" = xunknown; then
  m4_default([$3],[AC_MSG_ERROR([Fortran could not compile .$1 files])])
else
  FC_SRCEXT=$1
  if test "x$ac_cv_fc_srcext_$1" = xnone; then
    FCFLAGS_SRCEXT=""
    FCFLAGS_[]$1[]="$FCFLAGS_[]$1[]"
  else
    FCFLAGS_SRCEXT=$ac_cv_fc_srcext_$1
    FCFLAGS_[]$1[]="$FCFLAGS_[]$1[] $ac_cv_fc_srcext_$1"
  fi
  AC_SUBST(FCFLAGS_[]$1)
  $2
fi
])# AC_FC_SRCEXT


dnl # AC_FPP_SRCEXT(EXT, [ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
dnl # -----------------------------------------------------------
dnl # Set the source-code extension used in Fortran (FPP) tests to EXT (which
dnl # defaults to f).  Also, look for any necessary additional FPPFLAGS needed
dnl # to allow this extension, and store them in the output variable
dnl # FPPFLAGS_<EXT> (e.g. FPPFLAGS_f90 for EXT=f90).  If successful,
dnl # call ACTION-IF-SUCCESS.  If unable to compile source code with EXT,
dnl # call ACTION-IF-FAILURE, which defaults to failing with an error
dnl # message.
dnl #
dnl # (The flags for the current source-code extension, if any, are stored
dnl # in the FPPFLAGS_SRCEXT variable and are automatically used in subsequent
dnl # autoconf tests.)
dnl #
dnl # For ordinary extensions like f90, etcetera, the modified FPPFLAGS
dnl # are currently needed for IBM's xlf* and Intel's ifc (grrr).  Unfortunately,
dnl # xlf* will only take flags to recognize one extension at a time, so if the
dnl # user wants to compile multiple extensions (.f90 and .f95, say), she
dnl # will need to use the FPPFLAGS_F90 and FPPFLAGS_F95 individually rather
dnl # than just adding them all to FPPFLAGS, for example.
dnl #
dnl # Also, for Intel's ifc compiler (which does not accept .f95 by default in
dnl # some versions), the $FPPFLAGS_<EXT> variable *must* go immediately before
dnl # the source file on the command line, unlike other $FPPFLAGS.  Ugh.
dnl #
dnl #
dnl # The problem with this is ... we are doing this in order to find what
dnl # flag must be passed to the compiler to compile a preprocessed file;
dnl # but we haven't (or need not have) constructed the ac_fpp_compile 
dnl # variable yet. Thus we must check to see first if we know how to compile
dnl # preprocessed source, and otherwise pretend that we're compiling normal
dnl # Fortran.
dnl # 
dnl # The (other) problem is that flags may be different for preprocessed
dnl # free form, and preprocessed fixed form, source. This is true for
dnl # the xlf compiler, but this can be faked by treating the free/fixed
dnl # form question and the fpp/nofpp question separately, and adding both
dnl # to F(C|PP)FLAGS_SRCEXT. This requires good behaviour on the part of
dnl # automake or the package maintainer though. It really ought to be 
dnl # fixed at some point.
dnl #
dnl # For the moment, the flags we test for are:
dnl #       <none> : rely on current source extension
dnl # -qsuffix=cpp : xlf 
dnl #         -cpp : Tru64, (some versions of) Intel Fortran
dnl #         -fpp : NAG, (some versions of) Intel Fortran
dnl #  -lfe "-Cpp" : Lahey (undocumented)
dnl AC_DEFUN([AC_FPP_SRCEXT],
dnl [dnl
dnl AC_CACHE_CHECK([for Fortran flag to compile preprocessable .$1 files],
dnl                 ac_cv_fpp_srcext_$1,
dnl [ac_ext=$1
dnl ac_fpp_srcext_FPPFLAGS_SRCEXT_save=$FPPFLAGS_SRCEXT
dnl FPPFLAGS_SRCEXT=""
dnl ac_cv_fpp_srcext_$1=unknown
dnl for ac_flag in none -qsuffix=cpp=$1 -cpp -fpp '-lfe "-Cpp"'; do
dnl  test "x$ac_flag" != xnone && FPPFLAGS_SRCEXT="$ac_flag"
dnl   if test "x$ac_cv_fpp_build_rule" = x; then
dnl      # We're probably being called from AC_FPP_PROG. We don't know
dnl      # how to compile fpp programs yet.
dnl      AC_LANG_PUSH([Fortran])
dnl      AC_COMPILE_IFELSE([AC_LANG_PROGRAM()], [ac_cv_fpp_srcext_$1=$ac_flag; break])
dnl      AC_LANG_POP([Fortran])
dnl   else
dnl      # We already know what ac_fpp_compile etc. are:
dnl      AC_LANG_PUSH([Preprocessed Fortran])
dnl      ac_ext=$1
dnl      AC_COMPILE_IFELSE([
dnl       program main
dnl #if 1
dnl       end
dnl #endif
dnl ],   [ac_cv_fpp_srcext_$1=$ac_flag; break])
dnl      AC_LANG_POP([Preprocessed Fortran]) 
dnl   fi
dnl done
dnl rm -f conftest.$ac_objext conftest.$1
dnl FPPFLAGS_SRCEXT=$ac_fpp_srcext_FPPFLAGS_SRCEXT_save
dnl ])
dnl if test "x$ac_cv_fpp_srcext_$1" = xunknown; then
dnl   m4_default([$3],[AC_MSG_ERROR([Fortran could not compile .$1 files])])
dnl else
dnl   FPP_SRCEXT=$1
dnl   if test "x$ac_cv_fpp_srcext_$1" = xnone; then
dnl     FPPFLAGS_SRCEXT=""
dnl     FPPFLAGS_[]$1[]="$FPPFLAGS_[]$1[]"
dnl   else
dnl     FPPFLAGS_SRCEXT=$ac_cv_fpp_srcext_[]$1
dnl     FPPFLAGS_[]$1[]="$FPPFLAGS_[]$1[] $ac_cv_fpp_srcext_[]$1"
dnl   fi
dnl AC_SUBST(FPPFLAGS_[]$1)
dnl   $2
dnl fi
dnl ])# AC_FPP_SRCEXT

# -------------- #
# Utility macros #
# -------------- #

# AC_FC_FIXEDFORM([SRCEXT], [ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
# Look for compiler flags to make the Fortran (FC) compiler accept
# fixed-format source code, with a source extension of SRCEXT, 
# and puts any necessary flags in FCFLAGS_fixed_<SRCEXT>.  Call
# ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile fixed-format code using new extension) and ACTION-IF-FAILURE (defaults 
# to failing with an error message) if not.  
#
# The known flags are:
#                -FI: Intel compiler (icc, ecc)
#            -qfixed: IBM compiler (xlf)
#             -fixed: NAG compiler
#           -Mnofree: PGI compiler
# We try to test the "more popular" flags first, by some prejudiced
# notion of popularity.
AC_DEFUN([AC_FC_FIXEDFORM],
[AC_REQUIRE([AC_PROG_FC])
AC_CACHE_CHECK([for Fortran flag needed to allow fixed-form source for .$1 suffix],
                ac_cv_fc_fixedform_$1,
[ac_cv_fc_fixedform_$1=unknown
ac_ext=$1
ac_fc_fixedform_FCFLAGS_save=$FCFLAGS
for ac_flag in none -FI -qfixed -fixed --fix -Mnofree
do
  test "x$ac_flag" != xnone && FCFLAGS="$ac_fc_fixedform_FCFLAGS_save $ac_flag"
  AC_COMPILE_IFELSE([
      PROGRAM FIXEDFORM
C THIS COMMENT SHOULD CONFUSE FREEFORM COMPILERS
      PRI  NT*, 'HELLO '//
     .      'WORLD.'
      ENDP ROGRAM
],
                    [ac_cv_fc_fixedform_$1=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ac_fc_fixedform_FCFLAGS_save
])
if test "x$ac_cv_fc_fixedform_$1" = xunknown; then
  m4_default([$3],
             [AC_MSG_ERROR([Cannot compile fixed-form source with .$1 suffix], 77)])
else
  if test "x$ac_cv_fc_fixedform_$1" != xnone; then
    FCFLAGS_fixed_[]$1[]="$ac_cv_fc_fixedform_$1"
  fi
  $2
fi
])# AC_FC_FIXEDFORM


# AC_FC_FREEFORM([SRCEXT], [ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
# ------------------------------------------------------------------
# Look for compiler flags to make the Fortran (FC) compiler accept
# free-format source code, with a source extension of SRCEXT, 
# and puts any necessary flags in FCFLAGS_<SRCEXT>.  Call
# ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile fixed-format code using new extension) and ACTION-IF-FAILURE (defaults 
# to failing with an error message) if not.
#
# For backwards compatibility, this macro may be called without 
# specifying SRCEXT, in which case, a default extension of f90
# is used. This usage is deprecated.
#
# The known flags are:
#                        -ffree-form: GNU g77
#                                -FR: Intel compiler (icc, ecc)
#                              -free: Compaq compiler (fort), NAG compiler 
#         -qfree -qsuffix=f=<SRCEXT>: IBM compiler (xlf) (generates a warning with
#                                         recent versions)
#     -qfree=f90 -qsuffix=f=<SRCEXT>: Newer xlf versions 
#                             --nfix: Lahey compiler
#                 -Mfree, -Mfreeform: Portland Group compiler
#                          -freeform: SGI compiler
#                            -f free: Absoft Fortran
# We try to test the "more popular" flags first, by some prejudiced
# notion of popularity.
AC_DEFUN([AC_FC_FREEFORM],
[AC_REQUIRE([AC_PROG_FC])
AC_CACHE_CHECK([for Fortran flag needed to allow free-form source for .$1 suffix],
                ac_cv_fc_freeform_$1,
[ac_cv_fc_freeform_$1=unknown
ac_ext=$1
ac_fc_freeform_FCFLAGS_save=$FCFLAGS
for ac_flag in none -ffree-form -FR -free "-qfree=f90" "-qfree=f90 -qsuffix=f=$ac_ext"\
               -qfree "-qfree -qsuffix=f=$ac_ext" -Mfree -Mfreeform \
               -freeform "-f free" --nfix
do
  test "x$ac_flag" != xnone && FCFLAGS="$ac_fc_freeform_FCFLAGS_save $ac_flag"
  AC_COMPILE_IFELSE([
program freeform
! FIXME: how to best confuse non-freeform compilers?
print *, 'Hello ', &
'world.'
end program],
                    [ac_cv_fc_freeform_$1=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ac_fc_freeform_FCFLAGS_save
])
if test "x$ac_cv_fc_freeform_$1" = xunknown; then
  m4_default([$3],
             [AC_MSG_ERROR([Cannot compile free-form source with .$1 suffix], 77)])
else
  if test "x$ac_cv_fc_freeform_$1" != xnone; then
    FCFLAGS_free_[]$1[]="$ac_cv_fc_freeform_$1"
  fi
  $2
fi
])# AC_FC_FREEFORM


# _AC_FPP_FIXEDFORM_F
# Related to AC_FPP_FIXEDFORM, but used only from _AC_PROG_FC_FPP.
# How do we directly compile a preprocessable .F file?
# This should be a no-op on all systems except those with case-sensitive
# filenames, and those which can't do direct compilation anyway.
# Do not put result into cached variable if it fails.
AC_DEFUN([_AC_FPP_FIXEDFORM_F],[
ac_ext=F
ac_fpp_fixedform_FCFLAGS_save=$FCFLAGS
for ac_flag in none "-x f77-cpp-input" "-FI -cpp" "-qfixed -qsuffix=cpp=F" "-fixed -fpp" "-lfe \"-Cpp\" --fix"
do
  test "x$ac_flag" != xnone && FCFLAGS="$ac_fpp_fixedform_FCFLAGS_save $ac_flag"
    AC_COMPILE_IFELSE([
      PROGRAM FIXEDFORM
C THIS COMMENT SHOULD CONFUSE FREEFORM COMPILERS
      PRI  NT*, 'HELLO '//
     .      'WORLD.'
#ifndef OK
      ENDP ROGRAM
#endif
  ],
  [ac_cv_fpp_fixedform_F=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ac_fpp_fixedform_FCFLAGS_save
if test "x$ac_cv_fpp_fixedform_F" = x; then
  AC_MSG_WARN([Cannot compile fixed-form preprocessable Fortran with a .F extension.])
else
  if test "$ac_cv_fpp_fixedform_F" != none; then
    FPPFLAGS_fixed_F="$ac_cv_fpp_fixedform_F"
  fi
fi
])


# AC_FPP_FIXEDFORM([SRCEXT], [ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
# Look for compiler flags to make the Fortran (FC) compiler accept
# preprocessed fixed-format source code, with a source extension of SRCEXT, 
# and puts any necessary flags in FPPFLAGS_fixed_<SRCEXT>.  Call
# ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile fixed-format code using new extension) and ACTION-IF-FAILURE (defaults 
# to failing with an error message) if not.  
#
# Only applicable when using direct compilation.
#
# The known flags are:
#        -x f77-cpp-input: g77
#                -FI -cpp: Intel compiler (ifort)
# -qfixed -qsuffix=cpp=$1: IBM compiler (xlf)
#             -fixed -fpp: NAG compiler
#       -lfe "-Cpp" --fix: Lahey compiler
#                -Mnofree: PGI (no flag for preprocessing available)
# We try to test the "more popular" flags first, by some prejudiced
# notion of popularity.
# NB when updfating this list of flags, also update those of the previous
# macro.
AC_DEFUN([AC_FPP_FIXEDFORM],
[AC_REQUIRE([AC_PROG_FPP])
AC_CACHE_CHECK([for Fortran flag needed to allow preprocessed fixed-form source for .$1 suffix],
                ac_cv_fpp_fixedform_$1,
[
if test $ac_cv_fpp_build_rule = direct; then
  ac_cv_fpp_fixedform_$1=unknown
  ac_ext=$1
  ac_fpp_fixedform_FCFLAGS_save=$FCFLAGS
  for ac_flag in none "-x f77-cpp-input" "-FI -cpp" "-qfixed -qsuffix=cpp=$1" "-fixed -fpp" "-lfe \"-Cpp\" --fix"
  do
    test "x$ac_flag" != xnone && FCFLAGS="$ac_fpp_fixedform_FCFLAGS_save $ac_flag"
    AC_COMPILE_IFELSE([
      PROGRAM FIXEDFORM
C THIS COMMENT SHOULD CONFUSE FREEFORM COMPILERS
      PRI  NT*, 'HELLO '//
     .      'WORLD.'
#ifndef OK
      ENDP ROGRAM
#endif
  ],
                    [ac_cv_fpp_fixedform_$1=$ac_flag; break])
  done
  rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
  FCFLAGS=$ac_fpp_fixedform_FCFLAGS_save
else
  ac_cv_fpp_fixedform_$1=none
fi # test $ac_cv_fpp_build_rule = direct
])
if test "x$ac_cv_fpp_fixedform_$1" = xunknown; then
  m4_default([$3],
             [AC_MSG_ERROR([Cannot compile fixed-form source with .$1 suffix], 77)])
else
  if test "x$ac_cv_fpp_fixedform_$1" != xnone; then
    FPPFLAGS_fixed_[]$1[]="$ac_cv_fpp_fixedform_$1"
  fi
  $2
fi
if test $ac_cv_fpp_build_rule = direct; then
    FPP_FIXED_COMPILE_EXT=$1
    FPP_FIXED_PREPROCESS_EXT=DUMMY$1
dnl    FPP_FIXED_MAKE_FLAGS='$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS)'
else
    FPP_FIXED_COMPILE_EXT=DUMMYf
    FPP_FIXED_PREPROCESS_EXT=$1
    FPP_FIXED_MAKE_FLAGS=""
fi
])# AC_FPP_FIXEDFORM


# AC_FPP_FREEFORM([SRCEXT], [ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
# ------------------------------------------------------------------
# Look for compiler flags to make the Fortran (FC) compiler accept
# preprocessed free-format source code, with a source extension of SRCEXT, 
# and puts any necessary flags in FPPFLAGS_<SRCEXT>.  Call
# ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile fixed-format code using new extension) and ACTION-IF-FAILURE (defaults 
# to failing with an error message) if not.  
#
# Only applicable when using direct compilation.
#
# The known flags are:
#     -ffree-form -x f77-cpp-input: GNU g77
#                         -FR -cpp: Intel compiler (ifort)
#                       -free -cpp: Compaq compiler (fort), NAG compiler 
#     -qfree -qsuffix=cpp=<SRCEXT>: IBM compiler (xlf) (generates a warning with
#                                         recent versions)
# -qfree=f90 -qsuffix=cpp=<SRCEXT>: Newer xlf versions 
#                --nfix -lfe="Cpp": Lahey compiler
#               -Mfree, -Mfreeform: PGI (no flag for preprocessing available)
#          -freeform: SGI compiler
#            -f free: Absoft Fortran
# We try to test the "more popular" flags first, by some prejudiced
# notion of popularity.
AC_DEFUN([AC_FPP_FREEFORM],
[AC_REQUIRE([AC_PROG_FPP])
AC_CACHE_CHECK([for Fortran flag needed to allow free-form source for .$1 suffix],
                ac_cv_fpp_freeform_$1,
[
if test $ac_cv_fpp_build_rule = direct; then
  ac_cv_fpp_freeform_$1=unknown
  ac_ext=$1
  ac_fpp_freeform_FCFLAGS_save=$FCFLAGS
  for ac_flag in none -ffree-form -FR -free "-qfree=f90" "-qfree=f90 -qsuffix=cpp=$1"\
                 -qfree "-qfree -qsuffix=cpp=$1" -Mfree -Mfreeform \
                 -freeform "-f free" --nfix
  do
    test "x$ac_flag" != xnone && FCFLAGS="$ac_fpp_freeform_FCFLAGS_save $ac_flag"
    AC_COMPILE_IFELSE([
program freeform
! FIXME: how to best confuse non-freeform compilers?
print *, 'Hello ', &
'world.'
#ifndef OK
end program
#endif
  ],
                    [ac_cv_fpp_freeform_$1=$ac_flag; break])
  done
  rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
  FCFLAGS=$ac_fpp_freeform_FCFLAGS_save
else
  ac_cv_fpp_freeform_$1=none
fi # test $ac_cv_fpp_build_rule = direct 
])
if test "x$ac_cv_fpp_freeform_$1" = xunknown; then
  m4_default([$3],
             [AC_MSG_ERROR([Cannot compile free-form source with .$1 suffix], 77)])
else
  if test "x$ac_cv_fpp_freeform_$1" != xnone; then
    FPPFLAGS_free_[]$1[]="$ac_cv_fpp_freeform_$1"
  fi
  $2
fi
if test $ac_cv_fpp_build_rule = direct; then
    FPP_FREE_COMPILE_EXT=$1
    FPP_FREE_PREPROCESS_EXT=DUMMY$1
dnl     FPP_FREE_MAKE_FLAGS='$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS)'
else
    FPP_FREE_COMPILE_EXT=DUMMYf
    FPP_FREE_PREPROCESS_EXT=$1
    FPP_FREE_MAKE_FLAGS=""
fi
])# AC_FPP_FREEFORM


# AC_FC_OPEN_SPECIFIERS(specifier ...)
# ------------------------------------
#
# The Fortran OPEN statement is a rich source of portability problems,
# since there are numerous common extensions which consiste of extra
# specifiers, several of which are useful when they are available.
# For each of the specifiers in the (whitespace-separated) argument
# list, define HAVE_FC_OPEN_mungedspecifier if the specifier may be
# given as  argument to the OPEN statement.  The `mungedspecifier' is the
# `specifier' converted to uppercase and with all characters outside
# [a-zA-Z0-9_] deleted.  Note that this may include `specifiers' such
# as "access='append'" and "[access='sequential',recl=1]" (note quoting of
# comma) to check combinations of specifiers.  You may not include a
# space in the `specifier', even quoted.  Each argument must be a
# maximum of 65 characters in length (to abide by Fortran 77
# line-length limits). 
#
dnl Multiple m4_quote instances are necessary in case specifier includes comma.
dnl In the Fortran OPEN line, include status='scratch' unless status=???
dnl is in the specifier being tested.
dnl Put specifier on continuation line, in case it's long.
AC_DEFUN([AC_FC_OPEN_SPECIFIERS],
         [AC_REQUIRE([AC_PROG_FC])dnl
          AC_LANG_PUSH([Fortran])
          AC_FOREACH([Specifier],
                     m4_quote(m4_toupper([$1])),
                     [m4_define([mungedspec],
                                m4_bpatsubst(m4_quote(Specifier), [[^a-zA-Z0-9_]], []))
                      AC_CACHE_CHECK([whether ${FC} supports OPEN specifier ]m4_quote(Specifier),
                          [ac_cv_fc_spec_]mungedspec,
                          [AC_COMPILE_IFELSE(AC_LANG_PROGRAM([],
                                                             [      OPEN(UNIT=99,]m4_if(m4_bregexp(m4_quote(Specifier), [\<STATUS *=]), -1, [STATUS='SCRATCH'[,]], [])
     :m4_quote(Specifier)[)]),
                                             [ac_cv_fc_spec_]mungedspec=yes,
                                             [ac_cv_fc_spec_]mungedspec=no)])
                      if test $ac_cv_fc_spec_[]mungedspec = yes; then
                          AC_DEFINE([HAVE_FC_OPEN_]mungedspec, 1,
                                    [Define to 1 if the Fortran compiler supports OPEN specifier ]m4_quote(Specifier))
                      fi])
           AC_LANG_POP([Fortran])
])# AC_FC_OPEN_SPECIFIERS


# AC_FC_CHECK_INTRINSICS(intrinsic-function ...)
# ----------------------------------------------
#
# Like AC_CHECK_FUNCS, but instead determine the intrinsics available
# to the Fortran compiler.  For each intrinsic in the
# (whitespace-separated and case-insensitive) argument list, define
# HAVE_INTRINSIC_intrinsic-function (uppercased) if it is available.
# For example, AC_FC_CHECK_INTRINSICS(sin) would define
# HAVE_INTRINSIC_SIN if the `sin' intrinsic function were available
# (there are probably rather few Fortrans which don't have this
# function).  The macro works for both intrinsic functions and
# intrinsic subroutines.
AC_DEFUN([AC_FC_CHECK_INTRINSICS],
  [AC_REQUIRE([AC_PROG_FC])dnl
   AC_LANG_PUSH([Fortran])
   AC_FOREACH([IntrinsicName],
              dnl In case the user is mad, escape impossible names
              m4_bpatsubst(m4_toupper([$1]), [[^a-zA-Z0-9_ ]], [_]),
              [AC_CACHE_CHECK([whether ${FC} supports intrinsic ]IntrinsicName,
                              [ac_cv_fc_has_]IntrinsicName,
                              [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
dnl All we need do is attempt to declare the thing as an intrinsic, 
dnl since it is an error to thus declare a symbol which is not
dnl in fact an intrinsic function.
      intrinsic IntrinsicName
])],
                               ac_cv_fc_has_[]IntrinsicName=yes,
                               ac_cv_fc_has_[]IntrinsicName=no)])
               if test $ac_cv_fc_has_[]IntrinsicName = yes; then
                   AC_DEFINE([HAVE_INTRINSIC_]IntrinsicName, 1,
                             [Define to 1 if the Fortran compiler supports intrinsic ]IntrinsicName)
               fi
              ])
   AC_LANG_POP([Fortran])
])# AC_FC_CHECK_INTRINSICS


# AC_FC_RECL_UNIT
# ----------------
#
# When opening a file for direct access, you must specify
# the record length with the @samp{OPEN} specifier @samp{RECL};
# however in the case of unformatted direct access files, the
# @emph{units} of this specifier are processor dependent, and may be
# words or bytes.  This macro determines the units and defines
# @samp{FC_RECL_UNIT} to contain the number of bytes (1, 2, 4, 8, ...) in
# the processor's unit of measurement.
#
# Note that unformatted files are not themselves portable, and should
# only be used as either temporary files, or as data files which will
# be read by a program or library compiled with the same Fortran
# processor.  With this macro, however, you can read and write such
# files in a portable way.
AC_DEFUN([AC_FC_RECL_UNIT],
         [AC_REQUIRE([AC_PROG_FC])dnl
          AC_CACHE_CHECK([units for Fortran OPEN RECL],
                         [ac_cv_fc_recl_unit],
                         [AC_LANG_PUSH([Fortran])
                          AC_RUN_IFELSE([AC_LANG_SOURCE([dnl
      PROGRAM TESTRECL
      IMPLICIT NONE
      INTEGER NBYTES
*   Make sure these values agree
      PARAMETER ( NBYTES = 8 )
      REAL * 8 TOFILE, FROMFILE

      INTEGER RECLEN,UNITLEN,OUTUNIT

      TOFILE = 123456789D56
      OUTUNIT = 10

*   Record length to try
      RECLEN = 1
*   Unitlen is the result -- zero indicates that no value was successful
      UNITLEN = 0

*     Keep on increasing the record length until we hit a
*     size that allows us to write a number and read it back correctly.
      DO WHILE (RECLEN .LE. 8)

         OPEN(UNIT = OUTUNIT,
     :        FILE = 'conftest.rcl1',
     :        STATUS = 'REPLACE',
     :        FORM = 'UNFORMATTED',
     :        ACCESS = 'DIRECT',
     :        RECL = RECLEN,
     :        ERR = 101)

*      Write two records to the output file, so that the second will stomp
*      on the end of the first if the record length is too short.
         WRITE(UNIT=OUTUNIT,REC=1,ERR=101) TOFILE
         WRITE(UNIT=OUTUNIT,REC=2,ERR=101) TOFILE
         READ(UNIT=OUTUNIT,REC=1,ERR=101) FROMFILE
         IF (TOFILE .EQ. FROMFILE) THEN
            UNITLEN = NBYTES/RECLEN
            GOTO 102
         END IF

*      Error opening unit
 101     CONTINUE

         CLOSE(UNIT=OUTUNIT)

         RECLEN = RECLEN * 2
      END DO

*   Got a match
 102  CONTINUE

      OPEN(UNIT = OUTUNIT,
     :     FILE = 'conftest.rcl2',
     :     STATUS = 'REPLACE',
     :     ERR = 103)
      WRITE(OUTUNIT,'(I3)') UNITLEN
      CLOSE(UNIT = OUTUNIT)
 103  CONTINUE

      END
])],
                                        [],
                                        AC_MSG_FAILURE([Can't test for RECL length]),
                                        AC_MSG_FAILURE([Can't cross-compile: can't test for RECL length]))
                          AC_LANG_POP([Fortran])
                          if test -r conftest.rcl2; then
                              ac_cv_fc_recl_unit=`cat conftest.rcl2`
                          else
                              ac_cv_fc_recl_unit=0
                          fi
                          rm -f conftest.*])
          if test -n "$ac_cv_fc_recl_unit" -a $ac_cv_fc_recl_unit -gt 0; then
              AC_DEFINE_UNQUOTED([FC_RECL_UNIT], $ac_cv_fc_recl_unit,
                        [Define to the length in bytes of the unit that OPEN RECL expects])
          fi
])# AC_FC_RECL_UNIT


# AC_FC_MOD_PATH_FLAG
# -------------------------
# Check which flag is necessary to alter the compiler's search path
# for module files.
# This obviously requires that the compiler has some notion of 
# module files as separate from object files and some sensible 
# method of altering its search path. This will therefore not work 
# on early Cray F90 compilers, or on v5 (and 6?) of ifc.
#
# Nearly every compiler I have found uses -Ipath for this purpose;
# Sun F95 v7.1 (at least), uses -Mpath
# 
AC_DEFUN([AC_FC_CHECK_MOD_PATH_FLAG],[
          AC_REQUIRE([AC_PROG_FC])
          ac_cv_fc_mod_path_flag=no
          AC_MSG_CHECKING([for flag to alter module search path])
	  mkdir conftestdir
          cd conftestdir
          cat > conftest.$ac_ext << \_ACEOF
      module conftest
       implicit none
       integer :: i
      end module conftest
_ACEOF
          _AC_EVAL_STDERR($ac_compile)
          cd ..
          for i in -I -M; do
            if test $ac_cv_fc_mod_path_flag == "no"; then
               FCFLAGS_save=$FCFLAGS
               FCFLAGS="$FCFLAGS ${i}conftestdir"
               AC_COMPILE_IFELSE([
      subroutine test
       use conftest
       implicit none
       i = 0
      end subroutine test
],
              [FC_MOD_FLAG=$i; ac_cv_fc_mod_path_flag=$i],
              [:])
            fi
            FCFLAGS=$FCFLAGS_save
          done
          AC_MSG_RESULT([$ac_cv_fc_mod_path_flag])
          rm -rf conftestdir
          AS_IF([test $ac_cv_fc_mod_path_flag != "no"],
                [$1],
                [m4_default([$2],[AC_MSG_ERROR([Cannot find flag to alter module search path])])])
	  AC_SUBST(FC_MOD_FLAG)
])# AC_FC_MOD_PATH_FLAG

# -------------------------------------- #
# Feature tests for Preprocessed Fortran #
# -------------------------------------- #


# ------------------------------------------#
# Some test programs for different features #
# ------------------------------------------#

# _AC_LANG_PROGRAM_FPP_SIMPLE
# ---------------------------
# The minimum test program - any compiler supporting
# preprocessing should handle this
AC_DEFUN([_AC_LANG_PROGRAM_FPP_SIMPLE],
         [AC_LANG_PROGRAM([@%:@define OK], [dnl
#ifndef OK
      syntax error
#endif
])])#_AC_LANG_PROGRAM_FPP_SIMPLE


# _AC_LANG_PROGRAM_FPP_ONLY
# ---------------------------
# Test program for pure preprocessing
AC_DEFUN([_AC_LANG_PROGRAM_FPP_ONLY],
         [AC_LANG_PROGRAM([@%:@define OK], [dnl
#ifdef OK
      REAL A
#else
      syntax error
#endif
])])#_AC_LANG_PROGRAM_FPP_ONLY


# _AC_LANG_PROGRAM_FPP_D
# ---------------------------
# Like _AC_LANG_PROGRAM_FPP_SIMPLE, but OK is passed via -D switch
AC_DEFUN([_AC_LANG_PROGRAM_FPP_D],
[AC_LANG_PROGRAM([],[
#ifndef OK
      syntax error
#endif
])])#_AC_LANG_PROGRAM_FPP_D


# _AC_LANG_PROGRAM_FPP_I
# ---------------------------
# Test for #include statement
# If unsupported, this should give a type error
AC_DEFUN([_AC_LANG_PROGRAM_FPP_I],
[AC_LANG_PROGRAM([],[
      IMPLICIT CHARACTER (c)
!     Comments in test programs should be freeform compliant just in case.
!     conftest.inc contains the Fortran statement "REAL cc"
#include "conftest.inc"
      cc=1.
])])#_AC_LANG_PROGRAM_FPP_I


# _AC_LANG_PROGRAM_FPP_SUBS
# ---------------------------
# Test whether cpp symbols are expanded in Fortran code lines
# If not, this should give a type error
AC_DEFUN([_AC_LANG_PROGRAM_FPP_SUBS],
[AC_LANG_PROGRAM([
#define NM xxxx
], [      IMPLICIT CHARACTER (n)
      REAL xxxx
      NM=1.
])])#_AC_LANG_PROGRAM_FPP_SUBS


# _AC_LANG_PROGRAM_FPP_WRAP
# ---------------------------
# Test whether preprocessor breaks lines that become too long due
# to macro substitution. 
# If not, this gives an "unterminated character constant" error
AC_DEFUN([_AC_LANG_PROGRAM_FPP_WRAP],
[AC_LANG_PROGRAM([
#define LONG '901234567890123456789012345678901234567890123456789012345678901234567890'
],[      CHARACTER*80 A
      A=LONG
])])#_AC_LANG_PROGRAM_FPP_WRAP


# _AC_LANG_PROGRAM_FPP_CSTYLE
# ---------------------------
# Test program for C style comments
AC_DEFUN([_AC_LANG_PROGRAM_FPP_CSTYLE],
[AC_LANG_PROGRAM([],[
      A=1. /* C-style comment */
])])#_AC_LANG_PROGRAM_FPP_CSTYLE

# _AC_LANG_PROGRAM_FPP_CXXSTYLE
# ---------------------------
# Test program for C++ style comments
AC_DEFUN([_AC_LANG_PROGRAM_FPP_CXXSTYLE],
[
      PROGRAM MAIN
      CHARACTER*10 C
      C = "abcde" // "fghij"; END PROGRAM
]
)#_AC_LANG_PROGRAM_FPP_CXXSTYLE

# ----------------#
# Internal macros #
# ----------------#


# _AC_PROG_FPP_FEATURES ([feature list])
# --------------------------------------
# Parse the feature list from configure.in
AC_DEFUN([_AC_PROG_FPP_FEATURES],
[# defaults for needed features
ac_fpp_need_d=yes
ac_fpp_need_i=yes
ac_fpp_need_subs=yes
ac_fpp_need_wrap=no
ac_fpp_need_cstyle=no
ac_fpp_need_CSTYLE=no
ac_fpp_need_cxxstyle=yes
ac_fpp_need_CXXSTYLE=no
for _t in $1 nil
do
    case $_t in
        d*)    ac_fpp_need_d=yes    ;;
        nod*)  ac_fpp_need_d=no     ;;
        i*)    ac_fpp_need_i=yes    ;;   
        noi*)  ac_fpp_need_i=no     ;;
        s*)    ac_fpp_need_subs=yes    ;;   
        nos*)  ac_fpp_need_subs=no     ;;
        w*)    ac_fpp_need_wrap=yes   ;;   
        now*)  ac_fpp_need_wrap=no    ;;   
        cs*)    ac_fpp_need_cstyle=yes ;;   
        nocs*)  ac_fpp_need_cstyle=no  ;;   
        CS*)    ac_fpp_need_CSTYLE=yes ;;   
        noCS*)  ac_fpp_need_CSTYLE=no  ;;   
        cx*)    ac_fpp_need_cxxstyle=yes ;;   
        nocx*)  ac_fpp_need_cxxstyle=no  ;;   
        CX*)    ac_fpp_need_CXXSTYLE=yes ;;   
        noCX*)  ac_fpp_need_CXXSTYLE=no  ;;   
        nil)   ;;
    esac
done
# Wrapping requires substitution
test $ac_fpp_need_wrap = yes && ac_fpp_need_subs=yes
# Both CSTYLE and cstyle cannot be requested
# CSTYLE has precedence, since if it is not fulfilled,
# compile errors may arise
test $ac_fpp_need_CSTYLE = yes && ac_fpp_need_cstyle=no
dnl Similarly for cxxstyle
test $ac_fpp_need_CXXSTYLE = yes && ac_fpp_need_cxxstyle=no
])# _AC_PROG_FPP_FEATURES


# _AC_TEST_FPP ([command])
# ------------------------
# A helper macro to test correct fpp behaviour
# It sets ac_cv_prog_fpp and ac_fpp_out
AC_DEFUN([_AC_TEST_FPP],
[cat >conftest.$ac_ext << \_ACEOF
_AC_LANG_PROGRAM_FPP_ONLY
_ACEOF
ac_fpp_command=$1
if eval '$ac_fpp_command conftest.$ac_ext > conftest.log 2>/dev/null'; then
  test -f conftest.f &&
      cat conftest.f | grep '^      REAL A' >/dev/null 2>&1 &&
      cat conftest.f | grep 'syntax error' >/dev/null 2>&1  ||
    ac_cv_prog_fpp=$ac_fpp_command; ac_fpp_out=
  test -f conftest.log &&
      cat conftest.log | grep '^      REAL A' >/dev/null 2>&1 &&
      cat conftest.log | grep 'syntax error' >/dev/null 2>&1  ||
    ac_cv_prog_fpp=$ac_fpp_command; ac_fpp_out=' > conftest.f'
fi
rm -f conftest*
])# _AC_TEST_FPP

# _TEMP_AC_PROG_FPP
# -----------------
# Searches for the likely most suitable preprocessor for Fortran files.
# Tests first for the $FPP shell variable
# Then for a program called fpp
# Then for the fortran compiler (with appropriate flags)
# Then for g77 -E (in case we're using a better compiler than
# g77, but which doesn't preprocess).
# Then for a program called cpp.
# Then for gcc -E -x c
# And finally, the contents of $CPP,

# For each of these, we test if it will correctly preprocess
# a given file.

AC_DEFUN([_AC_PROG_FPP],
[AC_REQUIRE([_AC_PROG_FC_CPP])dnl
AC_CACHE_CHECK([how to preprocess Fortran files], ac_cv_prog_fpp,
[ac_cv_prog_fpp=
AC_LANG_ASSERT(Preprocessed Fortran)

# Let the user specify FPP
if test -n "$FPP"; then
  _AC_TEST_FPP([$FPP])
  if test -z "$ac_cv_prog_fpp"; then
    AC_MSG_ERROR([user-specified \$FPP ($FPP) does not work])
    FPP=
  fi
fi # test -n "$FPP"

if test -z "$ac_cv_prog_fpp"; then
  for ac_j in 'fpp' "$FC -F" "$FC -E" "g77 -E" 'cpp' '/lib/cpp' \
              '/usr/ccs/lib/cpp' "gcc -E -x c" "$CPP"; do
    _AC_TEST_FPP([$ac_j])
    test -n "$ac_cv_prog_fpp" && break;
  done
fi # test -z "$ac_cv_prog_fpp"

if test -z "$ac_cv_prog_fpp"; then
    AC_MSG_ERROR([cannot find a working Fortran preprocessor])
fi
])
AC_CACHE_CHECK([how to redirect $ac_cv_prog_fpp output],
  ac_cv_fpp_out, 
  [ac_cv_fpp_out=$ac_fpp_out])
FPP=$ac_cv_prog_fpp
ac_fpp_out=$ac_cv_fpp_out
])# TEMP_AC_PROG_FPP


# _AC_PROG_FPP
# ------------
# Try to figure out how to preprocess .fpp files for use with the selected
# Fortran compiler
#
# Must be run after _AC_PROG_FC_CPP (why?)
AC_DEFUN([_AC_PROG_FPP],
[AC_REQUIRE([_AC_PROG_FC_CPP])dnl
AC_CACHE_CHECK([how to preprocess Fortran files], ac_cv_prog_fpp,
[ac_cv_prog_fpp=
dnl AC_LANG_PUSH(Preprocessed Fortran)
AC_LANG_ASSERT(Preprocessed Fortran)

# Let the user specify FPP
if test -n "$FPP"; then
  _AC_TEST_FPP([$FPP])
  if test -z "$ac_cv_prog_fpp"; then
    AC_MSG_WARN([user-specified \$FPP ($FPP) does not work])
    FPP=
  fi
fi # test -n "$FPP"

# This next should never happen. We don't call this macro
# if we can already preprocess.
dnl if test -z "$ac_cv_prog_fpp" && test $ac_fpp_ok = yes; then
dnl # We know that the Fortran compiler can directly preprocess and compile.
dnl # All that remains is to find out how to run only the preprocessor.
dnl # We check that the compiler does not compile as well as
dnl # preprocessing, by looking for the file a.out (is this portable?  The
dnl # single-unix spec says that cc produces an a.out if no -o is given).
dnl # If we don't do this, and the compiler doesn't recognise, and
dnl # ignores, one of the options (for example g77 ignores -F
dnl # and returns without error), then the test appears to succeed.
dnl #
dnl # We only know the following methods of invocation: -F and -E
dnl # Tru64 Fortran uses -P, but puts the output in conftest.i 
dnl # It also returns success on f90 -E despite generating gibberish
dnl # What should we do with that?
dnl   for ac_j in "$FC -F" "$FC -E"; do
dnl     rm -f a.out
dnl     _AC_TEST_FPP([$ac_j])
dnl     if test -e a.out; then
dnl         rm -f a.out
dnl         ac_cv_prog_fpp=  # discard any value there
dnl     else
dnl         test -n "$ac_cv_prog_fpp" && break;
dnl     fi
dnl   done
dnl fi

# FIXME: should we bother testing for FC+flag preprocessing? 

if test -z "$ac_cv_prog_fpp"; then
# Either the Fortran compiler can't handle cpp, or doesn't have all the
# features, or can't be used for pure preprocessing.
# We must find another way for preprocessing.
# We try the "best" preprocessors first. We know that $FC can't preprocess
# by itself, but there is a small chance that F77 can be persuaded to
# preprocess, so we try that.
  for ac_j in 'fpp' "$CPP" "$CPP -x c" 'cpp' '/lib/cpp' '/usr/ccs/lib/cpp' \
              'g77 -E' '$CC -E' \
              "$FC -F" "$FC -E" "$F77 -F" "$F77 -E"; do
    _AC_TEST_FPP([$ac_j])
    test -n "$ac_cv_prog_fpp" && break;
  done
fi # test -z "$ac_cv_prog_fpp"

if test -z "$ac_cv_prog_fpp"; then
# This is only fatal if direct compilation doesn't work either
# but we're only here if direct compilation didn't work.
dnl   if test $ac_cv_prog_fc_cpp = no; then
    AC_MSG_ERROR([cannot find a working Fortran preprocessor])
dnl   else
dnl    AC_MSG_WARN([cannot find a working Fortran preprocessor])
dnl  fi
fi
dnl AC_LANG_POP()dnl
])
AC_CACHE_CHECK([how to redirect $ac_cv_prog_fpp output],
  ac_cv_fpp_out, 
  [ac_cv_fpp_out=$ac_fpp_out])
FPP=$ac_cv_prog_fpp
ac_fpp_out=$ac_cv_fpp_out
])# _AC_PROG_FPP


# _AC_PROG_FPP_P
# --------------
# Check whether we need to give FPP the -P option, to get it to
# produce code which FC can read.
AC_DEFUN([_AC_PROG_FPP_P],
[AC_CACHE_CHECK([whether $FPP needs the -P option],
ac_cv_prog_fpp_p,
[ac_cv_prog_fpp_p=unknown
dnl AC_LANG_PUSH(Preprocessed Fortran)
AC_LANG_ASSERT(Preprocessed Fortran)
# This will only be called from AC_PROG_FPP, and as such, the
# extension *will* be .F.
ac_ext=F
cat > conftest.$ac_ext << \_ACEOF
_AC_LANG_PROGRAM_FPP_ONLY
_ACEOF
dnl AC_LANG_POP(Preprocessed Fortran)dnl

AC_LANG_PUSH(Fortran)
ac_ext=F
ac_cmd='$FPP $FPPFLAGS conftest.$ac_ext '"$ac_fpp_out"
ac_link='$FC -o conftest$ac_exeext $FCFLAGS $LDFLAGS $FCFLAGS_SRCEXT conftest.f $LIBS'
if AC_TRY_EVAL(ac_cmd) &&
     AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
   ac_cv_prog_fpp_p=
else
   ac_save_FPPFLAGS=$FPPFLAGS
   FPPFLAGS="$FPPFLAGS -P"
   ac_cmd='$FPP $FPPFLAGS conftest.$ac_ext '"$ac_fpp_out"
   if AC_TRY_EVAL(ac_cmd) &&
       AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
     ac_cv_prog_fpp_p=-P
   fi
   FPPFLAGS=$ac_save_FPPFLAGS   
fi
rm -f conftest*
AC_LANG_POP(Fortran)dnl
])
if test "x$ac_cv_prog_fpp_p" = "xunknown"; then
   AC_MSG_WARN([$FPP cannot produce code that $FC compiles])
   FPP=
else
   FPPFLAGS="$FPPFLAGS $ac_cv_prog_fpp_p"
fi
])# _AC_PROG_FPP_P

# _AC_PROG_FPP_CSTYLE
# -------------------
# Check whether FPP lets C-style comments through to FC
AC_DEFUN([_AC_PROG_FPP_CSTYLE],
[AC_CACHE_CHECK([how to pass C-style comments to $FC], 
   ac_cv_prog_fpp_cstyle,
[ac_cv_prog_fpp_cstyle=unknown
dnl AC_LANG_PUSH(Preprocessed Fortran)
AC_LANG_ASSERT(Preprocessed Fortran)
cat > conftest.$ac_ext << \_ACEOF
_AC_LANG_PROGRAM_FPP_CSTYLE
_ACEOF
dnl AC_LANG_POP()dnl

AC_LANG_PUSH(Fortran)
ac_cmd='$FPP $FPPFLAGS conftest.$ac_ext '"$ac_fpp_out"
if AC_TRY_EVAL(ac_cmd) &&
   cat conftest.f | grep '[[*]]/.*[[*]]/' >/dev/null 2>&1; then
   ac_cv_prog_fpp_cstyle=
else
   ac_save_FPPFLAGS=$FPPFLAGS
   ac_name=`expr "x$FPP" : 'x\(fpp\)'` 
   if test "x$ac_name" = xfpp; then
     ac_flag="-c_com=no"
   else
     ac_flag="-C"
   fi
   FPPFLAGS="$FPPFLAGS $ac_flag"
   ac_cmd='$FPP $FPPFLAGS conftest.$ac_ext '"$ac_fpp_out"
   if AC_TRY_EVAL(ac_cmd) &&
     cat conftest.f | grep '/[[*]].*[[*]]/' >/dev/null 2>&1; then
     ac_cv_prog_fpp_cstyle=$ac_flag
   fi
   FPPFLAGS=$ac_save_FPPFLAGS   
fi
rm -f conftest*
AC_LANG_POP(Fortran)dnl
])
if test "x$ac_cv_prog_fpp_cstyle" = "xunknown"; then
  AC_MSG_WARN([cannot find a way to make $FPP pass C-style comments])
else
  FPPFLAGS="$FPPFLAGS $ac_cv_prog_fpp_cstyle"
fi
])# _AC_PROG_FPP_CSTYLE


# _AC_PROG_FC_CPP
# ----------------
# This macro checks whether the chosen preprocessing method
# has all the requested features.
#
#FIXME: this is only for fixed form code. Need a separate check for free-form.
#
# NB We are definitely using a suffix of .F in this case. If the filesystem
# is case-insensitive, we may need to force preprocessing.
# 
# Sets ac_fpp_ok to "no" if a requested feature is unavailable
#
AC_DEFUN([_AC_PROG_FC_CPP],
[ac_fpp_ok=yes
ac_prog_fc_cpp=no
ac_prog_fc_cpp_d=no
ac_prog_fc_cpp_i=no
ac_prog_fc_cpp_subs=no
ac_prog_fc_cpp_wrap=no
ac_prog_fc_cpp_CSTYLE=no
ac_prog_fc_cpp_cxxstyle=no

AC_LANG_ASSERT(Preprocessed Fortran)
AC_MSG_CHECKING([for fixed form Fortran preprocessor features])

if test $ac_fc_testing_fpp = direct; then
# On nearly all systems where direct compilation is possible, a .F file will 
# compile a preprocessable fixed-form file automatically. However, 
# case-insensitive filesystems (eg HFS+ on MacOSX) may get confused.
# Therefore, we must check for cpp flags.
  _AC_FPP_FIXEDFORM_F
  if test "x$ac_cv_fpp_fixedform_F" != x; then
    ac_prog_fc_cpp=yes
  else
    ac_fpp_ok=no
  fi

# It is possible we've failed the previous test because of a
# Tru64 bug where the compiler fails when called as 'f95' on 
# a .F file. It works when called as f90.
#FIXME: this does not protect the user's setting of FC, though
# we set it back if senesible.
 if test $ac_prog_fc_cpp = no && test $FC = f95; then
    FC=f90
      AC_LINK_IFELSE([_AC_LANG_PROGRAM_FPP_SIMPLE], 
        [ac_prog_fc_cpp=yes],
        [])
    if test $ac_prog_fc_cpp = no; then
      FC=f95
      ac_fpp_ok=no
    fi
  fi

  ac_first_save_FPPFLAGS=$FPPFLAGS
  FPPFLAGS="$FPPFLAGS $FPPFLAGS_F"
fi

# If we're doing direct compilation, we need only test if FC will
# do preprocessing. For indirect compilation, of course it will.
if test $ac_prog_fc_cpp = yes || test $ac_fc_testing_fpp = indirect; then

    if test $ac_fpp_need_d = yes; then
       ac_prog_fc_cpp_d=no
       ac_save_FPPFLAGS=$FPPFLAGS
       FPPFLAGS="$FPPFLAGS -DOK"
       AC_LINK_IFELSE([_AC_LANG_PROGRAM_FPP_D],
         [ac_prog_fc_cpp_d=yes; FPPFLAGS_DEF="-D"], 
         [:])
       FPPFLAGS=$ac_save_FPPFLAGS
       if test $ac_prog_fc_cpp_d = no; then
	  # stupid ibm compiler
          ac_save_FPPFLAGS=$FPPFLAGS
          FPPFLAGS="$FPPFLAGS -WF,-DOK"
          AC_LINK_IFELSE([_AC_LANG_PROGRAM_FPP_D],
            [ac_prog_fc_cpp_d=yes; FPPFLAGS_DEF="-WF,-D"], 
            [ac_fpp_ok=no])
          FPPFLAGS=$ac_save_FPPFLAGS
       fi
    fi
#FIXME we should probably do the AC_SUBST somewhere else.
    AC_SUBST(FPPFLAGS_DEF)

    if test $ac_fpp_need_i = yes; then
       mkdir conftst
       cat > conftst/conftest.inc << \_ACEOF
!     This statement overrides the IMPLICIT statement in the program
      REAL cc
_ACEOF
       ac_save_FPPFLAGS=$FPPFLAGS
       FPPFLAGS="$FPPFLAGS -Iconftst"
       AC_LINK_IFELSE([_AC_LANG_PROGRAM_FPP_I],
         [ac_prog_fc_cpp_i=yes],
         [ac_fpp_ok=no])
       rm -rf conftst
       FPPFLAGS=$ac_save_FPPFLAGS
    fi

    if test $ac_fpp_need_subs = yes; then
        AC_LINK_IFELSE([_AC_LANG_PROGRAM_FPP_SUBS],
           [ac_prog_fc_cpp_subs=yes], 
           [ac_fpp_ok=no])
    fi

    if test $ac_fpp_need_wrap = yes; then
        AC_LINK_IFELSE([_AC_LANG_PROGRAM_FPP_WRAP],
           [ac_prog_fc_cpp_wrap=yes], 
           [ac_fpp_ok=no])
    fi

    if test $ac_fpp_need_CSTYLE = yes; then
        AC_LINK_IFELSE([_AC_LANG_PROGRAM_FPP_CSTYLE],
           [ac_prog_fc_cpp_CSTYLE=yes], 
           [ac_fpp_ok=no])
    fi

    if test $ac_fpp_need_cxxstyle = yes; then
        AC_LINK_IFELSE([_AC_LANG_PROGRAM_FPP_CXXSTYLE],
           [ac_prog_fc_cpp_cxxstyle=yes], 
           [ac_fpp_ok=no])
    fi

fi
if test $ac_fc_testing_fpp = direct; then
  FPPFLAGS=$ac_first_save_FPPFLAGS
fi
AC_MSG_RESULT([done.])
])#_AC_PROG_FC_CPP


# _AC_FPP_BUILD_RULE
# ------------------
# Figure out how to build from cpp/Fortran sources
#
# If we need to use a separate preprocessor, we must override make's
# `direct' .fpp.o rule in order to do `indirect' compilation
# (.fpp -> .f then .f -> .o).
#
# Configure variables set here:
#
#   FPP_SRC_EXT
#     This is the file extension (without a dot) of preprocessable
#     Fortran files; by default this is `fpp', but `F' is also
#     common. [set in AC_PROG_FPP, but documented here for consistency]
#
#   FPP_COMPILE_EXT
#     This contains the file extension which the Fortran compiler will
#     accept.  It is $FPP_SRC_EXT if the compiler itself can do
#     preprocessing (`direct' compilation), or a dummy value (not .f,
#     since that would conflict with the pre-existing rule for
#     compiling Fortran) if it cannot, and the file must be
#     preprocessed separately (`indirect').
#
#   FPP_PREPROCESS_EXT
#     The partner of FPP_COMPILE_EXT.  This is $FPP_SRC_EXT for
#     indirect compilation, or a dummy value for direct compilation,
#     when the corresponding separate-processing generated rule should
#     be ignored.
#
#   FPP_MAKE_FLAGS
#     This is used to include CPP/FPP related flags into the compiler
#     call if we compile directly, and leave them out otherwise.
#
#   FPP_OUTPUT
#     This is used to redirect FPP output to the .f file in those
#     cases where FPP writes to stdout rather than to a file.
#
# These are used in Automake's lang_ppfc_finish subroutine.
#
# FIXME: Martin Wilck's original version of these macros noted that it
# was necessary to generate explicit rules for .fpp -> .o compilations
# in order to override make's builtin rules in a portable manner
# (i.e. without using make extensions).  Is this true?
# POSIX/Single-Unix states that inference rules can be redefined, and
# there's no warning against this in Autoconf's section on
# `Limitations of Make'.
#
AC_DEFUN([_AC_FPP_BUILD_RULE],
[AC_CACHE_CHECK([how to build from preprocessed Fortran sources],
  ac_cv_fpp_build_rule,
[ac_cv_fpp_build_rule=
if test $ac_cv_prog_fc_cpp_ok = yes; then
   ac_cv_fpp_build_rule=direct
   ac_fpp_status=ok
elif test $ac_cv_prog_fpp_ok = yes; then
   ac_cv_fpp_build_rule=indirect
   ac_fpp_status=ok
dnl elif test $ac_cv_prog_fc_cpp = no && $ac_cv_prog_fpp = no; then
dnl    ac_cv_fpp_build_rule=
dnl    ac_fpp_status=fatal
dnl elif test $ac_cv_prog_fc_cpp = no; then
dnl    ac_cv_fpp_build_rule=indirect
dnl    ac_fpp_status=fail
dnl elif test -z "$ac_cv_prog_fpp"; then
dnl    ac_cv_fpp_build_rule=direct
dnl    ac_fpp_status=fail
dnl elif test $ac_fpp_equals_fc = yes; then
dnl    ac_cv_fpp_build_rule=direct
dnl    ac_fpp_status=fail
dnl else   
dnl    ac_fpp_status=fail
dnl # Both methods - direct and indirect - fail to meet some requirements
dnl # We use a simple approach to choose the best alternative
dnl    ac_fc_score=0
dnl    ac_fpp_score=0
dnl    for ac_i in d i subs wrap cstyle CSTYLE; do
dnl       ac_tmp="\$ac_fpp_need_$ac_i"
dnl       ac_need=`eval echo $ac_tmp`
dnl       if test $ac_need = yes; then
dnl         ac_tmp="\$ac_cv_prog_fc_cpp_$ac_i"
dnl         ac_v_fc=`eval echo $ac_tmp`
dnl         ac_tmp="\$ac_cv_prog_fpp_$ac_i"
dnl         ac_v_fpp=`eval echo $ac_tmp`
dnl         test "x$ac_v_fc" = xyes && ac_fc_score=`expr $ac_fc_score + 1`
dnl         test "x$ac_v_fpp" = xyes && ac_fpp_score=`expr $ac_fpp_score + 1`
dnl       fi
dnl    done
dnl    if test $ac_fpp_score -gt $ac_fc_score; then
dnl      ac_cv_fpp_build_rule=indirect
dnl    else
dnl      ac_cv_fpp_build_rule=direct
dnl    fi
fi   
])

if test -z "$ac_fpp_status"; then
    # must be non-empty below
    ac_fpp_status=cached
fi

if test $ac_fpp_status = fatal; then

  AC_MSG_ERROR([cannot find a valid build rule for Fortran/cpp source: consider installing the free Fortran preprocessor fpp from ftp.netlib.org/fortran])
	
elif test $ac_fpp_status = fail; then

  AC_MSG_WARN([cannot find a build rule for Fortran/cpp source that fulfills all requirements.  The compilation may fail or run time errors may arise.  Consider installing the free Fortran preprocessor fpp from ftp.netlib.org/fortran])

fi

# FIXME: Are DEFAULT_INCLUDES and AM_CPPFLAGS the right set of
# variables to be including here?  I confess (NG) to being somewhat
# confused about which one's which (same applies to rules in lang-compile.am).
if test $ac_cv_fpp_build_rule = direct; then
    FPP_COMPILE_EXT=$FPP_SRC_EXT
    FPP_PREPROCESS_EXT=DUMMY$FPP_SRC_EXT
    FPP_MAKE_FLAGS='$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS)'
else
    FPP_COMPILE_EXT=DUMMYf
    FPP_PREPROCESS_EXT=$FPP_SRC_EXT
    FPP_MAKE_FLAGS=""
fi
if test -z "$ac_fpp_out"; then
   FPP_OUTPUT=" "
else
   FPP_OUTPUT=">\[$]@"
fi
])# _AC_FPP_BUILD_RULE


# -----------------------
# User macros
# -----------------------

# AC_PROG_FPP([required features])
# --------------------------------------------------
#
# [required features] is a space-separated list of features that the Fortran
# preprocessor must have for the code to compile.
# It is up to the package maintainer to properly set these requirements.
#
# This macro will find out how to compile a preprocessable fixed-form
# file, with a .F file extension. To the best of my knowledge, such a
# file is compilable everywhere (albeit flags may be needed on
# case-insensitive filesystems)
#
# This macro should be followed by calling AC_FPP_FIXEDFORM([.srcext])
# and AC_FPP_FREEFORM([.srcext]) as appropriate for whichever source
# extensions are used in the user's project.
#
# This will fail to give the correct result when fixed-format files may be
# preprocessed directly by the compiler, but free-format ones
# may not. 
#
# Supported features are:
#
# include   : correctly process #include directives and -I
# define    : correctly process -D
# substitute: substitute macros in Fortran code 
#             (some preprocessors touch only lines starting with #)
# wrap      : wrap lines that become too long through macro substitution  
#             fpp is probably the only preprocessor that does this.
# cstyle    : Do not suppress C style comments (-C option in cpp)
# CSTYLE    : *Do* suppress C style comments
#             (e.g. code contains C-style comments, and compiler may not
#             know how to handle them)
# cxxstyle  : Do not suppress C++ style comments (default)
# CXXSTYLE  : *Do* suppress C++ style comments (seems unlikely, but in here
#             for completeness
# 
# Features can be abbreviated: i, in, inc etc. are equivalent to include.
# Features can be deselected (feature not needed) by prepending "no", 
#   e.g. nodef (=nodefine), now (=nowrap).
#
# Default for the feature list is 
#       [include define substitute nowrap nocstyle noCSTYLE cxxstyle]
# Feature requirements corresponding to the defaults may be omitted
#
# Note that "wrap" implies "substitute", and CSTYLE and cstyle cannot
# be requested at the same time. The macro adjusts this automatically. 
#
# This macro sets and substitutes the variables FPP and FPPFLAGS, and
# causes to be set FPP_OUTPUT, FPP_MAKE_FLAGS, and FPP_COMPILE_EXT
# (actually set in macro _AC_FPP_BUILD_RULE)
#
# The macro depends on both FC and CPP, because we must possibly fall 
# back on CPP for preprocessing.
#
AC_DEFUN([AC_PROG_FPP],
[AC_REQUIRE([AC_PROG_FC])dnl
dnl We are not going to use AC_REQUIRE(AC_PROG_CPP) here for 
dnl two reasons:
dnl 1) we don't really need to if FC will preprocess itself
dnl 2) we can't pass in an optional parameter to change the
dnl    default CPP search order, which we need to. 
dnl AC_REQUIRE([AC_PROG_CPP([cpp])])dnl

# Prefer AC_PROG_FC to AC_PROG_F77
if test "X$F77" != X; then
    AC_MSG_WARN([Use A@&t@C_PROG_FC with A@&t@C_PROG_FPP, instead of A@&t@C_PROG_F77])
fi

AC_ARG_VAR([FPP], [Command to preprocess Fortran code])
AC_ARG_VAR([FPPFLAGS], [Flags for the Fortran preprocessor])
_AC_PROG_FPP_FEATURES([$1])

# We first try to use FC for compiling the source directly
# into object files
ac_fpp_compile='${FC-fc} -c $FPPFLAGS $FPPFLAGS_SRCEXT $FCFLAGS conftest.$ac_ext >&AS_MESSAGE_LOG_FD'
ac_fpp_link='${FC-fc} -o conftest${ac_exeext} $FPPFLAGS $FPPFLAGS_SRCEXT $FCFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&AS_MESSAGE_LOG_FD'

AC_LANG_PUSH(Preprocessed Fortran)
FPP_SRC_EXT=F

# _AC_PROG_FC_CPP stores results of the feature checks in non-cv variables,
# which we copy to cv variables afterwards.
# The reason for that is reusability of the macro for other cv variables (see below)
ac_fc_testing_fpp=direct
_AC_PROG_FC_CPP

AC_CACHE_CHECK([whether $FC compiles programs with cpp directives], 
   ac_cv_prog_fc_cpp, 
  [ac_cv_prog_fc_cpp=$ac_prog_fc_cpp])

if test $ac_prog_fc_cpp = yes; then

  if test $ac_fpp_need_d = yes; then
    AC_CACHE_CHECK([whether $FC accepts -D], 
       ac_cv_prog_fc_cpp_d, 
      [ac_cv_prog_fc_cpp_d=$ac_prog_fc_cpp_d])
  fi

  if test $ac_fpp_need_i = yes; then
    AC_CACHE_CHECK([whether $FC accepts -I], 
       ac_cv_prog_fc_cpp_i,
      [ac_cv_prog_fc_cpp_i=$ac_prog_fc_cpp_i])
  fi

  if test $ac_fpp_need_subs = yes; then
    AC_CACHE_CHECK([whether $FC substitutes macros in Fortran code], 
       ac_cv_prog_fc_cpp_subs,
      [ac_cv_prog_fc_cpp_subs=$ac_prog_fc_cpp_subs])
  fi

  if test $ac_fpp_need_wrap = yes; then 
    AC_CACHE_CHECK([whether $FC wraps long lines automatically], 
       ac_cv_prog_fc_cpp_wrap,
      [ac_cv_prog_fc_cpp_wrap=$ac_prog_fc_cpp_wrap])
  fi

# Don't need to test if $FC removes C++ comments - that 
# way madness lies.

fi # test $ac_prog_fc_cpp = yes

AC_CACHE_CHECK([whether $FC fulfills requested features],
  ac_cv_prog_fc_cpp_ok,
  [ac_cv_prog_fc_cpp_ok=$ac_fpp_ok])

# If so, we don't need to go any further.
if test $ac_fpp_ok = yes; then
  ac_cv_fpp_build_rule=direct
else

# Now we check how to invoke a preprocessor that outputs Fortran code
# that FC can understand
#FIXME: in a joint C/Fortran project, CPP might have already
# been defined. Here we are potentially (probably) redefining it.
# I don't think this matters. Not sure, though.
# In that case, AC_SUBST has already been called on CPP.
# We don't want to fail if we can't find cpp - we might be able
# to fall back on fpp.
#FIXME: actually , we should just prefer cpp to $CPP is all.
  AC_PROG_CPP([cpp],[],[])
# The next macro sets FPP (unless already set by the user)
  _AC_PROG_FPP
  _AC_PROG_FPP_P

# Redefine the compile and link commands for indirect compilation
  ac_fpp_compile='${FPP-fpp} $FPPFLAGS $FPPFLAGS_SRCEXT conftest.$ac_ext '"$ac_fpp_out"' && ${FC-fc} -c $FCFLAGS conftest.f >&AS_MESSAGE_LOG_FD'
  ac_fpp_link='${FPP-fpp} $FPPFLAGS conftest.$ac_ext $FPPFLAGS_SRCEXT '"$ac_fpp_out"' && ${FC-fc} -o conftest${ac_exeext} $FCFLAGS $LDFLAGS conftest.f $LIBS >&AS_MESSAGE_LOG_FD'
  
  ac_compile=$ac_fpp_compile
  ac_link=$ac_fpp_link
# Redo all the feature checks for indirect compilation.
  ac_fc_testing_fpp=indirect
  _AC_PROG_FC_CPP

  if test $ac_fpp_need_d = yes; then
    AC_CACHE_CHECK([whether $FPP accepts -D], 
       ac_cv_prog_fpp_d, 
      [ac_cv_prog_fpp_d=$ac_prog_fc_cpp_d])
  fi

  if test $ac_fpp_need_i = yes; then
    AC_CACHE_CHECK([whether $FPP accepts -I], 
       ac_cv_prog_fpp_i,
      [ac_cv_prog_fpp_i=$ac_prog_fc_cpp_i])
  fi

  if test $ac_fpp_need_subs = yes; then
    AC_CACHE_CHECK([whether $FPP substitutes macros in Fortran code], 
       ac_cv_prog_fpp_subs,
      [ac_cv_prog_fpp_subs=$ac_prog_fc_cpp_subs])
  fi
 
  if test $ac_fpp_need_wrap = yes; then 
    AC_CACHE_CHECK([whether $FPP wraps long lines automatically], 
       ac_cv_prog_fpp_wrap,
      [ac_cv_prog_fpp_wrap=$ac_prog_fc_cpp_wrap])
  fi

  if test $ac_fpp_need_CSTYLE = yes; then 
    AC_CACHE_CHECK([whether $FPP suppresses C-style comments], 
       ac_cv_prog_fpp_CSTYLE,
      [ac_cv_prog_fpp_CSTYLE=$ac_prog_fc_cpp_CSTYLE])

  elif test $ac_fpp_need_cstyle = yes; then 
  # It only makes sense to test this for indirect compilation, 
  # i.e., if .f files are generated
      _AC_PROG_FPP_CSTYLE
  fi

  if test $ac_fpp_need_cxxstyle = yes; then 
    AC_CACHE_CHECK([whether $FPP preserves C++-style comments], 
       ac_cv_prog_fpp_cxxstyle,
      [ac_cv_prog_fpp_cxxstyle=$ac_prog_fc_cpp_cxxstyle])
  fi

  AC_CACHE_CHECK([whether $FPP fulfills requested features],
    ac_cv_prog_fpp_ok,
   [ac_cv_prog_fpp_ok=$ac_fpp_ok])

  ac_cv_fpp_build_rule=indirect

fi # test ac_fpp_ok != yes

# We have all necessary information.
# It remains to construct optimal build rules 
# (direct: .fpp.o or indirect: .fpp.f)
# and carry out the substitutions.
dnl _AC_FPP_BUILD_RULE

if test -z "$ac_fpp_out"; then
   FPP_OUTPUT=" "
else
   FPP_OUTPUT=">\[$]@"
fi

AC_SUBST(FPP)
AC_SUBST(FPP_OUTPUT)
dnl AC_SUBST(FPPFLAGS)
dnl AC_SUBST(FPP_MAKE_FLAGS)
dnl AC_SUBST(FPP_COMPILE_EXT) 
dnl AC_SUBST(FPP_PREPROCESS_EXT) 
dnl AC_SUBST(FPP_SRC_EXT) 

AC_LANG_POP(Preprocessed Fortran)

])# AC_PROG_FPP