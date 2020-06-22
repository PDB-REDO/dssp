AC_INIT([dssp], 3.0, [m.hekkelman@nki.nl])

dnl Switch to a C++ compiler, and check if it works.
AC_LANG(C++)
AX_CXX_COMPILE_STDCXX_17([noext])

AC_CONFIG_SRCDIR([src/dssp.cpp])
AC_CONFIG_AUX_DIR(config)
dnl AC_CONFIG_MACRO_DIR([config/m4])

AC_PREFIX_DEFAULT(/usr/local)

AC_PROG_INSTALL

AC_ARG_VAR([DEBUG], [Build a debug version of the application])

AC_PATH_PROG([PKG_CONFIG], [pkg-config])

AX_PTHREAD

AX_BOOST_BASE([1.65.1], [], [AC_MSG_ERROR([Could not find a recent version of boost])])

AC_ARG_WITH([cif++],
	AS_HELP_STRING([--with-cif++=@<:@location@:>@],
		[Use the cif++ library as specified.]),
		[
			CPPFLAGS="$CPPFLAGS -I ${withval}/include"
			CXXFLAGS="$CXXFLAGS -I ${withval}/include"
			LDFLAGS="$LDFLAGS -L${withval}/.libs"
		])

AX_CHECK_LIBRARY([LIBCIFPP], [cif++/Config.h], [cif++],
	[],
	[AC_MSG_ERROR([libcif++ not found - cannot continue])])

dnl Process Makefile.in to create Makefile
AC_OUTPUT([GNUmakefile])
