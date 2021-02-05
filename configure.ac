AC_PREREQ([2.69])

AC_INIT([dssp], 4.0, [m.hekkelman@nki.nl])

dnl Switch to a C++ compiler, and check if it works.
AC_LANG(C++)
AX_CXX_COMPILE_STDCXX_17([noext])

AX_CHECK_COMPILE_FLAG([-fstandalone-debug],
	[
		CXXFLAGS="$CXXFLAGS -fstandalone-debug"
	] , , [-Werror])

AC_CONFIG_SRCDIR([src/dssp.cpp])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_HEADERS([src/config.hpp])

AC_PREFIX_DEFAULT(/usr/local)

AC_PROG_INSTALL
PKG_PROG_PKG_CONFIG

AC_ARG_VAR([DEBUG], [Build a debug version of the application])

AC_ARG_VAR([MRC], [Specify a location for the mrc executable])

PKG_PROG_PKG_CONFIG

dnl using resources?
USE_RSRC=0

if test "x$MRC" = "x"; then
	AC_PATH_PROG([MRC], [mrc])
fi

if test "x$MRC" = "x"; then
	AC_MSG_WARN([The mrc application was not found, not using resources.])
else
	AC_ARG_ENABLE(
		resources,
		[AS_HELP_STRING([--enable-resources], [Use mrc to store data in resources])])

	AS_IF([test "x$enable_resources" = "xyes" ], [
		USE_RSRC=1
	])
fi

AC_SUBST([USE_RSRC], [$USE_RSRC])

AC_DEFINE_UNQUOTED([USE_RSRC], [$USE_RSRC], [Use mrc to store resources])

dnl revision numbering is something used internally at the NKI
AC_ARG_ENABLE(
	revision,
	[AS_HELP_STRING([--enable-revision], [Create a build number as revision])])

AS_IF([test "x$enable_revision" = "xyes" ], [
	UPDATE_REVISION=1
])

AC_SUBST([UPDATE_REVISION], [$UPDATE_REVISION])

AX_PTHREAD

AC_CHECK_HEADER([filesystem], [], [AC_MSG_ERROR([The file <filesystem> is missing, perhaps you should install a more recent libstdc++ implementation.])])

dnl check if we need stdc++fs as library
AC_TRY_LINK(
	[#include <filesystem>],
	[(void)std::filesystem::current_path();],
	[],
	[
		LIBS="$LIBS -lstdc++fs"

		AC_TRY_LINK(
			[#include <filesystem>],
			[(void)std::filesystem::current_path();],
			[],
			[
				AC_MSG_ERROR([Could not link filesystem])
			]
		)
	]
)

AX_BOOST_BASE([1.65.1], [], [AC_MSG_ERROR([Could not find a recent version of boost])])
AX_BOOST_IOSTREAMS
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_DATE_TIME
AX_BOOST_REGEX

CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"

AC_ARG_WITH([cif++],
	AS_HELP_STRING([--with-cif++=@<:@location@:>@],
		[Use the cif++ library as specified.]),
		[
			AS_IF([test -d ${withval}/include], [], [
				AC_MSG_ERROR(['${withval}'' is not a valid directory for --with-cif++])
			])
			dnl AC_SUBST([CIFPP_CFLAGS], ["-I ${withval}/include"])
			dnl AC_SUBST([CIFPP_LIBS], ["-L${withval}/.libs -lcifpp"])

			CIFPP_CFLAGS="-I ${withval}/include"
			CIFPP_LIBS="-L${withval}/.libs -lcifpp"
			CIFPP_RSRC="${withval}/rsrc"

			AC_SUBST([CIFPP_CFLAGS], [$CIFPP_CFLAGS])
			AC_SUBST([CIFPP_LIBS], [$CIFPP_LIBS])
		])

AC_SUBST([CIFPP_RSRC], [$CIFPP_RSRC])

AS_IF([test "x$CIFPP_LIBS" = "x"], [
	if test -x "$PKG_CONFIG"
	then
		AX_PKG_CHECK_MODULES([CIFPP], [libcifpp], [], [], [AC_MSG_ERROR([the required package libcifpp is not installed])])
	else
		AC_CHECK_HEADER(
	        [cif++/Cif++.hpp],
	        [
				dnl CIFPP_CFLAGS="-I ${withval}/include"
			],
			[AC_MSG_ERROR([
Can't find the libcif++ header, Config.hpp.  Make sure that it
is installed, and either use the --with-cif++ option or install
pkg-config.])])

			AX_CHECK_LIBRARY([CIFPP], [cif++/Cif++.hpp], [cifpp],
				 [
					 LIBS="-lcifpp $LIBS"
				 ],
                 [AC_MSG_ERROR([libcif++ not found])])
	fi
])

dnl These are still needed outside the Debian environment
AX_CHECK_LIBRARY([LIBZ], [zlib.h], [z],
				 [ LIBS="$LIBS -lz" ],
                 [AC_MSG_ERROR([libz not found - compressed files not supported])])
AX_CHECK_LIBRARY([LIBBZ2], [bzlib.h], [bz2],
				 [ LIBS="$LIBS -lbz2"],
                 [AC_MSG_ERROR([libbz2 not found - compressed files not supported])])

AC_SUBST([LIBS], [$LIBS])

dnl Process Makefile.in to create Makefile
AC_OUTPUT([GNUmakefile])
