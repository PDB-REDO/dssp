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
dnl AC_CONFIG_MACRO_DIR([config/m4])

AC_PREFIX_DEFAULT(/usr/local)

AC_PROG_INSTALL

AC_ARG_VAR([DEBUG], [Build a debug version of the application])

AC_PATH_PROG([PKG_CONFIG], [pkg-config])

AC_ARG_VAR([MRC], [Specify a location for the mrc executable])

dnl We really want to use mrc
if test "x$MRC" = "x"; then
	AC_PATH_PROG([MRC], [mrc])
fi

if test "x$MRC" = "x"; then
	AC_MSG_ERROR([mrc not found, the application will be built without resources])
fi

AX_PTHREAD

AC_ARG_WITH([cif++],
	AS_HELP_STRING([--with-cif++=@<:@location@:>@],
		[Use the cif++ library as specified.]),
		[
			CPPFLAGS="$CPPFLAGS -I ${withval}/include"
			CXXFLAGS="$CXXFLAGS -I ${withval}/include"
			LDFLAGS="$LDFLAGS -L${withval}/.libs"
			AC_SUBST([CIFPP_RSRC], [${withval}/rsrc]) 
		])

AX_CHECK_LIBRARY([LIBCIFPP], [cif++/Config.hpp], [cif++],
	[],
	[AC_MSG_ERROR([libcif++ not found - cannot continue])])

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
AX_BOOST_THREAD
AX_BOOST_REGEX
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_DATE_TIME

AX_CHECK_LIBRARY([LIBZ], [zlib.h], [z],
				 [ LIBS="$LIBS -lz" ],
                 [AC_MSG_ERROR([libz not found - compressed files not supported])])
AX_CHECK_LIBRARY([LIBBZ2], [bzlib.h], [bz2],
				 [ LIBS="$LIBS -lbz2"],
                 [AC_MSG_ERROR([libbz2 not found - compressed files not supported])])

AC_SEARCH_LIBS([dladdr1], [dl],
	[
		AC_DEFINE([HAVE_DLADDR], 1, [Define if the system has dladdr])
	], [], [])

AC_SUBST([HAVE_DLADDR])

dnl Process Makefile.in to create Makefile
AC_OUTPUT([GNUmakefile])
