AC_INIT([dssp], 4.0, [m.hekkelman@nki.nl])

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

AX_BOOST_BASE([1.65.1], [], [AC_MSG_ERROR([Could not find a recent version of boost])])
AX_BOOST_IOSTREAMS
AX_BOOST_THREAD
AX_BOOST_REGEX
AX_BOOST_PROGRAM_OPTIONS

AC_ARG_VAR([CCP4], [The location where CCP4 is installed])

AS_IF([test x"$CCP4" != x""],
	  [
		CPPFLAGS="$CPPFLAGS -I ${CCP4}/include"
		CXXFLAGS="$CXXFLAGS -I ${CCP4}/include"
		LDFLAGS="$LDFLAGS -L${CCP4}/lib -Wl,-rpath=${CCP4}/lib"
	  ])


AC_ARG_WITH([clipper],
	AS_HELP_STRING([--with-clipper=@<:@location@:>@],
		[Use the clipper library as specified.
		@<:@location=$CCP4@:>@]),
		[
			CPPFLAGS="$CPPFLAGS -I ${withval}/include"
			CXXFLAGS="$CXXFLAGS -I ${withval}/include"
			LDFLAGS="$LDFLAGS -L${withval}/lib"
			LIBS="$LIBS -lclipper-ccp4 -lclipper-minimol -lclipper-mmdb -lclipper-cns -lclipper-phs -lclipper-contrib -lclipper-core"
		])

AC_CHECK_HEADER(
	[clipper/clipper.h],
	[],
	[AC_MSG_ERROR([
Can't find the main clipper include file clipper/clipper.h. Please
specify either the CCP4 location of the location of the installed
clipper using --with-clipper=<path>
])])

AX_CHECK_LIBRARY([LIBCLIPPER], [clipper/clipper.h], [clipper-core],
	[],
	[AC_MSG_ERROR([libclipper not found - cannot continue])])

AC_ARG_WITH([newuoa],
	AS_HELP_STRING([--with-newuoa=@<:@location@:>@],
		[Use the newuoa library as specified.]),
		[
			CPPFLAGS="$CPPFLAGS -I ${withval}/include"
			CXXFLAGS="$CXXFLAGS -I ${withval}/include"
			LDFLAGS="$LDFLAGS -L${withval}/lib"
			LIBS="$LIBS -lnewuoa"
		])

AC_CHECK_HEADER(
	[newuoa.h],
	[],
	[AC_MSG_ERROR([
Can't find the newuoa include file newuoa.h. Please install this
library and specify its location with --with-newuoa=<path>.
You can find newuoa at https://github.com/elsid/newuoa-cpp.git
])])

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

dnl AC_MSG_CHECKING([clipper version])
dnl AC_COMPILE_IFELSE(
dnl 	[read_test(clipper-test.cpp)],
dnl 	[],
dnl 	[AC_MSG_ERROR([The version of clipper is not up to date])])
dnl 
dnl AC_CHECK_HEADER(
dnl 	[newuoa.h],
dnl 	[],
dnl 	[AC_MSG_ERROR([
dnl Can't find the newuoa include file newuoa.h. Please install this
dnl library and specify its location with --with-newuoa=<path>.
dnl You can find newuoa at https://github.com/elsid/newuoa-cpp.git
dnl ])])

dnl Process Makefile.in to create Makefile
AC_OUTPUT([GNUmakefile])
