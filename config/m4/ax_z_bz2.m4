# SPDX-License-Identifier: BSD-2-Clause
# 
# Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Description: m4 macro to detect whether the z and bz2 libraries needs to be linked
# when using boost::iostreams

AC_DEFUN([read_test], [
	AC_LANG_SOURCE(esyscmd(config/tools/m4esc.sh config/test/$1))])

AC_DEFUN([AX_IOSTREAMS_Z],
[
	AC_MSG_CHECKING([need to link zlib])

	save_CPPFLAGS="$CPPFLAGS"
	CPPFLAGS="$BOOST_CPPFLAGS $CPPFLAGS "

	save_LDFLAGS="$LDFLAGS"
	LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"

	save_LIBS="$LIBS"
	LIBS="$LIBS $BOOST $BOOST_IOSTREAMS_LIB"

	NEED_ZLIB=0

	AC_LINK_IFELSE(
		[read_test(iostreams-zlib-test.cpp)],
		[AC_MSG_RESULT([no])],
		[
			LIBS="$LIBS -lz"

			AC_LINK_IFELSE(
				[read_test(iostreams-zlib-test.cpp)],
				[
					AC_MSG_RESULT([yes])
					NEED_ZLIB=1
				],
				[AC_MSG_ERROR([Could not link iostreams, even with zlib])])
		])
	
	CPPFLAGS="$save_CPPFLAGS"
	LDFLAGS="$save_LDFLAGS"
	LIBS="$save_LIBS"

	if [ test "x$NEED_ZLIB" != "x" ]; then
		LIBS="$LIBS -lz"
	fi
])

AC_DEFUN([AX_IOSTREAMS_BZ2],
[
	AC_MSG_CHECKING([need to link bzip2])

	save_CPPFLAGS="$CPPFLAGS"
	CPPFLAGS="$BOOST_CPPFLAGS $CPPFLAGS "

	save_LDFLAGS="$LDFLAGS"
	LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"

	save_LIBS="$LIBS"
	LIBS="$LIBS $BOOST $BOOST_IOSTREAMS_LIB"

	NEED_BZLIB=0

	AC_LINK_IFELSE(
		[read_test(iostreams-bzip2-test.cpp)],
		[AC_MSG_RESULT([no])],
		[
			LIBS="$LIBS -lbz2"

			AC_LINK_IFELSE(
				[read_test(iostreams-bzip2-test.cpp)],
				[
					AC_MSG_RESULT([yes])
					NEED_BZLIB=1
				],
				[AC_MSG_ERROR([Could not link iostreams, even with bzip2])])
		])
	
	CPPFLAGS="$save_CPPFLAGS"
	LDFLAGS="$save_LDFLAGS"
	LIBS="$save_LIBS"

	if [ test "x$NEED_BZLIB" != "x" ]; then
		LIBS="$LIBS -lbz2"
	fi
])