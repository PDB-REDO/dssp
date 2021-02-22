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
# Description: m4 macro to detect std::filesystem and optionally the linker flags to use it
# 
# Description: Check for libcifpp

AC_DEFUN([AX_LIBCIFPP],
[
	AC_ARG_WITH([cif++],
		AS_HELP_STRING([--with-cif++=@<:@location@:>@],
			[Use the cif++ library as specified.]),
			[
				AS_IF([test -d ${withval}/include], [], [
					AC_MSG_ERROR(['${withval}'' is not a valid directory for --with-cif++])
				])

				CIFPP_CFLAGS="-I ${withval}/include"
				CIFPP_LIBS="-L${withval}/.libs -lcifpp"
				LIBCIFPP_DATA_DIR="${withval}/rsrc"

				AC_SUBST([CIFPP_CFLAGS], [$CIFPP_CFLAGS])
				AC_SUBST([CIFPP_LIBS], [$CIFPP_LIBS])
			])

	AS_IF([test "x$CIFPP_LIBS" = "x"], [
		if test -x "$PKG_CONFIG"
		then
			AX_PKG_CHECK_MODULES([CIFPP], [libcifpp], [], [], [AC_MSG_ERROR([the required package libcif++ is not installed])])
			LIBCIFPP_DATA_DIR=$(pkg-config --variable=datalibdir libcifpp)
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
			AS_IF([ test -f /usr/local/share/libcifpp/mmcif_pdbx_v50.dic.gz ], [LIBCIFPP_DATA_DIR=/usr/local/share/libcifpp/ ])
			AS_IF([ test -f /var/cache/libcifpp/mmcif_pdbx_v50.dic.gz ], [LIBCIFPP_DATA_DIR=/var/cache/libcifpp ])
		fi
	])

	AC_ARG_VAR([LIBCIFPP_DATA_DIR], [Directory containing mmcif_pdbx_v50.dic file])
	AC_SUBST([LIBCIFPP_DATA_DIR], [$LIBCIFPP_DATA_DIR])
])