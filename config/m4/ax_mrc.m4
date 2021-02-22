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
# Description: m4 macro to detect mrc, the resource compiler

AC_DEFUN([AX_MRC],
[
	AC_ARG_VAR([MRC], [Specify a location for the mrc executable])

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
			[AS_HELP_STRING([--disable-resources], [Do not use mrc to store data in resources])])

		AS_IF([test "x$enable_resources" != "xno" ], [
			USE_RSRC=1
		])
	fi

	AC_SUBST([USE_RSRC], [$USE_RSRC])

	AC_DEFINE_UNQUOTED([USE_RSRC], [$USE_RSRC], [Use mrc to store resources])	
])