# SPDX-License-Identifier: BSD-2-Clause

# Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

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

cmake_minimum_required(VERSION 3.15)

# Create a revision file, containing the current git version info, if any
function(write_version_header)
	include(GetGitRevisionDescription)
	if(NOT(GIT-NOTFOUND OR HEAD-HASH-NOTFOUND))
		git_describe_working_tree(BUILD_VERSION_STRING --match=build --dirty)

		if(BUILD_VERSION_STRING MATCHES "build-([0-9]+)-g([0-9a-f]+)(-dirty)?")
			set(BUILD_GIT_TAGREF "${CMAKE_MATCH_2}")
			if(CMAKE_MATCH_3)
				set(BUILD_VERSION_STRING "${CMAKE_MATCH_1}*")
			else()
				set(BUILD_VERSION_STRING "${CMAKE_MATCH_1}")
			endif()
		endif()
	else()
		set(BUILD_VERSION_STRING "no git info available")
	endif()

	include_directories(${PROJECT_BINARY_DIR} PRIVATE)
	string(TIMESTAMP BUILD_DATE_TIME "%Y-%m-%dT%H:%M:%SZ" UTC)

	if(ARGC GREATER 0)
		set(VAR_PREFIX "${ARGV0}")
	endif()

	file(WRITE "${PROJECT_BINARY_DIR}/revision.hpp.in" [[// Generated revision file

#pragma once

#include <ostream>

const char k@VAR_PREFIX@ProjectName[] = "@PROJECT_NAME@";
const char k@VAR_PREFIX@VersionNumber[] = "@PROJECT_VERSION@";
const char k@VAR_PREFIX@VersionGitTag[] = "@BUILD_GIT_TAGREF@";
const char k@VAR_PREFIX@BuildInfo[] = "@BUILD_VERSION_STRING@";
const char k@VAR_PREFIX@BuildDate[] = "@BUILD_DATE_TIME@";

inline void write_version_string(std::ostream &os, bool verbose)
{
	os << k@VAR_PREFIX@ProjectName << " version " << k@VAR_PREFIX@VersionNumber << std::endl;
	if (verbose)
	{
		os << "build: " << k@VAR_PREFIX@BuildInfo << ' ' << k@VAR_PREFIX@BuildDate << std::endl;
		if (k@VAR_PREFIX@VersionGitTag[0] != 0)
			os << "git tag: " << k@VAR_PREFIX@VersionGitTag << std::endl;
	}
}
]])
	configure_file("${PROJECT_BINARY_DIR}/revision.hpp.in" "${PROJECT_BINARY_DIR}/revision.hpp" @ONLY)
endfunction()

