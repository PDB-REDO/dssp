# SPDX-License-Identifier: BSD-2-Clause

# Copyright (c) 2021-2023 NKI/AVL, Netherlands Cancer Institute

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

# This cmake extension writes out a revision.hpp file in a specified directory.
# The file will contain a C++ inline function that can be used to write out
# version information.

cmake_minimum_required(VERSION 3.15)

# We want the revision.hpp file to be updated whenever the status of the
# git repository changes. Use the same technique as in GetGitRevisionDescription.cmake


#[=======================================================================[.rst:
.. command:: write_version_header

  Write a file named revision.hpp containing version info::

	write_version_header(<destdir>
	                     [FILE_NAME <file-name>]
						 [LIB_NAME <library-name>]
	                    )
  
  This command will generate the code to write a file name
  revision.hpp in the directory ``<destdir>``.
  
  ``FILE_NAME``
	Specify the name of the file to create, default is ``revision.hpp``.

  ``LIB_NAME``
	Specify the library name which will be used as a prefix part for the
	variables contained in the revision file.
#]=======================================================================]

# First locate a .git file or directory.
function(_get_git_dir _start_dir _variable)

	set(cur_dir "${_start_dir}")
	set(git_dir "${_start_dir}/.git")

	while(NOT EXISTS "${git_dir}")
		# .git dir not found, search parent directories
		set(prev_dir "${cur_dir}")
		get_filename_component(cur_dir "${cur_dir}" DIRECTORY)
		if(cur_dir STREQUAL prev_dir OR cur_dir STREQUAL ${_start_dir})
			# we are not in git since we either hit root or
			# the CMAKE_SOURCE_DIR which should be the top
			set(${_variable}
				""
				PARENT_SCOPE)
			return()
		endif()
		set(git_dir "${cur_dir}/.git")
	endwhile()

	set(${_variable} "${git_dir}" PARENT_SCOPE)
endfunction()

# Locate the git refspec hash and load the hash
# This code locates the file containing the git refspec/hash
# and loads it. Doing it this way assures that each time the git
# repository changes the revision.hpp file gets out of date.
function(_get_git_hash _data_dir _variable)

	# Be pessimistic
	set(_variable "" PARENT_SCOPE)

	# Load git package if needed
	if(NOT GIT_FOUND)
		find_package(Git QUIET)
	endif()

	# And fail if not found
	if(NOT GIT_FOUND)
		return()
	endif()

	# Locate the nearest .git file or directory
	_get_git_dir(${CMAKE_CURRENT_SOURCE_DIR} GIT_DIR)

	# And fail if not found
	if("${GIT_DIR}" STREQUAL "")
        return()
    endif()

    # Check if the current source dir is a git submodule or a worktree.
    # In both cases .git is a file instead of a directory.
    #
    if(IS_DIRECTORY ${GIT_DIR})
		set(HEAD_SOURCE_FILE "${GIT_DIR}/HEAD")
	else()
		# The following git command will return a non empty string that
        # points to the super project working tree if the current
        # source dir is inside a git submodule.
        # Otherwise the command will return an empty string.
        #
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" rev-parse
                    --show-superproject-working-tree
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            OUTPUT_VARIABLE out
            ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(NOT "${out}" STREQUAL "")
            # If out is not empty, GIT_DIR/CMAKE_CURRENT_SOURCE_DIR is in a submodule
            file(READ ${GIT_DIR} submodule)
            string(REGEX REPLACE "gitdir: (.*)$" "\\1" GIT_DIR_RELATIVE
                                 ${submodule})
            string(STRIP ${GIT_DIR_RELATIVE} GIT_DIR_RELATIVE)
            get_filename_component(SUBMODULE_DIR ${GIT_DIR} PATH)
            get_filename_component(GIT_DIR ${SUBMODULE_DIR}/${GIT_DIR_RELATIVE}
                                   ABSOLUTE)
            set(HEAD_SOURCE_FILE "${GIT_DIR}/HEAD")
        else()
            # GIT_DIR/CMAKE_CURRENT_SOURCE_DIR is in a worktree
            file(READ ${GIT_DIR} worktree_ref)
            # The .git directory contains a path to the worktree information directory
            # inside the parent git repo of the worktree.
            #
            string(REGEX REPLACE "gitdir: (.*)$" "\\1" git_worktree_dir
                                 ${worktree_ref})
            string(STRIP ${git_worktree_dir} git_worktree_dir)
            _get_git_dir("${git_worktree_dir}" GIT_DIR)
            set(HEAD_SOURCE_FILE "${git_worktree_dir}/HEAD")
        endif()
	endif()

	# Fail if the 'head' file was not found
    if(NOT EXISTS "${HEAD_SOURCE_FILE}")
        return()
    endif()

	# Make a copy of the head file
    set(HEAD_FILE "${_data_dir}/HEAD")
    configure_file("${HEAD_SOURCE_FILE}" "${HEAD_FILE}" COPYONLY)

	# Now we create a cmake file that will read the contents of this
	# head file in the appropriate way
	file(WRITE "${_data_dir}/grab-ref.cmake.in" [[
set(HEAD_HASH)

file(READ "@HEAD_FILE@" HEAD_CONTENTS LIMIT 1024)

string(STRIP "${HEAD_CONTENTS}" HEAD_CONTENTS)
if(HEAD_CONTENTS MATCHES "ref")
	# named branch
	string(REPLACE "ref: " "" HEAD_REF "${HEAD_CONTENTS}")
	if(EXISTS "@GIT_DIR@/${HEAD_REF}")
		configure_file("@GIT_DIR@/${HEAD_REF}" "@VERSION_STRING_DATA@/head-ref" COPYONLY)
	else()
		configure_file("@GIT_DIR@/packed-refs" "@VERSION_STRING_DATA@/packed-refs" COPYONLY)
		file(READ "@VERSION_STRING_DATA@/packed-refs" PACKED_REFS)
		if(${PACKED_REFS} MATCHES "([0-9a-z]*) ${HEAD_REF}")
			set(HEAD_HASH "${CMAKE_MATCH_1}")
		endif()
	endif()
else()
	# detached HEAD
	configure_file("@GIT_DIR@/HEAD" "@VERSION_STRING_DATA@/head-ref" COPYONLY)
endif()

if(NOT HEAD_HASH)
	file(READ "@VERSION_STRING_DATA@/head-ref" HEAD_HASH LIMIT 1024)
	string(STRIP "${HEAD_HASH}" HEAD_HASH)
endif()
]])

    configure_file("${VERSION_STRING_DATA}/grab-ref.cmake.in"
                   "${VERSION_STRING_DATA}/grab-ref.cmake" @ONLY)
    
	# Include the aforementioned file, this will define
	# the HEAD_HASH variable we're looking for
	include("${VERSION_STRING_DATA}/grab-ref.cmake")

    set(${_variable} "${HEAD_HASH}" PARENT_SCOPE)
endfunction()

# Create a revision file, containing the current git version info, if any
function(write_version_header dir)

	set(flags )
	set(options LIB_NAME FILE_NAME)
	set(sources )
	cmake_parse_arguments(VERSION_STRING_OPTION "${flags}" "${options}" "${sources}" ${ARGN})

	# parameter check
	if(NOT IS_DIRECTORY ${dir})
		message(FATAL_ERROR "First parameter to write_version_header should be a directory where the final revision.hpp file will be placed")
	endif()

	if(VERSION_STRING_OPTION_FILE_NAME)
		set(file_name "${VERSION_STRING_OPTION_FILE_NAME}")
	else()
		set(file_name "revision.hpp")
	endif()

	# Where to store intermediate files
	set(VERSION_STRING_DATA "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/VersionString")
	if(NOT EXISTS "${VERSION_STRING_DATA}")
        file(MAKE_DIRECTORY "${VERSION_STRING_DATA}")
    endif()

	# Load the git hash using the wizzard-like code above.
	_get_git_hash("${VERSION_STRING_DATA}" GIT_HASH)

	# If git was found, fetch the git description string
	if(GIT_HASH)
		execute_process(
			COMMAND "${GIT_EXECUTABLE}" describe --dirty --match=build
			WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
			RESULT_VARIABLE res
			OUTPUT_VARIABLE out
			ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

		if(res EQUAL 0)
			set(REVISION_STRING "${out}")
		endif()
	endif()

	# Check the revision string, if it matches we fill in the required info
	if(REVISION_STRING MATCHES "build-([0-9]+)-g([0-9a-f]+)(-dirty)?")
		set(BUILD_NUMBER ${CMAKE_MATCH_1})
		if(CMAKE_MATCH_3)
			set(REVISION_GIT_TAGREF "${CMAKE_MATCH_2}*")
		else()
			set(REVISION_GIT_TAGREF "${CMAKE_MATCH_2}")
		endif()

		string(TIMESTAMP REVISION_DATE_TIME "%Y-%m-%dT%H:%M:%SZ" UTC)
	else()
		message(STATUS "no git info available, cannot update version string")

		set(REVISION_GIT_TAGREF "")
		set(BUILD_NUMBER 0)
		set(REVISION_DATE_TIME "")
	endif()

	if(VERSION_STRING_OPTION_LIB_NAME)
		set(VAR_PREFIX "${VERSION_STRING_OPTION_LIB_NAME}")
		set(IDENT_PREFIX "${VERSION_STRING_OPTION_LIB_NAME}_")
	else()
		set(VAR_PREFIX "")
		set(IDENT_PREFIX "")
	endif()

	file(WRITE "${VERSION_STRING_DATA}/${file_name}.in" [[// This file was generated by VersionString.cmake

#pragma once

#include <ostream>

constexpr const char k@VAR_PREFIX@ProjectName[] = "@PROJECT_NAME@";
constexpr const char k@VAR_PREFIX@VersionNumber[] = "@PROJECT_VERSION@";
constexpr int k@VAR_PREFIX@BuildNumber = @BUILD_NUMBER@;
constexpr const char k@VAR_PREFIX@RevisionGitTag[] = "@REVISION_GIT_TAGREF@";
constexpr const char k@VAR_PREFIX@RevisionDate[] = "@REVISION_DATE_TIME@";

#ifndef VERSION_INFO_DEFINED
#define VERSION_INFO_DEFINED 1

class version_info_base
{
  public:
	virtual ~version_info_base() = default;

	static void write(std::ostream &os, bool verbose)
	{
		auto s_head = head();
		if (s_head != nullptr)
			write(s_head, os, verbose);
	}

  protected:

	struct instance
	{
		const char *m_name;
		const char *m_version;
		int m_build;
		const char *m_git_tag;
		const char *m_revision_date;

		instance *m_next = nullptr;
	};

	static void write(const instance *inst, std::ostream &os, bool verbose)
	{
		if (inst->m_next != nullptr)
		{
			write(inst->m_next, os, verbose);
			if (not verbose)
				return;
			os << '-' << std::endl;
		}

		os << inst->m_name << " version " << inst->m_version << std::endl;

		if (verbose)
		{
			if (inst->m_build != 0)
			{
				os << "build: " << inst->m_build << ' ' << inst->m_revision_date << std::endl;
				if (inst->m_git_tag[0] != 0)
					os << "git tag: " << inst->m_git_tag << std::endl;
			}
			else
				os << "No revision information available" << std::endl;
		}
	}

	using instance_ptr = instance *;

	static instance_ptr &head()
	{
		static instance_ptr s_head = nullptr;
		return s_head;
	}
};

template<typename T>
class version_info : public version_info_base
{
  public:
	using implementation_type = T;

	struct register_object
	{
		register_object()
		{
			auto &s_head = version_info_base::head();
			static instance s_next{
				implementation_type::name(),
				implementation_type::version(),
				implementation_type::build_number(),
				implementation_type::git_tag(),
				implementation_type::revision_date(),
				s_head };
			s_head = &s_next;
		}
	};

	template<register_object&> struct referrence_object;

	static register_object s_registered_object;
	static referrence_object<s_registered_object> s_referrenced_object;
};

template<typename T> typename version_info<T>::register_object version_info<T>::s_registered_object;

inline void write_version_string(std::ostream &os, bool verbose)
{
	version_info_base::write(os, verbose);
}

#endif

class version_info_@IDENT_PREFIX@impl : public version_info<version_info_@IDENT_PREFIX@impl>
{
  public:
	static constexpr const char *name() { return k@VAR_PREFIX@ProjectName; }
	static constexpr const char *version() { return k@VAR_PREFIX@VersionNumber; }
	static constexpr int build_number() { return k@VAR_PREFIX@BuildNumber; }
	static constexpr const char *git_tag() { return k@VAR_PREFIX@RevisionGitTag; }
	static constexpr const char *revision_date() { return k@VAR_PREFIX@RevisionDate; }
};
]])
	configure_file("${VERSION_STRING_DATA}/${file_name}.in" "${dir}/${file_name}" @ONLY)
endfunction()

