# SPDX-License-Identifier: BSD-2-Clause

# Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.

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

# set the project name
project(mkdssp VERSION 4.4.6 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(CTest)
include(FetchContent)
include(VersionString)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers")
elseif(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

# Optionally build a version to be installed inside CCP4
option(BUILD_FOR_CCP4 "Build a version to be installed in CCP4" OFF)
option(BUILD_DOCUMENTATION "Generate the documentation files using pandoc" OFF)

if(BUILD_FOR_CCP4)
	if("$ENV{CCP4}" STREQUAL "" OR NOT EXISTS $ENV{CCP4})
		message(FATAL_ERROR "A CCP4 built was requested but CCP4 was not sourced")
	else()
		list(APPEND CMAKE_MODULE_PATH "$ENV{CCP4}")
		list(APPEND CMAKE_PREFIX_PATH "$ENV{CCP4}")
		set(CMAKE_INSTALL_PREFIX "$ENV{CCP4}")

		if(WIN32)
			set(BUILD_SHARED_LIBS ON)
		endif()
	endif()
endif()

if(MSVC)
	# make msvc standards compliant...
	add_compile_options(/permissive-)
	add_link_options(/NODEFAULTLIB:library)

	macro(get_WIN32_WINNT version)
		if(WIN32 AND CMAKE_SYSTEM_VERSION)
			set(ver ${CMAKE_SYSTEM_VERSION})
			string(REPLACE "." "" ver ${ver})
			string(REGEX REPLACE "([0-9])" "0\\1" ver ${ver})

			set(${version} "0x${ver}")
		endif()
	endmacro()

	get_WIN32_WINNT(ver)
	add_definitions(-D_WIN32_WINNT=${ver})

	if(BUILD_SHARED_LIBS)
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
	else()
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
	endif()
endif()

# Create a revision file, containing the current git version info
write_version_header("${CMAKE_CURRENT_SOURCE_DIR}/src")

# Optionally use mrc to create resources
find_package(Mrc QUIET)

if(MRC_FOUND)
	option(USE_RSRC "Use mrc to create resources" ON)
else()
	message(STATUS "Not using resources since mrc was not found")
endif()

set(CMAKE_THREAD_PREFER_PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG)
find_package(Threads)

if(NOT PDB_REDO_META)
	if(CMAKE_VERSION GREATER_EQUAL 3.28)
		set(EXLC "EXCLUDE_FROM_ALL")
	endif()

	FetchContent_Declare(
		libmcfp
		${EXLC}
		GIT_REPOSITORY https://github.com/mhekkel/libmcfp
		GIT_TAG v1.3.1)

	FetchContent_MakeAvailable(libmcfp)

	set(CIFPP_DOWNLOAD_CCD OFF)

	FetchContent_Declare(
		cifpp
		${EXLC}
		GIT_REPOSITORY https://github.com/pdb-redo/libcifpp.git
		GIT_TAG v7.0.3
	)

	FetchContent_MakeAvailable(cifpp)
endif()

add_subdirectory(libdssp)

add_executable(mkdssp ${CMAKE_CURRENT_SOURCE_DIR}/src/mkdssp.cpp)

target_link_libraries(mkdssp PRIVATE libmcfp::libmcfp dssp::dssp)

if(USE_RSRC)
	mrc_target_resources(mkdssp
		${CIFPP_SHARE_DIR}/mmcif_pdbx.dic
		${CIFPP_SHARE_DIR}/mmcif_ddl.dic
		${CMAKE_CURRENT_SOURCE_DIR}/libdssp/mmcif_pdbx/dssp-extension.dic)
endif()

# Install rules
install(TARGETS ${PROJECT_NAME}
	RUNTIME DESTINATION bin
)

if(BUILD_DOCUMENTATION)
	add_subdirectory(doc)
endif()

if(UNIX)
	install(FILES doc/mkdssp.1
		DESTINATION share/man/man1)
endif()

if(EXISTS "${CCP4}/html")
	install(FILES doc/mkdssp.html
		DESTINATION ${CCP4}/html)
endif()

# test
if(BUILD_TESTING)
	add_subdirectory(test)
endif()

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")

# Tarball generation
set(CPACK_SOURCE_TGZ ON)
set(CPACK_SOURCE_TBZ2 OFF)
set(CPACK_SOURCE_TXZ OFF)
set(CPACK_SOURCE_TZ OFF)
set(CPACK_SOURCE_IGNORE_FILES "/data/components.cif;/build;/.vscode;/.git")
set(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}-${PROJECT_VERSION}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME ${CPACK_PACKAGE_FILE_NAME})

# NSIS options
set(CPACK_NSIS_MODIFY_PATH ON)

# configuration done, include CPack
include(CPack)
