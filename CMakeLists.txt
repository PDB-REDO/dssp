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
project(mkdssp VERSION 4.1.0 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CheckFunctionExists)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(Dart)
include(FindFilesystem)
include(GenerateExportHeader)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Filesystem REQUIRED)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers")
elseif(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

if(NOT "$ENV{CCP4}" STREQUAL "")
	set(CCP4 $ENV{CCP4})
	list(PREPEND CMAKE_MODULE_PATH "${CCP4}/Lib")
	list(APPEND CMAKE_PREFIX_PATH ${CCP4})

	if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
		set(CMAKE_INSTALL_PREFIX ${CCP4})
	endif()
endif()

if(MSVC)
	# make msvc standards compliant...
	add_compile_options(/permissive-)

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
endif()

# Create a revision file, containing the current git version info
include(VersionString)
write_version_header()

# Optionally use mrc to create resources
find_package(Mrc)

if(MRC_FOUND)
	option(USE_RSRC "Use mrc to create resources" ON)
else()
	message(WARNING "Not using resources since mrc was not found")
endif()

if(USE_RSRC)
	set(USE_RSRC 1)

	message("Using resources compiled with ${MRC}")
	add_compile_definitions(USE_RSRC)
endif()

set(CMAKE_THREAD_PREFER_PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG)
find_package(Threads)

find_package(libconfig QUIET)
if(NOT libconfig_FOUND)
	add_subdirectory(libconfig EXCLUDE_FROM_ALL)
endif()

find_package(cifpp 5.0.0 REQUIRED)

add_subdirectory(date EXCLUDE_FROM_ALL)

# The DSSP code is in a separate library, optionally to be used by others
add_library(dssp_library OBJECT ${PROJECT_SOURCE_DIR}/src/DSSP.cpp)
add_library(dssp::dssp ALIAS dssp_library)

target_link_libraries(dssp_library cifpp::cifpp)

add_executable(mkdssp
	${PROJECT_SOURCE_DIR}/src/dssp_wrapper.cpp
	${PROJECT_SOURCE_DIR}/src/dssp_wrapper.hpp
	${PROJECT_SOURCE_DIR}/src/mkdssp.cpp
	$<TARGET_OBJECTS:dssp_library>)

target_include_directories(mkdssp PRIVATE ${CMAKE_BINARY_DIR})
target_link_libraries(mkdssp PRIVATE dssp_library cifpp::cifpp date::date libconfig::libconfig gxrio::gxrio)

if(USE_RSRC)
	mrc_target_resources(mkdssp ${CIFPP_SHARE_DIR}/mmcif_pdbx.dic)
endif()

install(TARGETS ${PROJECT_NAME}
	RUNTIME DESTINATION ${BIN_INSTALL_DIR}
)

# manual
if(UNIX)
	install(FILES doc/mkdssp.1
		DESTINATION ${CMAKE_INSTALL_DATADIR}/man/man1)
endif()

if(EXISTS "${CCP4}/html")
	install(FILES doc/mkdssp.html
		DESTINATION ${CCP4}/html)
endif()

# test
if(ENABLE_TESTING)
	# We still depend on boost for the testing headers-only library
	find_package(Boost REQUIRED)

	add_executable(dssp-unit-test ${PROJECT_SOURCE_DIR}/test/dssp-unit-test.cpp ${PROJECT_SOURCE_DIR}/src/dssp_wrapper.cpp)

	if(USE_RSRC)
		mrc_target_resources(dssp-unit-test ${CIFPP_SHARE_DIR}/mmcif_pdbx.dic)
	endif()

	target_include_directories(dssp-unit-test PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/src
		${CMAKE_CURRENT_SOURCE_DIR}/include
	)

	target_link_libraries(dssp-unit-test dssp_library cifpp::cifpp date::date Boost::boost)

	if(MSVC)
		# Specify unwind semantics so that MSVC knowns how to handle exceptions
		target_compile_options(dssp-unit-test PRIVATE /EHsc)
	endif()

	enable_testing()

	add_test(NAME dssp-unit-test COMMAND $<TARGET_FILE:dssp-unit-test> -- ${PROJECT_SOURCE_DIR}/test)
endif()

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")

# Tarball generation
set(CPACK_SOURCE_TGZ ON)
set(CPACK_SOURCE_TBZ2 OFF)
set(CPACK_SOURCE_TXZ OFF)
set(CPACK_SOURCE_TZ OFF)
set(CPACK_SOURCE_IGNORE_FILES "/data/components.cif;/build;/.vscode;/.git")
set (CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}-${PROJECT_VERSION}")
set (CPACK_SOURCE_PACKAGE_FILE_NAME ${CPACK_PACKAGE_FILE_NAME})

# NSIS options
set(CPACK_NSIS_MODIFY_PATH ON)

# configuration done, include CPack
include(CPack)
