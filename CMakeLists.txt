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
project(mkdssp VERSION 4.4.4 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CheckFunctionExists)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(FindFilesystem)
include(GenerateExportHeader)
include(CTest)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Filesystem REQUIRED)

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
write_version_header("${PROJECT_SOURCE_DIR}/src")

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
	find_package(libmcfp REQUIRED)
	find_package(cifpp 5.2.0 REQUIRED)
endif()

add_subdirectory(libdssp)

add_executable(mkdssp ${PROJECT_SOURCE_DIR}/src/mkdssp.cpp)

target_link_libraries(mkdssp PRIVATE libmcfp::libmcfp dssp::dssp)

if(USE_RSRC)
	mrc_target_resources(mkdssp
		${CIFPP_SHARE_DIR}/mmcif_pdbx.dic
		${CIFPP_SHARE_DIR}/mmcif_ddl.dic
		${PROJECT_SOURCE_DIR}/libdssp/mmcif_pdbx/dssp-extension.dic)
endif()

# Install rules
install(TARGETS ${PROJECT_NAME}
	RUNTIME DESTINATION ${BIN_INSTALL_DIR}
)

if(BUILD_DOCUMENTATION)
	# manual
	find_program(PANDOC pandoc)

	if(PANDOC)
		file(GENERATE OUTPUT ${PROJECT_BINARY_DIR}/pandoc-md.sh
			CONTENT "${PANDOC} ${PROJECT_SOURCE_DIR}/doc/mkdssp.1 -t markdown | sed -e \"s/\\\\\\' \\\\\\' (space)/  ' ' (space)/\" > ${PROJECT_SOURCE_DIR}/doc/mkdssp.md"
			FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)

		add_custom_command(OUTPUT ${PROJECT_SOURCE_DIR}/doc/mkdssp.md
			DEPENDS ${PROJECT_SOURCE_DIR}/doc/mkdssp.1
			COMMAND "${PROJECT_BINARY_DIR}/pandoc-md.sh"
			WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc/
			VERBATIM)

		add_custom_command(OUTPUT ${PROJECT_SOURCE_DIR}/doc/mkdssp.html
			DEPENDS ${PROJECT_SOURCE_DIR}/doc/mkdssp.1
			COMMAND ${PANDOC} -o mkdssp.html mkdssp.1
			WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc/
			VERBATIM)

		add_custom_command(OUTPUT ${PROJECT_SOURCE_DIR}/doc/mkdssp.pdf
			DEPENDS ${PROJECT_SOURCE_DIR}/doc/mkdssp.1
			COMMAND ${PANDOC} -o mkdssp.pdf mkdssp.1 -t html
			WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc/
			VERBATIM)

		add_custom_target(doc_files DEPENDS
			${PROJECT_SOURCE_DIR}/doc/mkdssp.html
			${PROJECT_SOURCE_DIR}/doc/mkdssp.pdf
			${PROJECT_SOURCE_DIR}/doc/mkdssp.md)
		
		add_dependencies(mkdssp doc_files)
	else()
		message(STATUS "Could not find pandoc, cannot recreate documentation files")
	endif()
endif()

if(UNIX)
	install(FILES doc/mkdssp.1
		DESTINATION ${CMAKE_INSTALL_DATADIR}/man/man1)
endif()

if(EXISTS "${CCP4}/html")
	install(FILES doc/mkdssp.html
		DESTINATION ${CCP4}/html)
endif()

# test
if(BUILD_TESTING)
	# We still depend on boost for the testing headers-only library
	find_package(Boost REQUIRED)

	add_executable(unit-test-dssp ${PROJECT_SOURCE_DIR}/test/unit-test-dssp.cpp ${PROJECT_SOURCE_DIR}/libdssp/src/dssp-io.cpp)

	if(USE_RSRC)
		mrc_target_resources(unit-test-dssp
			${CIFPP_SHARE_DIR}/mmcif_pdbx.dic ${CIFPP_SHARE_DIR}/mmcif_ddl.dic
			${CMAKE_CURRENT_SOURCE_DIR}/libdssp/mmcif_pdbx/dssp-extension.dic)
	endif()

	target_include_directories(unit-test-dssp PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/src
		${CMAKE_CURRENT_SOURCE_DIR}/include
	)

	target_link_libraries(unit-test-dssp PRIVATE dssp cifpp::cifpp Boost::boost)

	if(MSVC)
		# Specify unwind semantics so that MSVC knowns how to handle exceptions
		target_compile_options(unit-test-dssp PRIVATE /EHsc)
	endif()

	enable_testing()

	add_test(NAME unit-test-dssp COMMAND $<TARGET_FILE:unit-test-dssp> -- ${PROJECT_SOURCE_DIR}/test)
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
