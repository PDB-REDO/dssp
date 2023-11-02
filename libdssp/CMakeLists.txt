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
project(libdssp VERSION 4.4.4 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CheckFunctionExists)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(FindFilesystem)

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
endif()

# We do not want to write an export file for all our symbols...
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Create a revision file, containing the current git version info
include(VersionString)
write_version_header("${PROJECT_SOURCE_DIR}/src" LIB_NAME "libdssp")

# Optionally use mrc to create resources
find_package(Mrc QUIET)

set(CMAKE_THREAD_PREFER_PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG)
find_package(Threads)

if(NOT PDB_REDO_META)
	find_package(cifpp 5.1.0 REQUIRED)
endif()

# The DSSP code is in a separate library, optionally to be used by others
add_library(dssp
	${PROJECT_SOURCE_DIR}/src/dssp.cpp
	${PROJECT_SOURCE_DIR}/include/dssp.hpp
	${PROJECT_SOURCE_DIR}/src/dssp-io.cpp
	${PROJECT_SOURCE_DIR}/src/dssp-io.hpp
)

add_library(dssp::dssp ALIAS dssp)

target_link_libraries(dssp PUBLIC cifpp::cifpp)

set_target_properties(dssp PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(dssp
	PUBLIC
	"$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
	"$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

# Install rules
install(TARGETS dssp
	EXPORT dsspTargets
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(MSVC AND BUILD_SHARED_LIBS)
	install(
		FILES $<TARGET_PDB_FILE:dssp>
		DESTINATION ${CMAKE_INSTALL_LIBDIR}
		OPTIONAL)
endif()

install(EXPORT dsspTargets
	FILE "dsspTargets.cmake"
	NAMESPACE dssp::
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dssp
)

install(
	FILES include/dssp.hpp
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	COMPONENT Devel
)

set(CONFIG_TEMPLATE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/dsspConfig.cmake.in)

configure_package_config_file(
	${CONFIG_TEMPLATE_FILE}
	${CMAKE_CURRENT_BINARY_DIR}/dssp/dsspConfig.cmake
	INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dssp
)

install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/dssp/dsspConfig.cmake"
	"${CMAKE_CURRENT_BINARY_DIR}/dssp/dsspConfigVersion.cmake"
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dssp
	COMPONENT Devel
)

set(dssp_MAJOR_VERSION ${CMAKE_PROJECT_VERSION_MAJOR})
set_target_properties(dssp PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION ${dssp_MAJOR_VERSION}
	INTERFACE_dssp_MAJOR_VERSION ${dssp_MAJOR_VERSION})

set_property(TARGET dssp APPEND PROPERTY
	COMPATIBLE_INTERFACE_STRING dssp_MAJOR_VERSION
)

write_basic_package_version_file(
	"${CMAKE_CURRENT_BINARY_DIR}/dssp/dsspConfigVersion.cmake"
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY AnyNewerVersion
)

# pkgconfig support
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix ${CMAKE_INSTALL_PREFIX})
set(libdir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(includedir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libdssp.pc.in
	${CMAKE_CURRENT_BINARY_DIR}/libdssp.pc.in @ONLY)
file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libdssp.pc
	INPUT ${CMAKE_CURRENT_BINARY_DIR}/libdssp.pc.in)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libdssp.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

set(CIFPP_DATA_DIR "${CMAKE_INSTALL_FULL_DATADIR}/libcifpp")
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/mmcif_pdbx/dssp-extension.dic"
        DESTINATION ${CIFPP_DATA_DIR})
