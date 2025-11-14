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

cmake_minimum_required(VERSION 3.23)

# set the project name
project(mkdssp VERSION 4.5.5 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(CPM)
include(CTest)
include(VersionString)

# The default build type must be set before project()
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers")
elseif(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

# Optionally build a version to be installed inside CCP4
option(BUILD_FOR_CCP4 "Build a version to be installed in CCP4" OFF)
option(BUILD_DOCUMENTATION "Generate the documentation files using pandoc" OFF)
option(BUILD_PYTHON_MODULE "Build an experimental Python module" OFF)
option(INSTALL_LIBRARY "Install the libdssp library" OFF)

if(BUILD_FOR_CCP4)
	if("$ENV{CCP4}" STREQUAL "" OR NOT EXISTS $ENV{CCP4})
		message(FATAL_ERROR "dssp: A CCP4 built was requested but CCP4 was not sourced")
	else()
		list(APPEND CMAKE_MODULE_PATH "$ENV{CCP4}")
		list(APPEND CMAKE_PREFIX_PATH "$ENV{CCP4}")
		set(CMAKE_INSTALL_PREFIX "$ENV{CCP4}")

		if(WIN32)
			set(BUILD_SHARED_LIBS ON)
		endif()
	endif()
endif()

if(WIN32)
	if(${CMAKE_SYSTEM_VERSION} GREATER_EQUAL 10) # Windows 10
		add_definitions(-D _WIN32_WINNT=0x0A00)
	elseif(${CMAKE_SYSTEM_VERSION} EQUAL 6.3) # Windows 8.1
		add_definitions(-D _WIN32_WINNT=0x0603)
	elseif(${CMAKE_SYSTEM_VERSION} EQUAL 6.2) # Windows 8
		add_definitions(-D _WIN32_WINNT=0x0602)
	elseif(${CMAKE_SYSTEM_VERSION} EQUAL 6.1) # Windows 7
		add_definitions(-D _WIN32_WINNT=0x0601)
	elseif(${CMAKE_SYSTEM_VERSION} EQUAL 6.0) # Windows Vista
		add_definitions(-D _WIN32_WINNT=0x0600)
	else() # Windows XP (5.1)
		add_definitions(-D _WIN32_WINNT=0x0501)
	endif()

	# Man, this is 2024 we're living in...
	add_definitions(-DNOMINMAX)

	# We do not want to write an export file for all our symbols...
	set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

if(MSVC)
	# make msvc standards compliant...
	add_compile_options(/permissive- /bigobj)
	add_link_options(/NODEFAULTLIB:library)

	# This is dubious...
	if(BUILD_SHARED_LIBS)
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
	else()
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
	endif()
endif()

# Create a revision file, containing the current git version info
write_version_header("${CMAKE_CURRENT_SOURCE_DIR}/src")

# Optionally use mrc to create resources
find_package(mrc QUIET)

if(MRC_FOUND)
	option(USE_RSRC "Use mrc to create resources" ON)
else()
	message(STATUS "dssp: Not using resources since mrc was not found")
endif()

set(CMAKE_THREAD_PREFER_PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG)
find_package(Threads)

# Find fmt library for cifpp
find_package(fmt QUIET)
if(NOT fmt_FOUND)
        message(WARNING "fmt library not found - cifpp may fail to link")
endif()

CPMFindPackage(
	NAME mcfp
	GIT_REPOSITORY https://github.com/mhekkel/libmcfp
	GIT_TAG v1.4.2
	EXCLUDE_FROM_ALL YES)

# CPMFindPackage does not pass on variables set with set_and_check
# Working around this here
find_package(cifpp 9 QUIET)

if(NOT cifpp_FOUND)
	CPMAddPackage(
		NAME cifpp
		GIT_REPOSITORY https://github.com/pdb-redo/libcifpp.git
		GIT_TAG v9.0.0
		EXCLUDE_FROM_ALL YES)
endif()

if(TARGET cifpp)
	get_target_property(CIFPP_SOURCE_DIR cifpp SOURCE_DIR)
	set(CIFPP_SHARE_DIR share/libcifpp)
	set(CIFPP_DATA_DIR ${CIFPP_SOURCE_DIR}/rsrc)
elseif(DEFINED CIFPP_SHARE_DIR)
	set(CIFPP_DATA_DIR ${CIFPP_SHARE_DIR})
else()
	message(FATAL_ERROR "dssp: The CIFPP_SHARE_DIR variable is not found in the cifpp configuration files")
endif()

if(INSTALL_LIBRARY)
	add_subdirectory(libdssp)
else()
	add_subdirectory(libdssp EXCLUDE_FROM_ALL)
endif()

add_executable(mkdssp ${CMAKE_CURRENT_SOURCE_DIR}/src/mkdssp.cpp)

target_link_libraries(mkdssp PRIVATE mcfp::mcfp dssp::dssp)

if(USE_RSRC)
	mrc_target_resources(mkdssp
		${CIFPP_DATA_DIR}/mmcif_pdbx.dic
		${CMAKE_CURRENT_SOURCE_DIR}/libdssp/mmcif_pdbx/dssp-extension.dic)
else()
	get_target_property(LIBDSSP_SOURCE_DIR dssp SOURCE_DIR)

	install(FILES
		${CMAKE_CURRENT_SOURCE_DIR}/libdssp/mmcif_pdbx/dssp-extension.dic
		DESTINATION ${CIFPP_SHARE_DIR})

	if(TARGET cifpp)
		install(FILES ${CIFPP_DATA_DIR}/mmcif_pdbx.dic
			DESTINATION ${CIFPP_SHARE_DIR})
	endif()
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
if(BUILD_TESTING AND PROJECT_IS_TOP_LEVEL)
	add_subdirectory(test)
endif()

# Python
if(BUILD_PYTHON_MODULE)
	add_subdirectory(python-module)
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
