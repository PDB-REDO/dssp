# SPDX-License-Identifier: BSD-2-Clause
# 
# Copyright (c) 2025 NKI/AVL, Netherlands Cancer Institute
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

cmake_minimum_required(VERSION 3.23)

file(READ VERSION VER_FILE_CONTENT)
string(STRIP ${VER_FILE_CONTENT} VER_FILE_CONTENT)

project(libdssp VERSION ${VER_FILE_CONTENT} LANGUAGES CXX)

include(CMakePackageConfigHelpers)

# We do not want to write an export file for all our symbols...
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Create a revision file, containing the current git version info
write_version_header("${CMAKE_CURRENT_SOURCE_DIR}/src" LIB_NAME libdssp)

# The DSSP code is in a separate library, optionally to be used by others
add_library(dssp)
add_library(dssp::dssp ALIAS dssp)

target_sources(dssp
	PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/src/dssp.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/dssp-io.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/dssp-io.hpp
	PUBLIC
	FILE_SET dssp_header TYPE HEADERS
	BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
	FILES include/dssp.hpp
)

target_link_libraries(dssp PUBLIC cifpp::cifpp)

set_target_properties(dssp PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(dssp PUBLIC cxx_std_17)

target_include_directories(dssp
	PUBLIC
	"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include;${CMAKE_CURRENT_SOURCE_DIR}/../src>"
	"$<INSTALL_INTERFACE:include>"
)

# Install rules
install(TARGETS dssp
	EXPORT dssp
	FILE_SET dssp_header DESTINATION include)

if(MSVC AND BUILD_SHARED_LIBS)
	install(
		FILES $<TARGET_PDB_FILE:dssp>
		DESTINATION "lib"
		OPTIONAL)
endif()

install(EXPORT dssp
	FILE "dssp-targets.cmake"
	NAMESPACE dssp::
	DESTINATION lib/cmake/dssp
)

set(CONFIG_TEMPLATE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/dssp-config.cmake.in)

configure_package_config_file(
	${CONFIG_TEMPLATE_FILE}
	${CMAKE_CURRENT_BINARY_DIR}/dssp/dssp-config.cmake
	INSTALL_DESTINATION lib/cmake/dssp
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
	"${CMAKE_CURRENT_BINARY_DIR}/dssp/dssp-config-version.cmake"
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY AnyNewerVersion
)

install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/dssp/dssp-config.cmake"
	"${CMAKE_CURRENT_BINARY_DIR}/dssp/dssp-config-version.cmake"
	DESTINATION lib/cmake/dssp
	COMPONENT Devel
)


# set(CIFPP_DATA_DIR "${CMAKE_INSTALL_FULL_DATADIR}/libcifpp")
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/mmcif_pdbx/dssp-extension.dic"
        DESTINATION share/libcifpp)
