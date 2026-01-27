# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright (c) 2025 NKI/AVL, Netherlands Cancer Institute
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
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

include(CPM)

if(CMAKE_VERSION GREATER_EQUAL 3.30)
	cmake_policy(SET CMP0167 NEW)
endif()

find_package(Python REQUIRED COMPONENTS Interpreter Development)

if(NOT EXISTS ${Python_LIBRARIES})
	message(FATAL_ERROR "The python library ${Python_LIBRARIES} does not seem to exist?")
endif()

# Boost 1.84 is required since it contains code required to use numpy 2.x
if(NOT WIN32)	# On Windows we need to build boost ourselves, to find the correct python lib
	find_package(Boost 1.84 QUIET COMPONENTS python)
endif()

if(WIN32 OR NOT Boost_FOUND)
	# boost is a huge project and directly downloading the 'alternate release'
	# from github is much faster than recursively cloning the repo.
	CPMAddPackage(
		NAME Boost
		VERSION 1.84.0
		URL https://github.com/boostorg/boost/releases/download/boost-1.84.0/boost-1.84.0.tar.xz
		URL_HASH SHA256=2e64e5d79a738d0fa6fb546c6e5c2bd28f88d268a2a080546f74e5ff98f29d0e
		OPTIONS
		"BOOST_ENABLE_CMAKE ON"
		"BOOST_INCLUDE_LIBRARIES python"
		"BOOST_ENABLE_PYTHON ON"
		"CMAKE_POSITION_INDEPENDENT_CODE ON"
	)
endif()

# ---------
add_library(mkdssp_module SHARED dssp-python-plugin.cpp)

target_compile_features(mkdssp_module PUBLIC cxx_std_20)

target_include_directories(mkdssp_module PRIVATE ${Python_INCLUDE_DIRS})
target_link_libraries(mkdssp_module dssp::dssp Boost::python ${Python_LIBRARIES})

set_target_properties(mkdssp_module
	PROPERTIES
	PREFIX ""
	SUFFIX ".so"
	OUTPUT_NAME mkdssp
	LINKER_LANGUAGE CXX
	LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

if(WIN32)
	# python modules use this on Windows
	set_target_properties(
		mkdssp_module
		PROPERTIES
		SUFFIX ".pyd"
		RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
	)

	install(TARGETS mkdssp_module RUNTIME DESTINATION "${Python_SITELIB}")
else()
	install(TARGETS mkdssp_module LIBRARY DESTINATION "${Python_SITELIB}")
endif()