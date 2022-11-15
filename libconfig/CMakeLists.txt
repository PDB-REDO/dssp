# Copyright Maarten L. Hekkelman, 2022
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 3.16)

# set the project name
project(libconfig VERSION 1.2.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(Dart)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17 CACHE STRING "The minimum version of C++ required for this library")
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_TESTING "Build the unit test applications" OFF)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers")
elseif(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
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

add_library(libconfig INTERFACE)
add_library(libconfig::libconfig ALIAS libconfig)

target_include_directories(libconfig INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
)

# adding header sources just helps IDEs
target_sources(libconfig INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>$<INSTALL_INTERFACE:include>/cfg.hpp
)

set_target_properties(libconfig PROPERTIES PUBLIC_HEADER include/cfg.hpp)

# installation
set(version_config "${CMAKE_CURRENT_BINARY_DIR}/libconfigConfigVersion.cmake")

set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})

write_basic_package_version_file("${version_config}"
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY SameMajorVersion)

install(TARGETS libconfig
	EXPORT libconfigConfig
	PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(
	DIRECTORY include/cfg
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	COMPONENT Devel
)

export(TARGETS libconfig NAMESPACE libconfig:: FILE libconfigTargets.cmake)

if(WIN32 AND NOT CYGWIN)
	set(CONFIG_LOC CMake)
else()
	set(CONFIG_LOC "${CMAKE_INSTALL_LIBDIR}/cmake/libconfig")
endif()

configure_package_config_file(
	${PROJECT_SOURCE_DIR}/cmake/libconfigConfig.cmake.in
	${PROJECT_SOURCE_DIR}/cmake/libconfigConfig.cmake
	INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libconfig
	PATH_VARS INCLUDE_INSTALL_DIR
)

install(EXPORT libconfigConfig
	FILE libconfigTargets.cmake
	NAMESPACE libconfig::
	DESTINATION ${CONFIG_LOC})

install(
	FILES cmake/libconfigConfig.cmake "${version_config}"
	DESTINATION ${CONFIG_LOC})

if(ENABLE_TESTING)
	enable_testing()

	add_executable(libconfig-unit-test ${PROJECT_SOURCE_DIR}/test/unit-test.cpp)

	target_link_libraries(libconfig-unit-test libconfig::libconfig)

	if(MSVC)
		# Specify unwind semantics so that MSVC knowns how to handle exceptions
		target_compile_options(libconfig-unit-test PRIVATE /EHsc)
	endif()

	add_test(NAME libconfig-unit-test
		COMMAND $<TARGET_FILE:libconfig-unit-test> -- ${PROJECT_SOURCE_DIR}/test)
endif()
