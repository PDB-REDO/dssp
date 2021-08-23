cmake_minimum_required(VERSION 3.15)

# set the project name
project(mkdssp VERSION 4.0.0 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

enable_testing()

include(GNUInstallDirs)
include(CheckFunctionExists)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(Dart)
include(GenerateExportHeader)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W4")
endif()

if(MSVC)
    # make msvc standards compliant...
    add_compile_options(/permissive-)

	macro(get_WIN32_WINNT version)
		if (WIN32 AND CMAKE_SYSTEM_VERSION)
			set(ver ${CMAKE_SYSTEM_VERSION})
			string(REPLACE "." "" ver ${ver})
			string(REGEX REPLACE "([0-9])" "0\\1" ver ${ver})

			set(${version} "0x${ver}")
		endif()
	endmacro()

	get_WIN32_WINNT(ver)
	add_definitions(-D_WIN32_WINNT=${ver})

	# On Windows, do not install in the system location by default
	if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
		message(WARNING "The library and auxiliary files will be installed in $ENV{LOCALAPPDATA}/${PROJECT_NAME}")
		set(CMAKE_INSTALL_PREFIX "$ENV{LOCALAPPDATA}/${PROJECT_NAME}" CACHE PATH "..." FORCE)
	endif()

	# Find out the processor type for the target
	if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "AMD64")
		set(COFF_TYPE "x64")
	elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "i386")
		set(COFF_TYPE "x86")
	elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ARM64")
		set(COFF_TYPE "arm64")
	else()
		message(FATAL_ERROR "Unsupported or unknown processor type ${CMAKE_SYSTEM_PROCESSOR}")
	endif()	

	set(COFF_SPEC "--coff=${COFF_TYPE}")
endif()

# Create a revision file, containing the current git version info

find_package(Git)
if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
	include(GetGitRevisionDescription)
    get_git_head_revision(REFSPEC COMMITHASH)

	# Generate our own version string
	git_describe_working_tree(BUILD_VERSION_STRING --match=build --dirty)
else()
    message(WARNING "Git not found, cannot set version info")

    SET(BUILD_VERSION_STRING "unknown")
endif()

# generate version.h
include_directories(${CMAKE_BINARY_DIR} PRIVATE)
string(TIMESTAMP BUILD_DATE_TIME "%Y-%m-%d" UTC)
configure_file("${CMAKE_SOURCE_DIR}/src/revision.hpp.in" "${CMAKE_BINARY_DIR}/revision.hpp" @ONLY)

# Optionally use mrc to create resources
find_program(MRC mrc HINTS "$ENV{LOCALAPPDATA}/mrc" "${CMAKE_INSTALL_PREFIX}/../mrc" "/usr/local/bin")

if(MRC)
	option(USE_RSRC "Use mrc to create resources" ON)
else()
	message(WARNING "Not using resources since mrc was not found")
endif()

if(USE_RSRC STREQUAL "ON")
	set(USE_RSRC 1)

	message("Using resources compiled with ${MRC}")
	add_compile_definitions(USE_RSRC)
endif()

set(CMAKE_THREAD_PREFER_PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG)
find_package(Threads REQUIRED)

set (Boost_DETAILED_FAILURE_MSG ON)
find_package(Boost 1.71.0 REQUIRED COMPONENTS program_options system iostreams regex date_time)
# extra diagnostic -- helpful for problem with FindBoost.cmake
message(STATUS "Boost headers in: ${Boost_INCLUDE_DIR}")
message(STATUS "Boost libraries in: ${Boost_LIBRARY_DIRS}")

find_package(ZLIB)
find_package(BZip2)
find_package(cifpp 1.0 REQUIRED)

if(USE_RSRC)
	add_compile_definitions("DATA_DIR=\"${CIFPP_SHARE_DIR}\"")
endif()

include_directories(${Boost_INCLUDE_DIR} cifpp::cifpp ${CMAKE_SOURCE_DIR}/include)
link_libraries(${Boost_LIBRARIES} cifpp::cifpp ${CMAKE_THREAD_LIBS_INIT})

if (ZLIB_FOUND)
	link_libraries(ZLIB::ZLIB)
endif()

if (BZIP2_FOUND)
	link_libraries(BZip2::BZip2)
endif()

include_directories(
    ${PROJECT_SOURCE_DIR}/src
	PUBLIC ${cifpp_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS}
)

add_executable(mkdssp
	${PROJECT_SOURCE_DIR}/src/mkdssp.cpp)

install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION ${BIN_INSTALL_DIR}
)
