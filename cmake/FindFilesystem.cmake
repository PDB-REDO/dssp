# Simplistic reimplementation of https://github.com/vector-of-bool/CMakeCM/blob/master/modules/FindFilesystem.cmake

if(TARGET std::filesystem)
	return()
endif()

cmake_minimum_required(VERSION 3.10)

include(CMakePushCheckState)
include(CheckIncludeFileCXX)
include(CheckCXXSourceCompiles)

cmake_push_check_state()

set(CMAKE_CXX_STANDARD 17)

check_include_file_cxx("filesystem" _CXX_FILESYSTEM_HAVE_HEADER)
mark_as_advanced(_CXX_FILESYSTEM_HAVE_HEADER)

set(code [[
#include <cstdlib>
#include <filesystem>

int main() {
	auto cwd = std::filesystem::current_path();
	return EXIT_SUCCESS;
}
]])

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL 8.4.0)
	# >> https://stackoverflow.com/questions/63902528/program-crashes-when-filesystempath-is-destroyed
	set(CXX_FILESYSTEM_NO_LINK_NEEDED 0)
else()
	# Check a simple filesystem program without any linker flags
	check_cxx_source_compiles("${code}" CXX_FILESYSTEM_NO_LINK_NEEDED)
endif()

if(CXX_FILESYSTEM_NO_LINK_NEEDED)
	set(_found 1)
else()
	set(prev_libraries ${CMAKE_REQUIRED_LIBRARIES})
	# Add the libstdc++ flag
	set(CMAKE_REQUIRED_LIBRARIES ${prev_libraries} -lstdc++fs)
	check_cxx_source_compiles("${code}" CXX_FILESYSTEM_STDCPPFS_NEEDED)
	set(_found ${CXX_FILESYSTEM_STDCPPFS_NEEDED})
	if(NOT CXX_FILESYSTEM_STDCPPFS_NEEDED)
		# Try the libc++ flag
		set(CMAKE_REQUIRED_LIBRARIES ${prev_libraries} -lc++fs)
		check_cxx_source_compiles("${code}" CXX_FILESYSTEM_CPPFS_NEEDED)
		set(_found ${CXX_FILESYSTEM_CPPFS_NEEDED})
	endif()
endif()

if(_found)
	add_library(std::filesystem INTERFACE IMPORTED)
	set_property(TARGET std::filesystem APPEND PROPERTY INTERFACE_COMPILE_FEATURES cxx_std_17)

	if(CXX_FILESYSTEM_NO_LINK_NEEDED)
		# Nothing to add...
	elseif(CXX_FILESYSTEM_STDCPPFS_NEEDED)
		set_target_properties(std::filesystem PROPERTIES IMPORTED_LIBNAME stdc++fs)
	elseif(CXX_FILESYSTEM_CPPFS_NEEDED)
		set_target_properties(std::filesystem PROPERTIES IMPORTED_LIBNAME c++fs)
	endif()
endif()

cmake_pop_check_state()

set(Filesystem_FOUND ${_found} CACHE BOOL "TRUE if we can run a program using std::filesystem" FORCE)

if(Filesystem_FIND_REQUIRED AND NOT Filesystem_FOUND)
    message(FATAL_ERROR "Cannot run simple program using std::filesystem")
endif()

