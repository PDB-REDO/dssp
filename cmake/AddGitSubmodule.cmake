cmake_minimum_required(VERSION 3.16..3.19)

function(add_git_submodule dir)
	# add a Git submodule directory to CMake, assuming the
	# Git submodule directory is a CMake project.
	#
	# Usage: in CMakeLists.txt
	#
	# include(AddGitSubmodule.cmake)
	# add_git_submodule(mysubmod_dir)
	find_package(Git REQUIRED)

	if(NOT EXISTS ${dir}/CMakeLists.txt)
		if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.19)
			execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive -- ${dir}
				WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
				COMMAND_ERROR_IS_FATAL ANY)
		else()
			execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive -- ${dir}
				WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
		endif()
	endif()

	set(ENABLE_TESTING OFF)

	add_subdirectory(${dir} ${ARGV})
endfunction(add_git_submodule)