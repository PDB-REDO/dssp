CPMFindPackage(
	NAME Catch2 3
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG v3.4.0
	EXCLUDE_FROM_ALL YES)

add_executable(unit-test-dssp ${CMAKE_CURRENT_SOURCE_DIR}/unit-test-dssp.cpp ${PROJECT_SOURCE_DIR}/libdssp/src/dssp-io.cpp)

target_include_directories(unit-test-dssp PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/src
	${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(unit-test-dssp PRIVATE dssp cifpp::cifpp Catch2::Catch2)

if(MSVC)
	# Specify unwind semantics so that MSVC knowns how to handle exceptions
	target_compile_options(unit-test-dssp PRIVATE /EHsc)
endif()

add_test(NAME unit-test-dssp COMMAND $<TARGET_FILE:unit-test-dssp>
	--data-dir ${CMAKE_CURRENT_SOURCE_DIR}
	--rsrc-dir ${CMAKE_CURRENT_SOURCE_DIR}/../libdssp/mmcif_pdbx/)

if(BUILD_PYTHON_MODULE)
	find_package(Python REQUIRED Interpreter)

	add_test(NAME python_module COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/test-python.py")
	set_tests_properties(python_module PROPERTIES
		ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:mkdssp_module>"
		WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
endif()