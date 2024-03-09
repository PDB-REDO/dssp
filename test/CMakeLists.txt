find_package(Catch2 QUIET)

if(NOT Catch2_FOUND)
	include(FetchContent)

	FetchContent_Declare(
		Catch2
		GIT_REPOSITORY https://github.com/catchorg/Catch2.git
		GIT_TAG v2.13.9)

	FetchContent_MakeAvailable(Catch2)

	set(Catch2_VERSION "2.13.9")
endif()

add_executable(unit-test-dssp ${CMAKE_CURRENT_SOURCE_DIR}/unit-test-dssp.cpp ${PROJECT_SOURCE_DIR}/libdssp/src/dssp-io.cpp)

if(USE_RSRC)
	mrc_target_resources(unit-test-dssp
		${CIFPP_SHARE_DIR}/mmcif_pdbx.dic ${CIFPP_SHARE_DIR}/mmcif_ddl.dic
		${CMAKE_CURRENT_SOURCE_DIR}/../libdssp/mmcif_pdbx/dssp-extension.dic)
endif()

target_include_directories(unit-test-dssp PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/src
	${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(unit-test-dssp PRIVATE dssp cifpp::cifpp Catch2::Catch2)

if(MSVC)
	# Specify unwind semantics so that MSVC knowns how to handle exceptions
	target_compile_options(unit-test-dssp PRIVATE /EHsc)
endif()

if(${Catch2_VERSION} VERSION_GREATER_EQUAL 3.0.0)
	target_compile_definitions(unit-test-dssp PUBLIC CATCH22=0)
else()
	target_compile_definitions(unit-test-dssp PUBLIC CATCH22=1)
endif()

add_test(NAME unit-test-dssp COMMAND $<TARGET_FILE:unit-test-dssp> --data-dir ${CMAKE_CURRENT_SOURCE_DIR})
