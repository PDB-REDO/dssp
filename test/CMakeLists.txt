# We're using the older version 2 of Catch2
FetchContent_Declare(
	Catch2
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG v2.13.9)

FetchContent_MakeAvailable(Catch2)

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

target_compile_definitions(unit-test-dssp PUBLIC CATCH22=1)

add_test(NAME unit-test-dssp COMMAND $<TARGET_FILE:unit-test-dssp> --data-dir ${CMAKE_CURRENT_SOURCE_DIR})
