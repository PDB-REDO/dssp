# We do not want to write an export file for all our symbols...
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

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
	"$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

# Install rules
install(TARGETS dssp
	EXPORT dssp
	FILE_SET dssp_header DESTINATION "include")

if(MSVC AND BUILD_SHARED_LIBS)
	install(
		FILES $<TARGET_PDB_FILE:dssp>
		DESTINATION "lib"
		OPTIONAL)
endif()

install(EXPORT dssp
	FILE "dsspTargets.cmake"
	NAMESPACE dssp::
	DESTINATION lib/cmake/dssp
)

set(CONFIG_TEMPLATE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/dsspConfig.cmake.in)

configure_package_config_file(
	${CONFIG_TEMPLATE_FILE}
	${CMAKE_CURRENT_BINARY_DIR}/dssp/dsspConfig.cmake
	INSTALL_DESTINATION lib/cmake/dssp
)

install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/dssp/dsspConfig.cmake"
	DESTINATION lib/cmake/dssp
	COMPONENT Devel
)

set(CIFPP_DATA_DIR share/libcifpp)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/mmcif_pdbx/dssp-extension.dic"
        DESTINATION ${CIFPP_DATA_DIR})
