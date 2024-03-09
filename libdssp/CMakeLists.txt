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
	"${CMAKE_CURRENT_BINARY_DIR}/dssp/dsspConfigVersion.cmake"
	DESTINATION lib/cmake/dssp
	COMPONENT Devel
)

# set(dssp_MAJOR_VERSION ${CMAKE_PROJECT_VERSION_MAJOR})
# set_target_properties(dssp PROPERTIES
# 	VERSION ${PROJECT_VERSION}
# 	SOVERSION ${dssp_MAJOR_VERSION}
# 	INTERFACE_dssp_MAJOR_VERSION ${dssp_MAJOR_VERSION})

# set_property(TARGET dssp APPEND PROPERTY
# 	COMPATIBLE_INTERFACE_STRING dssp_MAJOR_VERSION
# )

# write_basic_package_version_file(
# 	"${CMAKE_CURRENT_BINARY_DIR}/dssp/dsspConfigVersion.cmake"
# 	VERSION ${PROJECT_VERSION}
# 	COMPATIBILITY AnyNewerVersion
# )

# # pkgconfig support
# set(prefix ${CMAKE_INSTALL_PREFIX})
# set(exec_prefix ${CMAKE_INSTALL_PREFIX})
# set(libdir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
# set(includedir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})

# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libdssp.pc.in
# 	${CMAKE_CURRENT_BINARY_DIR}/libdssp.pc.in @ONLY)
# file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libdssp.pc
# 	INPUT ${CMAKE_CURRENT_BINARY_DIR}/libdssp.pc.in)
# install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libdssp.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

set(CIFPP_DATA_DIR "${CMAKE_INSTALL_FULL_DATADIR}/libcifpp")
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/mmcif_pdbx/dssp-extension.dic"
        DESTINATION ${CIFPP_DATA_DIR})
