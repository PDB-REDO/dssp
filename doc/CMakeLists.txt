# manual 
find_program(PANDOC pandoc)

if(PANDOC)
	file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pandoc-md.sh
		CONTENT "${PANDOC} ${CMAKE_CURRENT_SOURCE_DIR}/mkdssp.1 -t markdown | sed -e \"s/\\\\\\' \\\\\\' (space)/  ' ' (space)/\" > ${CMAKE_CURRENT_SOURCE_DIR}/mkdssp.md"
		FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)

	add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/mkdssp.md
		DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/mkdssp.1
		COMMAND "${CMAKE_CURRENT_BINARY_DIR}/pandoc-md.sh"
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
		VERBATIM)

	add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/mkdssp.html
		DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/mkdssp.1
		COMMAND ${PANDOC} -o mkdssp.html mkdssp.1
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
		VERBATIM)

	add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/mkdssp.pdf
		DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/mkdssp.1
		COMMAND ${PANDOC} -o mkdssp.pdf mkdssp.1 -t html
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
		VERBATIM)

	add_custom_target(doc_files DEPENDS
		${CMAKE_CURRENT_SOURCE_DIR}/mkdssp.html
		${CMAKE_CURRENT_SOURCE_DIR}/mkdssp.pdf
		${CMAKE_CURRENT_SOURCE_DIR}/mkdssp.md)

	add_dependencies(mkdssp doc_files)
else()
	message(STATUS "Could not find pandoc, cannot recreate documentation files")
endif()
