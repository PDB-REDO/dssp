cmake_minimum_required(VERSION 3.10)

# set the project name
project(mkdssp VERSION 4.0.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set bindir, if not use -DBIN_INSTALL_DIR
if(NOT BIN_INSTALL_DIR)
  if(CMAKE_INSTALL_BINDIR)
    set(BIN_INSTALL_DIR ${CMAKE_INSTALL_BINDIR})
  else(CMAKE_INSTALL_BINDIR)
    set(BIN_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/bin")
  endif(CMAKE_INSTALL_BINDIR)
endif(NOT BIN_INSTALL_DIR)

# Set libdir, if not use -DLIB_INSTALL_DIR
if(NOT LIB_INSTALL_DIR)
  if(CMAKE_INSTALL_LIBDIR)
    set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})
  else(CMAKE_INSTALL_LIBDIR)
    set(LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")
  endif(CMAKE_INSTALL_LIBDIR)
endif(NOT LIB_INSTALL_DIR)

# Set includedir, if not use -DINCLUDE_INSTALL_DIR
if(NOT INCLUDE_INSTALL_DIR)
  if(CMAKE_INSTALL_INCLUDEDIR)
    set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
  else(CMAKE_INSTALL_INCLUDEDIR)
    set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include")
  endif(CMAKE_INSTALL_INCLUDEDIR)
endif(NOT INCLUDE_INSTALL_DIR)

# Set sharedir, if not use -DSHARE_INSTALL_DIR
if(NOT SHARE_INSTALL_DIR)
  if(CMAKE_INSTALL_DATADIR)
    set(SHARE_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}")
  else(CMAKE_INSTALL_DATADIR)
    set(SHARE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/share")
  endif(CMAKE_INSTALL_DATADIR)
endif(NOT SHARE_INSTALL_DIR)

set (Boost_DETAILED_FAILURE_MSG ON)
# set (BOOST_ROOT ${PROJECT_SOURCE_DIR}/../boost_1_75_0)
# set (Boost_COMPILER "-vc")
# set (Boost_USE_STATIC_RUNTIME ON)  

find_package(Boost 1.73.0 REQUIRED COMPONENTS system iostreams regex date_time program_options)

find_package(cifpp 1.0 REQUIRED)

set(CMAKE_THREAD_PREFER_PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG)
find_package(Threads)

include_directories(
    ${PROJECT_SOURCE_DIR}/src
	PUBLIC ${cifpp_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS}
)

add_executable(mkdssp ${PROJECT_SOURCE_DIR}/src/mkdssp.cpp)

install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION ${BIN_INSTALL_DIR}
)

target_link_libraries(${PROJECT_NAME} CifPP::cifpp ${Boost_LIBRARIES} Threads::Threads)

if(MSVC)
    # make msvc standards compliant...
    target_compile_options(${PROJECT_NAME} PRIVATE /permissive-)
else()
	find_library(Z z)
	find_library(BZ2 bz2)

	target_link_libraries(${PROJECT_NAME} ${Z} ${BZ2})
endif()

