#
# HO-CXX
# (c)2020-1 Applied Scientific Research, Inc.
#
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )
PROJECT( HO-CXX )
ENABLE_LANGUAGE (CXX)

# User-serviceable variables here
SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE "Release")
ENDIF ()
SET (CMAKE_BUILD_TYPE "Release" CACHE STRING "Select which configuration to build")
SET (USE_OMP FALSE CACHE BOOL "Use OpenMP multithreading")
SET (CMAKE_VERBOSE_MAKEFILE on)
SET (CMAKE_EXPORT_COMPILE_COMMANDS on)

# Define compiler-specific options
IF (CMAKE_COMPILER_IS_GNUCXX)
  SET (CMAKE_CXX_FLAGS "-Wall -Wformat -std=c++11")
  SET (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb3")
  SET (CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native")
  SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -march=native -g -ggdb3")
ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  SET (CMAKE_CXX_FLAGS "-Wall -Wformat -std=c++17")
  SET (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb")
  SET (CMAKE_CXX_FLAGS_RELEASE "-O3")
  SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -ggdb")
ELSEIF (MSVC)
  SET (CMAKE_CXX_FLAGS "/std:c++17 /EHsc /D_USE_MATH_DEFINES /DNOMINMAX")
  SET (CMAKE_CXX_FLAGS_DEBUG "/Zi")
  SET (CMAKE_CXX_FLAGS_RELEASE "/O2")
ENDIF ()

# suppress the FindOpenMP warnings
cmake_policy(SET CMP0054 NEW)

IF( USE_OMP )
  FIND_PACKAGE( OpenMP REQUIRED )
  IF( OPENMP_FOUND )
    #message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  ENDIF()
ELSE()
  # suppress the warning
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
ENDIF()

IF( WIN32 )
ELSEIF( APPLE )
ELSEIF( UNIX )
  INCLUDE_DIRECTORIES ( "/usr/include/eigen3" )
  INCLUDE_DIRECTORIES ( "/home/mstock/opt/amgcl-1.4.0" )
ENDIF ()

FILE(GLOB HEADERS "*.h" )
SET(SOURCES ${HEADERS}
            "HO_cpp_2D.cpp"
            "preprocess.cpp" )

# Define targets
ADD_LIBRARY ( ${PROJECT_NAME} STATIC ${SOURCES} )

#ADD_EXECUTABLE( "${PROJECT_NAME}.bin" ${SOURCES} "main.cpp" )
ADD_EXECUTABLE( "${PROJECT_NAME}.bin" "main.cpp" )
TARGET_LINK_LIBRARIES( "${PROJECT_NAME}.bin" LINK_PUBLIC ${PROJECT_NAME} ${LIBS} )
