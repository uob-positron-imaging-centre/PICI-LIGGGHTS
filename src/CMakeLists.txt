# CMake-based compilation
PROJECT(LIGGGHTS C CXX)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

add_compile_options(-Wall)
add_compile_options(-Wextra)
add_compile_options(-Wdisabled-optimization)

SET(LIGGGHTS_MAJOR_VERSION 3)
SET(LIGGGHTS_MINOR_VERSION 7)
SET(LIGGGHTS_PATCH_VERSION 0)
SET(LIGGGHTS_VERSION ${LIGGGHTS_MAJOR_VERSION}.${LIGGGHTS_MINOR_VERSION}.${LIGGGHTS_PATCH_VERSION})

IF (NOT CMAKE_CXX_FLAGS)
    SET(CMAKE_CXX_FLAGS "-O2 -std=c++11 -fPIC")
ENDIF()


# Inspired from the lethe-cfd library
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")


# Collect source and header files needed
FILE(GLOB SOURCES *.cpp)
FILE(GLOB HEADERS *.h)

#LIST(REMOVE_ITEM SOURCES main.cpp)

IF(WIN32)
  #include (GenerateExportHeader)
  INCLUDE_DIRECTORIES(WINDOWS/extra)
  ADD_DEFINITIONS(-D_USE_MATH_DEFINES -DNOMINMAX)

  IF(MSVC)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
  ENDIF()

  FIND_PACKAGE(Cygwin REQUIRED)

  ADD_CUSTOM_TARGET(GenerateHeaders
           COMMAND GenerateHeaders.bat ${CYGWIN_INSTALL_PATH}
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/WINDOWS)
  # enable parallelization
  foreach(flag_var
          CMAKE_C_FLAGS_RELEASE CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO
          CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
            set(${flag_var} "${${flag_var}} /MP4")
  endforeach()
ELSE()
  ADD_CUSTOM_TARGET(GenerateHeaders
           COMMAND sh Make.sh style
           COMMAND sh Make.sh models
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

ADD_LIBRARY(liggghts_obj OBJECT ${SOURCES})
ADD_DEPENDENCIES(liggghts_obj GenerateHeaders)

ADD_LIBRARY(liggghts_static STATIC $<TARGET_OBJECTS:liggghts_obj>)
SET_TARGET_PROPERTIES(liggghts_static PROPERTIES OUTPUT_NAME liggghts)
SET_PROPERTY(TARGET liggghts_static PROPERTY PUBLIC_HEADER ${HEADERS})

ADD_LIBRARY(liggghts_shared SHARED $<TARGET_OBJECTS:liggghts_obj>)
SET_TARGET_PROPERTIES(liggghts_shared PROPERTIES OUTPUT_NAME liggghts)
SET_PROPERTY(TARGET liggghts_shared PROPERTY PUBLIC_HEADER ${HEADERS})

#ADD_EXECUTABLE(liggghts_bin main.cpp)
#TARGET_LINK_LIBRARIES(liggghts_bin liggghts)
ADD_EXECUTABLE(liggghts_bin $<TARGET_OBJECTS:liggghts_obj>)
SET_TARGET_PROPERTIES(liggghts_bin PROPERTIES OUTPUT_NAME liggghts)

#=======================================
ADD_DEFINITIONS(-DSUPERQUADRIC_ACTIVE_FLAG)
ADD_DEFINITIONS(-DNONSPHERICAL_ACTIVE_FLAG)
FIND_PACKAGE(VTK NO_MODULE)

IF(VTK_FOUND)
  ADD_DEFINITIONS(-DLAMMPS_VTK)

  if (VTK_MAJOR_VERSION VERSION_LESS 9.0)
    INCLUDE(${VTK_USE_FILE})
    TARGET_LINK_LIBRARIES(liggghts_obj ${VTK_LIBRARIES})
    TARGET_LINK_LIBRARIES(liggghts_static ${VTK_LIBRARIES})
    TARGET_LINK_LIBRARIES(liggghts_shared ${VTK_LIBRARIES})
    TARGET_LINK_LIBRARIES(liggghts_bin ${VTK_LIBRARIES})
  else()

    # Only initialize minimum required VTK modules
    SET(VTK_MODULES
        "VTK::CommonDataModel" "VTK::IOLegacy" "VTK::IOParallelXML"
        "VTK::ParallelMPI" "VTK::FiltersCore")

    TARGET_LINK_LIBRARIES(liggghts_obj ${VTK_MODULES})
    TARGET_LINK_LIBRARIES(liggghts_static ${VTK_MODULES})
    TARGET_LINK_LIBRARIES(liggghts_shared ${VTK_MODULES})
    TARGET_LINK_LIBRARIES(liggghts_bin ${VTK_MODULES})

    vtk_module_autoinit(TARGETS liggghts_obj liggghts_static liggghts_shared liggghts_bin
                        MODULES ${VTK_MODULES})
  endif()

  MESSAGE(STATUS "Found VTK")
ELSE(VTK_FOUND)
  MESSAGE(STATUS "VTK NOT found!")
ENDIF(VTK_FOUND)

#=======================================


FIND_PACKAGE(JPEG)

IF(JPEG_FOUND)
  INCLUDE_DIRECTORIES(${JPEG_INCLUDE_DIR})
  ADD_DEFINITIONS(-DLAMMPS_JPEG)
  TARGET_LINK_LIBRARIES(liggghts_static ${JPEG_LIBRARIES})
  TARGET_LINK_LIBRARIES(liggghts_shared ${JPEG_LIBRARIES})
  TARGET_LINK_LIBRARIES(liggghts_bin ${JPEG_LIBRARIES})
ELSE(JPEG_FOUND)
  MESSAGE(STATUS "JPEG NOT found!")
ENDIF(JPEG_FOUND)

#=======================================

FIND_PACKAGE(MPI)

IF(MPI_FOUND)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
  TARGET_LINK_LIBRARIES(liggghts_static ${MPI_LIBRARIES})
  TARGET_LINK_LIBRARIES(liggghts_shared ${MPI_LIBRARIES})
  TARGET_LINK_LIBRARIES(liggghts_bin ${MPI_LIBRARIES})

  IF(MPI_COMPILE_FLAGS)
    SET_TARGET_PROPERTIES(liggghts_static PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(liggghts_shared PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(liggghts_bin PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  ENDIF()

  IF(MPI_LINK_FLAGS)
    SET_TARGET_PROPERTIES(liggghts_static PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    SET_TARGET_PROPERTIES(liggghts_shared PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    SET_TARGET_PROPERTIES(liggghts_bin PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  ENDIF()
ELSE(MPI_FOUND)
  MESSAGE(STATUS "MPI NOT found! Using MPI Stubs instead.")
  ADD_SUBDIRECTORY(STUBS)
  TARGET_LINK_LIBRARIES(liggghts_static mpi_stubs)
  TARGET_LINK_LIBRARIES(liggghts_shared mpi_stubs)
  TARGET_LINK_LIBRARIES(liggghts_bin mpi_stubs)
ENDIF(MPI_FOUND)

#======================================

#install(TARGETS liggghts liggghts_bin
#        RUNTIME DESTINATION bin
#        LIBRARY DESTINATION lib)

#install(TARGETS liggghts_bin
#        RUNTIME DESTINATION bin
#        LIBRARY DESTINATION lib)

install(TARGETS liggghts_bin liggghts_static liggghts_shared
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER DESTINATION include)
IF (WIN32)
    install(DIRECTORY WINDOWS DESTINATION include)
ENDIF()

