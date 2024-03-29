cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(SDOT)
set(${CMAKE_PROJECT_NAME}_MAJOR_VERSION 0)
set(${CMAKE_PROJECT_NAME}_MINOR_VERSION 1)
set(${CMAKE_PROJECT_NAME}_PATCH_VERSION 0)
set(${CMAKE_PROJECT_NAME}_VERSION ${${CMAKE_PROJECT_NAME}_MAJOR_VERSION}.${${CMAKE_PROJECT_NAME}_MINOR_VERSION}.${${CMAKE_PROJECT_NAME}_PATCH_VERSION})



# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set c++ compiler flags
set(CMAKE_CXX_FLAGS "-g")
set(CMAKE_BUILD_TYPE Release)

# Set the source files for the main SDOT library
set(LibrarySources
    src/BoundingBox.cpp
    src/LaguerreDiagram.cpp
    src/PolygonRasterize.cpp
    src/RegularGrid.cpp
    src/DiscretizedDistribution.cpp
    src/SemiDiscreteOT.cpp
    src/Polygon.cpp
    src/OptionUtilities.cpp
    src/Distances/Wasserstein2.cpp
    src/Distances/QuadraticRegularizationFunctions.cpp
    src/Distances/GHKFunctions.cpp
    src/Distances/QuadratureDistance.cpp
    src/Distances/TriangularQuadrature.cpp
    src/Distances/RectangularQuadrature.cpp
    src/Distances/LineQuadrature.cpp
)

# Set the source files for the python wrapper
set(WrapperSources
    src/PythonBindings.cpp
)

##################################
# Find CGAL
find_package(CGAL REQUIRED)
find_package(GMP)
find_package(MPFR)

include_directories(include ${CGAL_INCLUDE_DIRS})
set(ExternalLibraries ${CGAL_LIBRARY} ${GMP_LIBRARIES} ${MPFR_LIBRARIES})

# If we weren't able to find GMP or MPFR, just add the library names and hope for the best
if(NOT GMP_FOUND)
  set(ExternalLibraries ${ExternalLibraries} -lgmp)
endif()
if(NOT MPFR_FOUND)
  set(ExternalLibraries ${ExternalLibraries} -lmpfr)
endif()

add_library(sdot SHARED ${LibrarySources})
target_link_libraries(sdot ${ExternalLibraries})



###################################
## Add test executables
enable_testing()
find_package(GTest HINTS ${SDOT_GTEST_DIR})

if(GTest_FOUND)

  message(STATUS "GTest found.")

  set(unitTestSources
    tests/unit/DistanceTests.cpp
    tests/unit/LaguerreDiagramTests.cpp
    tests/unit/SemiDiscreteOTTests.cpp
    tests/unit/QuadratureTests.cpp
    tests/unit/RunUnitTests.cpp
  )

  add_executable(RunUnitTests ${unitTestSources})
  target_link_libraries(RunUnitTests sdot ${OtherLibraries} GTest::gtest)
  add_test(UnitTests RunUnitTests)

else()
  message(STATUS "GTest was not found and unit tests will be skipped.")
endif()

add_executable(ConstructLaguerre tests/ConstructLaguerre.cpp)
target_link_libraries(ConstructLaguerre sdot ${OtherLibraries})
add_test(BasicLaguerre ConstructLaguerre)

add_executable(ComputeSDOT tests/ComputeSDOT.cpp)
target_link_libraries(ComputeSDOT sdot ${OtherLibraries})
add_test(BasicOT ComputeSDOT)

add_executable(Convexify tests/Convexify.cpp)
target_link_libraries(Convexify sdot ${OtherLibraries})
add_test(Convexify Convexify)

add_executable(CircleProblem tests/CircleProblem.cpp)
target_link_libraries(CircleProblem sdot ${OtherLibraries})
add_test(SimpleCircles CircleProblem)

###################################
## Python Bindings
# Add pybind11 wrapper
set(PYBIND11_CPP_STANDARD -std=c++11)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/pybind11)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/pybind11/include)
pybind11_add_module(_pysdot SHARED NO_EXTRAS ${WrapperSources})
target_link_libraries(_pysdot PRIVATE sdot ${OtherLibraries})


# Add installation directives
install(TARGETS sdot _pysdot LIBRARY DESTINATION lib)

install(DIRECTORY src/python
        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        FILES_MATCHING PATTERN "*.py")

install(DIRECTORY include/SDOT
        DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
	      FILES_MATCHING PATTERN "*.h")



###############################################
# Add a target to build documentation
###############################################

find_package(Doxygen)
if(DOXYGEN_FOUND)
	message(STATUS "To build documentation, run 'make doc'")

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/sdot.doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/sdot.doxyfile @ONLY)
  add_custom_target(doc
      COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/sdot.doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM
  )

endif(DOXYGEN_FOUND)
