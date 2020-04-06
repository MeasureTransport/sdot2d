cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(SemiDiscreteOT VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set c++ compiler flags
set(CMAKE_CXX_FLAGS "-O2 -std=c++11")

# Set the source files for the main SDOT library
set(LibrarySources
    src/BoundingBox.cpp
    src/LaguerreDiagram.cpp
    src/PolygonRasterize.cpp
    src/RegularGrid.cpp
    src/DiscretizedDistribution.cpp
    src/SemiDiscreteOT.cpp
)

# Set the source files for the python wrapper
set(WrapperSources
    src/PythonBindings.cpp
)


include_directories(include)

set(ExternalLibraries
    gmp;CGAL;mpfr)

add_library(sdot SHARED ${LibrarySources})
target_link_libraries(sdot ${ExternalLibraries})

add_executable(ConstructLaguerre tests/ConstructLaguerre.cpp)
target_link_libraries(ConstructLaguerre sdot ${OtherLibraries})

add_executable(ComputeSDOT tests/ComputeSDOT.cpp)
target_link_libraries(ComputeSDOT sdot ${OtherLibraries})


# Add pybind11 wrapper
set(PYBIND11_CPP_STANDARD -std=c++11)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/pybind11)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/pybind11/include)
pybind11_add_module(pysdot SHARED NO_EXTRAS ${WrapperSources})
target_link_libraries(pysdot sdot ${OtherLibraries})
