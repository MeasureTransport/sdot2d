cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(SemiDiscreteOT VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set c++ compiler flags
set(CMAKE_CXX_FLAGS "-O2")

# Set the source files
set(LibrarySources
    src/LaguerreDiagram.cpp
    src/PolygonRasterize.cpp
    src/RegularGrid.cpp
    src/DiscretizedDistribution.cpp
    src/SemiDiscreteOT.cpp
    )

include_directories(include)

set(ExternalLibraries
    gmp;CGAL;mpfr)

add_library(sdot ${LibrarySources})
target_link_libraries(sdot ${ExternalLibraries})


add_executable(ConstructLaguerre tests/ConstructLaguerre.cpp)
target_link_libraries(ConstructLaguerre sdot ${OtherLibraries})

add_executable(ComputeSDOT tests/ComputeSDOT.cpp)
target_link_libraries(ComputeSDOT sdot ${OtherLibraries})
