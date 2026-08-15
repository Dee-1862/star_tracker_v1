cmake_minimum_required(VERSION 3.16)

project(
    star_tracker_suite
    VERSION 0.1.0
    DESCRIPTION "Allocation-free satellite star tracker flight software"
    LANGUAGES CXX
)

option(STAR_TRACKER_BUILD_TESTS "Build flight software unit tests" ON)
option(STAR_TRACKER_BUILD_BENCHMARKS "Build host-side benchmark tools" ON)

if(STAR_TRACKER_BUILD_TESTS)
    include(CTest)
endif()

add_subdirectory(flight_software)
