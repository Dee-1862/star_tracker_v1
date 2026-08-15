file(GLOB_RECURSE STAR_TRACKER_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

if(STAR_TRACKER_SOURCES)
    add_library(star_tracker_fsw STATIC ${STAR_TRACKER_SOURCES})
    target_include_directories(star_tracker_fsw PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
    )
    target_compile_features(star_tracker_fsw PUBLIC cxx_std_11)
    set_target_properties(star_tracker_fsw PROPERTIES
        CXX_EXTENSIONS OFF
        POSITION_INDEPENDENT_CODE OFF
    )

    if(MSVC)
        target_compile_options(star_tracker_fsw PRIVATE /W4 /permissive-)
    else()
        target_compile_options(star_tracker_fsw PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wshadow
        )
    endif()
else()
    add_library(star_tracker_fsw INTERFACE)
    target_include_directories(star_tracker_fsw INTERFACE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
    )
    target_compile_features(star_tracker_fsw INTERFACE cxx_std_11)
endif()

if(BUILD_TESTING AND STAR_TRACKER_BUILD_TESTS)
    file(GLOB STAR_TRACKER_TEST_SOURCES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/test/*_test.cpp"
    )

    foreach(TEST_SOURCE IN LISTS STAR_TRACKER_TEST_SOURCES)
        get_filename_component(TEST_NAME "${TEST_SOURCE}" NAME_WE)
        add_executable("${TEST_NAME}" "${TEST_SOURCE}")
        target_link_libraries("${TEST_NAME}" PRIVATE star_tracker_fsw)
        add_test(NAME "${TEST_NAME}" COMMAND "${TEST_NAME}")
    endforeach()
endif()

if(STAR_TRACKER_BUILD_BENCHMARKS)
    add_executable(
        star_tracker_benchmark_runner
        "${CMAKE_CURRENT_SOURCE_DIR}/benchmark_runner.cpp"
    )
    target_link_libraries(star_tracker_benchmark_runner PRIVATE star_tracker_fsw)
    if(WIN32)
        target_link_libraries(star_tracker_benchmark_runner PRIVATE Psapi)
    endif()

    # Frozen-centroid adapter: lets ours be compared against tetra3 and LOST
    # with centroiding removed as a variable.
    add_executable(
        star_tracker_centroid_runner
        "${CMAKE_CURRENT_SOURCE_DIR}/centroid_runner.cpp"
    )
    target_link_libraries(star_tracker_centroid_runner PRIVATE star_tracker_fsw)
endif()
