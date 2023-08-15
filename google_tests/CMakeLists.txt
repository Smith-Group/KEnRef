
project(Google_tests)

include(FetchContent)
message(STATUS "CMAKE_HOME_DIRECTORY is: ${CMAKE_HOME_DIRECTORY}")
set(FETCHCONTENT_BASE_DIR ${CMAKE_HOME_DIRECTORY}/google_tests/lib)
FetchContent_Declare(
        google_tests
        # Specify the commit you depend on and update it regularly.
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
        DOWNLOAD_DIR lib
        SOURCE_SUBDIR lib
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(google_tests_force_shared_crt ON CACHE BOOL "" FORCE)  #what should the first syllable be?

message(STATUS "google_tests_SOURCE_DIR is: ${google_tests_SOURCE_DIR}")
#set(google_tests_SOURCE_DIR Google_tests/lib)
message(STATUS "google_tests_SOURCE_DIR is: ${google_tests_SOURCE_DIR}")
message(STATUS "google_tests_SOURCE_DIR/include is: ${google_tests_SOURCE_DIR}/include")
message(STATUS "googletest_SOURCE_DIR is: ${google_tests_SOURCE_DIR}")

FetchContent_MakeAvailable(google_tests)

#add_subdirectory(${google_tests_SOURCE_DIR})
add_subdirectory(lib/google_tests-src)

message(STATUS "gtest_SOURCE_DIR is: ${gtest_SOURCE_DIR}")
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Now simply link against gtest or gtest_main as needed. Eg
add_executable(Google_tests_exe testKabsch.cpp)
target_link_libraries(Google_tests_exe gtest gtest_main)
##add_test(NAME Google_tests COMMAND example)