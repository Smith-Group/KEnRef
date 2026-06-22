cmake_minimum_required(VERSION 3.24)
project(Google_tests)

include(FetchContent)
message(STATUS "CMAKE_HOME_DIRECTORY is: ${CMAKE_HOME_DIRECTORY}")
# Pin googletest 1.17.0. FetchContent downloads it into THIS build dir's default _deps cache and builds
# it there, so the project configures cleanly in ANY build directory and is fully offline after the
# first fetch (per build dir). For an offline FIRST setup, point FETCHCONTENT_SOURCE_DIR_GOOGLETEST at a
# local extraction to skip the download. (Previously the FetchContent base was hardcoded to
# cmake-build-${build_type} with DOWNLOAD_DIR/SOURCE_SUBDIR + a manual add_subdirectory, which collided
# on the google_tests binary dir when configured in a separate build dir.)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/releases/download/v1.17.0/googletest-1.17.0.tar.gz
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

##TODO use the next line to find resources folder and test files (later). You may combine it with creating a header file
#cmake_path(ABSOLUTE_PATH "../res" BASE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" NORMALIZE OUTPUT_VARIABLE RESOURCES_ABSOLUTE_PATH)

FetchContent_MakeAvailable(googletest)
include_directories(${gtest_SOURCE_DIR})

message(STATUS "gtest_SOURCE_DIR is: ${gtest_SOURCE_DIR}")
include(GoogleTest)

# Now simply link against gtest or gtest_main as needed. Eg
add_executable(Google_tests_kenref_core_exe testKabsch.cpp testKEnRef.cpp testIoUtils.cpp testRestoreNoJump.cpp)
# Use namespace target for compatibility with both build and installed versions
if(TARGET KENREF::CORE)
	target_link_libraries(Google_tests_kenref_core_exe gtest gtest_main KENREF::CORE)
else()
	target_link_libraries(Google_tests_kenref_core_exe gtest gtest_main kenref_core)
endif()
target_include_directories(Google_tests_kenref_core_exe PRIVATE
	${gtest_SOURCE_DIR}/include
	${gtest_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/../include_core
)
gtest_discover_tests(Google_tests_kenref_core_exe)

if(BUILD_KENREF_GMX OR TARGET kenref_gmxinterface)
	add_executable(Google_tests_kenref_gmxinterface_exe testKEnRefForceProvider.cpp)
	if(TARGET KENREF::CORE)
		target_link_libraries(Google_tests_kenref_gmxinterface_exe gtest gtest_main KENREF::CORE kenref_gmxinterface)
	else()
		target_link_libraries(Google_tests_kenref_gmxinterface_exe gtest gtest_main kenref_core kenref_gmxinterface)
	endif()
	target_include_directories(Google_tests_kenref_gmxinterface_exe PRIVATE
		${gtest_SOURCE_DIR}/include
		${gtest_SOURCE_DIR}
		${CMAKE_CURRENT_SOURCE_DIR}/../include_core
	)
	gtest_discover_tests(Google_tests_kenref_gmxinterface_exe)
endif()
