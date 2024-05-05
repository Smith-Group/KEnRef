cmake_minimum_required(VERSION 3.24)
project(KEnRef)
set(CMAKE_CXX_STANDARD 17)

string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)
message(STATUS "the build type is: ${build_type}")

if (NOT DEFINED ACCEL)
	set(ACCEL "AVX2_256")
	message(STATUS "ACCEL variable not defined. Setting it to ${ACCEL}")
else ()
	message(STATUS "ACCEL variable is: ${ACCEL}")
endif ()
#TODO rearrange setting *_DIR variables after reading https://cmake.org/cmake/help/latest/command/find_package.html [then delete this line :)]
##change these if needed
set(SUFFIX_MPI _mpi)
set(GROMACS_SRC_DIR /home/amr/git/gromacs)
set(GROMACS_BUILD_DIR /home/amr/git/gromacs/cmake-build-${build_type}/mpi)
set(GROMACS_INSTALL_DIR /smithlab/opt/gromacs-dev/2023/2023.4/${build_type}/${ACCEL})
set(GROMACS_DIR ${GROMACS_INSTALL_DIR}/share/cmake/gromacs${SUFFIX_MPI}) #where libgromacs*.cmake are located

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	string(FIND ${GROMACS_INSTALL_DIR} "gromacs-dev" replaceable)
	if (replaceable GREATER_EQUAL -1)
		string(REPLACE "gromacs-dev" "kenref-dev" proposed_path ${GROMACS_INSTALL_DIR})
		set(CMAKE_INSTALL_PREFIX ${proposed_path} CACHE PATH "root directory to install KEnRef (default is a tree parallel to GROMACS_INSTALL_DIR)" FORCE)
	endif ()
endif()

#set(CMAKE_MODULE_PATH /usr/local/gromacs-dev/share/gromacs/template/cmake:$CMAKE_MODULE_PATH)
set(Eigen3_DIR /smithlab/opt/eigen/3.4.0/share/eigen3/cmake) #:/smithlab/opt/eigen/3.4.0/include/eigen3/unsupported
## comment the next line to make it float. Uncomment to make it double.
#add_compile_options(-DDOUBLE)

find_package(MPI REQUIRED)

find_package (Eigen3 3.4 REQUIRED)
if (Eigen3_FOUND)
	message("Eigen3 FOUND")
	message(STATUS "EIGEN3_INCLUDE_DIR is: ${EIGEN3_INCLUDE_DIR}")
	message(STATUS "EIGEN3_INCLUDE_DIRS is: ${EIGEN3_INCLUDE_DIRS}")
endif ()

message(STATUS "Searching for Gromacs in ${GROMACS_DIR}")
find_package(GROMACS NAMES gromacs_mpi gromacs REQUIRED)
if (GROMACS_FOUND)
	message("GROMACS FOUND")
	message(STATUS "GROMACS_DIR folder is: ${GROMACS_DIR}")
	message(STATUS "GROMACS_CONFIG file is: ${GROMACS_CONFIG}")
endif ()

find_package(OpenMP REQUIRED)
#if (OPENMP_FOUND)
#	message("OPENMP FOUND")
#	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif ()

FILE(GLOB core_sources "src/core/*.cpp")
FILE(GLOB gmxinterface_sources "src/gmxinterface/*.cpp")

include_directories(include)
include_directories(
#		"${GROMACS_INSTALL_DIR}/include"  # Not neeed here anymore, but needed for libgromacs
        "${GROMACS_BUILD_DIR}/src/include/"  #  Needed for "config.h"
#        "${GROMACS_SRC_DIR}/src/gromacs/gmxlib"
        "${GROMACS_SRC_DIR}/src/gromacs/mdtypes"
        "${GROMACS_SRC_DIR}/src/gromacs/utility/include"
        "${GROMACS_SRC_DIR}/src/programs"
        "${GROMACS_SRC_DIR}/src"
        "${GROMACS_SRC_DIR}/src/include"
        "${GROMACS_SRC_DIR}/src/gromacs/math/include"
)
include_directories("${EIGEN3_INCLUDE_DIR}")

#TODO Change the targets into FENREF::CORE and KENREF::GMXINTERFACE (and rearrange sources)
add_library(KENREF_CORE STATIC ${core_sources})
add_library(KENREF_GMXINTERFACE STATIC ${gmxinterface_sources})
add_executable(KEnRef "src/gmxinterface/gmxwrapper.cpp")

target_link_libraries(KENREF_CORE PUBLIC MPI::MPI_CXX Eigen3::Eigen
		#        Eigen3::unsupported # may need to restore it later
)
target_link_libraries(KENREF_GMXINTERFACE PUBLIC KENREF_CORE Gromacs::libgromacs)
target_link_libraries(KEnRef PRIVATE KENREF_CORE KENREF_GMXINTERFACE)

install(TARGETS KEnRef DESTINATION bin)
install(TARGETS KENREF_CORE KENREF_GMXINTERFACE DESTINATION lib)
###########################################
## Add gtest
add_subdirectory(google_tests)
###########################################

