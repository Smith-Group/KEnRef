cmake_minimum_required(VERSION 3.24)
project(KEnRef)

string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)
message(STATUS "the build type is: ${build_type}")

#change these if needed
set(GROMACS_SRC_DIR /home/amr/CLionProjects/gromacs)
set(GROMACS_BUILD_DIR /home/amr/CLionProjects/gromacs/cmake-build-${build_type})
set(GROMACS_INSTALL_DIR /smithlab/opt/gromacs-dev/${build_type})
## comment the next line to make it float. Uncomment to make it double.
add_compile_options(-DDOUBLE)

set(CMAKE_CXX_STANDARD 17)
set(GROMACS_DIR /smithlab/opt/gromacs-dev/${build_type}) #set by me then changed/masked (?) by CMake
#set(CMAKE_MODULE_PATH /usr/local/gromacs-dev/share/gromacs/template/cmake:$CMAKE_MODULE_PATH)
set(EIGEN3_INCLUDE_DIR /smithlab/opt/eigen/3.4.0/include/eigen3) #:/smithlab/opt/eigen/3.4.0/include/eigen3/unsupported

find_package(MPI REQUIRED)
find_package(Eigen3 REQUIRED)
#find_package (Eigen3 3.4 REQUIRED)

find_package(GROMACS NAMES gromacs_mpi)
#find_package(GROMACS)

include_directories(config core gmxinterface)
include_directories("${GROMACS_INSTALL_DIR}/include"
        "${GROMACS_BUILD_DIR}/src/include/"
        "${GROMACS_SRC_DIR}/src/gromacs/gmxlib"
        "${GROMACS_SRC_DIR}/src/gromacs/mdtypes"
        "${GROMACS_SRC_DIR}/src/gromacs/utility/include"
        "${GROMACS_SRC_DIR}/src/programs"
        "${GROMACS_SRC_DIR}/src"
        "${GROMACS_SRC_DIR}/src/include"
        "${GROMACS_SRC_DIR}/src/gromacs/math/include"
        )



add_executable(KEnRef
        config/KEnRefConfig.h
        core/ensembleutils.h
        core/IoUtils.cpp
        core/IoUtils.h
        core/kabsch.h
        core/KEnRef.cpp
        core/KEnRef.h
        core/kenreffinalizer.h
        core/kenrefinitializer.cpp
        core/kenrefinitializer.h
        core/kenrefrunner.h
        gmxinterface/gmxkenrefinitializer.cpp
        gmxinterface/gmxkenrefinitializer.h
        gmxinterface/gmxkenrefrunner.cpp
        gmxinterface/gmxkenrefrunner.h
        gmxinterface/gmxwrapper.cpp
        gmxinterface/gmxwrapper.h
        gmxinterface/KEnRefForceProvider.cpp
        gmxinterface/KEnRefForceProvider.h
        gmxinterface/KEnRefMDModule.cpp
        gmxinterface/KEnRefMDModule.h
        gmxinterface/KENRefRunner.cpp
        gmxinterface/mdrun.cpp)

target_link_libraries(KEnRef PRIVATE MPI::MPI_CXX
        Eigen3::Eigen
#        Eigen3::unsupported #TODO restore it later
	${GROMACS_INSTALL_DIR}/lib64/libgromacs_mpi.so)

###########################################
## Add gtest
add_subdirectory(google_tests)
###########################################

