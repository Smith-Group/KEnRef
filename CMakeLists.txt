cmake_minimum_required(VERSION 3.24)
project(KEnRef)
set(CMAKE_CXX_STANDARD 17)

string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)
message(STATUS "the build type is: ${build_type}")

#TODO rearrange setting *_DIR variables after reading https://cmake.org/cmake/help/latest/command/find_package.html [then delete this line :)]
##change these if needed
set(SUFFIX_MPI _mpi)
set(GROMACS_SRC_DIR /home/amr/CLionProjects/gromacs)
set(GROMACS_BUILD_DIR /home/amr/CLionProjects/gromacs/cmake-build-${build_type})

set(GROMACS_INSTALL_DIR /smithlab/opt/gromacs-dev1/${build_type})
set(GROMACS_DIR /smithlab/opt/gromacs-dev1/${build_type}/share/cmake/gromacs${SUFFIX_MPI}) #where libgromacs*.cmake are located
#set(GROMACS_INSTALL_DIR /usr/local/gromacs/${build_type})
#set(GROMACS_DIR /usr/local/gromacs/${build_type}/share/cmake/gromacs${SUFFIX_MPI}) #where libgromacs*.cmake are located

#set(CMAKE_MODULE_PATH /usr/local/gromacs-dev/share/gromacs/template/cmake:$CMAKE_MODULE_PATH)
set(Eigen3_DIR /smithlab/opt/eigen/3.4.0/share/eigen3/cmake) #:/smithlab/opt/eigen/3.4.0/include/eigen3/unsupported
## comment the next line to make it float. Uncomment to make it double.
#add_compile_options(-DDOUBLE)

find_package(MPI REQUIRED)
#find_package(Eigen3 REQUIRED)
find_package (Eigen3 3.4 REQUIRED)
if (Eigen3_FOUND)
	message("Eigen3 FOUND")
	message(STATUS "EIGEN3_INCLUDE_DIR is: ${EIGEN3_INCLUDE_DIR}")
	message(STATUS "EIGEN3_INCLUDE_DIRS is: ${EIGEN3_INCLUDE_DIRS}")
	message(STATUS "EIGEN3_LIBRARIES is: ${EIGEN3_LIBRARIES}")
endif ()


find_package(GROMACS NAMES gromacs_mpi gromacs)

# find and use OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
	message("OPENMP FOUND")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()


FILE(GLOB project_sources
		src/core/IoUtils.cpp
		src/core/KEnRef.cpp
		src/core/kenrefinitializer.cpp
		src/gmxinterface/gmxkenrefinitializer.cpp
		src/gmxinterface/gmxkenrefrunner.cpp
		src/gmxinterface/gmxwrapper.cpp
		src/gmxinterface/KEnRefForceProvider.cpp
		src/gmxinterface/KEnRefMDModule.cpp
		src/gmxinterface/KENRefRunner.cpp
		src/gmxinterface/mdrun.cpp
	)

FILE(GLOB project_includes
		config/KEnRefConfig.h
		core/ensembleutils.h
		core/IoUtils.h
		core/kabsch.h
		core/KEnRef.h
		core/kenreffinalizer.h
		core/kenrefinitializer.h
		core/kenrefrunner.h
		gmxinterface/gmxkenrefinitializer.h
		gmxinterface/gmxkenrefrunner.h
		gmxinterface/gmxwrapper.h
		gmxinterface/KEnRefForceProvider.h
		gmxinterface/KEnRefMDModule.h
	)

include_directories(include)
include_directories("${GROMACS_INSTALL_DIR}/include"
        "${GROMACS_BUILD_DIR}/src/include/"
        "${GROMACS_SRC_DIR}/src/gromacs/gmxlib"
        "${GROMACS_SRC_DIR}/src/gromacs/mdtypes"
        "${GROMACS_SRC_DIR}/src/gromacs/utility/include"
        "${GROMACS_SRC_DIR}/src/programs"
        "${GROMACS_SRC_DIR}/src"
        "${GROMACS_SRC_DIR}/src/include"
        "${GROMACS_SRC_DIR}/src/gromacs/math/include"
		"${EIGEN3_INCLUDE_DIR}"
        )

add_library(kenref ${project_sources})
add_executable(KEnRef ${project_sources})

target_link_libraries(KEnRef PRIVATE MPI::MPI_CXX Eigen3::Eigen
#        Eigen3::unsupported #TODO restore it later
	${GROMACS_INSTALL_DIR}/lib64/libgromacs_mpi.so)
target_link_libraries(kenref PRIVATE MPI::MPI_CXX Eigen3::Eigen
#        Eigen3::unsupported #TODO restore it later
	${GROMACS_INSTALL_DIR}/lib64/libgromacs_mpi.so)
###########################################
## Add gtest
add_subdirectory(google_tests)
###########################################

