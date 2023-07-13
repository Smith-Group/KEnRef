cmake_minimum_required(VERSION 3.25)
project(KEnRef)

set(CMAKE_CXX_STANDARD 17)
#set(GROMACS_DIR /usr/local/gromacs-dev/share/cmake/gromacs_mpi)
set(GROMACS_DIR /usr/local/gromacs-dev)
#set(CMAKE_MODULE_PATH /usr/local/gromacs-dev/share/gromacs/template/cmake:$CMAKE_MODULE_PATH)
set(EIGEN3_INCLUDE_DIR /smithlab/opt/eigen/3.4.0/include/eigen3) #:/smithlab/opt/eigen/3.4.0/include/eigen3/unsupported

find_package(MPI REQUIRED)
find_package(Eigen3 REQUIRED)
#find_package (Eigen3 3.4 REQUIRED)

find_package(GROMACS NAMES gromacs_mpi)
#find_package(GROMACS)


include_directories(config)
include_directories(core)
include_directories(gmxinterface)
include_directories(/usr/local/gromacs-dev/include)
include_directories(/home/amr/git/gromacs-build/eclipse/src/include/)
include_directories(/home/amr/git/gromacs/src/gromacs/gmxlib)
include_directories(/home/amr/git/gromacs/src/gromacs/mdtypes)
include_directories(/home/amr/git/gromacs/src/gromacs/utility/include)
include_directories(/home/amr/git/gromacs/src/programs)
include_directories(/home/amr/git/gromacs/src)
include_directories(/home/amr/git/gromacs/src/include)
include_directories(/home/amr/git/gromacs/src/gromacs/math/include)



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
        /usr/local/gromacs-dev/lib64/libgromacs_mpi.so)
