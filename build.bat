@echo off
rem ============================================================================
rem build.bat - KEnRef build orchestrator (Windows)
rem
rem Windows builds kenref_core ONLY. The MD engines KEnRef plugs into -- PLUMED
rem and GROMACS -- target Linux/macOS (autotools PLUMED, `plumed patch`, the
rem GROMACS force provider); they are not built here. For the full
rem kenref_core -> PLUMED -> GROMACS chain use build.sh on Linux/macOS/WSL.
rem
rem Toolchain specifics come from the environment (CXX/CC/CXXFLAGS) or your
rem CMake generator; this script keeps only portable defaults.
rem ============================================================================
setlocal enabledelayedexpansion

set "REPO_ROOT=%~dp0"
set "BUILD_TYPE=Release"
set "ACCEL=AVX2_256"
set "PREFIX=%REPO_ROOT%install\kenref"

:parse
if "%~1"=="" goto after
if /I "%~1"=="--build-type" ( set "BUILD_TYPE=%~2" & shift & shift & goto parse )
if /I "%~1"=="--accel"      ( set "ACCEL=%~2"      & shift & shift & goto parse )
if /I "%~1"=="--prefix"     ( set "PREFIX=%~2"     & shift & shift & goto parse )
if /I "%~1"=="--with-plumed"  goto engines_unsupported
if /I "%~1"=="--with-gromacs" goto engines_unsupported
if /I "%~1"=="-h"     goto usage
if /I "%~1"=="--help" goto usage
echo ERROR: unknown option "%~1" (try --help) 1>&2
exit /b 2

:engines_unsupported
echo ERROR: PLUMED/GROMACS integration is Linux/macOS-only. 1>&2
echo        Use build.sh (Linux/macOS/WSL) for the full chain; build.bat builds kenref_core only. 1>&2
exit /b 2

:usage
echo Usage: build.bat [--build-type Release^|Debug^|RelWithDebInfo] [--accel AVX_512^|AVX_256^|AVX2_256] [--prefix DIR]
echo(
echo Builds kenref_core only. PLUMED/GROMACS are Linux/macOS-only (see build.sh).
exit /b 0

:after
set "BUILD_DIR=%REPO_ROOT%cmake-build-orch-win"
echo ==^> kenref_core (%BUILD_TYPE% / %ACCEL%) -^> %PREFIX%

cmake -S "%REPO_ROOT%." -B "%BUILD_DIR%" ^
    -DCMAKE_BUILD_TYPE=%BUILD_TYPE% ^
    -DACCEL=%ACCEL% ^
    -DBUILD_KENREF_CORE=ON ^
    -DBUILD_KENREF_GMX=OFF ^
    -DBUILD_KENREF_PLUMED=OFF ^
    -DCMAKE_INSTALL_PREFIX="%PREFIX%"
if errorlevel 1 exit /b 1

cmake --build "%BUILD_DIR%" --config %BUILD_TYPE%
if errorlevel 1 exit /b 1

cmake --install "%BUILD_DIR%" --config %BUILD_TYPE%
if errorlevel 1 exit /b 1

echo ==^> DONE. kenref_core installed to %PREFIX%
endlocal
