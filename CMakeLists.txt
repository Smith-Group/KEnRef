cmake_minimum_required(VERSION 3.24)
project(KEnRef VERSION 2.0.0)
# Single source of truth for the KEnRef version. Flows to the pkg-config files (@PROJECT_VERSION@ in
# cmake/*.pc.in) and, via those, to PLUMED's configure floor (KENREF_MIN_VERSION in plumed configure.ac,
# checked as `kenref_core >= X`). Bump here when the public API / numerical contract changes; keep the
# PLUMED floor in lockstep. v2.0.0 == the model-abstraction restructure (branch v2.0 and its subbranches).
set(CMAKE_CXX_STANDARD 17)

# Emit compile_commands.json so clangd / IDEs resolve includes for every translation unit (including
# newly added core/model files). FORCE so it wins over a stale cache value. Harmless for the build.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Export compile_commands.json for clangd/IDEs" FORCE)

# Default to Release when the caller set nothing. build.sh always passes -DCMAKE_BUILD_TYPE, so this is for
# the raw-cmake path in INSTALL.md: an unset CMAKE_BUILD_TYPE means CMake adds NO optimisation flags at all,
# which silently produced an unoptimised install that still called itself one. GROMACS does exactly this
# (gromacs/CMakeLists.txt: `if(NOT CMAKE_BUILD_TYPE) set(... "Release" ...)`), so a KEnRef built beside it
# behaves the same way. Only applies with no cache value — an explicit -DCMAKE_BUILD_TYPE= always wins.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (Release/Debug/RelWithDebInfo/MinSizeRel)" FORCE)
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug RelWithDebInfo MinSizeRel)
	message(STATUS "CMAKE_BUILD_TYPE not set -> defaulting to Release (override with -DCMAKE_BUILD_TYPE=...)")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" build_type)   # quoted: tolerate an unset CMAKE_BUILD_TYPE
message(STATUS "the build type is: ${build_type}")

# ============================================================================
# BUILD OPTIONS
# ============================================================================
# Component switches are TRI-STATE STRINGs: ON / OFF / AUTO (default AUTO). See cmake/Components.cmake for
# the exact semantics — briefly: AUTO builds a component only if its heavy external engine (GROMACS / PLUMED)
# is ALREADY available, and NEVER downloads it; use =ON to force a build (fetching sub-deps per the download
# switch). The user's raw ON/OFF/AUTO stays in the cache; a resolved boolean of the SAME name shadows it for
# the rest of configure (so every `if(BUILD_KENREF_*)` below keeps working).
include(cmake/Components.cmake)
set(BUILD_KENREF_CORE   "AUTO" CACHE STRING "Build kenref_core (ON/OFF/AUTO)")
set(BUILD_KENREF_GMX    "AUTO" CACHE STRING "Build the kenref-gmx interface (ON/OFF/AUTO)")
# The PLUMED interface is OPT-IN — a plain kenref_core build neither COMPILES nor EXPORTS it (in production the
# PLUMED fork compiles src/plumedinterface itself, Option A).
#   * BUILD_KENREF_PLUMED          — COMPILE-VALIDATE plumedinterface against PLUMED's headers (needs
#                                    -DPLUMED_SRC_DIR). AUTO builds it iff those headers are present.
#   * KENREF_EXPORT_PLUMEDINTERFACE— INSTALL the plumedinterface headers + src/plumedinterface/*.cpp +
#                                    kenref_plumed.pc (the source a PLUMED build compiles in-tree, Option A).
#                                    Default OFF; forced ON below when a PLUMED build here needs it.
set(BUILD_KENREF_PLUMED "AUTO" CACHE STRING "Compile-validate the KENREF PLUMED interface (ON/OFF/AUTO; needs -DPLUMED_SRC_DIR)")
set_property(CACHE BUILD_KENREF_CORE BUILD_KENREF_GMX BUILD_KENREF_PLUMED PROPERTY STRINGS ON OFF AUTO)
option(KENREF_EXPORT_PLUMEDINTERFACE "Install the plumedinterface headers/source + kenref_plumed.pc (for a PLUMED build)" OFF)

# Download policy (single switch, also declared in cmake/Dependencies.cmake — declared here too so it is
# available for the component-resolution messages below; option() keeps the first/cache value): when ON
# (default), a COMMITTED component (ON, or AUTO that resolved ON) may FETCH a missing sub-dependency (Eigen;
# a stock GROMACS for the gmx path; the PLUMED fork for --with-plumed). When OFF, a missing sub-dependency is
# a clean FATAL instead of a download. It NEVER makes AUTO build an absent-engine component. build.sh exposes
# it as `--download ON|OFF`.
option(KENREF_FETCH_MISSING "Download & build missing dependencies (Eigen/GROMACS/PLUMED) instead of failing" ON)

# PLUMED: KEnRef does not REIMPLEMENT plumed's build; it DELEGATES to the plumed side (build-only.sh /
# build-and-batch.sh), and CMake itself invokes that script (via install(CODE) — after kenref_core is
# installed, so the plumed script reuses it). See cmake/InvokePlumed.cmake. build.sh just sets these.
# Tri-state (ON/OFF/AUTO): AUTO builds iff a PLUMED source is provided (KENREF_PLUMED_SRC_DIR). Resolved below.
set(KENREF_WITH_PLUMED "AUTO" CACHE STRING "Also build PLUMED (kenref module), delegated to the plumed side at install time (ON/OFF/AUTO)")
set_property(CACHE KENREF_WITH_PLUMED PROPERTY STRINGS ON OFF AUTO)
option(KENREF_WITH_PLUMED_GROMACS  "Also build PLUMED + a batched GROMACS (delegates to build-and-batch.sh)"        OFF)

# Per-model compile-time enable switches (Style-2 registration). Default ON; turning one OFF drops that
# model's registration from bootstrapModels() so it is unavailable at runtime. Passed to the core
# targets as KENREF_ENABLE_<MODEL> compile definitions (see cmake/KEnRefCore.cmake).
option(KENREF_ENABLE_SIGMA    "Enable the SIGMA energy model"    ON)
option(KENREF_ENABLE_PLATEAUS "Enable the PLATEAUS energy model" ON)
option(KENREF_ENABLE_RELAX    "Enable the RELAX energy model"    ON)

# Normalize the raw tri-state requests (upper-case) for the OFF-guards below; the final resolution to a
# boolean of the SAME variable name happens in COMPONENT RESOLUTION once engine availability is probed.
kenref_normalize_tristate("${BUILD_KENREF_CORE}"   _KN_CORE_REQ)
kenref_normalize_tristate("${BUILD_KENREF_GMX}"    _KN_GMX_REQ)
kenref_normalize_tristate("${BUILD_KENREF_PLUMED}" _KN_PLUMED_REQ)
kenref_normalize_tristate("${KENREF_WITH_PLUMED}"  _KN_WITHPLUMED_REQ)

### comment the next line to make it float. Uncomment to make it double.
add_compile_options(-DKENREF_DOUBLE)

if (NOT DEFINED ACCEL)
	# Auto-detect this machine's SIMD (like GROMACS): pick the highest of AVX_512 / AVX2_256 / AVX_256.
	# Override with -DACCEL=... for a fixed/portable build. Linux reads /proc/cpuinfo; elsewhere -> AVX2_256.
	set(ACCEL "AVX2_256")   # portable fallback
	if(EXISTS "/proc/cpuinfo")
		file(READ "/proc/cpuinfo" _kn_cpuinfo LIMIT 16384)
		if(_kn_cpuinfo MATCHES "avx512f")
			set(ACCEL "AVX_512")
		elseif(_kn_cpuinfo MATCHES "avx2")
			set(ACCEL "AVX2_256")
		elseif(_kn_cpuinfo MATCHES "[ \t]avx[ \t]")
			set(ACCEL "AVX_256")
		endif()
		unset(_kn_cpuinfo)
	endif()
	message(STATUS "ACCEL not set -> auto-detected ${ACCEL} for this machine (override with -DACCEL=...)")
else ()
	message(STATUS "ACCEL variable is: ${ACCEL}")
endif ()

# ACCEL -> compiler flags, matched to GROMACS's GMX_SIMD tier NAMES (authoritative: gromacs
# cmake/gmxSimdFlags.cmake). What the SIMD-match guard enforces is the Eigen ALIGN CLASS (16/32/64) the
# compiler width implies — NOT the tier name; the tier name is recorded for diagnostics. NB: a GROMACS tier
# name encodes GROMACS's internal kernel *datapath* width, so AVX_128_FMA / AVX2_128 (despite "128") still
# emit 256-bit AVX at the compiler level -> Eigen align 32, same class as the *_256 tiers.
# The three common x86 tiers keep their proven -march (the "incompatibility bug with PLUMED" fix).
# The flags are collected into KENREF_SIMD_FLAGS rather than passed straight to add_compile_options(),
# because they are part of kenref_core's ABI and must reach CONSUMERS too: they are re-emitted into the
# installed pkg-config files (see cmake/*.pc.in). A consumer that compiles Eigen with a different width
# gets a different EIGEN_MAX_ALIGN_BYTES, which the guard in core/EigenAbiCheck.h then rejects at compile
# time. Shipping the flag in the .pc means a plain `pkg-config --cflags` consumer is correct by default
# instead of having to learn this build's ACCEL out of band.
if(ACCEL STREQUAL "AVX_512")
	set(KENREF_SIMD_FLAGS -march=skylake-avx512)
elseif(ACCEL STREQUAL "AVX2_256")
	set(KENREF_SIMD_FLAGS -march=haswell)
elseif(ACCEL STREQUAL "AVX_256")
	set(KENREF_SIMD_FLAGS -march=sandybridge)
elseif(ACCEL STREQUAL "AVX2_128")
	set(KENREF_SIMD_FLAGS -mavx2 -mfma)          # AVX2, AMD Zen 128-bit datapath; compiler 256-bit -> align 32
elseif(ACCEL STREQUAL "AVX_128_FMA")
	set(KENREF_SIMD_FLAGS -mavx -mfma)           # AMD fam15h (true FMA4 parts need -mfma4 via CXXFLAGS) -> align 32
elseif(ACCEL STREQUAL "AVX_512_KNL")
	message(WARNING "ACCEL=AVX_512_KNL targets the DEPRECATED Knights Landing; prefer AVX_512 on modern x86.")
	set(KENREF_SIMD_FLAGS -mavx512f -mavx512er -mfma)
elseif(ACCEL STREQUAL "SSE4.1")
	set(KENREF_SIMD_FLAGS -msse4.1)
elseif(ACCEL STREQUAL "SSE2")
	set(KENREF_SIMD_FLAGS -msse2)
elseif(ACCEL STREQUAL "ARM_NEON_ASIMD")
	set(KENREF_SIMD_FLAGS -march=armv8-a+simd)
elseif(ACCEL STREQUAL "ARM_SVE")
	set(KENREF_SIMD_FLAGS -march=armv8-a+sve)    # SVE vector length is CPU / -msve-vector-bits dependent
elseif(ACCEL STREQUAL "IBM_VSX")
	set(KENREF_SIMD_FLAGS -mcpu=power8 -mvsx)
elseif(ACCEL STREQUAL "REFERENCE" OR ACCEL STREQUAL "NONE")
	set(KENREF_SIMD_FLAGS "")   # scalar reference build (Eigen may still use baseline SSE2 on x86-64)
else()
	message(FATAL_ERROR "Unknown ACCEL='${ACCEL}'. Valid (GROMACS GMX_SIMD tiers): AVX_512 AVX_512_KNL "
		"AVX2_256 AVX2_128 AVX_256 AVX_128_FMA SSE4.1 SSE2 ARM_NEON_ASIMD ARM_SVE IBM_VSX REFERENCE")
endif()
if(KENREF_SIMD_FLAGS)
	add_compile_options(${KENREF_SIMD_FLAGS})
endif()
# Space-separated form for the pkg-config templates (@KENREF_SIMD_FLAGS_STR@).
string(REPLACE ";" " " KENREF_SIMD_FLAGS_STR "${KENREF_SIMD_FLAGS}")

# Bake THIS build's ACCEL tier name into every TU compiled here (core, gmx, tests, the simd-check exes) as the
# diagnostic ISA string for the SIMD guard (see core/SimdSignature.h). Directory-scope (NOT a target INTERFACE
# definition), so a consumer built separately records its OWN tier rather than inheriting kenref_core's.
add_compile_definitions(KENREF_ACCEL="${ACCEL}")
# ============================================================================
# GROMACS (only for the kenref-gmx force provider)
# ============================================================================
# kenref-gmx is a GROMACS *plugin* (a ForceProvider / MDModule), so it needs THREE things from GROMACS, and
# each exists in exactly ONE place:
#   * the BUILD tree   -> the GENERATED src/include/config.h (GMX_LIB_MPI, GMX_DOUBLE, GMX_SIMD_*).
#                         Not present in the source, and never installed.
#   * the SOURCE tree  -> the INTERNAL headers GROMACS does not install (mdtypes/iforceprovider.h,
#                         mdtypes/imdmodule.h, mdrun/simulationcontext.h, mdrunutility/mdmodulesnotifiers.h,
#                         options/optionsection.h, utility/keyvaluetreebuilder.h, ...): the MDModules
#                         extension API is internal, not part of the public installed API.
#   * the INSTALL tree -> libgromacs*.so + its CMake package (the actual link target).
# So an INSTALLED GROMACS ALONE can never satisfy kenref-gmx.
#
# You normally pass just ONE parameter:
#
#   -DGROMACS_BUILD_DIR=<GROMACS's own CMake BUILD directory>
#        This is the directory you ran `cmake -B` into when you built GROMACS: it contains CMakeCache.txt and
#        src/include/config.h.  It is NOT the install prefix (the one with bin/gmx and share/cmake/gromacs*),
#        and it is NOT the source checkout.  From its CMakeCache.txt we read the matching source
#        (CMAKE_HOME_DIRECTORY) and install (CMAKE_INSTALL_PREFIX), and from its config.h the MPI flavor — so
#        the three trees are consistent BY CONSTRUCTION and cannot silently drift apart.
#
# Or, to have KEnRef build GROMACS for you from YOUR checkout:
#
#   -DGROMACS_SRC_DIR=<GROMACS source checkout>      (the directory containing GROMACS's top-level CMakeLists)
#
# Pass neither and KEnRef fetches + builds a stock GROMACS (when downloads are enabled).
# NOTHING is guessed: there is deliberately no sibling-directory ("../gromacs") or "/usr/local/gromacs"
# fallback — kenref and gromacs are not guaranteed to be co-located, and a wrong guess could silently pick a
# DIFFERENT GROMACS version. GROMACS_INSTALL_DIR is not a parameter: it is always derived.
if(NOT _KN_GMX_REQ STREQUAL "OFF")
	# env fallbacks for the two real parameters (an explicit -D always wins)
	foreach(_v GROMACS_BUILD_DIR GROMACS_SRC_DIR)
		if(NOT DEFINED ${_v} AND DEFINED ENV{${_v}})
			set(${_v} "$ENV{${_v}}")
		endif()
	endforeach()

	# GROMACS_INSTALL_DIR is derived, never accepted: honouring it is exactly how a config.h from one
	# configuration ends up paired with a libgromacs from another.
	if(DEFINED GROMACS_INSTALL_DIR)
		message(WARNING "GROMACS_INSTALL_DIR is no longer a KEnRef parameter — the install is derived from "
			"GROMACS_BUILD_DIR's CMakeCache (keeping build/source/install consistent). Ignoring the value "
			"'${GROMACS_INSTALL_DIR}'; pass -DGROMACS_BUILD_DIR=<gromacs build dir> instead.")
		unset(GROMACS_INSTALL_DIR)
		unset(GROMACS_INSTALL_DIR CACHE)
	endif()

	# Derive source + install + MPI flavor from the build dir GROMACS itself recorded.
	if(DEFINED GROMACS_BUILD_DIR AND EXISTS "${GROMACS_BUILD_DIR}/CMakeCache.txt")
		kenref_read_foreign_cache("${GROMACS_BUILD_DIR}" "CMAKE_HOME_DIRECTORY"  _kn_gmx_cache_src)
		kenref_read_foreign_cache("${GROMACS_BUILD_DIR}" "CMAKE_INSTALL_PREFIX"  _kn_gmx_cache_inst)
		if(_kn_gmx_cache_inst)
			set(GROMACS_INSTALL_DIR "${_kn_gmx_cache_inst}")
			message(STATUS "GROMACS: install derived from the build dir -> ${GROMACS_INSTALL_DIR}")
		endif()
		if(_kn_gmx_cache_src)
			if(NOT DEFINED GROMACS_SRC_DIR)
				set(GROMACS_SRC_DIR "${_kn_gmx_cache_src}")
				message(STATUS "GROMACS: source  derived from the build dir -> ${GROMACS_SRC_DIR}")
			elseif(NOT _kn_gmx_cache_src PATH_EQUAL GROMACS_SRC_DIR)
				message(WARNING "GROMACS_BUILD_DIR was configured from source '${_kn_gmx_cache_src}', but "
					"GROMACS_SRC_DIR is '${GROMACS_SRC_DIR}'. You would compile one tree's internal headers "
					"against another tree's config.h.")
			endif()
		endif()
	endif()

	# kenref-gmx compiles against GROMACS's public headers too, and GROMACS only INSTALLS them when it was
	# configured with GMX_INSTALL_LEGACY_API=ON (its default is OFF!). Without them the build dies much later
	# with a confusing "gromacs/math/vectypes.h file not found". The build dir's cache tells us up front.
	if(DEFINED GROMACS_BUILD_DIR AND EXISTS "${GROMACS_BUILD_DIR}/CMakeCache.txt")
		kenref_read_foreign_cache("${GROMACS_BUILD_DIR}" "GMX_INSTALL_LEGACY_API" _kn_gmx_legacy_api)
		if(_kn_gmx_legacy_api AND NOT _kn_gmx_legacy_api)
			message(WARNING
				"The GROMACS at '${GROMACS_BUILD_DIR}' was configured with GMX_INSTALL_LEGACY_API=OFF, so it does "
				"NOT install the public headers kenref-gmx compiles against (gromacs/math/vectypes.h, ...).\n"
				"  Reconfigure that GROMACS with -DGMX_INSTALL_LEGACY_API=ON and re-install it.")
		endif()
	endif()

	# MPI flavor comes from the GENERATED config.h — the authoritative source. NEVER inferred from the install's
	# directory/suffix name (`_mpi` is only a naming convention and can be customised).
	#
	# kenref-gmx REQUIRES a GROMACS built with a real external MPI library (GMX_LIB_MPI=1). This is not our
	# preference — it is GROMACS's own restriction: kenref-gmx drives MULTI-SIMULATION (replica) refinement, and
	# gromacs/mdrunutility/multisim.cpp throws
	#     "Multi-simulations are only supported when GROMACS has been configured with a proper external MPI
	#      library."
	# when !GMX_LIB_MPI. A thread-MPI or serial GROMACS therefore can never run kenref-gmx; detect that HERE and
	# say so, instead of failing later with a cryptic "thread_mpi/mpi_bindings.h not found" (thread-MPI) or
	# "typedef redefinition: 'ompi_communicator_t*' vs 'void*'" (serial, where gromacs shims MPI_Comm=void*).
	set(KENREF_GMX_LIB_MPI "")
	set(_kn_gmx_flavor "unknown")
	if(DEFINED GROMACS_BUILD_DIR)
		kenref_read_cpp_define("${GROMACS_BUILD_DIR}/src/include/config.h" "GMX_LIB_MPI"    KENREF_GMX_LIB_MPI)
		kenref_read_cpp_define("${GROMACS_BUILD_DIR}/src/include/config.h" "GMX_THREAD_MPI" _kn_gmx_thread_mpi)
		if(KENREF_GMX_LIB_MPI STREQUAL "1")
			set(_kn_gmx_flavor "library MPI")
			message(STATUS "GROMACS: flavor = library MPI (GMX_LIB_MPI=1) — usable by kenref-gmx")
		elseif(NOT KENREF_GMX_LIB_MPI STREQUAL "")
			if(_kn_gmx_thread_mpi STREQUAL "1")
				set(_kn_gmx_flavor "thread-MPI")
			else()
				set(_kn_gmx_flavor "serial")
			endif()
			message(STATUS "GROMACS: flavor = ${_kn_gmx_flavor} (GMX_LIB_MPI=0) — NOT usable by kenref-gmx")
		endif()
	endif()

	if(DEFINED GROMACS_SRC_DIR)
		message(STATUS "GROMACS_SRC_DIR: ${GROMACS_SRC_DIR}")
	endif()
	if(DEFINED GROMACS_BUILD_DIR)
		message(STATUS "GROMACS_BUILD_DIR: ${GROMACS_BUILD_DIR}")
	endif()
endif()

# The GROMACS CMake package dir. GROMACS names it share/cmake/gromacs<GMX_LIB_SUFFIX>: `_mpi` for a typical
# MPI build, EMPTY for a serial/thread-MPI one, or a custom suffix (e.g. `_d`). DETECT it rather than assuming
# `_mpi` — but note the name is only used to LOCATE the package, never to decide the MPI flavor (that comes
# from config.h above).
if(NOT DEFINED GROMACS_DIR AND DEFINED GROMACS_INSTALL_DIR)
	file(GLOB _kn_gmx_cfgdirs "${GROMACS_INSTALL_DIR}/share/cmake/gromacs*")
	foreach(_d IN LISTS _kn_gmx_cfgdirs)
		if(IS_DIRECTORY "${_d}")
			file(GLOB _kn_cfg "${_d}/gromacs*-config.cmake")
			if(_kn_cfg)
				set(GROMACS_DIR "${_d}")
				break()
			endif()
		endif()
	endforeach()
	unset(_kn_gmx_cfgdirs)
	unset(_kn_cfg)
endif()

# ============================================================================
# COMPONENT RESOLUTION (tri-state ON/OFF/AUTO -> boolean of the same name)
# ============================================================================
# Probe engine availability WITHOUT downloading, then resolve each tri-state. AUTO builds only if the engine
# is already present; the download switch (below, in the provide_* steps) only fetches sub-deps for a
# committed component. See cmake/Components.cmake.
message(STATUS "Resolving KEnRef components (ON/OFF/AUTO; downloads=${KENREF_FETCH_MISSING}):")

# CORE — the base library; always "available" (Eigen is a small sub-dep fetched under the download switch).
kenref_resolve_component("kenref_core" "${_KN_CORE_REQ}" TRUE BUILD_KENREF_CORE)

# GMX — "available" must mean USABLE, not merely findable. An INSTALLED GROMACS alone can never satisfy
# kenref-gmx (no internal MDModules headers, no generated config.h), so a bare find_package() success is not a
# sufficient signal: AUTO would resolve ON and then fail on the missing internal headers. Require all three:
#   generated config.h (build tree) + an internal extension-API header (source tree) + a findable CMake package.
set(_kn_gmx_available FALSE)
if(NOT _KN_GMX_REQ STREQUAL "OFF")
	if(DEFINED GROMACS_BUILD_DIR AND EXISTS "${GROMACS_BUILD_DIR}/src/include/config.h"
	   AND DEFINED GROMACS_SRC_DIR AND EXISTS "${GROMACS_SRC_DIR}/src/gromacs/mdtypes/iforceprovider.h"
	   AND DEFINED GROMACS_DIR AND IS_DIRECTORY "${GROMACS_DIR}"
	   AND KENREF_GMX_LIB_MPI STREQUAL "1")          # multi-sim needs real MPI — see the flavor note above
		find_package(GROMACS NAMES gromacs_mpi gromacs QUIET)
		if(GROMACS_FOUND)
			set(_kn_gmx_available TRUE)
		endif()
	endif()
endif()
# A non-lib-MPI GROMACS can never run kenref-gmx. Under AUTO that is just "not available" (skip quietly);
# when the user FORCED =ON, say exactly why rather than letting the compile fail cryptically later.
if(_KN_GMX_REQ STREQUAL "ON" AND NOT _kn_gmx_available
   AND NOT KENREF_GMX_LIB_MPI STREQUAL "" AND NOT KENREF_GMX_LIB_MPI STREQUAL "1")
	message(FATAL_ERROR
		"BUILD_KENREF_GMX=ON, but the GROMACS at '${GROMACS_BUILD_DIR}' is a ${_kn_gmx_flavor} build "
		"(GMX_LIB_MPI=0).\n"
		"  kenref-gmx runs MULTI-SIMULATION (replica) refinement, and GROMACS supports multi-simulations ONLY "
		"when built with a real external MPI library (see gromacs mdrunutility/multisim.cpp).\n"
		"  Rebuild that GROMACS with -DGMX_MPI=ON (and -DGMX_INSTALL_LEGACY_API=ON), then point "
		"-DGROMACS_BUILD_DIR at its build dir.")
endif()
kenref_resolve_component("kenref-gmx" "${_KN_GMX_REQ}" ${_kn_gmx_available} BUILD_KENREF_GMX)

# PLUMED interface (compile-validation) — available iff PLUMED headers are present at PLUMED_SRC_DIR.
set(_kn_plumed_iface_available FALSE)
if(DEFINED PLUMED_SRC_DIR AND NOT PLUMED_SRC_DIR STREQUAL "" AND EXISTS "${PLUMED_SRC_DIR}/bias/Bias.h")
	set(_kn_plumed_iface_available TRUE)
endif()
kenref_resolve_component("kenref_plumed (validate)" "${_KN_PLUMED_REQ}" ${_kn_plumed_iface_available} BUILD_KENREF_PLUMED)

# Delegated full PLUMED build — available (for AUTO) iff a PLUMED source tree is provided.
set(_kn_with_plumed_available FALSE)
if(DEFINED KENREF_PLUMED_SRC_DIR AND NOT KENREF_PLUMED_SRC_DIR STREQUAL "")
	set(_kn_with_plumed_available TRUE)
endif()
kenref_resolve_component("PLUMED (delegated build)" "${_KN_WITHPLUMED_REQ}" ${_kn_with_plumed_available} KENREF_WITH_PLUMED)

# Implications: batching a GROMACS implies building PLUMED; a PLUMED build here compiles plumedinterface from
# the kenref install, so it must be exported.
if(KENREF_WITH_PLUMED_GROMACS)
	set(KENREF_WITH_PLUMED ON)
endif()
if(KENREF_WITH_PLUMED)
	set(KENREF_EXPORT_PLUMEDINTERFACE ON CACHE BOOL "required by KENREF_WITH_PLUMED" FORCE)
endif()

# ============================================================================
# INSTALLATION PATHS
# ============================================================================
# ONE common install prefix for EVERYTHING from the KEnRef project (core libs + gmx executables +
# plumedinterface export) — flat, like GROMACS's /usr/local/gromacs. Default `${KENREF_INSTALL_BASE}/kenref`
# regardless of which components are built (CMake has one CMAKE_INSTALL_PREFIX per build; with a single prefix
# there is no duplication and no license/co-mingling warning to give). KENREF_INSTALL_BASE (default /usr/local)
# is the overridable system folder; or pass -DCMAKE_INSTALL_PREFIX to set the full path directly.
# Installing under /usr/local needs write permission (sudo) — build.sh handles that; for a no-sudo build pass
# -DCMAKE_INSTALL_PREFIX=<somewhere-writable> or -DKENREF_INSTALL_BASE=$HOME/.local.
set(KENREF_INSTALL_BASE "/usr/local" CACHE PATH "System folder for KEnRef's default install prefix")
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX "${KENREF_INSTALL_BASE}/kenref" CACHE PATH "root directory to install KEnRef" FORCE)
endif()

#set(CMAKE_MODULE_PATH /usr/local/gromacs-dev/share/gromacs/template/cmake:$CMAKE_MODULE_PATH)

# Dependency discovery (standard find -> fetch fallback; no hardcoded machine paths).
# Pass hints via -DCMAKE_PREFIX_PATH=... / -DEigen3_DIR=... when a dep lives outside default locations.
include(cmake/Dependencies.cmake)

# ============================================================================
# COMMON DEPENDENCIES
# ============================================================================
# MPI is required ONLY by the gmx component. kenref_core contains NO MPI code at all, and the PLUMED bias
# talks to replicas through PLUMED's own Communicator/multi_sim_comm (degrading to 1 replica when
# Communicator::plumedHasMPI() is false) and inherits PLUMED's MPI flags when PLUMED compiles it in-tree.
# Only src/gmxinterface + src/tools include <mpi.h>. So a core-only or core+plumed build needs no MPI.
if(BUILD_KENREF_GMX)
	kenref_provide_mpi()   # MPI::MPI_CXX / MPI::MPI_C — for the gmx interface + tools
endif()

if(BUILD_KENREF_CORE)
	kenref_provide_eigen()   # find Eigen3 (or fetch it if missing); sets Eigen3::Eigen + EIGEN3_INCLUDE_DIR
endif()

if(BUILD_KENREF_GMX)
	message(STATUS "Searching for Gromacs in ${GROMACS_DIR}")
	kenref_provide_fftw()      # detect system FFTW (informational; the gromacs sub-build can build its own)
	kenref_provide_gromacs()   # find_package(GROMACS) honoring GROMACS_DIR; fetch scaffold if missing
endif()

if(BUILD_KENREF_CORE)
	kenref_provide_openmp()   # optional: OpenMP::OpenMP_CXX if present, else kernels compile to serial
endif()

# ============================================================================
# INCLUDE COMPONENT MODULES
# ============================================================================

if(BUILD_KENREF_CORE)
	include(cmake/KEnRefCore.cmake)
endif()

if(BUILD_KENREF_GMX)
	# If building GMX without CORE, need to find installed CORE
	if(NOT BUILD_KENREF_CORE AND NOT TARGET KENREF::CORE)
		message(STATUS "Searching for installed KENREF CORE...")
		# Mirror the core install default (KENREF_INSTALL_BASE/kenref); override with
		# -DKENREF_CORE_CMAKE_PATH=<prefix>/lib/cmake/KEnRef or -DCMAKE_PREFIX_PATH=<core-prefix>.
		if(NOT DEFINED KENREF_CORE_CMAKE_PATH)
			set(KENREF_CORE_CMAKE_PATH "${KENREF_INSTALL_BASE}/kenref/lib/cmake/KEnRef")
		endif()
		message(STATUS "  Looking in: ${KENREF_CORE_CMAKE_PATH}")

		# KENREF CORE depends on Eigen3, so provide it first (find or fetch).
		kenref_provide_eigen()

		if(EXISTS "${KENREF_CORE_CMAKE_PATH}/KEnRefCoreConfig.cmake")
			message(STATUS "  Found KEnRefCoreConfig.cmake, including it...")
			include("${KENREF_CORE_CMAKE_PATH}/KEnRefCoreConfig.cmake")

			# The exported target is KENREF::kenref_core, create an alias to KENREF::CORE
			if(TARGET KENREF::kenref_core AND NOT TARGET KENREF::CORE)
				add_library(KENREF::CORE ALIAS KENREF::kenref_core)
			endif()

			message(STATUS "  ✓ Loaded KENREF CORE from ${KENREF_CORE_CMAKE_PATH}")
		else()
			message(FATAL_ERROR "KENREF CORE not found at ${KENREF_CORE_CMAKE_PATH}\nPlease build and install KENREF CORE first, or build with -DBUILD_KENREF_CORE=ON")
		endif()
	endif()

	include(cmake/KEnRefGMX.cmake)
endif()

if(BUILD_KENREF_PLUMED)
	# Needs KENREF::CORE (build it here, or it must be found as above for the GMX path).
	if(NOT TARGET KENREF::CORE)
		message(FATAL_ERROR "BUILD_KENREF_PLUMED requires KENREF::CORE (build with -DBUILD_KENREF_CORE=ON).")
	endif()
	include(cmake/KEnRefPlumed.cmake)
endif()

# Post-install: env.sh + TCL modulefile (always) + a clear summary of where things went.
if(BUILD_KENREF_CORE OR BUILD_KENREF_GMX)
	include(cmake/PostInstall.cmake)
endif()

# PLUMED delegation, invoked by CMake at install time (AFTER PostInstall so env.sh is installed first).
if(KENREF_WITH_PLUMED)
	if(NOT BUILD_KENREF_CORE)
		message(FATAL_ERROR "KENREF_WITH_PLUMED needs kenref_core (build with -DBUILD_KENREF_CORE=ON).")
	endif()
	include(cmake/InvokePlumed.cmake)
endif()

# ============================================================================
# TESTING
# ============================================================================
option(BUILD_TESTS "Build Google tests" OFF)
if(BUILD_TESTS)
	enable_testing()
	# KEnRefGMX.cmake also adds google_tests when the gmx interface is built; avoid a double add_subdirectory.
	if(NOT TARGET Google_tests_kenref_core_exe)
		add_subdirectory(google_tests)
	endif()
endif()

# ============================================================================
# PKG-CONFIG GENERATION FOR PLUMED INTEGRATION
# ============================================================================

# Get MPI flags for pkg-config
if(TARGET MPI::MPI_CXX)
    get_target_property(MPI_CXX_COMPILE_OPTIONS MPI::MPI_CXX INTERFACE_COMPILE_OPTIONS)
    get_target_property(MPI_CXX_LINK_LIBRARIES MPI::MPI_CXX INTERFACE_LINK_LIBRARIES)

    if(MPI_CXX_COMPILE_OPTIONS)
        string(JOIN " " MPI_CXX_COMPILE_OPTIONS_STR ${MPI_CXX_COMPILE_OPTIONS})
    else()
        set(MPI_CXX_COMPILE_OPTIONS_STR "")
    endif()
    if(MPI_CXX_LINK_LIBRARIES)
        string(JOIN " " MPI_CXX_LIBRARIES_STR ${MPI_CXX_LINK_LIBRARIES})
    else()
        set(MPI_CXX_LIBRARIES_STR "")
    endif()
endif()

# ============================================================================
# COMMON INSTALLATION SETTINGS
# ============================================================================

# Add install RPATH
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)