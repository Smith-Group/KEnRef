# ============================================================================
# RunEigenAbiGen.cmake — run the Eigen-ABI generator, but NEVER break the build if it cannot run.
#
# cmake/eigen_abi_gen.cpp is compiled with kenref_core's own flags and EXECUTED at build time to record the
# Eigen ABI kenref_core is built with. Executing a just-built binary is not always possible: with e.g.
# -stdlib=libc++ the generator needs the toolchain's runtime libs loadable at build time (no LD_LIBRARY_PATH
# -> exec fails with 127), and cross/sandboxed builds may refuse outright.
#
# The compile-time guard (core/EigenAbiCheck.h) is only a BELT — the RUNTIME check (core/SimdSignature.h,
# forced at install) is the primary defence and needs no generator. So if the generator cannot run we emit a
# header WITHOUT the KENREF_CORE_EIGEN_* defines: EigenAbiCheck.h then compiles to a no-op (it guards on
# `#ifdef KENREF_CORE_EIGEN_MAX_ALIGN_BYTES`) and the build proceeds with the runtime check intact.
#
# Invoked as: cmake -DKN_GEN=<generator exe> -DKN_OUT=<header path> -P RunEigenAbiGen.cmake
# ============================================================================
get_filename_component(_kn_out_dir "${KN_OUT}" DIRECTORY)
file(MAKE_DIRECTORY "${_kn_out_dir}")

execute_process(COMMAND "${KN_GEN}" "${KN_OUT}"
                RESULT_VARIABLE _kn_rc
                OUTPUT_QUIET ERROR_VARIABLE _kn_err)

if(NOT _kn_rc EQUAL 0)
    file(WRITE "${KN_OUT}"
        "// GENERATED: the Eigen-ABI generator could not be executed in this build environment, so the\n"
        "// COMPILE-TIME kenref<->Eigen guard is disabled here (core/EigenAbiCheck.h no-ops without these\n"
        "// defines). The RUNTIME guard (core/SimdSignature.h, enforced at install) is unaffected.\n"
        "#ifndef KENREF_CORE_EIGEN_ABI_H\n#define KENREF_CORE_EIGEN_ABI_H\n#endif\n")
    message(WARNING
        "KEnRef: could not run the Eigen-ABI generator (${KN_GEN}): ${_kn_err}\n"
        "  -> the COMPILE-TIME kenref<->Eigen ABI check is disabled for this build.\n"
        "     (The runtime SIMD/Eigen check still runs at install.) A common cause is the compiler's runtime\n"
        "     libraries not being loadable at build time (e.g. -stdlib=libc++ without LD_LIBRARY_PATH).")
endif()
