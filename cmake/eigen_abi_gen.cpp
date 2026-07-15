/*
 * eigen_abi_gen.cpp — build-time generator. Compiled with kenref_core's OWN flags (same -march + Eigen), it
 * writes core/KEnRefEigenAbi.h recording the Eigen ABI (alignment + version) kenref_core is built with. A
 * static_assert in core/EigenAbiCheck.h then fails ANY consumer TU whose Eigen produces a different value —
 * a COMPILE-TIME kenref<->Eigen guard complementing the runtime SIMD gate. See core/SimdSignature.h.
 *
 * Usage: eigen_abi_gen <output-header-path>
 */
#include <Eigen/Core>
#include <cstdio>

int main(int argc, char** argv) {
    std::FILE* f = (argc > 1) ? std::fopen(argv[1], "w") : stdout;
    if (!f) return 1;
    std::fprintf(f,
        "// GENERATED at kenref_core build time — records the Eigen ABI kenref_core was compiled with.\n"
        "#ifndef KENREF_CORE_EIGEN_ABI_H\n"
        "#define KENREF_CORE_EIGEN_ABI_H\n"
        "#define KENREF_CORE_EIGEN_MAX_ALIGN_BYTES %d\n"
        "#define KENREF_CORE_EIGEN_MAX_STATIC_ALIGN_BYTES %d\n"
        "#define KENREF_CORE_EIGEN_WORLD_VERSION %d\n"
        "#define KENREF_CORE_EIGEN_MAJOR_VERSION %d\n"
        "#endif\n",
        EIGEN_MAX_ALIGN_BYTES, EIGEN_MAX_STATIC_ALIGN_BYTES, EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION);
    if (f != stdout) std::fclose(f);
    return 0;
}
