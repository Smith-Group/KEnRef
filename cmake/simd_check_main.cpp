/*
 * simd_check_main.cpp — tiny gate executable. Compiled with a given build's flags and linked against a
 * kenref_core; returns 0 iff the caller's Eigen alignment matches the linked core's, non-zero otherwise.
 * Wired into an install(CODE) step so `cmake --install` FATALs on an ISA/ABI mismatch (see KEnRefCore.cmake
 * and KEnRefGMX.cmake). Trivially passes in a single-configure build; bites when a consumer links a
 * pre-installed core built with a different ACCEL.
 */
#include "core/SimdSignature.h"

#include <cstdio>
#include <string>

int main() {
    std::string why;
    if (kenref::simd_matches(&why)) {
        const kenref::SimdSignature s = kenref::core_simd_signature();
        std::printf("kenref SIMD/Eigen ABI OK (align=%d/%d eigen=%s isa=%s)\n",
                    s.max_align, s.max_static_align, s.eigen_ver.c_str(), s.isa.c_str());
        return 0;
    }
    std::fprintf(stderr, "kenref SIMD/Eigen ABI MISMATCH: %s\n", why.c_str());
    return 1;
}
