/*
 * SimdSignature.cpp — the ONE translation unit that bakes kenref_core's own SIMD/Eigen signature.
 * Compiled with the library's flags, so core_simd_signature() reports the alignment + ACCEL that
 * kenref_core was actually built with. See core/SimdSignature.h for the rationale.
 */
#include "core/SimdSignature.h"

#include <Eigen/Core>

namespace kenref {

SimdSignature core_simd_signature() {
    return SimdSignature{ EIGEN_MAX_ALIGN_BYTES, EIGEN_MAX_STATIC_ALIGN_BYTES, KENREF_EIGEN_ABI_EPOCH,
                          KENREF_ACCEL, KENREF_EIGEN_VERSION_STR };
}

} // namespace kenref
