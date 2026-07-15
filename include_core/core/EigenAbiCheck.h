/*
 * EigenAbiCheck.h — COMPILE-TIME kenref<->Eigen ABI guard.
 *
 * kenref_core stores Eigen objects directly in its own containers (e.g. std::vector<CoordsMatrixType>, and
 * Eigen members in structs). Eigen's alignment (EIGEN_MAX_ALIGN_BYTES) is therefore baked into the layout and
 * the aligned allocations inside libkenref_core — it is part of kenref_core's ABI. If a consumer compiles
 * Eigen (its own, or the bundled kenref_and_eigen3 headers) with a DIFFERENT alignment — from a different
 * -march/ACCEL, a different Eigen version, or an Eigen config macro — the two disagree on those objects and
 * memory is corrupted at runtime.
 *
 * This header static_asserts that the including TU's Eigen alignment (and ABI epoch) match the values
 * kenref_core was actually built with, recorded in the GENERATED core/KEnRefEigenAbi.h (see
 * cmake/eigen_abi_gen.cpp). It is included from core/KEnRef.h, so every consumer of the KEnRef API is checked
 * at compile time. Complements the runtime gate in core/SimdSignature.h (which also covers link-only
 * consumers and the cross-ACCEL install case). When the generated header is absent (a partial header copy) or
 * the build cross-compiles, the check degrades to a no-op rather than a false failure.
 */
#ifndef KENREF_CORE_EIGENABICHECK_H
#define KENREF_CORE_EIGENABICHECK_H

#include <Eigen/Core>   // EIGEN_MAX_ALIGN_BYTES / EIGEN_MAX_STATIC_ALIGN_BYTES / EIGEN_WORLD_VERSION

#if defined(__has_include)
#  if __has_include("core/KEnRefEigenAbi.h")
#    include "core/KEnRefEigenAbi.h"
#  endif
#endif

#if defined(KENREF_CORE_EIGEN_MAX_ALIGN_BYTES) && !defined(KENREF_NO_EIGEN_ABI_CHECK)
static_assert(EIGEN_MAX_ALIGN_BYTES == KENREF_CORE_EIGEN_MAX_ALIGN_BYTES,
    "KEnRef<->Eigen ABI mismatch: your Eigen build produces a different EIGEN_MAX_ALIGN_BYTES than kenref_core "
    "was built with. kenref_core stores Eigen objects in its containers, so mixing the two WILL corrupt memory. "
    "Rebuild your code with the SAME ACCEL/-march and the SAME Eigen as kenref_core (see INSTALL.md, SIMD/Eigen "
    "ABI). Set -DKENREF_NO_EIGEN_ABI_CHECK to bypass this at your own risk.");
static_assert(EIGEN_MAX_STATIC_ALIGN_BYTES == KENREF_CORE_EIGEN_MAX_STATIC_ALIGN_BYTES,
    "KEnRef<->Eigen ABI mismatch: EIGEN_MAX_STATIC_ALIGN_BYTES differs from kenref_core's. See INSTALL.md.");
// NB: compare WORLD *and* MAJOR. Eigen's release number and its version macros diverge — release "5.0.1"
// defines EIGEN_WORLD_VERSION=3 / EIGEN_MAJOR_VERSION=5 (macro version 3.5.0), exactly as 3.4.0 defines 3/4/0.
// So WORLD alone would NOT tell Eigen 3.4 from Eigen "5.x" (both WORLD=3); MAJOR is what discriminates.
static_assert(EIGEN_WORLD_VERSION == KENREF_CORE_EIGEN_WORLD_VERSION
              && EIGEN_MAJOR_VERSION == KENREF_CORE_EIGEN_MAJOR_VERSION,
    "KEnRef<->Eigen ABI mismatch: your Eigen version series differs from the one kenref_core was built with "
    "(Eigen keeps its ABI within a MINOR series, not across MAJOR). See INSTALL.md (SIMD/Eigen ABI).");
#endif

#endif // KENREF_CORE_EIGENABICHECK_H
