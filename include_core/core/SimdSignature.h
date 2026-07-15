/*
 * SimdSignature.h — runtime guard against a kenref_core / consumer SIMD (Eigen ABI) mismatch.
 *
 * Eigen is header-only and vectorizes per translation unit according to that TU's -march. If kenref_core is
 * built with one ACCEL/-march and a consumer (kenref-gmx, the PLUMED bias, …) with another, their
 * EIGEN_MAX_ALIGN_BYTES can differ — and passing Eigen objects across that boundary is an ABI mismatch that
 * corrupts alignment and crashes. Compilation alone cannot detect this (each TU compiles fine on its own);
 * only a RUNTIME comparison of a separately-built kenref_core against the consumer catches it.
 *
 *   core_simd_signature()  — DEFINED in SimdSignature.cpp, so it bakes kenref_core's OWN alignment + ISA.
 *   local_simd_signature() — inline, so it reflects the INCLUDING TU's alignment + ISA.
 *   simd_matches()         — true iff the two ALIGNMENT classes match (the ABI truth). The ISA string is
 *                            diagnostic only (a tier name may differ while the alignment class is identical).
 */
#ifndef KENREF_CORE_SIMDSIGNATURE_H
#define KENREF_CORE_SIMDSIGNATURE_H

#include <string>
#include <Eigen/Core>   // EIGEN_MAX_ALIGN_BYTES / EIGEN_MAX_STATIC_ALIGN_BYTES / EIGEN_*_VERSION

// The ACCEL/SIMD tier name, injected per-build via -DKENREF_ACCEL="..." (see CMakeLists.txt). Not exported as
// an INTERFACE definition, so a consumer TU records ITS OWN build's tier (or "unknown" if it set none).
#ifndef KENREF_ACCEL
#define KENREF_ACCEL "unknown"
#endif

// The Eigen ABI epoch. CAREFUL: Eigen's RELEASE number and its VERSION MACROS diverge — release "5.0.1"
// defines EIGEN_WORLD_VERSION=3, EIGEN_MAJOR_VERSION=5 (i.e. macro version 3.5.0), exactly like 3.4.0 defines
// 3/4/0. So WORLD alone does NOT distinguish Eigen 3.4 from Eigen "5.x" — both are WORLD=3. The meaningful
// epoch is WORLD.MAJOR (3.4 vs 3.5); Eigen keeps its ABI within a MINOR series, not across MAJOR.
#define KENREF_EIGEN_ABI_EPOCH  (EIGEN_WORLD_VERSION * 100 + EIGEN_MAJOR_VERSION)
#define KENREF_EIGEN_VERSION_STR \
    (std::to_string(EIGEN_WORLD_VERSION) + "." + std::to_string(EIGEN_MAJOR_VERSION) + "." + std::to_string(EIGEN_MINOR_VERSION))

namespace kenref {

struct SimdSignature {
    int         max_align;         // EIGEN_MAX_ALIGN_BYTES        — the Eigen ABI alignment (set by -march)
    int         max_static_align;  // EIGEN_MAX_STATIC_ALIGN_BYTES
    int         eigen_abi;         // WORLD*100+MAJOR              — Eigen ABI epoch (see the note above)
    std::string isa;               // ACCEL/SIMD tier name (diagnostic)
    std::string eigen_ver;         // full Eigen "W.M.m" (diagnostic)

    // ABI equality = same Eigen alignment AND same Eigen ABI epoch. isa / full version are diagnostic only.
    bool operator==(const SimdSignature& o) const {
        return max_align == o.max_align && max_static_align == o.max_static_align && eigen_abi == o.eigen_abi;
    }
    bool operator!=(const SimdSignature& o) const { return !(*this == o); }
};

// Baked into kenref_core (compiled with the library's flags + Eigen headers).
SimdSignature core_simd_signature();

// Evaluated with the including translation unit's flags + Eigen headers.
inline SimdSignature local_simd_signature() {
    return SimdSignature{ EIGEN_MAX_ALIGN_BYTES, EIGEN_MAX_STATIC_ALIGN_BYTES, KENREF_EIGEN_ABI_EPOCH,
                          KENREF_ACCEL, KENREF_EIGEN_VERSION_STR };
}

// True iff the caller's Eigen alignment matches the linked kenref_core's. On mismatch, *why (if non-null)
// gets a human-readable diff of both signatures.
inline bool simd_matches(std::string* why = nullptr) {
    const SimdSignature core = core_simd_signature();
    const SimdSignature loc  = local_simd_signature();
    if (core == loc) {
        return true;
    }
    if (why) {
        *why = "kenref_core[align=" + std::to_string(core.max_align) + "/" +
               std::to_string(core.max_static_align) + " eigen=" + core.eigen_ver + " isa=" + core.isa +
               "] != caller[align=" + std::to_string(loc.max_align) + "/" +
               std::to_string(loc.max_static_align) + " eigen=" + loc.eigen_ver + " isa=" + loc.isa + "]";
    }
    return false;
}

} // namespace kenref

#endif // KENREF_CORE_SIMDSIGNATURE_H
