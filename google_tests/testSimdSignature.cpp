/*
 * testSimdSignature.cpp — runtime assertion that this test TU's Eigen alignment matches the linked
 * kenref_core's. In the core test exe (same flags as the lib) this is trivially true; in the gmx test exe
 * (compiled with the gmx build's flags, linking core) it catches a cross-ACCEL mismatch.
 */
#include "gtest/gtest.h"
#include "core/SimdSignature.h"

TEST(SimdSignatureTest, CoreMatchesCaller) {
    std::string why;
    EXPECT_TRUE(kenref::simd_matches(&why)) << "kenref_core / consumer SIMD (Eigen ABI) mismatch: " << why;
}

TEST(SimdSignatureTest, CoreSignatureIsSane) {
    const kenref::SimdSignature s = kenref::core_simd_signature();
    EXPECT_GT(s.max_align, 0);
    EXPECT_GT(s.max_static_align, 0);
    EXPECT_GT(s.eigen_abi, 0);
    EXPECT_FALSE(s.isa.empty());
    EXPECT_FALSE(s.eigen_ver.empty());
}
