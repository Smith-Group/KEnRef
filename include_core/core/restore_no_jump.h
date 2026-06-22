/*
 * restore_no_jump.h
 *
 *  Engine-neutral periodic "no-jump" correction, moved into kenref_core so both the GROMACS and
 *  PLUMED consumers share one definition (they previously held byte-for-byte copies in
 *  KEnRefForceProvider::restoreNoJump and KEnRefBias::restoreNoJump).
 *
 *  The only change from the originals is the box representation: the GROMACS `matrix` is replaced
 *  by an Eigen 3x3 (box(i,j) = i-th box vector's j-th component, in raw engine units). The
 *  `toAngstrom` scale is applied here exactly as before (msmul(box, 10)), so the arithmetic — and
 *  thus the output — is identical. Each consumer's adapter converts its native box to this 3x3.
 */

#ifndef KENREF_RESTORE_NO_JUMP_H_
#define KENREF_RESTORE_NO_JUMP_H_

#include <Eigen/Dense>

namespace kenref {

/**
 * Shift each row of `atoms` by whole box vectors so it lies in the same periodic image as the
 * corresponding row of `reference` (typically the previous frame). `box_raw` is the simulation box
 * in raw engine units; `toAngstrom` scales it (and only it) by 10 to match Angstrom-based coords.
 *
 * Determinism: per-atom writes are disjoint, so coordinate output is bitwise-identical for any
 * numOmpThreads. With printStatistics, an extra reduction(||) computes only an "any atom moved"
 * flag for the log line; it does not affect coordinates.
 */
template<typename Real>
void restoreNoJump(Eigen::Matrix<Real, Eigen::Dynamic, 3, Eigen::RowMajor>& atoms,
                   const Eigen::Matrix<Real, Eigen::Dynamic, 3, Eigen::RowMajor>& reference,
                   const Eigen::Matrix<Real, 3, 3>& box_raw,
                   bool toAngstrom,
                   int numOmpThreads = 0,
                   bool printStatistics = false);

} // namespace kenref

#endif // KENREF_RESTORE_NO_JUMP_H_
