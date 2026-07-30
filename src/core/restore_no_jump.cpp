/*
 * restore_no_jump.cpp
 *
 *  Engine-neutral port of KEnRefForceProvider::restoreNoJump (which KEnRefBias::restoreNoJump
 *  mirrored byte-for-byte). The GROMACS `matrix`/`msmul`/XX,YY,ZZ,DIM are replaced by an Eigen 3x3
 *  and explicit indices; the algorithm is otherwise unchanged, so output is bitwise-identical.
 */

#include "core/restore_no_jump.h"

#include <vector>
#include <iostream>

namespace kenref {

template<typename Real>
void restoreNoJump(CoordsMatrixType<Real>& atoms,
                   const CoordsMatrixType<Real>& reference,
                   const Eigen::Matrix<Real, 3, 3>& box_raw,
                   const bool toAngstrom, int numOmpThreads, const bool printStatistics) {
    // Scale the box to the working units (Angstrom) exactly as the original msmul(box_, 10, box).
    const Real scale_factor = toAngstrom ? Real(10.0) : Real(1.0);
    const Eigen::Matrix<Real, 3, 3> box = box_raw * scale_factor;

    // Precompute box vectors and half the box length in each dimension
    const Eigen::RowVector3<Real> box_half =
        Real(0.5) * Eigen::RowVector3<Real>(box(0, 0), box(1, 1), box(2, 2));
    // Precompute non-zero dimensions to avoid branching in inner loop
    std::vector<int> active_dims;
    for (int m = 2 /*ZZ*/; m >= 0 /*XX*/; --m) {
        if (box_half[m] != 0) {
            active_dims.push_back(m);
        }
    }

    Eigen::RowVectorXi updatedLocations;
    bool updated = false;
    if (printStatistics) {
        updatedLocations = Eigen::RowVectorXi::Zero(atoms.rows());
    }

// default(shared), not default(none): the explicit shared() list below still documents every variable
// this region touches, but default(none) is not portable here. Eigen's templates reference their own
// constexpr `Eigen::Dynamic` inside the region, and under default(none) a pre-OpenMP-5.0 compiler
// demands an explicit data-sharing attribute for it -- which cannot be given, since it is a constexpr
// in another namespace, not one of our variables. OpenMP 5.0 made such constants predetermined shared,
// so newer clang accepts it and older clang does not:
//     Eigen/src/Core/util/ForwardDeclarations.h:88:45:
//     error: variable 'Dynamic' must have explicitly specified data sharing attributes
// Caught by the CI matrix (clang), not by the local build, whose clang is new enough to allow it.
#pragma omp parallel for num_threads(numOmpThreads) default(shared) reduction(||:updated) \
    shared(atoms, box, box_half, reference, updatedLocations, printStatistics, active_dims)
    for (int i = 0; i < atoms.rows(); ++i) {
        bool local_updated = false;
        // Get row references once per atom
        auto atom = atoms.row(i);
        const auto &refAtom = reference.row(i);

        //Iterate only over active dimensions
        for (int m : active_dims) {
            // Check if atom jumped across the box in this dimension
            while (atom[m] - refAtom[m] <= -box_half[m]) {
                // Jumped to negative image, correct by adding box size
                for (int d = 0; d <= m; ++d) {
                    atom[d] += box(m, d);
                }
                if (printStatistics) {
                    updatedLocations(i) = 1;
                    local_updated = true;
                }
            }
            while ((atom[m] - refAtom[m]) > box_half[m]) {
                // Jumped to positive image, correct by subtracting box size
                for (int d = 0; d <= m; ++d) {
                    atom[d] -= box(m, d);
                }
                if (printStatistics) {
                    updatedLocations(i) = 1;
                    local_updated = true;
                }
            }
        }
        updated = updated || local_updated; // this is faster than `updated |= local_updated`
    }

    if (printStatistics && updated) {
        const int updated_count = updatedLocations.sum();
        std::cout << "INFO: Restored NoJump in " << updated_count << " atoms (out of "
                  << updatedLocations.size() << ")" << std::endl;
        if (updated_count < 50) {  // Only show details for small numbers
            std::cout << "Updated atom indices:\n" << updatedLocations << std::endl;
        }
    }
}

// Explicit instantiations (the build fixes the scalar via KENREF_DOUBLE, but provide both like the
// rest of the core, e.g. Table.cpp).
template void restoreNoJump<double>(CoordsMatrixType<double>&, const CoordsMatrixType<double>&,
                                    const Eigen::Matrix<double, 3, 3>&, bool, int, bool);
template void restoreNoJump<float>(CoordsMatrixType<float>&, const CoordsMatrixType<float>&,
                                   const Eigen::Matrix<float, 3, 3>&, bool, int, bool);

} // namespace kenref
