//
// Created by amr on 5/16/24.
//

#ifndef KENREF_TESTHELPER_H
#define KENREF_TESTHELPER_H

#include "core/KEnRef.h"
#include <limits>
#include <numeric>
#include <type_traits>
#include <vector>
#include <Eigen/src/Core/Matrix.h>

template<typename KEnRef_Real>
class TestHelper{
public:
    static void EXPECT_MATRIX_NEAR(const Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> &leftSide,
                                   const Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> &rightSide,
                                   const KEnRef_Real epsilon = std::is_same_v<KEnRef_Real_t, float> ? 5e-6 :
                                                                 std::is_same_v<KEnRef_Real_t, double> ? 1e-14 :
                                                                 std::numeric_limits<KEnRef_Real>::min()
    ){
        EXPECT_EQ(leftSide.rows(), rightSide.rows()) << "Matrix rows do not match.";
        EXPECT_EQ(leftSide.cols(), rightSide.cols()) << "Matrix columns do not match.";
        for (auto i = 0; i < leftSide.rows(); i++) {
            for (int j = 0; j < leftSide.cols(); ++j) {
                if(abs(leftSide(i, j) - rightSide(i, j)) > epsilon)
                    EXPECT_NEAR(leftSide(i, j), rightSide(i, j), epsilon) << "Matrices differ at (" << i << ", " << j << ")";
            }
        }
    }

    // Central finite-difference gradient. energy_fn may modify its argument
    // (e.g. coord_array_to_sigma_energy scales in-place); two fresh copies are
    // made per element. atom_row_ids empty = perturb all rows.
    template<typename EnergyFn>
    static std::vector<CoordsMatrixType<KEnRef_Real>> finite_difference_grad(
        const std::vector<CoordsMatrixType<KEnRef_Real>>& coord_array,
        EnergyFn energy_fn,
        const std::vector<int>& atom_row_ids,
        KEnRef_Real delta = static_cast<KEnRef_Real>(1e-5)
    ) {
        const int numModels = static_cast<int>(coord_array.size());
        std::vector<CoordsMatrixType<KEnRef_Real>> fd_grad(numModels);
        for (int m = 0; m < numModels; ++m)
            fd_grad[m] = CoordsMatrixType<KEnRef_Real>::Zero(coord_array[m].rows(), coord_array[m].cols());

        std::vector<int> rows;
        if (atom_row_ids.empty()) {
            rows.resize(coord_array[0].rows());
            std::iota(rows.begin(), rows.end(), 0);
        } else {
            rows = atom_row_ids;
        }

        for (int m = 0; m < numModels; ++m) {
            const int nCols = static_cast<int>(coord_array[m].cols());
            for (int i : rows) {
                for (int j = 0; j < nCols; ++j) {
                    auto coords_plus = coord_array;
                    coords_plus[m](i, j) += delta;
                    const KEnRef_Real e_plus = energy_fn(coords_plus);

                    auto coords_minus = coord_array;
                    coords_minus[m](i, j) -= delta;
                    const KEnRef_Real e_minus = energy_fn(coords_minus);

                    fd_grad[m](i, j) = (e_plus - e_minus) / (2 * delta);
                }
            }
        }
        return fd_grad;
    }

    static void EXPECT_MATRIX_EQ(const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> &leftSide,
                                   const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> &rightSide){
    EXPECT_EQ(leftSide.rows(), rightSide.rows()) << "Matrix rows do not match.";
    EXPECT_EQ(leftSide.cols(), rightSide.cols()) << "Matrix columns do not match.";
    for (auto i = 0; i < leftSide.rows(); i++) {
        for (int j = 0; j < leftSide.cols(); ++j) {
            EXPECT_EQ(leftSide(i, j), rightSide(i, j)) << "Matrices differ at (" << i << ", " << j << ")";
        }
    }
    }
};

template class TestHelper<float>;
template class TestHelper<double>;

#endif //KENREF_TESTHELPER_H
