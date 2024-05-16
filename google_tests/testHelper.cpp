//
// Created by amr on 5/16/24.
//

#include "core/KEnRef.h"
#include <Eigen/src/Core/Matrix.h>
#include "gtest/gtest.h"

void EXPECT_MATRIX_NEAR(const Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, Eigen::Dynamic> &leftSide,
                        const Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, Eigen::Dynamic> &rightSide,
                        KEnRef_Real_t epsilon =
#ifdef DOUBLE
                        1.0e-14
#else
                        5.0e-6
#endif
) {
    EXPECT_EQ(leftSide.rows(), rightSide.rows()) << "Matrix rows do not match.";
    EXPECT_EQ(leftSide.cols(), rightSide.cols()) << "Matrix columns do not match.";

    for (auto i = 0; i < leftSide.rows(); i++) {
        for (int j = 0; j < leftSide.cols(); ++j) {
            if(abs(leftSide(i, j) - rightSide(i, j)) > epsilon)
                EXPECT_NEAR(leftSide(i, j), rightSide(i, j), epsilon) << "Matrices differ at (" << i << ", " << j << ")";
        }
    }
}
