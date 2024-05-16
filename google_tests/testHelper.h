//
// Created by amr on 5/16/24.
//

#ifndef KENREF_TESTHELPER_H
#define KENREF_TESTHELPER_H

#include "core/KEnRef.h"
#include <limits>
#include <type_traits>
#include <Eigen/src/Core/Matrix.h>

template<typename KEnRef_Real>
class TestHelper{
public:
    static void EXPECT_MATRIX_NEAR(const Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> &leftSide,
                                   const Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> &rightSide,
                                   const KEnRef_Real epsilon = std::is_same_v<KEnRef_Real_t, float> ? 5e-6 :
                                                                 std::is_same_v<KEnRef_Real_t, float> ? 1e-14 :
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

};

template class TestHelper<float>;
template class TestHelper<double>;

#endif //KENREF_TESTHELPER_H
