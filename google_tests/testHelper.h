//
// Created by amr on 5/16/24.
//

#ifndef KENREF_TESTHELPER_H
#define KENREF_TESTHELPER_H

#include "core/KEnRef.h"
#include <Eigen/src/Core/Matrix.h>

void EXPECT_MATRIX_NEAR(const Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, Eigen::Dynamic> &leftSide,
                        const Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, Eigen::Dynamic> &rightSide,
                        KEnRef_Real_t epsilon =
#ifdef DOUBLE
                        1.0e-14
#else
                        5.0e-6
#endif
);

#endif //KENREF_TESTHELPER_H
