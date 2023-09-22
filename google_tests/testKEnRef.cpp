#include "gtest/gtest.h"
#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include "KEnRef.h"

void EXPECT_MATRIX_NEAR(const Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, Eigen::Dynamic> &leftSide,
                        const Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, Eigen::Dynamic> &rightSide,
                        KEnRef_Real_t epsilon =
#ifdef DOUBLE
                                                1.0e-14
#else
                                                5.0e-6
#endif
                        ) {
    EXPECT_EQ(leftSide.rows(), rightSide.rows());
    EXPECT_EQ(leftSide.cols(), rightSide.cols());
    for (auto i = 0; i < leftSide.rows(); i++) {
        for (int j = 0; j < leftSide.cols(); ++j) {
            EXPECT_NEAR(leftSide(i, j), rightSide(i, j), epsilon);
        }
    }
}

Eigen::IOFormat fullPrecisionFmt(Eigen::FullPrecision);

TEST(KEnRefTestSuite, TestRArrayToDArray1) {

    std::vector<std::vector<std::vector<int>>> toy_grouping_list{{{0, 1, 2, 3}},
                                                                 {{0, 1}, {2, 3}},
                                                                 {{0},    {1}, {2}, {3}}};
    Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, 3> toy_r_mat(4, 3);
    toy_r_mat <<
              0.848351683690084, -0.529433112659379, 0,
            0.966177888683851, 0.257876496444355, 0,
            0.966177888683851, -0.257876496444355, 0,
            0.848351683690084, 0.529433112659379, 0;

    CoordsMatrixType<KEnRef_Real_t> expected_toy_r_mat(4, 3);
    expected_toy_r_mat <<
                       0.848351683690084, -0.529433112659379, 0,
            0.966177888683851, 0.257876496444355, 0,
            0.966177888683851, -0.257876496444355, 0,
            0.848351683690084, 0.529433112659379, 0;

    std::cout << "toy_r_mat" << std::endl << toy_r_mat << std::endl;
    std::cout << "toy_r_mat" << std::endl << toy_r_mat.format(fullPrecisionFmt) << std::endl;
    EXPECT_EQ(toy_r_mat, expected_toy_r_mat);

    auto [toy_d_array, toy_d_array_grad] = KEnRef<KEnRef_Real_t>::r_array_to_d_array(toy_r_mat, true);
    Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, 5> expected_toy_d_array(4, 5);
    Eigen::Matrix<KEnRef_Real_t, 5, Eigen::Dynamic> temp(5, 4);
    temp <<
         -0.500000000000001, -0.5, -0.5, -0.500000000000001,
            0.380532565661007, 0.750843527257811, 0.750843527257811, 0.380532565661007,
            0, 0, 0, 0,
            0, 0, 0, 0,
            -0.777942778404335, 0.431548372230798, -0.431548372230798, 0.777942778404335;
    std::cout << "toy_d_array" << std::endl << toy_d_array << std::endl;
    std::cout << "toy_d_array" << std::endl << toy_d_array.format(fullPrecisionFmt) << std::endl;
    expected_toy_d_array = temp.transpose();
    EXPECT_MATRIX_NEAR(toy_d_array, expected_toy_d_array);

    std::cout << "toy_d_array_grad" << std::endl << toy_d_array_grad << std::endl;
    std::cout << "toy_d_array_grad" << std::endl << toy_d_array_grad.format(fullPrecisionFmt) << std::endl;
    Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, Eigen::Dynamic> temp2(15, 4);
    temp2 <<
          1.27252752553513, 1.44926683302578, 1.44926683302578, 1.27252752553513,
            -0.144738995049283, -1.95377287713937, -1.95377287713937, -0.144738995049283,
            0, 0, 0, 0,
            0, 0, 0, 0,
            2.38284027903391, -1.63810728181505, 1.63810728181505, -2.38284027903391,
            -0.794149668989071, 0.386814744666533, -0.386814744666533, 0.794149668989071,
            1.92433775386622, -1.41477968485531, 1.41477968485531, -1.92433775386622,
            0, 0, 0, 0,
            0, 0, 0, 0,
            -0.589955114369631, 1.11703828096435, 1.11703828096435, -0.589955114369631,
            0, 0, 0, 0,
            0, 0, 0, 0,
            1.46938821883783, 1.67346919235006, 1.67346919235006, 1.46938821883783,
            -0.917005050335386, 0.446655193919479, -0.446655193919479, 0.917005050335386,
            0, 0, 0, 0;
    Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, Eigen::Dynamic> expected_toy_d_array_grad;
    int arr[15], idx = 0;
    for(int d = 0; d<5; d++)
        for(int xyz = 0; xyz < 3; xyz++)
            arr[idx++] = d + 5*xyz;
    expected_toy_d_array_grad = temp2.transpose()(Eigen::all, arr);
    std::cout << "expected_toy_d_array_grad" << std::endl << expected_toy_d_array_grad.format(fullPrecisionFmt) << std::endl;
    EXPECT_MATRIX_NEAR(toy_d_array_grad, expected_toy_d_array_grad);
}

TEST(KEnRefTestSuite, TestRArrayToDArray2){


}
