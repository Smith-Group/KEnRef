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

TEST(KEnRefTestSuite, restOfTestsToWrite){
    std::vector<std::vector<std::vector<int>>> toy_grouping_list {{{0, 1, 2, 3}}, {{0, 1}, {2, 3}}, {{0}, {1}, {2}, {3}}};
    Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, 3> toy_r_mat(4, 3);
    toy_r_mat <<
              0.848351683690084, -0.529433112659379, 0,
            0.966177888683851, 0.257876496444355, 0,
            0.966177888683851, -0.257876496444355, 0,
            0.848351683690084, 0.529433112659379, 0;

    std::cout << "toy_r_mat" << std::endl << toy_r_mat << std::endl;

    auto [toy_d_array, toy_d_array_grad] = KEnRef<KEnRef_Real_t>::r_array_to_d_array(toy_r_mat, true);

    std::vector<Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, 5>> toy_d_array_vec;
    toy_d_array_vec.reserve((toy_d_array.rows()));
//	std::cout << "toy_d_array_vec size\t" << toy_d_array_vec.size() << std::endl;
    for(int i = 0; i < toy_d_array.rows(); i++){
        toy_d_array_vec.emplace_back(toy_d_array.row(i));
    }

    for (int gg = 0; gg < toy_grouping_list.size(); gg++) {
        std::cout << "Calculating g" << gg+1 << std::endl;
        auto [toy_g_array, toy_g_array_grad] = KEnRef<KEnRef_Real_t>::d_array_to_g(toy_d_array_vec, toy_grouping_list[gg], true);
        std::cout << "toy_g_array" << std::endl << toy_g_array << std::endl;
        std::cout << "toy_g_array_grad" << std::endl;
        for(const auto& matrix: toy_g_array_grad){
            std::cout << matrix << /*std::endl <<*/ std::endl;
        }
        std::cout << "----------" << std::endl;
    }

    /////////////////////////////////////////////////////////////////////////

    Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, 3> model1(5, 3);
    model1 <<
           32.708, 53.484, 20.701,
            32.284, 52.123, 22.636,
            31.277, 51.654, 21.284,
            31.852, 49.646, 22.312,
            32.854, 49.716, 20.812;
    Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, 3> model2(5, 3);
    model2 <<
           32.733, 52.960, 22.152,
            33.130, 51.220, 23.736,
            31.878, 50.694, 22.613,
            33.471, 49.251, 21.415,
            34.819, 49.854, 22.481;
    std::vector<Eigen::Matrix<KEnRef_Real_t , Eigen::Dynamic, 3>> eros3_sub_coord{model1, model2};
    std::vector<std::tuple<int, int>> eros3_sub_atom_idPairs{{0, 1}, {0, 2}, {0, 3}, {0, 4}};
    std::vector<std::vector<std::vector<int>>> eros3_grouping_list {{{0}, {1}}, {{0, 1}}};

    auto eros3_sub_r_array = KEnRef<KEnRef_Real_t>::coord_array_to_r_array(eros3_sub_coord, eros3_sub_atom_idPairs);
    std::cout << "eros3_sub_r_array" << std::endl;
    for (int i = 0; i < eros3_sub_r_array.size(); i++) {
        auto r_array = eros3_sub_r_array[i];
        std::cout << "Model " << i+1 << std::endl << r_array << std::endl;
    }
    auto [eros3_sub_d_array, eros3_sub_d_array_grad] = KEnRef<KEnRef_Real_t>::r_array_to_d_array(eros3_sub_r_array, true);
    std::cout << "eros3_sub_d_array" << std::endl;
    for(int i = 0; i < eros3_sub_d_array.size(); i++){
        auto matrix = eros3_sub_d_array[i];
        std::cout << "eros3_sub_d_array " << i+1 << std::endl << matrix << std::endl;
    }
    std::cout << "eros3_sub_d_array_grad" << std::endl;
    for(int i = 0; i < eros3_sub_d_array_grad.size(); i++){
        auto matrix = eros3_sub_d_array_grad[i];
        std::cout << "eros3_sub_d_array_grad " << i+1 << std::endl << matrix << std::endl;
    }

    std::vector<Eigen::VectorX<KEnRef_Real_t>> g_list;
    for (int gg = 0; gg < eros3_grouping_list.size(); gg++) {
        std::cout << "Calculating eros3_sub g" << gg+1 << std::endl;
        auto [eros3_sub_g_list, eros3_sub_g_list_grad] = KEnRef<KEnRef_Real_t>::d_array_to_g(eros3_sub_d_array, eros3_grouping_list[gg], true); //TODO d_arrays_to_g ???
        std::cout << "eros3_sub_g_list" << std::endl << eros3_sub_g_list << std::endl;
        std::cout << "eros3_sub_g_list_grad" << std::endl;
        for(const auto& matrix: eros3_sub_g_list_grad){
            std::cout << matrix << /*std::endl <<*/ std::endl;
            std::cout << "----------" << std::endl;
        }
        std::cout << "=============" << std::endl;
        g_list.emplace_back(eros3_sub_g_list);
    }

    auto g_matrix = KEnRef<KEnRef_Real_t>::vectorOfVectors_to_Matrix(g_list);
    std::cout << "g_matrix" << std::endl << g_matrix << std::endl;


    auto eros3_sub_g = KEnRef<KEnRef_Real_t>::coord_array_to_g(eros3_sub_coord, eros3_sub_atom_idPairs, eros3_grouping_list);
    std::cout << "eros3_sub_g" << std::endl << eros3_sub_g << std::endl;

    std::vector<Eigen::MatrixX3<KEnRef_Real_t>> m1_twice = {eros3_sub_coord[0], eros3_sub_coord[0]};
    //get g values using M1 twice
    auto eros3_sub_1_g = KEnRef<KEnRef_Real_t>::coord_array_to_g(m1_twice, eros3_sub_atom_idPairs, eros3_grouping_list);
    std::cout << "eros3_sub_1_g" << std::endl << eros3_sub_1_g << std::endl;


    auto[eros3_sub_energy, eros3_sub_energy_grad] = KEnRef<KEnRef_Real_t>::coord_array_to_energy(eros3_sub_coord, eros3_sub_atom_idPairs, eros3_grouping_list, eros3_sub_1_g, 1.0, true);
    std::cout << "eros3_sub_energy" << std::endl << eros3_sub_energy << std::endl;
    for (const auto& mat: eros3_sub_energy_grad) {
        std::cout << "eros3_sub_energy_grad" << std::endl << mat << std::endl;
    }


    //	auto [g_list, eros3_sub_g_list_grad] = KEnRef::d_array_to_g(eros3_sub_d_array, eros3_grouping_list, true);
    //	KEnRef::g_to_energy(g_matrix, eros3_sub_1_g, 1.0, true);

    //	exit(0);


}
