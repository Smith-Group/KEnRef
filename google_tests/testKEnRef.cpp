#include "gtest/gtest.h"
#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <fstream>
#include "core/KEnRef.h"
#include "core/IoUtils.h"
#include "testHelper.h"

Eigen::IOFormat fullPrecisionFmt(Eigen::FullPrecision);

TEST(KEnRefTestSuite, TestRArrayToDArray1) {

    std::vector<std::vector<std::vector<int>>> toy_grouping_list{{{0, 1, 2, 3}},
                                                                 {{0, 1}, {2, 3}},
                                                                 {{0}, {1}, {2}, {3}}};
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

    const auto& [toy_d_array, toy_d_array_grad] = KEnRef<KEnRef_Real_t>::r_array_to_d_array(toy_r_mat, true);
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
    TestHelper<KEnRef_Real_t>::EXPECT_MATRIX_NEAR(toy_d_array, expected_toy_d_array);

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
    expected_toy_d_array_grad = temp2.transpose()(Eigen::indexing::all, arr);
    std::cout << "expected_toy_d_array_grad" << std::endl << expected_toy_d_array_grad.format(fullPrecisionFmt) << std::endl;
    TestHelper<KEnRef_Real_t>::EXPECT_MATRIX_NEAR(toy_d_array_grad, expected_toy_d_array_grad);
}

TEST(KEnRefTestSuite, testDArrayToG){
    std::vector<std::vector<std::vector<int>>> toy_grouping_list{{{0, 1, 2, 3}},
                                                                 {{0, 1}, {2, 3}},
                                                                 {{0}, {1}, {2}, {3}}};

    Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, 5> toy_d_array(4, 5);
    Eigen::Matrix<KEnRef_Real_t, 5, Eigen::Dynamic> temp(5, 4);
    temp <<
         -0.500000000000001, -0.5, -0.5, -0.500000000000001,
            0.380532565661007, 0.750843527257811, 0.750843527257811, 0.380532565661007,
            0, 0, 0, 0,
            0, 0, 0, 0,
            -0.777942778404335, 0.431548372230798, -0.431548372230798, 0.777942778404335;
    toy_d_array = temp.transpose();

//    auto [toy_d_array, toy_d_array_grad] = KEnRef<KEnRef_Real_t>::r_array_to_d_array(toy_r_mat, true);

    std::vector<Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, 5>> toy_d_array_vec;
    toy_d_array_vec.reserve((toy_d_array.rows()));
//	std::cout << "toy_d_array_vec size\t" << toy_d_array_vec.size() << std::endl;
    for(int i = 0; i < toy_d_array.rows(); i++){
        toy_d_array_vec.emplace_back(toy_d_array.row(i));
    }

    for (int gg = 0; gg < toy_grouping_list.size(); gg++) {
        std::cout << "Calculating g" << gg+1 << std::endl;
        auto [toy_g_array, toy_g_array_grad] = KEnRef<KEnRef_Real_t>::d_array_to_g(toy_d_array_vec,
                                                                                   toy_grouping_list[gg], true,
                                                                                   0);
        std::cout << "toy_g_array" << std::endl << toy_g_array << std::endl;
        std::cout << "toy_g_array_grad" << std::endl;
        for(const auto& matrix: toy_g_array_grad){
            std::cout << matrix << /*std::endl <<*/ std::endl;
        }
        std::cout << "----------" << std::endl;
    }
    //TODO calculate and validate
}

TEST(KEnRefTestSuite, testSaturate) {
    CoordsMatrixType<KEnRef_Real_t> testMatrix(1, 3);
    CoordsMatrixType<KEnRef_Real_t> expectedMatrix(1, 3);

    testMatrix << 0., 0., 0.;
    expectedMatrix << 0., 0., 0.;
    KEnRef<KEnRef_Real_t>::saturate(testMatrix, 0.0, 0);
    EXPECT_EQ(testMatrix, expectedMatrix);

    testMatrix << 1., 1., 1.;
    expectedMatrix << 1., 1., 1.;
    KEnRef<KEnRef_Real_t>::saturate(testMatrix, 3.0, 0);
    EXPECT_EQ(testMatrix, expectedMatrix);

    testMatrix << 0., 0., 900.;
    expectedMatrix << 0., 0., 900.;
    KEnRef<KEnRef_Real_t>::saturate(testMatrix, 1000. * 1000., 0);
    EXPECT_EQ(testMatrix, expectedMatrix);

    testMatrix << 1000., 0., 0.;
    expectedMatrix << 1000., 0., 0.;
    KEnRef<KEnRef_Real_t>::saturate(testMatrix, 1000.0 * 1000.0, 0);
    TestHelper<KEnRef_Real_t>::EXPECT_MATRIX_NEAR(testMatrix, expectedMatrix);

    testMatrix << 0., 10000., 0.;
    expectedMatrix << 0., 1000., 0.;
    KEnRef<KEnRef_Real_t>::saturate(testMatrix, 1000.0 * 1000.0, 0);
    TestHelper<KEnRef_Real_t>::EXPECT_MATRIX_NEAR(testMatrix, expectedMatrix, .0001);
}

TEST(KEnRefTestSuite, TestCoordArrayToEnergyFiniteDifferenceMethodTest){
    auto coordsArray_base = CoordsMatrixType<double>(409, 3);

    std::ifstream coordsFileStream("../../res/google_tests/coords.txt");
    auto experimentalData_table = IoUtils::readTable(
            "../../res/10nsstart+fitting/singleton_data_10nsstart+fit_3-5_1977pairs_80_A.csv", true);

    auto tempCoordsTable = IoUtils::read_uniform_table_of<double>(coordsFileStream);
    for (int i = 0; i < tempCoordsTable.size(); ++i) {
        coordsArray_base.row(i) = Eigen::RowVector3<double>{tempCoordsTable[i][0], tempCoordsTable[i][1], tempCoordsTable[i][2]};
    }

    const auto& atomIdPairs = IoUtils::readAtomIdPairs("../../res/google_tests/atomIdPairs.txt");
    auto atomIdPairsMatrix = Eigen::Matrix<int, Eigen::Dynamic, 2>(atomIdPairs.size(), 2);
    for (int i = 0; i < atomIdPairs.size(); ++i) {
        auto [lt, rt] = atomIdPairs[i];
        atomIdPairsMatrix(i, 0) = lt;
        atomIdPairsMatrix(i, 1) = rt;
    }
    std::cout << atomIdPairsMatrix.transpose() << std::endl;

    std::vector<std::vector<std::string> > data = std::get<1>(experimentalData_table);
    auto g0 = Eigen::Matrix<double, Eigen::Dynamic, 2>(data.size(), 2);
    for (int i = 0; i < data.size(); ++i) {
        const auto& record = data[i];
        std::istringstream temp1(record[5]), temp2(record[6]);
        temp1 >> g0(i, 0);
        temp2 >> g0(i, 1);
    }

    double k = 5e8;
    double n = 0.25;
    auto simulated_grouping_list = std::vector<std::vector<std::vector<int>>>{{{0}}, {{0}}};
    auto allSimulationsSubAtomsX_vector = std::vector<CoordsMatrixType<double>>{coordsArray_base};
    const auto &[energy_base, allDerivatives_vector_base] =
            KEnRef<double>::coord_array_to_energy(allSimulationsSubAtomsX_vector, atomIdPairs,
                                                         simulated_grouping_list, g0,
                                                         k, n, true, 20);

//    std::cout << "Energy " << energy_base << "\n first 20 lines of derivatives in TestKEnRef\n"
//              << allDerivatives_vector_base[0].topRows(20) << std::endl;

    EXPECT_NEAR(energy_base, 35.9427, 1e-2);

    double  delta = 1e-6;
    Eigen::IOFormat heavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
    for (int model = 0; model < allSimulationsSubAtomsX_vector.size(); ++model) {
        for (int i = 0; i < coordsArray_base.rows(); ++i) {
            for (int j = 0; j < 3; ++j) {
                CoordsMatrixType<double> coordsArray_derived(coordsArray_base);
                coordsArray_derived(i, j) += delta;
                allSimulationsSubAtomsX_vector = std::vector<CoordsMatrixType<double>>{coordsArray_derived};
                const auto &[energy_derived, allDerivatives_vector_derived] =
                        KEnRef<double>::coord_array_to_energy(allSimulationsSubAtomsX_vector, atomIdPairs,
                                                                     simulated_grouping_list, g0,
                                                                     k, n, true, 20);

                double E_delta = energy_derived - energy_base;
                auto diff_mat = allDerivatives_vector_derived[model] - allDerivatives_vector_base[model];
                std::cout << "changed (" << i << ", " << j << ")";
                std::cout << std::scientific;

                std::cout << "\tE_delta = " << E_delta /*<< std::endl*/;
                std::cout << "\tF_delta = " << diff_mat(i,j) /*<< std::endl*/;

                double d_FD = E_delta /delta;
                double d_Anal = allDerivatives_vector_base[model](i,j);
                double method_diff = abs((d_Anal - d_FD) / d_FD);
//                double method_diff = abs((d_Anal - d_FD) / d_Anal);

                std::cout << "\tdFD = " << d_FD;
                std::cout << "\td_Anal = " << d_Anal;
                std::cout << "\t(d_Anal - d_FD) / d_FD) = " << method_diff;

//                Eigen::Index maxRow, maxCol;
//                double max = diff_mat.cwiseAbs().maxCoeff(&maxRow, &maxCol);
//                std::cout << "\tMax F_delta found at ("<< maxRow << ", " << maxCol << "). Value = " << std::scientific << max;
//                std::cout << "\t(dE-df)/d " << ((E_delta - diff_mat(i,j))/delta);

                std::cout << std::endl;
//                std::cout << "\ndiff_mat=\n" << diff_mat.format(heavyFmt) << std::endl;
            }
        }
    }
}

TEST(KEnRefTestSuite, testRestOfTestsToWrite){
    CoordsMatrixType<KEnRef_Real_t> model1(5, 3);
    model1 <<
           32.708, 53.484, 20.701,
            32.284, 52.123, 22.636,
            31.277, 51.654, 21.284,
            31.852, 49.646, 22.312,
            32.854, 49.716, 20.812;
    CoordsMatrixType<KEnRef_Real_t> model2(5, 3);
    model2 <<
           32.733, 52.960, 22.152,
            33.130, 51.220, 23.736,
            31.878, 50.694, 22.613,
            33.471, 49.251, 21.415,
            34.819, 49.854, 22.481;
    std::vector<CoordsMatrixType<KEnRef_Real_t>> eros3_sub_coord{model1, model2};
    std::vector<std::tuple<int, int>> eros3_sub_atom_idPairs{{0, 1}, {0, 2}, {0, 3}, {0, 4}};
    std::vector<std::vector<std::vector<int>>> eros3_grouping_list {{{0}, {1}}, {{0, 1}}};

    auto eros3_sub_r_array = KEnRef<KEnRef_Real_t>::coord_array_to_r_array(eros3_sub_coord, eros3_sub_atom_idPairs);
    std::cout << "eros3_sub_r_array" << std::endl;
    for (int i = 0; i < eros3_sub_r_array.size(); i++) {
        const auto &r_array = eros3_sub_r_array[i];
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
        auto [eros3_sub_g_list, eros3_sub_g_list_grad] = KEnRef<KEnRef_Real_t>::d_array_to_g(eros3_sub_d_array,
                                                                                             eros3_grouping_list[gg],
                                                                                             true, 0); //TODO d_array_to_g_multiple_groupings ???
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

    std::vector<CoordsMatrixType<KEnRef_Real_t>> m1_twice = {eros3_sub_coord[0], eros3_sub_coord[0]};
    //get g values using M1 twice
    auto eros3_sub_1_g = KEnRef<KEnRef_Real_t>::coord_array_to_g(m1_twice, eros3_sub_atom_idPairs, eros3_grouping_list);
    std::cout << "eros3_sub_1_g" << std::endl << eros3_sub_1_g << std::endl;


    auto [eros3_sub_energy, eros3_sub_energy_grad] =
            KEnRef<KEnRef_Real_t>::coord_array_to_energy(eros3_sub_coord, eros3_sub_atom_idPairs, eros3_grouping_list,
                                                         eros3_sub_1_g, 1.0, 0.25, true);
    std::cout << "eros3_sub_energy" << std::endl << eros3_sub_energy << std::endl;
    for (const auto& mat: eros3_sub_energy_grad) {
        std::cout << "eros3_sub_energy_grad" << std::endl << mat << std::endl;
    }
    //	auto [g_list, eros3_sub_g_list_grad] = KEnRef::d_array_to_g(eros3_sub_d_array, eros3_grouping_list, true);
    //	KEnRef::g_to_energy(g_matrix, eros3_sub_1_g, 1.0, true);

    //	exit(0);
}

TEST(KEnRefTestSuite, TestS2OrderParameters){
    std::vector<std::string> files {"../../res/google_tests/ensemble_coord_moldel1.csv", "../../res/google_tests/ensemble_coord_moldel2.csv"};
//    std::vector<std::ifstream> ifStreams{};
//    ifStreams.emplace_back("../../res/google_tests/ensemble_coord_moldel1.csv");
//    ifStreams.emplace_back("../../res/google_tests/ensemble_coord_moldel2.csv");
    std::vector<CoordsMatrixType<KEnRef_Real_t>> coordsVector;
    coordsVector.reserve(files.size());
//    float f;
    for (int i = 0; i < files.size(); ++i) {
        auto tempCoordsData = std::get<1>(IoUtils::readTable(files[i], false, ","));
        coordsVector.emplace_back(tempCoordsData.size(), 3);
        for (int j = 0; j < tempCoordsData.size(); ++j) {
            for (int k = 0; k < 3; ++k) {
                std::istringstream temp(tempCoordsData[j][k]);
                temp >> coordsVector[i](j, k);
            }
//            coordsVector[i].row(j) = Eigen::RowVector3<double>{tempCoordsTable[j][0], tempCoordsTable[j][1], tempCoordsTable[j][2]};
        }
    }
    auto experimentalData_table = IoUtils::readTable(
            "../../res/google_tests/singleton_data_10nsstart+fit_0+10.csv", true, ",");

    std::vector<std::vector<std::string> > data = std::get<1>(experimentalData_table);
    std::vector<std::tuple<int, int>> atomIdPairs;
    auto expectedS2 = Eigen::VectorX<KEnRef_Real_t>(data.size());
    int i1, i2;
    KEnRef_Real_t f;
    for (int i = 0; i < data.size(); ++i) {
        const auto& record = data[i];
        std::istringstream temp1(record[3]), temp2(record[4]), temp3(record[7]);
        temp1 >> i1;
        temp2 >> i2;
        atomIdPairs.emplace_back(i1 - 1, i2 - 1);
        temp3 >> expectedS2(i);
    }

    KEnRef_Real_t epsilon;
    if constexpr (std::is_same_v<KEnRef_Real_t, float>){
        epsilon = 5e-6;
    }else{
        epsilon = 1e-13;
    }
    std::cout << "testing float S2 difference with epsilon = " << epsilon << "\n";
    const auto &experimentalS2Double = KEnRef<KEnRef_Real_t>::s2OrderParams(coordsVector, atomIdPairs,0);
    std::cout << "expectedS2    \t" << expectedS2.transpose() << std::endl;
    std::cout << "experimentalS2\t" << experimentalS2Double.transpose() << std::endl;
    TestHelper<KEnRef_Real_t>::EXPECT_MATRIX_NEAR(expectedS2, experimentalS2Double, epsilon);

}
