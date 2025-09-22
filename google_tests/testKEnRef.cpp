#include "gtest/gtest.h"
#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <fstream>
#include <vector>

#include "core/KEnRef.h"
#include "core/IoUtils.h"
#include "testHelper.h"

Eigen::IOFormat fullPrecisionFmt(Eigen::FullPrecision);

template <typename KEnRef_Real>
std::vector<SpecDenData<KEnRef_Real>> getSpecDenData(const bool handleNames) {
        std::vector<SpecDenData<double> > dipole_kinetic_data{};
    dipole_kinetic_data.emplace_back(SpecDenData<double>{
        {
            {
                IoUtils::normalizeName(" HA  MET A   1 ", handleNames),
                IoUtils::normalizeName(" H1  MET A   1 ", handleNames)
            }
        },
        {},
        {
            {{0, 1, 2}},
            {{0}, {1}, {2}}
        },
        NamedMatrix<double>{(Eigen::Matrix<double, 2, 2>() <<
         1, -1,
         0, 1
        ).finished(), {{"0", "kens"}}},
        Table{
            {{"0", "1"}},
            {{"0","kens"}},
            {{"kens"}}
        }
    });
    dipole_kinetic_data.emplace_back(SpecDenData<double>{
        {
            {
                IoUtils::normalizeName(" H1  MET A   1 ", handleNames),
                IoUtils::normalizeName(" HD1 TYR A   3 ", handleNames),
            },
            {
                IoUtils::normalizeName(" H1  MET A   1 ", handleNames),
                IoUtils::normalizeName(" HD2 TYR A   3 ", handleNames),
            },
        },
        {},
        {
            {{0, 1, 2, 3, 4, 5,},},
            {{0, 2, 4}, {1, 3, 5},},
            {{0, 1}, {2, 3}, {4, 5},},
            {{0}, {1}, {2}, {3}, {4}, {5},},
        },
        NamedMatrix<double>{(Eigen::Matrix<double, 4, 4>() <<
         1, -1, -1, 1,
         0, 1, 0, -1,
         0, 0, 1, -1,
         0, 0, 0, 1
        ).finished(), {{"0","karo","kens","kens+karo"}}},
        Table{
            {
                {"0", "0", "1", "1"},
                {"0", "1", "0", "1"}
            },
            {{"0", "karo", "kens", "kens+karo"}},
            {{"kens", "karo"}}
        }
    });
    dipole_kinetic_data.emplace_back(SpecDenData<double>{
        {
            {
                IoUtils::normalizeName(" H1  MET A   1 ", handleNames),
                IoUtils::normalizeName(" HE1 MET A   1 ", handleNames),
            },
            {
                IoUtils::normalizeName(" H1  MET A   1 ", handleNames),
                IoUtils::normalizeName(" HE2 MET A   1 ", handleNames),
            },
            {
                IoUtils::normalizeName(" H1  MET A   1 ", handleNames),
                IoUtils::normalizeName(" HE3 MET A   1 ", handleNames),
            },
        },
        {},
        {
            {{0, 1, 2, 3, 4, 5, 6, 7, 8}},
            {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}},
            {{0, 3, 6}, {1, 4, 7}, {2, 5, 8}},
            {{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}},
        },
        NamedMatrix<double>{(Eigen::Matrix<double, 4, 4>() <<
         1, -1, -1, 1,
         0, 1, 0, -1,
         0, 0, 1, -1,
         0, 0, 0, 1
        ).finished(), {{"0","kens","kmethyl","kens+kmethyl"}}},
        Table{
            {
                {"0", "1", "0", "1"},
                {"0", "0", "1", "1"}
            },
            {{"0", "kens", "kmethyl", "kens+kmethyl"}},
            {{"kens", "kmethyl"}}
        }
    });
    dipole_kinetic_data.emplace_back(SpecDenData<double>{
        {
            {
                IoUtils::normalizeName(" HD1 PHE A  30 ", handleNames),
                IoUtils::normalizeName(" HD1 TYR A   3 ", handleNames),
            },
            {
                IoUtils::normalizeName(" HD2 PHE A  30 ", handleNames),
                IoUtils::normalizeName(" HD1 TYR A   3 ", handleNames),
            },
            {
                IoUtils::normalizeName(" HD1 PHE A  30 ", handleNames),
                IoUtils::normalizeName(" HD2 TYR A   3 ", handleNames),
            },
            {
                IoUtils::normalizeName(" HD2 PHE A  30 ", handleNames),
                IoUtils::normalizeName(" HD2 TYR A   3 ", handleNames),
            },
        },
        {},
        {
            {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},},
            {{0, 1, 4, 5, 8, 9}, {2, 3, 6, 7, 10, 11},},
            {{0, 2, 4, 6, 8, 10}, {1, 3, 5, 7, 9, 11},},
            {{0, 4, 8}, {1, 5, 9}, {2, 6, 10}, {3, 7, 11}},
            {{0, 1, 2, 3}, {4, 5, 6, 7}, {8, 9, 10, 11}},
            {{0, 1}, {4, 5}, {8, 9}, {2, 3}, {6, 7}, {10, 11}},
            {{0, 2}, {4, 6}, {8, 10}, {1, 3}, {5, 7}, {9, 11}},
            {{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}},
        },
        NamedMatrix<double>{(Eigen::Matrix<double, 8, 6>() <<
         1, -2, 1, -1, 2, -1,
         0, 1, -1, 0, -1, 1,
         0, 1, -1, 0, -1, 1,
         0, 0, 1, 0, 0, -1,
         0, 0, 0, 1, -2, 1,
         0, 0, 0, 0, 1, -1,
         0, 0, 0, 0, 1, -1,
         0, 0, 0, 0, 0, 1
        ).finished(), {{"0","karo","2karo","kens","kens+karo","kens+2karo"}}},
        Table{
            {
                {"0","0","0","1","1","1"},
                {"0","1","2","0","1","2"}
            },
            {{"0","karo","2karo","kens","kens+karo","kens+2karo"}},
            {{"kens","karo"}}
        }
    });
    dipole_kinetic_data.emplace_back(SpecDenData<double>{
        {
            {
                IoUtils::normalizeName(" HD1 TYR A   3 ", handleNames),
                IoUtils::normalizeName(" HE1 MET A   1 ", handleNames),
            },
            {
                IoUtils::normalizeName(" HD2 TYR A   3 ", handleNames),
                IoUtils::normalizeName(" HE1 MET A   1 ", handleNames),
            },
            {
                IoUtils::normalizeName(" HD1 TYR A   3 ", handleNames),
                IoUtils::normalizeName(" HE2 MET A   1 ", handleNames),
            },
            {
                IoUtils::normalizeName(" HD2 TYR A   3 ", handleNames),
                IoUtils::normalizeName(" HE2 MET A   1 ", handleNames),
            },
            {
                IoUtils::normalizeName(" HD1 TYR A   3 ", handleNames),
                IoUtils::normalizeName(" HE3 MET A   1 ", handleNames),
            },
            {
                IoUtils::normalizeName(" HD2 TYR A   3 ", handleNames),
                IoUtils::normalizeName(" HE3 MET A   1 ", handleNames),
            },
        },
        {},
        {
            {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17},},
            {{0, 2, 4, 6, 8, 10, 12, 14, 16}, {1, 3, 5, 7, 9, 11, 13, 15, 17},},
            {{0, 1, 2, 3, 4, 5}, {6, 7, 8, 9, 10, 11}, {12, 13, 14, 15, 16, 17},},
            {{0, 2, 4}, {6, 8, 10}, {12, 14, 16}, {1, 3, 5}, {7, 9, 11}, {13, 15, 17},},
            {{0, 6, 12}, {2, 8, 14}, {4, 10, 16}, {1, 7, 13}, {3, 9, 15}, {5, 11, 17},},
            {{0}, {2}, {4}, {6}, {8}, {10}, {12}, {14}, {16}, {1}, {3}, {5}, {7}, {9}, {11}, {13}, {15}, {17},},
        },
        NamedMatrix<double>{(Eigen::Matrix<double, 6, 6>() <<
         1, -1, -1, 1, 0, 0,
         0,  1, 0, -1, -1, 1,
         0,  0, 1, -1, 0, 0,
         0,  0, 0, 1, 0, -1,
         0,  0, 0, 0, 1, -1,
         0,  0, 0, 0, 0, 1
        ).finished(), {{"0","karo","kens","kens+karo","kmethyl","kens+kmethyl"}}},
        Table{
            {
                {"0","0","1","1","0","1"},
                {"0","0","0","0","1","1"},
                {"0","1","0","1","0","0"},
            },
            std::optional<std::vector<std::string>>{{"0","karo","kens","kens+karo","kmethyl","kens+kmethyl"}},
            std::optional<std::vector<std::string>>{{"kens","kmethyl", "karo"}}
        }
    });
    dipole_kinetic_data.emplace_back(SpecDenData<double>{
        {
            {
                IoUtils::normalizeName(" HZ1 LYS A   4 ", handleNames),
                IoUtils::normalizeName(" HE1 MET A   1 ", handleNames),
            },
            {
                IoUtils::normalizeName(" HZ2 LYS A   4 ", handleNames),
                IoUtils::normalizeName(" HE1 MET A   1 ", handleNames),
            },
            {
                IoUtils::normalizeName(" HZ3 LYS A   4 ", handleNames),
                IoUtils::normalizeName(" HE1 MET A   1 ", handleNames),
            },
            {
                IoUtils::normalizeName(" HZ1 LYS A   4 ", handleNames),
                IoUtils::normalizeName(" HE2 MET A   1 ", handleNames),
            },
            {
                IoUtils::normalizeName(" HZ2 LYS A   4 ", handleNames),
                IoUtils::normalizeName(" HE2 MET A   1 ", handleNames),
            },
            {
                IoUtils::normalizeName(" HZ3 LYS A   4 ", handleNames),
                IoUtils::normalizeName(" HE2 MET A   1 ", handleNames),
            },
            {
                IoUtils::normalizeName(" HZ1 LYS A   4 ", handleNames),
                IoUtils::normalizeName(" HE3 MET A   1 ", handleNames),
            },
            {
                IoUtils::normalizeName(" HZ2 LYS A   4 ", handleNames),
                IoUtils::normalizeName(" HE3 MET A   1 ", handleNames),
            },
            {
                IoUtils::normalizeName(" HZ3 LYS A   4 ", handleNames),
                IoUtils::normalizeName(" HE3 MET A   1 ", handleNames),
            },
        },
        {},
        {
            {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26},},
            {{0, 1, 2, 3, 4, 5, 6, 7, 8}, {9, 10, 11, 12, 13, 14, 15, 16, 17}, {18, 19, 20, 21, 22, 23, 24, 25, 26},},
            {{0, 1, 2, 9, 10, 11, 18, 19, 20}, {3, 4, 5, 12, 13, 14, 21, 22, 23}, {6, 7, 8, 15, 16, 17, 24, 25, 26},},
            {
                {0, 1, 2}, {3, 4, 5}, {6, 7, 8}, {9, 10, 11}, {12, 13, 14}, {15, 16, 17}, {18, 19, 20}, {21, 22, 23},
                {24, 25, 26}
            },
            {{0, 3, 6, 9, 12, 15, 18, 21, 24}, {1, 4, 7, 10, 13, 16, 19, 22, 25}, {2, 5, 8, 11, 14, 17, 20, 23, 26},},
            {
                {0, 3, 6}, {1, 4, 7}, {2, 5, 8}, {9, 12, 15}, {10, 13, 16}, {11, 14, 17}, {18, 21, 24}, {19, 22, 25},
                {20, 23, 26},
            },
            {
                {0, 9, 18}, {1, 10, 19}, {2, 11, 20}, {3, 12, 21}, {4, 13, 22}, {5, 14, 23}, {6, 15, 24}, {7, 16, 25},
                {8, 17, 26}
            },
            {
                {0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}, {12}, {13}, {14}, {15}, {16}, {17}, {18},
                {19}, {20}, {21}, {22}, {23}, {24}, {25}, {26}
            },
        },
        NamedMatrix<double>{(Eigen::Matrix<double, 8, 6>() <<
         1, -1, -2,  2,  1, -1,
         0,  1,  0, -2,  0,  1,
         0,  0,  1, -1, -1,  1,
         0,  0,  0,  1,  0, -1,
         0,  0,  1, -1, -1,  1,
         0,  0,  0,  1,  0, -1,
         0,  0,  0,  0,  1, -1,
         0,  0,  0,  0,  0,  1
        ).finished(), {{"0","kens","kmethyl","kens+kmethyl","2kmethyl","kens+2kmethyl"}}},
        Table{
            {
                {"0","1","0","1","0","1"},
                {"0","0","1","1","2","2"}
            },
            {{"0","kens","kmethyl","kens+kmethyl","2kmethyl","kens+2kmethyl"}},
            {{"kens","kmethyl"}}
        }
    });
    return dipole_kinetic_data;
}

template <typename KEnRef_Real>
std::map<std::string, int> get_atomNameMapping(const std::map<std::string, int>& atomNameMapping_to1, bool handleNames) {
    std::map<std::string, int> atomNameMapping = {};
    int maxId0 = -1;
    for (const auto &[atomName, id1]: atomNameMapping_to1) {
        std::string tempName = atomName;
        std::string atomNameNormalized = IoUtils::normalizeName(tempName, handleNames);
        const int id0 = id1 - 1;
        if (id0 > maxId0) {
            maxId0 = id0;
        }
        atomNameMapping[atomNameNormalized] = id0;
    }
    return atomNameMapping;
}

template <typename KEnRef_Real>
std::vector<CoordsMatrixType<KEnRef_Real>> getAllModels_allAtomCoordsMatrix(const std::string& FILENAME, std::tuple<int, int> boundsInclusive) {

    auto[first, last] =  boundsInclusive;
    int numModels = static_cast<int>(last - first + 1);
    std::map<std::string, int> atomNameMapping_to1 = IoUtils::getAtomMappingFromPdb<std::string, int>( FILENAME, IoUtils::fill_atomId_to_index_Map);
    auto allAtomCoordsMap_raw_vector = IoUtils::getMultipleAtomMappingFromPdb<int, Eigen::RowVector3<double> >(
        FILENAME, IoUtils::fill_atomIndex1_to_coords_Map<double>);

    // allAtomCoordsMap_raw_vector.resize(numModels); //replaced with [first, last]
    // const bool handleNames = IoUtils::should_handleNames(atomNameMapping_to1);

    int maxId0 = -1;
    for (const auto &[atomName, id1]: atomNameMapping_to1) {
        std::string tempName = atomName;
        std::string atomNameNormalized = IoUtils::normalizeName(tempName, handleNames);
        const int id0 = id1 - 1;
        if (id0 > maxId0) {
            maxId0 = id0;
        }
    }
    //prepare coordsMatrix (0 based)
    std::vector<int> allModelAtomIdsVector = {};
    for (int i = 0; i <= maxId0; ++i) {
        allModelAtomIdsVector.push_back(i);
    }
    std::vector<CoordsMatrixType<double> > allModels_allAtomCoordsMap{};
    for (int i = 0; i < numModels; ++i) {
        auto &allAtomCoordsMap_raw = allAtomCoordsMap_raw_vector[i + first];
        const auto &allAtomCoordsMap = IoUtils::extractCoords(allModelAtomIdsVector, false, allAtomCoordsMap_raw, true);
        allModels_allAtomCoordsMap.emplace_back(allAtomCoordsMap);
    }
    return allModels_allAtomCoordsMap;
}

TEST(KEnRefTestSuite, TestRArrayToDArray1) {
    std::vector<std::vector<std::vector<int> > > toy_grouping_list{
        {{0, 1, 2, 3}},
        {{0, 1}, {2, 3}},
        {{0}, {1}, {2}, {3}}
    };
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

    const auto &[toy_d_array, toy_d_array_grad] = KEnRef<KEnRef_Real_t>::r_array_to_d_array(toy_r_mat, true);
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
    int arr[15], idx = 0;
    for (int d = 0; d < 5; d++)
        for (int xyz = 0; xyz < 3; xyz++)
            arr[idx++] = d + 5 * xyz;
    const Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, Eigen::Dynamic> expected_toy_d_array_grad =
        temp2.transpose()( Eigen::indexing::all, arr);
    std::cout << "expected_toy_d_array_grad" << std::endl << expected_toy_d_array_grad.format(fullPrecisionFmt) <<
            std::endl;
    TestHelper<KEnRef_Real_t>::EXPECT_MATRIX_NEAR(toy_d_array_grad, expected_toy_d_array_grad);
}

TEST(KEnRefTestSuite, testDArrayToG) {
    const std::vector<std::vector<std::vector<int>>> toy_grouping_list{
        {{0, 1, 2, 3}},
        {{0, 1}, {2, 3}},
        {{0}, {1}, {2}, {3}}
    };

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

    std::vector<Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, 5> > toy_d_array_vec;
    toy_d_array_vec.reserve(toy_d_array.rows());
    //	std::cout << "toy_d_array_vec size\t" << toy_d_array_vec.size() << std::endl;
    for (int i = 0; i < toy_d_array.rows(); i++) {
        toy_d_array_vec.emplace_back(toy_d_array.row(i));
    }

    for (int gg = 0; gg < toy_grouping_list.size(); gg++) {
        std::cout << "Calculating g" << gg + 1 << std::endl;
        auto [toy_g_array, toy_g_array_grad] = KEnRef<KEnRef_Real_t>::d_array_to_g(toy_d_array_vec,
            toy_grouping_list[gg],true, 0);
        std::cout << "toy_g_array" << std::endl << toy_g_array << std::endl;
        std::cout << "toy_g_array_grad" << std::endl;
        for (const auto &matrix: toy_g_array_grad) {
            std::cout << matrix << std::endl;
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

TEST(KEnRefTestSuite, TestCoordArrayToEnergyFiniteDifferenceMethodTest) {
    auto coordsArray_base = CoordsMatrixType<double>(409, 3);

    std::ifstream coordsFileStream("../../res/google_tests/coords.txt");
    auto experimentalData_table = IoUtils::readTable(
        "../../res/10nsstart+fitting/singleton_data_10nsstart+fit_3-5_1977pairs_80_A.csv", true);

    auto tempCoordsTable = IoUtils::read_uniform_table_of<double>(coordsFileStream);
    for (int i = 0; i < tempCoordsTable.size(); ++i) {
        coordsArray_base.row(i) = Eigen::RowVector3<double>{
            tempCoordsTable[i][0], tempCoordsTable[i][1], tempCoordsTable[i][2]
        };
    }

    const auto &atomIdPairs = IoUtils::readAtomIdPairs("../../res/google_tests/atomIdPairs.txt");
    auto atomIdPairsMatrix = Eigen::Matrix<int, Eigen::Dynamic, 2>(atomIdPairs.size(), 2);
    for (int i = 0; i < atomIdPairs.size(); ++i) {
        auto [lt, rt] = atomIdPairs[i];
        atomIdPairsMatrix(i, 0) = lt;
        atomIdPairsMatrix(i, 1) = rt;
    }
    std::cout << atomIdPairsMatrix.transpose() << std::endl;

    auto g0 = Eigen::Matrix<double, Eigen::Dynamic, 2>(experimentalData_table.rowCount(), 2);
    for (int i = 0; i < experimentalData_table.rowCount(); ++i) {
        std::istringstream temp1(experimentalData_table(i, "g1")), temp2(experimentalData_table(i, "g2"));
        temp1 >> g0(i, 0);
        temp2 >> g0(i, 1);
    }

    double k = 5e8;
    double n = 0.25;
    auto simulated_grouping_list = std::vector<std::vector<std::vector<int> > >{{{0}}, {{0}}};
    auto allSimulationsSubAtomsX_vector = std::vector<CoordsMatrixType<double> >{coordsArray_base};
    const auto &[energy_base, allDerivatives_vector_base] =
            KEnRef<double>::coord_array_to_g_energy(allSimulationsSubAtomsX_vector, atomIdPairs,
                                                    simulated_grouping_list, g0,
                                                    k, n, true, 20);

    //    std::cout << "Energy " << energy_base << "\n first 20 lines of derivatives in TestKEnRef\n"
    //              << allDerivatives_vector_base[0].topRows(20) << std::endl;

    EXPECT_NEAR(energy_base, 35.9427, 1e-2);

    Eigen::IOFormat heavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
    for (int model = 0; model < allSimulationsSubAtomsX_vector.size(); ++model) {
        for (int i = 0; i < coordsArray_base.rows(); ++i) {
            for (int j = 0; j < 3; ++j) {
                constexpr double delta = 1e-6;
                CoordsMatrixType<double> coordsArray_derived(coordsArray_base);
                coordsArray_derived(i, j) += delta;
                allSimulationsSubAtomsX_vector = std::vector<CoordsMatrixType<double> >{coordsArray_derived};
                const auto &[energy_derived, allDerivatives_vector_derived] =
                        KEnRef<double>::coord_array_to_g_energy(allSimulationsSubAtomsX_vector, atomIdPairs,
                                                                simulated_grouping_list, g0,
                                                                k, n, true, 20);

                double E_delta = energy_derived - energy_base;
                auto diff_mat = allDerivatives_vector_derived.value()[model] - allDerivatives_vector_base.value()[model];
                std::cout << "changed (" << i << ", " << j << ")";
                std::cout << std::scientific;

                std::cout << "\tE_delta = " << E_delta /*<< std::endl*/;
                std::cout << "\tF_delta = " << diff_mat(i, j) /*<< std::endl*/;

                double d_FD = E_delta / delta;
                double d_Anal = allDerivatives_vector_base.value()[model](i, j);
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

TEST(KEnRefTestSuite, testEROS3) {
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
    std::map<std::string, int> atomNameMapping0{
        {" HA  MET A   1 ", 0},
        {" HB2 MET A   1 ", 1},
        {" HB3 MET A   1 ", 2},
        {" HG2 MET A   1 ", 3},
        {" HG3 MET A   1 ", 4},
    };
    std::vector<std::tuple<std::string, std::string> > atomNamePairs = {
        {" HA  MET A   1 ", " HB2 MET A   1 ",},
        {" HA  MET A   1 ", " HB3 MET A   1 ",},
        {" HA  MET A   1 ", " HG2 MET A   1 ",},
        {" HA  MET A   1 ", " HG3 MET A   1 ",},
    };
    auto eros3_sub_atom_idPairs = *KEnRef<float>::atomNamePairs_2_atomIdPairs(atomNamePairs, atomNameMapping0);
    std::vector<std::tuple<int, int> > expectedEros3_sub_atom_idPairs{{0, 1}, {0, 2}, {0, 3}, {0, 4}};

    for (int i = 0; i < eros3_sub_atom_idPairs.size(); i++) {
        EXPECT_EQ(expectedEros3_sub_atom_idPairs[i], eros3_sub_atom_idPairs[i]);
    }


    std::vector<std::vector<std::vector<int> > > eros3_grouping_list{{{0}, {1}}, {{0, 1}}};

    auto eros3_sub_r_array = KEnRef<KEnRef_Real_t>::coord_array_to_r_array(eros3_sub_coord, eros3_sub_atom_idPairs);
    std::cout << "eros3_sub_r_array" << std::endl;
    for (int i = 0; i < eros3_sub_r_array.size(); i++) {
        const auto &r_array = eros3_sub_r_array[i];
        std::cout << "Model " << i + 1 << std::endl << r_array << std::endl;
    }
    auto [eros3_sub_d_array, eros3_sub_d_array_grad] = KEnRef<KEnRef_Real_t>::r_array_to_d_array(
        eros3_sub_r_array, true);
    std::cout << "eros3_sub_d_array" << std::endl;
    for (int i = 0; i < eros3_sub_d_array.size(); i++) {
        auto matrix = eros3_sub_d_array[i];
        std::cout << "eros3_sub_d_array " << i + 1 << std::endl << matrix << std::endl;
    }
    std::cout << "eros3_sub_d_array_grad" << std::endl;
    for (int i = 0; i < eros3_sub_d_array_grad.size(); i++) {
        auto matrix = eros3_sub_d_array_grad[i];
        std::cout << "eros3_sub_d_array_grad " << i + 1 << std::endl << matrix << std::endl;
    }

    auto [eros3_sub_g_list, eros3_sub_g_list_grad] = KEnRef<KEnRef_Real_t>::d_array_to_g_multiple_groupings(
        eros3_sub_d_array, eros3_grouping_list, true);
    std::vector<Eigen::VectorX<KEnRef_Real_t> > g_list;
    for (int gg = 0; gg < eros3_grouping_list.size(); gg++) {
        std::cout << "eros3_sub_g_list" << std::endl << eros3_sub_g_list[gg] << std::endl;
        std::cout << "eros3_sub_g_list_grad" << std::endl;
        for (const auto &matrix: eros3_sub_g_list_grad[gg]) {
            std::cout << matrix << /*std::endl <<*/ std::endl;
            std::cout << "----------" << std::endl;
        }
        std::cout << "=============" << std::endl;
    }

    auto eros3_sub_g =
            KEnRef<KEnRef_Real_t>::coord_array_to_g_matrix(eros3_sub_coord, eros3_sub_atom_idPairs,
                                                           eros3_grouping_list);
    std::cout << "eros3_sub_g" << std::endl << eros3_sub_g << std::endl;

    std::vector<CoordsMatrixType<KEnRef_Real_t> > m1_twice = {eros3_sub_coord[0], eros3_sub_coord[0]};
    //get g values using M1 twice
    auto eros3_sub_1_g =
            KEnRef<KEnRef_Real_t>::coord_array_to_g_matrix(m1_twice, eros3_sub_atom_idPairs, eros3_grouping_list);
    std::cout << "eros3_sub_1_g" << std::endl << eros3_sub_1_g << std::endl;


    auto [eros3_sub_energy, eros3_sub_energy_grad] =
            KEnRef<KEnRef_Real_t>::coord_array_to_g_energy(eros3_sub_coord, eros3_sub_atom_idPairs, eros3_grouping_list,
                                                           eros3_sub_1_g, 1.0, 0.25, true);
    std::cout << "eros3_sub_energy\n" << eros3_sub_energy << std::endl;
    for (const auto &mat: eros3_sub_energy_grad.value()) {
        std::cout << "eros3_sub_energy_grad" << std::endl << mat << std::endl;
    }

    auto [eros3_sub_energy_1, eros3_sub_energy_grad_1] =
            KEnRef<KEnRef_Real_t>::coord_array_to_g_energy_refactored(eros3_sub_coord, eros3_sub_atom_idPairs,
                                                                      eros3_grouping_list, eros3_sub_1_g, 1.0, 0.25,
                                                                      true);

    EXPECT_EQ(eros3_sub_energy, eros3_sub_energy_1);
    for (int i = 0; i < eros3_sub_energy_grad.value().size(); ++i) {
        TestHelper<KEnRef_Real_t>::EXPECT_MATRIX_NEAR(eros3_sub_energy_grad.value()[i], eros3_sub_energy_grad_1[i]);
    }
    //	auto [g_list, eros3_sub_g_list_grad] = KEnRef::d_array_to_g(eros3_sub_d_array, eros3_grouping_list, true);
    //	KEnRef::g_to_energy(g_matrix, eros3_sub_1_g, 1.0, true);

    //	exit(0);
}


TEST(KEnRefTestSuite, testGB3) {
    static const auto FILENAME = "../../res/google_tests/2lum.pdb";
    std::map<std::string, int> atomNameMapping_to1 = IoUtils::getAtomMappingFromPdb<std::string, int>(
        FILENAME, IoUtils::fill_atomId_to_index_Map);
    const std::vector<CoordsMatrixType<double> >& allModels_allAtomCoordsMap = getAllModels_allAtomCoordsMatrix<double>(FILENAME,{0, 2});
    int numModels = static_cast<int>(allModels_allAtomCoordsMap.size());

    auto allModels_allAtomCoordsMapMeters = std::vector<CoordsMatrixType<double> >(numModels);
    std::vector<std::vector<double> > expected_g_list_vector{
        {1.101178e+57, 2.440676e+57},
        {2.122924e+55, 2.952983e+55, 3.914569e+55, 6.73263e+55},
        {8.2869e+54, 3.089198e+55, 1.064407e+55, 3.56061e+55},
        {5.225709e+54, 5.667142e+54, 7.083363e+54, 7.737046e+54, 5.543656e+54, 7.648771e+54, 7.730481e+54, 1.103205e+55 },
        {5.121015e+53, 5.437515e+53, 6.015779e+53, 7.066374e+53, 5.541675e+53, 7.436357e+53},
        {6.060657e+52, 6.422185e+52, 6.089768e+52, 6.564926e+52, 6.097483e+52, 6.553809e+52, 6.127958e+52, 6.702502e+52 },
    };
    std::vector<double> expected_g_list_Epsilons_vector{1e52, 1e50, 1e49, 1e49, 1e48, 1e47};

    std::vector<Eigen::RowVectorX<double> > expected_a_matrix_vector{
        (Eigen::RowVectorX<double>(2)<< 1.101178e+57, 1.339498e+57).finished(),
        (Eigen::RowVectorX<double>(4)<< 2.122924e+55, 8.300591e+54, 1.791645e+55, 1.988002e+55).finished(),
        (Eigen::RowVectorX<double>(4)<< 8.2869e+54, 2.260508e+55, 2.357173e+54, 2.35694e+54).finished(),
        (Eigen::RowVectorX<double>(6)<< 5.225709e+54, 2.299086e+54, 2.122508e+53, 3.179466e+53, 1.992855e+54, 9.842026e+53).finished(),
        (Eigen::RowVectorX<double>(6)<< 5.121015e+53, 3.164998e+52, 8.947633e+52, 7.340958e+52, 1.041598e+52, 2.658234e+52).finished(),
        (Eigen::RowVectorX<double>(6)<< 6.060657e+52, 3.615275e+51, 6.59358e+50, 2.084292e+51, 1.364406e+49, 4.587508e+49).finished(),
    };
    std::vector<double> expected_a_matrix_Epsilons_vector{1e52, 1e49, 1e49, 1e48, 1e47, 1e46};

    std::vector<Eigen::RowVectorX<double> > expected_lambda_prime_vector{
        (Eigen::RowVectorX<double>(2)<< -2.5e+08, -7.5e+08).finished(),
        (Eigen::RowVectorX<double>(4)<< -250000000, -250010000, -750000000, -750010000).finished(),
        (Eigen::RowVectorX<double>(4)<< -2.50000e+08, -7.50000e+08, -1.00025e+12, -1.00075e+12).finished(),
        (Eigen::RowVectorX<double>(6)<< -250000000, -250010000, -250020000, -750000000, -750010000, -750020000).finished(),
        (Eigen::RowVectorX<double>(6)<< -2.50000e+08, -2.50010e+08, -7.50000e+08, -7.50010e+08, -1.00025e+12, -1.00075e+12).finished(),
        (Eigen::RowVectorX<double>(6)<< -2.50000e+08, -7.50000e+08, -1.00025e+12, -1.00075e+12, -2.00025e+12, -2.00075e+12).finished(),
    };
    std::vector<double> expected_lambda_prime_vect_Epsilons_vector{1e7, 1e-8, 1e3, 1e-8, 1e3, 1e3};

    std::vector<double> expected_sigma_vector{-0.3469721, -0.009440785, -0.003519733, -0.001993476, -0.0001351006, -1.400239e-05};
    std::vector<double> expected_sigma_vector_Epsilons_vector{1e-6, 1e-8, 1e-8, 1e-8, 1e-9, 1e-10};
    std::vector<NamedRowVector<double>> expected_sigma_grad{
        (Eigen::RowVectorX<double>(2) << -2.267168e-58, -7.265161e-59).finished(),
        (Eigen::RowVectorX<double>(4) << -2.267168e-58, -2.267076e-58, -7.265161e-59, -7.265055e-59).finished(),
        (Eigen::RowVectorX<double>(4) << -2.267168e-58, -7.265161e-59, 2.846774e-61, 2.845352e-61).finished(),
        (Eigen::RowVectorX<double>(6) << -2.267168e-58, -2.267076e-58, -2.266985e-58, -7.265161e-59, -7.265055e-59, -7.264949e-59).finished(),
        (Eigen::RowVectorX<double>(6) << -2.267168e-58, -2.267076e-58, -7.265161e-59, -7.265055e-59, 2.846774e-61, 2.845352e-61).finished(),
        (Eigen::RowVectorX<double>(6) << -2.267168e-58, -7.265161e-59, 2.846774e-61, 2.845352e-61, 1.423664e-61, 1.423308e-61).finished(),
    };
    std::vector<double> expected_sigma_grad_Epsilons_vector{1e-63, 1e-63, 1e-63, 1e-63, 1e-63, 1e-63};


    const auto& dipole_kinetic_data = getSpecDenData<double>(true); //handleNames);

    bool handleNames = IoUtils::should_handleNames(atomNameMapping_to1);
    const auto & atomNameMapping = get_atomNameMapping<double>(atomNameMapping_to1, handleNames);
    int maxId0 = -1;
    for (const auto &[atomName, id1]: atomNameMapping_to1) {
        std::string tempName = atomName;
        std::string atomNameNormalized = IoUtils::normalizeName(tempName, handleNames);
        int id0 = id1 - 1;
        if (id0 > maxId0) {
            maxId0 = id0;
        }
    }

    for (int interaction = 0; interaction < dipole_kinetic_data.size(); ++interaction) {
        std::cout << "======== new interaction ==========" << std::endl;
        auto atomIdPairs = IoUtils::atomNamePairs_2_atomIdPairs(dipole_kinetic_data[interaction].get_atom_pairs(), atomNameMapping);
        for (int i = 0; i < numModels; ++i)
            allModels_allAtomCoordsMapMeters[i] = allModels_allAtomCoordsMap[i] * 1e-10;
        auto r_arrays = KEnRef<double>::coord_array_to_r_array(allModels_allAtomCoordsMapMeters, atomIdPairs);
        // for (int i = 0; i < r_array.size(); ++i) {
        //     std::cout << "r_array[" << i << "] = " << std::endl;
        //     std::cout << r_array[i] << std::endl;
        // }

        const auto &[d_arrays, d_arrays_grad] = KEnRef<double>::r_array_to_d_array(r_arrays, true, false);
        // std::cout << "d_arrays\n" << d_arrays[0] << std::endl;
        const auto &grouping_matrix = IoUtils::subset_idx_to_grouping_mat(
            dipole_kinetic_data[interaction].get_multiple_grouping());
        uint nshift = grouping_matrix.cols() / numModels;
        const auto &d_arrays_shifted = KEnRef<double>::array_shift(d_arrays, nshift);
        // std::cout << "d_arrays:\n";
        // for (int i = 0; i < d_arrays.size(); ++i) {
        //     std::cout << "model "<<i << std::endl;
        //     std::cout << d_arrays[i] << std::endl;
        // }
        // std::cout << "d_arrays_shifted\n";
        // for (int i = 0; i < d_arrays_shifted.size(); ++i) {
        //     std::cout << "model - interactionComponent " << i << std::endl;
        //     std::cout << d_arrays_shifted[i] << std::endl;
        // }
        auto [g_list, g_list_grad] = KEnRef<double>::d_array_to_g_matrix(d_arrays_shifted, grouping_matrix/*, true*/);

        for (int i = 0; i < g_list.size(); i++) {
            std::cout << "g_list'\t" << g_list[i] << std::endl;
            if (g_list_grad.empty())
                continue;
            auto g_list_grad_i = g_list_grad[i];
            for (int j = 0; j < g_list_grad_i.size(); j++) {
                std::cout << "g_list_grad " << i + 1 << " " << j + 1 << std::endl;
                std::cout << g_list_grad_i[j] << std::endl;
            }
        }
        for (int i = 0; i < g_list.size(); ++i) {
            EXPECT_NEAR(expected_g_list_vector[interaction][i], g_list[i][0], expected_g_list_Epsilons_vector[interaction]);
        }

        //validation of a_matrix against the R code.
        auto g_matrix = KEnRef<double>::vectorOfVectors_to_Matrix(g_list);
        const auto &a_matrix = KEnRef<double>::g_matrix_to_a_matrix(g_matrix, dipole_kinetic_data[interaction].get_a_coef());
        std::cout << "a_matrix\n" << a_matrix << std::endl;
        for (int i = 0; i < a_matrix.cols(); ++i) {
            EXPECT_NEAR(expected_a_matrix_vector[interaction][i], a_matrix(0, i), expected_a_matrix_Epsilons_vector[interaction]);
        }

        Table rates_table({{"5.0e+08", "2.5e+08", "1.0e+12", "1.0e+04"}},
                          {{"kens", "kc", "kmethyl", "karo"}});
        // NamedRowVector<double> rates (rates_table.toNamedMatrix<double>()); //This line works, but I don't like it
        NamedRowVector<double> rates = rates_table.toNamedRowVector<double>();

        //calculate lambda eigenvalues
        //lambda_vec <- -colSums(rates[rownames(spec_den_data[["lambda_coef"]])]*spec_den_data[["lambda_coef"]])
        const auto &lambda_vector = KEnRef<double>::calculateLambdaVector(dipole_kinetic_data[interaction], rates);

        // update eigenvalues to account for molecular tumbling
        // NamedRowVector<double> lambda_prime_vec = (lambda_vector.array() - rates(KEnRef<double>::KC)).matrix();
        NamedRowVector<double> lambda_prime_vec(lambda_vector);
        lambda_prime_vec.array() -= rates(KEnRef<double>::KC);
        std::cout << lambda_prime_vec << std::endl;
        for (int i = 0; i < lambda_prime_vec.cols(); ++i) {
            EXPECT_NEAR(expected_lambda_prime_vector[interaction][i], lambda_prime_vec(0, i), expected_lambda_prime_vect_Epsilons_vector[interaction]);
        }
        const auto &[sigma, sigma_grad] = KEnRef<double>::a_matrix_to_sigma(a_matrix,lambda_prime_vec,700, true);
        EXPECT_NEAR(sigma(0,0), expected_sigma_vector[interaction], expected_sigma_vector_Epsilons_vector[interaction]);
        for (int i = 0; i < sigma_grad->size(); ++i) {
            EXPECT_NEAR(sigma_grad.value()(0, i), expected_sigma_grad[interaction][i], expected_sigma_grad_Epsilons_vector[interaction]);
        }
    }
}

TEST(KEnRefTestSuite, testCoordArrayToSigmaEnergy) {
    constexpr float n = 1; //0.25;
    constexpr float k = 1; //1e9;
    // proton field strength in MHz
    constexpr double proton_mhz = 700;
    const NamedRowVector<double> rates = Table(
        {{"5.0e+08", "2.5e+08", "1.0e+12", "1.0e+04"}},
        {{"kens", "kc", "kmethyl", "karo"}}
    ).toNamedRowVector<double>();

    static const auto FILENAME = "../../res/google_tests/2lum.pdb";

    const std::map<std::string, int> atomNameMapping_to1 = IoUtils::getAtomMappingFromPdb<std::string, int>( FILENAME, IoUtils::fill_atomId_to_index_Map);
    const bool handleNames = IoUtils::should_handleNames(atomNameMapping_to1);
    const std::map<std::string, int>& atomNameMapping = get_atomNameMapping<double>(atomNameMapping_to1, handleNames);

    auto spec_den_data_list = getSpecDenData<double>(handleNames);
    const std::vector<std::string> spec_den_data_prefixes  {"1-1", "1-2", "1-3", "2-2", "2-3", "3-3"};
    if (spec_den_data_list.size() != spec_den_data_prefixes.size()) {
        throw std::runtime_error("Data list and prefixes size mismatch");
    }
    std::vector<Table> tables{};
    // tables.reserve(spec_den_data_list.size());
    for (int i = 0; i < spec_den_data_list.size(); ++i) {
        std::string fileName = "../../res/google_tests/"+spec_den_data_prefixes[i]+"_atom_pairs.csv";
        std::cout << fileName << std::endl;
        const Table& table = IoUtils::readTable(fileName,true,false, "\\s*,\\s*", -1, true);
        std::vector<std::tuple<std::string, std::string>> atomPairs;
        atomPairs.reserve(table.rowCount());

        for (int j = 0; j < table.rowCount(); ++j) {
            // if (handleNames)
            //     atomPairs.emplace_back(IoUtils::normalizeName(table.at(j, 0), handleNames),
            //                            IoUtils::normalizeName(table.at(j, 1), handleNames));
            // else
            // atomPairs.emplace_back(IoUtils::normalizeName(table.at(j, 0), handleNames),
            //                        IoUtils::normalizeName(table.at(j, 1), handleNames));
            atomPairs.emplace_back(IoUtils::normalizeName(table.at(j, 0), handleNames),
                                    IoUtils::normalizeName(table.at(j, 1), handleNames));
        }
        auto& data = spec_den_data_list[i];
        data.setAtomPairs(atomPairs);
        NamedVector<double> sigmas(table.rowCount());
        for (int j = 0; j < sigmas.size(); ++j) {
            std::istringstream iss(table.at(j, 0));
            double value;
            iss >> value;
            sigmas(j, 0) = value;
        }
        data.set_sigmas(sigmas);
        tables.push_back(table);
    }

    //use models 4-6 to generate coords_synthetic
    const auto& coords_synthetic = getAllModels_allAtomCoordsMatrix<double>(FILENAME, {3,5});
    // pass coords_synthetic to coord_array_to_sigma to generate sigma0
    const auto &[sigma0, sigma0_grad] =
        KEnRef<double>::coord_array_to_sigma(coords_synthetic, rates, spec_den_data_list, proton_mhz, atomNameMapping, true,0);

    // verify that sigma0 is as expected from the files. No need to validate sigma0_grad.
    for (int i = 0; i < spec_den_data_list.size(); ++i) {
        const auto & specDenData = spec_den_data_list[i];
        const auto &atomPairs = specDenData.get_atom_pairs();
        for (int j = 0; j < sigma0[i].rows(); ++j) {
            // std::cout << i << ", " << j << "\t"<< sigma0[i](j,0) << "\t" << specDenData.getSigma(atomPairs.at(j)) << std::endl;
            double expectedValue = std::stod(tables[i].at(j,2));
            EXPECT_NEAR(sigma0[i](j,0), expectedValue, std::pow(10, static_cast<int>(log10(abs(expectedValue))) - 5));
        }
    }

    //set sigma0 in spec_den_data_list atom pairs
    for (int i = 0; i < spec_den_data_list.size(); ++i) {
        spec_den_data_list[i].set_sigmas(sigma0.at(i));
    }

    const auto& allModels_allAtomCoordsMap = getAllModels_allAtomCoordsMatrix<double>(FILENAME, {0,2});
    const auto &[sigma_energy, sigma_energy_grad] =
        KEnRef<double>::coord_array_to_sigma_energy(allModels_allAtomCoordsMap, rates, spec_den_data_list, proton_mhz,
        k, n, atomNameMapping, true, 0);
    // std::cout << "sigma_energy: " << sigma_energy << std::endl;
    // if (sigma_energy_grad.has_value()) {
    //     for (int i = 0; i < sigma_energy_grad->size(); ++i) {
    //         std::cout << "sigma_energy_grad["<<i<<"]: shape:(" <<sigma_energy_grad->at(i).rows()<<"," <<sigma_energy_grad->at(i).cols() << ")" << std::endl;
    //         std::cout << sigma_energy_grad->at(i) << std::endl;
    //     }
    // }
    EXPECT_NEAR(sigma_energy, 228.98, 1e-2);
    std::string fileName = "../../res/google_tests/2lum_sigma.csv";
    std::ifstream instream(fileName);
    if (!instream.is_open()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        throw std::runtime_error(std::string("Can't open file:").append(fileName));
    }
    std::vector<NamedMatrix<double>>modelsSigma(3);
    for (int model = 0; model < 3; ++model) {
        modelsSigma.at(model) = IoUtils::readTable(instream, true, true, "\\s*,\\s*", 423, false).toNamedMatrix<double>();
        // std::cout << modelsSigma.at(model) << std::endl ;
    }
    for (int m = 0; m < modelsSigma.size(); ++m) {
        auto &model_sigma = modelsSigma.at(m);
        auto &atomNames = model_sigma.rowNames();

        for (int i = 0; i < model_sigma.rows(); ++i) {
            for (int j = 0; j < model_sigma.cols(); ++j) {
                double expectedValue = model_sigma(i,j);
                auto normalizedAtomName = IoUtils::normalizeName(atomNames[i], true);
                double actualValue = sigma_energy_grad->at(m)(atomNameMapping.at(normalizedAtomName), j);
                EXPECT_NEAR(expectedValue, actualValue, std::pow(10, static_cast<int>(log10(abs(expectedValue))) - 5));
            }
        }
    }
}

TEST(KEnRefTestSuite, TestS2OrderParameters) {
    std::vector<std::string> files{
        "../../res/google_tests/ensemble_coord_model1.csv", "../../res/google_tests/ensemble_coord_model2.csv"
    };
    std::vector<CoordsMatrixType<KEnRef_Real_t> > coordsVector;
    coordsVector.reserve(files.size());

    for (int i = 0; i < files.size(); ++i) {
        auto tempCoordsData = IoUtils::readTable(files[i], false, false, ",");
        coordsVector.emplace_back(tempCoordsData.rowCount(), 3);
        for (int j = 0; j < tempCoordsData.rowCount(); ++j) {
            for (int k = 0; k < 3; ++k) {
                std::istringstream temp(tempCoordsData(j, k));
                temp >> coordsVector[i](j, k);
            }
        }
    }
    auto experimentalData_table = IoUtils::readTable(
        "../../res/google_tests/singleton_data_10nsstart+fit_0+10.csv", true, false, ",");

    std::vector<std::tuple<int, int> > atomIdPairs;
    auto expectedS2 = Eigen::VectorX<KEnRef_Real_t>(experimentalData_table.rowCount());
    int i1, i2;
    for (int i = 0; i < experimentalData_table.rowCount(); ++i) {
        std::istringstream temp1(experimentalData_table(i, "i1")), temp2(experimentalData_table(i, "i2")), temp3(experimentalData_table(i, "s2"));
        temp1 >> i1;
        temp2 >> i2;
        atomIdPairs.emplace_back(i1 - 1, i2 - 1);
        temp3 >> expectedS2(i);
    }

    KEnRef_Real_t epsilon;
    if constexpr (std::is_same_v<KEnRef_Real_t, float>) {
        epsilon = 5e-6;
    } else {
        epsilon = 1e-13;
    }
    std::cout << "testing float S2 difference with epsilon = " << epsilon << "\n";
    const auto &experimentalS2Double = KEnRef<KEnRef_Real_t>::s2OrderParams(coordsVector, atomIdPairs, 0);
    std::cout << "expectedS2    \t" << expectedS2.transpose() << std::endl;
    std::cout << "experimentalS2\t" << experimentalS2Double.transpose() << std::endl;
    TestHelper<KEnRef_Real_t>::EXPECT_MATRIX_NEAR(expectedS2, experimentalS2Double, epsilon);
}
