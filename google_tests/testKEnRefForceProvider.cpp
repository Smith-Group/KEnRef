#include <fstream>
#include "gtest/gtest.h"
#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include "testHelper.h"
#include "core/KEnRef.h"
#include "core/IoUtils.h"
#include "gmxinterface/KEnRefForceProvider.h"

Eigen::IOFormat fullPrecisionFmt(Eigen::FullPrecision);



TEST(KEnRefForceProviderTestSuite, TestRestoreNoJump){
    std::vector<std::string> files {"../../res/google_tests/ensemble_coord_10ns.csv", "../../res/google_tests/ensemble_coord_moldel2_PBC_Broken.csv"};
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
        }
    }

    std::vector<int> guideAtomsX0Indices = {4, 21, 38, 57, 77, 93, 115, 129, 148, 162, 169, 191, 205, 224, 238, 257, 272,
                                           288, 311, 317, 328, 340, 354, 373, 388, 402, 418, 440, 450, 472, 491, 508,
                                           520, 542, 557, 564, 591, 605, 611, 623, 640, 657, 681, 700, 719, 739, 749,
                                           756, 778, 795, 814, 829, 841, 848, 872, 886, 905, 916, 928, 949, 963, 982,
                                           999, 1021, 1036, 1047, 1061, 1080, 1097};
    std::vector<int> subAtomsX1Indices = {6, 8, 9, 11, 12, 21, 23, 25, 26, 28, 29, 33, 34, 38, 40, 42, 48, 49, 57, 59,
                                          61, 62, 69, 77, 79, 81, 93, 95, 97, 98, 100, 101, 103, 104, 106, 107, 115,
                                          117, 119, 129, 131, 133, 134, 136, 148, 150, 152, 162, 164, 165, 169, 171,
                                          173, 174, 176, 177, 179, 180, 182, 183, 191, 193, 195, 205, 207, 209, 215,
                                          216, 224, 226, 228, 238, 240, 242, 243, 245, 257, 259, 261, 262, 264, 265,
                                          272, 274, 276, 288, 290, 292, 293, 295, 296, 304, 305, 307, 308, 310, 311,
                                          313, 317, 319, 321, 322, 328, 330, 332, 333, 340, 342, 344, 354, 356, 358,
                                          364, 365, 373, 375, 377, 378, 380, 381, 388, 390, 392, 393, 397, 398, 402,
                                          404, 406, 418, 420, 422, 423, 425, 426, 428, 429, 431, 432, 440, 442, 450,
                                          452, 454, 455, 457, 458, 460, 461, 463, 464, 472, 474, 476, 482, 483, 491,
                                          493, 495, 496, 498, 499, 503, 504, 508, 510, 512, 513, 520, 522, 524, 525,
                                          527, 528, 530, 531, 533, 534, 542, 544, 546, 547, 549, 550, 557, 559, 560,
                                          564, 566, 568, 574, 575, 584, 585, 587, 588, 590, 591, 593, 598, 599, 601,
                                          602, 604, 605, 607, 611, 613, 615, 616, 623, 625, 627, 628, 630, 631, 635,
                                          636, 640, 642, 644, 645, 647, 648, 652, 653, 657, 659, 661, 662, 664, 665,
                                          667, 668, 681, 683, 685, 686, 688, 700, 702, 704, 710, 711, 719, 721, 723,
                                          724, 731, 739, 741, 749, 751, 752, 756, 758, 760, 761, 763, 764, 766, 767,
                                          769, 770, 778, 780, 782, 783, 785, 786, 790, 791, 795, 797, 799, 800, 802,
                                          814, 816, 818, 819, 821, 822, 829, 831, 833, 834, 841, 843, 844, 848, 850,
                                          852, 853, 855, 856, 858, 859, 872, 874, 876, 886, 888, 890, 891, 893, 905,
                                          907, 909, 910, 916, 918, 920, 921, 928, 930, 932, 933, 949, 951, 953, 954,
                                          958, 959, 963, 965, 967, 973, 974, 982, 984, 986, 987, 989, 990, 994, 995,
                                          999, 1001, 1003, 1004, 1006, 1007, 1009, 1010, 1012, 1013, 1021, 1023, 1025,
                                          1026, 1028, 1029, 1036, 1038, 1040, 1041, 1047, 1049, 1051, 1061, 1063, 1065,
                                          1066, 1068, 1080, 1082, 1084, 1085, 1089, 1093, 1097, 1099, 1101, 1102, 1104,
                                          1116, 1118, 1120, 1132, 1134, 1136, 1137, 1139, 1151, 1153, 1155, 1156, 1158,
                                          1159, 1161, 1162, 1175, 1177, 1179, 1180, 1182, 1194, 1196, 1198, 1199, 1201,
                                          1202, 1204, 1205, 1218, 1220, 1221, 1225};
    std::vector<CoordsMatrixType<KEnRef_Real_t>> all_guideAtomsX_ZEROIndexed{{guideAtomsX0Indices.size(), 3},
                                                                             {guideAtomsX0Indices.size(), 3}};
    CoordsMatrixType<KEnRef_Real_t> guideAtomsX_ZEROIndexed(guideAtomsX0Indices.size(), 3);
    for (int i = 0; i < guideAtomsX0Indices.size(); ++i) {
        for (int j = 0; j < 2; ++j) {
            all_guideAtomsX_ZEROIndexed[j].row(i) = coordsVector[j].row(guideAtomsX0Indices[i]);
        }
    }

    std::vector<CoordsMatrixType<KEnRef_Real_t>> all_subAtomsX{{subAtomsX1Indices.size(), 3},
                                                               {subAtomsX1Indices.size(), 3}};
    for (int i = 0; i < subAtomsX1Indices.size(); ++i) {
        for (int j = 0; j < 2; ++j) {
            all_subAtomsX[j].row(i) = coordsVector[j].row(subAtomsX1Indices[i]-1);
        }
    }
//    typedef KEnRef_Real_t matrix[3][3];
//    matrix box_ = {{6.40261984, 0,          0},
//                   {0,          6.40261984, 0},
//                   {0,          0,          6.40261984}};

    matrix box_ = {{6.40261993, 0,          0},
                   {0,          6.40261993, 0},
                   {0,          0,          6.40261993}};

//    matrix box_ = {{6.4026, 0,          0},
//                   {0,          6.4026, 0},
//                   {0,          0,          6.4026}};


    KEnRef_Real_t epsilon;
    if constexpr (std::is_same_v<KEnRef_Real_t, float>){
//        epsilon = 5e-6;
        epsilon = 0.119;
//        epsilon = 0.10;
    }else{
        epsilon = 1e-13;
    }
    std::cout << "Testing no jump on C-Alpha guide atoms" << std::endl;
    CoordsMatrixType<KEnRef_Real_t> guideAtomsDiffB4NoJump = all_guideAtomsX_ZEROIndexed[1] - all_guideAtomsX_ZEROIndexed[0];
    std::cout << "differences before restoreNoJump:\n" << guideAtomsDiffB4NoJump.transpose() << std::endl;
    KEnRefForceProvider::restoreNoJump(all_guideAtomsX_ZEROIndexed[1], all_guideAtomsX_ZEROIndexed[0], box_, true, 1);
    TestHelper<KEnRef_Real_t>::EXPECT_MATRIX_NEAR(all_guideAtomsX_ZEROIndexed[1], all_guideAtomsX_ZEROIndexed[0], epsilon);
    CoordsMatrixType<KEnRef_Real_t> guideAtomsDiffAfterNoJump = all_guideAtomsX_ZEROIndexed[1] - all_guideAtomsX_ZEROIndexed[0];
    std::cout << "differences after  restoreNoJump:\n" << guideAtomsDiffAfterNoJump.transpose() << std::endl;


    std::cout << "\n\nTesting no jump on Hydrogen subset of atoms" << std::endl;
    CoordsMatrixType<KEnRef_Real_t> subAtomsDiffB4NoJump = all_subAtomsX[1] - all_subAtomsX[0];
    std::cout << "differences before restoreNoJump:\n" << subAtomsDiffB4NoJump.transpose() << std::endl;
    KEnRefForceProvider::restoreNoJump(all_subAtomsX[1], all_subAtomsX[0], box_, true, 1);
    TestHelper<KEnRef_Real_t>::EXPECT_MATRIX_NEAR(all_subAtomsX[1], all_subAtomsX[0], epsilon);
    CoordsMatrixType<KEnRef_Real_t> subAtomsDiffAfterNoJump = all_subAtomsX[1] - all_subAtomsX[0];
    std::cout << "differences after  restoreNoJump:\n" << subAtomsDiffAfterNoJump.transpose() << std::endl;

//    std::cout << "testing float restoreNoJump with epsilon = " << epsilon << "\n";
//    const auto &experimentalS2Double = KEnRef<KEnRef_Real_t>::s2OrderParams(coordsVector, atomIdPairs,0);
//    std::cout << "expectedS2    \t" << expectedS2.transpose() << std::endl;
//    std::cout << "experimentalS2\t" << experimentalS2Double.transpose() << std::endl;
//    TestHelper<KEnRef_Real_t>::EXPECT_MATRIX_NEAR(expectedS2, experimentalS2Double, epsilon);

}

