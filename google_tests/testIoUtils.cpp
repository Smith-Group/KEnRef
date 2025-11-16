#include <iostream>
#include <filesystem>
#include <fstream>
#include "gtest/gtest.h"
#include <vector>

#include "testHelper.h"
#include "core/IoUtils.h"

TEST(IoUtilsTestSuite, TestGetAtomNameMappingFromPdb){
    static const char *const ATOMNAME_MAPPING_FILENAME = "../../res/google_tests/00ns.pdb";
    std::cout << std::filesystem::current_path() << std::endl;
    std::map<std::string, int> atomNameMapping = IoUtils::getAtomMappingFromPdb<std::string, int>(ATOMNAME_MAPPING_FILENAME, IoUtils::fill_atomId_to_index_Map);
    EXPECT_TRUE(!atomNameMapping.empty());
    EXPECT_EQ(atomNameMapping.size(), 25951); // was 1231 before adding solvent
    EXPECT_EQ(atomNameMapping[" HA  MET     1 "], 6);
    EXPECT_EQ(atomNameMapping[" HB1 MET     1 "], 8);
    EXPECT_EQ(atomNameMapping[" HA  GLN     2 "], 23);
    EXPECT_EQ(atomNameMapping[" OC2 GLY    76 "], 1231);
}

TEST(IoUtilsTestSuite, testReadUniformTable){
    std::istringstream input(
            "1 2 3\n"
            "4 5 6\n"
            "7 8 9\n");
    const std::vector<std::vector<int>> expected {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    auto table = IoUtils::read_uniform_table_of<int>(input);
    for (int i = 0; i< 3; i++) {
        EXPECT_EQ(expected[i], table[i]);
        IoUtils::printVector(expected[i]);
    }
}

TEST(IoUtilsTestSuite, testStripEnclosingQuotes){
    std::string strin, strout, strexp;
    strin = "\"doublequote in the beginning";
    strexp = "\"doublequote in the beginning";
    strout = IoUtils::strip_enclosing_quotes(strin);
    std::cout << "Testing: " << strin << "-->" << strout << std::endl;
    EXPECT_EQ(strout, strexp);
    
    strin = "doublequote at the end\"";
    strexp = "doublequote at the end\"";
    strout = IoUtils::strip_enclosing_quotes(strin);
    std::cout << "Testing: " << strin << "-->" << strout << std::endl;
    EXPECT_EQ(strout, strexp);

    strin = "doublequote in the (\") middle";
    strexp = "doublequote in the (\") middle";
    strout = IoUtils::strip_enclosing_quotes(strin);
    std::cout << "Testing: " << strin << "-->" << strout << std::endl;
    EXPECT_EQ(strout, strexp);

    strin = "escaped doublequote at the end\\\"";
    strexp = "escaped doublequote at the end\\\"";
    strout = IoUtils::strip_enclosing_quotes(strin);
    std::cout << "Testing: " << strin << "-->" << strout << std::endl;
    EXPECT_EQ(strout, strexp);

    strin = "\"doublequote at both ends (should be stripped out)\"";
    strexp = "doublequote at both ends (should be stripped out)";
    strout = IoUtils::strip_enclosing_quotes(strin);
    std::cout << "Testing: " << strin << "-->" << strout << std::endl;
    EXPECT_EQ(strout, strexp);

    strin = "`escaped back quote at the end\\`";
    strexp = "`escaped back quote at the end\\`";
    strout = IoUtils::strip_enclosing_quotes(strin, '`');
    std::cout << "Testing: " << strin << "-->" << strout << std::endl;
    EXPECT_EQ(strout, strexp);

    strin = "`back quote at both ends (should be stripped out)`";
    strexp = "back quote at both ends (should be stripped out)";
    strout = IoUtils::strip_enclosing_quotes(strin, '`');
    std::cout << "Testing: " << strin << "-->" << strout << std::endl;
    EXPECT_EQ(strout, strexp);
}

TEST(IoUtilsTestSuite, testSubset_idx_to_grouping_mat) {
    std::vector<std::vector<std::vector<int>>> multipleGroupings{
        {{0,1,2,3,4,5,}},
        {{0,2,4,}, {1,3,5,}},
        {{0,1,}, {2,3,}, {4,5,}},
        {{0},{2},{4},{1},{3},{5},},
    };
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> expected(4,6);
    expected <<
    1, 1, 1, 1, 1, 1,
    1, 2, 1, 2, 1, 2,
    1, 1, 2, 2, 3, 3,
    1, 4, 2, 5, 3, 6;
    expected.array() -= 1;
    std::cout << "expeted\n" << expected << std::endl;
    std::cout << "calculated\n" << IoUtils::subset_idx_to_grouping_mat(multipleGroupings) << std::endl;

    TestHelper<KEnRef_Real_t>::EXPECT_MATRIX_EQ(expected, IoUtils::subset_idx_to_grouping_mat(multipleGroupings));
}
