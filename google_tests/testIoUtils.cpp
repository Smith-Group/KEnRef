#include <cstdio>
#include <iostream>
#include <filesystem>
#include <fstream>
#include "gtest/gtest.h"
#include <vector>
#include "IoUtils.h"

TEST(IoUtilsTestSuite, TestGetAtomNameMappingFromPdb){
    static const char *const ATOMNAME_MAPPING_FILENAME = "../../res/cleanstart/6v5d_step0_for_atomname_mapping.pdb";
    std::cout << std::filesystem::current_path() << std::endl;
    std::map<std::string, int> atomNameMapping = IoUtils::getAtomNameMappingFromPdb(ATOMNAME_MAPPING_FILENAME);
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
    strout = IoUtils::strip_enclosing_quotoes(strin);
    std::cout << "Testing: " << strin << "-->" << strout << std::endl;
    EXPECT_EQ(strout, strexp);
    
    strin = "doublequote at the end\"";
    strexp = "doublequote at the end\"";
    strout = IoUtils::strip_enclosing_quotoes(strin);
    std::cout << "Testing: " << strin << "-->" << strout << std::endl;
    EXPECT_EQ(strout, strexp);

    strin = "doublequote in the (\") middle";
    strexp = "doublequote in the (\") middle";
    strout = IoUtils::strip_enclosing_quotoes(strin);
    std::cout << "Testing: " << strin << "-->" << strout << std::endl;
    EXPECT_EQ(strout, strexp);

    strin = "escaped doublequote at the end\\\"";
    strexp = "escaped doublequote at the end\\\"";
    strout = IoUtils::strip_enclosing_quotoes(strin);
    std::cout << "Testing: " << strin << "-->" << strout << std::endl;
    EXPECT_EQ(strout, strexp);

    strin = "\"doublequote at both ends (should be stripped out)\"";
    strexp = "doublequote at both ends (should be stripped out)";
    strout = IoUtils::strip_enclosing_quotoes(strin);
    std::cout << "Testing: " << strin << "-->" << strout << std::endl;
    EXPECT_EQ(strout, strexp);

    strin = "`escaped back quote at the end\\`";
    strexp = "`escaped back quote at the end\\`";
    strout = IoUtils::strip_enclosing_quotoes(strin, '`');
    std::cout << "Testing: " << strin << "-->" << strout << std::endl;
    EXPECT_EQ(strout, strexp);

    strin = "`back quote at both ends (should be stripped out)`";
    strexp = "back quote at both ends (should be stripped out)";
    strout = IoUtils::strip_enclosing_quotoes(strin, '`');
    std::cout << "Testing: " << strin << "-->" << strout << std::endl;
    EXPECT_EQ(strout, strexp);

}
TEST(IoUtilsTestSuite, restOfTests){
    //TODO make them tests
    std::ifstream in("../../res/noe_1.tsv");
    const auto& [group1names, group2names, values] = IoUtils::read_noe_table(in);
    for(int i = 0; i < group1names.size(); ++i){
        std::cout<< group1names[i] << "\t| " << group2names[i] << "\t| " << values[i] << std::endl;
    }

    std::ifstream noe_groups_file("../../res/noe_groups.R");
    auto noe_groups = IoUtils::read_noe_groups(noe_groups_file);
    for(auto [key, val] : noe_groups){
        std::cout << "<" << key<< ">"  << " " << val.size() << ":";
        for(const auto& str: val){
            std::cout << " [" << str << ']';
        }
        std::cout << std::endl;
    }
}