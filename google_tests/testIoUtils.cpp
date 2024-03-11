#include <cstdio>
#include <iostream>
#include <filesystem>
#include <fstream>
#include "gtest/gtest.h"
#include <vector>
#include "core/IoUtils.h"

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

TEST(IoUtilsTestSuite, testReadParams) {
    std::istringstream input(
            "#note line\n"
            "# key = value\n"
            "noSpace=noSpace\n"
            " headingSpace= headingSpace\n"
            "trailingSpace  =trailingSpace    \n"
            " bothSpaces =  bothSpaces   \n"
            "internal Space=internal Space\n"
            " internal and external Space = internal and external Space \n"
            " spaces and a comment = spaces and a comment    # this is an inline comment\n"
            "non matching line\n"
            );
    auto map = IoUtils::readParams(input);
    EXPECT_EQ(map.size(), 7);
    EXPECT_EQ(map["noSpace"], "noSpace");
    EXPECT_EQ(map["headingSpace"], "headingSpace");
    EXPECT_EQ(map["trailingSpace"], "trailingSpace");
    EXPECT_EQ(map["bothSpaces"], "bothSpaces");
    EXPECT_EQ(map["internal Space"], "internal Space");
    EXPECT_EQ(map["internal and external Space"], "internal and external Space");
    EXPECT_EQ(map["spaces and a comment"], "spaces and a comment");
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