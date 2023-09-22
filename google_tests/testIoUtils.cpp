#include <iostream>
#include <filesystem>
#include "gtest/gtest.h"
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
    //for a larger file (later), test the last atom as well.
}


