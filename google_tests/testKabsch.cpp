#include "gtest/gtest.h"

// Just some dummy tests
TEST(KabschTestSuite, TestTranslate){
    EXPECT_EQ(737761, 737761);
}


TEST(KabschTestSuite, Fail){
    EXPECT_EQ(737761, 0);
}
