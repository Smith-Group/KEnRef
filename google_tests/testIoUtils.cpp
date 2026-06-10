#include <iostream>
#include <filesystem>
#include <fstream>
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <stdexcept>

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

// ---------------------------------------------------------------------------
// load_spec_den_relax_data: C++ port of read_spec_den_relax_data() in ke.R.
// Self-contained: each test writes its own *_atom_relax/_groupings/_a_coef/
// _lambda_coef CSVs into a temp dir that TearDown removes.
// ---------------------------------------------------------------------------
class SpecDenRelaxLoaderTest : public ::testing::Test {
protected:
    std::filesystem::path dir;

    void SetUp() override {
        const auto *info = ::testing::UnitTest::GetInstance()->current_test_info();
        dir = std::filesystem::temp_directory_path() / (std::string("kenref_relax_") + info->name());
        std::filesystem::remove_all(dir);
        std::filesystem::create_directories(dir);
    }
    void TearDown() override { std::filesystem::remove_all(dir); }

    static void wr(const std::filesystem::path &p, const std::string &content) {
        std::ofstream o(p);
        o << content;
    }
    // groupings / a_coef / lambda_coef — readable but otherwise immaterial to these parsing tests
    void writeAuxFiles(const std::string &prefix) const {
        wr(dir / (prefix + "_groupings.csv"), "1,1,2,2\n1,2,1,2\n");
        wr(dir / (prefix + "_a_coef.csv"), "kens,kc\n0.5,0.5\n0.25,0.75\n");
        wr(dir / (prefix + "_lambda_coef.csv"), ",kens,kc\nmode1,1.0,2.0\nmode2,3.0,4.0\n");
    }
};

// unit=TRUE headers, two rates (r1_400 = value+k+two terms; r2_500 = value-only+one term),
// N==M==2 (no aliasing). Also exercises the generalized, non-`\d+-\d+` prefix "gb3_15n".
TEST_F(SpecDenRelaxLoaderTest, LoadsMultipleRates) {
    const std::string prefix = "gb3_15n";
    const std::string header =
        "a1_unit,a2_unit,r1_400_value,r1_400_k,r1_400_wHmwN_coef,r1_400_wHmwN_freq,"
        "r1_400_wN_coef,r1_400_wN_freq,r2_500_value,r2_500_0_coef,r2_500_0_freq";
    wr(dir / (prefix + "_atom_relax.csv"),
       header + "\n"
       "\" N   MET A   3 \",\" H   MET A   3 \",1.5,4.0,0.1,1000,0.3,2000,5.5,0.7,0\n"
       "\" N   GLN A   4 \",\" H   GLN A   4 \",2.5,9.0,0.11,1100,0.33,2200,6.5,0.77,0\n");
    writeAuxFiles(prefix);

    // generalized prefix discovery
    const auto prefixes = IoUtils::find_spec_den_relax_data_prefixes(dir.string());
    ASSERT_EQ(prefixes.size(), 1u);
    EXPECT_EQ(prefixes[0], prefix);

    const auto data = IoUtils::load_spec_den_relax_data(dir.string(), false);
    ASSERT_EQ(data.size(), 1u);
    const auto &d = data.front();

    EXPECT_TRUE(d.is_unit());
    EXPECT_EQ(d.n_data_rows(), 2u);
    EXPECT_EQ(d.get_atom_pairs().size(), 2u);
    EXPECT_EQ(d.get_atom_relax_header(), header);   // verbatim round-trip

    const auto &rates = d.get_relax_data_list();
    ASSERT_EQ(rates.size(), 2u);
    EXPECT_EQ(rates[0].rate_name, "r1_400");        // CSV order preserved
    EXPECT_EQ(rates[1].rate_name, "r2_500");
    EXPECT_TRUE(d.has_rate("r1_400"));
    EXPECT_FALSE(d.has_rate("nope"));

    // r1_400: value + k + two terms (underscore-containing rate name parsed correctly)
    const auto &r1 = d.rate("r1_400");
    ASSERT_TRUE(r1.value.has_value());
    ASSERT_TRUE(r1.k.has_value());
    EXPECT_EQ(r1.coef().rows(), 2);
    EXPECT_EQ(r1.coef().cols(), 2);
    EXPECT_EQ(r1.coef().colNames(), (std::vector<std::string>{"wHmwN", "wN"}));
    EXPECT_EQ(r1.freq().colNames(), (std::vector<std::string>{"wHmwN", "wN"}));
    EXPECT_NEAR(r1.coef().coeff(0, 0), 0.1, 1e-9);    // wHmwN, row0
    EXPECT_NEAR(r1.coef().coeff(0, 1), 0.3, 1e-9);    // wN, row0
    EXPECT_NEAR(r1.coef().coeff(1, 1), 0.33, 1e-9);   // wN, row1
    EXPECT_NEAR(r1.freq().coeff(1, 0), 1100.0, 1e-9); // wHmwN, row1
    EXPECT_NEAR(r1.value->coeff(0), 1.5, 1e-9);
    EXPECT_NEAR(r1.k->coeff(1), 9.0, 1e-9);

    // r2_500: value only (no k), single term "0"
    const auto &r2 = d.rate("r2_500");
    EXPECT_TRUE(r2.value.has_value());
    EXPECT_FALSE(r2.k.has_value());
    EXPECT_EQ(r2.coef().cols(), 1);
    EXPECT_EQ(r2.coef().colNames(), (std::vector<std::string>{"0"}));
    EXPECT_NEAR(r2.coef().coeff(1, 0), 0.77, 1e-9);
    EXPECT_NEAR(r2.value->coeff(1), 6.5, 1e-9);
}

// No `_unit` suffix -> unit=FALSE; trailing atom-pair-only rows alias back into the
// first N=2 data rows (M=4), exactly as sigma wraps via `% N`.
TEST_F(SpecDenRelaxLoaderTest, WrapsAliasRows) {
    const std::string prefix = "3-5";
    wr(dir / (prefix + "_atom_relax.csv"),
       "a1,a2,r1_400_value,r1_400_wN_coef,r1_400_wN_freq\n"
       "\" N   MET A   3 \",\" H   MET A   3 \",1.5,0.3,2000\n"
       "\" N   GLN A   4 \",\" H   GLN A   4 \",2.5,0.33,2200\n"
       "\" N   MET A   3 \",\" H   MET A   3 \",,,\n"
       "\" N   GLN A   4 \",\" H   GLN A   4 \",,,\n");
    writeAuxFiles(prefix);

    const auto data = IoUtils::load_spec_den_relax_data(dir.string(), false);
    ASSERT_EQ(data.size(), 1u);
    const auto &d = data.front();

    EXPECT_FALSE(d.is_unit());
    EXPECT_EQ(d.n_data_rows(), 2u);              // only the first two rows carry relax data
    EXPECT_EQ(d.get_atom_pairs().size(), 4u);    // all rows are atom pairs (M = 2N)

    const auto &r1 = d.rate("r1_400");
    ASSERT_TRUE(r1.value.has_value());
    EXPECT_EQ(r1.value->rows(), 2);              // compact N-row block, not M
    EXPECT_EQ(r1.coef().rows(), 2);
    EXPECT_NEAR(r1.value->coeff(1), 2.5, 1e-9);

    // wrapping: trailing rows alias to the first N
    EXPECT_EQ(d.blockRow(static_cast<size_t>(0)), 0u);
    EXPECT_EQ(d.blockRow(static_cast<size_t>(2)), 0u);
    EXPECT_EQ(d.blockRow(static_cast<size_t>(3)), 1u);
}

// Exactly one of the first two columns ends in `_unit` -> must throw.
TEST_F(SpecDenRelaxLoaderTest, ThrowsOnUnitSuffixMismatch) {
    const std::string prefix = "3-5";
    wr(dir / (prefix + "_atom_relax.csv"),
       "a1_unit,a2,r1_400_wN_coef,r1_400_wN_freq\n"
       "\" N   MET A   3 \",\" H   MET A   3 \",0.3,2000\n");
    EXPECT_THROW(IoUtils::load_spec_den_relax_data(dir.string(), false), std::runtime_error);
}
