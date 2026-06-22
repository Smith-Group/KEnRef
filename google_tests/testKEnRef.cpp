#include "gtest/gtest.h"
#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "core/KEnRef.h"
#include "core/IoUtils.h"
#include "core/EnergyModel.h"
#include "core/EngineAdapter.h"
#include "core/ModelRegistry.h"
#include "testHelper.h"

// Minimal EngineAdapter for driving an EnergyModel's init in unit tests: only getRawParam() is
// exercised (during buildCache); the per-step methods are unused here and assert if called.
class TestEngineAdapter : public kenref::EngineAdapter<KEnRef_Real_t> {
public:
    explicit TestEngineAdapter(std::map<std::string, std::string> params) : params_(std::move(params)) {}
    [[nodiscard]] std::optional<std::string> getRawParam(const std::string& key) const override {
        const auto it = params_.find(key);
        return it == params_.end() ? std::nullopt : std::optional<std::string>(it->second);
    }
    [[nodiscard]] int numModelsInThisProcess() const override { return 0; }
    void getLocalModelX(int, CoordsMatrixType<KEnRef_Real_t>&, CoordsMatrixType<KEnRef_Real_t>&,
                        Eigen::Matrix<KEnRef_Real_t, 3, 3>&) const override { ADD_FAILURE(); }
    void addLocalModelDerivatives(int, const CoordsMatrixType<KEnRef_Real_t>&) override { ADD_FAILURE(); }
    [[nodiscard]] int numModelsTotal() const override { return 0; }
    [[nodiscard]] int simulationIndex() const override { return 0; }
    void gatherFittedSubAtomsX(const std::vector<CoordsMatrixType<KEnRef_Real_t>>&,
                               std::vector<CoordsMatrixType<KEnRef_Real_t>>&) const override { ADD_FAILURE(); }
    void scatterModelDerivatives(const std::vector<CoordsMatrixType<KEnRef_Real_t>>&,
                                 std::vector<CoordsMatrixType<KEnRef_Real_t>>&) const override { ADD_FAILURE(); }
private:
    std::map<std::string, std::string> params_;
};

// Benchmark helper: OpenMP thread counts to sweep, overridable via KENREF_BENCH_THREADS
// (comma-separated, e.g. "1" to pin a clean single-thread perf profile, or "1,8,0").
static std::vector<int> benchThreadList(const std::vector<int>& dflt) {
    if (const char* e = std::getenv("KENREF_BENCH_THREADS")) {
        std::vector<int> v; std::stringstream ss(e); std::string tok;
        while (std::getline(ss, tok, ',')) if (!tok.empty()) v.push_back(std::stoi(tok));
        if (!v.empty()) return v;
    }
    return dflt;
}

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
        NamedMatrix<double>{(Eigen::Matrix<double, 1, 2>() <<
        0, 1
        ).finished(),
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
        NamedMatrix<double> { (Eigen::Matrix<double, 2, 4>() <<
                0, 0, 1, 1,
                0, 1, 0, 1).finished(),
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
        NamedMatrix<double>{ (Eigen::Matrix<double, 2, 4>() <<
             0, 1, 0, 1,
             0, 0, 1, 1
            ).finished(),
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
        NamedMatrix<double> { (Eigen::Matrix<double, 2, 6>() <<
            0, 0, 0, 1, 1, 1,
            0, 1, 2, 0, 1, 2
            ).finished(),
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
        NamedMatrix<double> { (Eigen::Matrix<double, 3, 6>() <<
                0, 0, 1, 1, 0, 1,
                0, 0, 0, 0, 1, 1,
                0, 1, 0, 1, 0, 0
            ).finished(),
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
        NamedMatrix<double> { (Eigen::Matrix<double, 2, 6>() <<
                0, 1, 0, 1, 0, 1,
                0, 0, 1, 1, 2, 2
            ).finished(),
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
std::vector<int> generate_contiguous_indixes(std::tuple<int, int> boundsInclusive) {
    auto [first, last] = boundsInclusive;
    std::vector<int> indices(last - first + 1);
    std::iota(indices.begin(), indices.end(), first);
    return indices;
}

template <typename KEnRef_Real>
std::vector<CoordsMatrixType<KEnRef_Real>> getAllModels_allAtomCoordsMatrix(
        const std::string& FILENAME, const std::vector<int>& modelIndices) {
    const auto atomNameMapping_to1 = IoUtils::getAtomMappingFromPdb<std::string, int>(
        FILENAME, IoUtils::fill_atomId_to_index_Map);
    auto allAtomCoordsMap_raw_vector = IoUtils::getMultipleAtomMappingFromPdb<int, Eigen::RowVector3<double>>(
        FILENAME, IoUtils::fill_atomIndex1_to_coords_Map<double>);

    int maxId0 = -1;
    for (const auto& [atomName, id1] : atomNameMapping_to1)
        maxId0 = std::max(maxId0, id1 - 1);
    std::vector<int> allModelAtomIdsVector(maxId0 + 1);
    std::iota(allModelAtomIdsVector.begin(), allModelAtomIdsVector.end(), 0);
    std::vector<CoordsMatrixType<double>> result;
    result.reserve(modelIndices.size());
    for (int idx : modelIndices)
        result.emplace_back(IoUtils::extractCoords(allModelAtomIdsVector, false, allAtomCoordsMap_raw_vector[idx], true));
    return result;
}

template<typename KEnRef_Real>
std::vector<CoordsMatrixType<KEnRef_Real> > getAllModels_allAtomCoordsMatrix(
    const std::string &FILENAME, const std::tuple<int, int> boundsInclusive) {
    return getAllModels_allAtomCoordsMatrix<KEnRef_Real>(
        FILENAME, generate_contiguous_indixes<KEnRef_Real>(boundsInclusive));
}

// ── Shared GB3 test fixture ───────────────────────────────────────────────────

static constexpr const char*  GB3_FILENAME       = "../../res/google_tests/2lum.pdb";
static constexpr const char*  GB3_PROTON_FILENAME = "../../res/google_tests/2lum_subset_proton.pdb";
static constexpr double       GB3_PROTON_MHZ  = 700.0;
static constexpr const char*  GB3_SPEC_DEN_DIR   = "../../res/google_tests/";
static const std::vector<std::string> SPEC_DEN_PREFIXES{"1-1","1-2","1-3","2-2","2-3","3-3"};

struct GB3SigmaEnergySetup {
    NamedRowVector<double>              rates;
    std::map<std::string, int>          atomNameMapping;
    bool                                handleNames;
    std::vector<SpecDenData<double>>    spec_den_data_list;  // sigma0 loaded from CSV
    std::vector<CoordsMatrixType<double>> coord_array;       // models {0,2,4} (0-based), unscaled
};

static GB3SigmaEnergySetup makeGB3SigmaEnergySetup() {
    const NamedRowVector<double> rates = Table(
        {{"5.0e+08", "2.5e+08", "1.0e+12", "1.0e+04"}},
        {{"kens", "kc", "kmethyl", "karo"}}
    ).toNamedRowVector<double>();

    const auto atomNameMapping_to1 = IoUtils::getAtomMappingFromPdb<std::string, int>(
        GB3_PROTON_FILENAME, IoUtils::fill_atomId_to_index_Map);
    const bool handleNames = IoUtils::should_handleNames(atomNameMapping_to1);
    const auto atomNameMapping = get_atomNameMapping<double>(atomNameMapping_to1, handleNames);

    // Load spec_den_data from CSVs — atom pairs, groupings, a_coef, lambda_coef, sigma all from R
    auto spec_den_data_list = IoUtils::load_spec_den_data(GB3_SPEC_DEN_DIR, handleNames);

    // coord_array: models {0,2,4} (0-based) = R's c(1,3,5) (1-based)
    return {rates, atomNameMapping, handleNames, std::move(spec_den_data_list),
            getAllModels_allAtomCoordsMatrix<double>(GB3_PROTON_FILENAME, {0, 2, 4})};
}

// ─────────────────────────────────────────────────────────────────────────────

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

    const auto &[toy_d_array, toy_d_array_grad] = KEnRef<KEnRef_Real_t>::r_array_to_d_array(toy_r_mat, false, true);
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
    auto eros3_sub_atom_idPairs = KEnRef<float>::atomNamePairs_2_atomIdPairs(atomNamePairs, atomNameMapping0);
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
        eros3_sub_r_array, false, true);
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
    const std::vector<CoordsMatrixType<double> >& allModels_allAtomCoordsMap = getAllModels_allAtomCoordsMatrix<double>(FILENAME, std::make_tuple(0, 2));
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
    // calculateLambdaVector is now a deterministic GEMV; lambda_prime matches the R reference exactly
    // (verified at tol=0). Tightened from the old 1e7/1e3 to 1e-2 on the large-magnitude interactions
    // (~1e12-2e12, where double granularity is ~5e-4, so 1e-2 stays cross-platform-safe); the already
    // tight 1e-8 on interactions 1 and 3 (max ~7.5e8) is kept.
    std::vector<double> expected_lambda_prime_vect_Epsilons_vector{1e-2, 1e-8, 1e-2, 1e-8, 1e-2, 1e-2};

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

        const auto &[d_arrays, d_arrays_grad] = KEnRef<double>::r_array_to_d_array(r_arrays, false, true, false);
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
    constexpr double k = 1;
    constexpr double n = 1;
    auto s = makeGB3SigmaEnergySetup();

    // Validate sigma0: C++ from synthetic models {1,3,5} must match R-generated CSV values
    auto coords_synthetic = getAllModels_allAtomCoordsMatrix<double>(GB3_PROTON_FILENAME, {1, 3, 5});
    auto [sigma0, unused_grad] = KEnRef<double>::coord_array_to_sigma(
        coords_synthetic, s.rates, s.spec_den_data_list, GB3_PROTON_MHZ, s.atomNameMapping, false, 0);
    (void)unused_grad;
    for (int i = 0; i < (int)s.spec_den_data_list.size(); ++i)
        for (int j = 0; j < sigma0[i].rows(); ++j) {
            const double expected = s.spec_den_data_list[i].sigmas()(j, 0);
            EXPECT_NEAR(sigma0[i](j, 0), expected,
                        std::pow(10, static_cast<int>(log10(std::abs(expected))) - 5));
        }

    // Run coord_array_to_sigma_energy on models {0,2,4} (R's c(1,3,5))
    auto coord_array = s.coord_array;  // copy: function scales in-place
    const auto& [sigma_energy, sigma_energy_grad] = KEnRef<double>::coord_array_to_sigma_energy(
        coord_array, s.rates, s.spec_den_data_list, GB3_PROTON_MHZ, k, n, s.atomNameMapping, true, 0);

    EXPECT_NEAR(sigma_energy, 124.8785, 1e-4);

    // Validate gradient against expected values from 2lum_sigma.csv
    std::ifstream instream("../../res/google_tests/2lum_sigma.csv");
    ASSERT_TRUE(instream.is_open()) << "Cannot open 2lum_sigma.csv";
    std::vector<NamedMatrix<double>> modelsSigma(3);
    for (int model = 0; model < 3; ++model)
        modelsSigma[model] = IoUtils::readTable(instream, true, true, "\\s*,\\s*", 423, false).toNamedMatrix<double>();
    for (int m = 0; m < (int)modelsSigma.size(); ++m) {
        const auto& model_sigma = modelsSigma[m];
        const auto& atomNames   = model_sigma.rowNames();
        for (int i = 0; i < model_sigma.rows(); ++i)
            for (int j = 0; j < model_sigma.cols(); ++j) {
                const double expected = model_sigma(i, j);
                const auto normName   = IoUtils::normalizeName(atomNames[i], true);
                const double actual   = sigma_energy_grad->at(m)(s.atomNameMapping.at(normName), j);
                EXPECT_NEAR(expected, actual, 1e-9);
            }
    }
}

// Drive the SIGMA model through the registry (buildCache -> finalizeIndexing -> compute) and confirm
// it reproduces the direct-kernel energy 124.8785 from testCoordArrayToSigmaEnergy. The model is a
// thin wrapper over coord_array_to_sigma_energy; here coord_array is the full-atom set and the
// name->id map doubles as the (identity) sub-indexing, so the wrapped call matches the direct one.
TEST(KEnRefTestSuite, testSigmaModelViaRegistry) {
    kenref::bootstrapModels();
    auto model = kenref::ModelRegistry<KEnRef_Real_t>::create("SIGMA");
    ASSERT_NE(model, nullptr);

    auto s = makeGB3SigmaEnergySetup();

    TestEngineAdapter adapter({
        {"EXP_DATA_FOLDER", GB3_SPEC_DEN_DIR},
        {"PROTON_MHZ", std::to_string(GB3_PROTON_MHZ)},
    });
    // Empty schema is fine: the adapter supplies every param the model reads, so no defaults needed.
    kenref::ParamSchema emptySchema;
    kenref::InitContext<KEnRef_Real_t> initCtx{adapter, emptySchema, s.atomNameMapping, s.handleNames, 3};
    model->buildCache(initCtx);

    kenref::IndexingContext<KEnRef_Real_t> indexCtx{s.atomNameMapping};
    model->finalizeIndexing(indexCtx);

    auto coord_array = s.coord_array;  // copy: the kernel scales in place
    kenref::StepContext<KEnRef_Real_t> stepCtx{coord_array, /*k*/ 1.0, /*n*/ 1.0, /*gradient*/ true, 0};
    auto [energy, grad] = model->compute(stepCtx);
    EXPECT_NEAR(energy, 124.8785, 1e-4);
    EXPECT_TRUE(grad.has_value());
    EXPECT_EQ(model->forceUnitScale(), 10.0);
}

// Drive the PLATEAUS model through the registry and confirm it forwards correctly to
// coord_array_to_g_energy. We build a small EXP_DATA_FILE from real atom-name pairs (valid keys of
// the 2lum mapping) with synthetic g1/g2, then compare model->compute() against a direct kernel call
// using the model's own atomIdPairs() + the same g0/grouping. (The kernel is R-validated elsewhere;
// here we test the wrapper: data-file reading, grouping, g0, and dispatch.)
TEST(KEnRefTestSuite, testPlateausModelViaRegistry) {
    kenref::bootstrapModels();
    auto model = kenref::ModelRegistry<KEnRef_Real_t>::create("PLATEAUS");
    ASSERT_NE(model, nullptr);

    auto s = makeGB3SigmaEnergySetup();  // atomNameMapping + coord_array (3 models) + handleNames

    // Gather up to 30 real atom-name pairs (already valid keys of s.atomNameMapping) and synthetic g.
    std::vector<std::tuple<std::string, std::string>> namePairs;
    std::vector<std::pair<double, double>> gvals;
    for (const auto& d : s.spec_den_data_list) {
        for (const auto& ap : d.get_atom_pairs()) {
            const int i = static_cast<int>(namePairs.size());
            namePairs.push_back(ap);
            gvals.emplace_back(0.01 + 0.001 * i, 0.02 + 0.001 * i);
            if (namePairs.size() >= 30) break;
        }
        if (namePairs.size() >= 30) break;
    }
    ASSERT_FALSE(namePairs.empty());

    const std::string path = "/tmp/kenref_plateaus_test.csv";
    {
        std::ofstream o(path);
        o << "\"\",\"atom1\",\"atom2\",\"g1\",\"g2\"\n";
        for (size_t i = 0; i < namePairs.size(); ++i)
            o << '"' << i << "\",\"" << std::get<0>(namePairs[i]) << "\",\"" << std::get<1>(namePairs[i])
              << "\"," << gvals[i].first << ',' << gvals[i].second << '\n';
    }

    TestEngineAdapter adapter({{"EXP_DATA_FILE", path}});
    kenref::ParamSchema emptySchema;
    kenref::InitContext<KEnRef_Real_t> initCtx{adapter, emptySchema, s.atomNameMapping, s.handleNames, 3};
    model->buildCache(initCtx);
    kenref::IndexingContext<KEnRef_Real_t> indexCtx{s.atomNameMapping};
    model->finalizeIndexing(indexCtx);

    // Reconstruct g0 + grouping the same way the model did, for the direct comparison.
    Eigen::MatrixX<double> g0(namePairs.size(), 2);
    for (int i = 0; i < (int)namePairs.size(); ++i) { g0(i, 0) = gvals[i].first; g0(i, 1) = gvals[i].second; }
    const std::vector<std::vector<std::vector<int>>> grouping{{{0, 1, 2}}, {{0}, {1}, {2}}};

    auto coords_model = s.coord_array;   // copies (g_energy does not scale in place, but keep parallel)
    auto coords_direct = s.coord_array;
    kenref::StepContext<KEnRef_Real_t> stepCtx{coords_model, /*k*/ 5e8, /*n*/ 0.25, /*gradient*/ true, 0};
    auto [energy_model, grad_model] = model->compute(stepCtx);

    auto [energy_direct, grad_direct] = KEnRef<double>::coord_array_to_g_energy(
        coords_direct, model->atomIdPairs(), grouping, g0, 5e8, 0.25, true, 0);

    EXPECT_NEAR(energy_model, energy_direct, 1e-9);
    EXPECT_EQ(model->forceUnitScale(), 1.0);
    ASSERT_TRUE(grad_model.has_value());
    ASSERT_TRUE(grad_direct.has_value());
    for (size_t m = 0; m < grad_model->size(); ++m)
        TestHelper<double>::EXPECT_MATRIX_NEAR((*grad_model)[m], (*grad_direct)[m], 1e-9);
}

// Relaxation-rate restraint energy + gradient end-to-end check against R's
// coord_array_to_relax_energy() (demo/relax_deriv_check.R, isotropic diffusion).
TEST(KEnRefTestSuite, testCoordArrayToRelaxEnergy) {
    constexpr double k = 1;     // R power_scaled_loss default
    constexpr double n = 1;     // R power_scaled_loss default p=1
    constexpr const char* RELAX_DIR = "../../res/google_tests/relax/";

    // rates: kens, kmethyl, karo, Dx=Dy=Dz (isotropic), from relax_deriv_check.R
    const NamedRowVector<double> rates = Table(
        {{"5.0e+08", "1.0e+12", "1.0e+04", "2.5e+08", "2.5e+08", "2.5e+08"}},
        {{"kens", "kmethyl", "karo", "Dx", "Dy", "Dz"}}
    ).toNamedRowVector<double>();

    const auto atomNameMapping_to1 = IoUtils::getAtomMappingFromPdb<std::string, int>(
        GB3_PROTON_FILENAME, IoUtils::fill_atomId_to_index_Map);
    const bool handleNames = IoUtils::should_handleNames(atomNameMapping_to1);
    const auto atomNameMapping = get_atomNameMapping<double>(atomNameMapping_to1, handleNames);

    auto relax_data_list = IoUtils::load_spec_den_relax_data(RELAX_DIR, handleNames);
    ASSERT_EQ(relax_data_list.size(), 6u);

    // coord_array: models {0,2,4} (R's c(1,3,5))
    auto coord_array = getAllModels_allAtomCoordsMatrix<double>(GB3_PROTON_FILENAME, {0, 2, 4});

    const auto& [relax_energy, relax_energy_grad] = KEnRef<double>::coord_array_to_relax_energy(
        coord_array, rates, relax_data_list, k, n, atomNameMapping, true, 0);

    // expected energy from res/google_tests/relax/relax_energy.txt (R reference)
    EXPECT_NEAR(relax_energy, 1840.2342988885609, 1e-2);

    // gradient vs 2lum_relax_grad.csv (3 model-blocks of 423 rows: name, X, Y, Z)
    std::ifstream instream(std::string(RELAX_DIR) + "2lum_relax_grad.csv");
    ASSERT_TRUE(instream.is_open()) << "Cannot open 2lum_relax_grad.csv";
    ASSERT_TRUE(relax_energy_grad.has_value());
    for (int m = 0; m < 3; ++m) {
        const auto model_grad = IoUtils::readTable(instream, true, true, "\\s*,\\s*", 423, false).toNamedMatrix<double>();
        const auto& atomNames = model_grad.rowNames();
        for (int i = 0; i < model_grad.rows(); ++i)
            for (int j = 0; j < model_grad.cols(); ++j) {
                const double expected = model_grad(i, j);
                const auto normName   = IoUtils::normalizeName(atomNames[i], handleNames);
                const double actual   = relax_energy_grad->at(m)(atomNameMapping.at(normName), j);
                EXPECT_NEAR(expected, actual, 1e-4) << "model " << m << " atom " << atomNames[i] << " xyz " << j;
            }
    }
}

// coord_array_to_relax(): predicted per-pair relaxation rates vs R reference
// (relax_pred.csv). Compared as a sorted multiset so the test is independent of
// substructure-discovery order.
TEST(KEnRefTestSuite, testCoordArrayToRelax) {
    constexpr const char* RELAX_DIR = "../../res/google_tests/relax/";
    const NamedRowVector<double> rates = Table(
        {{"5.0e+08", "1.0e+12", "1.0e+04", "2.5e+08", "2.5e+08", "2.5e+08"}},
        {{"kens", "kmethyl", "karo", "Dx", "Dy", "Dz"}}
    ).toNamedRowVector<double>();

    const auto atomNameMapping_to1 = IoUtils::getAtomMappingFromPdb<std::string, int>(
        GB3_PROTON_FILENAME, IoUtils::fill_atomId_to_index_Map);
    const bool handleNames = IoUtils::should_handleNames(atomNameMapping_to1);
    const auto atomNameMapping = get_atomNameMapping<double>(atomNameMapping_to1, handleNames);

    auto relax_data_list = IoUtils::load_spec_den_relax_data(RELAX_DIR, handleNames);
    auto coord_array = getAllModels_allAtomCoordsMatrix<double>(GB3_PROTON_FILENAME, {0, 2, 4});

    const auto relax = KEnRef<double>::coord_array_to_relax(
        coord_array, rates, relax_data_list, atomNameMapping, 0);

    std::vector<double> actual;
    for (const auto& m : relax)
        for (int j = 0; j < m.cols(); ++j)
            for (int i = 0; i < m.rows(); ++i)
                actual.push_back(m(i, j));

    // expected values from R's coord_array_to_relax (relax_pred.csv, column "value")
    auto predTable = IoUtils::readTable(std::string(RELAX_DIR) + "relax_pred.csv", true, false, "\\s*,\\s*", -1, false);
    std::vector<double> expected;
    expected.reserve(predTable.rowCount());
    for (int r = 0; r < predTable.rowCount(); ++r)
        expected.push_back(std::stod(predTable.at(r, "value")));

    ASSERT_EQ(actual.size(), expected.size());
    std::sort(actual.begin(), actual.end());
    std::sort(expected.begin(), expected.end());
    for (size_t i = 0; i < actual.size(); ++i)
        EXPECT_NEAR(actual[i], expected[i], 1e-6) << "sorted index " << i;
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

TEST(KEnRefTestSuite, TestCoordArrayToSigmaEnergyFD) {
    constexpr double k = 1;
    constexpr double n = 1;
    auto s = makeGB3SigmaEnergySetup();

    // Per-interaction FD check (mirrors R's gradient_list loop over spec_den_data_ind_list)
    constexpr double delta = 1e-5;
    double overall_max_analytical = 0.0, overall_max_diff = 0.0;

    for (int item_idx = 0; item_idx < (int)s.spec_den_data_list.size(); ++item_idx) {
        const std::string& prefix = SPEC_DEN_PREFIXES[item_idx];
        const std::vector<SpecDenData<double>> item_list{s.spec_den_data_list[item_idx]};

        // Active atoms for this interaction
        std::set<int> active_ids_set;
        for (const auto& [a1, a2] : s.spec_den_data_list[item_idx].get_atom_pairs()) {
            auto it1 = s.atomNameMapping.find(a1), it2 = s.atomNameMapping.find(a2);
            if (it1 != s.atomNameMapping.end()) active_ids_set.insert(it1->second);
            if (it2 != s.atomNameMapping.end()) active_ids_set.insert(it2->second);
        }
        const std::vector<int> active_atom_ids(active_ids_set.begin(), active_ids_set.end());

        // Analytical gradient for this interaction
        auto coord_copy = s.coord_array;
        const auto& [energy, grad] = KEnRef<double>::coord_array_to_sigma_energy(
            coord_copy, s.rates, item_list, GB3_PROTON_MHZ, k, n, s.atomNameMapping, true, 0);
        ASSERT_TRUE(grad.has_value()) << "No gradient for " << prefix;

        // FD gradient for this interaction
        const auto fd = TestHelper<double>::finite_difference_grad(
            s.coord_array,
            [&](std::vector<CoordsMatrixType<double>> coords) -> double {
                return std::get<0>(KEnRef<double>::coord_array_to_sigma_energy(
                    coords, s.rates, item_list, GB3_PROTON_MHZ, k, n,
                    s.atomNameMapping, false, 0));
            },
            active_atom_ids, delta);

        // Collect (fd, analytical) pairs, compute stats
        std::vector<std::pair<double, double>> pairs;
        double max_anal = 0.0, max_diff = 0.0;
        for (int m = 0; m < (int)grad->size(); ++m)
            for (int atomId : active_atom_ids)
                for (int j = 0; j < 3; ++j) {
                    const double av = grad->at(m)(atomId, j);
                    const double fv = fd[m](atomId, j);
                    pairs.emplace_back(fv, av);
                    max_anal = std::max(max_anal, std::abs(av));
                    max_diff = std::max(max_diff, std::abs(av - fv));
                }

        overall_max_analytical = std::max(overall_max_analytical, max_anal);
        overall_max_diff       = std::max(overall_max_diff,       max_diff);

        std::cout << prefix
                  << "  Max|anal|=" << max_anal
                  << "  Max|diff|=" << max_diff;
        if (max_anal > 0) {
            std::cout << "  rel="
            << std::setprecision(std::numeric_limits<double>::max_digits10) << std::defaultfloat
            << max_diff / max_anal;
        }
        std::cout << "\n";

        // CSV
        {
            std::ofstream csv(prefix + "_fd_check.csv");
            csv << "fd,analytical\n";
            // Set maximum precision for the CSV stream before the loop
            csv << std::setprecision(std::numeric_limits<double>::max_digits10)
                << std::defaultfloat;  // or std::fixed for fixed-point notation
            for (const auto& [fv, av] : pairs) csv << fv << "," << av << "\n";
            std::cout << "  written " << prefix << "_fd_check.csv\n";
        }

        // Gnuplot PNG — write a script file then invoke gnuplot (avoids SIGPIPE)
        {
            const std::string gp_file = prefix + "_fd_check.gp";
            std::ofstream gp(gp_file);
            gp << "set terminal pngcairo size 600,600\n"
               << "set output '" << prefix << "_fd_check.png'\n"
               << "set datafile separator ','\n"
               << "set xlabel 'Finite Difference Gradient'\n"
               << "set ylabel 'Analytical Gradient'\n"
               << "set title '" << prefix << "'\n"
               << "plot '" << prefix << "_fd_check.csv' skip 1 using 1:2 "
                  "with points pt 7 ps 0.8 notitle, "
                  "x lc rgb 'red' title 'y=x'\n";
        }
        if (std::system(("gnuplot " + prefix + "_fd_check.gp ").c_str()) == 0)
            std::cout << "  written " << prefix << "_fd_check.png\n";

        EXPECT_LT(max_diff / std::max(max_anal, 1.0), 1e-4) << "FD check failed for " << prefix;
    }

    std::cout << "Overall Max|anal|=" << overall_max_analytical
              << "  Max|diff|=" << overall_max_diff << "\n";
}

// a_matrix_to_relax: port of a_matrix_to_relax() in ke.R. Hand-computed reference for
// 1 pair, 2 terms, 2 internal eigenvalues, 1 overall mode.
//   a_int = [2, 3], lambda_int = [-3, -5], a_overall = [5], lambda_overall = [-1]
//   coef  = [0.5, 0.2], freq = [10, 20]
//   j=0: lp=-4, rowsum=0.5/116+0.2/416=0.00479111406, term=-5*-4*rowsum=0.09582228117, +=2*term
//   j=1: lp=-6, rowsum=0.5/136+0.2/436=0.00413518618, term=-5*-6*rowsum=0.12405558554, +=3*term
//   value = 2*0.09582228117 + 3*0.12405558554 = 0.56381131895
TEST(KEnRefTestSuite, AMatrixToRelaxKernel) {
    using Real = double;
    Eigen::MatrixX<Real> a_int(1, 2);            a_int << 2, 3;
    Eigen::RowVectorX<Real> lambda_int(2);       lambda_int << -3, -5;
    Eigen::MatrixX<Real> a_overall(1, 1);        a_overall << 5;
    Eigen::RowVectorX<Real> lambda_overall(1);   lambda_overall << -1;
    NamedMatrix<Real> coef(1, 2);                coef << 0.5, 0.2;
    NamedMatrix<Real> freq(1, 2);                freq << 10, 20;
    const SpecDenTermArray<Real> sdt{coef, freq};

    auto [value, grad] = KEnRef<Real>::a_matrix_to_relax(a_int, lambda_int, a_overall, lambda_overall, sdt, true);
    ASSERT_EQ(value.size(), 1);
    EXPECT_NEAR(value(0), 0.5638113189451191, 1e-9);
    ASSERT_TRUE(grad.has_value());
    EXPECT_EQ(grad->rows(), 1);
    EXPECT_EQ(grad->cols(), 2);
    EXPECT_NEAR((*grad)(0, 0), 0.0958222811671087, 1e-9);   // d value / d a_int[:,0]
    EXPECT_NEAR((*grad)(0, 1), 0.1240555855369672, 1e-9);   // d value / d a_int[:,1]

    // value-only path returns no gradient
    auto [value2, grad2] = KEnRef<Real>::a_matrix_to_relax(a_int, lambda_int, a_overall, lambda_overall, sdt, false);
    EXPECT_NEAR(value2(0), 0.5638113189451191, 1e-9);
    EXPECT_FALSE(grad2.has_value());
}

// ---------------------------------------------------------------------------
// Profiling harness for the vector r_array_to_d_array (gradient path) — DISABLED
// by default. Run with:
//   for B in 100000000 64 128 256 512; do \
//     KENREF_DARRAY_BLOCK=$B ./Google_tests_kenref_core_exe \
//       --gtest_also_run_disabled_tests \
//       --gtest_filter='*RArrayToDArrayBlockThreadSweep*'; done
// The block size is read once per process, so sweep it across process invocations
// (B >= N == one block per model == model-axis-only parallelism: the A/B baseline).
// Within a process it sweeps OpenMP thread counts and prints ms/call + speedup.
// ---------------------------------------------------------------------------
TEST(KEnRefBench, DISABLED_RArrayToDArrayBlockThreadSweep) {
    using R = double;
    const char* blkEnv = std::getenv("KENREF_DARRAY_BLOCK");
    const std::string blkStr = blkEnv ? blkEnv : "256 (default)";
    // KENREF_DARRAY_REUSE=1 -> reuse persistent ret1/ret2 buffers via r_array_to_d_array_into (no
    // per-call allocation/first-touch); default -> fresh-allocating returning overload.
    const bool reuse = std::getenv("KENREF_DARRAY_REUSE") != nullptr;

    struct Cfg { int M; long N; int iters; };
    const std::vector<Cfg> cfgs = {{3, 2000, 200}, {3, 50000, 50}, {8, 50000, 30}, {3, 300000, 20}};
    const std::vector<int> threads = benchThreadList({1, 2, 4, 8, 16, 0}); // 0 == "all" (project convention)

    std::printf("\n# r_array_to_d_array (gradient)  |  KENREF_DARRAY_BLOCK = %s  |  buffers = %s\n",
                blkStr.c_str(), reuse ? "REUSED (_into)" : "fresh-alloc");
    std::printf("# %-7s %-8s %-8s %-12s %-9s\n", "models", "pairs", "threads", "ms/call", "vs 1-thread");

    for (const auto& cfg : cfgs) {
        std::mt19937 gen(20260617u);
        std::uniform_real_distribution<R> d(0.15e-9, 0.6e-9); // internuclear vectors, ~Å in meters
        std::vector<CoordsMatrixType<R>> models(cfg.M);
        for (int m = 0; m < cfg.M; ++m)
            models[m] = CoordsMatrixType<R>::NullaryExpr(cfg.N, 3, [&]() { return d(gen); });

        // persistent buffers for the reuse path (allocated once here, refilled in place each call)
        std::vector<Eigen::Matrix<R, Eigen::Dynamic, 5>> ret1;
        std::vector<Eigen::Matrix<R, Eigen::Dynamic, 15>> ret2;

        double baseMs = 0.0;
        for (int t : threads) {
            volatile R sink = 0;
            { // warm-up (allocation / first-touch not timed)
                if (reuse) KEnRef<R>::r_array_to_d_array_into(models, false, true, false, t, ret1, ret2);
                else { auto w = KEnRef<R>::r_array_to_d_array(models, false, true, false, t); ret1 = std::move(std::get<0>(w)); ret2 = std::move(std::get<1>(w)); }
                sink += ret1[0](0, 0);
            }
            const auto t0 = std::chrono::steady_clock::now();
            for (int it = 0; it < cfg.iters; ++it) {
                if (reuse) {
                    KEnRef<R>::r_array_to_d_array_into(models, false, true, false, t, ret1, ret2);
                    sink += ret1[0](0, 0) + ret2[0](0, 0);
                } else {
                    auto r = KEnRef<R>::r_array_to_d_array(models, false, true, false, t);
                    sink += std::get<0>(r)[0](0, 0) + std::get<1>(r)[0](0, 0); // defeat DCE
                }
            }
            const auto t1 = std::chrono::steady_clock::now();
            const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / cfg.iters;
            if (t == 1) baseMs = ms;
            const std::string tlabel = (t == 0) ? "all" : std::to_string(t);
            std::printf("  %-7d %-8ld %-8s %-12.3f %-9.2f\n",
                        cfg.M, cfg.N, tlabel.c_str(), ms, baseMs / ms);
            (void) sink;
        }
    }
    SUCCEED();
}

// ---------------------------------------------------------------------------
// End-to-end benchmark of the relax-energy pathway (gradient), DISABLED by
// default. Unlike testCoordArrayToRelaxEnergy, all fixtures (PDB atom mapping,
// spec_den_relax_data, coordinates) are loaded ONCE up front so the one-time
// std::regex name-normalization / file parsing does NOT pollute the timing or a
// perf profile of this test. The timed loop averages many compute-only calls and
// sweeps OpenMP thread counts (incl. 16). Run with:
//   ./Google_tests_kenref_core_exe --gtest_also_run_disabled_tests \
//     --gtest_filter='*RelaxEnergyThreadSweep*'
// (must run from a build's google_tests/ dir so RELAX_DIR resolves).
// ---------------------------------------------------------------------------
TEST(KEnRefBench, DISABLED_RelaxEnergyThreadSweep) {
    using R = double;
    constexpr R k = 1, n = 1;
    constexpr const char* RELAX_DIR = "../../res/google_tests/relax/";

    // ---- load all data ONCE (untimed) ----
    const NamedRowVector<R> rates = Table(
        {{"5.0e+08", "1.0e+12", "1.0e+04", "2.5e+08", "2.5e+08", "2.5e+08"}},
        {{"kens", "kmethyl", "karo", "Dx", "Dy", "Dz"}}
    ).toNamedRowVector<R>();
    const auto atomNameMapping_to1 = IoUtils::getAtomMappingFromPdb<std::string, int>(
        GB3_PROTON_FILENAME, IoUtils::fill_atomId_to_index_Map);
    const bool handleNames = IoUtils::should_handleNames(atomNameMapping_to1);
    const auto atomNameMapping = get_atomNameMapping<R>(atomNameMapping_to1, handleNames);
    auto relax_data_list = IoUtils::load_spec_den_relax_data(RELAX_DIR, handleNames);
    ASSERT_EQ(relax_data_list.size(), 6u);
    // Prime the atom-id-pairs cache ONCE, exactly as coord_array_to_sigma_energy relies on (the GMX
    // force provider does this for the SIGMA model, KEnRefForceProvider.cpp ~691). Relax is not wired
    // into the force provider yet, so we prime it here; otherwise coord_array_to_relax_energy re-resolves
    // atom names via std::map on every call (the ~22% hotspot seen when the cache is empty). When relax
    // is wired into the force provider, populate the cache there in a relax branch like the sigma one.
    for (auto& rd : relax_data_list)
        rd.set_atomIdPairs_to_sub0Atom_id_pairs_cache(
            {KEnRef<R>::atomNamePairs_2_atomIdPairs(rd.get_atom_pairs(), atomNameMapping)});
    const auto coord_master = getAllModels_allAtomCoordsMatrix<R>(GB3_PROTON_FILENAME, {0, 2, 4});

    const std::vector<int> threads = benchThreadList({1, 2, 4, 8, 16, 0}); // 0 == "all"
    const int iters = 300;

    std::printf("\n# coord_array_to_relax_energy (gradient)  |  data loaded once, %d-iter average\n", iters);
    std::printf("# %-8s %-12s %-11s\n", "threads", "ms/call", "vs 1-thread");

    double baseMs = 0.0;
    for (int t : threads) {
        // coord_array_to_relax_energy scales coords in place (Å->m), so each call needs a fresh
        // Å copy; the copy is a few hundred atoms x 3 models (~tens of KB) and is included in the
        // timing (it is part of feeding the kernel, not parsing).
        { auto cc = coord_master; // warm-up (not timed)
          volatile R s = std::get<0>(KEnRef<R>::coord_array_to_relax_energy(
              cc, rates, relax_data_list, k, n, atomNameMapping, true, t)); (void) s; }

        volatile R sink = 0;
        const auto t0 = std::chrono::steady_clock::now();
        for (int it = 0; it < iters; ++it) {
            auto cc = coord_master;
            auto res = KEnRef<R>::coord_array_to_relax_energy(
                cc, rates, relax_data_list, k, n, atomNameMapping, true, t);
            sink += std::get<0>(res);
        }
        const auto t1 = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / iters;
        if (t == 1) baseMs = ms;
        const std::string tlabel = (t == 0) ? "all" : std::to_string(t);
        std::printf("  %-8s %-12.4f %-11.2f\n", tlabel.c_str(), ms, baseMs / ms);
        (void) sink;
    }
    SUCCEED();
}
