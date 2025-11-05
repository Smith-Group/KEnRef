#include "tools/AbstractCalculator.h"


class EnergyCalculator : public AbstractCalculator{
    const NamedRowVector<KEnRef_Real_t> rates = Table(		//TODO provide a way to change it
        {{"5.0e+08", "2.5e+08", "1.0e+12", "1.0e+04"}},
        {{"kens", "kc", "kmethyl", "karo"}}
        ).toNamedRowVector<KEnRef_Real_t>();
    KEnRef_Real_t k=0, n=0;
    KEnRef_Real_t proton_mhz = 700.0;
    std::vector<std::vector<std::vector<int> > > simulated_grouping_list;
    std::vector<SpecDenData<KEnRef_Real_t> > spec_den_data_vector;
    Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, Eigen::Dynamic> g0;

public:
    EnergyCalculator() : AbstractCalculator("energy calculator", "energy") {
    }
    void fill_spec_den_data_vector(const std::string &spec_den_data_prefix, const Table &atomPairAndSigmaTable, std::vector<std::tuple<std::string, std::string>> iterationAtomPairs) override {
        // sigma
        std::vector<KEnRef_Real_t> sigmasVec{};
        for (int row = 0; row < atomPairAndSigmaTable.rowCount(); ++row) {
            if (!atomPairAndSigmaTable.isRowComplete(row))
                break;
            const auto &valueStr = atomPairAndSigmaTable.at(row, "sigma");
            std::istringstream iss(valueStr);
            KEnRef_Real_t value;
            iss >> value;
            sigmasVec.emplace_back(value);
        }
        std::optional<NamedVector<KEnRef_Real_t> > sigma(sigmasVec.size());
        for (int j = 0; j < sigma->rows(); ++j) {
            sigma.value()(j, 0) = sigmasVec[j];
        }
        //multiple_grouping
        std::string multiple_grouping_fileName = experimentalDataFolder + spec_den_data_prefix + "_groupings.csv";
        const auto &grouping_matrix = NamedMatrix<int>(IoUtils::readTable(
                                                           multiple_grouping_fileName, false, false).toNamedMatrix<int>().array() - 1);
        const auto &multiple_grouping = IoUtils::grouping_mat_to_subset_idx(grouping_matrix);
        //a_coef
        std::string aCoefFileName = experimentalDataFolder + spec_den_data_prefix + "_a_coef.csv";
        const auto &a_coef = IoUtils::readTable(aCoefFileName, true, false, "\\s*,\\s*", -1, false).
                toNamedMatrix<KEnRef_Real_t>();
        //lambda_coef
        std::string lambdaCoefFileName =
                experimentalDataFolder + spec_den_data_prefix + "_lambda_coef.csv";
        const auto &lambda_coef = IoUtils::readTable(lambdaCoefFileName, true, true, "\\s*,\\s*", -1, false).
                toNamedMatrix<KEnRef_Real_t>();

        spec_den_data_vector.emplace_back(std::move(SpecDenData{
            iterationAtomPairs, sigma, multiple_grouping, a_coef, lambda_coef })); //you can leave it `emplace_back`, but if you want to use at(), then you have to use reserve()
    }

    void addSpecificParameters(CLI::App& app) override {
        app.add_option("-k,--k", k, "K force constant")->envname(toolParamPrefix + "_K");
        app.add_option("-n,--n", n, "power scaling")->envname(toolParamPrefix + "_N");
        app.add_option("-z,--proton-mhz", proton_mhz, "spectrometer proton field strength in MHz")
                ->envname(toolParamPrefix + "_PROTON_MHZ");
    }

    void handle_plateaus_energy_model() override {
        switch (num_models) {
            case 1:
                simulated_grouping_list = std::vector<std::vector<std::vector<int> > >{{{0}}, {{0}}};
                break;
            case 2:
                simulated_grouping_list = std::vector<std::vector<std::vector<int> > >{{{0, 1}}, {{0}, {1}}};
                break;
            case 3:
                simulated_grouping_list = std::vector<std::vector<std::vector<int> > > {{{0, 1, 2}}, {{0}, {1}, {2}}};
                break;
            default:
                assert(num_models <= 3 && "I don't know how to handle more than 3 simulations yet");
        }
        auto experimentalData_table = IoUtils::readTable(experimentalDataFileName, true, false, ",");
        g0.resize(static_cast<int>(experimentalData_table.rowCount()), 2);
        for (int i = 0; i < experimentalData_table.rowCount(); ++i) {
            std::istringstream temp1(experimentalData_table(i, "i1")), temp2(experimentalData_table(i, "i2")), temp3(experimentalData_table(i, "g1")), temp4(experimentalData_table(i, "g2"));
            int i1, i2;
            temp1 >> i1;
            temp2 >> i2;
            atomIdPairs.emplace_back(i1 - 1, i2 - 1);
            temp3 >> g0(i, 0);
            temp4 >> g0(i, 1);
        }
    }

    void fill_atomName_to_atomSub0Id_map_if_needed(const std::map<std::string, int>& atomName_to_atomGlobalId_map, int maxAtomIdOfInterest, std::vector<int> global0Id_to_sub0Id) override {
        for (const auto &[name, global1Id]: atomName_to_atomGlobalId_map){
            if (global1Id > maxAtomIdOfInterest + 1)
                continue;
            atomName_to_atomSub0Id_map[name] = global0Id_to_sub0Id.at(global1Id - 1);
        }
    }

    void calculate_and_report(std::vector<t_file_state>::value_type &fst) override {
        if (fst.step % dt == 0) {
            KEnRef_Real_t energy = 0;
            std::optional<std::vector<CoordsMatrixType<KEnRef_Real_t> > > allDerivatives_vector;
            if (selected_energy_model == KEnRef<double>::energyModel::PLATEAUS) {
                std::tie(energy, allDerivatives_vector) =
                        KEnRef<KEnRef_Real_t>::coord_array_to_g_energy(
                            allSimulationsSubAtomsXVector, subAtomIdPairs, simulated_grouping_list, g0, k, n, false, 1);
            }else if (selected_energy_model == KEnRef<double>::energyModel::SIGMA) {
                std::tie(energy, allDerivatives_vector) =
                        KEnRef<KEnRef_Real_t>::coord_array_to_sigma_energy(
                            allSimulationsSubAtomsXVector, rates, spec_den_data_vector, proton_mhz, k, n,
                            atomName_to_atomSub0Id_map, false, 1);
            }
            out_file_stream << fst.step << '\t' << energy << std::endl;
        }
    }
};

int main(int argc, char** argv){
    MPI_Init(&argc, &argv);
    EnergyCalculator().calc(argc, argv);
}

