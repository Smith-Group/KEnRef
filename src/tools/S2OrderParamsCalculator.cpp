#include "tools/AbstractCalculator.h"

Eigen::IOFormat insideCsvLineFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "", "", "", "", "");

class S2OrderParamsCalculator : public AbstractCalculator{

public:
    S2OrderParamsCalculator() : AbstractCalculator("S2 Order Parameters calculator", "s2") { }

    void fill_spec_den_data_vector(const std::string &spec_den_data_prefix, const Table &atomPairAndSigmaTable, std::vector<std::tuple<std::string, std::string>> iterationAtomPairs) override { }

    void addSpecificParameters(CLI::App &app) override { }

    void handle_plateaus_energy_model() override {
        auto experimentalData_table = IoUtils::readTable(experimentalDataFileName, true, false, ",");
        for (int i = 0; i < experimentalData_table.rowCount(); ++i) {
            std::istringstream temp1(experimentalData_table(i, "i1")), temp2(experimentalData_table(i, "i2")), temp3(experimentalData_table(i, "g1")), temp4(experimentalData_table(i, "g2"));
            int i1, i2;
            temp1 >> i1;
            temp2 >> i2;
            atomIdPairs.emplace_back(i1 - 1, i2 - 1);
        }
    }

    void fill_atomName_to_atomSub0Id_map_if_needed(const std::map<std::string, int>& atomName_to_atomGlobalId_map, int maxAtomIdOfInterest, std::vector<int> global0Id_to_sub0Id) override {
    }

    void calculate_and_report(std::vector<t_file_state>::value_type &fst) override {
        if (fst.step % dt == 0) {
            const auto &frameS2OrderParams = KEnRef<KEnRef_Real_t>::s2OrderParams( allSimulationsSubAtomsXVector, subAtomIdPairs, 0);
            if (fst.nframe == 0) {
                std::cout << "First S2OrderParams of first step\n" << frameS2OrderParams.topRows(25).transpose() << std::endl;
            }
            out_file_stream << fst.step << ", " << frameS2OrderParams.transpose().format(insideCsvLineFormat) << std::endl;
        }
    }
};

int main(int argc, char** argv){
    MPI_Init(&argc, &argv);
    S2OrderParamsCalculator().calc(argc, argv);
}

