#ifndef KENREF_ABSTRACTCALCULATOR_H
#define KENREF_ABSTRACTCALCULATOR_H

#include <filesystem>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <regex>
#include <utility>
#include "CLI11/CLI11.hpp"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/gmxlib/network.h"
#include "core/kabsch.h"
#include "core/KEnRef.h"
#include "core/IoUtils.h"

struct t_file_state{
    t_fileio *xd;
    rvec     *x;
    matrix    box;
    int       nframe, natoms;
    int64_t   step;
    real      prec, time;
    gmx_bool  bOK;
    real     lambda;
};

inline void fillX(CoordsMatrixType<KEnRef_Real_t> &targetAtomsX, const std::vector<int> &idxes0, const rvec *x, const bool toAngstrom) {
    for (int i = 0; i < targetAtomsX.rows(); i++) {
        //        const int *piGlobal = new int{sub0Id_to_global1Id[i] - 1};
        //        const int *piLocal = forceProviderInput.cr_.dd->ga2la->findHome(*piGlobal);
        //        GMX_ASSERT(piLocal, "ERROR: Can't find local index of atom");
        const gmx::RVec atom_x = x[idxes0[i]];
#if VERBOSE
        std::cout << sub0Id_to_global1Id[i] << "\t" << *piGlobal << "\t" << *piLocal << "\t x: " << atom_x[0] << ", " << atom_x[1] << ", " << atom_x[2] << std::endl;
#endif
        if constexpr (std::is_same_v<KEnRef_Real_t, real>) {
            auto subAtomsX_buffer = targetAtomsX.data();
            const auto rvec = atom_x.as_vec();
            std::copy_n(rvec, 3, &subAtomsX_buffer[i * 3]);
        } else {
            for (int j = 0; j < 3; ++j) {
                targetAtomsX(i, j) = static_cast<KEnRef_Real_t>(atom_x[j]);
            }
        }
    }
    if (toAngstrom)
        targetAtomsX *= 10;
}


// Quick one-liners if you don't need a full class
inline std::string to_upper(std::string s) {  // must be copied
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    return s;
}

inline std::string capitalize(std::string s) { // must be copied
    if (!s.empty()) {
        s[0] = std::toupper(s[0]);
        std::transform(s.begin() + 1, s.end(), s.begin() + 1, ::tolower);
    }
    return s;
}

class AbstractCalculator {
protected:
    std::string toolName;
    std::string metricToolCalculates;
    std::string toolParamPrefix, metricToolCalculatesCapitalized;
    CoordsMatrixType<KEnRef_Real_t> lastFrameSubAtomsX_; //Used only for proper NoJump algorithm
    CoordsMatrixType<KEnRef_Real_t> lastFrameGuideAtomsX_ZEROIndexed_; //Used only for proper NoJump algorithm
    enum class InputFileType { xtc, trr, UNKNOWN };
    InputFileType inputFileType = InputFileType::UNKNOWN;
    std::string outputPathName;

    bool debug = false;
    std::string GUIDE_C_ALPHA = "guideC-alpha";
    std::vector<std::string> inputFiles;
    std::string indexFileName = "index.ndx";
    std::string refFileName = "ref.pdb";
    KEnRef<KEnRef_Real_t>::energyModel selected_energy_model = KEnRef<KEnRef_Real_t>::energyModel::UNKNOWN;
    std::string experimentalDataFileName, experimentalDataFolder;
    int max_frame = -1;
    uint dt = 10;

    std::vector<CoordsMatrixType<KEnRef_Real_t> > allSimulationsSubAtomsXVector;
    int num_models = 0;
    std::vector<std::tuple<int, int> > atomIdPairs; // Zero based
    std::vector<std::tuple<int, int> > subAtomIdPairs;
    std::map<std::string, int> atomName_to_atomSub0Id_map;


    std::ofstream out_file_stream;


public:
    AbstractCalculator(const std::string& toolName, const std::string& metricToolCalculates) {
        this->toolName = toolName;
        this->metricToolCalculates = metricToolCalculates;
        metricToolCalculatesCapitalized = capitalize(metricToolCalculates);
        toolParamPrefix = to_upper(metricToolCalculates);
        outputPathName = metricToolCalculates+".out";
    }
    virtual ~AbstractCalculator() = default;

    //remember that the data must be PBC corrected (in every step)
    int calc(int argc, char** argv);

    void add_common_parameters(CLI::App& app);

    virtual void addSpecificParameters(CLI::App &app) = 0;

    void handle_sigma_energy_model(std::map<std::string, int> atomName_to_atomGlobalId_map, bool handleNames);

    virtual void fill_spec_den_data_vector(const std::string &spec_den_data_prefix, const Table &atomPairAndSigmaTable, std::vector<std::tuple<std::string, std::string>> iterationAtomPairs) = 0;
    virtual void handle_plateaus_energy_model() = 0;
    virtual void fill_atomName_to_atomSub0Id_map_if_needed(const std::map<std::string, int>& atomName_to_atomGlobalId_map, int maxAtomIdOfInterest, std::vector<int> global0Id_to_sub0Id) = 0;
    virtual void calculate_and_report(std::vector<t_file_state>::value_type &fst) = 0;

};



#endif //KENREF_ABSTRACTCALCULATOR_H
