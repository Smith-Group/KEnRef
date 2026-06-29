#ifndef KENREF_ABSTRACTCALCULATOR_H
#define KENREF_ABSTRACTCALCULATOR_H

#include <filesystem>
#include <iostream>
#include <iomanip>
#include <memory>
#include <unistd.h>
#include <regex>
#include <utility>
#include "CLI11/CLI11.hpp"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/gmxlib/network.h"
#include "core/kabsch.h"
#include "core/KEnRef.h"
#include "core/IoUtils.h"
#include "core/EnergyModel.h"
#include "core/ModelRegistry.h"
#include "gmxinterface/TrajectoryEngineAdapter.h"

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

/*
 * AbstractCalculator — the shared driver for the offline trajectory tools (energycalc / s2calc).
 *
 * After the model-abstraction restructure it no longer branches on a per-tool energyModel enum nor
 * inlines the trr/xtc open/fit pipeline. calc() now:
 *   1. reads the reference PDB atom-name -> global-id map,
 *   2. builds the selected EnergyModel from the registry (by name) and lets it load its inputs,
 *   3. derives the sub-atom indexing from the model's atom-name pairs (identical to the GROMACS force
 *      provider's fillParamsStep0), and finalizes the model's indexing,
 *   4. hands the prepared model to the concrete tool via prepareConsumer() (energycalc builds a
 *      KEnRefDriver; s2calc just keeps the sub-atom pairs), then
 *   5. steps every trajectory in lock-step through a TrajectoryEngineAdapter, calling perFrameReport().
 */
class AbstractCalculator {
protected:
    std::string toolName;
    std::string metricToolCalculates;
    std::string toolParamPrefix, metricToolCalculatesCapitalized;
    std::string outputPathName;

    bool debug = false;
    std::string GUIDE_C_ALPHA = "guideC-alpha";
    std::vector<std::string> inputFiles;
    std::string indexFileName = "index.ndx";
    std::string refFileName = "ref.pdb";
    std::string modelName;                       // registry model name (e.g. "SIGMA" / "PLATEAUS")
    std::string experimentalDataFileName, experimentalDataFolder;
    int max_frame = -1;
    uint dt = 10;

    int num_models = 0;
    // Compact (0-based sub-atom) pairs the selected model restrains; sourced from the model after
    // finalizeIndexing. Used by s2calc (and available to any tool) at the tool level.
    std::vector<std::tuple<int, int> > subAtomIdPairs;
    // Centered guide-atom reference coords (Angstrom) for Kabsch fitting.
    CoordsMatrixType<KEnRef_Real_t> guideAtomsReferenceCoordsCentered;

    std::unique_ptr<TrajectoryEngineAdapter> adapter;
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

    // ---- per-tool customization -----------------------------------------------------------------
    // Add tool-specific CLI options (energycalc: K/N/PROTON_MHZ; s2calc: none).
    virtual void addSpecificParameters(CLI::App &app) = 0;
    // Inject any extra model-tier params into the adapter's param map (energycalc: PROTON_MHZ/RATES_FILE).
    virtual void addModelParams(std::map<std::string, std::string>& /*params*/) {}
    // Number of OpenMP threads the adapter forwards to the kernels (energycalc historically ran serial).
    [[nodiscard]] virtual int kernelNumOmpThreads() const { return 1; }
    // Take ownership of the prepared model (energycalc builds a KEnRefDriver; s2calc ignores it, having
    // already captured subAtomIdPairs). Receives the centered guide reference for fitting.
    virtual void prepareConsumer(std::unique_ptr<kenref::EnergyModel<KEnRef_Real_t>> model) = 0;
    // Process the adapter's CURRENT frame (already advanced by calc()) and write any output.
    virtual void perFrameReport() = 0;
};

#endif //KENREF_ABSTRACTCALCULATOR_H
