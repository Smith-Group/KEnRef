#include <limits>

#include "tools/AbstractCalculator.h"
#include "core/KEnRefDriver.h"

/*
 * energycalc — reports the KEnRef restraint energy of a (multi-)trajectory, frame by frame.
 *
 * After the restructure it owns NO per-model logic: AbstractCalculator builds the selected EnergyModel
 * from the registry and the sub-atom indexing; this tool only constructs a KEnRefDriver around the
 * prepared model and runs it once per reported frame (energy only — addLocalModelDerivatives is a no-op
 * in the offline adapter, so the gradient the driver computes is simply not applied).
 */
class EnergyCalculator : public AbstractCalculator {
    KEnRef_Real_t k = 0, n = 0;
    KEnRef_Real_t proton_mhz = 700.0;
    std::string ratesFile;
    std::unique_ptr<kenref::KEnRefDriver<KEnRef_Real_t>> driver_;

public:
    EnergyCalculator() : AbstractCalculator("energy calculator", "energy") {}

    void addSpecificParameters(CLI::App& app) override {
        app.add_option("-k,--k", k, "K force constant")->envname(toolParamPrefix + "_K");
        app.add_option("-n,--n", n, "power scaling")->envname(toolParamPrefix + "_N");
        app.add_option("-z,--proton-mhz", proton_mhz, "spectrometer proton field strength in MHz")
                ->envname(toolParamPrefix + "_PROTON_MHZ");
        app.add_option("--rates-file", ratesFile, "CSV file overriding the default SIGMA relaxation rates")
                ->envname(toolParamPrefix + "_RATES_FILE");
    }

    void addModelParams(std::map<std::string, std::string>& params) override {
        params["PROTON_MHZ"] = std::to_string(proton_mhz);
        if (!ratesFile.empty())
            params["RATES_FILE"] = ratesFile;
    }

    void prepareConsumer(std::unique_ptr<kenref::EnergyModel<KEnRef_Real_t>> model) override {
        // Energy is saturation-independent (saturation only clamps the derivatives, which we never
        // apply), so no max-force is needed: an infinite threshold leaves the energy untouched.
        constexpr KEnRef_Real_t noSaturation = std::numeric_limits<KEnRef_Real_t>::infinity();
        driver_ = std::make_unique<kenref::KEnRefDriver<KEnRef_Real_t>>(
            std::move(model), k, n, noSaturation, guideAtomsReferenceCoordsCentered);
    }

    void perFrameReport() override {
        if (adapter->currentStep() % dt == 0) {
            const KEnRef_Real_t energy = driver_->step(*adapter);
            out_file_stream << adapter->currentStep() << '\t' << energy << std::endl;
        }
    }
};

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    EnergyCalculator().calc(argc, argv);
}
