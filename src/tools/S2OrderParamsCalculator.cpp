#include "tools/AbstractCalculator.h"

/*
 * s2calc — reports S2 order parameters of a (multi-)trajectory, frame by frame.
 *
 * S2 is a metric, NOT an energy model: it does not run the KEnRefDriver. AbstractCalculator still
 * builds the selected EnergyModel (only to source the restrained atom pairs via finalizeIndexing ->
 * atomIdPairs, captured in subAtomIdPairs); this tool fits each model's sub-atoms to the reference and
 * calls KEnRef::s2OrderParams at the tool level. Offline trajectories are pre-PBC-corrected, so the fit
 * is Kabsch-only (no no-jump), matching the pre-restructure tool exactly.
 */
Eigen::IOFormat insideCsvLineFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "", "", "", "", "");

class S2OrderParamsCalculator : public AbstractCalculator {
public:
    S2OrderParamsCalculator() : AbstractCalculator("S2 Order Parameters calculator", "s2") {}

    void addSpecificParameters(CLI::App& /*app*/) override {}

    // s2OrderParams historically ran on all available threads.
    [[nodiscard]] int kernelNumOmpThreads() const override { return 0; }

    void prepareConsumer(std::unique_ptr<kenref::EnergyModel<KEnRef_Real_t>> /*model*/) override {
        // Not an energy model: subAtomIdPairs was already captured by AbstractCalculator::calc().
    }

    void perFrameReport() override {
        if (adapter->currentStep() % dt != 0)
            return;
        const int numModels = adapter->numModelsInThisProcess();
        std::vector<CoordsMatrixType<KEnRef_Real_t> > allModels(numModels);
        CoordsMatrixType<KEnRef_Real_t> guideX, subX;
        Eigen::Matrix<KEnRef_Real_t, 3, 3> box;
        for (int i = 0; i < numModels; ++i) {
            adapter->getLocalModelX(i, guideX, subX, box);
            const auto affine = Kabsch_Umeyama<KEnRef_Real_t>::find3DAffineTransform(
                guideX, guideAtomsReferenceCoordsCentered, false, false, false);
            allModels[i] = Kabsch_Umeyama<KEnRef_Real_t>::applyTransform(affine, subX);
        }
        const auto frameS2OrderParams = KEnRef<KEnRef_Real_t>::s2OrderParams(allModels, subAtomIdPairs, 0);
        if (adapter->currentFrame() == 0) {
            std::cout << "First S2OrderParams of first step\n"
                      << frameS2OrderParams.topRows(25).transpose() << std::endl;
        }
        out_file_stream << adapter->currentStep() << ", "
                        << frameS2OrderParams.transpose().format(insideCsvLineFormat) << std::endl;
    }
};

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    S2OrderParamsCalculator().calc(argc, argv);
}
