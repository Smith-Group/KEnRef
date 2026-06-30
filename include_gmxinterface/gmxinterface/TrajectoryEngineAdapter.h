/*
 * TrajectoryEngineAdapter.h
 *
 *  ROLE: the offline (trajectory-file) implementation of kenref::EngineAdapter, letting the
 *  command-line tools (energycalc / s2calc) drive the SAME KEnRefDriver + model registry the live MD
 *  engines use, instead of their own duplicated open/fit/dispatch pipeline.
 *
 *  One process owns ALL N models (one per input trajectory) — there is no MPI, so the cross-process
 *  gather/scatter are trivial moves and there are no forces to apply (energy/metric only). It
 *  encapsulates the trr/xtc frame stepping that AbstractCalculator used to inline: openAll() reads the
 *  first frame of every trajectory, advanceAll() steps them in lock-step, and getLocalModelX()
 *  presents the current frame's guide + sub atoms (Angstrom) and box (raw nm) exactly as the GROMACS
 *  force provider does — so the driver's no-jump + Kabsch fit run byte-for-byte the same math.
 *
 *  No-jump note: the driver runs restoreNoJump (skipped on the first frame). Offline trajectories are
 *  expected to be PBC-corrected already (the tools always assumed this), so it is effectively a no-op.
 */
#ifndef KENREF_TRAJECTORYENGINEADAPTER_H
#define KENREF_TRAJECTORYENGINEADAPTER_H

#include <map>
#include <optional>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "gromacs/fileio/xtcio.h"
#include "gromacs/fileio/trrio.h"

#include "core/DefaultEngineAdapter.h"
#include "core/KEnRef.h"   // CoordsMatrixType<>, KEnRef_Real_t

class TrajectoryEngineAdapter final : public kenref::DefaultEngineAdapter<KEnRef_Real_t> {
public:
    enum class InputFileType { xtc, trr, UNKNOWN };

    // params: model-tier keyword -> raw string value (EXP_DATA_FOLDER, EXP_DATA_FILE, PROTON_MHZ,
    // RATES_FILE, ...) — backs getRawParam() so the registry model reads its config the same way the
    // live engines do. guideAtom0Indices: 0-based. numOmpThreads is forwarded to every kernel
    // (energycalc historically ran serial => default 1, preserving its bit-for-bit output).
    TrajectoryEngineAdapter(std::vector<std::string> inputFiles,
                            std::vector<int> guideAtom0Indices,
                            std::map<std::string, std::string> params,
                            int numOmpThreads = 1);
    ~TrajectoryEngineAdapter() override;

    // The compact sub-atom set (0-based global ids) the selected model restrains; set by the tool
    // after it derives the indexing from the model's atom-name pairs (mirrors the force provider).
    void setSubAtoms0Ids(std::vector<int> subAtoms0Ids) { subAtoms0Ids_ = std::move(subAtoms0Ids); }

    // ---- trajectory lifecycle -------------------------------------------------------------------
    void openAll();      // open every file and read its first frame; throws on failure
    bool advanceAll();   // read the next frame of every file; false once any file is exhausted
    void closeAll();
    [[nodiscard]] int64_t currentStep()  const { return models_.empty() ? 0 : models_.front().step; }
    [[nodiscard]] int     currentFrame() const { return nframe_; }
    [[nodiscard]] InputFileType inputFileType() const { return inputFileType_; }

    // ---- kenref::EngineAdapter ------------------------------------------------------------------
    [[nodiscard]] std::optional<std::string> getRawParam(const std::string& key) const override;
    [[nodiscard]] int numOmpThreads() const override { return numOmpThreads_; }
    [[nodiscard]] int numModelsInThisProcess() const override { return static_cast<int>(models_.size()); }
    void getLocalModelX(int localModel, CoordsMatrixType<KEnRef_Real_t>& guideX,
                        CoordsMatrixType<KEnRef_Real_t>& subX,
                        Eigen::Matrix<KEnRef_Real_t, 3, 3>& box) const override;
    void addLocalModelDerivatives(int /*localModel*/,
                                  const CoordsMatrixType<KEnRef_Real_t>& /*derivs*/) override {
        // energy/metric-only: the offline tools never apply forces.
    }
    [[nodiscard]] int numModelsTotal()  const override { return static_cast<int>(models_.size()); }
    // simulationIndex() (=0), gatherFittedSubAtomsX / scatterModelDerivatives (trivial moves) and the
    // numModelsInThisProcess/Total split: this single offline process owns the whole ensemble, so the
    // DefaultEngineAdapter defaults apply (numModelsInThisProcess/Total are overridden above to N).

private:
    // Per-trajectory open file + current frame (the old AbstractCalculator::t_file_state, owned here).
    struct ModelFile {
        t_fileio* xd = nullptr;
        rvec*     x  = nullptr;
        matrix    box{};
        int       natoms = 0;
        int64_t   step   = 0;
        real      prec = 0, time = 0, lambda = 0;
        gmx_bool  bOK = false;
    };

    // Fill `target` (resized [idxes0.size() x 3]) from frame `mf`'s coords at the given 0-based atom
    // indices, converting nm -> Angstrom (matching the live engines' getLocalModelX).
    static void fillFromFrame(CoordsMatrixType<KEnRef_Real_t>& target,
                              const std::vector<int>& idxes0, const ModelFile& mf);

    std::vector<std::string> inputFiles_;
    std::vector<int>         guideAtom0Indices_;
    std::vector<int>         subAtoms0Ids_;
    std::map<std::string, std::string> params_;
    int                      numOmpThreads_;
    InputFileType            inputFileType_ = InputFileType::UNKNOWN;
    std::vector<ModelFile>   models_;
    int                      nframe_ = 0;
};

#endif // KENREF_TRAJECTORYENGINEADAPTER_H
