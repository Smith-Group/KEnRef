/*
 * TrajectoryEngineAdapter.cpp
 *
 *  Offline EngineAdapter implementation — see TrajectoryEngineAdapter.h. The trr/xtc frame stepping is
 *  ported verbatim from the old AbstractCalculator::calc() open/read loop (now owned here, once).
 */
#include "gmxinterface/TrajectoryEngineAdapter.h"

#include <filesystem>
#include <stdexcept>

#include "gromacs/utility/smalloc.h"

TrajectoryEngineAdapter::TrajectoryEngineAdapter(std::vector<std::string> inputFiles,
                                                 std::vector<int> guideAtom0Indices,
                                                 std::map<std::string, std::string> params,
                                                 int numOmpThreads)
    : inputFiles_(std::move(inputFiles)), guideAtom0Indices_(std::move(guideAtom0Indices)),
      params_(std::move(params)), numOmpThreads_(numOmpThreads) {
    const std::string ext = std::filesystem::path(inputFiles_.front()).extension().string();
    inputFileType_ = ext == ".trr" ? InputFileType::trr
                   : ext == ".xtc" ? InputFileType::xtc
                                   : InputFileType::UNKNOWN;
    models_.resize(inputFiles_.size());
}

TrajectoryEngineAdapter::~TrajectoryEngineAdapter() {
    closeAll();
}

std::optional<std::string> TrajectoryEngineAdapter::getRawParam(const std::string& key) const {
    const auto it = params_.find(key);
    if (it == params_.end() || it->second.empty())
        return std::nullopt;
    return it->second;
}

void TrajectoryEngineAdapter::fillFromFrame(CoordsMatrixType<KEnRef_Real_t>& target,
                                            const std::vector<int>& idxes0, const ModelFile& mf) {
    target.resize(static_cast<Eigen::Index>(idxes0.size()), 3);
    for (int i = 0; i < target.rows(); i++) {
        const gmx::RVec atom_x = mf.x[idxes0[i]];
        if constexpr (std::is_same_v<KEnRef_Real_t, real>) {
            const auto rvec = atom_x.as_vec();
            std::copy_n(rvec, 3, &target.data()[i * 3]);
        } else {
            for (int j = 0; j < 3; ++j)
                target(i, j) = static_cast<KEnRef_Real_t>(atom_x[j]);
        }
    }
    target *= 10; // nm -> Angstrom (matches the live engines' getLocalModelX)
}

void TrajectoryEngineAdapter::getLocalModelX(int localModel, CoordsMatrixType<KEnRef_Real_t>& guideX,
                                             CoordsMatrixType<KEnRef_Real_t>& subX,
                                             Eigen::Matrix<KEnRef_Real_t, 3, 3>& box) const {
    const ModelFile& mf = models_.at(localModel);
    fillFromFrame(guideX, guideAtom0Indices_, mf);
    fillFromFrame(subX,   subAtoms0Ids_,      mf);
    // GROMACS box (nm) -> Eigen 3x3 (raw; the driver scales it to Angstrom via toAngstrom).
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            box(i, j) = static_cast<KEnRef_Real_t>(mf.box[i][j]);
}

void TrajectoryEngineAdapter::openAll() {
    for (size_t m = 0; m < inputFiles_.size(); ++m) {
        ModelFile& mf = models_[m];
        const std::string& path = inputFiles_[m];
        switch (inputFileType_) {
            case InputFileType::xtc:
                mf.xd = open_xtc(path, "r");
                read_first_xtc(mf.xd, &mf.natoms, &mf.step, &mf.time, mf.box, &mf.x, &mf.prec, &mf.bOK);
                break;
            case InputFileType::trr: {
                gmx_trr_header_t header;
                gmx_trr_read_single_header(path, &header);
                snew(mf.x, header.x_size);
                mf.xd  = gmx_trr_open(path, "r");
                mf.bOK = gmx_trr_read_frame(mf.xd, &mf.step, &mf.time, &mf.lambda, mf.box,
                                            &mf.natoms, mf.x, nullptr, nullptr);
                break;
            }
            default:
                throw std::runtime_error("TrajectoryEngineAdapter: unrecognized input file type for '" + path + "'");
        }
    }
    nframe_ = 0;
}

bool TrajectoryEngineAdapter::advanceAll() {
    bool allOk = true;
    for (ModelFile& mf : models_) {
        int ret = 0;
        switch (inputFileType_) {
            case InputFileType::xtc:
                ret = read_next_xtc(mf.xd, mf.natoms, &mf.step, &mf.time, mf.box, mf.x, &mf.prec, &mf.bOK);
                break;
            case InputFileType::trr:
                ret = gmx_trr_read_frame(mf.xd, &mf.step, &mf.time, &mf.lambda, mf.box,
                                         &mf.natoms, mf.x, nullptr, nullptr);
                mf.bOK = ret; // not exactly correct, but it is fine (mirrors the old loop)
                break;
            default:
                throw std::runtime_error("TrajectoryEngineAdapter: unrecognized input file type");
        }
        if (ret == 0)
            allOk = false;
    }
    if (allOk)
        ++nframe_;
    return allOk;
}

void TrajectoryEngineAdapter::closeAll() {
    for (ModelFile& mf : models_) {
        if (mf.x) {
            sfree(mf.x);
            mf.x = nullptr;
        }
        if (mf.xd) {
            switch (inputFileType_) {
                case InputFileType::xtc: close_xtc(mf.xd); break;
                case InputFileType::trr: gmx_trr_close(mf.xd); break;
                default: break;
            }
            mf.xd = nullptr;
        }
    }
}
