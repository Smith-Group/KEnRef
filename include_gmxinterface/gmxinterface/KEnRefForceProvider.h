/*
 * KEnRefForceProvider.h
 *
 *      Author: amr
 *
 *  GROMACS consumer of kenref_core. After the model-abstraction restructure this is a thin
 *  EngineAdapter: it provides the per-step GROMACS I/O (atom coords, force application, the box, MPI
 *  gather/scatter, and parameter lookup from Settings) and delegates the whole per-step refinement
 *  pipeline to a kenref::KEnRefDriver that holds the user-selected EnergyModel. There is no per-model
 *  switch here any more — adding a model touches only kenref_core.
 */

#ifndef KENREFFORCEPROVIDER_H_
#define KENREFFORCEPROVIDER_H_

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <gromacs/mdtypes/iforceprovider.h>
#include <gromacs/mdrun/simulationcontext.h>
#include <gromacs/selection/selection.h>
#include <chrono>

#include "gmxwrapper.h"
#include "core/KEnRef.h"
#include "core/DefaultEngineAdapter.h"
#include "core/KEnRefDriver.h"


class KEnRefForceProvider final : public gmx::IForceProvider,
                                  public kenref::DefaultEngineAdapter<KEnRef_Real_t> {

private:
	//force the user to select one value
	std::string selectedEnergyModel; // registry model name (SIGMA / PLATEAUS / RELAX / ...)

	gmx::SimulationContext *simulationContext_ = nullptr;
	bool paramsInitialized = false; // We use it instead of (step == 0), because the simulation may be continuing after 0
	std::shared_ptr<std::vector<int> const> guideAtom0Indices_; //ZERO indexed
	std::shared_ptr<CoordsMatrixType<KEnRef_Real_t> const> guideAtomsReferenceCoords_; //ZERO indexed
	std::shared_ptr<CoordsMatrixType<KEnRef_Real_t> const> guideAtomsReferenceCoordsCentered_; //ZERO indexed. Cashed for faster Kabsch Algorithm
	std::shared_ptr<std::map<std::string, int> const> atomName_to_atomGlobalId_map_; //built from the mapping PDB, passed to buildModelIndexing
	std::shared_ptr<CoordsMatrixType<KEnRef_Real_t> > subAtomsX_; //per-step scratch for this replica's restrained atoms
	std::shared_ptr<std::vector<int> > sub0Id_to_global1Id_; //Global ID is ONE based, subId is a small subset and is ZERO based (from buildModelIndexing)
	// (atomName_to_atomSub0Id_map_ + global1Id_to_sub0Id_ removed: the sub-indexing now lives in
	//  kenref::buildModelIndexing, and neither was used at step time.)
	std::shared_ptr<CoordsMatrixType<KEnRef_Real_t> > allSimulationsSubAtomsX_; //allocated once, used every step (MPI gather buffer)
	std::shared_ptr<KEnRef_Real_t[]> allDerivatives_buffer_;
	std::shared_ptr<KEnRef_Real_t[]> derivatives_buffer_;
	long long calculateForces_time = 0;

	// The shared per-step pipeline + the selected energy model live here now.
	std::unique_ptr<kenref::KEnRefDriver<KEnRef_Real_t> > driver_;

	// ---- per-step state, stashed at the top of calculateForces so the EngineAdapter callbacks
	//      (invoked by driver_->step(*this)) can reach the current GROMACS input/output + MPI info ----
	const gmx::ForceProviderInput *currentInput_ = nullptr;
	gmx::ForceProviderOutput *currentOutput_ = nullptr;
	bool isMultiSimulation_ = false;
	int numSimulations_ = 1;
	int simulationIndex_ = 0;
	MPI_Comm mainRanksComm_ = MPI_COMM_NULL;

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	KEnRefForceProvider();

	~KEnRefForceProvider() override;

    [[maybe_unused]] [[maybe_unused]] KEnRefForceProvider(const KEnRefForceProvider &other) = delete;

    [[maybe_unused]] [[maybe_unused]] KEnRefForceProvider(KEnRefForceProvider &&other) noexcept = delete;

	KEnRefForceProvider &operator=(const KEnRefForceProvider &other) = delete;

	KEnRefForceProvider &operator=(KEnRefForceProvider &&other) noexcept = delete;

	void calculateForces(const gmx::ForceProviderInput &forceProviderInput,
	                     gmx::ForceProviderOutput *forceProviderOutput) override;

    [[maybe_unused]]void setSimulationContext(gmx::SimulationContext *simulationContext);

	void setGuideAtom0Indices(std::shared_ptr<std::vector<int> const> targetAtoms0Indices);

	void setGuideAtomsReferenceCoords(
		std::shared_ptr<const CoordsMatrixType<KEnRef_Real_t>> &guideAtomsReferenceCoords);

    void fillParamsStep0(size_t homenr, int numSimulations, const gmx::ForceProviderInput &forceProviderInput);

	static CoordsMatrixType<KEnRef_Real_t>
	getGuideAtomsX(const std::vector<int> &guideAtom0Indices,
                    const gmx::ForceProviderInput &forceProviderInput, bool toAngstrom);

    static void fillSubAtomsX(CoordsMatrixType<KEnRef_Real_t> &subAtomsX, const std::vector<int> &sub0Id_to_global1Id,
                              const gmx::ForceProviderInput &forceProviderInput, bool toAngstrom);

	void set_selected_energy_model(std::string selected_energy_model) {
		selectedEnergyModel = std::move(selected_energy_model);
	}

	// ----------------------------------------------------------------------------------------------
	//  EngineAdapter implementation (GROMACS, live MD: exactly one model in this process/replica)
	// ----------------------------------------------------------------------------------------------
	[[nodiscard]] std::optional<std::string> getRawParam(const std::string &key) const override;
	[[nodiscard]] int numOmpThreads() const override;
	// numModelsInThisProcess() == 1 (one replica per process) comes from DefaultEngineAdapter.
	void getLocalModelX(int localModel, CoordsMatrixType<KEnRef_Real_t> &guideX,
	                    CoordsMatrixType<KEnRef_Real_t> &subX, Eigen::Matrix<KEnRef_Real_t, 3, 3> &box) const override;
	void addLocalModelDerivatives(int localModel, const CoordsMatrixType<KEnRef_Real_t> &derivs) override;
	[[nodiscard]] int numModelsTotal() const override { return numSimulations_; }
	[[nodiscard]] int simulationIndex() const override { return simulationIndex_; }
	void gatherFittedSubAtomsX(const std::vector<CoordsMatrixType<KEnRef_Real_t>> &localFitted,
	                           std::vector<CoordsMatrixType<KEnRef_Real_t>> &all) const override;
	void scatterModelDerivatives(const std::vector<CoordsMatrixType<KEnRef_Real_t>> &allPerModel,
	                             std::vector<CoordsMatrixType<KEnRef_Real_t>> &localPerModel) const override;
};

#endif /* KENREFFORCEPROVIDER_H_ */
