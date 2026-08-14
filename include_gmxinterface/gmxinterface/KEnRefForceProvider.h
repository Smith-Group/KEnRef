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
#include <gromacs/domdec/localatomset.h>
#include <gromacs/domdec/localatomsetmanager.h>
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
	/*! \brief The same sub-atom list as sub0Id_to_global1Id_, but ZERO-based.
	 *
	 * GROMACS numbers atoms from 0 everywhere (ga2la, LocalAtomSet), while the mapping PDB gives
	 * 1-based serials. Converting once at setup keeps the "-1" out of the per-step loops and gives
	 * both the LocalAtomSet registration and the ga2la fallback a single source of truth. */
	std::vector<int> subGlobal0Ids_;
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
	/*! \brief Whether this rank is the main rank of ITS OWN simulation.
	 *
	 * True for every rank when there is one rank per replica, which is why nothing depended on it
	 * before. Under domain decomposition it distinguishes the one rank per replica that may use
	 * mainRanksComm_ (MPI_COMM_NULL elsewhere) from the ranks that may not. Refreshed each step from
	 * the intra-simulation communicator rather than cached at setup, because it costs one integer
	 * and cannot then go stale. */
	bool isSimMainRank_ = true;
	MPI_Comm mainRanksComm_ = MPI_COMM_NULL;

	/*! \brief Owns every LocalAtomSet; handed over at setup time, null if the notifier never fired.
	 *
	 * Null is a supported state, not an error: it means this build did not receive the
	 * simulation-setup notification, and the fill loops fall back to their ga2la path. */
	gmx::LocalAtomSetManager *localAtomSetManager_ = nullptr;
	/*! \brief The restrained (sub) atoms and the guide atoms, maintained by DD across repartitions.
	 *
	 * Registered once at setup. `localIndex()[k]` is the local index of the k-th atom this rank owns
	 * and `collectiveIndex()[k]` is that atom's slot in the set's own (global) ordering -- exactly the
	 * row index the fill loops need. Empty optionals when no manager was delivered. */
	std::optional<gmx::LocalAtomSet> subAtomSet_;
	std::optional<gmx::LocalAtomSet> guideAtomSet_;

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

	/*! \brief Hand over the manager that owns every LocalAtomSet, before any atom set is registered.
	 *
	 * Must be called before initParamsAtSetup(), which does the registering. Passing nullptr is
	 * allowed and simply leaves KEnRef on its ga2la fallback. */
	void setLocalAtomSetManager(gmx::LocalAtomSetManager *manager);

	void setGuideAtom0Indices(std::shared_ptr<std::vector<int> const> targetAtoms0Indices);

	void setGuideAtomsReferenceCoords(
		std::shared_ptr<const CoordsMatrixType<KEnRef_Real_t>> &guideAtomsReferenceCoords);

	/*! \brief One-time setup: atom mapping, model indexing, per-step buffers and the driver.
	 *
	 * Runs at SETUP time (from KEnRefMDModule::initForceProviders), not at step 0. That is a hard
	 * requirement rather than a tidy-up: a LocalAtomSet may only be registered before the first
	 * dd_partition_system(), because LocalAtomSetData's constructor initialises localIndex_ to the
	 * GLOBAL indices and only setIndicesInDomainDecomposition() (called from each repartition) makes
	 * them local. A set registered at step 0 would therefore hand out global indices as if they were
	 * local until the next repartition -- silently wrong under DD. The sub-atom index list comes from
	 * buildModelIndexing(), so that has to happen here too.
	 *
	 * Safe to call this early because it reads only Settings, the mapping files and the simulation
	 * context: the old signature took `homenr` and `ForceProviderInput` but used neither. */
	void initParamsAtSetup();

	/*! \brief Fill the rows of \p dest owned by THIS rank, leave the rest zero, then reduce.
	 *
	 * The zero-fill convention that makes this exact is described on kenrefSumOverSimRanks(). Two
	 * equivalent ways to decide ownership:
	 *  - \p atomSet present: iterate the set's localIndex()/collectiveIndex(), which DD maintains.
	 *    O(atoms this rank owns), and no ga2la lookups at all.
	 *  - \p atomSet absent: look every global index up in ga2la and skip the ones not home here.
	 *    O(all atoms in the set) per rank, and the only path available when no manager was delivered.
	 * Both write each row on exactly one rank, which KENREF_DD_SELFCHECK=1 verifies at runtime.
	 *
	 * \param[in]  global0Ids ZERO-based global index of each row, used only by the fallback path. */
	void fillOwnedRows(CoordsMatrixType<KEnRef_Real_t> &dest,
	                   const std::optional<gmx::LocalAtomSet> &atomSet,
	                   const std::vector<int> &global0Ids,
	                   const gmx::ForceProviderInput &forceProviderInput,
	                   const char *what) const;

	[[nodiscard]] CoordsMatrixType<KEnRef_Real_t>
	getGuideAtomsX(const gmx::ForceProviderInput &forceProviderInput, bool toAngstrom) const;

    void fillSubAtomsX(CoordsMatrixType<KEnRef_Real_t> &subAtomsX,
                       const gmx::ForceProviderInput &forceProviderInput, bool toAngstrom) const;

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
	/*! \brief Which model this process contributes; 0 identifies the MASTER, per EngineAdapter.
	 *
	 * The driver's only use of this is `if (adapter.simulationIndex() == 0)` -- the test for "am I the
	 * process that assembles the ensemble and computes the energy". Under domain decomposition EVERY
	 * rank of replica 0 has simulationIndex_ == 0, so the bare member would make all of them run the
	 * model and drive the cross-replica collectives. Reporting a non-master value on a non-main rank
	 * keeps that test meaning exactly what EngineAdapter documents it to mean. Unchanged at one rank
	 * per replica, where every rank is its simulation's main rank. */
	[[nodiscard]] int simulationIndex() const override { return isSimMainRank_ ? simulationIndex_ : -1; }

	/*! \brief Whether this is THE ONE rank in the whole run that speaks for the ensemble.
	 *
	 * True on exactly one rank of the entire job: the main rank of replica 0. That rank is the one
	 * that assembles every replica's coordinates, runs the energy model, and therefore holds the only
	 * meaningful energy value -- every other rank gets 0 back from the driver.
	 *
	 * BOTH halves are needed, and each fails differently on its own:
	 *  - `simulationIndex_ == 0` alone is true on EVERY rank of replica 0 under domain decomposition,
	 *    so per-domain duplicates appear (and, in the driver, the model would run once per domain).
	 *  - `isSimMainRank_` alone is true on the main rank of EVERY replica, so each replica reports.
	 * Keeping the pair in one named place stops the two copies drifting apart. */
	[[nodiscard]] bool isEnsembleMasterRank() const { return simulationIndex_ == 0 && isSimMainRank_; }
	void gatherFittedSubAtomsX(const std::vector<CoordsMatrixType<KEnRef_Real_t>> &localFitted,
	                           std::vector<CoordsMatrixType<KEnRef_Real_t>> &all) const override;
	void scatterModelDerivatives(const std::vector<CoordsMatrixType<KEnRef_Real_t>> &allPerModel,
	                             std::vector<CoordsMatrixType<KEnRef_Real_t>> &localPerModel) const override;
};

#endif /* KENREFFORCEPROVIDER_H_ */
