/*
 * KENRefMDModule.cpp
 *
 *      Author: amr
 */

#include <iostream>
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/topology/topology.h"
#include "gmxinterface/gmxkenrefinitializer.h"
#include "gmxinterface/KEnRefMDModule.h"

#include "gmxinterface/gmxwrapper.h"


KEnRefMDModule::KEnRefMDModule() {

	std::vector<int>const& indices = GmxKEnRefInitializer::loadGmxIndexGroup(Settings::GUIDE_C_ALPHA, Settings::indexFileName);
	this->guideAtoms0Indexed = std::make_shared<std::vector<int> const>(indices);
    std::cout << "Guide atoms indices (0 indexed):\n";
    IoUtils::printVector(indices);

    auto allAtomReferenceCoords = IoUtils::getAtomMappingFromPdb<int, Eigen::RowVector3<KEnRef_Real_t>>(
            Settings::refFileName,
            IoUtils::fill_atomIndex1_to_coords_Map<KEnRef_Real_t>);
    //Read and save guideAtomsReferenceCoords_
    CoordsMatrixType<KEnRef_Real_t> guideAtomsReferenceCoords(guideAtoms0Indexed->size(), 3);
    for (int index0 = 0; index0 < guideAtoms0Indexed->size(); ++index0) {
        guideAtomsReferenceCoords(index0, Eigen::placeholders::all) = allAtomReferenceCoords[guideAtoms0Indexed->at(index0) + 1];
    }
//    std::cout << "Reference Atoms"<< std::endl << guideAtomsReferenceCoords_ << std::endl;
//    guideAtomsReferenceCoords *= 10; // Don't multiply. they are already in Angstrom
    this->guideAtomsReferenceCoords_ = std::make_shared<const CoordsMatrixType<KEnRef_Real_t >>(guideAtomsReferenceCoords);

    //TODO: build another matrix for subAtomsXReferenceCoords

}

KEnRefMDModule::~KEnRefMDModule() = default;
[[maybe_unused]] KEnRefMDModule::KEnRefMDModule(const KEnRefMDModule &other) = default;
[[maybe_unused]] KEnRefMDModule::KEnRefMDModule(KEnRefMDModule &&other) noexcept = default;
KEnRefMDModule& KEnRefMDModule::operator=(const KEnRefMDModule &other) = default;
KEnRefMDModule& KEnRefMDModule::operator=(KEnRefMDModule &&other) noexcept = default;

//! Returns an interface for handling mdp input (and tpr I/O).
gmx::IMdpOptionProvider* KEnRefMDModule::mdpOptionProvider() {
    std::cout << "====> KEnRefMDModule::mdpOptionProvider() called" << std::endl;
	return &kEnRefOptions;
}

//! Returns an interface for handling output files during simulation.
gmx::IMDOutputProvider* KEnRefMDModule::outputProvider() {
    std::cout << "====> KEnRefMDModule::outputProvider() called" << std::endl;
	return &kEnRefOutputProvider;
}

//! Initializes force provider(s) from this module (only one now).
void KEnRefMDModule::initForceProviders(gmx::ForceProviders* forceProviders) {
	std::cout << "====> KEnRefMDModule::initForceProviders() called" << std::endl;
	forceProvider_ = std::make_shared<KEnRefForceProvider>();
	forceProvider_->setSimulationContext(simulationContext_);
    forceProvider_->setGuideAtom0Indices(this->guideAtoms0Indexed);
    forceProvider_->setGuideAtomsReferenceCoords(this->guideAtomsReferenceCoords_);
	forceProvider_->set_selected_energy_model(Settings::selected_energy_model);
forceProviders->addForceProvider(forceProvider_.get());
}

//! Subscribe to simulation setup notifications
//
// IMPORTANT: GROMACS never calls this for us on its own. `MDModules::add()` stores the module in
// `impl_->modules_`, and that vector is iterated in exactly ONE place -- `initForceProviders`. Every
// `subscribe*` dispatcher in mdmodules.cpp hard-codes the five built-in modules (densityFitting_,
// qmmm_, colvars_, plumed_, nnpot_), so an externally added module silently receives a subset of the
// IMDModule interface it implements. KEnRef therefore calls this itself from its own mdrun.cpp, on
// the notifier object returned by the public MDModules::notifiers() -- which is the same object the
// runner emits on (runner.cpp: `mdModules_->notifiers().simulationSetupNotifier_`), so registering
// here is enough to receive everything the built-ins receive.
void KEnRefMDModule::subscribeToSimulationSetupNotifications(gmx::MDModulesNotifiers* notifiers) {
	if (notifiers == nullptr)
		return;

	// The manager that owns every LocalAtomSet. GROMACS emits it once, after the domain decomposition
	// has been built but before the first repartition (runner.cpp: `setupNotifier.notify(&atomSets)`),
	// which is exactly when a set may still be registered. Sets added here are then kept up to date by
	// DD across every repartition, which is what makes them usable under domain decomposition.
	const auto storeLocalAtomSetManager = [this](gmx::LocalAtomSetManager* manager) {
		this->localAtomSetManager_ = manager;
		// Proof that the subscription fired. GROMACS does not deliver this to externally added modules
		// on its own, so if this line never appears the self-subscription in mdrun.cpp has regressed
		// (or GROMACS changed how notifications are dispatched) and LocalAtomSet is unavailable.
		std::cout << "KEnRef: LocalAtomSetManager received (" << (manager ? "non-null" : "NULL")
		          << ") — simulation-setup notification reached the module" << std::endl;
	};
	notifiers->simulationSetupNotifier_.subscribe(storeLocalAtomSetManager);
}

//! Subscribe to pre processing notifications
#if GMX_VERSION >= 20260000
void KEnRefMDModule::subscribeToSimulationRunNotifications(gmx::MDModulesNotifiers* /*notifiers*/) {
    // Nothing to subscribe to: KEnRef reads everything it needs at setup time and per force call.
}
#endif

void KEnRefMDModule::subscribeToPreProcessingNotifications(gmx::MDModulesNotifiers* notifier) {

    std::cerr << "subscribeToPreProcessingNotifications() called +++++++++" << std::endl;
	//TODO Implement. May need Coordinates, PBC and box
//    // Notification of the Coordinates, box and pbc during pre-processing
//    const auto processCoordinatesFunction = [this](const CoordinatesAndBoxPreprocessed& coord) { //needs importing "gromacs/mdrunutility/mdmodulesnotifiers.h"
//        qmmmOptions_.processCoordinates(coord);
//    };
//    notifier->preProcessingNotifier_.subscribe(processCoordinatesFunction);

    // Modification of the topology during pre-processing
    const auto modifyKEnRefTopologyFunction = [/*this*/](gmx_mtop_t* mtop) {
        std::cout << "Topology changed +++++++++" << std::endl;
        std::cout << mtop->name << std::endl <<std::endl;
        // qmmmOptions_.modifyQMMMTopology(mtop);
        // TODO Handle topology
    };
    // Modification of the topology during pre-processing
    notifier->preProcessingNotifier_.subscribe(modifyKEnRefTopologyFunction);
}

void KEnRefMDModule::setSimulationContext(gmx::SimulationContext* simulationContext){
	this->simulationContext_ = simulationContext;
}
