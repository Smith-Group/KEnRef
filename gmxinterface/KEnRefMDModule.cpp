/*
 * KENRefMDModule.cpp
 *
 *      Author: amr
 */

#include <iostream>
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/topology/topology.h"
#include "../core/IoUtils.h"
#include "../gmxinterface/gmxkenrefinitializer.h"
#include "KEnRefMDModule.h"


KEnRefMDModule::KEnRefMDModule() {
	std::vector<int>const& indices = GmxKEnRefInitializer::loadGmxIndexGroup(KEnRefMDModule::GUIDE_C_ALPHA, KEnRefMDModule::INDEX_FILE_LOCATION);
	this->guideAtoms0Indexed = std::make_shared<std::vector<int> const>(indices);
    IoUtils::printVector(indices);
}

KEnRefMDModule::~KEnRefMDModule() = default;
KEnRefMDModule::KEnRefMDModule(const KEnRefMDModule &other) = default;
KEnRefMDModule::KEnRefMDModule(KEnRefMDModule &&other)  noexcept = default;
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
	forceProviders->addForceProvider(forceProvider_.get());
}

//! Subscribe to simulation setup notifications
void KEnRefMDModule::subscribeToSimulationSetupNotifications(gmx::MDModulesNotifiers* notifiers) {
	//TODO Implement
}

//! Subscribe to pre processing notifications
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
