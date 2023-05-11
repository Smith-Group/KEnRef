/*
 * NENRefMDModule.cpp
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
	const char* H_HA_HA1_HA2 = "H_HA_HA1_HA2";
	const char* indexFileLocation = "../../res/josh/KEnRefAtomIndex.ndx";
//	auto indexGroups = GmxKEnRefInitializer::loadGmxIndexFile(indexFileLocation);
//	for(auto group: indexGroups){
//		std::cout << "'" << group.first << "':" << std::endl;
//		IoUtils::printVector(group.second);
//	}
//	const std::vector<int>& indices = indexGroups[H_HA_HA1_HA2];
//	std::cout << "[" << H_HA_HA1_HA2 << "]:" << std::endl;
//	IoUtils::printVector(indices);
	std::vector<int>const& indices2 = GmxKEnRefInitializer::loadGmxIndexGroup(H_HA_HA1_HA2, indexFileLocation);
	IoUtils::printVector(indices2);
//	const auto& indices3 = IoUtils::getGmxNdxGroup(indexFileLocation, H_HA_HA1_HA2);
//	IoUtils::printVector(indices3);

	this->targetAtoms = std::make_shared<std::vector<int> const>(indices2.begin(), indices2.end());
//	for(auto v : *vv1){std::cout << v << " ";}std::cout << std::endl;
//	std::cout << this->targetAtoms << std::endl;

}

KEnRefMDModule::~KEnRefMDModule() {
	// TODO Auto-generated destructor stub
}

KEnRefMDModule::KEnRefMDModule(const KEnRefMDModule &other) {
	// TODO Auto-generated constructor stub
}

KEnRefMDModule::KEnRefMDModule(KEnRefMDModule &&other) {
	// TODO Auto-generated constructor stub
}

//KEnRefMDModule& NENRefMDModule::operator=(const NENRefMDModule &other) {
//	// TODO Auto-generated method stub
//}
//KEnRefMDModule& NENRefMDModule::operator=(NENRefMDModule &&other) {
//	// TODO Auto-generated method stub
//}

//! Returns an interface for handling mdp input (and tpr I/O).
gmx::IMdpOptionProvider* KEnRefMDModule::mdpOptionProvider() {
	return &kEnRefOptions;
}

//! Returns an interface for handling output files during simulation.
gmx::IMDOutputProvider* KEnRefMDModule::outputProvider() {
	return &kEnRefOutputProvider;
}

//! Initializes force provider(s) from this module (only one now).
void KEnRefMDModule::initForceProviders(gmx::ForceProviders* forceProviders) {
	std::cout << "KEnRefMDModule::initForceProviders()" << std::endl;
	forceProvider_ = std::make_unique<KEnRefForceProvider>();
	forceProvider_->setSimulationContext(simulationContext);
	forceProvider_->setTargetAtomIndices(this->targetAtoms);
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
	this->simulationContext = simulationContext;
}

const std::string KEnRefModuleInfo::name_ = "Kinetic-Ensemble-Refinement";

