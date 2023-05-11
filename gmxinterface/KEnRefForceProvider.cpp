/*
 * KEnRefForceProvider.cpp
 *
 *      Author: amr
 */

#include <iostream>
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/gmxlib/network.h"
#include "KEnRefForceProvider.h"
#include "../core/IoUtils.h"
//#include <gromacs/selection.h>
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include <Eigen/Dense>
#include <Eigen/Core>


KEnRefForceProvider::KEnRefForceProvider() {
	std::cout << "KEnRefForceProvider constructed [cout]" << std::endl;
	std::cerr << "KEnRefForceProvider constructed [cerr]" << std::endl;
}

KEnRefForceProvider::~KEnRefForceProvider() {
	// TODO Auto-generated destructor stub
}

KEnRefForceProvider::KEnRefForceProvider(const KEnRefForceProvider &other) {
	// TODO Auto-generated constructor stub

}

KEnRefForceProvider::KEnRefForceProvider(KEnRefForceProvider &&other) {
	// TODO Auto-generated constructor stub

}

void KEnRefForceProvider::setSimulationContext(gmx::SimulationContext* simulationContext){
	this->simulationContext = simulationContext;
}
void KEnRefForceProvider::setTargetAtomIndices(std::shared_ptr<std::vector<int> const> targetAtomIndices){
	this->targetAtomIndices = targetAtomIndices;
}

void KEnRefForceProvider::calculateForces(const gmx::ForceProviderInput& forceProviderInput, gmx::ForceProviderOutput* forceProviderOutput){
	std::cout << "calculateForces() called" << std::endl;
    const int homenr = forceProviderInput.homenr_;
    GMX_ASSERT(homenr >= 0, "number of home atoms must be non-negative.");

//    const auto& box = forceProviderInput.box_;
    matrix box = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    copy_mat(forceProviderInput.box_, box);

    GMX_ASSERT(check_box(PbcType::Unset, box) == nullptr, "Invalid box.");
    t_pbc temp{};
    t_pbc *pbc = &temp;
    set_pbc(pbc, PbcType::Unset, box);

    const auto& x  = forceProviderInput.x_;
    const auto& cr = forceProviderInput.cr_;
    const auto& step = forceProviderInput.step_;
//    const auto& t  = forceProviderInput.t_; // Not needed (at least yet)

    std::cout << "--> rankInDefaultCommunicator " << cr.rankInDefaultCommunicator;
    if(this->simulationContext->multiSimulation_){
    	std::cout<< " " << this->simulationContext->multiSimulation_->simulationIndex_;
    }
    std::cout << std::endl;
//	std::cerr << "calculateForces() called" << std::endl;
    const auto& force = forceProviderOutput->forceWithVirial_.force_;

    if(this->simulationContext->multiSimulation_){
    	std::cout<< "--> " << this->simulationContext->multiSimulation_->simulationIndex_ << "\t";
    }

//    ///////////// selection ////////////////////////////////////
//    gmx::Selection selection = *(this->selection);
////    std::unique_ptr<gmx::SelectionCollection> selectionCollection = std::make_unique<gmx::SelectionCollection>();
////    gmx::SelectionList selectionList = selectionCollection->parseFromString(selectionString);
////    selectionCollection->compile();
////    //TODO move the above declarations to constructor and use options
////    for (gmx::Selection selection: selectionList) {
//        std::cout << selection.name() << ":" << selection.posCount() << std::endl;
//        if(selection.atomIndices().empty())
//        	std::cout << "selection empty" << std::endl;
//        else{
//			for (int atomIdx: selection.atomIndices()){
//				std::cout << atomIdx << " ";
//			}
//        }
//		std::cout << std::endl;
////	}
//    ////////////// End Selection ////////////////////////////////////////////////
    //Using gmx::Selection failed. using targetAtoms: a vector<int*>
    std::vector<int> const& targetAtomIndices = *this->targetAtomIndices;
    if(step == 0){
    	std::cout << "Number of atoms = " << homenr << std::endl;
    	this->selectionMask = new bool[homenr];
    	std::fill(selectionMask, selectionMask+homenr, false);
    	IoUtils::printVector(targetAtomIndices);
    	for(int idx: targetAtomIndices){
    		selectionMask[idx] = true;
    	}
    	std::cout << "havePPDomainDecomposition(cr): " << havePPDomainDecomposition(&cr) << std::endl;
    	std::cout << "haveDDAtomOrdering(cr): " <<  haveDDAtomOrdering(cr) << std::endl;
    	std::cout << "cr.dd->nnodes: " << cr.dd->nnodes << std::endl;

		if (false) {
			std::cout<< "Atom number mapping:" << std::endl;
			for(int i = 0; i < homenr; i++){
				const int* aLocal = &i;
				if ((cr.dd == nullptr) || (aLocal = cr.dd->ga2la->findHome(i))){
					std::cout<< i<< " : " << static_cast<size_t>(*aLocal) << std::endl;
				}
			}
		}
    }


	if (false) {
		int iZero = 0;
		const int *piZero = cr.dd->ga2la->findHome(iZero);
		std::cout << x[*piZero][0] << " " << x[*piZero][1] << " " << x[*piZero][2] << "\t\t" << std::endl << std::endl;
	}


    auto && selectionMask = &this->selectionMask;
	//TODO extract targetAtomsPositions (maybe using Eigen selction expression)
//    for(size_t i = 0; i < homenr; i++){
//    	bool b = selectionMask[i];
//    	std::cout << b << " ";
//    }
//    std::cout << std::endl;

    size_t targetAtomIndicesSize = this->targetAtomIndices.get()->size();
    float targetAtomsX_buffer[targetAtomIndicesSize * 3];
    for(auto i =0; i < this->targetAtomIndices.get()->size(); i++){
    	const int *pi = &i;
    	const int *pilocal = cr.dd->ga2la->findHome(*pi);
    	GMX_ASSERT(pilocal, "ERROR: Can't find local index of atom");
		const gmx::RVec atom_x = x[*pilocal];
//		const gmx::BasicVector<float>::RawArray& xi_asVec = x[*pilocal].toRVec();
//		const gmx::BasicVector<float>::RawArray& xi_asVec = x[*pilocal].as_vec();
//    	float* xi_asVec = x[*pilocal].as_vec();

		//TODO Is there a faster way than this for direct mapping of coordinates?
		for(auto j = 0; j < 3; j++){
			targetAtomsX_buffer[i * 3 + j] = atom_x[j];
			//FIXME test only
			targetAtomsX_buffer[i * 3 + j] = (this->simulationContext->multiSimulation_? 100000 * this->simulationContext->multiSimulation_->simulationIndex_ : 0) + 100 * i + j;
		}
    }
	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> targetAtomsX = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>>(targetAtomsX_buffer, targetAtomIndicesSize, 3); //(targetAtomIndicesSize, 3);
    std::cout << targetAtomsX << std::endl;

    if (haveDDAtomOrdering(cr)){
    	//TODO handle Domain Decomposition
    }

    //TODO restore atom coordinates without the PBC
    // pbc_dx(pbc, *box_const, atom_x, atoms_nopbc[*pii]);

	// Broadcast targetAtomsPositions of model 0
    float targetAtoms_model0_X_buffer[targetAtomIndicesSize * 3];
    if(this->simulationContext->multiSimulation_ && this->simulationContext->multiSimulation_->simulationIndex_ == 0){
    	std::copy(targetAtomsX_buffer, targetAtomsX_buffer + targetAtomIndicesSize * 3, targetAtoms_model0_X_buffer);
    }else{
    	//FIXME leave unintialized. this initialization is for testing only
    	std::copy(targetAtomsX_buffer, targetAtomsX_buffer + targetAtomIndicesSize * 3, targetAtoms_model0_X_buffer);
    }
	if (this->simulationContext->multiSimulation_ ) {
//		//FIXME for testing only
//		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> tempMatrix = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>>(targetAtomsX_buffer, 5, 3); //(targetAtomIndicesSize, 3);
//		std::cout << "from simulation ((" << this->simulationContext->multiSimulation_->simulationIndex_ << ")) before bCast" << std::endl << tempMatrix << std::endl;

		gmx_bcast(targetAtomIndicesSize * sizeof(float), targetAtomsX_buffer, this->simulationContext->multiSimulation_->mainRanksComm_);
		//I don't think this line is important. only for easy printing
		gmx_barrier(this->simulationContext->multiSimulation_->mainRanksComm_);
		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> model0TargetAtomsX = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>>(targetAtomsX_buffer, targetAtomIndicesSize, 3);

//		//FIXME for testing only
//		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> tempMatrix1 = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>>(targetAtomsX_buffer, 5, 3); //(targetAtomIndicesSize, 3);
//		std::cout << "from simulation ((" << this->simulationContext->multiSimulation_->simulationIndex_ << ")) after bCast" << std::endl << tempMatrix1 << std::endl;
	}

	//TODO fit every replica (except 0) to model 0
	//TODO calculate position distances locally

	//Either this block
	//TODO scatter and gather the calculated distances as(into) average(s) at ranks of every model.

	//or this block
	//TODO Reduce the fitted positions as(into) average(s) at rank of model 0
	//TODO calculate average position distances at rank of model 0

	//in both cases, you can alternatively scatter and gather the passed-in data in case you want to have the average on ranks of every model.

	//TODO do force calculations (call KEnRef::calc_derivative() method).


	//////////////////////////////////////////////////////////////////////////////////////////////////
    //	std::cout << forceProviderInput.x_ << std::endl << std::endl;
    std::cout << forceProviderOutput << std::endl;
    //	std::cout << forceProviderOutput->forceWithVirial_ << std::endl;
    std::cout << forceProviderOutput->forceWithVirial_.getVirial() << std::endl;
    std::cout << "computeVirial_ = " << std::boolalpha  << forceProviderOutput->forceWithVirial_.computeVirial_ << std::endl;
    ////	const gmx::ArrayRef<gmx::BasicVector<float> > 	force = forceProviderOutput->forceWithVirial_.force_;
    //	for (int i = 0; i < 5 /*force.size()/100*/; ++i) {
    //		std::cout  << "Sum forces on Atom # " << i << ":" << std::endl;
    ////		auto forceitem = force[i];
    //		for(int j = 0; j < 3; j++){
    //			std::cout << forceProviderInput.x_[i][j] << " ";
    //			//			std::cout << forceitem[j] << " ";
    //		}
    //		std::cout << std::endl;
    //	}
    //	std::cout << std::endl << std::endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////




    if (haveDDAtomOrdering(cr)){
        // Note: this assumes that all ranks are hitting this line, which is not generally true.
        // I need to find the right subcommunicator. What I really want is a _scoped_ communicator...
        gmx_barrier(cr.mpi_comm_mygroup);
    }

//    // From /Gromacs-Source@gromacs/src/gromacs/restraint/restraintmdmodule.cpp
//    const int  site1  = static_cast<int>(sites_.front().index());
//    const int* aLocal = &site1;
//    // Set forces using index `site1` if no domain decomposition, otherwise set with local index if available.
//    const auto& force = forceProviderOutput->forceWithVirial_.force_;
//    if ((cr.dd == nullptr) || (aLocal = cr.dd->ga2la->findHome(site1))){
//        force[static_cast<size_t>(*aLocal)] += result.force;
//    }


//    for(int i = 0; i< forceProviderInput.x_.size(); i++)
//    	force[0] += {0.0, 10.0, 0.0}; //THIS is the prototype line


}


//KEnRefForceProvider& KEnRefForceProvider::operator=(
//		const KEnRefForceProvider &other) {
//	// TODO Auto-generated method stub
//}
//KEnRefForceProvider& KEnRefForceProvider::operator=(
//		KEnRefForceProvider &&other) {
//	// TODO Auto-generated method stub
//}

