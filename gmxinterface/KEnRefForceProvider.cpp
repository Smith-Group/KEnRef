/*
 * KEnRefForceProvider.cpp
 *
 *      Author: amr
 */

#include <iostream>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/gmxlib/network.h"
#include "KEnRefForceProvider.h"
#include "../core/IoUtils.h"
#include "../core/kabsch.h"
#include "../core/KEnRef.h"
//#include <gromacs/selection.h>
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include <fstream>// needed only during simulated data


#define VERBOSE false
#define KEnRef_Real float //notice that it is defined in KEnRefForceProvider.h as well. TODO Remove this duplication later

KEnRefForceProvider::KEnRefForceProvider() {}

KEnRefForceProvider::~KEnRefForceProvider() {}

KEnRefForceProvider::KEnRefForceProvider(const KEnRefForceProvider &other) {}

KEnRefForceProvider::KEnRefForceProvider(KEnRefForceProvider &&other) {}

//KEnRefForceProvider& KEnRefForceProvider::operator=(const KEnRefForceProvider &other) {}
//KEnRefForceProvider& KEnRefForceProvider::operator=(KEnRefForceProvider &&other) {}

void KEnRefForceProvider::setSimulationContext(gmx::SimulationContext* simulationContext){
	this->simulationContext = simulationContext;
}
void KEnRefForceProvider::setGuideAtomIndices(std::shared_ptr<std::vector<int> const> targetAtomIndices){
	this->guideAtomIndices = targetAtomIndices;
}

void KEnRefForceProvider::calculateForces(const gmx::ForceProviderInput& forceProviderInput, gmx::ForceProviderOutput* forceProviderOutput){
	std::cout << "calculateForces() called" << std::endl;
//	std::cerr << "calculateForces() called" << std::endl;
    const size_t homenr = forceProviderInput.homenr_; // total number of atoms in the system (or domain dec ?)
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
    const auto& force = forceProviderOutput->forceWithVirial_.force_;

	bool isMultiSimulation = this->simulationContext->multiSimulation_ != nullptr;
	int numSimulations = isMultiSimulation ? this->simulationContext->multiSimulation_->numSimulations_ : 1;
	int simulationIndex = isMultiSimulation ? this->simulationContext->multiSimulation_->simulationIndex_: 0;
	MPI_Comm mainRanksComm = isMultiSimulation ? this->simulationContext->multiSimulation_->mainRanksComm_ : MPI_COMM_NULL;
	std::cout 	<< "--> isMultiSimulation: " << std::boolalpha << isMultiSimulation << "\n"
				<< "--> numSimulations " << numSimulations << "\n"
    			<< "--> rankInDefaultCommunicator " << cr.rankInDefaultCommunicator << " " << (isMultiSimulation? simulationIndex : -1) << "\n"
    			<< "--> simulationIndex " << simulationIndex << "\tstep " << step << std::endl;

//    ///////////// selection ////////////////////////////////////
//    gmx::Selection selection = *(this->selection);
////    std::unique_ptr<gmx::SelectionCollection> selectionCollection = std::make_unique<gmx::SelectionCollection>();
////    gmx::SelectionList selectionList = selectionCollection->parseFromString(selectionString);
////    selectionCollection->compile();
////    //TO DO move the above declarations to constructor and use options
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
    //Using gmx::Selection failed. using guideAtoms: a vector<int*>

    if(step == 0){
    	std::cout << "Number of atoms = " << homenr << std::endl;
    	std::cout << "havePPDomainDecomposition(cr): " << havePPDomainDecomposition(&cr) << std::endl;
    	std::cout << "haveDDAtomOrdering(cr): " <<  haveDDAtomOrdering(cr) << std::endl;
    	std::cout << "cr.dd->nnodes: " << cr.dd->nnodes << std::endl;

		if (VERBOSE) {
			std::cout<< "Global to Local Atom number mapping:" << std::endl;
			for(int i = 0; i < homenr; i++){
				const int* aLocal = &i;
				if ((cr.dd == nullptr) || (aLocal = cr.dd->ga2la->findHome(i))){
					std::cout<< i<< " : " << static_cast<size_t>(*aLocal) << std::endl;
				}
			}
		}

		this->atomName_to_atomGlobalId_map = std::make_shared<std::map<std::string, int>>(IoUtils::getAtomNameMappingFromPdb("../6v5d_for_atomname_mapping.pdb"));
		GMX_ASSERT(atomName_to_atomGlobalId_map->size() > 0, "No atom mapping found");
		auto& atomName_to_atomGlobalId_map = *this->atomName_to_atomGlobalId_map;
		if (VERBOSE) {
			for (const auto& [name, globalId] : atomName_to_atomGlobalId_map){
				std::cout << "[" << name << "]\t:" << globalId << std::endl;
			}
//			for(auto& entry: atomName_to_atomGlobalId_map){
//				auto[name, globalId] = entry;
//				std::cout << "[" << name << "]\t:" << globalId << std::endl;
//			}
		}
		this->simulatedData_table = std::make_shared<std::tuple<std::vector<std::string>, std::vector<std::vector<std::string>>>>(IoUtils::readTable("../singleton_data.csv"));
		GMX_ASSERT(simulatedData_table && std::get<1>(*simulatedData_table).size() > 0, "No simulated data found");
		if (VERBOSE) {
			auto [header, data] = *simulatedData_table;
			IoUtils::printVector(header);
			for(auto record: data){
				IoUtils::printVector(record);
			}
		}
		std::vector<std::vector<std::string>> data = std::get<1>(*simulatedData_table);
		this->atomName_pairs = new std::vector<std::tuple<std::string, std::string>>();
		for (auto record : data) {
			std::string atom1 = IoUtils::normalizeName(record[1], true);
			std::string atom2 = IoUtils::normalizeName(record[2], true);
			atomName_pairs->emplace_back(std::make_tuple(atom1, atom2));
		}
		if(VERBOSE){
			for(auto [atom1, atom2]: *atomName_pairs){
				std::cout << "[" << atom1 << "], [" << atom2 << "]" << std::endl;
			}
		}
		this->g0 = new Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>(data.size(), 2);
		for (int i = 0; i < data.size(); ++i) {
			auto record = data[i];
        	std::istringstream temp1(record[3]), temp2(record[4]);
        	temp1 >> (*g0)(i,0);
        	temp2 >> (*g0)(i,1);
		}
		if(VERBOSE)
			std::cout << *g0 << std::endl;

		int maxAtomIdOfInterest = -1;// If you want to use size_t, then you can NOT use -1 as an initial value
		this->globalAtomIdFlags = std::make_shared<std::vector<bool>>(homenr, false);
		auto& globalAtomIdFlags = *this->globalAtomIdFlags;
		//scan the atom pairs to do:
		//1) find the highest globalAtomId number of interest
		//2) fill in the subAtomsFilter
		//3) do a quick sanity scan on the availability of all atomname atomID maping
		int tempI;
		for(auto [a1, a2]: *this->atomName_pairs){
			if(VERBOSE){
				std::cout << "[" << a1 << "]\t" << atomName_to_atomGlobalId_map.at(a1) << "\t";
				std::cout << "[" << a2 << "]\t" << atomName_to_atomGlobalId_map.at(a2) << std::endl;
			}
			//In the next lines I use .at() instead of [] deliberately; to throw an exception if unexpected name found
			if((tempI = atomName_to_atomGlobalId_map.at(a1)) > maxAtomIdOfInterest) maxAtomIdOfInterest = tempI;
			globalAtomIdFlags[tempI] = true;
			if((tempI = atomName_to_atomGlobalId_map.at(a2)) > maxAtomIdOfInterest) maxAtomIdOfInterest = tempI;
			globalAtomIdFlags[tempI] = true;
		}
		if (VERBOSE) {
			IoUtils::printVector(globalAtomIdFlags);
		}
		globalAtomIdFlags.resize(maxAtomIdOfInterest +1);

		this->globalId_to_subId = std::make_shared<std::vector<int>>(globalAtomIdFlags.size(), -1);
		this->subId_to_globalId = std::make_shared<std::vector<int>>(globalAtomIdFlags.size(), -1);
		auto& globalId_to_subId = *this->globalId_to_subId;
		auto& subId_to_globalId = *this->subId_to_globalId;
		int localId = 0;
		for(int i = 0; i < globalAtomIdFlags.size(); i++){
			if(globalAtomIdFlags[i]){
				globalId_to_subId[i] = localId;
				subId_to_globalId[localId] = i;
				localId++;
			}
		}
		subId_to_globalId.resize(localId);

		this->atomName_to_atomSubId_map = std::make_shared<std::map<std::string, int>>();
		auto& atomName_to_atomSubId_map = *this->atomName_to_atomSubId_map;
		for(auto [name, globalId]: atomName_to_atomGlobalId_map){
			atomName_to_atomSubId_map[name] = globalId_to_subId[globalId];
		}
		if (VERBOSE) {
			for(auto entry: atomName_to_atomSubId_map){
				auto[name, localId] = entry;
				std::cout << "[" << name << "]\t:" << localId << std::endl;
			}
		}

		this->subAtomsX = std::make_shared<CoordsMatrixType>(subId_to_globalId.size(), 3);//contains needed atoms only
		if(VERBOSE){auto subAtomsX = *this->subAtomsX; std::cout << "subAtomsX shape is (" << subAtomsX.rows() << ", " << subAtomsX.cols() << ")" << std::endl;}
		this->allSimulationsSubAtomsX = std::make_shared<CoordsMatrixType>(numSimulations * subAtomsX->rows(), 3);
		if(VERBOSE){auto allSimulationsSubAtomsX = *this->allSimulationsSubAtomsX; std::cout << "allSimulationsSubAtomsX shape is (" << allSimulationsSubAtomsX.rows() << ", " << allSimulationsSubAtomsX.cols() << ")" << std::endl;}
    }

    std::vector<int> const& guideAtomIndices = *this->guideAtomIndices;
	auto& atomName_to_atomLocalId_map = *this->atomName_to_atomSubId_map;
//	auto& globalAtomIdFlags = *this->globalAtomIdFlags;
//	auto& atomName_to_atomGlobalId_map = *this->atomName_to_atomGlobalId_map;
//	auto& globalId_to_subId = *this->globalId_to_subId;
	auto& subId_to_globalId = *this->subId_to_globalId;
	auto simulatedData_table = *this->simulatedData_table;
	auto atomName_pairs = *this->atomName_pairs;
	auto g0 = *this->g0;
	CoordsMatrixType& subAtomsX = *this->subAtomsX;
	CoordsMatrixType& allSimulationsSubAtomsX = *this->allSimulationsSubAtomsX;
	//setting the coordinate values to ZERO (or ONE) is dangerous because it causes Invalid floating point operation
	std::vector<std::vector<std::vector<int>>> simulated_grouping_list {{{0}, {1}, {2}}, {{0, 1, 2}}};
	if(!isMultiSimulation){
		simulated_grouping_list = {{{0}}, {{0}}};
	}

	if (VERBOSE) {
		int iZero = 0;
		const int *piZero = cr.dd->ga2la->findHome(iZero);
		std::cout << x[*piZero][0] << " " << x[*piZero][1] << " " << x[*piZero][2] << "\t\t" << std::endl << std::endl;
		//Note that the first atom of guideAtomsX (i.e. guideAtomsX[0]) is not the same subAtomsX[0], and even subAtomsX[0] ((may)) later not be the first atom in the system.
	}

    int guideAtomIndicesSize = guideAtomIndices.size();
    float guideAtomsX_buffer[guideAtomIndicesSize * 3];
    for(auto i =0; i < guideAtomIndicesSize; i++){
    	const int *pi = &guideAtomIndices[i];
    	const int *piLocal = cr.dd->ga2la->findHome(*pi);
    	GMX_ASSERT(piLocal, "ERROR: Can't find local index of atom");
		const gmx::RVec atom_x = x[*piLocal];

		auto rvec = atom_x.as_vec();
		std::copy(rvec, rvec + 3, &guideAtomsX_buffer[i*3]);
//		// for test only
//		for(auto j = 0; j < 3; j++){
//			guideAtomsX_buffer[i * 3 + j] = (isMultiSimulation? 100000 * simulationIndex : 0) + 100 * i + j;
//		}

    }
	CoordsMapType guideAtomsX = CoordsMapType(guideAtomsX_buffer, guideAtomIndicesSize, 3);//TODO CoordsMapType or CoordsMatrixType?
	if(VERBOSE){
		std::cout << "guideAtomsX shape is (" << guideAtomsX.rows() << ", " << guideAtomsX.cols() << ")" << std::endl;
		std::cout << "guideAtomsX" << std::endl << guideAtomsX << std::endl;
	}

    if (haveDDAtomOrdering(cr)){
    	//TODO handle Domain Decomposition
    }

    //TODO restore atom coordinates without the PBC
    // pbc_dx(pbc, *box_const, atom_x, atoms_nopbc[*pii]);

	// Broadcast targetAtomsPositions of model 0
    float guideAtoms_model0_X_buffer[guideAtomIndicesSize * 3]; //TODO Do it the other way around: create the matrix and get its data in buffer
	Eigen::MatrixX3<float> subAtomsXAfterFitting;

	if(simulationIndex == 0){
		std::copy(guideAtomsX_buffer, guideAtomsX_buffer + guideAtomIndicesSize * 3, guideAtoms_model0_X_buffer);
	}else{
//		//Leave unintialized.
		//This initialization is for testing only
//		std::copy(guideAtomsX_buffer, guideAtomsX_buffer + targetAtomIndicesSize * 3, targetAtoms_model0_X_buffer);
	}

//	//For testing only
//	CoordsMapType tempMap = CoordsMapType(guideAtomsX_buffer, 5, 3); //guideAtomIndicesSize, 3);
//	std::cout << "from simulation ((" << simulationIndex << ")) before bCast" << std::endl << tempMap << std::endl;

	if (isMultiSimulation) {
		//Broadcast guideAtomsX from rank 0 to all other ranks
		gmx_bcast(guideAtomIndicesSize * 3 * sizeof(float), guideAtomsX_buffer, mainRanksComm);
		//I don't think this line is important. Only for easy printing
		gmx_barrier(mainRanksComm);
	}

	CoordsMapType model0guideAtomsX = CoordsMapType(guideAtomsX_buffer, guideAtomIndicesSize, 3);

//	//For testing only
//	CoordsMapType tempMatrix1 = CoordsMapType(guideAtomsX_buffer, 5, 3); //guideAtomIndicesSize, 3);
//	std::cout << "from simulation ((" << simulationIndex << ")) after bCast" << std::endl << tempMatrix1 << std::endl;

	Eigen::Transform<KEnRef_Real, 3, Eigen::Affine> affine = (
			simulationIndex == 0 ?
					Eigen::Transform<KEnRef_Real, 3, Eigen::Affine>::Identity() :
					Kabsch<KEnRef_Real>::Find3DAffineTransform(guideAtomsX, model0guideAtomsX));
	if (VERBOSE) {
		std::cout << "Affine Matrix" << std::endl << affine.matrix() << std::endl;
	}

	auto subAtomsX_buffer = subAtomsX.data();
	// Fill needed atoms of subAtomsX with atoms (in the original order). The rest were set to zero earlier
	for(int i = 0; i < subAtomsX.rows(); i++){
		const int *pi = &i;
		const int *piLocal = cr.dd->ga2la->findHome(*pi);
//		GMX_ASSERT(piLocal, "ERROR: Can't find local index of atom");
		const gmx::RVec atom_x = x[*piLocal];
		const auto rvec = atom_x.as_vec();
		std::copy(rvec, rvec + 3, &subAtomsX_buffer[i*3]);
	}

	// Fit every replica (except 0) to model 0
	if(simulationIndex == 0){
		subAtomsXAfterFitting = subAtomsX;
	}else{
		subAtomsXAfterFitting = (subAtomsX.cast<KEnRef_Real>().rowwise().homogeneous() * affine.matrix().transpose()).leftCols(3).cast<float>();
	}

//	// Copy all atomsX into its corresponding section of allSimulationsSubAtomsX after fitting
//	allSimulationsSubAtomsX.setZero();
//	allSimulationsSubAtomsX.block(simulationIndex*subAtomsXAfterFitting.rows(), 0, subAtomsXAfterFitting.rows(), 3);

	// Reduce allSimulationsSubAtomsX to rank 0
	if (isMultiSimulation) {
		MPI_Gather(subAtomsXAfterFitting.data(), subAtomsXAfterFitting.size(), MPI_FLOAT, allSimulationsSubAtomsX.data(), subAtomsXAfterFitting.size(), MPI_FLOAT, 0, mainRanksComm);
		//I don't think this line is important. Only for easy printing
		gmx_barrier(mainRanksComm);
	}else{
		allSimulationsSubAtomsX=subAtomsXAfterFitting;
	}

	if (VERBOSE) {
		if(simulationIndex == 0){
			std::cout << "allSimulationsSubAtomsX shape is (" << allSimulationsSubAtomsX.rows() << ", " << allSimulationsSubAtomsX.cols() << ")" << std::endl;
//			std::cout << "allSimulationsSubAtomsX" << std::endl << allSimulationsSubAtomsX << std::endl;
			for (int r = 0; r < allSimulationsSubAtomsX.size()/10; ++r) {
				for(int i = 0; i< numSimulations; i++){
					for(int j = 0; j < 3; j++){
						std::cout << std::setw(6) << allSimulationsSubAtomsX(subAtomsX.rows() * i + r , j) << " ";
					}
					std::cout << "\t| ";
				}
				std::cout << std::endl;
			}
		}
	}

	float derivatives_buffer[subAtomsX.size()], allDerivatives_buffer[subAtomsX.size()*numSimulations];// TODO try to avoid repeated memory allocation
	CoordsMatrixType derivatives_rectified;
	//in rank 0
	if(simulationIndex == 0){
		//collect all matrices of all replica into a vector of atom coordinates.
		std::vector<Eigen::MatrixX3<float>> allSimulationsSubAtomsX_vector = std::vector<Eigen::MatrixX3<float>>(numSimulations);
		for(int i = 0; i < numSimulations; i++){
			//TODO Review and Simplify this
			CoordsMatrixType temp1Matrix = CoordsMapType(&allSimulationsSubAtomsX.data()[i*subAtomsX.size()], subAtomsX.rows(), 3);
			CoordsMatrixType temp2Matrix = temp1Matrix;
			allSimulationsSubAtomsX_vector.at(i) = temp2Matrix;
		}

//		float energy;
//		std::vector<Eigen::MatrixX3<float>> allDerivatives;
		//do force calculations
		auto [energy, allDerivatives] = KEnRef::coord_array_to_energy(allSimulationsSubAtomsX_vector, atomName_pairs, simulated_grouping_list, g0, 1e-20, atomName_to_atomLocalId_map, true);
		if(VERBOSE){
			std::cout << "energy = " << energy << ", allDerivatives:" << std::endl;
			for(int i = 0; i < allDerivatives.size(); i++){
				std::cout << "model "<< i << " shape (" << allDerivatives[i].rows() << " x " << allDerivatives[i].cols() << ")" << std::endl << allDerivatives[i] << std::endl;
			}
		}
		//I will use the slow method of copying data now, as it is less error prone
		for (int i = 0; i < allDerivatives.size(); ++i) {
			auto matrix = allDerivatives[i];
			std::copy(matrix.data(), matrix.data() + subAtomsX.size(), &allDerivatives_buffer[i*subAtomsX.size()]);
		}

		derivatives_rectified = allDerivatives[0];

	}

	CoordsMapType derivatives_matrix(nullptr, 0, 3);
	std::cout << "Before MPI_Scatter in Thread " << simulationIndex << std::endl;
	if (isMultiSimulation) {
		// Distribute all derivatives
		MPI_Scatter(allDerivatives_buffer, subAtomsX.size(), MPI_FLOAT, derivatives_buffer, subAtomsX.size(), MPI_FLOAT, 0, mainRanksComm);
		std::cout << "After MPI_Scatter in Thread " << simulationIndex << std::endl;
		//once you have the derivatives, retrieve them from the buffer
		new (&derivatives_matrix) CoordsMapType(derivatives_buffer, subAtomsX.rows(), 3);
		std::cout << "After \"CoordsMapType derivatives_matrix\" in Thread [[" << simulationIndex << "]]. derivatives_matrix shape (" << derivatives_matrix.rows() << " x " << derivatives_matrix.cols() << ")"<< std::endl;
	}

	// Transform them back
	if(simulationIndex != 0){
		std::cout << "Before derivatives_rectified in Thread " << simulationIndex << std::endl;
		derivatives_rectified = (derivatives_matrix.cast<KEnRef_Real>().rowwise().homogeneous() * affine.matrix().inverse().transpose()).leftCols(3).cast<float>();
		std::cout << "derivatives_rectified # " << simulationIndex << " shape (" << derivatives_rectified.rows() << " x " << derivatives_rectified.cols() << ")" /*<< std::endl << derivatives_rectified*/ << std::endl;
		std::cout << "After derivatives_rectified in Thread " << simulationIndex << std::endl;
	}

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

    //Finally, add them to corresponding atoms
    std::cout << "Before adding forces in thread " << simulationIndex << std::endl;
    for(int i = 0; i< subAtomsX.rows(); i++){
    	const int *pGlobalId = &subId_to_globalId[i];
//    	//TODO Check whether it use global or local ID?
    	const int *pLocalId = cr.dd->ga2la->findHome(*pGlobalId);
    	force[*pLocalId] += {derivatives_rectified(i, 0), derivatives_rectified(i, 1), derivatives_rectified(i, 2)};// TODO optimize this line/process
    }

    std::cout << "final force values of simulation # " << simulationIndex << std::endl;
    for(int i = 0; i < globalAtomIdFlags->size()/10; i++){
    	std::cout << ((*globalAtomIdFlags).at(i) ? "*" : " ") << "\t" <<force[i][0] << "\t" << force[i][1] << "\t" << force[i][2] << ((*globalAtomIdFlags).at(i) ? "\t*" : "") << std::endl;
    }

	//I don't think this line is important. Only for easy printing
	if(isMultiSimulation) gmx_barrier(mainRanksComm);
    std::cout << "=================" << std::endl;
	gmx_barrier(mainRanksComm);
}

#undef KEnRef_Real
#undef VERBOSE


