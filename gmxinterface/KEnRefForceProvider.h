/*
 * KEnRefForceProvider.h
 *
 *      Author: amr
 */

#ifndef KENREFFORCEPROVIDER_H_
#define KENREFFORCEPROVIDER_H_

#include <map>
#include <memory>
#include <Eigen/Dense>
#include <gromacs/mdtypes/iforceprovider.h>
#include <gromacs/mdtypes/forceoutput.h>
#include <gromacs/mdrun/simulationcontext.h>
#include <gromacs/selection/selection.h>



class KEnRefForceProvider: public gmx::IForceProvider {

#define KEnRef_Real float //notice that it is defined in KEnRefForceProvider.h as well. TODO Remove this duplication later
typedef Eigen::Matrix<KEnRef_Real,Eigen::Dynamic, 3, Eigen::RowMajor> CoordsMatrixType;
typedef Eigen::Map<CoordsMatrixType> CoordsMapType;
typedef Eigen::Map<const CoordsMatrixType> CoordsMapTypeConst;   // a read-only map


	gmx::SimulationContext* simulationContext = nullptr;
//	gmx::Selection* selection = nullptr;
//	std::string selectionString = "resid 16 to 20 and pdbname CB"; // "distance from [3., 3., 3.] < 5.0";//"atomnr 6";
	std::shared_ptr<std::vector<int> const> guideAtomIndices;
//	bool *selectionMask = nullptr;
	std::shared_ptr<std::map<std::string, int> const> atomName_to_atomGlobalId_map; //TODO later you may remove this and keep atomName_to_atomSubId_map, or update it and delete atomName_to_atomSubId_map
	std::shared_ptr<std::map<std::string, int>> atomName_to_atomSubId_map;
	std::shared_ptr<std::tuple<std::vector<std::string>, std::vector<std::vector<std::string>>>> simulatedData_table = nullptr; //TODO remove this pointer when it is no longer needed
	std::vector<std::tuple<std::string, std::string>> *atomName_pairs = nullptr;
	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> *g0 = nullptr;
	std::shared_ptr<CoordsMatrixType> subAtomsX;
	std::shared_ptr<std::vector<bool>> globalAtomIdFlags;
	std::shared_ptr<std::vector<int>> subId_to_globalId;
	std::shared_ptr<std::vector<int>> globalId_to_subId;
	std::shared_ptr<CoordsMatrixType> allSimulationsSubAtomsX;

public:
	KEnRefForceProvider();
	virtual ~KEnRefForceProvider();
	KEnRefForceProvider(const KEnRefForceProvider &other);
	KEnRefForceProvider(KEnRefForceProvider &&other);
//	KEnRefForceProvider& operator=(const KEnRefForceProvider &other);
//	KEnRefForceProvider& operator=(KEnRefForceProvider &&other);
    virtual void calculateForces(const gmx::ForceProviderInput& forceProviderInput, gmx::ForceProviderOutput* forceProviderOutput);
    virtual void setSimulationContext(gmx::SimulationContext* simulationContext);
    virtual void setGuideAtomIndices(std::shared_ptr<std::vector<int> const> targetAtomIndices);
	void setAtomNameAtomIdMap(std::map<std::string, int> *atomNameAtomIdMap);
	void setSimulatedDataTable(std::tuple<std::vector<std::string>, std::vector<std::vector<std::string> > > *simulatedDataTable);
#undef KEnRef_Real
};

#endif /* KENREFFORCEPROVIDER_H_ */
