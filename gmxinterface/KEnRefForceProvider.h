/*
 * KEnRefForceProvider.h
 *
 *      Author: amr
 */

#ifndef KENREFFORCEPROVIDER_H_
#define KENREFFORCEPROVIDER_H_

#include <gromacs/mdtypes/iforceprovider.h>
#include <gromacs/mdtypes/forceoutput.h>
#include <gromacs/mdrun/simulationcontext.h>
#include <gromacs/selection.h>
#include <map>
#include <Eigen/Dense>

class KEnRefForceProvider: public gmx::IForceProvider {
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
	std::shared_ptr<Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>> subAtomsX;
	std::shared_ptr<std::vector<bool>> globalAtomIdFlags;
	std::shared_ptr<std::vector<int>> subId_to_globalId;
	std::shared_ptr<std::vector<int>> globalId_to_subId;
	std::shared_ptr<Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>> allSimulationsSubAtomsX;

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
};

#endif /* KENREFFORCEPROVIDER_H_ */
