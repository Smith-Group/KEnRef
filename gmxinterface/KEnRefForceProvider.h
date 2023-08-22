/*
 * KEnRefForceProvider.h
 *
 *      Author: amr
 */

#ifndef KENREFFORCEPROVIDER_H_
#define KENREFFORCEPROVIDER_H_

#include <map>
#include <memory>
#include <gromacs/mdtypes/iforceprovider.h>
#include <gromacs/mdtypes/forceoutput.h>
#include <gromacs/mdrun/simulationcontext.h>
#include <gromacs/selection/selection.h>
#include "../core/KEnRef.h"



class KEnRefForceProvider: public gmx::IForceProvider {

	gmx::SimulationContext* simulationContext_ = nullptr;
    KEnRef_Real maxForce_ = 100.0;
//	gmx::Selection* selection = nullptr;
//	std::string selectionString = "resid 16 to 20 and pdbname CB"; // "distance from [3., 3., 3.] < 5.0";//"atomnr 6";
	std::shared_ptr<std::vector<int> const> guideAtomIndices_;
//	bool *selectionMask = nullptr;
	std::shared_ptr<std::map<std::string, int> const> atomName_to_atomGlobalId_map_; //TODO later you may remove this and keep atomName_to_atomSubId_map_, or update it and delete atomName_to_atomSubId_map_
	std::shared_ptr<std::map<std::string, int>> atomName_to_atomSubId_map_;
	std::shared_ptr<std::tuple<std::vector<std::string>, std::vector<std::vector<std::string>>>> simulatedData_table_ = nullptr; //TODO remove this pointer when it is no longer needed
	std::vector<std::tuple<std::string, std::string>> *atomName_pairs_ = nullptr;
	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> *g0_ = nullptr;
	std::shared_ptr<CoordsMatrixType> subAtomsX_{};
	std::shared_ptr<std::vector<bool>> globalAtomIdFlags_;
	std::shared_ptr<std::vector<int>> subId_to_globalId_;
	std::shared_ptr<std::vector<int>> globalId_to_subId_;
	std::shared_ptr<CoordsMatrixType> allSimulationsSubAtomsX_{};

public:
	KEnRefForceProvider();
	virtual ~KEnRefForceProvider();
//	KEnRefForceProvider(const KEnRefForceProvider &other);
//	KEnRefForceProvider(KEnRefForceProvider &&other) noexcept ;
//	KEnRefForceProvider& operator=(const KEnRefForceProvider &other);
//	KEnRefForceProvider& operator=(KEnRefForceProvider &&other);
    void calculateForces(const gmx::ForceProviderInput& forceProviderInput, gmx::ForceProviderOutput* forceProviderOutput) override;
    virtual void setSimulationContext(gmx::SimulationContext* simulationContext);
    virtual void setGuideAtomIndices(std::shared_ptr<std::vector<int> const> targetAtomIndices);
    void fillParamsStep0(size_t homenr, int numSimulations);
};

#endif /* KENREFFORCEPROVIDER_H_ */
