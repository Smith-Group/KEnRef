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
#include <chrono>
#include "core/KEnRef.h"


class KEnRefForceProvider : public gmx::IForceProvider {
	//	gmx::Selection* selection = nullptr;
	//	std::string selectionString = "resid 16 to 20 and pdbname CB"; // "distance from [3., 3., 3.] < 5.0";//"atomnr 6";
	//	bool *selectionMask = nullptr;

	gmx::SimulationContext *simulationContext_ = nullptr;
	KEnRef_Real_t maxForceSquared_ = 200.0 * 200.0;
	KEnRef_Real_t k_ = 1.0;
	KEnRef_Real_t n_ = 0.25;
	bool paramsInitialized = false;
	std::shared_ptr<std::vector<int> const> guideAtom0Indices_; //ZERO indexed
	std::shared_ptr<CoordsMatrixType<KEnRef_Real_t> const> guideAtomsReferenceCoords_; //ZERO indexed
	std::shared_ptr<std::map<std::string, int> const> atomName_to_atomGlobalId_map_; //TODO later you may remove this and keep atomName_to_atomSubId_map_, or update it and delete atomName_to_atomSubId_map_
	std::shared_ptr<std::map<std::string, int> > atomName_to_atomSub0Id_map_; //atomName is normalized string. SubId is a small subset and is ZERO based
	std::shared_ptr<std::tuple<std::vector<std::string>, std::vector<std::vector<std::string> > > > experimentalData_table_ = nullptr; //TODO remove this pointer when it is no longer needed
	std::vector<std::tuple<std::string, std::string> > *atomName_pairs_ = nullptr;
	Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, Eigen::Dynamic> *g0_ = nullptr;
	std::shared_ptr<std::vector<std::tuple<int, int> > > atomId_pairs_;
	std::shared_ptr<CoordsMatrixType<KEnRef_Real_t> > subAtomsX_;
	std::shared_ptr<std::vector<bool> > globalAtomIdFlags_; //ONE based
	std::shared_ptr<std::vector<int> > sub0Id_to_global1Id_; //Global ID is ONE based, subId is a small subset and is ZERO based
	std::shared_ptr<std::vector<int> > global1Id_to_sub0Id_; //Global ID is ONE based, subId is a small subset and is ZERO based
	std::shared_ptr<CoordsMatrixType<KEnRef_Real_t> > allSimulationsSubAtomsX_; //allocated once, used every step
	std::shared_ptr<KEnRef_Real_t[]> allDerivatives_buffer_;
	std::shared_ptr<KEnRef_Real_t[]> derivatives_buffer_;
	long long calculateForces_time = 0;

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	KEnRefForceProvider();

	virtual ~KEnRefForceProvider();

	KEnRefForceProvider(const KEnRefForceProvider &other);

	KEnRefForceProvider(KEnRefForceProvider &&other) noexcept;

	KEnRefForceProvider &operator=(const KEnRefForceProvider &other);

	KEnRefForceProvider &operator=(KEnRefForceProvider &&other) noexcept;

	void calculateForces(const gmx::ForceProviderInput &forceProviderInput,
	                     gmx::ForceProviderOutput *forceProviderOutput) override;

	virtual void setSimulationContext(gmx::SimulationContext *simulationContext);

	virtual void setGuideAtom0Indices(std::shared_ptr<std::vector<int> const> targetAtoms0Indices);

	virtual void setGuideAtomsReferenceCoords(
		std::shared_ptr<const CoordsMatrixType<KEnRef_Real_t>> &guideAtomsReferenceCoords);

	void fillParamsStep0(size_t homenr, int numSimulations);

	static CoordsMatrixType<KEnRef_Real_t>
	getGuideAtomsX(const gmx::ArrayRef<const gmx::RVec> &x, const t_commrec &cr,
	               const std::vector<int> &guideAtom0Indices);
};

#endif /* KENREFFORCEPROVIDER_H_ */
