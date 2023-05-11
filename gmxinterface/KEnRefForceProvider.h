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

class KEnRefForceProvider: public gmx::IForceProvider {
	gmx::SimulationContext* simulationContext = nullptr;
//	gmx::Selection* selection = nullptr;
//	std::string selectionString = "resid 16 to 20 and pdbname CB"; // "distance from [3., 3., 3.] < 5.0";//"atomnr 6";
	std::shared_ptr<std::vector<int> const> targetAtomIndices;
	bool *selectionMask = nullptr;
public:
	KEnRefForceProvider();
	virtual ~KEnRefForceProvider();
	KEnRefForceProvider(const KEnRefForceProvider &other);
	KEnRefForceProvider(KEnRefForceProvider &&other);
//	KEnRefForceProvider& operator=(const KEnRefForceProvider &other);
//	KEnRefForceProvider& operator=(KEnRefForceProvider &&other);
    virtual void calculateForces(const gmx::ForceProviderInput& forceProviderInput, gmx::ForceProviderOutput* forceProviderOutput);
    virtual void setSimulationContext(gmx::SimulationContext* simulationContext);
    virtual void setTargetAtomIndices(std::shared_ptr<std::vector<int> const> targetAtomIndices);
};

#endif /* KENREFFORCEPROVIDER_H_ */
