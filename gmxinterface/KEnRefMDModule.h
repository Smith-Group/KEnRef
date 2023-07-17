/*
 * NENRefMDModule.h
 *
 *      Author: amr
 */

#ifndef KENREFMDMODULE_H_
#define KENREFMDMODULE_H_

#include <iostream>
#include <gromacs/utility/keyvaluetreebuilder.h>//TODO you may remove it later
//#include <gromacs/options.h>
#include <gromacs/options/optionsection.h>
#include <gromacs/options/optionfiletype.h>
#include <gromacs/selection/selection.h>
//#include <gromacs/options/ioptionscontainerwithsections.h>
#include <gromacs/mdtypes/imdmodule.h>
#include <gromacs/mdtypes/imdoutputprovider.h>
#include <gromacs/mdtypes/imdpoptionprovider.h>
#include <gromacs/mdrun/simulationcontext.h>
#include "KEnRefForceProvider.h"


/*
 * \brief Handle file output for KEnRef simulations.
 * empty implementation as KEnRef does not use that yet. TODO It might move to a separate file later.
 */
class KEnRefOutputProvider final : public gmx::IMDOutputProvider{
public:
    //! Initialize output
    void initOutput(FILE* /*fplog*/,
                    int /*nfile*/,
                    const t_filenm /*fnm*/[],
                    bool /*bAppendFiles*/,
                    const gmx_output_env_t* /*oenv*/) override{}
    //! Finalizes output from a simulation run.
    void finishOutput() override {}
};

/*
 * \brief Handle the MDP options for KEnRef simulations.
 * empty implementation as KEnRef does not use that yet. TODO It might move to a separate file later.
 */
class KEnRefOptions final : public gmx::IMdpOptionProvider{
    virtual void initMdpTransform(gmx::IKeyValueTreeTransformRules* transform) override {}
    /*! \brief
     * Initializes options that declare input (mdp) parameters for this
     * module.
     */
    void initMdpOptions(gmx::IOptionsContainerWithSections* options) override {
//    	std::cout << "====> initMdpOptions() called" << std::endl;
//        auto section = options->addSection(gmx::OptionSection("KEnRef Options"));
//         section.addOption(gmx::FileNameOption("top")
//                               .filetype(gmx::OptionFileType::Topology)
//                               .inputFile()
//                               .store(&topfile_)
//                               .defaultBasename("topol")
//                               .description("Topology file"));
//         section.addOption(gmx::SelectionOption("select")
////                               .storeVector(&selections_).multiValue()
//                               .store(&selection)
//							   .defaultSelectionText("resid 16 to 20 and pdbname CB")
//                               .description("Selection string for atoms used for calculations"));
    }
    //! Prepares to write a flat key-value tree like an mdp file.
    void buildMdpOutput(gmx::KeyValueTreeObjectBuilder* builder) const override {
        builder->addValue("top", "topol.tpr");
        builder->addValue("select", selection.name()/*"resnum 5"*/);
    }

public:
	gmx::Selection& getSelection() /*const*/ {
		return selection;
	}
private:
    std::string topfile_;
//    gmx::SelectionList selections_;
    gmx::Selection selection;
};

class KEnRefMDModule final: public gmx::IMDModule {

	KEnRefOptions kEnRefOptions;
	KEnRefOutputProvider kEnRefOutputProvider;
    //! Object that evaluates the forces
    std::unique_ptr<KEnRefForceProvider> forceProvider_;
    gmx::SimulationContext* simulationContext = nullptr;
    std::shared_ptr<std::vector<int> const> guideAtoms ; //= nullptr; // std::make_shared<std::vector<int>>();

public:
	KEnRefMDModule();
	virtual ~KEnRefMDModule();
	KEnRefMDModule(const KEnRefMDModule &other);
	KEnRefMDModule(KEnRefMDModule &&other) noexcept;
//	KEnRefMDModule& operator=(const NENRefMDModule &other);
//	KEnRefMDModule& operator=(NENRefMDModule &&other);

    //! Returns an interface for handling mdp input (and tpr I/O).
    virtual gmx::IMdpOptionProvider* mdpOptionProvider() override;
    //! Returns an interface for handling output files during simulation.
    virtual gmx::IMDOutputProvider* outputProvider() override;
    //! Initializes force providers from this module.
    virtual void initForceProviders(gmx::ForceProviders* forceProviders) override;
    //! Subscribe to simulation setup notifications
    virtual void subscribeToSimulationSetupNotifications(gmx::MDModulesNotifiers* notifiers) override;
    //! Subscribe to pre processing notifications
    virtual void subscribeToPreProcessingNotifications(gmx::MDModulesNotifiers* notifiers) override;
    virtual void setSimulationContext(gmx::SimulationContext* simulationContext);

};

/*! \libinternal \brief Information about the Kinetic Ensemble Refinement module.
 *
 * Provides name and method to create a Kinetic Ensemble Refinement module.
 */
struct KEnRefModuleInfo
{
    /*! \brief
     * Creates a module for applying forces to refine NMR structures.
     */
    static std::unique_ptr<gmx::IMDModule> create(){
        return std::make_unique<KEnRefMDModule>();
    }
    //! The name of the module
    static const std::string name_;
};


#endif /* KENREFMDMODULE_H_ */
