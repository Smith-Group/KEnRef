/*
 * NENRefMDModule.h
 *
 *      Author: amr
 */

#ifndef KENREFMDMODULE_H_
#define KENREFMDMODULE_H_

#include <iostream>
#include <gromacs/utility/keyvaluetreebuilder.h>//TODO you may remove it later
//#include <gromacs/options/ioptionscontainerwithsections.h>
//#include <gromacs/options.h>
#include <gromacs/options/optionsection.h>
#include <gromacs/options/optionfiletype.h>
#include <gromacs/selection/selection.h>
#include <gromacs/mdtypes/imdmodule.h>
#include <gromacs/version.h>
#include <gromacs/mdtypes/imdoutputprovider.h>
#include <gromacs/mdtypes/imdpoptionprovider.h>
#include <gromacs/mdrun/simulationcontext.h>
#include "KEnRefForceProvider.h"
#include "core/IoUtils.h"


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
                    const gmx_output_env_t* /*oenv*/) override{
        std::cout << "====> initOutput() called" << std::endl;
    }
    //! Finalizes output from a simulation run.
    void finishOutput() override {
        std::cout << "====> finishOutput() called" << std::endl;
    }
};

/*
 * \brief Handle the MDP options for KEnRef simulations.
 * empty implementation as KEnRef does not use that yet. TODO It might move to a separate file later.
 */
class KEnRefOptions final : public gmx::IMdpOptionProvider{
    void initMdpTransform(gmx::IKeyValueTreeTransformRules* transform) override {
        std::cout << "====> initMdpTransform() called" << std::endl;
    }
    /*! \brief
     * Initializes options that declare input (mdp) parameters for this
     * module.
     */
    void initMdpOptions(gmx::IOptionsContainerWithSections* options) override {
    	std::cout << "====> initMdpOptions() called" << std::endl;
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
        std::cout << "====> buildMdpOutput() called" << std::endl;
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
    std::shared_ptr<KEnRefForceProvider> forceProvider_;
    gmx::SimulationContext* simulationContext_ = nullptr;
    std::shared_ptr<std::vector<int> const> guideAtoms0Indexed; //ZERO indexed
    std::shared_ptr<const CoordsMatrixType<KEnRef_Real_t>> guideAtomsReferenceCoords_;

    static std::string get_param_value(const std::map<std::string, std::string> &params, const std::string &key) {
        GMX_ASSERT(params.find(key) != params.end(), ("ERROR: Parameter '" + key + "' is not found!").c_str());
        std::string value = params.at(key);
        std::cout << "Parameter '" + key + "' = " << value << std::endl;
        return value;
    }

public:
    KEnRefMDModule();
	~KEnRefMDModule() override;
    [[maybe_unused]] KEnRefMDModule(const KEnRefMDModule &other);
    [[maybe_unused]] KEnRefMDModule(KEnRefMDModule &&other) noexcept;
	KEnRefMDModule& operator=(const KEnRefMDModule &other);
	KEnRefMDModule& operator=(KEnRefMDModule &&other) noexcept;

    //! Returns an interface for handling mdp input (and tpr I/O).
    gmx::IMdpOptionProvider* mdpOptionProvider() override;
    //! Returns an interface for handling output files during simulation.
    gmx::IMDOutputProvider* outputProvider() override;
    //! Initializes force providers from this module.
    void initForceProviders(gmx::ForceProviders* forceProviders) override;
    //! Subscribe to simulation setup notifications
    void subscribeToSimulationSetupNotifications(gmx::MDModulesNotifiers* notifiers) override;
    //! Subscribe to pre processing notifications
    void subscribeToPreProcessingNotifications(gmx::MDModulesNotifiers* notifiers) override;
#if GMX_VERSION >= 20260000
    /*! \brief Subscribe to simulation run notifications.
     *
     * GROMACS 2026 added this as a pure virtual on IMDModule, so it must be present to instantiate
     * the module there; on 2025 the method does not exist and `override` would not compile. KEnRef
     * needs nothing from the run-time notifiers, so the implementation is empty. */
    void subscribeToSimulationRunNotifications(gmx::MDModulesNotifiers* notifiers) override;
#endif
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
    inline static const std::string name_ = "Kinetic-Ensemble-Refinement";
};


#endif /* KENREFMDMODULE_H_ */
