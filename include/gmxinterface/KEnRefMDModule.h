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

    static void readParams(const char *kenref_params) {
        const std::map<std::string, std::string> &params = IoUtils::readParams(kenref_params);
        KEnRefMDModule::GUIDE_C_ALPHA = params.at("GUIDE_C_ALPHA");
        KEnRefMDModule::INDEX_FILE_LOCATION = params.at("INDEX_FILE_LOCATION");
        KEnRefMDModule::ATOMNAME_MAPPING_FILENAME = params.at("ATOMNAME_MAPPING_FILENAME");
        KEnRefMDModule::EXPERIMENTAL_DATA_FILENAME = params.at("EXPERIMENTAL_DATA_FILENAME");
        auto refFileLocation = params.find("REFERENCE_FILENAME");
        if (refFileLocation == params.end()){ //NOT Found
            KEnRefMDModule::REFERENCE_FILENAME = ATOMNAME_MAPPING_FILENAME;
            std::cout << "REFERENCE_FILENAME not found, using ATOMNAME_MAPPING_FILENAME (" << ATOMNAME_MAPPING_FILENAME << ")" << std::endl;
        }else{ //Found
            KEnRefMDModule::REFERENCE_FILENAME = refFileLocation->second;
            std::cout << "REFERENCE_FILENAME found, (" << REFERENCE_FILENAME << ")" << std::endl;
        }
    }
public:
    inline static std::string GUIDE_C_ALPHA; // = "guideC-alpha";
    inline static std::string INDEX_FILE_LOCATION; // = "../../res/cleanstart/KEnRefAtomIndex.ndx";
    inline static std::string REFERENCE_FILENAME; // = "../../res/cleanstart/6v5d_step0_for_atomname_mapping.pdb";
    inline static std::string ATOMNAME_MAPPING_FILENAME; // = "../../res/cleanstart/6v5d_step0_for_atomname_mapping.pdb";
    inline static std::string EXPERIMENTAL_DATA_FILENAME; // = "../../res/cleanstart/singleton_data_step0_model01.csv";

    KEnRefMDModule();
	~KEnRefMDModule() override;
	KEnRefMDModule(const KEnRefMDModule &other);
	KEnRefMDModule(KEnRefMDModule &&other) noexcept;
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
