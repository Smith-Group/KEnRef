/*
 * KEnRefSettingsOptions.cpp
 *
 * KEnRef's settings, shared by bin/KEnRef and the GROMACS plugin.
 */
#include "gmxinterface/KEnRefSettingsOptions.h"

#include <iostream>

#include "gmxinterface/gmxwrapper.h"
#include "core/ModelRegistry.h"

namespace kenref {

void addSettingsOptions(CLI::App *app) {
    app->add_flag("--debug", Settings::debug, "enable debugging (holds for debugging)")-> group("");

    app->add_option("--k", Settings::k, "K force constant")->envname("KENREF_K");
    app->add_option("--n", Settings::n, "power scaling")->envname("KENREF_N");
    app->add_option("--max-force", Settings::max_force, "maximum force")->envname("KENREF_MAX_FORCE");

    app->add_option("-g,--guide",  Settings::GUIDE_C_ALPHA, "name of guide group")->envname("KENREF_GUIDE");
    app->add_option("-d,--index", Settings::indexFileName, "Index file name")->envname("KENREF_INDEX");
    app->add_option("-r,--ref", Settings::refFileName, "Reference file")->envname("KENREF_REF");

    app->add_option("-m,--model", Settings::selected_energy_model, "energy model name (e.g. SIGMA, PLATEAUS, RELAX)")
            ->required() ->envname("KENREF_MODEL")
            ->check([](const std::string &name) -> std::string {
                return kenref::ModelRegistry<KEnRef_Real_t>::has(name)
                    ? std::string{} : ("'" + name + "' is not a registered KEnRef energy model");
            });
    app->add_option("-x,--exp-data-folder", Settings::experimentalDataFolder, "experimental data folder for sigma")
        ->envname("KENREF_EXP_DATA_FOLDER"); // ->check(CLI::ExistingDirectory);
    app->add_option("-X,--exp-data-file", Settings::experimentalDataFileName, "experimental data file for plateaus")
        ->envname("KENREF_EXP_DATA_FILE"); //->check(CLI::ExistingFile);

    app->add_option("--atomname-mapping", Settings::atomName_mapping_fileName, "atom name mapping file name")
        ->envname("KENREF_ATOMNAME_MAPPING"); // ->check(CLI::ExistingFile);

    app->add_option("-z,--proton-mhz", Settings::proton_mhz, "spectrometer proton field strength in MHz")->envname("KENREF_PROTON_MHZ");
    // app->add_option("--max-atom-to-read", Settings::max_atom_to_read, "maximum number of atoms to read")->envname("KENREF_MAX_ATOM_TO_READ");
    app->add_option("--alt-out-path", Settings::alt_out_path, "path to redirect output stream")->envname("KENREF_ALT_OUT_PATH");

    // Load from config file
    app->set_config("--params", "params.toml", "Read a TOML config file", false)->envname("KENREF_PARAMS");
}

void useAtomNameMappingAsDefaultReference() {
    if (Settings::refFileName.empty()) {
        std::cout << "Reference file parameter (--ref) not found, using atomname-mapping (" << Settings::atomName_mapping_fileName << ")" << std::endl;
        Settings::refFileName = Settings::atomName_mapping_fileName;
    }else {
        std::cout << "Reference file parameter (--ref) set to , (" << Settings::refFileName << ")" << std::endl;
    }
}

} // namespace kenref
