/*
 * gmxwrapper.cpp
 *
 *  Created on: Nov 29, 2022
 *      Author: amr
 */
#include <string>
#include <map>
#include <unistd.h>
#include <CLI11/CLI11.hpp>
#include "gmxpre.h"
#include "mdrun/mdrun_main.h"
#include "gromacs/version.h"
#include "gromacs/commandline/cmdlinemodule.h"
#if GMX_VERSION >= 20260000
// GROMACS 2026 split CommandLineModuleSettings out of cmdlinemodule.h, which now only
// forward-declares it; the definition is needed for setDefaultNiceLevel() below.
#    include "gromacs/commandline/cmdlinemodulesettings.h"
#endif
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gmxinterface/gmxwrapper.h"
#include "core/KEnRef.h"
// Remember to put all local imports (especially those which reference Eigen) after
// all Gromacs imports; to avoid compiler error "reference to 'real' is ambiguous"

GMXWrapper::GMXWrapper() = default;
GMXWrapper::~GMXWrapper() = default;
GMXWrapper::GMXWrapper(const GMXWrapper &other) = default;
GMXWrapper::GMXWrapper(GMXWrapper &&other) noexcept = default;
GMXWrapper &GMXWrapper::operator=(const GMXWrapper &other) = default;
GMXWrapper &GMXWrapper::operator=(GMXWrapper &&other) noexcept = default;

//! Initializer for a module that defaults to nice level zero. //copied from mdrun_main.cpp
void initSettingsNoNice(gmx::CommandLineModuleSettings *settings) {
    settings->setDefaultNiceLevel(0);
}

struct ArgvDeleter {
    void operator()(char **argv) const {
        if (argv) {
            for (int i = 0; argv[i] != nullptr; ++i)
                free(argv[i]);
            delete[] argv;
        }
    }
};

using ArgvPtr = std::unique_ptr<char *[], ArgvDeleter>;

ArgvPtr create_argv(int argc, char *original_argv[], const std::vector<std::string> &remaining) {
    // Allocate array of char*
    // ReSharper disable once CppDFAMemoryLeak
    char **new_argv = new char *[remaining.size() + 2](); // +2 for program name and nullptr
    // Duplicate program name
    new_argv[0] = strdup(original_argv[0]);
    // Duplicate remaining arguments
    for (size_t i = 0; i < remaining.size(); ++i)
        new_argv[i + 1] = strdup(remaining[i].c_str());
    // Last element is nullptr (conventional)
    new_argv[remaining.size() + 1] = nullptr;
    return ArgvPtr(new_argv);
}


int main(int argc, char *argv[]) {
    // This should be the main function that loads GMXWrapper (which should replace GROMACS main class)
    CLI::App app{"KEnRef"};

    app.add_flag("--debug", Settings::debug, "enable debugging (holds for debugging)")-> group("");

    app.add_option("--k", Settings::k, "K force constant")->envname("KENREF_K");
    app.add_option("--n", Settings::n, "power scaling")->envname("KENREF_N");
    app.add_option("--max-force", Settings::max_force, "maximum force")->envname("KENREF_MAX_FORCE");

    app.add_option("-g,--guide",  Settings::GUIDE_C_ALPHA, "name of guide group")->envname("KENREF_GUIDE");
    app.add_option("-d,--index", Settings::indexFileName, "Index file name")->envname("KENREF_INDEX");
    app.add_option("-r,--ref", Settings::refFileName, "Reference file")->envname("KENREF_REF");

    app.add_option("-m,--model", Settings::selected_energy_model, "energy model")
            ->required() ->envname("KENREF_MODEL")
            ->transform(CLI::CheckedTransformer(KEnRef<KEnRef_Real_t>::energyModelMap, CLI::ignore_case));
    app.add_option("-x,--exp-data-folder", Settings::experimentalDataFolder, "experimental data folder for sigma")
        ->envname("KENREF_EXP_DATA_FOLDER"); // ->check(CLI::ExistingDirectory);
    app.add_option("-X,--exp-data-file", Settings::experimentalDataFileName, "experimental data file for plateaus")
        ->envname("KENREF_EXP_DATA_FILE"); //->check(CLI::ExistingFile);

    app.add_option("--atomname-mapping", Settings::atomName_mapping_fileName, "atom name mapping file name")
        ->envname("KENREF_ATOMNAME_MAPPING"); // ->check(CLI::ExistingFile);

    app.add_option("-z,--proton-mhz", Settings::proton_mhz, "spectrometer proton field strength in MHz")->envname("KENREF_PROTON_MHZ");
    // app.add_option("--max-atom-to-read", Settings::max_atom_to_read, "maximum number of atoms to read")->envname("KENREF_MAX_ATOM_TO_READ");
    app.add_option("--alt-out-path", Settings::alt_out_path, "path to redirect output stream")->envname("KENREF_ALT_OUT_PATH");


    // Load from config file
    app.set_config("--params", "params.toml", "Read a TOML config file", false)->envname("KENREF_PARAMS");
    app.allow_extras();

    app.usage("Kinetic Ensemble Refinement library.");
    app.footer(
        "\nUse '--' to separate KEnRef options from GROMACS options.\n"
        "This is very useful if KEnRef parameter(s) colloides with Gromacs parameter(s)\n"
        "Example: KEnRef -g CA -d index.ndx -- -deffnm simulation\n"
        "or when you want to pass a parameter to Gromacs that could be consumed by KEnRef otherwise\n"
        "Example: KEnRef -m SIGMA -- -h");
    CLI11_PARSE(app, argc, argv);

    if (Settings::debug) {
        volatile bool holdToDebug = Settings::debug;
        while (/*simulationIndex > 0 &&*/ holdToDebug) {
            sleep(1);
        }
    }

    if (Settings::refFileName.empty()) {
        std::cout << "Reference file parameter (--ref) not found, using atomname-mapping (" << Settings::atomName_mapping_fileName << ")" << std::endl;
        Settings::refFileName = Settings::atomName_mapping_fileName;
    }else {
        std::cout << "Reference file parameter (--ref) set to , (" << Settings::refFileName << ")" << std::endl;
    }

    //============================================================
    auto remaining = app.remaining();
    for (auto i = remaining.begin(); i != remaining.end(); ++i) {
        if (std::string(*i) == "--") {
            remaining.erase(i);
            break;
        }
    }
    auto new_argv = create_argv(argc, argv, remaining);
    int new_argc = static_cast<int>(remaining.size() + 1); // +1 for program name

    //next line is copied from mdrun_main.cpp
    return gmx::CommandLineModuleManager::runAsMainCMainWithSettings(
        new_argc, new_argv.get(), &gmx::/*KEnRef_*/gmx_mdrun, &initSettingsNoNice);
}
