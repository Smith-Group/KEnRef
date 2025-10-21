/*
 * gmxwrapper.cpp
 *
 *  Created on: Nov 29, 2022
 *      Author: amr
 */
#include <string>
#include <map>
#include <CLI11/CLI11.hpp>
#include "gmxpre.h"
#include "mdrun/mdrun_main.h"
#include "gromacs/commandline/cmdlinemodule.h"
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

    // KENREF_MAXFORCE=999
    app.add_option("--k", Settings::k, "K force constant");
    app.add_option("--n", Settings::n, "power scaling");
    app.add_option("--max-force", Settings::max_force, "maximum force");

    app.add_option("-g,--guide",  Settings::GUIDE_C_ALPHA, "name of guide group");
    app.add_option("-d,--index", Settings::indexFileName, "Index file name");
    app.add_option("-r,--ref", Settings::refFileName, "Reference file");

    app.add_option("-m,--model", Settings::selected_energy_model, "energy model")
            ->required()
            ->transform(CLI::CheckedTransformer(KEnRef<KEnRef_Real_t>::energyModelMap, CLI::ignore_case));
    app.add_option("-x,--exp-data-folder", Settings::experimentalDataFolder, "experimental data folder for sigma")
        ->check(CLI::ExistingDirectory);
    app.add_option("-X,--exp-data-file", Settings::experimentalDataFileName, "experimental data file for plateaus")
        ->check(CLI::ExistingFile);

    app.add_option("--atomname-mapping", Settings::atomName_mapping_fileName, "atom name mapping file name")
        ->check(CLI::ExistingFile);

    app.add_option("-z,--proton-mhz", Settings::proton_mhz, "spectrometer proton field strength in MHz");
    app.add_option("--max-atom-to-read", Settings::max_atom_to_read, "maximum number of atoms to read");
    app.add_option("--alt-out-path", Settings::alt_out_path, "path to redirect output stream");


    // Load from config file
    app.set_config("--params", "params.toml", "Read a TOML config file", false);
    app.allow_extras();
    CLI11_PARSE(app, argc, argv);

    if (Settings::refFileName.empty()) {
        std::cout << "Reference file (--ref) not found, using ATOMNAME_MAPPING_FILENAME (" << Settings::atomName_mapping_fileName << ")" << std::endl;
        Settings::refFileName = Settings::atomName_mapping_fileName;
    }else {
        std::cout << "Reference file (--ref) found, (" << Settings::refFileName << ")" << std::endl;
    }

    //============================================================
    const auto remaining = app.remaining();
    auto new_argv = create_argv(argc, argv, remaining);
    int new_argc = static_cast<int>(remaining.size() + 1); // +1 for program name

    //next line is copied from mdrun_main.cpp
    return gmx::CommandLineModuleManager::runAsMainCMainWithSettings(
        new_argc, new_argv.get(), &gmx::/*KEnRef_*/gmx_mdrun, &initSettingsNoNice);
}
