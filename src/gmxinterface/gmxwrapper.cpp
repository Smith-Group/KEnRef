/*
 * gmxwrapper.cpp
 *
 *  Created on: Nov 29, 2022
 *      Author: amr
 */
#include <string>
#include <iostream>
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
#include "gromacs/utility/baseversion.h"
#include "gmxinterface/gmxwrapper.h"
#include "gmxinterface/KEnRefGmxVersion.h"
#include "gmxinterface/KEnRefSettingsOptions.h"
#include "core/KEnRef.h"
#include "core/ModelRegistry.h"
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
    // Refuse to run inside a GROMACS other than the one KEnRef was compiled against: the internal headers we
    // use have no ABI promise, and a mismatch does not crash -- it misroutes notifier subscriptions without
    // a word. First thing, before anything reaches into libgromacs. See KEnRefGmxVersion.h.
    if (const long loaded = kenref::gmxVersionNumber(gmx_version()); loaded != GMX_VERSION) {
        std::cerr << "KEnRef: compiled against GROMACS " << GMX_VERSION / 10000 << "." << GMX_VERSION % 10000
                  << ", but the libgromacs loaded at run time reports '" << gmx_version() << "'.\n"
                  << "  KEnRef uses GROMACS internals that change between releases, so it must run with the "
                     "GROMACS it was built against.\n"
                  << "  Check LD_LIBRARY_PATH and any loaded GROMACS module, or rebuild KEnRef against this "
                     "GROMACS." << std::endl;
        return 1;
    }
    // This should be the main function that loads GMXWrapper (which should replace GROMACS main class)
    // Register every compile-time-enabled energy model NOW so --model can be validated against the
    // registry at CLI-parse time (buildModelIndexing bootstraps again later; the call is idempotent).
    kenref::bootstrapModels();
    CLI::App app{"KEnRef"};

    kenref::addSettingsOptions(&app);
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

    kenref::useAtomNameMappingAsDefaultReference();

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
