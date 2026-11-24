/*
 * KEnRefPlugin.cpp
 *
 * KEnRef as an MDModules plugin, which a GROMACS built with GMX_MDMODULE_PLUGINS=ON loads when the
 * mdp names it:
 *
 *     mdmodule-plugins = kenref
 *
 * GROMACS looks for kenref.so in GMX_MDMODULE_PLUGIN_PATH. The plugin reads its settings from the
 * same KENREF_* environment variables, and KENREF_PARAMS, as bin/KEnRef.
 */
#include <memory>
#include <stdexcept>
#include <string>

#include "gromacs/mdtypes/mdmodulepluginentrypoints.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/version.h"
#include "gmxinterface/KEnRefGmxVersion.h"
#include "gmxinterface/KEnRefMDModule.h"

namespace {

std::unique_ptr<gmx::IMDModule> createKEnRefModule() {
    // KEnRef uses GROMACS headers that the plugin interface does not version, so it must run in the
    // GROMACS it was built against, as bin/KEnRef checks too
    if (kenref::gmxVersionNumber(gmx_version()) != GMX_VERSION) {
        throw std::runtime_error(std::string("KEnRef was built against GROMACS ") + std::to_string(GMX_VERSION / 10000)
                                 + "." + std::to_string(GMX_VERSION % 10000) + ", but runs in '" + gmx_version()
                                 + "'. Rebuild KEnRef against this GROMACS.");
    }
    return std::make_unique<KEnRefMDModule>(true);
}

} // namespace

GMX_MDMODULE_PLUGIN("kenref", KENREF_PLUGIN_VERSION, createKEnRefModule)
