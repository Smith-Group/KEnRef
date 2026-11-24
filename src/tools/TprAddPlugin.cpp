/*
 * TprAddPlugin.cpp
 *
 * Test tool: writes a copy of a tpr that names an MDModules plugin, as grompp writes it for
 * "mdmodule-plugins = <name>" in the mdp. plugin_validate.sh uses it to run rt_validate's tprs with
 * the plugin, since the GB3 fixture has no topology to run grompp on.
 *
 *     kenref_tpr_add_plugin <in.tpr> <out.tpr> <plugin name>
 */
#include <iostream>
#include <memory>
#include <utility>

#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"

int main(int argc, char *argv[]) {
    if (argc != 4) {
        std::cerr << "usage: " << argv[0] << " <in.tpr> <out.tpr> <plugin name>" << std::endl;
        return 2;
    }
    // The key under which grompp records the plugins (MDModules::writePluginRecord)
    const char *const recordKey = "mdmodule-plugins";

    t_inputrec ir;
    t_state state;
    gmx_mtop_t mtop;
    read_tpx_state(argv[1], &ir, &state, &mtop);

    gmx::KeyValueTreeBuilder builder;
    gmx::KeyValueTreeObjectBuilder root = builder.rootObject();
    if (ir.internalParameters != nullptr) {
        if (ir.internalParameters->keyExists(recordKey)) {
            std::cerr << argv[1] << " already names plugins" << std::endl;
            return 1;
        }
        root.mergeObject(std::move(*ir.internalParameters));
    }
    // mdrun reads only the names; the plugin is checked when it is loaded
    root.addObject(recordKey).addObject(argv[3]);
    ir.internalParameters = std::make_unique<gmx::KeyValueTreeObject>(builder.build());

    write_tpx_state(argv[2], &ir, &state, mtop);
    return 0;
}
