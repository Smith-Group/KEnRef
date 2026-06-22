/*
 * ModelRegistry.h
 *
 *  Name -> {factory, schema} table so the user selects an energy model at runtime by string, and
 *  adding a model never edits a consumer. Selection is two-stage (see the restructure design):
 *    - link time : bootstrapModels() registers every COMPILE-TIME-ENABLED model (KENREF_ENABLE_*).
 *    - run  time : create(name) instantiates the ONE model the user picked.
 *
 *  Linking (Style 2): bootstrapModels() is defined in a CMake-generated TU that calls each model's
 *  register*Model() free function under `#if KENREF_ENABLE_<MODEL>`. Because the driver references
 *  bootstrapModels(), the linker keeps those TUs — no `--whole-archive`, works for .a and .so, and
 *  `-DKENREF_ENABLE_<MODEL>=OFF` drops a model at compile time.
 */

#ifndef KENREF_MODELREGISTRY_H_
#define KENREF_MODELREGISTRY_H_

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "EnergyModel.h"
#include "ParamSchema.h"

namespace kenref {

template<typename Real>
class ModelRegistry {
public:
    using Factory  = std::function<std::unique_ptr<EnergyModel<Real>>()>;
    using SchemaFn = ParamSchema (*)();   // the model's Model-tier ParamSchema (no instance needed)

    struct Entry { Factory factory; SchemaFn schemaFn; };

    // Register a model under `name` (case-sensitive, e.g. "SIGMA"). Called from each model's
    // register*Model() free function, which bootstrapModels() invokes under #if KENREF_ENABLE_*.
    static void add(const std::string& name, Factory factory, SchemaFn schemaFn) {
        table()[name] = Entry{std::move(factory), schemaFn};
    }

    // Instantiate the runtime-selected model, or nullptr if `name` is unknown/disabled.
    [[nodiscard]] static std::unique_ptr<EnergyModel<Real>> create(const std::string& name) {
        const auto it = table().find(name);
        return it == table().end() ? nullptr : it->second.factory();
    }

    [[nodiscard]] static bool has(const std::string& name) { return table().count(name) != 0; }

    // Model-tier schema of every registered model, paired with its name (for union keyword
    // registration on PLUMED / CLI11 option registration on GROMACS, and for diagnostics).
    [[nodiscard]] static std::vector<std::pair<std::string, ParamSchema>> allSchemas() {
        std::vector<std::pair<std::string, ParamSchema>> out;
        out.reserve(table().size());
        for (const auto& [name, e] : table()) out.emplace_back(name, e.schemaFn());
        return out;
    }

    [[nodiscard]] static std::vector<std::string> names() {
        std::vector<std::string> out;
        out.reserve(table().size());
        for (const auto& [name, e] : table()) out.push_back(name);
        return out;
    }

private:
    // Function-static singleton: avoids static-init-order issues across translation units.
    static std::map<std::string, Entry>& table() {
        static std::map<std::string, Entry> t;
        return t;
    }
};

/**
 * Register every compile-time-enabled model. DEFINED in the CMake-generated ModelBootstrap.cpp,
 * which calls each register*Model() under `#if KENREF_ENABLE_<MODEL>`. KEnRefDriver calls this once
 * before first use. Non-template (operates on ModelRegistry<KEnRef_Real_t>) since the build fixes
 * the scalar type via KENREF_DOUBLE.
 */
void bootstrapModels();

} // namespace kenref

#endif // KENREF_MODELREGISTRY_H_
