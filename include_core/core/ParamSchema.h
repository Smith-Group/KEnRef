/*
 * ParamSchema.h
 *
 *  Declarative parameter specification for KEnRef energy models (restructure concern #2).
 *
 *  Each model declares WHAT configuration it needs (key, value-kind, default, required-ness)
 *  as plain data; each engine binds that schema to its own parameter mechanism — PLUMED native
 *  keywords (registerKeywords + parse), or the GROMACS CLI11/TOML Settings. The schema is the
 *  single source of truth so that adding a model never edits a consumer's parameter code.
 */

#ifndef KENREF_PARAMSCHEMA_H_
#define KENREF_PARAMSCHEMA_H_

#include <string>
#include <vector>
#include <optional>

namespace kenref {

/**
 * Which engine-level tier a parameter belongs to.
 *  - General : framework-owned, identical for every model (e.g. K, N, MAX_FORCE). On PLUMED these
 *              are registered as `compulsory` keywords.
 *  - Model   : declared by one specific EnergyModel; only meaningful when that model is selected.
 *              On PLUMED these are registered as `optional` keywords (since registerKeywords runs
 *              before MODEL is parsed, the union of all models' Model-tier keys is registered), and
 *              their required-ness is enforced per-selected-model in EnergyModel::buildCache().
 */
enum class ParamTier { General, Model };

/**
 * The value-kind of a parameter, so each engine can route it to its own typed parse call
 * (PLUMED parse()/parseFlag(); GROMACS CLI11 add_option/add_flag). `Path` is a String that
 * additionally denotes a filesystem path (an engine may add existence checks).
 */
enum class ParamType { Real, Int, String, Bool, Path };

/**
 * One declared parameter. `defaultValue` is kept as a string (parsed by the consumer to the
 * ParamType) so the schema stays a pure-data description independent of any engine's typed parser;
 * std::nullopt means "no default". `doc` becomes the keyword help text (e.g. PLUMED registration).
 */
struct ParamSpec {
    std::string key;
    ParamType   type;
    ParamTier   tier;
    bool        required = false;
    std::optional<std::string> defaultValue = std::nullopt;
    std::string doc;
};

/**
 * An ordered collection of ParamSpecs (a model's, or the framework's General set). Ordering is
 * preserved so engines can register keywords deterministically.
 */
class ParamSchema {
public:
    ParamSchema& add(ParamSpec spec) { specs_.push_back(std::move(spec)); return *this; }

    [[nodiscard]] const std::vector<ParamSpec>& specs() const { return specs_; }

    // Convenience: append every spec of another schema (used to merge General + Model tiers).
    ParamSchema& merge(const ParamSchema& other) {
        for (const auto& s : other.specs_) specs_.push_back(s);
        return *this;
    }

private:
    std::vector<ParamSpec> specs_;
};

} // namespace kenref

#endif // KENREF_PARAMSCHEMA_H_
