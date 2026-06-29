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
#include <stdexcept>

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

    // Is a spec with this key already present?
    [[nodiscard]] bool has(const std::string& key) const { return findSpec(key) != nullptr; }

    // Append every spec of `other`, with strict conflict detection (used to build the General + Model
    // tier per model, and the union of all models' Model-tier specs for keyword registration):
    //   - key not present            -> add it
    //   - key present, SAME spec     -> skip (idempotent dedup; never produces a duplicate keyword,
    //                                    which e.g. PLUMED registerKeywords would reject)
    //   - key present, DIFFERENT spec -> throw std::invalid_argument (model-tier keyword names must be
    //                                    globally unique, so this is a caller bug, not a silent merge)
    // "Same spec" compares type/tier/required/default; the doc/help text is ignored as benign.
    // (Schemas are tiny, so the O(n*m) scan is irrelevant.)
    ParamSchema& merge(const ParamSchema& other) {
        for (const auto& s : other.specs_) {
            const ParamSpec* existing = findSpec(s.key);
            if (existing == nullptr) {
                specs_.push_back(s);
            } else if (!sameDeclaration(*existing, s)) {
                throw std::invalid_argument("ParamSchema::merge: conflicting declarations for key '" + s.key + "'");
            }
            // else: identical declaration already present -> dedup (skip)
        }
        return *this;
    }

private:
    [[nodiscard]] const ParamSpec* findSpec(const std::string& key) const {
        for (const auto& s : specs_) if (s.key == key) return &s;
        return nullptr;
    }

    // Two specs for the same key are "the same declaration" iff their semantically meaningful fields
    // match. doc (help text) is intentionally excluded — differing help text for the same parameter is
    // benign and should not be treated as a conflict.
    static bool sameDeclaration(const ParamSpec& a, const ParamSpec& b) {
        return a.type == b.type && a.tier == b.tier && a.required == b.required
               && a.defaultValue == b.defaultValue;
    }

    std::vector<ParamSpec> specs_;
};

} // namespace kenref

#endif // KENREF_PARAMSCHEMA_H_
