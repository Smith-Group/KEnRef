/*
 * EnergyModel.h
 *
 *  The per-model abstraction at the heart of the restructure. A concrete EnergyModel (SigmaModel,
 *  PlateausModel, RelaxModel, …) is ONE self-contained unit declaring the three things that vary
 *  between models — and nothing else:
 *    (1) INPUTS it consumes               -> loaded in buildCache()
 *    (2) PARAM SPECIFICATION              -> the static ParamSchema (registered with the model)
 *    (3) CACHED values reused across steps -> built in buildCache()/finalizeIndexing(), held as
 *                                             members, consumed by compute()
 *
 *  Everything else (param binding, the per-step pipeline, engine glue) is shared. Each consumer
 *  holds a unique_ptr<EnergyModel> obtained from ModelRegistry by the user-selected name.
 *
 *  Init is two-phase to resolve the chicken-and-egg in the current consumers (the atom-name pairs
 *  come FROM the model's data files, but the sub-atom indexing that turns names into compact ids
 *  is built by the driver from those pairs, and the engine must request atoms before computing):
 *     1. buildCache(InitContext)      : load data files; expose atomNamePairs().
 *     2. (driver) build sub-indexing from atomNamePairs(); request atoms from the engine.
 *     3. finalizeIndexing(IndexingContext) : convert names->sub0Id pairs; prime per-data caches.
 */

#ifndef KENREF_ENERGYMODEL_H_
#define KENREF_ENERGYMODEL_H_

#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "KEnRef.h"          // CoordsMatrixType<>, KEnRef_Real_t, KEnRef<>
#include "ParamSchema.h"
#include "EngineAdapter.h"

namespace kenref {

/**
 * Resolved-parameter + init-data view passed to buildCache(). Typed getters read the user's value
 * through the EngineAdapter and fall back to the schema default; the conversion mirrors how the old
 * consumers parsed each keyword. Non-owning references — valid only for the duration of buildCache().
 */
template<typename Real>
struct InitContext {
    const EngineAdapter<Real>& adapter;
    const ParamSchema&         schema;            // merged General + this model's Model-tier specs
    const std::map<std::string, int>& atomNameToGlobalId; // normalised name -> 1-based serial
    bool handleNames = false;                     // IoUtils::should_handleNames(...) result
    int  numModelsTotal = 1;                      // for model-count-dependent setup (e.g. grouping)

    [[nodiscard]] std::optional<std::string> raw(const std::string& key) const {
        if (auto v = adapter.getRawParam(key)) return v;
        for (const auto& s : schema.specs())
            if (s.key == key) return s.defaultValue;  // may be nullopt
        return std::nullopt;
    }
    [[nodiscard]] std::string getString(const std::string& key, const std::string& fallback = "") const {
        auto v = raw(key); return v ? *v : fallback;
    }
    [[nodiscard]] Real getReal(const std::string& key, Real fallback = Real(0)) const {
        auto v = raw(key); if (!v) return fallback; std::istringstream is(*v); Real r; is >> r; return r;
    }
    [[nodiscard]] int getInt(const std::string& key, int fallback = 0) const {
        auto v = raw(key); if (!v) return fallback; std::istringstream is(*v); int r; is >> r; return r;
    }
    [[nodiscard]] bool getFlag(const std::string& key, bool fallback = false) const {
        auto v = raw(key); if (!v) return fallback;
        return *v == "1" || *v == "true" || *v == "TRUE" || *v == "on" || *v == "ON";
    }
};

/**
 * View passed to finalizeIndexing() once the driver has built the sub-atom indexing from the
 * model's atomNamePairs(). The model uses it to convert its name pairs to compact (sub0Id) ids and
 * prime any per-data caches (e.g. SpecDenData's atomId-pair cache).
 */
template<typename Real>
struct IndexingContext {
    const std::map<std::string, int>& atomNameToSub0Id; // normalised name -> 0-based compact id
};

/**
 * Per-step inputs to compute(). `coord_array` is the assembled, already-fitted coordinates of ALL
 * models (one [nSub x 3] per model), valid on the master only. It is NON-const because some kernels
 * (SIGMA, RELAX) scale it in place (Å→m); the driver passes a fresh per-step buffer so that is safe.
 * `k`/`n` are the general force-constant / power params (framework-owned). Return matches the
 * existing kernels: (energy, optional per-model derivatives [nSub x 3]).
 */
template<typename Real>
struct StepContext {
    std::vector<CoordsMatrixType<Real>>& coord_array;
    Real k = Real(1.0);
    Real n = Real(0.25);
    bool gradient = false;
    int  numOmpThreads = 0;
};

template<typename Real>
class EnergyModel {
public:
    virtual ~EnergyModel() = default;

    // (1)+(3): load input data files and build the between-step caches.
    virtual void buildCache(const InitContext<Real>& ctx) = 0;

    // The model's atom-name pairs (available after buildCache); the driver builds sub-indexing from
    // their union, and S2OrderParamsCalculator reads pairs from the selected model.
    [[nodiscard]] virtual const std::vector<std::tuple<std::string, std::string>>&
        atomNamePairs() const = 0;

    // (3 cont.): convert name pairs -> sub0Id pairs and prime per-data caches.
    virtual void finalizeIndexing(const IndexingContext<Real>& ctx) = 0;

    // Compact (sub0Id) atom-id pairs (available after finalizeIndexing); used by S2 and by models
    // (e.g. PLATEAUS) that pass atomId_pairs straight to their kernel.
    [[nodiscard]] virtual const std::vector<std::tuple<int, int>>& atomIdPairs() const = 0;

    // The per-step compute — forwards to an existing coord_array_to_*_energy kernel (unchanged).
    virtual std::tuple<Real, std::optional<std::vector<CoordsMatrixType<Real>>>>
        compute(const StepContext<Real>& ctx) = 0;

    // Force/derivative unit scale applied by the driver after inverse-fit (replaces the per-model
    // "*= 10" branch): SIGMA/RELAX = 10 (Å⁻¹ -> nm⁻¹), PLATEAUS = 1 (manuscript back-compat).
    [[nodiscard]] virtual Real forceUnitScale() const = 0;
};

} // namespace kenref

#endif // KENREF_ENERGYMODEL_H_
