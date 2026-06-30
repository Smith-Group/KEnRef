/*
 * buildModelIndexing.h
 *
 *  ROLE: the engine-agnostic replacement for the SHARED half of every consumer's "fillParamsStep0" —
 *  the sub-atom indexing that GROMACS KEnRefForceProvider, PLUMED KEnRefBias, and the offline tools
 *  used to each duplicate. It creates the named EnergyModel, loads its inputs (buildCache via the
 *  adapter), derives the compact sub-atom indexing from the model's atom-name pairs, and finalizes the
 *  model's indexing — then hands the prepared model back.
 *
 *  It deliberately does NOT construct the KEnRefDriver: driver construction needs consumer-specific
 *  values (k / n / maxForceSquared / the centered reference) and not every consumer wants one (s2calc
 *  is a metric, not an energy model — it only needs the model's atomIdPairs()). So each consumer does
 *  the 1-line `make_unique<KEnRefDriver>(std::move(mi.model), ...)` itself; that line is not duplication.
 *
 *  The engine supplies its atom-name -> global-id map (ANY consistent base: 1-based serials for the live
 *  engines, or 0-based for the tools — the returned subAtomGlobalIds are in that SAME base) and then
 *  builds its own per-step position-fetch mapping from subAtomGlobalIds (GROMACS ga2la, PLUMED
 *  serial->local, the tools' trajectory index).
 *
 *  Header-only (it only orchestrates other core calls). No engine headers — pure kenref_core.
 */

#ifndef KENREF_BUILDMODELINDEXING_H_
#define KENREF_BUILDMODELINDEXING_H_

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "EnergyModel.h"       // EnergyModel<>, InitContext / IndexingContext, ParamSchema
#include "EngineAdapter.h"
#include "ModelRegistry.h"

namespace kenref {

template<typename Real>
struct ModelIndexing {
    std::unique_ptr<EnergyModel<Real>> model;         // built + finalized; ready for compute()/atomIdPairs()
    std::vector<int> subAtomGlobalIds;                // sorted-unique global-ids (engine's base) of the
                                                      // restrained sub-atoms; sub0Id == index into this
    std::map<std::string, int> atomName_to_sub0Id;    // normalised name -> compact 0-based sub-id
};

/**
 * Build the model + sub-indexing. The caller must have validated `modelName` against
 * ModelRegistry::has() and reported an engine-appropriate error; here `create` is assumed to succeed.
 */
template<typename Real>
ModelIndexing<Real> buildModelIndexing(
    const std::string& modelName,
    const EngineAdapter<Real>& adapter,
    const std::map<std::string, int>& atomName_to_globalId,
    bool handleNames, int numModelsTotal, int numOmpThreads)
{
    bootstrapModels();
    ModelIndexing<Real> out;
    out.model = ModelRegistry<Real>::create(modelName);

    // (1)+(3) load the model's input data files + build its between-step caches.
    {
        ParamSchema emptySchema; // the adapter supplies every param the model reads via getRawParam
        InitContext<Real> initCtx{adapter, emptySchema, atomName_to_globalId, handleNames, numModelsTotal};
        out.model->buildCache(initCtx);
    }

    // Sub-atom set = the sorted-unique global-ids appearing in any atom-name pair (std::set gives both
    // uniqueness and ascending order, matching every consumer's increasing-global-id sub0Id assignment).
    {
        std::set<int> ids;
        for (const auto& [a1, a2] : out.model->atomNamePairs()) {
            ids.insert(atomName_to_globalId.at(a1)); // .at(): throw on an unexpected name
            ids.insert(atomName_to_globalId.at(a2));
        }
        out.subAtomGlobalIds.assign(ids.begin(), ids.end());
    }
    std::map<int, int> globalId_to_sub0Id;
    for (int sub0Id = 0; sub0Id < static_cast<int>(out.subAtomGlobalIds.size()); ++sub0Id)
        globalId_to_sub0Id[out.subAtomGlobalIds[sub0Id]] = sub0Id;
    for (const auto& [name, gid] : atomName_to_globalId) {
        const auto it = globalId_to_sub0Id.find(gid);
        if (it != globalId_to_sub0Id.end())
            out.atomName_to_sub0Id[name] = it->second;
    }

    // (3 cont.) convert the model's name pairs -> sub0Id pairs + prime its per-data atomId caches.
    {
        IndexingContext<Real> idxCtx{out.atomName_to_sub0Id, numOmpThreads};
        out.model->finalizeIndexing(idxCtx);
    }

    return out;
}

} // namespace kenref

#endif // KENREF_BUILDMODELINDEXING_H_
