/*
 * PlateausModel.cpp
 *
 *  PLATEAUS (g-value) energy model wrapped as a registry EnergyModel. Concerns it owns:
 *    (1) inputs   : the experimental table (EXP_DATA_FILE) giving atom-name pairs + target g values;
 *                   the model-count-dependent grouping list.
 *    (2) params   : the Model-tier ParamSchema below.
 *    (3) caches   : g0_, grouping_, atom-name pairs and the compact atomId pairs.
 *  compute() forwards to the unchanged KEnRef::coord_array_to_g_energy kernel.
 *
 *  Faithful to BOTH production consumers (KEnRefForceProvider::fillParamsStep0 PLATEAUS branch and
 *  KEnRefBias): pairs come from atom1/atom2 NAMES (normalised, with the per-table "unprepared" check)
 *  and g0 from g1/g2 — NOT the i1/i2 indices the EnergyCalculator tool used (unified away here).
 */

#include "core/EnergyModel.h"
#include "core/ModelRegistry.h"
#include "core/KEnRef.h"
#include "core/IoUtils.h"
#include "core/Table.h"

namespace kenref {

namespace {

class PlateausModel : public EnergyModel<KEnRef_Real_t> {
public:
    static ParamSchema schema() {
        ParamSchema s;
        s.add({"EXP_DATA_FILE", ParamType::Path, ParamTier::Model, /*required*/ true, std::nullopt,
               "Experimental data file with atom1/atom2 pairs and g1/g2 targets (PLATEAUS model)"});
        return s;
    }

    void buildCache(const InitContext<KEnRef_Real_t>& ctx) override {
        const std::string exp_data_file = ctx.getString("EXP_DATA_FILE");
        const Table table = IoUtils::readTable(exp_data_file, true, false, "\\s*,\\s*", -1);

        // Mirror the consumers: confirm whether we must handle unprepared atom names (per-table).
        bool handleUnpreparedAtomNames = IoUtils::should_handleNames(table);

        atomName_pairs_.clear();
        for (int row = 0; row < table.rowCount(); row++) {
            atomName_pairs_.emplace_back(
                IoUtils::normalizeName(table(row, "atom1"), handleUnpreparedAtomNames),
                IoUtils::normalizeName(table(row, "atom2"), handleUnpreparedAtomNames));
        }

        // g0_ (numPairs x 2) from the g1/g2 columns.
        g0_.resize(static_cast<Eigen::Index>(table.rowCount()), 2);
        for (int i = 0; i < table.rowCount(); ++i) {
            std::istringstream t1(table(i, "g1")), t2(table(i, "g2"));
            t1 >> g0_(i, 0);
            t2 >> g0_(i, 1);
        }

        // Grouping list keyed on the total number of models (the consumers' hardcoded switch).
        grouping_ = groupingForModelCount(ctx.numModelsTotal);
    }

    [[nodiscard]] const std::vector<std::tuple<std::string, std::string>>& atomNamePairs() const override {
        return atomName_pairs_;
    }

    void finalizeIndexing(const IndexingContext<KEnRef_Real_t>& ctx) override {
        atomId_pairs_ = IoUtils::atomNamePairs_2_atomIdPairs(atomName_pairs_, ctx.atomNameToSub0Id, ctx.numOmpThreads);
    }

    [[nodiscard]] const std::vector<std::tuple<int, int>>& atomIdPairs() const override {
        return atomId_pairs_;
    }

    std::tuple<KEnRef_Real_t, std::optional<std::vector<CoordsMatrixType<KEnRef_Real_t>>>>
    compute(const StepContext<KEnRef_Real_t>& ctx) override {
        return KEnRef<KEnRef_Real_t>::coord_array_to_g_energy(
            ctx.coord_array, atomId_pairs_, grouping_, g0_, ctx.k, ctx.n, ctx.gradient, ctx.numOmpThreads);
    }

    [[nodiscard]] KEnRef_Real_t forceUnitScale() const override {
        // PLATEAUS keeps scale 1 (manuscript back-compat; the consumers' "*= 10 unless PLATEAUS").
        // TODO Consider changing it to 10.0 in next releases
        return KEnRef_Real_t(1.0);
    }

private:
    static std::vector<std::vector<std::vector<int>>> groupingForModelCount(int numModels) {
        switch (numModels) {
            case 1:  return {{{0}}, {{0}}};
            case 2:  return {{{0, 1}}, {{0}, {1}}};
            case 3:  return {{{0, 1, 2}}, {{0}, {1}, {2}}};
            default:
                throw std::runtime_error("PLATEAUS: I don't know how to handle more than 3 simulations yet");
        }
    }

    Eigen::MatrixX<KEnRef_Real_t> g0_;
    std::vector<std::vector<std::vector<int>>> grouping_;
    std::vector<std::tuple<std::string, std::string>> atomName_pairs_;
    std::vector<std::tuple<int, int>> atomId_pairs_;
};

} // namespace

// Called by the generated ModelBootstrap::bootstrapModels() under #if KENREF_ENABLE_PLATEAUS.
void registerPlateausModel() {
    ModelRegistry<KEnRef_Real_t>::add(
        "PLATEAUS",
        [] { return std::make_unique<PlateausModel>(); },
        &PlateausModel::schema);
}

} // namespace kenref
