/*
 * SigmaModel.cpp
 *
 *  SIGMA energy model wrapped as a registry EnergyModel. Concerns it owns:
 *    (1) inputs   : spectral-density data (load_spec_den_data from EXP_DATA_FOLDER), relaxation
 *                   rates (default, or from RATES_FILE), proton MHz.
 *    (2) params   : the Model-tier ParamSchema below.
 *    (3) caches   : spec_den_data_list_ + its primed atomId-pair caches; the name->id map and the
 *                   compact atomId pairs.
 *  compute() forwards to the unchanged KEnRef::coord_array_to_sigma_energy kernel.
 *
 *  The RATES_FILE override mirrors Step 0's KEnRefBias change, now shared by every engine.
 */

#include "core/EnergyModel.h"
#include "core/ModelRegistry.h"
#include "core/KEnRef.h"
#include "core/IoUtils.h"
#include "core/Table.h"

namespace kenref {

namespace {

class SigmaModel : public EnergyModel<KEnRef_Real_t> {
public:
    static ParamSchema schema() {
        ParamSchema s;
        s.add({"EXP_DATA_FOLDER", ParamType::Path,   ParamTier::Model, /*required*/ true,  std::nullopt,
               "Folder with spectral-density experimental data (SIGMA model)"});
        s.add({"PROTON_MHZ",      ParamType::Real,   ParamTier::Model, /*required*/ false, std::string("700.0"),
               "Spectrometer proton field strength in MHz (SIGMA model)"});
        s.add({"RATES_FILE",      ParamType::Path,   ParamTier::Model, /*required*/ false, std::nullopt,
               "CSV file (a_coef.csv style: header row of rate names + one data row of values) "
               "overriding the default SIGMA relaxation rates"});
        return s;
    }

    void buildCache(const InitContext<KEnRef_Real_t>& ctx) override {
        proton_mhz_      = ctx.getReal("PROTON_MHZ", KEnRef_Real_t(700.0));
        exp_data_folder_ = ctx.getString("EXP_DATA_FOLDER");

        spec_den_data_list_ = IoUtils::load_spec_den_data(exp_data_folder_, ctx.handleNames);

        // Parameter-validation: num_models must match the ensemble-member count the exp-data was built
        // for (InitContext::numModelsTotal is provided for exactly this model-count-dependent setup).
        // Done once here, NOT per-step in the hot coord_array_to_sigma_energy kernel.
        IoUtils::assertEnsembleMemberCount(spec_den_data_list_, ctx.numModelsTotal, "SIGMA");

        // Optionally override the default rates_ from a CSV file (header = rate names, one data row
        // of values) read as a 1-row table -> NamedRowVector; otherwise keep the built-in default.
        const std::string rates_file = ctx.getString("RATES_FILE");
        if (!rates_file.empty()) {
            rates_ = IoUtils::readTable(rates_file, true, false).toNamedRowVector<KEnRef_Real_t>();
        } else {
            rates_ = Table({{"5.0e+08", "2.5e+08", "1.0e+12", "1.0e+04"}},
                           {{"kens", "kc", "kmethyl", "karo"}}).toNamedRowVector<KEnRef_Real_t>();
        }

        // Collect all atom-name pairs across every SpecDenData object (the driver builds sub-indexing
        // from their union).
        atomName_pairs_.clear();
        for (const auto& sdd : spec_den_data_list_)
            for (const auto& ap : sdd.get_atom_pairs())
                atomName_pairs_.emplace_back(ap);
    }

    [[nodiscard]] const std::vector<std::tuple<std::string, std::string>>& atomNamePairs() const override {
        return atomName_pairs_;
    }

    void finalizeIndexing(const IndexingContext<KEnRef_Real_t>& ctx) override {
        atomNames_2_atomIds_ = ctx.atomNameToSub0Id; // the kernel maps SpecDenData name pairs through this
        atomId_pairs_ = IoUtils::atomNamePairs_2_atomIdPairs(atomName_pairs_, atomNames_2_atomIds_, ctx.numOmpThreads);

        // Prime each SpecDenData's atomId-pair cache (the SIGMA-only init the consumers did by hand).
        for (auto& sdd : spec_den_data_list_)
            sdd.set_atomIdPairs_to_sub0Atom_id_pairs_cache(
                IoUtils::atomNamePairs_2_atomIdPairs(sdd.get_atom_pairs(), atomNames_2_atomIds_, ctx.numOmpThreads));
    }

    [[nodiscard]] const std::vector<std::tuple<int, int>>& atomIdPairs() const override {
        return atomId_pairs_;
    }

    std::tuple<KEnRef_Real_t, std::optional<std::vector<CoordsMatrixType<KEnRef_Real_t>>>>
    compute(const StepContext<KEnRef_Real_t>& ctx) override {
        return KEnRef<KEnRef_Real_t>::coord_array_to_sigma_energy(
            ctx.coord_array, rates_, spec_den_data_list_, proton_mhz_,
            ctx.k, ctx.n, atomNames_2_atomIds_, ctx.gradient, ctx.numOmpThreads);
    }

    [[nodiscard]] KEnRef_Real_t forceUnitScale() const override { return static_cast<KEnRef_Real_t>(10.0); }

private:
    KEnRef_Real_t proton_mhz_ = static_cast<KEnRef_Real_t>(700.0);
    std::string   exp_data_folder_;
    NamedRowVector<KEnRef_Real_t> rates_;
    std::vector<SpecDenData<KEnRef_Real_t>> spec_den_data_list_;
    std::vector<std::tuple<std::string, std::string>> atomName_pairs_;
    std::vector<std::tuple<int, int>> atomId_pairs_;
    std::map<std::string, int> atomNames_2_atomIds_;
};

} // namespace

// Called by the generated ModelBootstrap::bootstrapModels() under #if KENREF_ENABLE_SIGMA.
void registerSigmaModel() {
    ModelRegistry<KEnRef_Real_t>::add(
        "SIGMA",
        [] { return std::make_unique<SigmaModel>(); },
        &SigmaModel::schema);
}

} // namespace kenref
