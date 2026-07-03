/*
 * RelaxModel.cpp
 *
 *  RELAX (relaxation-rate) energy model wrapped as a registry EnergyModel — the payoff of the
 *  model-abstraction restructure: wiring a new model is ONE self-registering TU, with NO edits to any
 *  engine consumer (GROMACS / PLUMED / the offline tools all dispatch it through the registry). Concerns:
 *    (1) inputs   : relaxation spectral-density data (load_spec_den_relax_data from EXP_DATA_FOLDER),
 *                   relaxation rates (default, or from RATES_FILE).
 *    (2) params   : the Model-tier ParamSchema below.
 *    (3) caches   : relax_data_list_ + its primed atomId-pair caches; the name->id map and the compact
 *                   atomId pairs.
 *  compute() forwards to the unchanged, R-validated KEnRef::coord_array_to_relax_energy kernel
 *  (energy 1840.2342988885609 + gradient against the R reference; see testCoordArrayToRelaxEnergy).
 *
 *  Mirrors SigmaModel exactly, minus PROTON_MHZ: RELAX has no spectrometer field and its rates set
 *  carries the diffusion-tensor components (Dx/Dy/Dz) in addition to kens/kmethyl/karo.
 */

#include "core/EnergyModel.h"
#include "core/ModelRegistry.h"
#include "core/KEnRef.h"
#include "core/IoUtils.h"
#include "core/Table.h"

namespace kenref {

namespace {

class RelaxModel : public EnergyModel<KEnRef_Real_t> {
public:
    static ParamSchema schema() {
        ParamSchema s;
        s.add({"EXP_DATA_FOLDER", ParamType::Path, ParamTier::Model, /*required*/ true, std::nullopt,
               "Folder with relaxation spectral-density experimental data (RELAX model)"});
        s.add({"RATES_FILE",      ParamType::Path, ParamTier::Model, /*required*/ false, std::nullopt,
               "CSV file (a_coef.csv style: header row of rate names + one data row of values) "
               "overriding the default RELAX relaxation rates"});
        return s;
    }

    void buildCache(const InitContext<KEnRef_Real_t>& ctx) override {
        exp_data_folder_ = ctx.getString("EXP_DATA_FOLDER");
        relax_data_list_ = IoUtils::load_spec_den_relax_data(exp_data_folder_, ctx.handleNames);

        // Parameter-validation: num_models must match the ensemble-member count the exp-data was built
        // for. RELAX shares SIGMA's file-based grouping (grouping_matrix.cols() / numModels inside the hot
        // coord_array_to_relax_energy kernel), so a mismatch would likewise corrupt the heap. Check once
        // here at setup, not per-step. numModelsTotal is provided for this model-count-dependent setup.
        IoUtils::assertEnsembleMemberCount(relax_data_list_, ctx.numModelsTotal, "RELAX");

        // Optionally override the default rates_ from a CSV file (header = rate names, one data row of
        // values) read as a 1-row table -> NamedRowVector; otherwise keep the built-in default. The
        // RELAX default carries the diffusion-tensor components Dx/Dy/Dz alongside kens/kmethyl/karo.
        const std::string rates_file = ctx.getString("RATES_FILE");
        if (!rates_file.empty()) {
            rates_ = IoUtils::readTable(rates_file, true, false).toNamedRowVector<KEnRef_Real_t>();
        } else {
            rates_ = Table({{"5.0e+08", "1.0e+12", "1.0e+04", "2.5e+08", "2.5e+08", "2.5e+08"}},
                           {{"kens", "kmethyl", "karo", "Dx", "Dy", "Dz"}}).toNamedRowVector<KEnRef_Real_t>();
        }

        // Collect all atom-name pairs across every SpecDenRelaxData object (the driver builds
        // sub-indexing from their union).
        atomName_pairs_.clear();
        for (const auto& rd : relax_data_list_)
            for (const auto& ap : rd.get_atom_pairs())
                atomName_pairs_.emplace_back(ap);
    }

    [[nodiscard]] const std::vector<std::tuple<std::string, std::string>>& atomNamePairs() const override {
        return atomName_pairs_;
    }

    void finalizeIndexing(const IndexingContext<KEnRef_Real_t>& ctx) override {
        atomNames_2_atomIds_ = ctx.atomNameToSub0Id; // the kernel maps SpecDenRelaxData name pairs through this
        atomId_pairs_ = IoUtils::atomNamePairs_2_atomIdPairs(atomName_pairs_, atomNames_2_atomIds_, ctx.numOmpThreads);

        // Prime each SpecDenRelaxData's atomId-pair cache (the RELAX equivalent of the SIGMA priming
        // the consumers did by hand; otherwise the kernel re-resolves names via std::map every call).
        for (auto& rd : relax_data_list_)
            rd.set_atomIdPairs_to_sub0Atom_id_pairs_cache(
                IoUtils::atomNamePairs_2_atomIdPairs(rd.get_atom_pairs(), atomNames_2_atomIds_, ctx.numOmpThreads));
    }

    [[nodiscard]] const std::vector<std::tuple<int, int>>& atomIdPairs() const override {
        return atomId_pairs_;
    }

    std::tuple<KEnRef_Real_t, std::optional<std::vector<CoordsMatrixType<KEnRef_Real_t>>>>
    compute(const StepContext<KEnRef_Real_t>& ctx) override {
        return KEnRef<KEnRef_Real_t>::coord_array_to_relax_energy(
            ctx.coord_array, rates_, relax_data_list_,
            ctx.k, ctx.n, atomNames_2_atomIds_, ctx.gradient, ctx.numOmpThreads);
    }

    [[nodiscard]] KEnRef_Real_t forceUnitScale() const override { return static_cast<KEnRef_Real_t>(10.0); }

private:
    std::string   exp_data_folder_;
    NamedRowVector<KEnRef_Real_t> rates_;
    std::vector<SpecDenRelaxData<KEnRef_Real_t>> relax_data_list_;
    std::vector<std::tuple<std::string, std::string>> atomName_pairs_;
    std::vector<std::tuple<int, int>> atomId_pairs_;
    std::map<std::string, int> atomNames_2_atomIds_;
};

} // namespace

// Called by ModelBootstrap::bootstrapModels() under #if KENREF_ENABLE_RELAX.
void registerRelaxModel() {
    ModelRegistry<KEnRef_Real_t>::add(
        "RELAX",
        [] { return std::make_unique<RelaxModel>(); },
        &RelaxModel::schema);
}

} // namespace kenref
