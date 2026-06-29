#include "tools/AbstractCalculator.h"
#include "gmxinterface/gmxkenrefinitializer.h"

int AbstractCalculator::calc(int argc, char **argv) {
    CLI::App app{toolName};
    argv = app.ensure_utf8(argv); // new argv memory is held by app

    add_common_parameters(app);
    addSpecificParameters(app);

    CLI11_PARSE(app, argc, argv);
    if (debug) {
        volatile bool holdToDebug = true;
        while (holdToDebug) {
            sleep(1);
        }
    }
    ///////////////////////////////////////////////////////
    num_models = static_cast<int>(inputFiles.size());
    modelName = to_upper(modelName); // registry keys are upper-case (e.g. SIGMA, PLATEAUS)
    std::cout << "Number of Models = " << num_models << std::endl;
    std::cout << "Energy model: " << modelName << "\n";

    // ---- 1. reference PDB atom-name -> 1-based global id ----
    auto atomName_to_atomGlobalId_map /*(1-based)*/ = IoUtils::getAtomMappingFromPdb<std::string, int>(
        refFileName, IoUtils::fill_atomId_to_index_Map);
    GMX_ASSERT(!atomName_to_atomGlobalId_map.empty(), "No atom mapping found in reference file");
    const bool handleNames = IoUtils::should_handleNames(atomName_to_atomGlobalId_map);

    // ---- 2. guide indices + System size from the index file ----
    std::cout << "Current path is " << std::filesystem::current_path() << std::endl;
    const std::vector<int> guideAtom0Indices = GmxKEnRefInitializer::loadGmxIndexGroup(GUIDE_C_ALPHA, indexFileName);
    IoUtils::printVector(guideAtom0Indices);
    const size_t homenr = GmxKEnRefInitializer::loadGmxIndexGroup("System", indexFileName).size();
    assert(homenr > 0 && "No group named \"System\" found in index file.");

    // ---- 3. param map + the offline EngineAdapter ----
    std::map<std::string, std::string> params;
    if (!experimentalDataFolder.empty())   params["EXP_DATA_FOLDER"] = experimentalDataFolder;
    if (!experimentalDataFileName.empty())  params["EXP_DATA_FILE"]   = experimentalDataFileName;
    addModelParams(params); // tool-specific (energycalc: PROTON_MHZ / RATES_FILE)
    adapter = std::make_unique<TrajectoryEngineAdapter>(inputFiles, guideAtom0Indices, params, kernelNumOmpThreads());

    // ---- 4. build the selected model from the registry and let it load its inputs ----
    kenref::bootstrapModels();
    auto model = kenref::ModelRegistry<KEnRef_Real_t>::create(modelName);
    GMX_ASSERT(model != nullptr,
               (toolParamPrefix + "_ENERGY_MODEL must name a registered model (e.g. SIGMA or PLATEAUS).").c_str());
    {
        kenref::ParamSchema schema; // empty: the adapter supplies every param the model reads
        kenref::InitContext<KEnRef_Real_t> initCtx{*adapter, schema, atomName_to_atomGlobalId_map, handleNames,
                                                   num_models};
        model->buildCache(initCtx);
    }

    // ---- 5. derive the sub-atom indexing from the model's atom-name pairs (mirrors the GROMACS force
    //         provider's fillParamsStep0). Everything here is 0-based to index the trajectory arrays. ----
    const auto &atomName_pairs = model->atomNamePairs();
    int maxAtomId0OfInterest = -1;
    std::vector<bool> globalAtom0IdFlags(homenr, false);
    int tmp;
    for (const auto &[a1, a2]: atomName_pairs) {
        //In the next lines I use .at() instead of [] deliberately; to throw an exception if unexpected name found
        if ((tmp = atomName_to_atomGlobalId_map.at(a1) - 1) > maxAtomId0OfInterest) maxAtomId0OfInterest = tmp;
        globalAtom0IdFlags[tmp] = true;
        if ((tmp = atomName_to_atomGlobalId_map.at(a2) - 1) > maxAtomId0OfInterest) maxAtomId0OfInterest = tmp;
        globalAtom0IdFlags[tmp] = true;
    }
    globalAtom0IdFlags.resize(maxAtomId0OfInterest + 1);

    std::vector<int> global0Id_to_sub0Id(globalAtom0IdFlags.size(), -1);
    std::vector<int> subAtoms0Ids;
    {
        int sub0Id = 0;
        for (int i = 0; i < static_cast<int>(globalAtom0IdFlags.size()); ++i) {
            if (globalAtom0IdFlags[i]) {
                global0Id_to_sub0Id[i] = sub0Id++;
                subAtoms0Ids.emplace_back(i);
            }
        }
    }
    adapter->setSubAtoms0Ids(subAtoms0Ids);

    std::map<std::string, int> atomName_to_atomSub0Id_map;
    for (const auto &[name, global1Id]: atomName_to_atomGlobalId_map) {
        if (global1Id - 1 > maxAtomId0OfInterest)
            continue;
        const int sub = global0Id_to_sub0Id.at(global1Id - 1);
        if (sub >= 0)
            atomName_to_atomSub0Id_map[name] = sub;
    }

    // ---- 6. finalize the model's indexing and capture the compact sub-atom pairs ----
    {
        kenref::IndexingContext<KEnRef_Real_t> idxCtx{atomName_to_atomSub0Id_map, kernelNumOmpThreads()};
        model->finalizeIndexing(idxCtx);
    }
    subAtomIdPairs = model->atomIdPairs();

    // ---- 7. centered guide-atom reference coords (Angstrom) for Kabsch fitting ----
    auto referenceStructureAllAtomCoordsMap = IoUtils::getAtomMappingFromPdb<int, Eigen::RowVector3<KEnRef_Real_t> >(
        refFileName, IoUtils::fill_atomIndex1_to_coords_Map<KEnRef_Real_t>);
    const auto guideAtomsReferenceCoords = IoUtils::extractCoords(guideAtom0Indices, false,
                                                                  referenceStructureAllAtomCoordsMap, true);
    guideAtomsReferenceCoordsCentered = Kabsch_Umeyama<KEnRef_Real_t>::translateCenterOfMassToOrigin(guideAtomsReferenceCoords);

    // ---- 8. let the concrete tool take the prepared model (build a driver, etc.) ----
    prepareConsumer(std::move(model));

    // ---- 9. output file ----
    std::cout << metricToolCalculatesCapitalized+" output file path: " << outputPathName << std::endl;
    std::filesystem::path file_path(outputPathName);
    if (auto parent = file_path.parent_path(); !parent.empty() && !exists(parent)) {
        std::error_code ec;
        std::filesystem::create_directories(parent, ec);
        if (ec) {
            std::cerr << "Error creating directories: " << ec.message() << std::endl;
        }
    }
    out_file_stream.open(outputPathName);
    if (out_file_stream.is_open()) {
        std::cout << metricToolCalculatesCapitalized+" Output file open successfully\n";
    } else {
        std::cerr << "FATAL ERROR: Can't open file [" << outputPathName << "] for writing.\n";
        exit(-1);
    }

    // ---- 10. step every trajectory in lock-step; the concrete tool processes the current frame ----
    try {
        adapter->openAll();
        do {
            perFrameReport();
        } while (adapter->advanceAll() && (max_frame < 0 || adapter->currentFrame() <= max_frame));
        adapter->closeAll();
        out_file_stream.close();
    } catch (...) {
        adapter->closeAll();
        out_file_stream.flush();
        out_file_stream.close();
        std::rethrow_exception(std::current_exception());
    }
    return 0;
}

void AbstractCalculator::add_common_parameters(CLI::App &app) {
    app.add_flag("--debug", debug, "enable debugging (holds for debugging)");
    app.add_option("-g,--guide", GUIDE_C_ALPHA, "name of guide group")->envname(toolParamPrefix + "_GUIDE");
    app.add_option("-i,--input", inputFiles, "Input files")
            ->required()->envname(toolParamPrefix + "_INPUT")->check(CLI::ExistingFile);
    app.add_option("-d,--index", indexFileName, "Index file name")->envname(toolParamPrefix + "_INDEX");
    app.add_option("-r,--ref", refFileName, "Reference file")->envname(toolParamPrefix + "_REF");
    app.add_option("-o,--output", outputPathName, "output " + metricToolCalculates + " file")
            ->envname(toolParamPrefix + "_OUTPUT");
    app.add_option("-m,--model", modelName, "energy model name (e.g. SIGMA, PLATEAUS)")
            ->required()->envname(toolParamPrefix + "_ENERGY_MODEL");
    app.add_option("-x,--exp-data-folder", experimentalDataFolder, "experimental data folder for sigma")
            ->envname(toolParamPrefix + "_EXP_DATA_FOLDER")->check(CLI::ExistingDirectory);
    app.add_option("-X,--exp-data-file", experimentalDataFileName, "experimental data file for plateaus")
            ->envname(toolParamPrefix + "_EXP_DATA_FILE")->check(CLI::ExistingFile);
    app.add_option("--max-frame", max_frame, "maximum number of frames to read")
            ->envname(toolParamPrefix + "_MAX_FRAME");
    app.add_option("--dt", dt, "dt time step to report energy")->envname(toolParamPrefix + "_DT");
    // Load from config file
    app.set_config("--params", "params.toml", "Read a TOML config file", false)
            ->envname(toolParamPrefix + "_PARAMS");
}
