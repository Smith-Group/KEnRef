#include "tools/AbstractCalculator.h"
#include "gmxinterface/gmxkenrefinitializer.h"
#include "core/buildModelIndexing.h"

#include <filesystem>   // std::filesystem::{current_path,path,create_directories} — libc++ pulls this in
                        // transitively, libstdc++ does NOT; include it explicitly so both toolchains build.

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

    // ---- 4-6. shared model + sub-atom indexing (registry create + buildCache + sub0Id indexing +
    //          finalizeIndexing), via kenref::buildModelIndexing — the SAME path the GROMACS and PLUMED
    //          engines use. We pass the 1-based PDB serials and convert the returned 1-based
    //          subAtomGlobalIds to the 0-based indices the trajectory arrays use. ----
    GMX_ASSERT(kenref::ModelRegistry<KEnRef_Real_t>::has(modelName),
               (toolParamPrefix + "_ENERGY_MODEL must name a registered model (e.g. SIGMA or PLATEAUS).").c_str());
    auto mi = kenref::buildModelIndexing<KEnRef_Real_t>(
        modelName, *adapter, atomName_to_atomGlobalId_map, handleNames, num_models, kernelNumOmpThreads());

    std::vector<int> subAtoms0Ids(mi.subAtomGlobalIds.size());
    for (size_t i = 0; i < subAtoms0Ids.size(); ++i)
        subAtoms0Ids[i] = mi.subAtomGlobalIds[i] - 1; // 1-based PDB serial -> 0-based trajectory index
    adapter->setSubAtoms0Ids(std::move(subAtoms0Ids));
    subAtomIdPairs = mi.model->atomIdPairs();

    // ---- 7. centered guide-atom reference coords (Angstrom) for Kabsch fitting ----
    auto referenceStructureAllAtomCoordsMap = IoUtils::getAtomMappingFromPdb<int, Eigen::RowVector3<KEnRef_Real_t> >(
        refFileName, IoUtils::fill_atomIndex1_to_coords_Map<KEnRef_Real_t>);
    const auto guideAtomsReferenceCoords = IoUtils::extractCoords(guideAtom0Indices, false,
                                                                  referenceStructureAllAtomCoordsMap, true);
    guideAtomsReferenceCoordsCentered = Kabsch_Umeyama<KEnRef_Real_t>::translateCenterOfMassToOrigin(guideAtomsReferenceCoords);

    // ---- 8. let the concrete tool take the prepared model (build a driver, etc.) ----
    prepareConsumer(std::move(mi.model));

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
