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
    std::string ext = std::filesystem::path(inputFiles.front()).extension().string();
    inputFileType = ext == ".trr" ? InputFileType::trr : ext == ".xtc" ? InputFileType::xtc : InputFileType::UNKNOWN;
    using CLI::enums::operator<<;
    std::cout << "Energy model: " << selected_energy_model << "\n";

    auto atomName_to_atomGlobalId_map /*(1-based)*/ = IoUtils::getAtomMappingFromPdb<std::string, int>(refFileName, IoUtils::fill_atomId_to_index_Map);
    const bool handleNames = IoUtils::should_handleNames(atomName_to_atomGlobalId_map);

    if (selected_energy_model == KEnRef<double>::energyModel::SIGMA) {
        handle_sigma_energy_model(atomName_to_atomGlobalId_map, handleNames);
    } else if (selected_energy_model == KEnRef<KEnRef_Real_t>::energyModel::PLATEAUS) {
        handle_plateaus_energy_model();
    } else {
        GMX_ASSERT(selected_energy_model != KEnRef<KEnRef_Real_t>::energyModel::UNKNOWN, (toolParamPrefix + "_ENERGY_MODEL must be set to either 'SIGMA' or 'PLATEAUS'.").c_str());
    }

    std::cout << "Current path is " << std::filesystem::current_path() << std::endl;
    //Guide atom indices
    const std::vector<int> &guideAtom0Indices = GmxKEnRefInitializer::loadGmxIndexGroup(GUIDE_C_ALPHA, indexFileName);
    IoUtils::printVector(guideAtom0Indices);
    //Total number of atoms in the system
    size_t homenr = GmxKEnRefInitializer::loadGmxIndexGroup("System", indexFileName).size();
    assert(homenr > 0 && "No group named \"System\" found in index file.");

    //Guide atoms X
    //load all reference coordinates (including both guide atoms and reference atoms)
    auto referenceStructureAllAtomCoordsMap = IoUtils::getAtomMappingFromPdb<int, Eigen::RowVector3<KEnRef_Real_t> >(
        refFileName, IoUtils::fill_atomIndex1_to_coords_Map<KEnRef_Real_t>);

    //calculate guideAtomsReferenceCoordsCentered (to use for fitting)
    //const auto &guideAtomsReferenceCoordsCentered = Kabsch_Umeyama<KEnRef_Real_t>::translateCenterOfMassToOrigin(
    //        IoUtils::extractCoords(guideAtom0Indices, false, referenceStructureAllAtomCoordsMap, true));
    const auto &guideAtomsReferenceCoords = IoUtils::extractCoords(guideAtom0Indices, false, referenceStructureAllAtomCoordsMap, true);

    //N.B. globalAtomIdFlags is ZERO based, in contrast to its corresponding KEnRefForceProvider one
    std::vector globalAtomIdFlags(homenr, false);
    int maxAtomIdOfInterest = -1;
    {
        for (const auto &[a1, a2]: atomIdPairs) {
            //In the next lines I use .at() instead of [] deliberately; to throw an exception if unexpected name found
            if (a1 > maxAtomIdOfInterest) maxAtomIdOfInterest = a1;
            if (a2 > maxAtomIdOfInterest) maxAtomIdOfInterest = a2;
            globalAtomIdFlags.at(a1) = true;
            globalAtomIdFlags.at(a2) = true;
        }
        globalAtomIdFlags.resize(maxAtomIdOfInterest + 1);
    }
    //prepare sub0Id_to_global0Id
    auto global0Id_to_sub0Id = std::vector(globalAtomIdFlags.size(), -1);
    // auto sub0Id_to_global0Id = std::vector(globalAtomIdFlags.size(), -1);
    std::vector<int> subAtoms0Ids;
    {
        int sub0Id = 0;
        for (int i = 0; i < globalAtomIdFlags.size(); i++) {
            if (globalAtomIdFlags[i]) {
                global0Id_to_sub0Id[i] = sub0Id;
                // sub0Id_to_global0Id[sub0Id] = i;
                sub0Id++;
                subAtoms0Ids.emplace_back(i);
            }
        }
        // sub0Id_to_global0Id.resize(sub0Id);
    }
    // by the end of the above block, sub0Id_to_global0Id should not have any -1 items
    for (auto &[atom1, atom2]: atomIdPairs) {
        subAtomIdPairs.emplace_back(global0Id_to_sub0Id[atom1], global0Id_to_sub0Id[atom2]);
    }

    fill_atomName_to_atomSub0Id_map_if_needed(atomName_to_atomGlobalId_map, maxAtomIdOfInterest, global0Id_to_sub0Id);

    std::cout << metricToolCalculatesCapitalized+" output file path: " << outputPathName << std::endl;
    std::filesystem::path file_path(outputPathName);
    // Create the directories if they do not exist
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

    std::vector<t_file_state> fsts(num_models);
    Eigen::VectorX<int> returns(num_models);
    Eigen::VectorX<int> oks(num_models);


    //I need
    // Data Once: atomIdPairs, guideAtomsX_ZEROIndexed
    // Actions Once: translating guideAtomsX_ZEROIndexed to the point of origin
    // Data Every frame: coordsVector
    // Actions Every frame: fitting coordsVector to guideAtomsX_ZEROIndexed
    try {
        for (int modelId = 0; modelId < num_models; ++modelId) {
            auto &fst = fsts[modelId];
            const std::string &modelPathName = inputFiles[modelId];

            switch (inputFileType) {
                case InputFileType::xtc:
                    fst.xd = open_xtc(modelPathName, "r");
                    read_first_xtc(fst.xd, &fst.natoms, &fst.step, &fst.time, fst.box, &fst.x, &fst.prec, &fst.bOK);
                    break;
                case InputFileType::trr:
                    gmx_trr_header_t header;
                    gmx_trr_read_single_header(modelPathName, &header);
                    // fst.x = new rvec[header.x_size];
                    snew(fst.x, header.x_size);
                    fst.xd = gmx_trr_open(modelPathName, "r");
                    fst.bOK = gmx_trr_read_frame(fst.xd, &fst.step, &fst.time, &fst.lambda, fst.box, &fst.natoms, fst.x, nullptr, nullptr);
                    break;
                default:
                    std::cerr << "FATAL ERROR: Unrecognized input file type.\n";
            }
            fst.nframe = 0;
        }

        auto guideAtomsX_ZEROIndexed = CoordsMatrixType<KEnRef_Real_t>(guideAtom0Indices.size(), 3);
        allSimulationsSubAtomsXVector.resize(num_models);
        do {

            //Use the data from the 2 models
            for (int modelIdx = 0; modelIdx < num_models; ++modelIdx) {
                auto &fst = fsts[modelIdx];
                //calculate guideAtomsX_ZEROIndexed (every model every frame)
                fillX(guideAtomsX_ZEROIndexed, guideAtom0Indices, fst.x, true);

                //remember that the data must be PBC corrected (in every step)

                //calculate the transformation matrix
                //const auto &affineForProcessing = Kabsch_Umeyama<KEnRef_Real_t>::find3DAffineTransform(guideAtomsX_ZEROIndexed, guideAtomsReferenceCoordsCentered, true, true);
                //another way to achieve the above line
                const auto &affineForProcessing = Kabsch_Umeyama<KEnRef_Real_t>::find3DAffineTransform(
                    guideAtomsX_ZEROIndexed, guideAtomsReferenceCoords, false, true);

                //and calculate subAtomsXAfterTransform
                auto subAtomsX = CoordsMatrixType<KEnRef_Real_t>(subAtoms0Ids.size(), 3);
                fillX(subAtomsX, subAtoms0Ids, fst.x, true);
                allSimulationsSubAtomsXVector.at(modelIdx) = Kabsch_Umeyama<KEnRef_Real_t>::applyTransform(affineForProcessing, subAtomsX);

                if (modelIdx == num_models - 1) {
                    calculate_and_report(fst);
                }

                fst.nframe++;
                switch (inputFileType) {
                    case InputFileType::xtc:
                        returns[modelIdx] = read_next_xtc(fst.xd, fst.natoms, &fst.step, &fst.time, fst.box, fst.x, &fst.prec, &fst.bOK);
                        break;
                    case InputFileType::trr:
                        returns[modelIdx] = gmx_trr_read_frame(fst.xd, &fst.step, &fst.time, &fst.lambda, fst.box, &fst.natoms, fst.x, nullptr, nullptr);
                        fst.bOK = returns[modelIdx]; //this is not exactly correct, but it is fine
                        break;
                    default:
                        std::cerr << "FATAL ERROR: Unrecognized input file type.\n";
                }
                oks[modelIdx] = fst.bOK;
            }
        } while ((returns.array() != 0).all() && (max_frame < 0 || fsts[0].nframe <= max_frame));
        if ((inputFileType == InputFileType::xtc && !oks.all()) ||
            (inputFileType == InputFileType::trr && oks.rows() > 1 && !(oks.array() == 0).all())) {
            fprintf(stderr, "\nWARNING: Incomplete frame.\n");
        }
        for (auto &fst: fsts) {
            sfree(fst.x);
            switch (inputFileType) {
                case InputFileType::xtc:
                    close_xtc(fst.xd);
                    break;
                case InputFileType::trr:
                    gmx_trr_close(fst.xd);
                    break;
                default:
                    std::cerr << "FATAL ERROR: Unrecognized input file type.\n";
            }
        }
        out_file_stream.close();
    } catch (...) {
        for (auto &fst: fsts) {
            sfree(fst.x);
            switch (inputFileType) {
                case InputFileType::xtc:
                    close_xtc(fst.xd);
                    break;
                case InputFileType::trr:
                    gmx_trr_close(fst.xd);
                    break;
                default:
                    std::cerr << "FATAL ERROR: Unrecognized input file type.\n";
            }
        }
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
            ->required()->envname(toolParamPrefix + "_INPUT") ->check(CLI::ExistingFile);
    app.add_option("-d,--index", indexFileName, "Index file name")->envname(toolParamPrefix + "_INDEX");
    app.add_option("-r,--ref", refFileName, "Reference file")->envname(toolParamPrefix + "_REF");
    app.add_option("-o,--output", outputPathName, "output " + metricToolCalculates + " file")
            ->envname(toolParamPrefix + "_OUTPUT");
    app.add_option("-m,--model", selected_energy_model, "energy model")
            ->required()->envname(toolParamPrefix + "_ENERGY_MODEL")
            ->transform(CLI::CheckedTransformer(KEnRef<KEnRef_Real_t>::energyModelMap, CLI::ignore_case));
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

void AbstractCalculator::handle_sigma_energy_model(std::map<std::string, int> atomName_to_atomGlobalId_map, const bool handleNames) {
    const std::vector<std::string> &spec_den_data_prefixes = KEnRef<KEnRef_Real_t>::spec_den_data_prefixes;
    // spec_den_data_vector.reserve(spec_den_data_prefixes.size());
    for (const auto &spec_den_data_prefix : spec_den_data_prefixes) {
        //AtomPairs
        std::string atomPairAndSigmaFileName = experimentalDataFolder + spec_den_data_prefix + "_atom_pairs.csv";
        std::cout << atomPairAndSigmaFileName << std::endl;
        const Table &atomPairAndSigmaTable =
                IoUtils::readTable(atomPairAndSigmaFileName, true, false, "\\s*,\\s*", -1, true);
        std::vector<std::tuple<std::string, std::string> > iterationAtomPairs(atomPairAndSigmaTable.rowCount());
        for (int j = 0; j < atomPairAndSigmaTable.rowCount(); ++j) {
            auto atom1 = IoUtils::normalizeName(atomPairAndSigmaTable.at(j, "atom1"), handleNames);
            auto atom2 = IoUtils::normalizeName(atomPairAndSigmaTable.at(j, "atom2"), handleNames);
            iterationAtomPairs.at(j) = std::move(std::tuple{atom1, atom2});
            int i1 = atomName_to_atomGlobalId_map[atom1] - 1;
            int i2 = atomName_to_atomGlobalId_map[atom2] - 1;
            atomIdPairs.emplace_back(i1, i2);
        }
        fill_spec_den_data_vector(spec_den_data_prefix, atomPairAndSigmaTable, iterationAtomPairs);
    }
}