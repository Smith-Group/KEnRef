#include <filesystem>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <regex>
#include "CLI11/CLI11.hpp"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/gmxlib/network.h"
#include "core/kabsch.h"
#include "core/KEnRef.h"
#include "core/IoUtils.h"
#include "gmxinterface/gmxkenrefinitializer.h"


struct t_file_state{
    t_fileio *xd;
    rvec     *x;
    matrix    box;
    int       nframe, natoms;
    int64_t   step;
    real      prec, time;
    gmx_bool  bOK;
};

void fillX(CoordsMatrixType<KEnRef_Real_t> &targetAtomsX, const std::vector<int> &idxes0, const rvec *x, const bool toAngstrom) {
    for (int i = 0; i < targetAtomsX.rows(); i++) {
//        const int *piGlobal = new int{sub0Id_to_global1Id[i] - 1};
//        const int *piLocal = forceProviderInput.cr_.dd->ga2la->findHome(*piGlobal);
//        GMX_ASSERT(piLocal, "ERROR: Can't find local index of atom");
        const gmx::RVec atom_x = x[idxes0[i]];
#if VERBOSE
        std::cout << sub0Id_to_global1Id[i] << "\t" << *piGlobal << "\t" << *piLocal << "\t x: " << atom_x[0] << ", " << atom_x[1] << ", " << atom_x[2] << std::endl;
#endif
        if constexpr (std::is_same_v<KEnRef_Real_t, real>) {
            auto subAtomsX_buffer = targetAtomsX.data();
            const auto rvec = atom_x.as_vec();
            std::copy_n(rvec, 3, &subAtomsX_buffer[i * 3]);
        } else {
            for (int j = 0; j < 3; ++j) {
                targetAtomsX(i, j) = static_cast<KEnRef_Real_t>(atom_x[j]);
            }
        }
    }
    if (toAngstrom)
        targetAtomsX *= 10;
}

class EnergyCalculator {

    CoordsMatrixType<KEnRef_Real_t> lastFrameSubAtomsX_; //Used only for proper NoJump algorithm
    CoordsMatrixType<KEnRef_Real_t> lastFrameGuideAtomsX_ZEROIndexed_; //Used only for proper NoJump algorithm
    const NamedRowVector<KEnRef_Real_t> rates = Table(		//TODO provide a way to change it
        {{"5.0e+08", "2.5e+08", "1.0e+12", "1.0e+04"}},
        {{"kens", "kc", "kmethyl", "karo"}}
        ).toNamedRowVector<KEnRef_Real_t>();
public:

//remember that the data must be PBC corrected (in every step)
    int calc(int argc, char** argv) const {
        CLI::App app{"energy calculator"};
        argv = app.ensure_utf8(argv);  // new argv memory is held by app
        bool debug = false;
        app.add_flag("--debug", debug, "enable debugging (holds for debugging)");

        std::string GUIDE_C_ALPHA = "guideC-alpha";
        app.add_option("-g,--guide", GUIDE_C_ALPHA, "name of guide group")->envname("ENERGY_GUIDE");

        std::vector<std::string> inputFiles;
        app.add_option("-i,--input", inputFiles, "Input files")
            ->required()->envname("ENERGY_INPUT") ->check(CLI::ExistingFile);

        std::string indexFileName = "index.ndx";
        app.add_option("-d,--index", indexFileName, "Index file name")->envname("ENERGY_INDEX");

        std::string refFileName = "ref.pdb";
        app.add_option("-r,--ref", refFileName, "Reference file")->envname("ENERGY_REF");

        KEnRef_Real_t k, n;
        app.add_option("-k,--k", k, "K force constant")->envname("ENERGY_K");
        app.add_option("-n,--n", n, "power scaling")->envname("ENERGY_N");

        std::string energyOutputFileName = "energy.out";
        app.add_option("-o,--output", energyOutputFileName, "output energy file")
            ->envname("ENERGY_OUTPUT");

        KEnRef<KEnRef_Real_t>::energyModel selected_energy_model;
        app.add_option("-m,--model", selected_energy_model, "energy model")
            ->required()->envname("ENERGY_MODEL")
            ->transform(CLI::CheckedTransformer(KEnRef<KEnRef_Real_t>::energyModelMap, CLI::ignore_case));

        std::string experimentalDataFileName, experimentalDataFolder;
        app.add_option("-x,--exp-data-folder", experimentalDataFolder, "experimental data folder for sigma")
            ->envname("ENERGY_EXP_DATA_FOLDER")->check(CLI::ExistingDirectory);
        app.add_option("-X,--exp-data-file", experimentalDataFileName, "experimental data file for plateaus")
            ->envname("ENERGY_EXP_DATA_FILE")->check(CLI::ExistingFile);

        KEnRef_Real_t proton_mhz = 700.0;
        app.add_option("-z,--proton-mhz", proton_mhz, "spectrometer proton field strength in MHz")
            ->envname("ENERGY_PROTON_MHZ");

        int max_frame = -1;
        app.add_option("--max-frame", max_frame, "maximum number of frames to read")
            ->envname("ENERGY_MAX_FRAME");

        uint dt = 10;
        app.add_option("--dt", dt, "dt time step to report energy")->envname("ENERGY_DT");

        // Load from config file
        app.set_config("--params", "params.toml", "Read a TOML config file", false)
            ->envname("ENERGY_PARAMS");
        CLI11_PARSE(app, argc, argv);
        if (debug) {
            volatile bool holdToDebug = true;
            while (holdToDebug) {
                sleep(1);
            }
        }

        ///////////////////////////////////////////////////////

        int numModels = static_cast<int>(inputFiles.size());
        using CLI::enums::operator<<;
        std::cout << "Energy model: " << selected_energy_model << "\n";

        std::vector<std::vector<std::vector<int>>>simulated_grouping_list;
        std::vector<SpecDenData<KEnRef_Real_t> > spec_den_data_vector;
        std::vector<std::tuple<int, int>> atomIdPairs; // Zero based
        Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, Eigen::Dynamic> g0;

        auto atomName_to_atomGlobalId_map /*(1-based)*/= IoUtils::getAtomMappingFromPdb<std::string, int>(refFileName,IoUtils::fill_atomId_to_index_Map);
        const bool handleNames = IoUtils::should_handleNames(atomName_to_atomGlobalId_map);

        if (selected_energy_model == KEnRef<double>::energyModel::SIGMA) {
            const std::vector<std::string> &spec_den_data_prefixes = KEnRef<KEnRef_Real_t>::spec_den_data_prefixes;
            spec_den_data_vector.reserve(spec_den_data_prefixes.size());
            for (const auto & spec_den_data_prefix : spec_den_data_prefixes) {
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
                    int i1 = atomName_to_atomGlobalId_map[atom1] -1;
                    int i2 = atomName_to_atomGlobalId_map[atom2] -1;
                    atomIdPairs.emplace_back(i1, i2);
                }
                // sigma
                std::vector<KEnRef_Real_t> sigmasVec{};
                for (int row = 0; row < atomPairAndSigmaTable.rowCount(); ++row) {
                    if (!atomPairAndSigmaTable.isRowComplete(row))
                        break;
                    const auto &valueStr = atomPairAndSigmaTable.at(row, "sigma");
                    std::istringstream iss(valueStr);
                    KEnRef_Real_t value;
                    iss >> value;
                    sigmasVec.emplace_back(value);
                }
                std::optional<NamedVector<KEnRef_Real_t> > sigma(sigmasVec.size());
                for (int j = 0; j < sigma->rows(); ++j) {
                    sigma.value()(j, 0) = sigmasVec[j];
                }
                //multiple_grouping
                std::string multiple_grouping_fileName = experimentalDataFolder + spec_den_data_prefix + "_groupings.csv";
                const auto &grouping_matrix = NamedMatrix<int>(IoUtils::readTable(
                    multiple_grouping_fileName, false, false).toNamedMatrix<int>().array() - 1);
                const auto &multiple_grouping = IoUtils::grouping_mat_to_subset_idx(grouping_matrix);
                //a_coef
                std::string aCoefFileName = experimentalDataFolder + spec_den_data_prefix + "_a_coef.csv";
                const auto &a_coef = IoUtils::readTable(aCoefFileName, true, false, "\\s*,\\s*", -1, false).
                        toNamedMatrix<KEnRef_Real_t>();
                //lambda_coef
                std::string lambdaCoefFileName =
                        experimentalDataFolder + spec_den_data_prefix + "_lambda_coef.csv";
                const auto &lambda_coef = IoUtils::readTable(lambdaCoefFileName, true, true, "\\s*,\\s*", -1, false).
                        toNamedMatrix<KEnRef_Real_t>();

                spec_den_data_vector.emplace_back(std::move(SpecDenData{
                    iterationAtomPairs, sigma, multiple_grouping, a_coef, lambda_coef }));
            }
        }else if (selected_energy_model == KEnRef<KEnRef_Real_t>::energyModel::PLATEAUS) {
            switch (numModels) {
                case 1:
                    simulated_grouping_list = std::vector<std::vector<std::vector<int> > >{{{0}}, {{0}}};
                    break;
                case 2:
                    simulated_grouping_list = std::vector<std::vector<std::vector<int> > >{{{0, 1}}, {{0}, {1}}};
                    break;
                case 3:
                    simulated_grouping_list = std::vector<std::vector<std::vector<int> > > {{{0, 1, 2}}, {{0}, {1}, {2}}};
                    break;
                default:
                    assert(numModels <= 3 && "I don't know how to handle more than 3 simulations yet");
            }
            auto experimentalData_table = IoUtils::readTable(experimentalDataFileName, true, false, ",");
            g0.resize(static_cast<int>(experimentalData_table.rowCount()), 2);
            for (int i = 0; i < experimentalData_table.rowCount(); ++i) {
                std::istringstream temp1(experimentalData_table(i, "i1")), temp2(experimentalData_table(i, "i2")), temp3(experimentalData_table(i, "g1")), temp4(experimentalData_table(i, "g2"));
                int i1, i2;
                temp1 >> i1;
                temp2 >> i2;
                atomIdPairs.emplace_back(i1 - 1, i2 - 1);
                temp3 >> g0(i, 0);
                temp4 >> g0(i, 1);
            }
        }else {
            GMX_ASSERT(selected_energy_model != KEnRef<KEnRef_Real_t>::energyModel::UNKNOWN, "ENER_MODEL must be set to either 'SIGMA' or 'PLATEAUS'.");
        }

        std::cout << "Current path is "<< std::filesystem::current_path()<< std::endl;
        //Guide atom indices
        const std::vector<int> &guideAtom0Indices = GmxKEnRefInitializer::loadGmxIndexGroup(GUIDE_C_ALPHA, indexFileName);
        IoUtils::printVector(guideAtom0Indices);
        //Total number of atoms in the system
        size_t homenr = GmxKEnRefInitializer::loadGmxIndexGroup("System", indexFileName).size();
        assert(homenr > 0 && "No group named \"System\" found in index file.");

        //Guide atoms X
        //load all reference coordinates (including both guide atoms and reference atoms)
        auto referenceStructureAllAtomCoordsMap = IoUtils::getAtomMappingFromPdb<int, Eigen::RowVector3<KEnRef_Real_t>>(
                refFileName, IoUtils::fill_atomIndex1_to_coords_Map < KEnRef_Real_t > );

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
        auto sub0Id_to_global0Id = std::vector(globalAtomIdFlags.size(), -1);
        std::vector<int> subAtoms0Ids;
        {
            int sub0Id = 0;
            for (int i = 0; i < globalAtomIdFlags.size(); i++) {
                if (globalAtomIdFlags[i]) {
                    global0Id_to_sub0Id[i] = sub0Id;
                    sub0Id_to_global0Id[sub0Id] = i;
                    sub0Id++;
                    subAtoms0Ids.emplace_back(i);
                }
            }
            sub0Id_to_global0Id.resize(sub0Id);
        }
        // by the end of the above block, sub0Id_to_global0Id should not have any -1 items
        std::vector<std::tuple<int, int>> subAtomIdPairs;
        for (auto &[atom1, atom2]:atomIdPairs) {
            subAtomIdPairs.emplace_back(global0Id_to_sub0Id[atom1], global0Id_to_sub0Id[atom2]);
        }

        std::map<std::string, int> atomName_to_atomSub0Id_map;
        for (const auto &[name, global1Id]: atomName_to_atomGlobalId_map){
            if (global1Id > maxAtomIdOfInterest + 1)
                continue;
            atomName_to_atomSub0Id_map[name] = global0Id_to_sub0Id.at(global1Id - 1);
        }

        std::cout << "Energy output file path: " << energyOutputFileName << std::endl;
        std::filesystem::path file_path(energyOutputFileName);
        // Create the directories if they do not exist
        if (auto parent = file_path.parent_path(); !parent.empty() && !exists(parent)){
            std::error_code ec;
            std::filesystem::create_directories(parent, ec);
            if (ec) {
                std::cerr << "Error creating directories: " << ec.message() << std::endl;
            }
        }
        std::ofstream energyOutFileStream;
        energyOutFileStream.open(energyOutputFileName);
        if(energyOutFileStream.is_open()){
            std::cout << "Energy Output file open successfully\n";
        }else{
            std::cerr << "FATAL ERROR: Can't open file [" << energyOutputFileName << "] for writing.\n";
            exit(-1);
        }

        std::vector<t_file_state> fsts(numModels);
        Eigen::VectorX<int> returns(numModels);
        Eigen::VectorX<int> oks(numModels);


        //I need
        // Data Once: atomIdPairs, guideAtomsX_ZEROIndexed
        // Actions Once: translating guideAtomsX_ZEROIndexed to the point of origin
        // Data Every frame: coordsVector
        // Actions Every frame: fitting coordsVector to guideAtomsX_ZEROIndexed
        try {
            for (int modelId = 0; modelId < numModels; ++modelId) {
                auto& fst = fsts[modelId];
                const std::string& modelPathName = inputFiles[modelId];

                fst.xd = open_xtc(modelPathName, "r");
                read_first_xtc(fst.xd, &fst.natoms, &fst.step, &fst.time, fst.box, &fst.x, &fst.prec, &fst.bOK);
                fst.nframe = 0;
            }

            auto guideAtomsX_ZEROIndexed = CoordsMatrixType<KEnRef_Real_t>(guideAtom0Indices.size(), 3);
            do {
                std::vector<CoordsMatrixType<KEnRef_Real_t>> allSimulationsSubAtomsXVector(numModels);

                //Use the data from the 2 models
                for (int modelIdx = 0; modelIdx < numModels; ++modelIdx) {
                    auto &fst = fsts[modelIdx];
                    //calculate guideAtomsX_ZEROIndexed (every model every frame)
                    fillX(guideAtomsX_ZEROIndexed, guideAtom0Indices, fst.x, true);

                    //remember that the data must be PBC corrected (in every step)

                    //calculate the transformation matrix
                    //const auto &affineForEnergy = Kabsch_Umeyama<KEnRef_Real_t>::find3DAffineTransform(guideAtomsX_ZEROIndexed, guideAtomsReferenceCoordsCentered, true, true);
                    //another way to achieve the above line
                    const auto &affineForEnergy = Kabsch_Umeyama<KEnRef_Real_t>::find3DAffineTransform(
                            guideAtomsX_ZEROIndexed, guideAtomsReferenceCoords, false, true);

                    //and calculate subAtomsXAfterTransform
                    auto subAtomsX = CoordsMatrixType<KEnRef_Real_t>(subAtoms0Ids.size(), 3);
                    fillX(subAtomsX, subAtoms0Ids, fst.x, true);
                    allSimulationsSubAtomsXVector.at(modelIdx) = Kabsch_Umeyama<KEnRef_Real_t>::applyTransform(affineForEnergy, subAtomsX);

                    if (modelIdx == numModels - 1){
                        KEnRef_Real_t energy = 0;
                        std::optional<std::vector<CoordsMatrixType<KEnRef_Real_t> > > allDerivatives_vector;
                        if (selected_energy_model == KEnRef<double>::energyModel::PLATEAUS) {
                            std::tie(energy, allDerivatives_vector) =
                                KEnRef<KEnRef_Real_t>::coord_array_to_g_energy(allSimulationsSubAtomsXVector, subAtomIdPairs,
                                    simulated_grouping_list, g0, k, n, false, 1);
                        }else if (selected_energy_model == KEnRef<double>::energyModel::SIGMA) {
                            std::tie(energy, allDerivatives_vector) =
                                    KEnRef<KEnRef_Real_t>::coord_array_to_sigma_energy(
                                        allSimulationsSubAtomsXVector, rates, spec_den_data_vector, proton_mhz, k, n,
                                        atomName_to_atomSub0Id_map, false, 1);
                        }
                        if (fst.step % dt == 0) {
                            // std::cout << "Step: " << fst.step << " Energy: " << energy << std::endl;
                            // energyOutFileStream << "Step: " << fst.step << " Energy: " << energy << std::endl;
                            energyOutFileStream << fst.step << '\t' << energy << std::endl;
                        }
                    }

                    fst.nframe++;
                    returns[modelIdx] = read_next_xtc(fst.xd, fst.natoms, &fst.step, &fst.time, fst.box, fst.x, &fst.prec, &fst.bOK);
                    oks[modelIdx] = fst.bOK;
                }
            } while ((returns.array() != 0).all() && (max_frame < 0 || fsts[0].nframe <= max_frame));
            if (! oks.all()) {
                fprintf(stderr, "\nWARNING: Incomplete frame.\n");
            }
            for (auto & fst:fsts) {
                sfree(fst.x);
                close_xtc(fst.xd);
            }
            energyOutFileStream.close();
        } catch (...) {
            for (auto & fst:fsts) {
                sfree(fst.x);
                close_xtc(fst.xd);
            }
            energyOutFileStream.flush();
            energyOutFileStream.close();
            std::rethrow_exception(std::current_exception());
        }
        return 0;
    }

    static std::string constructFileNamePath(const std::string &pathTemplate, const std::string &fileNameTemplate,
        const std::unordered_map<std::string, std::string> &replacements) {

        std::string currentModelPathName = pathTemplate;
        if (currentModelPathName.back() != '/')
            currentModelPathName.append("/");
        currentModelPathName.append(fileNameTemplate);
//        std::cout << currentModelPathName << std::endl;

        for (const auto &[fst, snd]: replacements) {
            std::regex placeholder("\\$\\{"+ fst + "\\}");
            currentModelPathName = std::regex_replace(currentModelPathName, placeholder, snd);
        }
        return currentModelPathName;
    }
};

int main(int argc, char** argv){
    MPI_Init(&argc, &argv);
    EnergyCalculator().calc(argc, argv);
}

