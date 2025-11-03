#include <iostream>
#include <unistd.h>
#include "CLI11/CLI11.hpp"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/fileio/trrio.h"
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
    real     lambda;
};

KEnRef_Real_t pearsonCorrelation(const Eigen::VectorX<KEnRef_Real_t>& x, const Eigen::VectorX<KEnRef_Real_t>& y) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("Vectors must be of the same length.");
    }

    const KEnRef_Real_t mean_x = x.mean();
    const KEnRef_Real_t mean_y = y.mean();

    Eigen::VectorX<KEnRef_Real_t> diff_x = x.array() - mean_x;
    Eigen::VectorX<KEnRef_Real_t> diff_y = y.array() - mean_y;

    const KEnRef_Real_t numerator = (diff_x.array() * diff_y.array()).sum();
    const KEnRef_Real_t denominator = std::sqrt((diff_x.array().square().sum()) * (diff_y.array().square().sum()));

    return numerator / denominator;
}

Eigen::IOFormat insideCsvLineFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "", "", "", "", "");

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

class S2OrderParamsCalculator {

    CoordsMatrixType<KEnRef_Real_t> lastFrameSubAtomsX_; //Used only for proper NoJump algorithm
    CoordsMatrixType<KEnRef_Real_t> lastFrameGuideAtomsX_ZEROIndexed_; //Used only for proper NoJump algorithm

    enum class InputFileType { xtc, trr, UNKNOWN };
public:

//remember that the data must be PBC corrected (in every step)
    int calc(int argc, char** argv) const {
        CLI::App app{"energy calculator"};
        argv = app.ensure_utf8(argv); // new argv memory is held by app
        bool debug = false;
        app.add_flag("--debug", debug, "enable debugging (holds for debugging)");

        std::string GUIDE_C_ALPHA = "guideC-alpha";
        app.add_option("-g,--guide", GUIDE_C_ALPHA, "name of guide group")->envname("s2_GUIDE");

        std::vector<std::string> inputFiles;
        app.add_option("-i,--input", inputFiles, "Input files")
                ->required()->envname("S2_INPUT") ->check(CLI::ExistingFile);

        std::string indexFileName = "index.ndx";
        app.add_option("-d,--index", indexFileName, "Index file name")->envname("S2_INDEX");

        std::string refFileName = "ref.pdb";
        app.add_option("-r,--ref", refFileName, "Reference file")->envname("S2_REF");

        std::string s2OutputPathName = "s2.out";
        app.add_option("-o,--output", s2OutputPathName, "output S2 file")
                ->envname("S2_OUTPUT");

        KEnRef<KEnRef_Real_t>::energyModel selected_energy_model;
        app.add_option("-m,--model", selected_energy_model, "energy model")
                ->required()->envname("S2_ENERGY_MODEL")
                ->transform(CLI::CheckedTransformer(KEnRef<KEnRef_Real_t>::energyModelMap, CLI::ignore_case));

        std::string experimentalDataFileName, experimentalDataFolder;
        app.add_option("-x,--exp-data-folder", experimentalDataFolder, "experimental data folder for sigma")
                ->envname("S2_EXP_DATA_FOLDER")->check(CLI::ExistingDirectory);
        app.add_option("-X,--exp-data-file", experimentalDataFileName, "experimental data file for plateaus")
                ->envname("S2_EXP_DATA_FILE")->check(CLI::ExistingFile);

        int max_frame = -1;
        app.add_option("--max-frame", max_frame, "maximum number of frames to read")
                ->envname("S2_MAX_FRAME");

        uint dt = 10;
        app.add_option("--dt", dt, "dt time step to report energy")->envname("S2_DT");

        // Load from config file
        app.set_config("--params", "params.toml", "Read a TOML config file", false)
                ->envname("S2_PARAMS");
        CLI11_PARSE(app, argc, argv);
        if (debug) {
            volatile bool holdToDebug = true;
            while (holdToDebug) {
                sleep(1);
            }
        }

        ///////////////////////////////////////////////////////

        int numModels = static_cast<int>(inputFiles.size());
        std::string ext = std::filesystem::path(inputFiles.front()).extension().string();
        InputFileType inputFileType = ext == ".trr" ? InputFileType::trr : ext == ".xtc" ? InputFileType::xtc : InputFileType::UNKNOWN;
        using CLI::enums::operator<<;
        std::cout << "Energy model: " << selected_energy_model << "\n";

        std::vector<std::vector<std::vector<int> > > simulated_grouping_list;
        std::vector<SpecDenData<KEnRef_Real_t> > spec_den_data_vector;
        std::vector<std::tuple<int, int> > atomIdPairs; // Zero based

        auto atomName_to_atomGlobalId_map /*(1-based)*/ = IoUtils::getAtomMappingFromPdb<std::string, int>(refFileName, IoUtils::fill_atomId_to_index_Map);
        const bool handleNames = IoUtils::should_handleNames(atomName_to_atomGlobalId_map);

        if (selected_energy_model == KEnRef<double>::energyModel::SIGMA) {
            const std::vector<std::string> &spec_den_data_prefixes = KEnRef<KEnRef_Real_t>::spec_den_data_prefixes;
            spec_den_data_vector.reserve(spec_den_data_prefixes.size());
            for (const auto &spec_den_data_prefix: spec_den_data_prefixes) {
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
            }
        } else if (selected_energy_model == KEnRef<KEnRef_Real_t>::energyModel::PLATEAUS) {
            auto experimentalData_table = IoUtils::readTable(experimentalDataFileName, true, false, ",");
            for (int i = 0; i < experimentalData_table.rowCount(); ++i) {
                std::istringstream temp1(experimentalData_table(i, "i1")), temp2(experimentalData_table(i, "i2")), temp3(experimentalData_table(i, "g1")), temp4(experimentalData_table(i, "g2"));
                int i1, i2;
                temp1 >> i1;
                temp2 >> i2;
                atomIdPairs.emplace_back(i1 - 1, i2 - 1);
            }
        } else {
            GMX_ASSERT(selected_energy_model != KEnRef<KEnRef_Real_t>::energyModel::UNKNOWN, "ENER_MODEL must be set to either 'SIGMA' or 'PLATEAUS'.");
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
        // auto sub0Id_to_global0Id = std::vector(globalAtomIdFlags.size(), -1); // NOT NEEDED?
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
        std::vector<std::tuple<int, int> > subAtomIdPairs;
        for (auto &[atom1, atom2]: atomIdPairs) {
            subAtomIdPairs.emplace_back(global0Id_to_sub0Id[atom1], global0Id_to_sub0Id[atom2]);
        }


        std::cout << "S2 output file path: " << s2OutputPathName << std::endl;
        std::filesystem::path file_path(s2OutputPathName);
        // Create the directories if they do not exist
        if (auto parent = file_path.parent_path(); !parent.empty() && !exists(parent)) {
            std::error_code ec;
            std::filesystem::create_directories(parent, ec);
            if (ec) {
                std::cerr << "Error creating directories: " << ec.message() << std::endl;
            }
        }
        std::ofstream s2OutFileStream;
        s2OutFileStream.open(s2OutputPathName);
        if (s2OutFileStream.is_open()) {
            std::cout << "S2 Output file open successfully\n";
        } else {
            std::cerr << "FATAL ERROR: Can't open file [" << s2OutputPathName << "] for writing.\n";
            exit(-1);
        }

        std::vector<t_file_state> fsts(numModels);
        Eigen::VectorX<int> returns(numModels);
        Eigen::VectorX<int> oks(numModels);


        //        //Calculate reference atoms (to use to calculate reference S2 order params later)
        //        // To do so, we need 1) subAtomsIds (0?) and 2) global coords from the reference PDB file.
        //        const auto& referenceStructureSubAtomCoords = extractCoords(subAtoms0Ids, false, referenceStructureAllAtomCoordsMap, true);
        //          No need. We already have S2 from the experimental table csv file

        //I need
        // Data Once: atomIdPairs, guideAtomsX_ZEROIndexed
        // Actions Once: translating guideAtomsX_ZEROIndexed to the point of origin
        // Data Every frame: coordsVector
        // Actions Every frame: fitting coordsVector to guideAtomsX_ZEROIndexed
        try {
            for (int modelId = 0; modelId < numModels; ++modelId) {
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
            do {
                std::vector<CoordsMatrixType<KEnRef_Real_t> > allSimulationsSubAtomsXVector(numModels);

                //Use the data from the 2 models
                for (int modelIdx = 0; modelIdx < numModels; ++modelIdx) {
                    auto &fst = fsts[modelIdx];
                    //calculate guideAtomsX_ZEROIndexed (every model every frame)
                    fillX(guideAtomsX_ZEROIndexed, guideAtom0Indices, fst.x, true);

                    //remember that the data must be PBC corrected (in every step)

                    //calculate the transformation matrix
                    //const auto &affineForS2 = Kabsch_Umeyama<KEnRef_Real_t>::find3DAffineTransform(guideAtomsX_ZEROIndexed, guideAtomsReferenceCoordsCentered, true, true);
                    //another way to achieve the above line
                    const auto &affineForS2 = Kabsch_Umeyama<KEnRef_Real_t>::find3DAffineTransform(
                        guideAtomsX_ZEROIndexed, guideAtomsReferenceCoords, false, true);

                    //and calculate subAtomsXAfterTransform
                    auto subAtomsX = CoordsMatrixType<KEnRef_Real_t>(subAtoms0Ids.size(), 3);
                    fillX(subAtomsX, subAtoms0Ids, fst.x, true);
                    allSimulationsSubAtomsXVector.at(modelIdx) = Kabsch_Umeyama<KEnRef_Real_t>::applyTransform(affineForS2, subAtomsX);

                    if (modelIdx == numModels - 1) {
                        if (fst.step % dt == 0) {
                            const auto &frameS2OrderParams = KEnRef<KEnRef_Real_t>::s2OrderParams( allSimulationsSubAtomsXVector, subAtomIdPairs, 0);
                            if (fst.nframe == 0) {
                                std::cout << "First S2OrderParams of first step\n" << frameS2OrderParams.topRows(25).transpose() << std::endl;
                            }
                            s2OutFileStream << fst.step << ", " << frameS2OrderParams.transpose().format(insideCsvLineFormat) << std::endl;
                        }
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
            s2OutFileStream.close();
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
            s2OutFileStream.flush();
            s2OutFileStream.close();
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
    S2OrderParamsCalculator().calc(argc, argv);
}

