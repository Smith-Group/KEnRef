#include <filesystem>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <regex>
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

Eigen::IOFormat insideCsvLineFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "", "", "", "", "");

void fillX(CoordsMatrixType<KEnRef_Real_t> &targetAtomsX, const std::vector<int> &idxes0,
           const rvec *x, bool toAngstrom) {
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

    int numModels=0;
    CoordsMatrixType<KEnRef_Real_t> lastFrameSubAtomsX_; //Used only for proper NoJump algorithm
    CoordsMatrixType<KEnRef_Real_t> lastFrameGuideAtomsX_ZEROIndexed_; //Used only for proper NoJump algorithm
public:

//remember that the data must be PBC corrected (in every step)
    void calc() {
        const std::string GUIDE_C_ALPHA = IoUtils::getEnvParam("ENER_GUIDE_C_ALPHA" ,"guideC-alpha");
        const std::string INDEX_FILE_LOCATION = IoUtils::getEnvParam("ENER_INDEX_FILE_LOCATION", "../res/10nsstart+fitting/KEnRefAtomIndex.ndx");
        const std::string REFERENCE_FILENAME = IoUtils::getEnvParam("ENER_REFERENCE_FILENAME", "../res/10nsstart+fitting/00ns.pdb");
//        const std::string inPathTemplate = IoUtils::getEnvParam("ENER_IN_PATH", "/smithlab/home/aalhossary/temp/res/out_restrained/out_multi_data${dataDir}_structure${structure}/${repl}");
        const std::string inPathTemplate = IoUtils::getEnvParam("ENER_IN_PATH", "/smithlab/home/aalhossary/temp/temp_xtc/data${dataDir}_structure${structure}");
        const std::string outPathTemplate = IoUtils::getEnvParam("ENER_OUT_PATH", "/smithlab/home/aalhossary/temp/tempout");
        const std::string experimentalDataFileName = IoUtils::getEnvParam("ENER_experimentalDataFileName", "../res/10nsstart+fitting/singleton_data_10nsstart+fit_3-5_1977pairs_80_A.csv");
//        const std::string inputFileNameTemplate = IoUtils::getEnvParam("ENER_fileNameTemplate", "KEnRef_multi_${K}_${N}_${max}_${variant}_${trial}.xtc");
        const std::string inputFileNameTemplate = IoUtils::getEnvParam("ENER_fileNameTemplate", "${repl}_KEnRef_multi_${K}_${N}_${max}_${variant}_${trial}.xtc");
        const std::string enerOutFileTemplate = IoUtils::getEnvParam("ENER_OUTPUT_FILE", "KEnRef_multi_${K}_${N}_${max}_${variant}_${trial}-1.out");
        const std::string dataDir = IoUtils::getEnvParam("ENER_DATA_DIR", "alef+alef");
        const std::string structure = IoUtils::getEnvParam("ENER_STRUCTURE", "alef+alef");
        const std::string k = IoUtils::getEnvParam("ENER_K", "1e9");
        const std::string max = IoUtils::getEnvParam("ENER_MAX", "800");
        const std::string subSet = IoUtils::getEnvParam("ENER_SUBSET_", "80");
        const std::string n = IoUtils::getEnvParam("ENER_N", ".25");
        const std::string variant = IoUtils::getEnvParam("ENER_VARIANT", "A");
        const int trial = IoUtils::getEnvParam("ENER_TRIAL", 1);
        const int MAX_FRAME = IoUtils::getEnvParam("ENER_MAX_FRAME", 5000);

        const std::string replicatesIn = IoUtils::getEnvParam("ENER_REPLICATES", "repl_01,repl_02");
        std::vector<std::string> replicates = IoUtils::split(replicatesIn, ",");
        numModels = static_cast<int>(replicates.size());

        std::vector<std::vector<std::vector<int>>>simulated_grouping_list;
        switch (numModels) {
            case 1:
                simulated_grouping_list = std::vector<std::vector<std::vector<int>>>{{{0}}, {{0}}};
                break;
            case 2:
                simulated_grouping_list = std::vector<std::vector<std::vector<int>>>{{{0, 1}}, {{0}, {1}}};
                break;
            case 3:
                simulated_grouping_list = std::vector<std::vector<std::vector<int>>>{{{0, 1, 2}}, {{0}, {1}, {2}}};
                break;
            default:
                assert(numModels <= 3 && "I don't know how to handle more than 3 simulations yet");
        }

        std::cout << "Current path is "<< std::filesystem::current_path()<< std::endl;
        //Guide atom indices
        const std::vector<int> &guideAtom0Indices = GmxKEnRefInitializer::loadGmxIndexGroup(GUIDE_C_ALPHA, INDEX_FILE_LOCATION);
        IoUtils::printVector(guideAtom0Indices);
        //Total number of atoms in the system
        long homenr = GmxKEnRefInitializer::loadGmxIndexGroup("System", INDEX_FILE_LOCATION).size();
        assert(homenr > 0 && "No group named \"System\" found in index file.");

        //Guide atoms X
        //load all reference coordinates (including both guide atoms and reference atoms)
        auto referenceStructureAllAtomCoordsMap = IoUtils::getAtomMappingFromPdb<int, Eigen::RowVector3<KEnRef_Real_t>>(
                REFERENCE_FILENAME, IoUtils::fill_atomIndex1_to_coords_Map < KEnRef_Real_t > );

        //calculate guideAtomsReferenceCoordsCentered (to use for fitting)
        //const auto &guideAtomsReferenceCoordsCentered = Kabsch_Umeyama<KEnRef_Real_t>::translateCenterOfMassToOrigin(
        //        IoUtils::extractCoords(guideAtom0Indices, false, referenceStructureAllAtomCoordsMap, true));
        const auto &guideAtomsReferenceCoords = IoUtils::extractCoords(guideAtom0Indices, false, referenceStructureAllAtomCoordsMap, true);
        auto experimentalData_table = IoUtils::readTable(experimentalDataFileName, true, ",");
        std::vector<std::vector<std::string> > experimentalData_tableData = std::get<1>(experimentalData_table);

        std::vector<std::tuple<int, int>> atomIdPairs;
        Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, Eigen::Dynamic> g0(experimentalData_tableData.size(), 2);

        for (int i = 0; i < experimentalData_tableData.size(); ++i) {
            const auto & record = experimentalData_tableData[i];
            std::istringstream temp1(record[3]), temp2(record[4]), temp3(record[5]), temp4(record[6]);
            int i1, i2;
            temp1 >> i1;
            temp2 >> i2;
            atomIdPairs.emplace_back(i1 - 1, i2 - 1);
            temp3 >> g0(i, 0);
            temp4 >> g0(i, 1);
        }
        //N.B. globalAtomIdFlags is ZERO based, in contrast to its corresponding KEnRefForceProvider one
        std::vector<bool> globalAtomIdFlags(homenr, false);
        {
            int maxAtomIdOfInterest = -1;
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
        auto global0Id_to_sub0Id = std::vector<int>(globalAtomIdFlags.size(), -1);
        auto sub0Id_to_global0Id = std::vector<int>(globalAtomIdFlags.size(), -1);
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

        std::unordered_map<std::string, std::string> replacements{
                {"dataDir",   dataDir},
                {"structure", structure},
                {"K",         k},
                {"N",         n},
                {"max",       max},
                {"subSet",    subSet},
                {"variant",   variant},
                {"trial",     std::to_string(trial)}
        };
        std::string enerOutputPathName = constructFileNamePath(outPathTemplate, enerOutFileTemplate, replacements);
        std::cout << "Energy output file path: " << enerOutputPathName << std::endl;
        std::filesystem::path file_path(enerOutputPathName);
        // Create the directories if they do not exist
        if (!exists(file_path.parent_path())){
            std::error_code ec;
            std::filesystem::create_directories(file_path.parent_path(), ec);
            if (ec) {
                std::cerr << "Error creating directories: " << ec.message() << std::endl;
            }
        }
        std::ofstream enerOutFileStream;
        enerOutFileStream.open(enerOutputPathName);
        if(enerOutFileStream.is_open()){
            std::cout << "Energy Output file open successfully\n";
        }else{
            std::cerr << "FATAL ERROR: Can't open file [" << enerOutputPathName << "] for writing.\n";
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
                replacements.insert_or_assign("repl", replicates[modelId]);
                std::string modelPathName = constructFileNamePath(inPathTemplate, inputFileNameTemplate, replacements);

                fst.xd = open_xtc(modelPathName, "r");
                read_first_xtc(fst.xd, &fst.natoms, &fst.step, &fst.time, fst.box, &fst.x, &fst.prec, &fst.bOK);
                fst.nframe = 0;
            }

            CoordsMatrixType<KEnRef_Real_t> guideAtomsX_ZEROIndexed = CoordsMatrixType<KEnRef_Real_t>(guideAtom0Indices.size(), 3);
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
                    CoordsMatrixType<KEnRef_Real_t> subAtomsX = CoordsMatrixType<KEnRef_Real_t>(subAtoms0Ids.size(), 3);
                    fillX(subAtomsX, subAtoms0Ids, fst.x, true);
                    allSimulationsSubAtomsXVector.at(modelIdx) = Kabsch_Umeyama<KEnRef_Real_t>::applyTransform(affineForEnergy, subAtomsX);

                    if (modelIdx == numModels - 1){
                        KEnRef_Real_t energy = 0;
                        std::vector<CoordsMatrixType<KEnRef_Real_t> > allDerivatives_vector;
                        if constexpr(std::is_same_v<KEnRef_Real_t, float>) {
                            std::tie(energy, allDerivatives_vector) =
                                    KEnRef<KEnRef_Real_t>::coord_array_to_energy(allSimulationsSubAtomsXVector,
                                                                                 subAtomIdPairs,
                                                                                 simulated_grouping_list, g0,
                                                                                 std::stof(k), std::stof(n),
                                                                                 false, 1, (fst.step % 10 == 0));
                        } else {
                            std::tie(energy, allDerivatives_vector) =
                                    KEnRef<KEnRef_Real_t>::coord_array_to_energy(allSimulationsSubAtomsXVector,
                                                                                 subAtomIdPairs,
                                                                                 simulated_grouping_list, g0,
                                                                                 std::stod(k), std::stod(n),
                                                                                 false, 1, (fst.step % 10 == 0));
                        }
                        if (fst.step % 10 == 0) {
//                            std::cout << "Step: " << fst.step << " Energy: " << energy << std::endl;
                            enerOutFileStream << "Step: " << fst.step << " Energy: " << energy << std::endl;
                        }
                    }

                    fst.nframe++;
                    returns[modelIdx] = read_next_xtc(fst.xd, fst.natoms, &fst.step, &fst.time, fst.box, fst.x, &fst.prec, &fst.bOK);
                    oks[modelIdx] = fst.bOK;
                }
            } while ((returns.array() != 0).all() && fsts[0].nframe <= MAX_FRAME);
            if (! oks.all()) {
                fprintf(stderr, "\nWARNING: Incomplete frame.\n");
            }
            for (auto & fst:fsts) {
                sfree(fst.x);
                close_xtc(fst.xd);
            }
            enerOutFileStream.close();
        } catch (...) {
            for (auto & fst:fsts) {
                sfree(fst.x);
                close_xtc(fst.xd);
            }
            enerOutFileStream.flush();
            enerOutFileStream.close();
            std::rethrow_exception(std::current_exception());
        }
    }

    static std::string constructFileNamePath(const std::string &pathTemplate, const std::string &fileNameTemplate,
                                             const std::unordered_map<std::string, std::string> &replacements) {
        std::string currentModelPathName = pathTemplate;
        if (currentModelPathName.back() != '/')
            currentModelPathName.append("/");
        currentModelPathName.append(fileNameTemplate);
//        std::cout << currentModelPathName << std::endl;

        for (const auto &pair: replacements) {
            std::regex placeholder("\\$\\{"+ pair.first + "\\}");
            currentModelPathName = std::regex_replace(currentModelPathName, placeholder, pair.second);
        }
        return currentModelPathName;
    }
};


int main(int argc, char** argv){
    MPI_Init(&argc, &argv);
    if (argc >= 2 && std::strcmp(argv[1], "debug") == 0) {
        volatile bool holdToDebug = true;
        while (holdToDebug) {
            sleep(1);
        }
    }
    EnergyCalculator().calc();
}

