#include <iostream>
#include <iomanip>
#include <unistd.h>
#include "gromacs/fileio/xtcio.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
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

KEnRef_Real_t pearsonCorrelation(const Eigen::VectorX<KEnRef_Real_t>& x, const Eigen::VectorX<KEnRef_Real_t>& y) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("Vectors must be of the same length.");
    }

    KEnRef_Real_t mean_x = x.mean();
    KEnRef_Real_t mean_y = y.mean();

    Eigen::VectorX<KEnRef_Real_t> diff_x = x.array() - mean_x;
    Eigen::VectorX<KEnRef_Real_t> diff_y = y.array() - mean_y;

    KEnRef_Real_t numerator = (diff_x.array() * diff_y.array()).sum();
    KEnRef_Real_t denominator = std::sqrt((diff_x.array().square().sum()) * (diff_y.array().square().sum()));

    return numerator / denominator;
}

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

template<typename KEnRef_Real>
class S2OrderParamsCalculator {

    int numModels;

public:

    //remember that the data must be PBC corrected (in every step)
    void calc() {
        const std::string GUIDE_C_ALPHA = IoUtils::getEnvParam("S2_GUIDE_C_ALPHA" ,"guideC-alpha");
        const std::string INDEX_FILE_LOCATION = IoUtils::getEnvParam("S2_INDEX_FILE_LOCATION", "../res/10nsstart+fitting/KEnRefAtomIndex.ndx");
        const std::string REFERENCE_FILENAME = IoUtils::getEnvParam("S2_REFERENCE_FILENAME", "../res/10nsstart+fitting/00ns.pdb");
        const std::string inPathTemplate = IoUtils::getEnvParam("S2_IN_PATH", "/smithlab/home/aalhossary/temp/temp_xtc/data${dataDir}_structure${structure}");
        const std::string outPathTemplate = IoUtils::getEnvParam("S2_OUT_PATH", "/smithlab/home/aalhossary/temp/s2-order-params");
        const std::string experimentalDataFileName = IoUtils::getEnvParam("S2_experimentalDataFileName", "../res/10nsstart+fitting/singleton_data_10nsstart+fit_alef+baaa_3-5_1814pairs_80_A.csv");
        const std::string fileNameTemplate = IoUtils::getEnvParam("S2_fileNameTemplate", "repl_${repl}_KEnRef_multi_${K}_.25_${max}_${variant}_${trial}.xtc");
        const std::string s2OutFileTemplate = IoUtils::getEnvParam("S2_OUTPUT_FILE", "S2_data${dataDir}_structure${structure}_KEnRef_multi_${K}_.25_${max}_${variant}_${trial}.csv");
//        const std::string s2rOutFileTemplate = IoUtils::getEnvParam("S2R_OUTPUT_FILE", "S2r_data${dataDir}_structure${structure}_KEnRef_multi_${K}_.25_${max}_${variant}_${trial}.csv");

        const std::string dataDir = IoUtils::getEnvParam("S2_DATA_DIR", "alef+baaa");
        const std::string structure = IoUtils::getEnvParam("S2_STRUCTURE", "alef+baaa");
        const std::string k = IoUtils::getEnvParam("S2_K", "1e9");
        const std::string max = IoUtils::getEnvParam("S2_MAX", "800");
//        const std::string subSet = IoUtils::getEnvParam("S2_SUBSET", "80");
        const int subSet = IoUtils::getEnvParam("S2_SUBSET_", 80);
        const std::string n = IoUtils::getEnvParam("S2_N", ".25");
        const std::string variant = IoUtils::getEnvParam("S2_VARIANT", "A");
        const int trial = IoUtils::getEnvParam("S2_TRIAL", 2);
        const int MAX_FRAME = IoUtils::getEnvParam("S2_MAX_FRAME", 1000);

        const std::string replicatesIn = IoUtils::getEnvParam("S2_REPLICATES", "repl_01,repl_02");
        std::vector<std::string> replicates = IoUtils::split(replicatesIn, ",");
        numModels = (int) replicates.size();

        //Guide atom indices
        const std::vector<int> &guideAtom0Indices = GmxKEnRefInitializer::loadGmxIndexGroup(GUIDE_C_ALPHA, INDEX_FILE_LOCATION);
        IoUtils::printVector(guideAtom0Indices);
        //Total number of atoms in the system
        long homenr = GmxKEnRefInitializer::loadGmxIndexGroup("System", INDEX_FILE_LOCATION).size();
        assert(homenr > 0);

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
        Eigen::VectorX<KEnRef_Real_t> referenceS2OrderParams(experimentalData_tableData.size());

        for (int i = 0; i < experimentalData_tableData.size(); ++i) {
            const auto & record = experimentalData_tableData[i];
            std::istringstream temp1(record[3]), temp2(record[4]), temp3(record[7]);
            int i1, i2; KEnRef_Real_t f;
            temp1 >> i1;
            temp2 >> i2;
            temp3 >> f;
            atomIdPairs.emplace_back(i1 - 1, i2 - 1);
            referenceS2OrderParams(i) = f;
        }
#if VERBOSE
//        std::cout << "referenceS2OrderParams\n" << referenceS2OrderParams.transpose() << std::endl;
#endif
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

        std::unordered_map<std::string, std::string> replacements{};
//                {"dataDir",   dataDir},
//                {"structure", structure},
//                {"K",         k},
//                {"N",         n},
//                {"max",       max},
//                {"subSet",    subSet}, //TODO FIXME Avoiding a strange bug
//                {"variant",   variant},
//                {"trial",     std::to_string(trial)}
//        };
        replacements["dataDir"] = dataDir;
        replacements["subSet"] = std::to_string(subSet);
        replacements["structure"] = structure;
        replacements["K"] = k;
        replacements["N"] = n;
        replacements["variant"] = variant;
        replacements["max"] = max;
        replacements["trial"] = std::to_string(trial);
        std::string s2OutputPathName = constructFileNamePath(outPathTemplate, s2OutFileTemplate, replacements);
        std::cout << "S2 output file path: " << s2OutputPathName << std::endl;
        std::ofstream s2OutFileStream;
        s2OutFileStream.open(s2OutputPathName);
        if(s2OutFileStream.is_open()){
            std::cout << "S2 Output file open successfully\n";
        }else{
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
                auto& fst = fsts[modelId];
//                replacements.insert_or_assign("repl", IoUtils::padWithZeros(modelId + 1, 2));
                replacements.insert_or_assign("repl", replicates[modelId]);
                std::string modelPathName = constructFileNamePath(inPathTemplate, fileNameTemplate, replacements);

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
                    //const auto &affineForS2 = Kabsch_Umeyama<KEnRef_Real_t>::find3DAffineTransform(guideAtomsX_ZEROIndexed, guideAtomsReferenceCoordsCentered, true, true);
                    //another way to achieve the above line
                    const auto &affineForS2 = Kabsch_Umeyama<KEnRef_Real_t>::find3DAffineTransform(
                            guideAtomsX_ZEROIndexed, guideAtomsReferenceCoords, false, true);

                    //and calculate subAtomsXAfterTransform
                    CoordsMatrixType<KEnRef_Real_t> subAtomsX = CoordsMatrixType<KEnRef_Real_t>(subAtoms0Ids.size(), 3);
                    fillX(subAtomsX, subAtoms0Ids, fst.x, true);
                    allSimulationsSubAtomsXVector.at(modelIdx) = Kabsch_Umeyama<KEnRef_Real_t>::applyTransform(affineForS2, subAtomsX);

                    if (modelIdx == numModels - 1){
                        const auto &frameS2OrderParams = KEnRef<KEnRef_Real_t>::s2OrderParams(allSimulationsSubAtomsXVector, subAtomIdPairs, 0);
                        if (fst.nframe == 0){
                            std::cout << "referenceS2OrderParams\n" << referenceS2OrderParams.topRows(25).transpose() << "\n";
                            std::cout << "frameS2OrderParams\n" << frameS2OrderParams.topRows(25).transpose() << "\n";
                        }
                        s2OutFileStream << fst.nframe << ", " << frameS2OrderParams.transpose().format(insideCsvLineFormat) << std::endl;
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
            s2OutFileStream.close();
        } catch (const std::runtime_error& e) {
            for (auto & fst:fsts) {
                sfree(fst.x);
                close_xtc(fst.xd);
            }
            s2OutFileStream.close();
            throw e;
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
            std::string str = "\\$\\{"+ pair.first + "\\}";
//            std::string str = ""; //"\\$\\{";
//            str.append(pair.first);
            std::regex placeholder(str);
            currentModelPathName = std::regex_replace(currentModelPathName, placeholder, pair.second);
        }
        return currentModelPathName;
    }
};

//void testRead() {
//    const char* fn = "/smithlab/home/aalhossary/temp/res/out_restrained/out_multi_dataalef+alef_structurealef+alef/repl_01/KEnRef_multi_1e8_.25_800_A_1.xtc";
//    t_fileio* xd;
//    int       indent;
//    char      buf[256];
//    rvec*     x;
//    matrix    box;
//    int       nframe, natoms;
//    int64_t   step;
//    real      prec, time;
//    gmx_bool  bOK;
//
//    xd = open_xtc(fn, "r");
//    read_first_xtc(xd, &natoms, &step, &time, box, &x, &prec, &bOK);
//
//    nframe = 0;
//    do {
//        sprintf(buf, "%s frame %d", fn, nframe);
//        indent = 0;
////        indent = pr_title(stdout, indent, buf);
////        pr_indent(stdout, indent);
//        fprintf(stdout, "natoms=%10d  step=%10" PRId64 "  time=%12.7e  prec=%10g\n", natoms, step, time, prec);
//        pr_rvecs(stdout, indent, "box", box, DIM);
//        pr_rvecs(stdout, indent, "x", x, natoms);
//        nframe++;
//    } while (read_next_xtc(xd, natoms, &step, &time, box, x, &prec, &bOK) != 0);
//    if (!bOK) {
//        fprintf(stderr, "\nWARNING: Incomplete frame at time %g\n", time);
//    }
//    sfree(x);
//    close_xtc(xd);
//}

int main(int argc, char** argv){
    MPI_Init(&argc, &argv);
    if (argc >= 2 && std::strcmp(argv[1], "debug") == 0) {
        volatile bool holdToDebug = true;
        while (holdToDebug) {
            sleep(1);
        }
    }
//    testRead();
    S2OrderParamsCalculator<KEnRef_Real_t>().calc();
}

