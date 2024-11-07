#include <iostream>
#include <iomanip>
#include "mpi.h"
#include "core/KEnRef.h"
#include "core/IoUtils.h"
#include "core/kabsch.h"
#include "gmxinterface/gmxkenrefinitializer.h"

const int numModels = 2;


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


template<typename KEnRef_Real>
class S2OrderParamsCalculator {

public:

    //remember that the data must be PBC corrected (in every step)
    static void calc() {
        const std::string GUIDE_C_ALPHA = IoUtils::getEnvParam("S2_GUIDE_C_ALPHA" ,"guideC-alpha");
        const std::string INDEX_FILE_LOCATION = IoUtils::getEnvParam("S2_INDEX_FILE_LOCATION", "../res/10nsstart+fitting/KEnRefAtomIndex.ndx");
        const std::string REFERENCE_FILENAME = IoUtils::getEnvParam("S2_REFERENCE_FILENAME", "../res/10nsstart+fitting/00ns.pdb");
        const std::string inPath = IoUtils::getEnvParam("S2_IN_PATH", "/smithlab/home/aalhossary/temp/test");
        const std::string outPath = IoUtils::getEnvParam("S2_OUT_PATH", "/smithlab/home/aalhossary/temp/s2-order-params");
        //const std::string experimentalDataFileName = IoUtils::getEnvParam("S2_experimentalDataFileName", "../res/google_tests/singleton_data_10nsstart+fit_0+10.csv");
        //const std::string experimentalDataFileName = IoUtils::getEnvParam("S2_experimentalDataFileName", "../res/10nsstart+fitting/singleton_data_10nsstart+fit_3-5_1977pairs_all.csv");
        //const std::string experimentalDataFileName = IoUtils::getEnvParam("S2_experimentalDataFileName", "../res/10nsstart+fitting/singleton_data_10nsstart+fit_0+10_3-5_1812pairs_all.csv");
        const std::string experimentalDataFileName = IoUtils::getEnvParam("S2_experimentalDataFileName", "../res/10nsstart+fitting/singleton_data_10nsstart+fit_0+10_3-5_1812pairs_80_A.csv");
        const std::string fileNameTemplate = IoUtils::getEnvParam("S2_fileNameTemplate", "data${dataDir}_structure${structure}_repl_${repl}_KEnRef_multi_${K}_.25_${max}_${variant}_${trial}.${frame}.pdb");
        const std::string s2OutFileTemplate = IoUtils::getEnvParam("S2_OUTPUT_FILE", "S2_data${dataDir}_structure${structure}_KEnRef_multi_${K}_.25_${max}_${variant}_${trial}.csv");
        const std::string s2rOutFileTemplate = IoUtils::getEnvParam("S2R_OUTPUT_FILE", "S2r_data${dataDir}_structure${structure}_KEnRef_multi_${K}_.25_${max}_${variant}_${trial}.csv");

        std::string dataDir = IoUtils::getEnvParam("S2_DATA_DIR", "0+10");
        std::string structure = IoUtils::getEnvParam("S2_STRUCTURE", "0+10");
        std::string k = IoUtils::getEnvParam("S2_K", "1e9");
        std::string max = IoUtils::getEnvParam("S2_MAX", "999");
        std::string subSet = IoUtils::getEnvParam("S2_SUBSET", "80");
        std::string variant = IoUtils::getEnvParam("S2_VARIANT", "A");
        int trial = IoUtils::getEnvParam("S2_TRIAL", 2);
        const int MAX_FRAME = IoUtils::getEnvParam("S2_MAX_FRAME", 1000);

        //Guide atom indices
        const std::vector<int> &guideAtom0Indices = GmxKEnRefInitializer::loadGmxIndexGroup(GUIDE_C_ALPHA, INDEX_FILE_LOCATION);
        IoUtils::printVector(guideAtom0Indices);

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
        std::vector<bool> globalAtomIdFlags(experimentalData_tableData.size(), false);
        {
            int maxAtomIdOfInterest = -1;
            for (const auto &[a1, a2]: atomIdPairs) {
                //In the next lines I use .at() instead of [] deliberately; to throw an exception if unexpected name found
                if (a1 > maxAtomIdOfInterest) maxAtomIdOfInterest = a1;
                if (a2 > maxAtomIdOfInterest) maxAtomIdOfInterest = a2;
                globalAtomIdFlags[a1] = true;
                globalAtomIdFlags[a2] = true;
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

        std::map<std::string, std::string> replacements{
                {"dataDir",   dataDir},
                {"structure", structure},
                {"K",         k},
                {"max",       max},
                {"subSet",    subSet},
                {"variant",   variant},
                {"trial",     std::to_string(trial)},
        };
        std::string s2OutputPathName = constructFileNamePath(outPath, s2OutFileTemplate, replacements);
        std::string s2rOutputPathName = constructFileNamePath(outPath, s2rOutFileTemplate, replacements);
        std::cout << "S2 output file path: " << s2OutputPathName << std::endl;
        std::cout << "S2r output file path: " << s2rOutputPathName << std::endl;
        std::ofstream s2OutFileStream, s2rOutFileStream;
        s2OutFileStream.open(s2OutputPathName);
        s2rOutFileStream.open(s2rOutputPathName);
        if(s2OutFileStream.is_open()){
            std::cout << "S2 Output file open successfully\n";
        }else{
            std::cerr << "FATAL ERROR: Can't open file [" << s2OutputPathName << "] for writing.\n";
            exit(-1);
        }
        if(s2rOutFileStream.is_open()){
            std::cout << "S2r Output file open successfully\n";
        }else{
            std::cerr << "FATAL ERROR: Can't open file [" << s2rOutputPathName << "] for writing.\n";
            exit(-1);
        }


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
//            KEnRef_Real_t correlationCoefficient;
            for (int frame = 0; frame <= MAX_FRAME; frame++) {
                CoordsMatrixType<KEnRef_Real_t> guideAtomsX_ZEROIndexed;
                std::vector<CoordsMatrixType<KEnRef_Real_t>> allSimulationsSubAtomsXVector(numModels);

                for (int modelIdx = 0; modelIdx < numModels; modelIdx++) {

                    std::map<std::string, std::string> replacements_ext = {
                            {"frame",     std::to_string(frame)},
                            {"repl",      IoUtils::padWithZeros(modelIdx + 1, 2)},
                    };
                    for (auto const& node:replacements_ext) {
                        replacements.insert_or_assign(node.first, node.second);
                    }

                    std::string currentModelPathName = constructFileNamePath(inPath, fileNameTemplate, replacements);

                    //calculate guideAtomsX_ZEROIndexed (every model every frame)
                    auto modelAllAtomCoordsMap = IoUtils::getAtomMappingFromPdb<int, Eigen::RowVector3<KEnRef_Real_t>>(
                            currentModelPathName, IoUtils::fill_atomIndex1_to_coords_Map < KEnRef_Real_t > );
                    guideAtomsX_ZEROIndexed = IoUtils::extractCoords(guideAtom0Indices, false, modelAllAtomCoordsMap, true);

                    //remember that the data must be PBC corrected (in every step)

                    //calculate the transformation matrix
                    //const auto &affineForS2 = Kabsch_Umeyama<KEnRef_Real_t>::find3DAffineTransform(guideAtomsX_ZEROIndexed, guideAtomsReferenceCoordsCentered, true, true);
                    //another way to achieve the above line
                    const auto &affineForS2 = Kabsch_Umeyama<KEnRef_Real_t>::find3DAffineTransform(guideAtomsX_ZEROIndexed, guideAtomsReferenceCoords, false, true);

                    //and calculate subAtomsXAfterTransform
                    CoordsMatrixType<KEnRef_Real_t> subAtomsXAfterTransform;
                    {
                        CoordsMatrixType<KEnRef_Real_t> subAtomsX = IoUtils::extractCoords(subAtoms0Ids, false,
                                                                                           modelAllAtomCoordsMap, true);
                        if (modelIdx == 0) {
                            Eigen::RowVector3<KEnRef_Real_t> center_of_mass = guideAtomsX_ZEROIndexed.colwise().mean();
                            subAtomsXAfterTransform = subAtomsX.rowwise() - center_of_mass;
                        } else
                            subAtomsXAfterTransform = Kabsch_Umeyama<KEnRef_Real_t>::applyTransform(affineForS2, subAtomsX);
                    }

                    allSimulationsSubAtomsXVector.at(modelIdx) = subAtomsXAfterTransform;
                } // end of for(model)

                const auto &frameS2OrderParams = KEnRef<KEnRef_Real_t>::s2OrderParams(allSimulationsSubAtomsXVector, subAtomIdPairs, 0);
//                correlationCoefficient = pearsonCorrelation(frameS2OrderParams, referenceS2OrderParams);
#if false
                std::cout << frameS2OrderParams.transpose() << std::endl;
#endif
                if (frame == 0){
                    std::cout << "referenceS2OrderParams\n" << referenceS2OrderParams.topRows(25).transpose() << "\n";
                    std::cout << "frameS2OrderParams\n" << frameS2OrderParams.topRows(25).transpose() << "\n";
                }
//                std::cout << "pearsonCorrelation = " << correlationCoefficient << std::endl;
//            if (frame == 0){
//                std::cout << "referenceS2OrderParams = " << referenceS2OrderParams.transpose() << std::endl;
//                std::cout << "difference = " << (frameS2OrderParams - referenceS2OrderParams).transpose() << std::endl;
//            }
                s2OutFileStream << frame << ", " << frameS2OrderParams.transpose().format(insideCsvLineFormat) << "\n";
//                s2rOutFileStream << frame << ", " << correlationCoefficient << "\n";
            }
            s2OutFileStream.close();
            s2rOutFileStream.close();
//            return correlationCoefficient;
        } catch (const std::runtime_error& e) {
            s2OutFileStream.close();
            s2rOutFileStream.close();
            throw e;
        }
    }

    static std::string constructFileNamePath(const std::string &path, const std::string &fileNameTemplate,
                                             const std::map<std::string, std::string> &replacements) {
        std::string inputFileName = fileNameTemplate;
        for (const auto &pair: replacements) {
            std::regex placeholder("\\$\\{" + pair.first + "\\}");
            inputFileName = std::regex_replace(inputFileName, placeholder, pair.second);
        }

        std::string currentModelPathName = path;
        if (currentModelPathName.back() != '/')
            currentModelPathName.append("/");
        currentModelPathName.append(inputFileName);
//        std::cout << currentModelPathName << std::endl;

        return currentModelPathName;
    }
};

int main(int argc, char** argv){
    MPI_Init(&argc, &argv);
    S2OrderParamsCalculator<KEnRef_Real_t>::calc();
}

