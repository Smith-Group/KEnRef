/*
 * IoUtils.h
 *
 *  Created on: Oct 13, 2022
 *      Author: aalhossary
 */

#ifndef IOUTILS_H_
#define IOUTILS_H_

#include <regex>
#include <vector>
#include <string>
#include <map>
#include <functional>
#include <fstream>

//#include "../config/KEnRefConfig.h"
#include "core/KEnRef.h"
#include "core/Table.h"

class IoUtils final {
    const inline static std::map<std::string, std::string> aa1_to_aa3{
        {"TRP", "W"}, {"PHE", "F"}, {"TYR", "Y"}, {"MET", "M"}, {"LEU", "L"}, {"ILE", "I"}, {"VAL", "V"},
        {"ALA", "A"}, {"GLY", "G"}, {"SER", "S"}, {"THR", "T"}, {"ARG", "R"}, {"LYS", "K"}, {"HIS", "H"},
        {"ASN", "N"}, {"GLN", "Q"}, {"ASP", "D"}, {"GLU", "E"}, {"PRO", "P"}, {"CYS", "C"},
    };
    const inline static std::map<std::string, std::string> aa3_to_aa1{
        {"W", "TRP"}, {"F", "PHE"}, {"Y", "TYR"}, {"M", "MET"}, {"L", "LEU"}, {"I", "ILE"}, {"V", "VAL"},
        {"A", "ALA"}, {"G", "GLY"}, {"S", "SER"}, {"T", "THR"}, {"R", "ARG"}, {"K", "LYS"}, {"H", "HIS"},
        {"N", "ASN"}, {"Q", "GLN"}, {"D", "ASP"}, {"E", "GLU"}, {"P", "PRO"}, {"C", "CYS"},
    };

    static std::vector<std::vector<int> >
    grouping_mat_row_to_subset_idx(const Eigen::RowVectorX<int> &grouping_matrix);

public:
    IoUtils() = delete;

    template<typename T>
    static std::vector<std::vector<T> >
    read_uniform_table_of(std::istream &ins);


    static Table
    readTable(const std::string &fileName, bool has_columnNames = true, bool has_rowNames = false,
              const std::string &delimiter = "\\s*,\\s*", int max_rows = -1, bool tolerateMissingLastColumn = false);

    static Table
    readTable(std::ifstream &instream, bool has_columnNames, bool has_rowNames = false,
              const std::string &delimiter = "\\s*,\\s*", int max_rows = -1, bool tolerateMissingLastColumn = false);

    static std::vector<std::string> splitCSVLine(const std::string &line, const std::string &delimiter_pattern,
                                                 char quote_char = '"');

    static std::map<std::string, std::string> readParams(const std::string &fileName);

    static Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>
    subset_idx_to_grouping_mat(const std::vector<std::vector<std::vector<int> > > &multipleGroupings);

    static std::vector<std::vector<std::vector<int> > >
    grouping_mat_to_subset_idx(const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> &grouping_matrix);

    static std::map<std::string, std::string>
    readParams(std::istream &paramsFileStream);

    static int getEnvParam(const std::string &paramName, int defaultValue);

    static long getEnvParam(const std::string &paramName, long defaultValue);

    static std::string getEnvParam(const char *paramName, const char *defaultValue);

    static std::string getEnvParam(const std::string &paramName, const char *defaultValue);

    static std::string getEnvParam(const std::string &paramName, const std::string &defaultValue);

    template<typename KEnRef_Real_t>
    static KEnRef_Real_t getEnvParam(const std::string &paramName, KEnRef_Real_t defaultValue);


    static std::string
    strip_enclosing_quotes(const std::string &str, char delim = '"');

    static std::tuple<
        std::vector<std::string>, std::vector<std::string>, std::vector<int> >
    read_noe_table(std::istream &instream, bool has_header = true);

    std::vector<std::string>
    static split(const std::string &str, const std::string &delim = "\\s+");

    static std::map<std::string, std::vector<std::string> > read_noe_groups(std::istream &ins);

    static std::vector<int> getGmxNdxGroup(const std::string &filename, const std::string &groupName);

    static std::map<std::string, std::vector<int> > getAllGmxNdxGroups(const std::string &filename);

    template<typename TYPE>
    static void printVector(const std::vector<TYPE> &vec);

    static std::string padWithZeros(int value, int width);

    template<typename KEnRef_Real_t>
    static CoordsMatrixType<KEnRef_Real_t>
    extractCoords(const std::vector<int> &atomIndices, bool indicesOneBased,
                  std::map<int, Eigen::RowVector3<KEnRef_Real_t> > &allAtomCoords, bool mapOneBased);

    template<typename retMapKey, typename retMapValue>
    static std::map<retMapKey, retMapValue>
    getAtomMappingFromPdb(const std::string &pdbFilename,
                          const std::function<void(std::map<retMapKey, retMapValue> &ret, const std::smatch &sm)> &
                          mappingFunc);

    template<typename retMapKey, typename retMapValue>
    static std::vector<std::map<retMapKey, retMapValue> >
    getMultipleAtomMappingFromPdb(const std::string &pdbFilename,
                                  const std::function<void (std::map<retMapKey, retMapValue> &ret, const std::smatch &sm)> &mappingFunc);

    static std::string normalizeName(std::string atomId, bool lowerNameRanks = true);

    static bool isNotPrepared(const std::string &atomName);

    static bool should_handleNames(const std::map<std::string, int> &atomNameMapping);

    inline const static auto HB2_MET = std::regex("HB2.MET");
    inline const static auto HB3_MET = std::regex("HB3.MET");
    inline const static auto HG2_MET = std::regex("HG2.MET");
    inline const static auto HG3_MET = std::regex("HG3.MET");
    inline const static auto HB2_GLN = std::regex("HB2.GLN");
    inline const static auto HB3_GLN = std::regex("HB3.GLN");
    inline const static auto HG2_GLN = std::regex("HG2.GLN");
    inline const static auto HG3_GLN = std::regex("HG3.GLN");
    inline const static auto HB2_GLU = std::regex("HB2.GLU");
    inline const static auto HB3_GLU = std::regex("HB3.GLU");
    inline const static auto HG2_GLU = std::regex("HG2.GLU");
    inline const static auto HG3_GLU = std::regex("HG3.GLU");
    inline const static auto HB2_PHE = std::regex("HB2.PHE");
    inline const static auto HB3_PHE = std::regex("HB3.PHE");
    inline const static auto HB2_LYS = std::regex("HB2.LYS");
    inline const static auto HB3_LYS = std::regex("HB3.LYS");
    inline const static auto HG2_LYS = std::regex("HG2.LYS");
    inline const static auto HG3_LYS = std::regex("HG3.LYS");
    inline const static auto HD2_LYS = std::regex("HD2.LYS");
    inline const static auto HD3_LYS = std::regex("HD3.LYS");
    inline const static auto HE2_LYS = std::regex("HE2.LYS");
    inline const static auto HE3_LYS = std::regex("HE3.LYS");
    inline const static auto HB2_LEU = std::regex("HB2.LEU");
    inline const static auto HB3_LEU = std::regex("HB3.LEU");
    inline const static auto HA2_GLY = std::regex("HA2.GLY");
    inline const static auto HA3_GLY = std::regex("HA3.GLY");
    inline const static auto HB2_PRO = std::regex("HB2.PRO");
    inline const static auto HB3_PRO = std::regex("HB3.PRO");
    inline const static auto HG2_PRO = std::regex("HG2.PRO");
    inline const static auto HG3_PRO = std::regex("HG3.PRO");
    inline const static auto HD2_PRO = std::regex("HD2.PRO");
    inline const static auto HD3_PRO = std::regex("HD3.PRO");
    inline const static auto HB2_SER = std::regex("HB2.SER");
    inline const static auto HB3_SER = std::regex("HB3.SER");
    inline const static auto HB2_ASP = std::regex("HB2.ASP");
    inline const static auto HB3_ASP = std::regex("HB3.ASP");
    inline const static auto HB2_ASN = std::regex("HB2.ASN");
    inline const static auto HB3_ASN = std::regex("HB3.ASN");
    inline const static auto HB2_ARG = std::regex("HB2.ARG");
    inline const static auto HB3_ARG = std::regex("HB3.ARG");
    inline const static auto HG2_ARG = std::regex("HG2.ARG");
    inline const static auto HG3_ARG = std::regex("HG3.ARG");
    inline const static auto HD2_ARG = std::regex("HD2.ARG");
    inline const static auto HD3_ARG = std::regex("HD3.ARG");
    inline const static auto HB2_TYR = std::regex("HB2.TYR");
    inline const static auto HB3_TYR = std::regex("HB3.TYR");
    inline const static auto HB2_HIS = std::regex("HB2.HIS");
    inline const static auto HB3_HIS = std::regex("HB3.HIS");
    inline const static auto HG12_ILE = std::regex("HG12.ILE");
    inline const static auto HG13_ILE = std::regex("HG13.ILE");

    inline const static auto UNPREPARED_NAMES_MASK = std::regex("(HA3.GLY)|"
        "(HB3.(PHE|LEU|SER|ASP|ASN|TYR|HIS))|"
        "((HB3|HG3).(MET|GLN|GLU))|"
        "((HB3|HG3|HD3).(ARG|PRO))|"
        "((HB3|HG3|HD3|HE3).(LYS))|"
        "(HG13.ILE)");

    //This method removes AltLoc and ChainID and keeps only the last one of each
    static void fill_atomId_to_index_Map(std::map<std::string, int> &ret, const std::smatch &sm);

    template<typename KEnRef_Real>
    static void fill_atomIndex1_to_coords_Map(std::map<int, Eigen::RowVector3<KEnRef_Real> > &ret,
                                              const std::smatch &sm);

    static std::vector<std::tuple<int, int> > atomNamePairs_2_atomIdPairs(
        const std::vector<std::tuple<std::string, std::string>> &atomName_pairs,
        const std::map<std::string, int> &atomNames_2_atomIds);


    static std::vector<std::tuple<int, int> >
    readAtomIdPairs(const std::string &fileName);

    inline static const std::string &aa3(const std::string &aa1);

    inline static const std::string &aa1(const std::string &aa3);

    template<typename Scalar> static Scalar convertValue(const std::string& str);
    template<typename Scalar> static bool tryConvertValue(const std::string& str, Scalar& out) noexcept;
};
#endif /* IOUTILS_H_ */
