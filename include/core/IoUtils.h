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
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <functional>
#include <Eigen/Core>

#include "../config/KEnRefConfig.h"

class IoUtils final {
private:

public:
    IoUtils()=delete;
    template<typename T>
	static std::vector<std::vector<T>>
	read_uniform_table_of(std::istream &ins);

	static std::tuple<std::vector<std::string>, std::vector<std::vector<std::string>>>
	readTable(const std::string& fileName, bool has_header=true, int max_rows=-1);

    static std::map<std::string, std::string>
    readParams(const std::string& fileName);

    static std::map<std::string, std::string>
    readParams(std::istream &paramsFileStream);

	static std::string
	strip_enclosing_quotoes(const std::string& str, char delim = '\"');

	static std::tuple<
		std::vector<std::string>,
		std::vector<std::string>,
		std::vector<int>>
	read_noe_table(std::istream &instream, bool has_header=true);

	std::vector<std::string>
	static split(const std::string& str, const std::string& delim = "\\s+");

	static std::map<std::string, std::vector<std::string>>
	read_noe_groups(std::istream &ins);
	static std::vector<int> getGmxNdxGroup(const std::string& filename, const std::string& groupName);
	static std::map<std::string, std::vector<int>> getAllGmxNdxGroups(const std::string& filename);

	template<typename TYPE>
	static void printVector(const std::vector<TYPE>& vec);

    template<typename retMapKey, typename retMapValue>
    static std::map<retMapKey, retMapValue>
    getAtomMappingFromPdb(const std::string &pdbFilename, std::function<void(std::map<retMapKey, retMapValue> &ret, const std::smatch &sm)> mappingFunc);
    static std::string &normalizeName(std::string &atomId, bool lowerNameRanks = true);

    static bool isNotPrepared(const std::string &atomName);

    inline const static std::regex HB2_MET = std::regex("HB2.MET");
    inline const static std::regex HB3_MET = std::regex("HB3.MET");
    inline const static std::regex HG2_MET = std::regex("HG2.MET");
    inline const static std::regex HG3_MET = std::regex("HG3.MET");
    inline const static std::regex HB2_GLN = std::regex("HB2.GLN");
    inline const static std::regex HB3_GLN = std::regex("HB3.GLN");
    inline const static std::regex HG2_GLN = std::regex("HG2.GLN");
    inline const static std::regex HG3_GLN = std::regex("HG3.GLN");
    inline const static std::regex HB2_GLU = std::regex("HB2.GLU");
    inline const static std::regex HB3_GLU = std::regex("HB3.GLU");
    inline const static std::regex HG2_GLU = std::regex("HG2.GLU");
    inline const static std::regex HG3_GLU = std::regex("HG3.GLU");
    inline const static std::regex HB2_PHE = std::regex("HB2.PHE");
    inline const static std::regex HB3_PHE = std::regex("HB3.PHE");
    inline const static std::regex HB2_LYS = std::regex("HB2.LYS");
    inline const static std::regex HB3_LYS = std::regex("HB3.LYS");
    inline const static std::regex HG2_LYS = std::regex("HG2.LYS");
    inline const static std::regex HG3_LYS = std::regex("HG3.LYS");
    inline const static std::regex HD2_LYS = std::regex("HD2.LYS");
    inline const static std::regex HD3_LYS = std::regex("HD3.LYS");
    inline const static std::regex HE2_LYS = std::regex("HE2.LYS");
    inline const static std::regex HE3_LYS = std::regex("HE3.LYS");
    inline const static std::regex HB2_LEU = std::regex("HB2.LEU");
    inline const static std::regex HB3_LEU = std::regex("HB3.LEU");
    inline const static std::regex HA2_GLY = std::regex("HA2.GLY");
    inline const static std::regex HA3_GLY = std::regex("HA3.GLY");
    inline const static std::regex HB2_PRO = std::regex("HB2.PRO");
    inline const static std::regex HB3_PRO = std::regex("HB3.PRO");
    inline const static std::regex HG2_PRO = std::regex("HG2.PRO");
    inline const static std::regex HG3_PRO = std::regex("HG3.PRO");
    inline const static std::regex HD2_PRO = std::regex("HD2.PRO");
    inline const static std::regex HD3_PRO = std::regex("HD3.PRO");
    inline const static std::regex HB2_SER = std::regex("HB2.SER");
    inline const static std::regex HB3_SER = std::regex("HB3.SER");
    inline const static std::regex HB2_ASP = std::regex("HB2.ASP");
    inline const static std::regex HB3_ASP = std::regex("HB3.ASP");
    inline const static std::regex HB2_ASN = std::regex("HB2.ASN");
    inline const static std::regex HB3_ASN = std::regex("HB3.ASN");
    inline const static std::regex HB2_ARG = std::regex("HB2.ARG");
    inline const static std::regex HB3_ARG = std::regex("HB3.ARG");
    inline const static std::regex HG2_ARG = std::regex("HG2.ARG");
    inline const static std::regex HG3_ARG = std::regex("HG3.ARG");
    inline const static std::regex HD2_ARG = std::regex("HD2.ARG");
    inline const static std::regex HD3_ARG = std::regex("HD3.ARG");
    inline const static std::regex HB2_TYR = std::regex("HB2.TYR");
    inline const static std::regex HB3_TYR = std::regex("HB3.TYR");
    inline const static std::regex HB2_HIS = std::regex("HB2.HIS");
    inline const static std::regex HB3_HIS = std::regex("HB3.HIS");
    inline const static std::regex HG12_ILE = std::regex("HG12.ILE");
    inline const static std::regex HG13_ILE = std::regex("HG13.ILE");

    inline const static std::regex UNPREPARED_NAMES_MASK = std::regex("(HA3.GLY)|"
                                                                 "(HB3.(PHE|LEU|SER|ASP|ASN|TYR|HIS))|"
                                                                 "((HB3|HG3).(MET|GLN|GLU))|"
                                                                 "((HB3|HG3|HD3).(ARG|PRO))|"
                                                                 "((HB3|HG3|HD3|HE3).(LYS))|"
                                                                 "(HG13.ILE)");

    //This method removes AltLoc and ChainID and keeps only the last one of each
    static void fill_atomId_to_index_Map(std::map<std::string, int> &ret, const std::smatch &sm);

    template<typename KEnRef_Real>
    static void fill_atomIndex1_to_coords_Map(std::map<int, Eigen::RowVector3<KEnRef_Real>> &ret, const std::smatch &sm);
};

#endif /* IOUTILS_H_ */
