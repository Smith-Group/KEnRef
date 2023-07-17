/*
 * IoUtils.h
 *
 *  Created on: Oct 13, 2022
 *      Author: aalhossary
 */

#ifndef IOUTILS_H_
#define IOUTILS_H_

#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <map>

#include "../config/KEnRefConfig.h"

class IoUtils final {
private:

public:
    IoUtils()=delete;
    template<typename T>
	static std::vector<std::vector<T>>
	read_uniform_table_of(std::istream &ins);

	static std::tuple<std::vector<std::string>, std::vector<std::vector<std::string>>>
	readTable(const std::string& fileName, bool has_header=true);

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

	static void printVector(const std::vector<int>& vec);
	static void printVector(const std::vector<bool>& vec);
	static void printVector(const std::vector<std::string>& vec);
	template<typename TYPE>
	static void printVector(const std::vector<TYPE>& vec);

	static std::map<std::string, int>
	getAtomNameMappingFromPdb(const std::string& filename); //This method removes AltLoc and ChainID and keeps only the last one of each
	static std::string& normalizeName(std::string &atomId, bool lowerMet=true);
};

#endif /* IOUTILS_H_ */
