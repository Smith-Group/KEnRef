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
	IoUtils();


public:
	template<typename T>
	static std::vector<std::vector<T>>
	read_uniform_table_of(std::istream &ins);

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
	static std::vector<int> getGmxNdxGroup(const std::string filename, const std::string groupName);
	static std::map<std::string, std::vector<int>> getAllGmxNdxGroups(const std::string filename);

	static void printVector(const std::vector<int>& vec);
	template<typename T>
	static void printVector(const std::vector<T, std::allocator<T>>& vec);
};

#endif /* IOUTILS_H_ */
