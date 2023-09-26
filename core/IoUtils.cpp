#include <filesystem>
#include <regex>
#include <fstream>
#include <string>
#include <system_error>
#include "IoUtils.h"

template<typename T>
std::vector<std::vector<T>> IoUtils::read_uniform_table_of(std::istream &ins) {
	std::vector<std::vector<T>> result;
	std::string s;
	while (getline(ins, s)) {
		std::vector<T> row;
		std::istringstream ss(s);
		T value;
		while (ss >> value)
			row.emplace_back(value);
		result.emplace_back(row);
	}
	return result;
}

//caution: this method does not handle the case when there is a comma inside a string
std::tuple<std::vector<std::string>, std::vector<std::vector<std::string>>> IoUtils::readTable(const std::string& fileName, bool has_header){
	std::vector<std::string> header{};
	std::vector<std::vector<std::string>> data{};
	bool header_consumed = false;
	std::string line;

    std::ifstream instream(fileName);
    if (!instream.is_open()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return {header, data};
    }
	while(instream.peek() != EOF){
		std::getline(instream, line);
		auto tokens = split(line, ",\\s*");
		for(size_t i = 0; i < tokens.size(); ++i){ // NOLINT(modernize-loop-convert)
			tokens[i] = strip_enclosing_quotoes(tokens[i], '\"');
		}
		if(has_header && ! header_consumed){
			//use consumed header
			header = tokens;
//			for(auto colName: header){std::cout << colName << '\t';} std::cout << std::endl;
			header_consumed = true;
		}else{
			data.emplace_back(tokens);
//			for(std::string value: data.back()){std::cout << value << '\t';} std::cout << std::endl;
		}
	}

//    std::cerr << "------- before table tuple" << std::endl;
//    std::tuple<std::vector<std::string>, std::vector<std::vector<std::string>>> ret = std::make_tuple(header, data);
//    std::cerr << "------- after table tuple" << std::endl;
//    return std::move(& ret);
    return {header, data};
}

std::string IoUtils::strip_enclosing_quotoes(const std::string &str, char delim) {
	// unfortunately, lookbehind does not work
	// std::string ret = std::regex_replace(str, std::regex("((?=^)\")|(?<!\\)(\"$)"), "");
	if (str.empty() || str.front() != delim) //no starting double quote
		return str;
	if ((str.back() != delim) || (str.length() > 1 && (str[str.length() - 2] == '\\'))) //no matching end double quote (or escaped)
		return str;
	else
		return str.substr(1, str.length()-2);
}

std::tuple<std::vector<std::string>, std::vector<std::string>,
		std::vector<int>> IoUtils::read_noe_table(std::istream &instream, bool skip_header) {
	std::vector<std::tuple<std::string, std::string, int>> temp { };

	std::string line;
	std::regex lineTemplate {"(.+)\t(.+)\t(.+)" };
	std::smatch sm;
	bool header_consumed = false;
	while (std::getline(instream, line)) {
		std::regex_match(line, sm, lineTemplate);
		std::string g1 = sm[1];
		std::string g2 = sm[2];
		std:: string val_str = sm[3];
		int val_int = -1;
		if(skip_header && ! header_consumed){
			//TODO use consumed header
//			std::cout << g1 << '\t' << g2 << '\t' << val_str << std::endl;
			header_consumed = true;
		}else{
			g1 = strip_enclosing_quotoes(g1);
			g2 = strip_enclosing_quotoes(g2);
			val_int = std::stoi(val_str);
//			std::cout << g1 << '\t' << g2 << '\t' << val_int << std::endl;
			temp.emplace_back(g1, g2, val_int);
		}
	}
	std::vector<std::string> group1{}, group2{};
	std::vector<int> values{};
	for(const auto& [g1, g2, v] : temp){
		group1.emplace_back(g1);
		group2.emplace_back(g2);
		values.emplace_back(v);
	}
	//TODO change this tuple of vectors into a vector of tuples
	std::tuple<std::vector<std::string>, std::vector<std::string>,
			std::vector<int>> ret {group1, group2, values};
	return ret;
}


std::vector<std::string>
IoUtils::split(const std::string& str, const std::string& delim){
	std::regex delim_re(delim);
	std::vector<std::string> ret;
	std::sregex_token_iterator iter(str.begin(), str.end(), delim_re, -1);
	std::sregex_token_iterator end;

	while (iter != end) {
//		ret.emplace_back(*iter);
		ret.emplace_back(iter->str());
		++iter;
	}

//    while(next_delim_pos != std::string::npos){
//        ret.emplace_back(str.begin() + startIndex, str.begin() + next_delim_pos);
//        startIndex = next_delim_pos + delim.size();
//        next_delim_pos = str.find(delim, startIndex);
//    }
//    if(startIndex != str.size())
//        ret.emplace_back(str.begin()+startIndex, str.end());
    return ret;
}


std::map<std::string, std::vector<std::string>>
IoUtils::read_noe_groups(std::istream &instream){

	std::map<std::string, std::vector<std::string>> ret{};
	std::string line;
	std::regex lineTemplate { R"(^\s*(.*?)\s*=\s*(.*?)(\s*,?\s*)$)" };
	std::smatch sm;

	while(std::getline(instream, line)){
		if(! std::regex_match(line, sm, lineTemplate))
			continue;
		std::string key = strip_enclosing_quotoes(sm[1], '`');
		std::string temp_val = sm[2];
		std::vector<std::string> val;
		if(temp_val.rfind("c(", 0) == 0){ // str.startwith() is available in C++20 at least
			auto tokens = split(temp_val.substr(2, temp_val.length()-3), ",\\s*");
			for(auto i = 0; i < tokens.size(); ++i){ // NOLINT(modernize-loop-convert)
				tokens[i] = strip_enclosing_quotoes(tokens[i], '\"');
			}
			val = tokens;
		}else{
			val = {strip_enclosing_quotoes(temp_val, '\"')};
		}
		ret.insert({key, val});
	}
	return ret;
}


std::vector<int>
IoUtils::getGmxNdxGroup(const std::string& filename, const std::string& groupName){
    std::vector<int> ret;
    std::ifstream indexFile(filename);
    if (!indexFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return ret;
    }
    std::string line;
    bool groupFound = false;
    while (std::getline(indexFile, line)) {
        if (line.empty()) continue;
        if (line[0] == '[') {
            groupFound = (line.find(groupName) != std::string::npos);
            continue;
        }
        if (!groupFound) continue;
        int index;
        std::istringstream iss(line);
        while (iss >> index) {
            ret.push_back(index - 1);
        }
    }
    indexFile.close();
    return ret;
}

std::map<std::string, std::vector<int>>
IoUtils::getAllGmxNdxGroups(const std::string& filename){
    std::map<std::string, std::vector<int>> ret;
    std::ifstream indexFile(filename);
    if (!indexFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return ret;
    }
    std::string line;
	std::vector<int> indices;//FIXME is it safe to declare it this way?
    while (std::getline(indexFile, line)) {
        if (line.empty()) continue;
        if (line[0] == '[') {
        	auto closing = line.find(']');
        	if (closing != std::string::npos && closing > 0) {
					std::string groupName;
        	        groupName =line.substr(2, closing-2); // TODO trim the whitespace better than hard coding (2 , closing -2)
        	        ret[groupName] = indices;
//        	        std::cout << "--------------------{{" << groupName << "}}" << std::endl;
					continue;
        	}else {
                std::cerr << "error parsing line: " << line << std::endl;
                return ret;
			}
        }

        int index;
        std::istringstream iss(line);
        while (iss >> index) {
            indices.push_back(index - 1);
        }
    }
    indexFile.close();
    return ret;
}

bool IoUtils::isNotPrepared(const std::string& atomName) {
    return std::regex_search(atomName, IoUtils::UNPREPARED_NAMES_MASK);
}

/**
 * IMPORTANT: Notice that we assume for this function to lower name ranks correctly that the function is called
 * sequentially with atom names in lexical order (e.g. "HB2 MET" is called before "HB3 MET" and not vise versa).
 */
std::string& IoUtils::normalizeName(std::string &atomId, bool lowerNameRanks) {
	atomId[4] = ' '; //remove alternate location
	atomId[9] = ' '; //remove chain ID
	if (lowerNameRanks) {
		atomId = std::regex_replace(atomId, HB2_MET, "HB1 MET");
		atomId = std::regex_replace(atomId, HB3_MET, "HB2 MET");
		atomId = std::regex_replace(atomId, HG2_MET, "HG1 MET");
		atomId = std::regex_replace(atomId, HG3_MET, "HG2 MET");

		atomId = std::regex_replace(atomId, HB2_GLN, "HB1 GLN");
		atomId = std::regex_replace(atomId, HB3_GLN, "HB2 GLN");
		atomId = std::regex_replace(atomId, HG2_GLN, "HG1 GLN");
		atomId = std::regex_replace(atomId, HG3_GLN, "HG2 GLN");

		atomId = std::regex_replace(atomId, HB2_GLU, "HB1 GLU");
		atomId = std::regex_replace(atomId, HB3_GLU, "HB2 GLU");
		atomId = std::regex_replace(atomId, HG2_GLU, "HG1 GLU");
		atomId = std::regex_replace(atomId, HG3_GLU, "HG2 GLU");

		atomId = std::regex_replace(atomId, HB2_PHE, "HB1 PHE");
		atomId = std::regex_replace(atomId, HB3_PHE, "HB2 PHE");

		atomId = std::regex_replace(atomId, HB2_LYS, "HB1 LYS");
		atomId = std::regex_replace(atomId, HB3_LYS, "HB2 LYS");
		atomId = std::regex_replace(atomId, HG2_LYS, "HG1 LYS");
		atomId = std::regex_replace(atomId, HG3_LYS, "HG2 LYS");
		atomId = std::regex_replace(atomId, HD2_LYS, "HD1 LYS");
		atomId = std::regex_replace(atomId, HD3_LYS, "HD2 LYS");
		atomId = std::regex_replace(atomId, HE2_LYS, "HE1 LYS");
		atomId = std::regex_replace(atomId, HE3_LYS, "HE2 LYS");

		atomId = std::regex_replace(atomId, HB2_LEU, "HB1 LEU");
		atomId = std::regex_replace(atomId, HB3_LEU, "HB2 LEU");

		atomId = std::regex_replace(atomId, HA2_GLY, "HA1 GLY");
		atomId = std::regex_replace(atomId, HA3_GLY, "HA2 GLY");

		atomId = std::regex_replace(atomId, HB2_PRO, "HB1 PRO");
		atomId = std::regex_replace(atomId, HB3_PRO, "HB2 PRO");
		atomId = std::regex_replace(atomId, HG2_PRO, "HG1 PRO");
		atomId = std::regex_replace(atomId, HG3_PRO, "HG2 PRO");
		atomId = std::regex_replace(atomId, HD2_PRO, "HD1 PRO");
		atomId = std::regex_replace(atomId, HD3_PRO, "HD2 PRO");

		atomId = std::regex_replace(atomId, HB2_SER, "HB1 SER");
		atomId = std::regex_replace(atomId, HB3_SER, "HB2 SER");

		atomId = std::regex_replace(atomId, HB2_ASP, "HB1 ASP");
		atomId = std::regex_replace(atomId, HB3_ASP, "HB2 ASP");

		atomId = std::regex_replace(atomId, HB2_ASN, "HB1 ASN");
		atomId = std::regex_replace(atomId, HB3_ASN, "HB2 ASN");

		atomId = std::regex_replace(atomId, HB2_ARG, "HB1 ARG");
		atomId = std::regex_replace(atomId, HB3_ARG, "HB2 ARG");
		atomId = std::regex_replace(atomId, HG2_ARG, "HG1 ARG");
		atomId = std::regex_replace(atomId, HG3_ARG, "HG2 ARG");
		atomId = std::regex_replace(atomId, HD2_ARG, "HD1 ARG");
		atomId = std::regex_replace(atomId, HD3_ARG, "HD2 ARG");

		atomId = std::regex_replace(atomId, HB2_TYR, "HB1 TYR");
		atomId = std::regex_replace(atomId, HB3_TYR, "HB2 TYR");

		atomId = std::regex_replace(atomId, HB2_HIS, "HB1 HIS");
		atomId = std::regex_replace(atomId, HB3_HIS, "HB2 HIS");

		atomId = std::regex_replace(atomId, HG12_ILE, "HG11 ILE");
		atomId = std::regex_replace(atomId, HG13_ILE, "HG12 ILE");
	}
	return atomId;
}

std::map<std::string, int>
IoUtils::getAtomNameMappingFromPdb(const std::string& pdbFilename){
	std::map<std::string, int> ret = {};

//    std::filesystem::path pdbFilePath(pdbFilename);
//    std::ifstream pdbFileStream;
//    try {
//        pdbFileStream = std::ifstream (pdbFilename);
//    } catch (std::filesystem::filesystem_error& e) {
//        throw e;
//    }
    std::ifstream pdbFileStream(pdbFilename); //no need to try and catch here.

    std::string line;
    std::regex atomRecordTemplate {"^((ATOM  )|(HETATM))([0-9 ]{5}) (.{15}).+$" };
    std::smatch sm;
    int index;
    while (std::getline(pdbFileStream, line)) {
        if (line.empty()) continue;
        if(std::regex_match(line, sm, atomRecordTemplate)){
        	std::string atomId = sm[5];
			normalizeName(atomId, false);
        	std::istringstream iss(sm[4]);
        	iss >> index;
        	ret[atomId] = index;
        }
    }
    pdbFileStream.close();
    std::cerr << "number of items in the map = " << ret.size() << std::endl;
    return ret; //std::move(&ret);
}


void IoUtils::printVector(const std::vector<int>& vec){
	for(int val : vec)
		std::cout << val << " ";
	std::cout << std::endl;
}
void IoUtils::printVector(const std::vector<bool>& vec){
	for(bool val : vec)
		std::cout << std::boolalpha << val << " ";
	std::cout << std::endl;
}
void IoUtils::printVector(const std::vector<std::string>& vec){
	for(const std::string& val : vec)
		std::cout << val << " ";
	std::cout << std::endl;
}
template<typename TYPE>
void IoUtils::printVector(const std::vector<TYPE>& vec){
    for (const TYPE &val: vec)
        std::cout << val << " ";
    std::cout << std::endl;
}

///////////////////////////////////////////////////////
void test() {
	std::istringstream input(
			"1 2 3\n"
			"4 5 6\n"
			"7 8 9\n");
	auto table = IoUtils::read_uniform_table_of<int>(input);
	for (const auto& record : table) {
		for (auto field : record)
			std::cout << field << " ";
		std::cout << "\n";
	}
	std::cout << std::endl;
	///////////////////////////////////
	std::string strin, strout;
	strin = "\"doublequote in the beginning";
	strout = IoUtils::strip_enclosing_quotoes(strin);
	std::cout << strout << std::endl;
	strin = "doublequote at the end\"";
	strout = IoUtils::strip_enclosing_quotoes(strin);
	std::cout << strout << std::endl;
	strin = "doublequote in the (\") middle";
	strout = IoUtils::strip_enclosing_quotoes(strin);
	std::cout << strout << std::endl;
	strin = R"("escaped doublequote at the end\")";
	strout = IoUtils::strip_enclosing_quotoes(strin);
	std::cout << strout << std::endl;
	strin = "\"doublequote at both ends\"";
	strout = IoUtils::strip_enclosing_quotoes(strin);
	std::cout << strout << std::endl;
	strin = "`escaped back quote at the end\\`";
	strout = IoUtils::strip_enclosing_quotoes(strin);
	std::cout << strout << std::endl;
	strin = "`back quote at both ends`";
	strout = IoUtils::strip_enclosing_quotoes(strin, '`');
	std::cout << strout << std::endl;
	///////////////////////////////////
	std::ifstream in("res/noe_1.tsv");
	if (in) {
		const auto& [group1names, group2names, values] = IoUtils::read_noe_table(in);
		for(int i = 0; i < group1names.size(); ++i){
			std::cout<< group1names[i] << "\t| " << group2names[i] << "\t| " << values[i] << std::endl;
		}
	}

	std::ifstream noe_groups_file("../res/noe_groups.R");
	auto noe_groups = IoUtils::read_noe_groups(noe_groups_file);
	for(auto [key, val] : noe_groups){
		std::cout << "<" << key<< ">"  << " " << val.size() << ":";
		for(const auto& str: val){
			std::cout << " [" << str << ']';
		}
		std::cout << std::endl;
	}
}

#if _TEST_
int main() {
	test();
}
#endif
