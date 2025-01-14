#include <filesystem>
#include <regex>
#include <fstream>
#include <string>
#include <system_error>
#include <Eigen/Core>
#include "core/IoUtils.h"

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
std::tuple<std::vector<std::string>, std::vector<std::vector<std::string>>> IoUtils::readTable(const std::string& fileName, bool has_header, const std::string& delimiter, int max_rows){
    std::ifstream instream(fileName);
    if (!instream.is_open()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return {{}, {}};
    }
    return readTable(instream, has_header, delimiter, max_rows);
}

std::tuple<std::vector<std::string>, std::vector<std::vector<std::string>>>
IoUtils::readTable(std::ifstream &instream, bool has_header, const std::string &delimiter, int max_rows) {
    std::vector<std::string> header{};
    std::vector<std::vector<std::string>> data{};
    bool header_consumed = false;
    std::string line;
    int counter = 0;
    //TODO create a unit test to validate the effect when max_rows is 0, -1, or a positive int
    while(instream.peek() != EOF && max_rows == -1 || counter < max_rows){
        std::getline(instream, line);
        auto tokens = split(line, delimiter);
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
            counter++;
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

std::regex atomRecordTemplate {"^((ATOM  )|(HETATM))([0-9 ]{5}) (.{15})   ([0-9 .-]{8})([0-9 .-]{8})([0-9 .-]{8}).+$" };

template<typename retMapKey, typename retMapValue>
std::map<retMapKey, retMapValue>
IoUtils::getAtomMappingFromPdb(const std::string &pdbFilename, std::function<void(std::map<retMapKey, retMapValue> &ret, const std::smatch &sm)> mappingFunc){
    std::map<retMapKey, retMapValue> ret = {};

    std::ifstream pdbFileStream(pdbFilename);
    if (pdbFileStream.fail()) {
        std::cerr << "ERROR: Could not open " << pdbFilename << "\n";
        throw std::runtime_error("File not found: " + pdbFilename);
    }
    std::string line;
    std::smatch sm;
    while (std::getline(pdbFileStream, line)) {
        if (line.empty()) continue;
        if(std::regex_match(line, sm, atomRecordTemplate)){
//            fill_atomId_to_index_Map(ret, sm);
            mappingFunc(ret, sm);
        }
    }
    pdbFileStream.close();
//    std::cerr << "number of items in the map = " << ret.size() << std::endl;
    return ret; //std::move(&ret);
}

void IoUtils::fill_atomId_to_index_Map(std::map<std::string, int> &ret, const std::smatch &sm){
    int atomIndex;
    std::string atomId = sm[5];
    normalizeName(atomId, false);
    std::istringstream iss(sm[4]);
    iss >> atomIndex;
    ret[atomId] = atomIndex;
}

template<typename KEnRef_Real>
void IoUtils::fill_atomIndex1_to_coords_Map(std::map<int, Eigen::RowVector3<KEnRef_Real>> &ret, const std::smatch &sm){
    int atomIndex1;
    KEnRef_Real x, y, z;
    std::istringstream iss1(sm[4]);
    iss1 >> atomIndex1;
    std::istringstream issX(sm[6]);
    issX >> x;
    std::istringstream issY(sm[7]);
    issY >> y;
    std::istringstream issZ(sm[8]);
    issZ >> z;
    ret[atomIndex1] = std::move(Eigen::RowVector3<KEnRef_Real>(x, y, z));
}


std::map<std::string, std::string> IoUtils::readParams(const std::string &fileName) {
    std::ifstream paramsFileStream(fileName);
    return readParams(paramsFileStream);
}

std::map<std::string, std::string> IoUtils::readParams(std::istream &paramsFileStream) {
    std::string line;
    std::regex recordTemplate {R"(^\s*(.+?)\s*=\s*(\S.*?)\s*(#.*)?)" };
    std::smatch sm;
    std::map<std::string, std::string> ret{};
    while (std::getline(paramsFileStream, line)) {
        if (line.empty()) continue;
        if (std::regex_match(line, sm, recordTemplate)) {
            const std::string &key = sm[1].str();
            if(key[0] == '#') continue;
            const std::string &value = strip_enclosing_quotoes(sm[2].str());
            ret[key] = value;
        }
    }
    return ret;
}

int IoUtils::getEnvParam(const std::string& paramName, int defaultValue){
    if (const char *pSEnvParam = std::getenv(paramName.c_str())) {
        int retValue = std::atoi(pSEnvParam);
        std::cout << paramName << " is: " << retValue << '\n';
        return retValue;
    } else {
        std::cout << "No "<< paramName << " identified. Will use default value of " << defaultValue << std::endl;
        return defaultValue;
    }
}
long IoUtils::getEnvParam(const std::string& paramName, long defaultValue){
    if (const char *pSEnvParam = std::getenv(paramName.c_str())) {
        long retValue = std::atol(pSEnvParam);
        std::cout << paramName << " is: " << retValue << '\n';
        return retValue;
    } else {
        std::cout << "No "<< paramName << " identified. Will use default value of " << defaultValue << std::endl;
        return defaultValue;
    }
}

std::string IoUtils::getEnvParam(const char *paramName, const char *defaultValue){
    return getEnvParam(std::string(paramName), std::string(defaultValue));
}
std::string IoUtils::getEnvParam(const std::string& paramName, const char *defaultValue){
    return getEnvParam(paramName, std::string(defaultValue));
}
std::string IoUtils::getEnvParam(const std::string& paramName, const std::string& defaultValue){
    std::string retValue = defaultValue;
    if (const char *pSEnvParam = std::getenv(paramName.c_str())) {
        retValue = pSEnvParam;
        std::cout << paramName << " is: " << retValue << '\n';
    } else {
        std::cout << "No "<< paramName << " identified. Will use default value of " << defaultValue << std::endl;
    }
    return retValue;
}
template<typename KEnRef_Real_t>
KEnRef_Real_t IoUtils::getEnvParam(const std::string& paramName, KEnRef_Real_t defaultValue){
    if (const char *pSEnvParam = std::getenv(paramName.c_str())) {
        std::stringstream sstream(pSEnvParam);
        KEnRef_Real_t retValue;
        sstream >> retValue;
        std::cout << paramName << " is: " << retValue << '\n';
        return retValue;
    } else {
        std::cout << "No "<< paramName << " identified. Will use default value of " << defaultValue << std::endl;
        return defaultValue;
    }
}


std::vector<std::tuple<int, int>> IoUtils::readAtomIdPairs(const std::string &fileName) {
    std::ifstream atomIdPairsFileStream(fileName);
    auto tempAtomIdPairsTable = IoUtils::read_uniform_table_of<int>(atomIdPairsFileStream);
    std::vector<std::tuple<int, int>> atomIdPairs;
    for (auto row : tempAtomIdPairsTable) {
        atomIdPairs.emplace_back(row[0], row[1]);
    }
    std::cout << "Atom ID Pairs (" << atomIdPairs.size() << "):";
    std::cout << "<< \n";
    return atomIdPairs;
}

std::string IoUtils::padWithZeros(int value, int width) {
    std::ostringstream oss;
    oss << std::setw(width) << std::setfill('0') << value;
    return oss.str();
}

template<typename KEnRef_Real_t>
CoordsMatrixType<KEnRef_Real_t>
IoUtils::extractCoords(const std::vector<int> &atomIndices, bool indicesOneBased,
                       std::map<int, Eigen::RowVector3<KEnRef_Real_t>> &allAtomCoords, bool mapOneBased) {
    auto coords = CoordsMatrixType<KEnRef_Real_t>(atomIndices.size(), 3);
    int delta = 0;
    if(indicesOneBased ^ mapOneBased)
        indicesOneBased ? delta-- : delta++;


    for (int i = 0; i < atomIndices.size(); ++i) {
        int key = atomIndices[i] + delta;
        coords(i, Eigen::all) = allAtomCoords[key];
    }
    return coords;
}

template<typename TYPE>
void IoUtils::printVector(const std::vector<TYPE>& vec){
    for (const TYPE &val: vec)
        std::cout << val << " ";
    std::cout << std::endl;
}

///////// Template Declarations /////////////////////////

template std::vector<std::vector<int>>IoUtils::read_uniform_table_of(std::istream&);
template std::vector<std::vector<double>>IoUtils::read_uniform_table_of(std::istream&);

template float IoUtils::getEnvParam(const std::string& paramName, float defaultValue);
template double IoUtils::getEnvParam(const std::string& paramName, double defaultValue);

template CoordsMatrixType<float>
IoUtils::extractCoords(const std::vector<int> &atomIndices, bool indicesOneBased,
                       std::map<int, Eigen::RowVector3<float>> &allAtomCoords, bool mapOneBased);
template CoordsMatrixType<double>
IoUtils::extractCoords(const std::vector<int> &atomIndices, bool indicesOneBased,
                       std::map<int, Eigen::RowVector3<double>> &allAtomCoords, bool mapOneBased);

template void IoUtils::printVector(const std::vector<bool> &vec);
template void IoUtils::printVector(const std::vector<int> &vec);
template void IoUtils::printVector(const std::vector<float> &vec);
template void IoUtils::printVector(const std::vector<double> &vec);
template void IoUtils::printVector(const std::vector<std::string> &vec);

template std::map<std::string, int> IoUtils::getAtomMappingFromPdb(const std::string& pdbFilename, std::function<void(std::map<std::string, int> &, const std::smatch &)> mappingFunc);
template void IoUtils::fill_atomIndex1_to_coords_Map<float>(std::map<int, Eigen::RowVector3<float>> &ret, const std::smatch &sm);
template void IoUtils::fill_atomIndex1_to_coords_Map<double>(std::map<int, Eigen::RowVector3<double>> &ret, const std::smatch &sm);
template std::map<int, Eigen::RowVector3<float>> IoUtils::getAtomMappingFromPdb<int, Eigen::RowVector3<float>>(const std::string& pdbFilename, std::function<void(std::map<int, Eigen::RowVector3<float>> &, const std::smatch &)> mappingFunc);
template std::map<int, Eigen::RowVector3<double>> IoUtils::getAtomMappingFromPdb<int, Eigen::RowVector3<double>>(const std::string& pdbFilename, std::function<void(std::map<int, Eigen::RowVector3<double>> &, const std::smatch &)> mappingFunc);
