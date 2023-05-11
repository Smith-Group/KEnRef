#include <regex>
#include <fstream>
#include "IoUtils.h"

IoUtils::IoUtils() {
}

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
	std::regex tempelate { "\(.+)\t\(.+)\t\(.+)" };
	std::smatch sm;
	bool header_consumed = false;
	while (std::getline(instream, line)) {
		std::regex_match(line, sm, tempelate);
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
			temp.emplace_back(std::tuple<std::string, std::string, int> {g1, g2, val_int});
		}
	}
	std::vector<std::string> group1{}, group2{};
	std::vector<int> values{};
	for(auto entry : temp){
		auto[g1, g2, v] = entry;
		group1.emplace_back(g1);
		group2.emplace_back(g2);
		values.emplace_back(v);
	}
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
	std::regex tempelate { "^\\s*(.*?)\\s*=\\s*(.*?)(\\s*,?\\s*)$" };
	std::smatch sm;

	while(std::getline(instream, line)){
		if(! std::regex_match(line, sm, tempelate))
			continue;
		std::string key = strip_enclosing_quotoes(sm[1], '`');
		std::string temp_val = sm[2];
		std::vector<std::string> val;
		if(temp_val.rfind("c(", 0) == 0){ // str.startwith() is available in C++20 at least
			auto tokens = split(temp_val.substr(2, temp_val.length()-3), ",\\s*");
			for(size_t i = 0; i < tokens.size(); ++i){
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
IoUtils::getGmxNdxGroup(const std::string filename, const std::string groupName){
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
IoUtils::getAllGmxNdxGroups(const std::string filename){
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
        	int closing = line.find(']');
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

void IoUtils::printVector(const std::vector<int>& vec){
	for(int val : vec)
		std::cout << val << " ";
	std::cout << std::endl;
}
template<typename T>
void IoUtils::printVector(const std::vector<T, std::allocator<T>>& vec) {
	for (T val : vec) {
		std::cout << val << " ";
	}
	std::cout << std::endl;
}



///////////////////////////////////////////////////////
void test() {
	std::istringstream input(
			"1 2 3\n"
			"4 5 6\n"
			"7 8 9\n");
	auto table = IoUtils::read_uniform_table_of<int>(input);
	for (auto record : table) {
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
	strin = "\"escaped doublequote at the end\\\"";
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
		auto table = IoUtils::read_noe_table(in);
		auto[group1names, group2names, values] = table;
		for(int i = 0; i < group1names.size(); ++i){
			std::cout<< group1names[i] << "\t| " << group2names[i] << "\t| " << values[i] << std::endl;
		}
	}

	std::ifstream noe_groups_file("../res/noe_groups.R");
	auto noe_groups = IoUtils::read_noe_groups(noe_groups_file);
	for(auto [key, val] : noe_groups){
		std::cout << "<" << key<< ">"  << " " << val.size() << ":";
		for(auto str: val){
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
