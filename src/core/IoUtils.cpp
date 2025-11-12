//#include <filesystem>
//#include <regex>
//#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
//#include <Eigen/Core>
#include "core/IoUtils.h"

template<typename T>
std::vector<std::vector<T> > IoUtils::read_uniform_table_of(std::istream &ins) {
    std::vector<std::vector<T> > result;
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


std::vector<std::tuple<int, int>> IoUtils::atomNamePairs_2_atomIdPairs(
    const std::vector<std::tuple<std::string, std::string>> &atomName_pairs,
    const std::map<std::string, int> &atomNames_2_atomIds){

    auto atomId_pairs = std::vector<std::tuple<int, int>>(atomName_pairs.size());
    // Fill the vector using atomNames_2_atomIds
    for (int i = 0; i < atomName_pairs.size(); ++i) {
        auto [left, right] = atomName_pairs.at(i);
        // I use at() instead of operator[] to force an exception to be thrown
        atomId_pairs.at(i) = std::move(std::tuple{atomNames_2_atomIds.at(left), atomNames_2_atomIds.at(right)});
    }
    return atomId_pairs;
}



std::string IoUtils::strip_enclosing_quotes(const std::string &str, char delim) {
    // unfortunately, lookbehind does not work
    // std::string ret = std::regex_replace(str, std::regex("((?=^)\")|(?<!\\)(\"$)"), "");
    if (str.empty() || str.front() != delim) //no starting double quote
        return str;
    if ((str.back() != delim) || (str.length() > 1 && str[str.length() - 2] == '\\'))
        //no matching end double quote (or escaped)
        return str;
    else
        return str.substr(1, str.length() - 2);
}

std::vector<std::string>
IoUtils::splitCSVLine(const std::string &line, const std::string &delimiter_pattern, const char quote_char) {
    std::vector<std::string> tokens;
    std::string currentToken;
    bool inQuotes = false;
    size_t i = 0;

    const std::regex delimiter_regex(delimiter_pattern);

    while (i < line.size()) {
        char c = line[i];

        if (c == quote_char) {
            // Toggle quote state unless it's an escaped quote
            if (i > 0 && line[i - 1] == '\\') {
                currentToken.back() = quote_char; // replace \ with "
            } else {
                inQuotes = !inQuotes;
            }
            currentToken += c;
            ++i;
        }
        else if (!inQuotes) {
            std::smatch match;
            std::string remaining = line.substr(i);
            if (std::regex_search(remaining, match, delimiter_regex, std::regex_constants::match_continuous)) {
                tokens.push_back(strip_enclosing_quotes(currentToken, quote_char));
                currentToken.clear();
                i += match.length(); // skip delimiter
            } else {
                currentToken += c;
                ++i;
            }
        }
        else {
            currentToken += c;
            ++i;
        }
    }

    // Add the last token
    if (!currentToken.empty()) {
        tokens.push_back(strip_enclosing_quotes(currentToken, quote_char));
    }

    return tokens;
}

bool IoUtils::ends_with(const std::string &s, const std::string &suffix) {
    const std::size_t n = s.size();
    const std::size_t m = suffix.size();
    return n >= m && s.rfind(suffix, n - m) == n - m;
}

bool IoUtils::ext_xtc_or_trr(const std::string &path) {
    const std::string ext = std::filesystem::path(path).extension().string(); // e.g., ".xtc"
    return ext == ".xtc" || ext == ".trr";
}

KEnRef_Real_t IoUtils::pearsonCorrelation(const Eigen::VectorX<KEnRef_Real_t> &x, const Eigen::VectorX<KEnRef_Real_t> &y) {
    if (x.size() != y.size())
        throw std::invalid_argument("Vectors must be of the same length.");
    const KEnRef_Real_t mean_x = x.mean();
    const KEnRef_Real_t mean_y = y.mean();
    Eigen::VectorX<KEnRef_Real_t> diff_x = x.array() - mean_x;
    Eigen::VectorX<KEnRef_Real_t> diff_y = y.array() - mean_y;
    const KEnRef_Real_t numerator = (diff_x.array() * diff_y.array()).sum();
    const KEnRef_Real_t denominator = std::sqrt((diff_x.array().square().sum()) * (diff_y.array().square().sum()));
    return numerator / denominator;
}

std::vector<SpecDenData<KEnRef_Real_t>> IoUtils::load_spec_den_data(const std::string &experimentalDataFolder, const bool handleNames) {
    const std::vector<std::string> &spec_den_data_prefixes = KEnRef<KEnRef_Real_t>::spec_den_data_prefixes;
    std::vector<SpecDenData<KEnRef_Real_t>> spec_den_data_vector;
    spec_den_data_vector.reserve(spec_den_data_prefixes.size());
    for (const auto & spec_den_data_prefix : spec_den_data_prefixes) {
        //AtomPairs
        std::string atomPairAndSigmaFileName = experimentalDataFolder + spec_den_data_prefix + "_atom_pairs.csv";
        std::cout << atomPairAndSigmaFileName << std::endl;
        const Table& atomPairAndSigmaTable = IoUtils::readTable(atomPairAndSigmaFileName,true,false, "\\s*,\\s*", -1, true);
        std::vector<std::tuple<std::string, std::string>> atomPairs(atomPairAndSigmaTable.rowCount());
        for (int j = 0; j < atomPairAndSigmaTable.rowCount(); ++j) {
            auto atom1 = IoUtils::normalizeName(atomPairAndSigmaTable.at(j, 0), handleNames);
            auto atom2 = IoUtils::normalizeName(atomPairAndSigmaTable.at(j, 1), handleNames);
            atomPairs.at(j) = std::move(std::tuple<std::string, std::string>{ atom1, atom2 });
        }
        // sigma
        std::vector<KEnRef_Real_t> sigmasVec{};
        for (int row = 0; row < atomPairAndSigmaTable.rowCount(); ++row) {
            if (! atomPairAndSigmaTable.isRowComplete(row))
                break;
            const auto &valueStr = atomPairAndSigmaTable.at(row, 2);
            std::istringstream iss(valueStr);
            KEnRef_Real_t value;
            iss >> value;
            sigmasVec.emplace_back(value);
        }
        std::optional<NamedVector<KEnRef_Real_t>> sigma = NamedVector<KEnRef_Real_t>(sigmasVec.size());
        for (int j = 0; j < sigma->rows(); ++j) {
            sigma.value()(j, 0) = sigmasVec[j];
        }
        //multiple_grouping
        std::string multiple_grouping_fileName = experimentalDataFolder + spec_den_data_prefix+"_groupings.csv";
        const auto &grouping_matrix = NamedMatrix<int>(IoUtils::readTable(multiple_grouping_fileName, false, false).toNamedMatrix<int>().array() - 1);
        const auto &multiple_grouping = IoUtils::grouping_mat_to_subset_idx(grouping_matrix);
        //a_coef
        std::string aCoefFileName = experimentalDataFolder + spec_den_data_prefix+"_a_coef.csv";
        const auto &a_coef = IoUtils::readTable(aCoefFileName, true,false, "\\s*,\\s*", -1, false).toNamedMatrix<KEnRef_Real_t>();
        //lambda_coef
        std::string lambdaCoefFileName = experimentalDataFolder + spec_den_data_prefix+"_lambda_coef.csv";
        const auto &lambda_coef = IoUtils::readTable(lambdaCoefFileName, true,true, "\\s*,\\s*", -1, false).toNamedMatrix<KEnRef_Real_t>();

        // spec_den_data_vector.emplace_back(atomPairs, sigma, multiple_grouping, a_coef, lambda_coef);
        spec_den_data_vector.emplace_back(std::move(SpecDenData<KEnRef_Real_t>{atomPairs, sigma, multiple_grouping, a_coef, lambda_coef}));
    }
    return spec_den_data_vector;
}

Table
IoUtils::readTable(const std::string &fileName, const bool has_columnNames, const bool has_rowNames,
    const std::string &delimiter, const int max_rows, bool tolerateMissingLastColumn) {
    std::ifstream instream(fileName);
    if (!instream.is_open()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        throw std::runtime_error(std::string("Can't open file:").append(fileName));
    }
    return readTable(instream, has_columnNames, has_rowNames, delimiter, max_rows, tolerateMissingLastColumn);

}

Table
IoUtils::readTable(std::ifstream &instream, const bool has_columnNames, const bool has_rowNames,
    const std::string &delimiter, const int max_rows, bool tolerateMissingLastColumn) {

    std::vector<std::vector<std::string>> data;
    std::optional<std::vector<std::string>> colNames;
    std::optional<std::vector<std::string>> rowNames;
    std::string line;

    // Read column names if present
    if (has_columnNames && std::getline(instream, line)) {
        colNames = splitCSVLine(line, delimiter);
        if (has_rowNames) {
            colNames.value().erase(colNames.value().begin());
        }
    }

    // Read data rows
    int row_count = 0;
    while ((max_rows < 0 || row_count < max_rows) && std::getline(instream, line)) {
        auto row = splitCSVLine(line, delimiter);

        if (row.empty()) continue;  // Skip empty lines

        // Handle row names if present
        if (has_rowNames) {
            if (!rowNames) {
                rowNames.emplace();
            }
            rowNames->push_back(row.front());
            row.erase(row.begin());
        }
        data.push_back(row);
        row_count++;
    }
    return Table(data, colNames, rowNames, tolerateMissingLastColumn);
}
std::tuple<std::vector<std::string>, std::vector<std::string>,
    std::vector<int> > IoUtils::read_noe_table(std::istream &instream, bool has_header) {
    std::vector<std::tuple<std::string, std::string, int> > temp{};

    std::string line;
    std::regex lineTemplate{"(.+)\t(.+)\t(.+)"};
    std::smatch sm;
    bool header_consumed = false;
    while (std::getline(instream, line)) {
        std::regex_match(line, sm, lineTemplate);
        std::string g1 = sm[1];
        std::string g2 = sm[2];
        std::string val_str = sm[3];
        if (has_header && !header_consumed) {
            //TODO use consumed header
            //			std::cout << g1 << '\t' << g2 << '\t' << val_str << std::endl;
            header_consumed = true;
        } else {
            g1 = strip_enclosing_quotes(g1);
            g2 = strip_enclosing_quotes(g2);
            int val_int = std::stoi(val_str);
            //			std::cout << g1 << '\t' << g2 << '\t' << val_int << std::endl;
            temp.emplace_back(g1, g2, val_int);
        }
    }
    std::vector<std::string> group1{}, group2{};
    std::vector<int> values{};
    for (const auto &[g1, g2, v]: temp) {
        group1.emplace_back(g1);
        group2.emplace_back(g2);
        values.emplace_back(v);
    }
    //TODO change this tuple of vectors into a vector of tuples
    std::tuple<std::vector<std::string>, std::vector<std::string>,
        std::vector<int> > ret{group1, group2, values};
    return ret;
}


std::vector<std::string>
IoUtils::split(const std::string &str, const std::string &delim) {
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


std::map<std::string, std::vector<std::string> >
IoUtils::read_noe_groups(std::istream &instream) {
    std::map<std::string, std::vector<std::string> > ret{};
    std::string line;
    std::regex lineTemplate{R"(^\s*(.*?)\s*=\s*(.*?)(\s*,?\s*)$)"};
    std::smatch sm;

    while (std::getline(instream, line)) {
        if (!std::regex_match(line, sm, lineTemplate))
            continue;
        std::string key = strip_enclosing_quotes(sm[1], '`');
        std::string temp_val = sm[2];
        std::vector<std::string> val;
        if (temp_val.rfind("c(", 0) == 0) {
            // str.startwith() is available in C++20 at least
            auto tokens = split(temp_val.substr(2, temp_val.length() - 3), ",\\s*");
            for (auto i = 0; i < tokens.size(); ++i) { // NOLINT(modernize-loop-convert)
                tokens[i] = strip_enclosing_quotes(tokens[i], '\"');
            }
            val = tokens;
        } else {
            val = {strip_enclosing_quotes(temp_val, '\"')};
        }
        ret.insert({key, val});
    }
    return ret;
}


std::vector<int>
IoUtils::getGmxNdxGroup(const std::string &filename, const std::string &groupName) {
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

std::map<std::string, std::vector<int> >
IoUtils::getAllGmxNdxGroups(const std::string &filename) {
    std::map<std::string, std::vector<int> > ret;
    std::ifstream indexFile(filename);
    if (!indexFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return ret;
    }
    std::string line;
    std::vector<int> indices; //FIXME is it safe to declare it this way?
    while (std::getline(indexFile, line)) {
        if (line.empty()) continue;
        if (line[0] == '[') {
            auto closing = line.find(']');
            if (closing != std::string::npos && closing > 0) {
                std::string groupName;
                groupName = line.substr(2, closing - 2);
                // TODO trim the whitespace better than hard coding (2 , closing -2)
                ret[groupName] = indices;
                //        	        std::cout << "--------------------{{" << groupName << "}}" << std::endl;
                continue;
            } else {
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

bool IoUtils::isNotPrepared(const std::string &atomName) {
    return std::regex_search(atomName, IoUtils::UNPREPARED_NAMES_MASK);
}

bool IoUtils::should_handleNames(const std::map<std::string, int> &atomNameMapping) {
    for (const auto &[atomName, id]: atomNameMapping)
        if (isNotPrepared(atomName)) {
            std::cerr << "WARNING: It seems that your data is from an unprepared file. "
                    "We will try to handle it, but we can not guarantee the results." << std::endl;
            return true;
        }
    return false;
}

/**
 * IMPORTANT: Notice that we assume for this function to lower name ranks correctly that the function is called
 * sequentially with atom names in lexical order (e.g. "HB2 MET" is called before "HB3 MET" and not vise versa).
 */
std::string IoUtils::normalizeName(std::string atomId, const bool lowerNameRanks) {
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

std::regex atomRecordTemplate{"^((ATOM  )|(HETATM))([0-9 ]{5}) (.{15})   ([0-9 .-]{8})([0-9 .-]{8})([0-9 .-]{8}).+$"};
std::regex modelRecordTemplate{"^MODEL \\s*([0-9]+).*$"};

template<typename retMapKey, typename retMapValue>
std::map<retMapKey, retMapValue>
IoUtils::getAtomMappingFromPdb(const std::string &pdbFilename,
                               const std::function<void(std::map<retMapKey, retMapValue> &ret, const std::smatch &sm)> &
                               mappingFunc) {
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
        if (std::regex_match(line, sm, atomRecordTemplate)) {
            //fill_atomId_to_index_Map(ret, sm);
            mappingFunc(ret, sm);
        }
        else if (std::regex_match(line, sm, modelRecordTemplate)) {
            int modelNo = static_cast<int>(strtol(sm[1].str().c_str(), nullptr, 10));
            if (modelNo ==1) {
                //NOTHING
            }
        } else if (line.find("ENDMDL") == 0) {
            break; //TODO Handle other cases of ENDMDL
        }
    }
    pdbFileStream.close();
    //    std::cerr << "number of items in the map = " << ret.size() << std::endl;
    return ret; //std::move(&ret);
}

template<typename retMapKey, typename retMapValue>
std::vector<std::map<retMapKey, retMapValue> >
IoUtils::getMultipleAtomMappingFromPdb(const std::string &pdbFilename,
                                       const std::function<void(std::map<retMapKey, retMapValue> &ret,
                                                                const std::smatch &sm)> &mappingFunc) {
    std::vector<std::map<retMapKey, retMapValue> > ret_vector = {};

    std::ifstream pdbFileStream(pdbFilename);
    if (pdbFileStream.fail()) {
        std::cerr << "ERROR: Could not open " << pdbFilename << "\n";
        throw std::runtime_error("File not found: " + pdbFilename);
    }
    std::string line;
    std::smatch sm;
    bool addNextModel = true;
    std::map<retMapKey, retMapValue> modelMap{};

    while (std::getline(pdbFileStream, line)) {
        if (line.empty()) continue;
        if (std::regex_match(line, sm, atomRecordTemplate)) {
            //fill_atomId_to_index_Map(ret, sm);
            mappingFunc(modelMap, sm);
        } else if (std::regex_match(line, sm, modelRecordTemplate)) {
            addNextModel = true;
        } else if (line.find("ENDMDL") == 0) {
            ret_vector.emplace_back(std::move(modelMap));
            modelMap = {};
            addNextModel = false;
        }
    }
    pdbFileStream.close();
    if (addNextModel) {
        ret_vector.emplace_back(std::move(modelMap));
    }
    return ret_vector; //std::move(&ret);
}

void IoUtils::fill_atomId_to_index_Map(std::map<std::string, int> &ret, const std::smatch &sm) {
    int atomIndex;
    std::string atomId = sm[5];
    normalizeName(atomId, false);
    std::istringstream iss(sm[4]);
    iss >> atomIndex;
    ret[atomId] = atomIndex;
}

template<typename KEnRef_Real>
void IoUtils::fill_atomIndex1_to_coords_Map(std::map<int, Eigen::RowVector3<KEnRef_Real> > &ret, const std::smatch &sm) {
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
    if (paramsFileStream.fail()) {
        std::cerr << "ERROR: Could not open " << fileName << "\n";
        throw std::runtime_error("File not found: " + fileName);
    }
    return readParams(paramsFileStream);
}

std::map<std::string, std::string> IoUtils::readParams(std::istream &paramsFileStream) {
    std::string line;
    const std::regex recordTemplate{R"(^\s*(.+?)\s*=\s*(\S.*?)\s*(#.*)?)"};
    std::smatch sm;
    std::map<std::string, std::string> ret{};
    while (std::getline(paramsFileStream, line)) {
        if (line.empty()) continue;
        if (std::regex_match(line, sm, recordTemplate)) {
            const std::string &key = sm[1].str();
            if (key[0] == '#') continue;
            const std::string &value = strip_enclosing_quotes(sm[2].str());
            ret[key] = value;
        }
    }
    return ret;
}

Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>
IoUtils::subset_idx_to_grouping_mat(const std::vector<std::vector<std::vector<int>>>& multipleGroupings) {
    int maxRowLength = -1;
    for (const auto& grouping : multipleGroupings) {
        int rowLenSum = 0;
        for (const auto& group: grouping) {
            rowLenSum += static_cast<int>(group.size());
        }
        if (maxRowLength < rowLenSum) {
            maxRowLength = rowLenSum;
        }
    }
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> ret(multipleGroupings.size(), maxRowLength);

    for (int i = 0; i < multipleGroupings.size(); ++i) {
        for (int j = 0; j < multipleGroupings[i].size(); ++j) {
            auto subsetWithSameIndex = multipleGroupings[i][j];
            for (int index: subsetWithSameIndex) {
                ret(i, index) = j;
            }
        }
    }
    return ret;
}

std::vector<std::vector<std::vector<int>>>
IoUtils::grouping_mat_to_subset_idx(const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& grouping_matrix) {
    auto rows = grouping_matrix.rows();
    std::vector<std::vector<std::vector<int>>> ret(rows);
    for (int row = 0; row < rows; ++row) {
        Eigen::RowVectorX<int> grouping_row = grouping_matrix.row(row);
        // std::cout << "grouping_row " << grouping << std::endl;
        ret.at(row) = std::move(grouping_mat_row_to_subset_idx(grouping_row));
    }
    return ret;
}

std::vector<std::vector<int>>
IoUtils::grouping_mat_row_to_subset_idx(const Eigen::RowVectorX<int>& grouping_matrix) {
    std::map<int, std::vector<int>> group_map;
    // Iterate through the grouping vector (0-based or 1-based)
    for (int i = 0; i < grouping_matrix.size(); ++i) {
        int group = grouping_matrix(i);
        group_map[group].push_back(i);  // i+1 for 1-based indexing (like R)
    }

    // Convert the map into a vector of vectors
    std::vector<std::vector<int>> ret{};
    for (const auto& pair : group_map) {
        ret.push_back(pair.second);
    }
    return ret;
}

int IoUtils::getEnvParam(const std::string &paramName, int defaultValue) {
    return static_cast<int>(getEnvParam(paramName, static_cast<long>(defaultValue)));
}

long IoUtils::getEnvParam(const std::string &paramName, const long defaultValue) {
    if (const char *pSEnvParam = std::getenv(paramName.c_str())) {
        long retValue = std::strtol(pSEnvParam, nullptr, 10);
        std::cout << paramName << " is: " << retValue << '\n';
        return retValue;
    } else {
        std::cout << "No " << paramName << " identified. Will use default value of " << defaultValue << std::endl;
        return defaultValue;
    }
}

std::string IoUtils::getEnvParam(const char *paramName, const char *defaultValue) {
    return getEnvParam(std::string(paramName), std::string(defaultValue));
}
std::string IoUtils::getEnvParam(const std::string &paramName, const char *defaultValue) {
    return getEnvParam(paramName, std::string(defaultValue));
}
std::string IoUtils::getEnvParam(const std::string &paramName, const std::string &defaultValue) {
    std::string retValue = defaultValue;
    if (const char *pSEnvParam = std::getenv(paramName.c_str())) {
        retValue = pSEnvParam;
        std::cout << paramName << " is: " << retValue << '\n';
    } else {
        std::cout << "No " << paramName << " identified. Will use default value of " << defaultValue << std::endl;
    }
    return retValue;
}

template<typename KEnRef_Real>
KEnRef_Real IoUtils::getEnvParam(const std::string &paramName, KEnRef_Real defaultValue) {
    if (const char *pSEnvParam = std::getenv(paramName.c_str())) {
        std::stringstream sStream(pSEnvParam);
        KEnRef_Real retValue;
        sStream >> retValue;
        std::cout << paramName << " is: " << retValue << '\n';
        return retValue;
    } else {
        std::cout << "No " << paramName << " identified. Will use default value of " << defaultValue << std::endl;
        return defaultValue;
    }
}


std::vector<std::tuple<int, int> > IoUtils::readAtomIdPairs(const std::string &fileName) {
    std::ifstream atomIdPairsFileStream(fileName);
    auto tempAtomIdPairsTable = IoUtils::read_uniform_table_of<int>(atomIdPairsFileStream);
    std::vector<std::tuple<int, int> > atomIdPairs;
    for (auto row: tempAtomIdPairsTable) {
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

template<typename KEnRef_Real>
CoordsMatrixType<KEnRef_Real>
IoUtils::extractCoords(const std::vector<int> &atomIndices, bool indicesOneBased,
                       std::map<int, Eigen::RowVector3<KEnRef_Real> > &allAtomCoords, bool mapOneBased) {
    auto coords = CoordsMatrixType<KEnRef_Real>(atomIndices.size(), 3);
    int delta = 0;
    if (indicesOneBased ^ mapOneBased)
        indicesOneBased ? delta-- : delta++;


    for (int i = 0; i < atomIndices.size(); ++i) {
        int key = atomIndices[i] + delta;
        coords(i, Eigen::indexing::all) = allAtomCoords[key];
    }
    return coords;
}

template<typename TYPE>
void IoUtils::printVector(const std::vector<TYPE> &vec) {
    for (const auto &val: vec)
        std::cout << val << " ";
    std::cout << std::endl;
}

inline const std::string &IoUtils::aa3(const std::string &aa1) {
    return aa1_to_aa3.at(aa1);
}

inline const std::string &IoUtils::aa1(const std::string &aa3) {
    return aa3_to_aa1.at(aa3);
}

// Type conversion helper
template<typename Scalar>
Scalar IoUtils::convertValue(const std::string& str) {
    if constexpr (std::is_same_v<Scalar, std::string>) {
        return str;
    }
    else if constexpr (std::is_integral_v<Scalar>) {
        try {
            size_t pos = 0;
            long long val = std::stoll(str, &pos);
            if (pos != str.length()) {
                throw std::invalid_argument("Invalid characters in numeric conversion");
            }
            return static_cast<Scalar>(val);
        } catch (const std::exception& e) {
            throw std::invalid_argument("Failed to convert '" + str + "' to integer: " + e.what());
        }
    }
    else if constexpr (std::is_floating_point_v<Scalar>) {
        try {
            size_t pos = 0;
            double val = std::stod(str, &pos);
            if (pos != str.length()) {
                throw std::invalid_argument("Invalid characters in numeric conversion");
            }
            return static_cast<Scalar>(val);
        } catch (const std::exception& e) {
            throw std::invalid_argument("Failed to convert '" + str + "' to floating point: " + e.what());
        }
    }
    else {
        // static_assert(sizeof(Scalar) == 0, "Unsupported scalar type for conversion");
        throw std::invalid_argument("Unsupported scalar type for conversion");
    }
}

// Safe conversion alternative (no-throw version)
template<typename Scalar>
bool IoUtils::tryConvertValue(const std::string& str, Scalar& out) noexcept {
    try {
        out = convertValue<Scalar>(str);
        return true;
    } catch (...) {
        return false;
    }
}

///////// Template Declarations /////////////////////////

template std::vector<std::vector<int> > IoUtils::read_uniform_table_of(std::istream &);
template std::vector<std::vector<double> > IoUtils::read_uniform_table_of(std::istream &);

template float IoUtils::getEnvParam(const std::string &paramName, float defaultValue);
template double IoUtils::getEnvParam(const std::string &paramName, double defaultValue);

template CoordsMatrixType<float> IoUtils::extractCoords(const std::vector<int> &atomIndices, bool indicesOneBased, std::map<int, Eigen::RowVector3<float> > &allAtomCoords, bool mapOneBased);
template CoordsMatrixType<double> IoUtils::extractCoords(const std::vector<int> &atomIndices, bool indicesOneBased, std::map<int, Eigen::RowVector3<double> > &allAtomCoords, bool mapOneBased);

template void IoUtils::printVector(const std::vector<bool> &vec);
template void IoUtils::printVector(const std::vector<int> &vec);
template void IoUtils::printVector(const std::vector<float> &vec);
template void IoUtils::printVector(const std::vector<double> &vec);
template void IoUtils::printVector(const std::vector<std::string> &vec);

template void IoUtils::fill_atomIndex1_to_coords_Map<float>(std::map<int, Eigen::RowVector3<float> > &ret, const std::smatch &sm);
template void IoUtils::fill_atomIndex1_to_coords_Map<double>(std::map<int, Eigen::RowVector3<double> > &ret, const std::smatch &sm);
template std::map<std::string, int> IoUtils::getAtomMappingFromPdb(const std::string &pdbFilename, const std::function<void( std::map<std::string, int> &, const std::smatch &)> &mappingFunc);
template std::map<int, Eigen::RowVector3<float> > IoUtils::getAtomMappingFromPdb<int, Eigen::RowVector3<float> >( const std::string &pdbFilename, const std::function<void(std::map<int, Eigen::RowVector3<float> > &, const std::smatch &)> &mappingFunc);
template std::map<int, Eigen::RowVector3<double> > IoUtils::getAtomMappingFromPdb<int, Eigen::RowVector3<double> >( const std::string &pdbFilename, const std::function<void(std::map<int, Eigen::RowVector3<double> > &, const std::smatch &)> &mappingFunc);
template std::vector<std::map<int, Eigen::RowVector3<float>>>  IoUtils::getMultipleAtomMappingFromPdb<int, Eigen::RowVector3<float> > ( const std::string &pdbFilename, const std::function<void(std::map<int, Eigen::RowVector3<float> >  &, const std::smatch &)> &mappingFunc);
template std::vector<std::map<int, Eigen::RowVector3<double>>> IoUtils::getMultipleAtomMappingFromPdb<int, Eigen::RowVector3<double> >( const std::string &pdbFilename, const std::function<void(std::map<int, Eigen::RowVector3<double> > &, const std::smatch &)> &mappingFunc);

template int IoUtils::convertValue<int>(const std::string& str);
template long IoUtils::convertValue<long>(const std::string& str);
template float IoUtils::convertValue<float>(const std::string& str);
template double IoUtils::convertValue<double>(const std::string& str);
template bool IoUtils::tryConvertValue<int>(const std::string& str, int& out) noexcept;
template bool IoUtils::tryConvertValue<long>(const std::string& str, long& out) noexcept;
template bool IoUtils::tryConvertValue<float>(const std::string& str, float& out) noexcept;
template bool IoUtils::tryConvertValue<double>(const std::string& str, double& out) noexcept;