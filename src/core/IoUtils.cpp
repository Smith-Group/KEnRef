#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "core/IoUtils.h"

std::regex atomRecordTemplate{"^((ATOM  )|(HETATM))([0-9 ]{5}) (.{15})   ([0-9 .-]{8})([0-9 .-]{8})([0-9 .-]{8}).+$"};
std::regex modelRecordTemplate{"^MODEL \\s*([0-9]+).*$"};

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
    const std::map<std::string, int> &atomNames_2_atomIds, int numOmpThreads){

    auto atomId_pairs = std::vector<std::tuple<int, int>>(atomName_pairs.size());
    // Fill the vector using atomNames_2_atomIds (map lookups are read-only => safe to parallelize)
#pragma omp parallel for num_threads(numOmpThreads)
    for (int i = 0; i < atomName_pairs.size(); ++i) {
        auto [left, right] = atomName_pairs.at(i);
        // I use at() instead of operator[] to force an exception to be thrown
        atomId_pairs.at(i) = std::move(std::tuple{atomNames_2_atomIds.at(left), atomNames_2_atomIds.at(right)});
    }
    return atomId_pairs;
}

std::vector<std::string>
IoUtils::find_matching_files(const std::string &folder_path, const std::string &pattern_string) {
    std::vector<std::string> matching_files;
    const std::regex pattern(pattern_string);
    try {
        for (const auto& entry : std::filesystem::directory_iterator(folder_path)) {
            if (entry.is_regular_file()) {
                std::string filename = entry.path().filename().string();
                if (std::regex_match(filename, pattern)) {
                    matching_files.push_back(entry.path().filename().string());
                }
            }
        }
    } catch (const std::filesystem::filesystem_error& ex) {
        std::cerr << "Filesystem error: " << ex.what() << std::endl;
    }
    return matching_files;
}

// Function that also extracts the captured group
std::vector<std::string>
IoUtils::find_spec_den_data_prefixes(const std::string& folder_path) {
    std::vector<std::string> results;
    const std::regex pattern(R"((\d+-\d+)_atom_pairs\.csv)");

    try {
        for (const auto& entry : std::filesystem::directory_iterator(folder_path)) {
            if (entry.is_regular_file()) {
                std::string filename = entry.path().filename().string();
                std::smatch matches;

                if (std::regex_match(filename, matches, pattern) /*&& matches.size() >= 1*/) {
                    results.emplace_back(matches[1].str());
                }
            }
        }
    } catch (const std::filesystem::filesystem_error& ex) {
        std::cerr << "Filesystem error: " << ex.what() << std::endl;
    }
    return results;
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
    const std::vector<std::string> &spec_den_data_prefixes = find_spec_den_data_prefixes(experimentalDataFolder);
    std::vector<SpecDenData<KEnRef_Real_t>> spec_den_data_vector;
    spec_den_data_vector.reserve(spec_den_data_prefixes.size());
    for (const auto & spec_den_data_prefix : spec_den_data_prefixes) {
        //AtomPairs
        std::string atomPairAndSigmaFileName = std::filesystem::path(experimentalDataFolder) / (spec_den_data_prefix + "_atom_pairs.csv");
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
            if (! atomPairAndSigmaTable.isRowComplete(row) || atomPairAndSigmaTable.at(row, "sigma").empty())
                break;
            const auto &valueStr = atomPairAndSigmaTable.at(row, "sigma");
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
        std::string multiple_grouping_fileName = std::filesystem::path(experimentalDataFolder) / (spec_den_data_prefix +"_groupings.csv");
        const auto &grouping_matrix = NamedMatrix<int>(IoUtils::readTable(multiple_grouping_fileName, false, false).toNamedMatrix<int>().array() - 1);
        const auto &multiple_grouping = IoUtils::grouping_mat_to_subset_idx(grouping_matrix);
        //a_coef
        std::string aCoefFileName = std::filesystem::path(experimentalDataFolder) / ( spec_den_data_prefix +"_a_coef.csv");
        const auto &a_coef = IoUtils::readTable(aCoefFileName, true,false, "\\s*,\\s*", -1, false).toNamedMatrix<KEnRef_Real_t>();
        //lambda_coef
        std::string lambdaCoefFileName = std::filesystem::path(experimentalDataFolder) / (spec_den_data_prefix +"_lambda_coef.csv");
        const auto &lambda_coef = IoUtils::readTable(lambdaCoefFileName, true,true, "\\s*,\\s*", -1, false).toNamedMatrix<KEnRef_Real_t>();

        // spec_den_data_vector.emplace_back(atomPairs, sigma, multiple_grouping, a_coef, lambda_coef);
        spec_den_data_vector.emplace_back(std::move(SpecDenData<KEnRef_Real_t>{atomPairs, sigma, multiple_grouping, a_coef, lambda_coef}));
    }
    return spec_den_data_vector;
}

namespace {
// Read the first line of a file verbatim (sans line terminator) so the *_atom_relax.csv
// header can be kept for byte-exact debugging / backtracking.
std::string readHeaderLine(const std::string &fileName) {
    std::ifstream in(fileName);
    if (!in.is_open()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        throw std::runtime_error(std::string("Can't open file:").append(fileName));
    }
    std::string line;
    std::getline(in, line);
    if (!line.empty() && line.back() == '\r') line.pop_back();   // drop CR from CRLF endings
    return line;
}

// A column carries relaxation data (vs. an atom-pair / metadata column) iff its name
// ends in one of these suffixes. Mirrors `!grepl("_(value|k|coef|freq)$", ...)` in ke.R.
bool is_relax_col(const std::string &c) {
    return IoUtils::ends_with(c, "_value") || IoUtils::ends_with(c, "_k")
        || IoUtils::ends_with(c, "_coef") || IoUtils::ends_with(c, "_freq");
}

size_t index_of(const std::vector<std::string> &v, const std::string &s) {
    return static_cast<size_t>(std::find(v.begin(), v.end(), s) - v.begin());
}

// Mirror of atom_relax_df_to_spec_den_term_array_list() + atom_relax_columns_to_spec_den_term_array():
// reconstruct one RelaxEntry per rate from the first `n_data_rows` rows of the table.
std::vector<RelaxEntry<KEnRef_Real_t> >
build_relax_data_list(const Table &table, const size_t n_data_rows) {
    static const std::regex term_full_re(R"(^(.+)_([^_]+)_(coef|freq)$)");   // <rate>_<freq>_<coef|freq>
    static const std::regex term_suffix_re(R"(^([^_]+)_(coef|freq)$)");      // <freq>_<coef|freq>
    const auto &colNames = table.getColNames();

    // 1) discover rate names, in first-appearance (column) order. Parsing is right-anchored,
    //    so rate names may themselves contain underscores (e.g. "r1_400").
    std::vector<std::string> rate_names;
    auto addRate = [&](const std::string &r) {
        if (std::find(rate_names.begin(), rate_names.end(), r) == rate_names.end())
            rate_names.push_back(r);
    };
    for (const auto &col : colNames) {
        if (IoUtils::ends_with(col, "_value")) { addRate(col.substr(0, col.size() - 6)); continue; }
        if (IoUtils::ends_with(col, "_k"))     { addRate(col.substr(0, col.size() - 2)); continue; }
        std::smatch m;
        if (std::regex_match(col, m, term_full_re)) addRate(m[1].str());
        // otherwise: atom-pair / metadata column -> ignored
    }

    std::vector<RelaxEntry<KEnRef_Real_t> > out;
    out.reserve(rate_names.size());
    const auto N = static_cast<Eigen::Index>(n_data_rows);

    for (const auto &rate_name : rate_names) {
        const std::string value_col = rate_name + "_value";
        const std::string k_col = rate_name + "_k";
        const std::string rate_prefix = rate_name + "_";

        // collect this rate's term columns; parse <freq>,<component>; gather unique freq labels in order
        std::vector<std::string> term_cols, freq_labels, components, freq_levels;
        for (const auto &col : colNames) {
            if (col == value_col || col == k_col) continue;
            if (col.size() <= rate_prefix.size() || col.compare(0, rate_prefix.size(), rate_prefix) != 0) continue;
            const std::string suffix = col.substr(rate_prefix.size());
            std::smatch m;
            if (!std::regex_match(suffix, m, term_suffix_re))
                throw std::runtime_error("All term columns for rate `" + rate_name + "` must match `" +
                    rate_name + "_<freq_name>_<coef|freq>` with no underscores in <freq_name>");
            term_cols.push_back(col);
            freq_labels.push_back(m[1].str());
            components.push_back(m[2].str());
            if (std::find(freq_levels.begin(), freq_levels.end(), m[1].str()) == freq_levels.end())
                freq_levels.push_back(m[1].str());
        }
        if (term_cols.empty())
            throw std::runtime_error("No spectral-density term columns found for rate `" + rate_name + "`");

        // each freq term must carry both a coef and a freq column
        for (const auto &fl : freq_levels) {
            bool hasCoef = false, hasFreq = false;
            for (size_t j = 0; j < freq_labels.size(); ++j)
                if (freq_labels[j] == fl) (components[j] == "coef" ? hasCoef : hasFreq) = true;
            if (!(hasCoef && hasFreq))
                throw std::runtime_error("Rate `" + rate_name + "` term `" + fl + "` must have both `coef` and `freq` columns");
        }

        const auto T = static_cast<Eigen::Index>(freq_levels.size());
        NamedMatrix<KEnRef_Real_t> coef(N, T);
        NamedMatrix<KEnRef_Real_t> freq(N, T);
        for (size_t j = 0; j < term_cols.size(); ++j) {
            const size_t c = index_of(freq_levels, freq_labels[j]);
            const size_t colIdx = table.getColIndex(term_cols[j]);
            NamedMatrix<KEnRef_Real_t> &target = (components[j] == "coef") ? coef : freq;
            for (size_t r = 0; r < n_data_rows; ++r)
                target(r, c) = IoUtils::convertValue<KEnRef_Real_t>(table.at(r, colIdx));
        }
        coef.setColNames(freq_levels);
        freq.setColNames(freq_levels);

        // optional target values / weights (first N rows)
        std::optional<NamedVector<KEnRef_Real_t> > value;
        if (std::find(colNames.begin(), colNames.end(), value_col) != colNames.end()) {
            const size_t colIdx = table.getColIndex(value_col);
            NamedVector<KEnRef_Real_t> v(N);
            for (size_t r = 0; r < n_data_rows; ++r) v(r, 0) = IoUtils::convertValue<KEnRef_Real_t>(table.at(r, colIdx));
            value = std::move(v);
        }
        std::optional<NamedVector<KEnRef_Real_t> > k;
        if (std::find(colNames.begin(), colNames.end(), k_col) != colNames.end()) {
            const size_t colIdx = table.getColIndex(k_col);
            NamedVector<KEnRef_Real_t> kk(N);
            for (size_t r = 0; r < n_data_rows; ++r) kk(r, 0) = IoUtils::convertValue<KEnRef_Real_t>(table.at(r, colIdx));
            k = std::move(kk);
        }

        out.push_back(RelaxEntry<KEnRef_Real_t>{rate_name, std::move(value), std::move(k),
            SpecDenTermArray<KEnRef_Real_t>{std::move(coef), std::move(freq)}});
    }
    return out;
}
}   // anonymous namespace

std::vector<std::string>
IoUtils::find_spec_den_relax_data_prefixes(const std::string &folder_path) {
    std::vector<std::string> results;
    // Generalized prefix: capture anything preceding the `_atom_relax.csv` suffix
    // (e.g. "3-5", "gb3_15n"), not only the `\d+-\d+` interaction convention.
    const std::regex pattern(R"((.+)_atom_relax\.csv)");
    try {
        for (const auto &entry : std::filesystem::directory_iterator(folder_path)) {
            if (entry.is_regular_file()) {
                std::string filename = entry.path().filename().string();
                std::smatch matches;
                if (std::regex_match(filename, matches, pattern)) {
                    results.emplace_back(matches[1].str());
                }
            }
        }
    } catch (const std::filesystem::filesystem_error &ex) {
        std::cerr << "Filesystem error: " << ex.what() << std::endl;
    }
    return results;
}

std::vector<SpecDenRelaxData<KEnRef_Real_t> >
IoUtils::load_spec_den_relax_data(const std::string &experimentalDataFolder, const bool handleNames) {
    const std::vector<std::string> &prefixes = find_spec_den_relax_data_prefixes(experimentalDataFolder);
    std::vector<SpecDenRelaxData<KEnRef_Real_t> > result;
    result.reserve(prefixes.size());
    for (const auto &prefix : prefixes) {
        // ---- *_atom_relax.csv ----
        std::string atomRelaxFileName = std::filesystem::path(experimentalDataFolder) / (prefix + "_atom_relax.csv");
        std::cout << atomRelaxFileName << std::endl;
        std::string headerLine = readHeaderLine(atomRelaxFileName);
        const Table &atomRelaxTable = IoUtils::readTable(atomRelaxFileName, true, false, "\\s*,\\s*", -1, true);
        const auto &colNames = atomRelaxTable.getColNames();
        const size_t numRows = atomRelaxTable.rowCount();
        if (colNames.size() < 2)
            throw std::runtime_error("`" + atomRelaxFileName + "` must have at least two atom-pair columns");

        // unit flag: the first two column headers must both end in `_unit`, or neither may.
        // This is independent of `handleNames` (the identifiers are still atom names).
        const bool col0Unit = ends_with(colNames[0], "_unit");
        const bool col1Unit = ends_with(colNames[1], "_unit");
        if (col0Unit != col1Unit)
            throw std::runtime_error("The first two columns of `*_atom_relax.csv` must either both end in `_unit` or neither may");
        const bool unit = col0Unit && col1Unit;

        // relaxation columns (everything that is not an atom-pair / metadata column)
        std::vector<size_t> relaxColIdx;
        for (size_t c = 0; c < colNames.size(); ++c)
            if (is_relax_col(colNames[c])) relaxColIdx.push_back(c);

        // n_data_rows N: leading rows carrying a full set of (non-empty) relaxation columns. Trailing
        // rows hold atom pairs only and alias back into the first N via blockRow (idx % N), as sigma does.
        size_t n_data_rows = 0;
        if (relaxColIdx.empty()) {
            n_data_rows = numRows;   // no rates -> every row is just an atom pair
        } else {
            for (size_t r = 0; r < numRows; ++r) {
                if (!atomRelaxTable.isRowComplete(r)) break;   // guards against touching short rows' cells
                bool anyEmpty = false;
                for (const size_t c : relaxColIdx)
                    if (atomRelaxTable.at(r, c).empty()) { anyEmpty = true; break; }
                if (anyEmpty) break;
                ++n_data_rows;
            }
        }

        // atom pairs: the first two columns are the identifiers (for all rows), normalized as in the
        // sigma loader. `unit` does not change this — handleNames still applies.
        std::vector<std::tuple<std::string, std::string> > atomPairs(numRows);
        for (size_t r = 0; r < numRows; ++r) {
            auto a1 = IoUtils::normalizeName(atomRelaxTable.at(r, static_cast<size_t>(0)), handleNames);
            auto a2 = IoUtils::normalizeName(atomRelaxTable.at(r, static_cast<size_t>(1)), handleNames);
            atomPairs[r] = std::tuple<std::string, std::string>{a1, a2};
        }

        // relax_data_list from the first N rows
        std::vector<RelaxEntry<KEnRef_Real_t> > relax_data_list = build_relax_data_list(atomRelaxTable, n_data_rows);

        // ---- *_groupings.csv ----
        std::string groupingsFileName = std::filesystem::path(experimentalDataFolder) / (prefix + "_groupings.csv");
        const auto &grouping_matrix = NamedMatrix<int>(IoUtils::readTable(groupingsFileName, false, false).toNamedMatrix<int>().array() - 1);
        const auto &multiple_grouping = IoUtils::grouping_mat_to_subset_idx(grouping_matrix);
        // ---- *_a_coef.csv  (a_int_coef) ----
        std::string aCoefFileName = std::filesystem::path(experimentalDataFolder) / (prefix + "_a_coef.csv");
        const auto &a_int_coef = IoUtils::readTable(aCoefFileName, true, false, "\\s*,\\s*", -1, false).toNamedMatrix<KEnRef_Real_t>();
        // ---- *_lambda_coef.csv  (lambda_int_coef) ----
        std::string lambdaCoefFileName = std::filesystem::path(experimentalDataFolder) / (prefix + "_lambda_coef.csv");
        const auto &lambda_int_coef = IoUtils::readTable(lambdaCoefFileName, true, true, "\\s*,\\s*", -1, false).toNamedMatrix<KEnRef_Real_t>();

        result.emplace_back(SpecDenRelaxData<KEnRef_Real_t>{
            atomPairs, unit, std::move(relax_data_list), n_data_rows, std::move(headerLine),
            multiple_grouping, a_int_coef, lambda_int_coef});
    }
    return result;
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

std::vector<std::string>
IoUtils::split(const std::string &str, const std::string &delim) {
    std::regex delim_re(delim);
    std::vector<std::string> ret;
    std::sregex_token_iterator iter(str.begin(), str.end(), delim_re, -1);
    std::sregex_token_iterator end;

    while (iter != end) {
        ret.emplace_back(iter->str());
        ++iter;
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
    std::vector<int> indices;
    while (std::getline(indexFile, line)) {
        if (line.empty()) continue;
        if (line[0] == '[') {
            if (auto closing = line.find(']'); closing != std::string::npos && closing > 0) {
                std::string groupName;
                groupName = line.substr(2, closing - 2);
                // TODO trim the whitespace better than hard coding (2 , closing -2)
                ret[groupName] = indices;
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

bool IoUtils::should_handleNames(const Table &table) {
    bool handleUnpreparedAtomNames = false;
    for (int row = 0; row < table.rowCount(); row++) {
        if (IoUtils::isNotPrepared(table(row, "atom1")) || IoUtils::isNotPrepared(table(row, "atom2"))) {
            std::cerr << "WARNING: It seems that your data is from an unprepared file. We will try "
                    "to handle it, but we can not guarantee the results." << std::endl;
            handleUnpreparedAtomNames = true;
            break;
        }
    }
    return handleUnpreparedAtomNames;
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
    const std::function<void(std::map<retMapKey, retMapValue> &ret, const std::smatch &sm)> &mappingFunc) {
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

int IoUtils::grouping_width(const std::vector<std::vector<std::vector<int>>>& multipleGroupings) {
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
    return maxRowLength;
}

Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>
IoUtils::subset_idx_to_grouping_mat(const std::vector<std::vector<std::vector<int>>>& multipleGroupings) {
    const int maxRowLength = grouping_width(multipleGroupings);
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

// template float IoUtils::getEnvParam(const std::string &paramName, float defaultValue);
// template double IoUtils::getEnvParam(const std::string &paramName, double defaultValue);

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