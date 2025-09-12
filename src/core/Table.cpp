#include <algorithm>
#include "core/Table.h"
#include "core/IoUtils.h"

Table::Table() = default;

Table::Table(const std::vector<std::vector<std::string> > &data,
             const std::optional<std::vector<std::string> > &colNames,
             const std::optional<std::vector<std::string> > &rowNames,
             bool tolerateMissingLastColumn)
    : data_(data), colNames_(colNames), rowNames_(rowNames) {
    validateDimensions(tolerateMissingLastColumn);
}

// Accessors
std::string &Table::at(size_t row, size_t col) {
    checkBounds(row, col);
    return data_.at(row).at(col);
}

[[nodiscard]] const std::string &Table::at(const size_t row, const size_t col) const {
    checkBounds(row, col);
    return data_.at(row).at(col);
}

std::string &Table::operator()(const size_t row, const size_t col) { return at(row, col); }
const std::string &Table::operator()(const size_t row, const size_t col) const { return at(row, col); }

// Named access (throws if names not available)
std::string &Table::at(const std::string &rowName, const std::string &colName) {
    return at(getRowIndex(rowName), getColIndex(colName));
}

[[nodiscard]] const std::string &Table::at(const std::string &rowName, const std::string &colName) const {
    return at(getRowIndex(rowName), getColIndex(colName));
}

// Dimension accessors
[[nodiscard]] size_t Table::rowCount() const { return data_.size(); }
[[nodiscard]] size_t Table::colCount() const { return rowCount() > 0 ? data_[0].size() : 0; }

// Name management
[[nodiscard]] bool Table::hasColNames() const { return colNames_.has_value(); }
[[nodiscard]] bool Table::hasRowNames() const { return rowNames_.has_value() && rowNames_->size() == rowCount(); }

[[nodiscard]] const std::vector<std::string> &Table::getColNames() const {
    if (!colNames_) throw std::runtime_error("Column names not set");
    return *colNames_;
}

[[nodiscard]] const std::vector<std::string> &Table::getRowNames() const {
    if (!hasRowNames()) throw std::runtime_error("Row names not properly set");
    return *rowNames_;
}

void Table::setColNames(const std::vector<std::string> &names) {
    if (names.size() != colCount()) {
        throw std::invalid_argument("Column names size doesn't match table dimensions");
    }
    colNames_ = names;
}

void Table::setRowNames(const std::vector<std::string> &names) {
    if (names.size() != rowCount()) {
        throw std::invalid_argument("Row names size doesn't match table dimensions");
    }
    rowNames_ = names;
}

[[nodiscard]] std::vector<size_t> Table::colIndices(const std::vector<std::string> &colNames) const {
    if (!colNames_) throw std::runtime_error("Column names not set");
    std::vector<size_t> ret(colNames.size());
    for (int i = 0; i < colNames.size(); ++i)
        ret.at(i) = getColIndex(colNames[i]);
    return ret;
}

[[nodiscard]] std::vector<size_t> Table::rowIndices(const std::vector<std::string> &rowNames) const {
    if (!hasRowNames()) throw std::runtime_error("Row names not properly set");
    std::vector<size_t> ret(rowNames.size());
    for (int i = 0; i < rowNames.size(); ++i)
        ret.at(i) = getRowIndex(rowNames[i]);
    return ret;
}

// Data modification
void Table::addRow(const std::vector<std::string> &row,
                   const std::optional<std::string> &rowName) {
    if (row.size() != colCount() && colCount() > 0) {
        throw std::invalid_argument("Row size doesn't match table column count");
    }

    data_.push_back(row);

    if (rowName) {
        if (!rowNames_) {
            // Initialize with empty names for previous rows
            rowNames_ = std::vector<std::string>(rowCount() - 1, "");
        }
        rowNames_->push_back(*rowName);
    } else if (rowNames_) {
        // Add empty name if other rows have names
        rowNames_->push_back("");
    }
}

void Table::addColumn(const std::vector<std::string> &col,
                      const std::optional<std::string> &colName) {
    if (col.size() != rowCount()) {
        throw std::invalid_argument("Column size doesn't match table row count");
    }

    for (size_t i = 0; i < rowCount(); ++i) {
        data_[i].push_back(col[i]);
    }

    if (colName) {
        if (!colNames_) colNames_.emplace();
        colNames_->push_back(*colName);
    }
}

// Printing
void Table::print(std::ostream &os) const {
    // Print column names if available
    if (colNames_) {
        if (hasRowNames()) os << ",";
        for (size_t i = 0; i < colNames_->size(); ++i) {
            if (i > 0) os << ",";
            os << colNames_->at(i);
        }
        os << "\n";
    }

    // Print each row
    for (size_t i = 0; i < rowCount(); ++i) {
        // Print row name if available
        if (hasRowNames()) {
            os << rowNames_->at(i) << ",";
        }

        // Print row data
        for (size_t j = 0; j < colCount(); ++j) {
            if (j > 0) os << ",";
            os << data_[i][j];
        }
        os << "\n";
    }
}

void Table::validateDimensions(bool tolerateMissingLastColumn) const {
    if (data_.empty()) return;

    const size_t cols = data_[0].size();
    for (const auto &row: data_) {
        if (row.size() != cols) {
            if (!tolerateMissingLastColumn && row.size() != cols- 1)
                throw std::invalid_argument("All rows must have the same number of columns");
        }
    }

    if (rowNames_ && rowNames_->size() != rowCount()) {
        throw std::invalid_argument("Row names size doesn't match row count");
    }

    if (colNames_ && colNames_->size() != colCount()) {
        throw std::invalid_argument("Column names size doesn't match column count");
    }
}

void Table::checkBounds(const size_t row, const size_t col) const {
    if (row >= rowCount() || col >= colCount()) {
        throw std::out_of_range("Table index out of range");
    }
}

[[nodiscard]] size_t Table::getRowIndex(const std::string &name) const {
    if (!hasRowNames()) throw std::runtime_error("Row names not available");
    const auto it = std::find(rowNames_->begin(), rowNames_->end(), name);
    if (it == rowNames_->end()) throw std::runtime_error("Row name not found");
    return std::distance(rowNames_->begin(), it);
}

[[nodiscard]] size_t Table::getColIndex(const std::string &name) const {
    if (!colNames_) throw std::runtime_error("Column names not available");
    const auto it = std::find(colNames_->begin(), colNames_->end(), name);
    if (it == colNames_->end()) throw std::runtime_error("Column name not found");
    return std::distance(colNames_->begin(), it);
}

// template<typename Scalar, int Options>
// [[nodiscard]] Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> Table::toEigenMatrix() const {
//     const size_t rows = rowCount();
//     const size_t cols = colCount();
//     Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> matrix(rows, cols);
//
//     for (size_t i = 0; i < rows; ++i) {
//         for (size_t j = 0; j < cols; ++j) {
//             matrix(i, j) = IoUtils::convertValue<Scalar>(data_[i][j]);
//         }
//     }
//     return matrix;
// }

// template<typename KEnRef_Real>
template<typename Scalar, int  Options>
NamedMatrix<Scalar, Options> Table::toNamedMatrix() const {
    const size_t rows = rowCount();
    const size_t cols = colCount();
    NamedMatrix<Scalar, Options> matrix(rows, cols);

    // Convert data
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            matrix(i,j) = IoUtils::convertValue<Scalar>(data_[i][j]);
        }
    }

    // Set names if available
    if (colNames_) {
        matrix.setColNames(*colNames_);
    }
    if (rowNames_) {
        matrix.setRowNames(*rowNames_);
    }
    return matrix;
}
template<typename Scalar, int Options>
NamedVector<Scalar> Table::toNamedVector() const {
    if (colCount() != 1) {
        throw std::invalid_argument("Table must have exactly 1 column");
    }
    return NamedVector<Scalar, Options>(toNamedMatrix<Scalar>());
}

template<typename Scalar, int Options>
NamedRowVector<Scalar> Table::toNamedRowVector() const {
    if (rowCount() != 1) {
        throw std::invalid_argument("Table must have exactly 1 row");
    }
    return NamedRowVector<Scalar, Options>(toNamedMatrix<Scalar>());
}


// template Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> Table::toEigenMatrix() const;
// template Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Table::toEigenMatrix() const;
// template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Table::toEigenMatrix() const;
// template Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Table::toEigenMatrix() const;
// template Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Table::toEigenMatrix() const;
// template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Table::toEigenMatrix() const;
template NamedMatrix<int> Table::toNamedMatrix() const;
template NamedMatrix<float> Table::toNamedMatrix() const;
template NamedMatrix<double> Table::toNamedMatrix() const;
template NamedVector<int> Table::toNamedVector() const;
template NamedVector<float> Table::toNamedVector() const;
template NamedVector<double> Table::toNamedVector() const;
template NamedRowVector<int> Table::toNamedRowVector() const;
template NamedRowVector<float> Table::toNamedRowVector() const;
template NamedRowVector<double> Table::toNamedRowVector() const;
