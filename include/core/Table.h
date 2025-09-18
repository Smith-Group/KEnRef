#ifndef TABLE_H
#define TABLE_H

#include <vector>
#include <string>
#include <optional>
#include <iostream>
#include "core/NamedMatrix.h"


class Table {
public:
    // Constructors
    Table();

    explicit Table(const std::vector<std::vector<std::string>>& data,
                   const std::optional<std::vector<std::string>>& colNames = std::nullopt,
                   const std::optional<std::vector<std::string>>& rowNames = std::nullopt,
                   bool tolerateMissingLastColumn = false);

    // Accessors
    std::string& operator()(size_t row, size_t col);
    const std::string& operator()(size_t row, size_t col) const;
    std::string& operator()(size_t row, const std::string &colName);
    const std::string& operator()(size_t row, const std::string &colName) const;
    std::string& operator()(const std::string &rowName, size_t col);
    const std::string& operator()(const std::string &rowName, size_t col) const;
    std::string& at(size_t row, size_t col);
    [[nodiscard]] const std::string& at(size_t row, size_t col) const;
    std::string& at(const std::string& rowName, const std::string& colName);
    [[nodiscard]] const std::string& at(const std::string& rowName, const std::string& colName) const;
    std::string& at(size_t row, const std::string& colName);
    [[nodiscard]] const std::string& at(size_t row, const std::string& colName) const;
    std::string& at(const std::string& rowName, size_t col);
    [[nodiscard]] const std::string& at(const std::string& rowName, size_t col) const;

    // Dimension accessors
    [[nodiscard]] size_t rowCount() const;
    [[nodiscard]] size_t colCount() const;

    // Name management
    [[nodiscard]] bool hasColNames() const;
    [[nodiscard]] bool hasRowNames() const;
    [[nodiscard]] const std::vector<std::string>& getColNames() const;
    [[nodiscard]] const std::vector<std::string>& getRowNames() const;
    void setColNames(const std::vector<std::string>& names);
    void setRowNames(const std::vector<std::string>& names);
    [[nodiscard]] std::vector<size_t> colIndices(const std::vector<std::string> &colNames) const;
    [[nodiscard]] std::vector<size_t> rowIndices(const std::vector<std::string> &rowNames) const;
    [[nodiscard]] size_t getRowIndex(const std::string& name) const;
    [[nodiscard]] size_t getColIndex(const std::string& name) const;

    // Data modification
    void addRow(const std::vector<std::string>& row,
                const std::optional<std::string>& rowName = std::nullopt);
    void addColumn(const std::vector<std::string>& col,
                   const std::optional<std::string>& colName = std::nullopt);



    // // Conversion to Matrix
    // template<typename Scalar, int Options = Eigen::ColMajor>
    // [[nodiscard]] Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> toEigenMatrix() const;
    // Conversion to NamedMatrix
    template<typename Scalar, int Options = Eigen::ColMajor>
    [[nodiscard]] NamedMatrix<Scalar, Options> toNamedMatrix() const;
    template<typename Scalar, int Options = Eigen::ColMajor>
    NamedVector<Scalar> toNamedVector() const;
    template<typename Scalar, int Options = Eigen::RowMajor>
    NamedRowVector<Scalar> toNamedRowVector() const;

    // Printing
    void print(std::ostream& os = std::cout) const;

private:
    std::vector<std::vector<std::string>> data_;
    std::optional<std::vector<std::string>> colNames_;
    std::optional<std::vector<std::string>> rowNames_;

    void validateDimensions(bool tolerateMissingLastColumn=false) const;
    void checkBounds(size_t row, size_t col) const;
};

#endif // TABLE_H
