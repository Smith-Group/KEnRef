#ifndef NAMED_MATRIX_H
#define NAMED_MATRIX_H

#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <stdexcept>
#include <type_traits>
#include <optional>   // std::optional (libstdc++ needs it explicitly; libc++ pulled it in transitively)

template<typename Scalar, int Rows, int Cols, int Options = Eigen::ColMajor>
class NamedMatrixImpl;

// User-facing aliases
template<typename Scalar, int Options = Eigen::ColMajor>
using NamedMatrix = NamedMatrixImpl<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options>;

template<typename Scalar, int Options = Eigen::ColMajor>
using NamedVector = NamedMatrixImpl<Scalar, Eigen::Dynamic, 1, Options>;

template<typename Scalar, int Options = Eigen::RowMajor>
using NamedRowVector = NamedMatrixImpl<Scalar, 1, Eigen::Dynamic, Options>;

// Core implementation
template<typename Scalar, int Rows, int Cols, int Options>
class NamedMatrixImpl : public Eigen::Matrix<Scalar, Rows, Cols, Options> {
public:
    using Base = Eigen::Matrix<Scalar, Rows, Cols, Options>;
    using Base::Base;

    // Construct from Eigen expressions
    template<typename OtherDerived>
    explicit NamedMatrixImpl(const Eigen::MatrixBase<OtherDerived>& other) : Base(other) {
        // std::cout << "Constructor from Eigen::MatrixBase" << std::endl;
        colNames_ = std::nullopt;
        rowNames_ = std::nullopt;
    }

    template<typename OtherDerived>
    explicit NamedMatrixImpl(const Eigen::MatrixBase<OtherDerived>& other,
        const std::optional<std::vector<std::string>>& colNames = std::nullopt,
        const std::optional<std::vector<std::string>>& rowNames = std::nullopt) : Base(other) {
        // std::cout << "Constructor from Eigen::MatrixBase with names" << std::endl;
        colNames_ = colNames;
        rowNames_ = rowNames;
    }

    template<typename OtherScalar, int OtherRows, int OtherCols, int OtherOptions>
    explicit NamedMatrixImpl(const NamedMatrixImpl<OtherScalar, OtherRows, OtherCols, OtherOptions> &other)
        : Base(other) {
        // std::cout << "(old) conversion constructor from another NamedMatrixImpl" << std::endl;
        colNames_ = other.hasColNames() ? std::optional(other.colNames()) : std::nullopt;
        rowNames_ = other.hasRowNames() ? std::optional(other.rowNames()) : std::nullopt;
    }

    // template<typename OtherScalar, int OtherRows, int OtherCols, int OtherOptions>
    // explicit NamedMatrixImpl(const NamedMatrixImpl<OtherScalar, OtherRows, OtherCols, OtherOptions>& other)
    //     : Base(other),
    //       colNames_(other.colNames()),
    //       rowNames_(other.rowNames())
    // {
    //     std::cout << "(new) conversion constructor from another NamedMatrixImpl" << std::endl;
    //     // static_assert(Rows == OtherRows || Rows == Eigen::Dynamic, "Row count mismatch in conversion");
    //     // static_assert(Cols == OtherCols || Cols == Eigen::Dynamic, "Column count mismatch in conversion");
    // }

    // Copy assignment from any compatible Eigen expression
    template<typename OtherDerived>
    NamedMatrixImpl& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
        // std::cout << "operator=(const Eigen::MatrixBase&)" << std::endl;
        Base::operator=(other);
        return *this;
    }

    // Special assignment from another NamedMatrixImpl
    template<typename OtherScalar, int OtherRows, int OtherCols, int OtherOptions>
    NamedMatrixImpl& operator=(const NamedMatrixImpl<OtherScalar, OtherRows, OtherCols, OtherOptions>& other) {
        // std::cout << "operator=(const const NamedMatrixImpl&)" << std::endl;
        Base::operator=(other);
        // Copy names if dimensions match
        if (Rows == OtherRows || Rows == Eigen::Dynamic) {
            rowNames_ = other.rowNames();
        }
        if (Cols == OtherCols || Cols == Eigen::Dynamic) {
            colNames_ = other.colNames();
        }
        return *this;
    }

    template<int R=Rows, int C=Cols>
    std::enable_if_t<R != 1 && C != 1, NamedMatrix<Scalar>>
    transpose() const {
        NamedMatrix<Scalar> result = Base::transpose();
        if (hasColNames()) result.setRowNames(colNames_.value());
        if (hasRowNames()) result.setColNames(rowNames_.value());
        return result;
    }
    template<int R=Rows>
    std::enable_if_t<R == 1, NamedVector<Scalar>>
    transpose() const {
        NamedVector<Scalar> result = Base::transpose();
        if (hasColNames()) result.setRowNames(colNames_.value());
        if (hasRowNames()) result.setColNames(rowNames_.value());
        return result;
    }
    template<int C=Cols>
    std::enable_if_t<C == 1, NamedRowVector<Scalar>>
    transpose() const {
        NamedRowVector<Scalar> result = Base::transpose();
        if (hasColNames()) result.setRowNames(colNames_.value());
        if (hasRowNames()) result.setColNames(rowNames_.value());
        return result;
    }


    // Name management
    void setColNames(const std::vector<std::string>& names) {
        if constexpr (Cols == 1) {
            throw std::logic_error("Cannot set column names on column vectors");
        }
        if (names.size() != this->cols()) {
            throw std::invalid_argument("Column count mismatch");
        }
        colNames_ = names;
    }

    void setRowNames(const std::vector<std::string>& names) {
        if constexpr (Rows == 1) {
            throw std::logic_error("Cannot set row names on row vectors");
        }
        if (names.size() != this->rows()) {
            throw std::invalid_argument("Row count mismatch");
        }
        rowNames_ = names;
    }

    // Vector-specific name management
    template<int R = Rows, int C = Cols>
    std::enable_if_t<R == 1 || C == 1>
    setNames(const std::vector<std::string>& names) {
        if constexpr (C == 1) { // Column vector
            if (names.size() != this->rows()) {
                throw std::invalid_argument("Name count mismatch");
            }
            rowNames_ = names;
        } else { // Row vector
            if (names.size() != this->cols()) {
                throw std::invalid_argument("Name count mismatch");
            }
            colNames_ = names;
        }
    }

    // Matrix element access - only enabled for non-vector types
    template<int R = Rows, int C = Cols>
    std::enable_if_t<R != 1 && C != 1, Scalar&>
    operator()(const std::string& row, const std::string& col) {
        return Base::operator()(getRowIndex(row), getColIndex(col));
    }

    // Const matrix access - only enabled for non-vector types
    template<int R = Rows, int C = Cols>
    std::enable_if_t<R != 1 && C != 1, const Scalar&>
    operator()(const std::string& row, const std::string& col) const {
        return Base::operator()(getRowIndex(row), getColIndex(col));
    }

    // Vector element access
    // const Scalar& operator()(const std::string& name) const;
    template<int R = Rows, int C = Cols>
    std::enable_if_t<R == 1 || C == 1, const Scalar&>
    operator()(const std::string& name) const{
        if constexpr (C == 1) {  // Column vector
            return Base::operator()(getRowIndex(name), 0);
        } else {  // Row vector
            return Base::operator()(0, getColIndex(name));
        }
    }
    // non-const version of Scalar& operator()(const std::string& name)
    template<int R = Rows, int C = Cols>
    std::enable_if_t<R == 1 || C == 1, Scalar &>
    operator()(const std::string &name) {
        // Cast to const to call the const version, then remove constness
        return const_cast<Scalar &>(static_cast<const NamedMatrixImpl *>(this)->operator()(name));
    }


    Scalar& operator()(const size_t row, const size_t col){
        return Base::operator()(row, col);
    }

    const Scalar& operator()(const size_t row, const size_t col) const{
        return Base::operator()(row, col);
    }

    // Const access methods
    const Scalar& operator()(const std::string& row, const std::string& col) const;

    // Name accessors
    [[nodiscard]] bool hasColNames() const { return colNames_.has_value(); }
    [[nodiscard]] bool hasRowNames() const { return rowNames_.has_value(); }
    [[nodiscard]] const std::vector<std::string>& colNames() const { return *colNames_; }
    [[nodiscard]] const std::vector<std::string>& rowNames() const { return *rowNames_; }
    [[nodiscard]] const std::vector<std::string>& names() const {
        if constexpr (Cols == 1) return *rowNames_;
        else return *colNames_;
    }


    // Output operator
    friend std::ostream& operator<<(std::ostream& os, const NamedMatrixImpl& matrix) {
        // Calculate maximum widths for each column
        std::vector<size_t> colWidths(matrix.cols(), 0);

        // Check column names first
        if (matrix.hasColNames()) {
            for (int j = 0; j < matrix.cols(); ++j) {
                colWidths[j] = std::max(colWidths[j], matrix.colNames()[j].size());
            }
        }

        // Check all data cells
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                std::ostringstream oss;
                oss << matrix(i, j);
                colWidths[j] = std::max(colWidths[j], oss.str().size());
            }
        }

        // Calculate row name width if needed
        size_t rowNameWidth = 0;
        if (matrix.hasRowNames()) {
            for (const auto& name : matrix.rowNames()) {
                rowNameWidth = std::max(rowNameWidth, name.size());
            }
            rowNameWidth += 1; // Add padding
        }

        // Print column headers if they exist
        if (matrix.hasColNames()) {
            if (matrix.hasRowNames()) {
                os << std::setw(static_cast<int>(rowNameWidth)) << " ";
            }

            for (int j = 0; j < matrix.cols(); ++j) {
                os << std::setw(static_cast<int>(colWidths[j]) + 2) << std::left << matrix.colNames()[j];
            }
            os << "\n";
        }

        // Print each row
        for (int i = 0; i < matrix.rows(); ++i) {
            // Print row name if it exists
            if (matrix.hasRowNames()) {
                os << std::setw(static_cast<int>(rowNameWidth)) << std::left << matrix.rowNames()[i];
            }

            // Print matrix values
            for (int j = 0; j < matrix.cols(); ++j) {
                os << std::setw(static_cast<int>(colWidths[j]) + 2) << std::left << matrix(i, j);
            }
            os << "\n";
        }

        return os;
    }

private:
    std::optional<std::vector<std::string>> colNames_{};
    std::optional<std::vector<std::string>> rowNames_{};

    // [[nodiscard]] size_t getRowIndex(const std::string& name);
    // [[nodiscard]] size_t getColIndex(const std::string& name);
    [[nodiscard]] size_t getRowIndex(const std::string& name) const;
    [[nodiscard]] size_t getColIndex(const std::string& name) const;
    // [[nodiscard]] const size_t getRowIndex(const std::string &name) const;
    // [[nodiscard]] const size_t getColIndex(const std::string &name) const;
};


#endif // NAMED_MATRIX_H