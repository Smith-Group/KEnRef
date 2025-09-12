#include "core/NamedMatrix.h"

// Name lookup implementations
template<typename Scalar, int Rows, int Cols, int Options>
size_t NamedMatrixImpl<Scalar, Rows, Cols, Options>::getRowIndex(const std::string &name) const {
    if constexpr (Rows != 1) { // Not a row vector
        if (!rowNames_) throw std::runtime_error("No row names set");
        const auto it = std::find(rowNames_->begin(), rowNames_->end(), name);
        if (it == rowNames_->end()) throw std::runtime_error("Row name not found: " + name);
        return std::distance(rowNames_->begin(), it);
    } else { // Row vector - use column names
        if (!colNames_) throw std::runtime_error("No names set");
        const auto it = std::find(colNames_->begin(), colNames_->end(), name);
        if (it == colNames_->end()) throw std::runtime_error("Name not found: " + name);
        return std::distance(colNames_->begin(), it);
    }
}

template<typename Scalar, int Rows, int Cols, int Options>
size_t NamedMatrixImpl<Scalar, Rows, Cols, Options>::getColIndex(const std::string& name) const {
    if constexpr (Cols == 1) { // Column vector - use row names
        if (!rowNames_) throw std::runtime_error("No names set");
        auto it = std::find(rowNames_->begin(), rowNames_->end(), name);
        if (it == rowNames_->end()) throw std::runtime_error("Name not found: " + name);
        return std::distance(rowNames_->begin(), it);
    } else { // Row vector or matrix - use column names
        if (!colNames_) throw std::runtime_error("No column names set");
        auto it = std::find(colNames_->begin(), colNames_->end(), name);
        if (it == colNames_->end()) throw std::runtime_error("Column name not found: " + name);
        return std::distance(colNames_->begin(), it);
    }
}


// Const access implementations
template<typename Scalar, int Rows, int Cols, int Options>
const Scalar& NamedMatrixImpl<Scalar, Rows, Cols, Options>::operator()(
    const std::string& row, const std::string& col) const{
    return Base::operator()(getRowIndex(row), getColIndex(col));
}


// Explicit instantiation for common types
template class NamedMatrixImpl<int, Eigen::Dynamic, Eigen::Dynamic>;
template class NamedMatrixImpl<float, Eigen::Dynamic, Eigen::Dynamic>;
template class NamedMatrixImpl<double, Eigen::Dynamic, Eigen::Dynamic>;
template class NamedMatrixImpl<int, Eigen::Dynamic, 1>;
template class NamedMatrixImpl<float, Eigen::Dynamic, 1>;
template class NamedMatrixImpl<double, Eigen::Dynamic, 1>;
template class NamedMatrixImpl<int, 1, Eigen::Dynamic, Eigen::RowMajor>;
template class NamedMatrixImpl<float, 1, Eigen::Dynamic, Eigen::RowMajor>;
template class NamedMatrixImpl<double, 1, Eigen::Dynamic, Eigen::RowMajor>;
