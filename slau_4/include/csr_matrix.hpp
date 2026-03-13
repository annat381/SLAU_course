#pragma once
#include <vector>
#include "vector.hpp"

class CSRMatrix {
private:
    size_t num_rows;
    size_t num_cols;

public:
    std::vector<double> values;
    std::vector<size_t> col_indices;
    std::vector<size_t> row_ptr;

    CSRMatrix(size_t rows, size_t cols);
    
    void add_element(size_t row, size_t col, double value);
    void finalize();
    double get(size_t i, size_t j) const;
    size_t get_rows() const;
    size_t get_cols() const;
    
    Vector operator*(const Vector& vec) const;
};
