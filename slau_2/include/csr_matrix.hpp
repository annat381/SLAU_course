#ifndef CSR_MATRIX_HPP
#define CSR_MATRIX_HPP

#include <vector>
#include "vector.hpp"

class CSRMatrix {
private:
    size_t num_rows;
    size_t num_cols;
    
    std::vector<double> values;        // v - значения
    std::vector<size_t> col_indices;   // c - столбцы
    std::vector<size_t> row_ptr;       // r - указатели на строки

public:
    CSRMatrix(size_t rows, size_t cols)
        : num_rows(rows), num_cols(cols) {
        row_ptr.resize(num_rows + 1, 0);
    }

   
    void add_element(size_t row, size_t col, double value) {
        values.push_back(value);
        col_indices.push_back(col);
        row_ptr[row + 1]++;
    }

   
    void finalize() {
        for (size_t i = 0; i < num_rows; ++i) {
            row_ptr[i + 1] += row_ptr[i];
        }
    }

    
    double get(size_t i, size_t j) const {
        size_t start = row_ptr[i];
        size_t end = row_ptr[i + 1];
        
        for (size_t k = start; k < end; ++k) {
            if (col_indices[k] == j) {
                return values[k];
            }
        }
        return 0.0;
    }

    size_t get_rows() const { return num_rows; }
    size_t get_cols() const { return num_cols; }

    
    Vector operator*(const Vector& vec) const {
        Vector result(num_rows);
        
        for (size_t i = 0; i < num_rows; ++i) {
            double sum = 0.0;
            for (size_t k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                size_t col = col_indices[k];
                double val = values[k];
                sum += val * vec[col];
            }
            result[i] = sum;
        }
        
        return result;
    }
};

#endif 
