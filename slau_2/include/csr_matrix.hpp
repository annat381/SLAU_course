#ifndef CSR_MATRIX_HPP
#define CSR_MATRIX_HPP

#include <vector>
#include <map>
#include <utility>
#include <stdexcept>
#include <algorithm>
#include "vector.hpp"

class CSRMatrix {
private:
    size_t num_rows;
    size_t num_cols;
    
  
    std::vector<double> values;        // v - значения ненулевых элементов
    std::vector<size_t> col_indices;   // c - индексы столбцов
    std::vector<size_t> row_ptr;       // r - указатели на начало строк
    
    
    std::map<std::pair<size_t, size_t>, double> temp_data;
    bool is_built;

    void build_csr() {
        if (is_built) return;
        
       
        values.clear();
        col_indices.clear();
        row_ptr.clear();
        row_ptr.resize(num_rows + 1, 0);
        
       
        for (const auto& entry : temp_data) {
            size_t row = entry.first.first;
            row_ptr[row + 1]++;
        }
        
        
        for (size_t i = 0; i < num_rows; ++i) {
            row_ptr[i + 1] += row_ptr[i];
        }
        
        
        values.resize(temp_data.size());
        col_indices.resize(temp_data.size());
        
     
        std::vector<size_t> row_pos = row_ptr;
        
        for (const auto& entry : temp_data) {
            size_t row = entry.first.first;
            size_t col = entry.first.second;
            double val = entry.second;
            
            size_t idx = row_pos[row];
            values[idx] = val;
            col_indices[idx] = col;
            row_pos[row]++;
        }
        
        
        for (size_t i = 0; i < num_rows; ++i) {
            size_t start = row_ptr[i];
            size_t end = row_ptr[i + 1];
            
            
            for (size_t j = start; j < end - 1; ++j) {
                for (size_t k = j + 1; k < end; ++k) {
                    if (col_indices[j] > col_indices[k]) {
                        std::swap(col_indices[j], col_indices[k]);
                        std::swap(values[j], values[k]);
                    }
                }
            }
        }
        
        is_built = true;
    }

public:
    CSRMatrix(size_t rows, size_t cols)
        : num_rows(rows), num_cols(cols), is_built(false) {}

    void set(size_t i, size_t j, double val) {
        if (i >= num_rows || j >= num_cols) {
            throw std::out_of_range("Index out of range");
        }
        temp_data[{i, j}] = val;
        is_built = false;  
    }

    double get(size_t i, size_t j) const {
        if (i >= num_rows || j >= num_cols) {
            throw std::out_of_range("Index out of range");
        }

        CSRMatrix* self = const_cast<CSRMatrix*>(this);
        self->build_csr();
        
       
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
        if (num_cols != vec.size()) {
            throw std::runtime_error("size bad");
        }
        
        CSRMatrix* self = const_cast<CSRMatrix*>(this);
        self->build_csr();
        
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
