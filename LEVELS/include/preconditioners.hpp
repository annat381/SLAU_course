#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "csr_matrix.hpp"

using LevelMatrix = std::vector<std::vector<size_t>>;


inline LevelMatrix computeLevelMatrix(const CSRMatrix& A) {
    size_t n = A.get_rows();
    LevelMatrix lvl(n, std::vector<size_t>(n, n - 1));
    
    
    for (size_t i = 0; i < n; ++i)
        for (size_t k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k)
            lvl[i][A.col_indices[k]] = 0;
    
    
    for (size_t k = 0; k < n; ++k)
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                if (lvl[i][k] < n-1 && lvl[k][j] < n-1)
                    lvl[i][j] = std::min(lvl[i][j], lvl[i][k] + lvl[k][j] + 1);
    
    return lvl;
}


inline CSRMatrix incompleteCholesky0(const CSRMatrix& A) {
    size_t n = A.get_rows();
    CSRMatrix L(n, n);
    L.row_ptr = A.row_ptr;
    L.col_indices = A.col_indices;
    L.values.resize(A.values.size(), 0.0);
    
    std::vector<double> L_row(n, 0.0);
    
    for (size_t i = 0; i < n; i++) {
        std::fill(L_row.begin(), L_row.end(), 0.0);
        
        for (size_t ptr = A.row_ptr[i]; ptr < A.row_ptr[i + 1]; ptr++) {
            size_t j = A.col_indices[ptr];
            if (j > i) continue; 
            
            double sum = A.values[ptr];
            
            
            size_t ptr_i = A.row_ptr[i];
            while (ptr_i < ptr) {
                size_t k = A.col_indices[ptr_i];
                if (k >= j) break;
                
                size_t ptr_j = A.row_ptr[j];
                size_t end_j = A.row_ptr[j + 1];
                while (ptr_j < end_j && A.col_indices[ptr_j] < k) ++ptr_j;
                
                if (ptr_j < end_j && A.col_indices[ptr_j] == k) {
                    sum -= L_row[k] * L.values[ptr_j];
                }
                ptr_i++;
            }
            
            if (i == j) {
                L.values[ptr] = std::sqrt(sum);
                L_row[i] = L.values[ptr];
            } else {
                double l_jj = 0.0;
                for (size_t p = A.row_ptr[j]; p < A.row_ptr[j+1]; p++)
                    if (A.col_indices[p] == j) l_jj = L.values[p];
                
                L.values[ptr] = sum / l_jj;
                L_row[j] = L.values[ptr];
            }
        }
    }
    return L;
}
