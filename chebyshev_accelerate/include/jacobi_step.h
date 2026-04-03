#pragma once
#include "csr_matrix.h"
#include <functional>

inline Vector jacobi_step(const CSRMatrix& A, const Vector& b, const Vector& x) {
    size_t n = A.get_rows();
    Vector xn(n);
    
    for (size_t i = 0; i < n; ++i) {
        double diag = A.get(i, i);
        double sum = b[i];
        
        for (size_t k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
            size_t j = A.col_indices[k];
            if (j != i) {
                sum -= A.values[k] * x[j];
            }
        }
        
        xn[i] = sum / diag;
    }
    return xn;
}

using StepFunction = std::function<Vector(const CSRMatrix&, const Vector&, const Vector&)>;
