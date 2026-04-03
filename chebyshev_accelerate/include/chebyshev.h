#pragma once
#include "jacobi_step.h"

inline size_t chebyshev_accelerate(StepFunction step, const CSRMatrix& A, const Vector& b,
                                   Vector& x, double tol, double rho, size_t max_iter) {
    size_t n = A.get_rows();
    double sigma = rho * rho;
    
    Vector x_old = x;
    Vector x_curr = step(A, b, x);
    
    double err = (A * x_curr - b).norm();
    if (err < tol) { 
        x = x_curr; 
        return 1; 
    }
    
    double omega = 2.0 / (2.0 - sigma);
    size_t iter = 1;
    
    while (iter < max_iter) {
        Vector x_prev = x_curr;
        Vector x_new = step(A, b, x_curr);
        
        if (iter > 1) {
            omega = 1.0 / (1.0 - sigma * omega / 4.0);
        }
        
        for (size_t i = 0; i < n; ++i) {
            x_curr[i] = x_old[i] + omega * (x_new[i] - x_old[i]);
        }
        
        x_old = x_prev;
        err = (A * x_curr - b).norm();
        
        if (err < tol) { 
            x = x_curr; 
            return iter + 1; 
        }
        ++iter;
    }
    
    x = x_curr;
    return max_iter;
}

inline size_t run_plain(StepFunction step, const CSRMatrix& A, const Vector& b,
                        Vector& x, double tol, size_t max_iter) {
    size_t iter = 0;
    while (iter < max_iter) {
        x = step(A, b, x);
        double err = (A * x - b).norm();
        if (err < tol) break;
        ++iter;
    }
    return iter + 1;
}
