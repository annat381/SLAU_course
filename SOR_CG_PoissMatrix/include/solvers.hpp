#pragma once
#include "csr_matrix.hpp"
#include <cmath>
#include <iostream>

inline Vector SOR(const CSRMatrix& A, const Vector& b, double epsilon,
                  size_t max_iter, double omega = 1.5) {
    size_t n = A.get_rows();
    Vector x(n);
    size_t iter = 0;

    while (iter < max_iter) {
        Vector x_old = x;

        for (size_t i = 0; i < n; i++) {
            double diag = A.get(i, i);
            double sum = 0.0;
            size_t start = A.row_ptr[i];
            size_t stop = A.row_ptr[i + 1];

            for (size_t k = start; k < stop; k++) {
                size_t j = A.col_indices[k];
                double val = A.values[k];
                if (j != i) {
                    sum += val * x[j];
                }
            }

            double gauss_seidel_value = (b[i] - sum) / diag;
            x[i] = (1.0 - omega) * x[i] + omega * gauss_seidel_value;
        }

        double error = (x - x_old).norm();
        if (error < epsilon) {
            return x;
        }
        iter++;
    }
    return x;
}

inline Vector FastDown(const CSRMatrix& A, const Vector& b,
                       double epsilon, size_t max_iter) {
    size_t n = A.get_rows();
    Vector x(n);

    for (size_t iter = 0; iter < max_iter; ++iter) {
        Vector Ax = A * x;
        Vector r(n);
        for (size_t i = 0; i < n; ++i) {
            r[i] = b[i] - Ax[i];
        }

        double r_norm = r.norm();
        if (r_norm < epsilon) {
            return x;
        }

        Vector Ar = A * r;
        double r_dot_r = 0.0;
        double r_dot_Ar = 0.0;
        for (size_t i = 0; i < n; ++i) {
            r_dot_r += r[i] * r[i];
            r_dot_Ar += r[i] * Ar[i];
        }

        double alpha = r_dot_r / r_dot_Ar;
        for (size_t i = 0; i < n; ++i) {
            x[i] += alpha * r[i];
        }
    }
    return x;
}

inline Vector CG(const CSRMatrix& A, const Vector& b,
                 double epsilon, size_t max_iter) {
    size_t n = A.get_rows();
    Vector x(n);
    Vector r = b;
    Vector p = r;
    double r_dot_r = r * r;

    for (size_t iter = 0; iter < max_iter; ++iter) {
        if (std::sqrt(r_dot_r) < epsilon) {
            return x;
        }

        Vector Ap = A * p;
        double p_dot_Ap = p * Ap;
        double alpha = r_dot_r / p_dot_Ap;

        for (size_t i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
        }

        for (size_t i = 0; i < n; ++i) {
            r[i] -= alpha * Ap[i];
        }

        double r_new_dot_r = r * r;
        if (std::sqrt(r_new_dot_r) < epsilon) {
            return x;
        }

        double beta = r_new_dot_r / r_dot_r;
        for (size_t i = 0; i < n; ++i) {
            p[i] = r[i] + beta * p[i];
        }

        r_dot_r = r_new_dot_r;
    }
    return x;
}
