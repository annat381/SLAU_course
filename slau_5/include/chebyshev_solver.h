#pragma once
#include "csr_matrix.h"
#include <vector>
#include <cmath>
#include <iostream>

inline double powerMethodMaxEigenvalue(const CSRMatrix& A, size_t max_iter = 1000, double tol = 1e-10) {
    size_t n = A.get_rows();
    Vector x(n);

    for (size_t i = 0; i < n; ++i) {
        x[i] = 1.0;
    }

    double lambda = 0.0;

    for (size_t iter = 0; iter < max_iter; iter++) {
        Vector Ax = A * x;

        double new_lambda = (x * Ax) / (x * x);

        double norm = Ax.norm();
        if (norm < 1e-15) break;

        for (size_t i = 0; i < n; ++i) {
            x[i] = Ax[i] / norm;
        }

        if (std::abs(new_lambda - lambda) < tol) {
            return new_lambda;
        }
        lambda = new_lambda;
    }

    return lambda;
}

inline double powerMethodMinEigenvalue(const CSRMatrix& A, double lambda_max, size_t max_iter = 1000, double tol = 1e-10) {
    size_t n = A.get_rows();

    CSRMatrix shifted(n, n);
    shifted.row_ptr.resize(n + 1, 0);

    for (size_t i = 0; i < n; ++i) {
        for (size_t k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
            size_t j = A.col_indices[k];
            double val = A.values[k];
            if (i == j) {
                shifted.add_element(i, j, lambda_max - val);
            } else {
                shifted.add_element(i, j, -val);
            }
        }
    }
    shifted.finalize();

    double lambda_shifted_max = powerMethodMaxEigenvalue(shifted, max_iter, tol);

    return lambda_max - lambda_shifted_max;
}

inline std::vector<double> computeChebyshevRoots(size_t n, double lambda_min, double lambda_max) {
    std::vector<double> roots(n);

    for (size_t k = 0; k < n; ++k) {
        double theta = M_PI * (2 * k + 1) / (2 * n);
        roots[k] = 0.5 * (lambda_max + lambda_min) +
                   0.5 * (lambda_max - lambda_min) * std::cos(theta);
    }

    std::vector<double> permuted_roots(n);
    size_t left = 0, right = n - 1;
    for (size_t i = 0; i < n; ++i) {
        if (i % 2 == 0) {
            permuted_roots[i] = roots[left++];
        } else {
            permuted_roots[i] = roots[right--];
        }
    }

    return permuted_roots;
}

inline Vector ChebyshevAcceleration(const CSRMatrix& A, const Vector& b,
                             double epsilon, size_t max_cycles,
                             size_t n_roots = 5,
                             double tau_init = 0.0) {
    size_t n = A.get_rows();

    double lambda_max = powerMethodMaxEigenvalue(A);
    double lambda_min = powerMethodMinEigenvalue(A, lambda_max);

    double tau = tau_init;
    if (tau <= 0) {
        tau = 2.0 / (lambda_min + lambda_max);
    }

    Vector x(n);
    size_t total_iterations = 0;

    for (size_t cycle = 0; cycle < max_cycles; ++cycle) {
        std::vector<double> cheb_roots = computeChebyshevRoots(n_roots, lambda_min, lambda_max);

        for (size_t i = 0; i < n_roots; ++i) {
            double root = cheb_roots[i];
            double tau_i = 1.0 / root; 

            Vector Ax = A * x;
            Vector residual = b - Ax;
            Vector x_new = x + residual * tau_i;

            double error = (x_new - x).norm();
            x = x_new;
            total_iterations++;

            if (error < epsilon) {
                std::cout << "нужная точность " << std::endl;
                return x;
            }
        }

        Vector Ax = A * x;
        Vector residual = b - Ax;
        double residual_norm = residual.norm();

        if (residual_norm < epsilon) {
            std::cout << "есть точность" << std::endl;
            return x;
        }
    }

    std::cerr << "Нет сходимости " << std::endl;
    return x;
}
