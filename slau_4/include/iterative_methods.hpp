#pragma once
#include "csr_matrix.hpp"
#include "vector.hpp"

Vector Jacobi(const CSRMatrix& A, const Vector& b, double epsilon, size_t max_iter);

Vector Simple(const CSRMatrix& A, const Vector& b, const Vector& x0
                       double tau, double tolerance);

Vector GaussZeidel(const CSRMatrix& A, const Vector& b, double epsilon, size_t max_iter);
