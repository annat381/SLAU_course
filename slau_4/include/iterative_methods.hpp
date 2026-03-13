#pragma once
#include "csr_matrix.hpp"
#include "vector.hpp"

Vector Jacobi(const CSRMatrix& A, const Vector& b, double epsilon, size_t max_iter);

Vector SimpleIteration(const CSRMatrix& A, const Vector& b, 
                       double tau, double epsilon, size_t max_iter);

Vector GaussZeidel(const CSRMatrix& A, const Vector& b, double epsilon, size_t max_iter);
