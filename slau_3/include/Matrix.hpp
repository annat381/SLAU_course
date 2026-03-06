#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "Vector.hpp"
#include <vector>
#include <cstddef>
#include <utility>

class Matrix {
private:
    size_t rows_;
    size_t cols_;
    std::vector<double> data_;

public:
    Matrix(size_t rows, size_t cols);
    
    double operator()(size_t i, size_t j) const;
    double& operator()(size_t i, size_t j);
    
    size_t rows() const;
    size_t cols() const;
    
    Vector operator*(const Vector& vec) const;
    Matrix multiply(const Matrix& other) const;
    Matrix transpose() const;
    
};


std::pair<Matrix, Matrix> householderQR(Matrix A);


Vector solveQR(const Matrix& A, const Vector& b);

#endif // MATRIX_HPP
