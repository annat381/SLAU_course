#include "Matrix.hpp"
#include <iostream>
#include <iomanip>
#include <stdexcept>

Matrix::Matrix(size_t rows, size_t cols)
    : rows_(rows), cols_(cols), data_(rows * cols, 0.0) {}

double Matrix::operator()(size_t i, size_t j) const {
    if (i >= rows_ || j >= cols_)
        throw std::out_of_range("index out of range");
    return data_[i * cols_ + j];
}

double& Matrix::operator()(size_t i, size_t j) {
    if (i >= rows_ || j >= cols_)
        throw std::out_of_range("index out of range");
    return data_[i * cols_ + j];
}

size_t Matrix::rows() const { return rows_; }
size_t Matrix::cols() const { return cols_; }

Vector Matrix::operator*(const Vector& vec) const {
    if (cols_ != vec.size())
        throw std::runtime_error("size mismatch");
    Vector result(rows_);
    for (size_t i = 0; i < rows_; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < cols_; ++j)
            sum += data_[i * cols_ + j] * vec[j];
        result[i] = sum;
    }
    return result;
}

Matrix Matrix::multiply(const Matrix& other) const {
    if (cols_ != other.rows_)
        throw std::runtime_error("size mismatch");
    Matrix result(rows_, other.cols_);
    for (size_t i = 0; i < rows_; ++i)
        for (size_t j = 0; j < other.cols_; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < cols_; ++k)
                sum += data_[i * cols_ + k] * other.data_[k * other.cols_ + j];
            result(i, j) = sum;
        }
    return result;
}

Matrix Matrix::transpose() const {
    Matrix result(cols_, rows_);
    for (size_t i = 0; i < rows_; ++i)
        for (size_t j = 0; j < cols_; ++j)
            result.data_[j * rows_ + i] = data_[i * cols_ + j];
    return result;
}

