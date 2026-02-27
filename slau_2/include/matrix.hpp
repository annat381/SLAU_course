#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <stdexcept>
#include "vector.hpp"

class Matrix {
private:
    size_t rows_;
    size_t cols_;
    std::vector<double> data_;

public:
    Matrix(size_t rows, size_t cols)
        : rows_(rows), cols_(cols), data_(rows * cols, 0.0) {}

    double operator()(size_t i, size_t j) const {
        return data_[i * cols_ + j];
    }

    double& operator()(size_t i, size_t j) {
        return data_[i * cols_ + j];
    }

    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }

    Vector operator*(const Vector& vec) const {
        if (cols_ != vec.size()) {
            throw std::runtime_error("Matrix-Vector bad size ");
        }

        Vector result(rows_);
        for (size_t i = 0; i < rows_; ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < cols_; ++j) {
                sum += data_[i * cols_ + j] * vec[j];
            }
            result[i] = sum;
        }
        return result;
    }
};

#endif
