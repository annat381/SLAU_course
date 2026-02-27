#ifndef CSR_MATRIX_HPP
#define CSR_MATRIX_HPP

#include <map>
#include <utility>
#include <vector>
#include "vector.hpp"

class CSRMatrix {
private:
    size_t num_rows;
    size_t num_cols;
    std::map<std::pair<size_t, size_t>, double> data;

public:
    CSRMatrix(size_t rows, size_t cols)
        : num_rows(rows), num_cols(cols) {}

    void set(size_t i, size_t j, double val) {
        data[{i, j}] = val;
    }

    double get(size_t i, size_t j) const {
        auto it = data.find({i, j});
        if (it != data.end()) {
            return it->second;
        }
        return double{};
    }

    size_t get_rows() const { return num_rows; }
    size_t get_cols() const { return num_cols; }

    Vector operator*(const Vector& vec) const {
        Vector result(num_rows);
        for (auto it = data.begin(); it != data.end(); ++it) {
            std::pair<size_t, size_t> coords = it->first;
            size_t row = coords.first;
            size_t col = coords.second;
            double value = it->second;
            result[row] += value * vec[col];
        }
        return result;
    }
};

#endif
