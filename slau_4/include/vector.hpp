#pragma once
#include <vector>
#include <stdexcept>
#include <cmath>

class Vector {
private:
    size_t size_;
    std::vector<double> data_;

public:
    Vector(size_t size) : size_(size), data_(size, 0.0) {}

    double& operator[](size_t i) { return data_[i]; }
    const double& operator[](size_t i) const { return data_[i]; }
    size_t size() const { return size_; }

    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;
    Vector operator*(double scalar) const;
    double operator*(const Vector& other) const;
    double norm() const;
};
