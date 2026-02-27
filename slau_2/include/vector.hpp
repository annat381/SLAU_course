#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>
#include <stdexcept>

class Vector {
private:
    size_t size_;
    std::vector<double> data_;

public:
    Vector(size_t size) : size_(size), data_(size, 0.0) {}

    double& operator[](size_t i) { return data_[i]; }
    const double& operator[](size_t i) const { return data_[i]; }

    size_t size() const { return size_; }

    Vector operator+(const Vector& other) const {
        if (size_ != other.size_) throw std::runtime_error("Bad sizes");
        Vector result(size_);
        for (size_t i = 0; i < size_; ++i) {
            result[i] = data_[i] + other[i];
        }
        return result;
    }

    Vector operator*(double scalar) const {
        Vector result(size_);
        for (size_t i = 0; i < size_; ++i) {
            result[i] = data_[i] * scalar;
        }
        return result;
    }

    double operator*(const Vector& other) const {
        if (size_ != other.size_) throw std::runtime_error("Bad size!");
        double result = 0.0;
        for (size_t i = 0; i < size_; ++i) {
            result += data_[i] * other[i];
        }
        return result;
    }
};

#endif 
