#include "vector.hpp"

Vector Vector::operator+(const Vector& other) const {
    if (size_ != other.size_) throw std::runtime_error("Bad sizes");
    Vector result(size_);
    for (size_t i = 0; i < size_; ++i)
        result[i] = data_[i] + other[i];
    return result;
}

Vector Vector::operator-(const Vector& other) const {
    if (size_ != other.size_)
        throw std::runtime_error("sizes mismatch");
    Vector result(size_);
    for (size_t i = 0; i < size_; ++i)
        result[i] = data_[i] - other[i];
    return result;
}

Vector Vector::operator*(double scalar) const {
    Vector result(size_);
    for (size_t i = 0; i < size_; ++i)
        result[i] = data_[i] * scalar;
    return result;
}

double Vector::operator*(const Vector& other) const {
    if (size_ != other.size_) throw std::runtime_error("Bad size!");
    double result = 0.0;
    for (size_t i = 0; i < size_; ++i)
        result += data_[i] * other[i];
    return result;
}

double Vector::norm() const {
    double sum = 0.0;
    for (size_t i = 0; i < size_; ++i)
        sum += data_[i] * data_[i];
    return std::sqrt(sum);
}
