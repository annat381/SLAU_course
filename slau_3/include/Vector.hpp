#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>
#include <cstddef>
#include <cmath>

class Vector {
private:
    size_t size_;
    std::vector<double> data_;

public:
    Vector(size_t size);
    
    double& operator[](size_t i);
    const double& operator[](size_t i) const;
    
    size_t size() const;
    
    Vector operator+(const Vector& other) const;
    Vector operator*(double scalar) const;
    double operator*(const Vector& other) const; 
    
    double norm() const;  
};

#endif // VECTOR_HPP
