#include <gtest/gtest.h>
#include "Matrix.hpp"
#include <cmath>

const double EPS = 1e-10;

TEST(QRTest, SimpleSystem2x2) {
    
    Matrix A(2, 2);
    A(0,0) = 1; A(0,1) = 1;
    A(1,0) = 2; A(1,1) = 1;
    
    Vector b(2);
    b[0] = 5;
    b[1] = 8;
    
    Vector x = solveQR(A, b);
    
    EXPECT_NEAR(x[0], 3.0, EPS);
    EXPECT_NEAR(x[1], 2.0, EPS);
    
    // A*x = b
    Vector Ax = A * x;
    EXPECT_NEAR(Ax[0], 5.0, EPS);
    EXPECT_NEAR(Ax[1], 8.0, EPS);
}

TEST(QRTest, ComplexSystem3x3) {
    
    Matrix A(3, 3);
    A(0,0) = 4; A(0,1) = 1; A(0,2) = 2;
    A(1,0) = 1; A(1,1) = 5; A(1,2) = 3;
    A(2,0) = 2; A(2,1) = 3; A(2,2) = 6;
    
    Vector b(3);
    b[0] = 1;
    b[1] = 2;
    b[2] = 3;
    
    Vector x = solveQR(A, b);
    
    Vector Ax = A * x;
    EXPECT_NEAR(Ax[0], b[0], EPS);
    EXPECT_NEAR(Ax[1], b[1], EPS);
    EXPECT_NEAR(Ax[2], b[2], EPS);
}

TEST(QRTest, DiagonalMatrix) {
    
    Matrix A(2, 2);
    A(0,0) = 2; A(0,1) = 0;
    A(1,0) = 0; A(1,1) = 3;
    
    Vector b(2);
    b[0] = 6;
    b[1] = 12;
    
    Vector x = solveQR(A, b);
    
    EXPECT_NEAR(x[0], 3.0, EPS);
    EXPECT_NEAR(x[1], 4.0, EPS);
}

