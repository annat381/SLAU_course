#include <gtest/gtest.h>
#include "iterative_methods.hpp"


CSRMatrix createTestMatrix() {
    CSRMatrix A(4, 4);
    
   
    A.add_element(0, 0, 10.0); A.add_element(0, 1, 1.0);
    A.add_element(0, 2, 2.0);  A.add_element(0, 3, 1.0);
    
    A.add_element(1, 0, 1.0);  A.add_element(1, 1, 10.0);
    A.add_element(1, 2, 1.0);  A.add_element(1, 3, 2.0);
    
    A.add_element(2, 0, 2.0);  A.add_element(2, 1, 1.0);
    A.add_element(2, 2, 10.0); A.add_element(2, 3, 1.0);
    
    A.add_element(3, 0, 1.0);  A.add_element(3, 1, 2.0);
    A.add_element(3, 2, 1.0);  A.add_element(3, 3, 10.0);
    
    A.finalize();
    return A;
}

Vector createTestVector() {
    Vector b(4);
    b[0] = 22.0;
    b[1] = 32.0;
    b[2] = 38.0;
    b[3] = 48.0;
    return b;
}


TEST(IterativeMethodsTest, JacobiMethod) {
    CSRMatrix A = createTestMatrix();
    Vector b = createTestVector();
    
    double epsilon = 1e-6;
    size_t max_iter = 1000;
    
    Vector x = Jacobi(A, b, epsilon, max_iter);
    
    // решение [1, 2, 3, 4]
    EXPECT_NEAR(x[0], 1.0, 1e-5);
    EXPECT_NEAR(x[1], 2.0, 1e-5);
    EXPECT_NEAR(x[2], 3.0, 1e-5);
    EXPECT_NEAR(x[3], 4.0, 1e-5);
    
   
    Vector residual = A * x - b;
    EXPECT_LT(residual.norm(), 1e-5);
}


TEST(IterativeMethodsTest, SimpleIterationMethod) {
    CSRMatrix A = createTestMatrix();
    Vector b = createTestVector();
    
    double tau = 0.1;
    double epsilon = 1e-6;
    size_t max_iter = 1000;
    
    Vector x = SimpleIteration(A, b, tau, epsilon, max_iter);
    
    EXPECT_NEAR(x[0], 1.0, 1e-5);
    EXPECT_NEAR(x[1], 2.0, 1e-5);
    EXPECT_NEAR(x[2], 3.0, 1e-5);
    EXPECT_NEAR(x[3], 4.0, 1e-5);
    
    Vector residual = A * x - b;
    EXPECT_LT(residual.norm(), 1e-5);
}


TEST(IterativeMethodsTest, GaussSeidelMethod) {
    CSRMatrix A = createTestMatrix();
    Vector b = createTestVector();
    
    double epsilon = 1e-6;
    size_t max_iter = 1000;
    
    Vector x = GaussZeidel(A, b, epsilon, max_iter);
    
    EXPECT_NEAR(x[0], 1.0, 1e-5);
    EXPECT_NEAR(x[1], 2.0, 1e-5);
    EXPECT_NEAR(x[2], 3.0, 1e-5);
    EXPECT_NEAR(x[3], 4.0, 1e-5);
    
    Vector residual = A * x - b;
    EXPECT_LT(residual.norm(), epsilon);
}


