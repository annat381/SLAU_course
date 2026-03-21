#include <gtest/gtest.h>
#include "chebyshev_solver.h"   



TEST(ChebyshevSolverTest, SolveTridiagonal) {
    size_t n = 5;
    CSRMatrix A(n, n);

    for (size_t i = 0; i < n; ++i) {
        A.add_element(i, i, 4.0);
        if (i > 0) {
            A.add_element(i, i-1, -1.0);
        }
        if (i < n-1) {
            A.add_element(i, i+1, -1.0);
        }
    }
    A.finalize();

    Vector b(n);
    b[0] = 3.0;
    b[1] = 2.0;
    b[2] = 2.0;
    b[3] = 2.0;
    b[4] = 3.0;

    double epsilon = 1e-5;
    size_t max_cycles = 20;
    size_t n_roots = 3;

    Vector x = ChebyshevAcceleration(A, b, epsilon, max_cycles, n_roots);

    
    for (size_t i = 0; i < n; ++i) {
        EXPECT_NEAR(x[i], 1.0, 1e-5);
    }

   
    Vector Ax = A * x;
    Vector residual = b - Ax;
    EXPECT_LT(residual.norm(), epsilon);
}






int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
