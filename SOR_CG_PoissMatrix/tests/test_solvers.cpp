#include <gtest/gtest.h>
#include "../include/vector.hpp"
#include "../include/csr_matrix.hpp"
#include "../include/poisson_matrix.hpp"
#include "../include/solvers.hpp"


TEST(SolverTests, SOR_Method) {
    const size_t N = 3;
    CSRMatrix A(N, N);

    A.add_element(0, 0,  4.0); A.add_element(0, 1, -1.0);
    A.add_element(1, 0, -1.0); A.add_element(1, 1,  4.0); A.add_element(1, 2, -1.0);
    A.add_element(2, 1, -1.0); A.add_element(2, 2,  4.0);
    A.finalize();

    Vector b(N);
    b[0] = 3.0; b[1] = 2.0; b[2] = 3.0;

    double omega = 1.5;
    double eps = 1e-8;
    size_t max_iter = 1000;

    Vector x = SOR(A, b, eps, max_iter, omega);

   
    EXPECT_NEAR(x[0], 1.0, 1e-6);
    EXPECT_NEAR(x[1], 1.0, 1e-6);
    EXPECT_NEAR(x[2], 1.0, 1e-6);

    
    Vector residual = A * x - b;
    EXPECT_LT(residual.norm(), 1e-6);
}


TEST(SolverTests, FastDown_Method) {
    const size_t H = 3;
    CSRMatrix B(H, H);

    B.add_element(0, 0, 4.0); B.add_element(0, 1, -1.0);
    B.add_element(1, 0, -1.0); B.add_element(1, 1, 4.0); B.add_element(1, 2, -1.0);
    B.add_element(2, 1, -1.0); B.add_element(2, 2, 4.0);
    B.finalize();

    Vector b_(H);
    b_[0] = 3.0; b_[1] = 2.0; b_[2] = 3.0;

    double eps_ = 1e-8;
    size_t max_iter_ = 1000;

    Vector x_ = FastDown(B, b_, eps_, max_iter_);

   
    EXPECT_NEAR(x_[0], 1.0, 1e-6);
    EXPECT_NEAR(x_[1], 1.0, 1e-6);
    EXPECT_NEAR(x_[2], 1.0, 1e-6);

   
    Vector residual = B * x_ - b_;
    EXPECT_LT(residual.norm(), 1e-6);
}


TEST(SolverTests, CG_Method) {
    const size_t DIM = 3;
    CSRMatrix sparseMatrix(DIM, DIM);

    sparseMatrix.add_element(0, 0, 4.0); sparseMatrix.add_element(0, 1, -1.0);
    sparseMatrix.add_element(1, 0, -1.0); sparseMatrix.add_element(1, 1, 4.0); sparseMatrix.add_element(1, 2, -1.0);
    sparseMatrix.add_element(2, 1, -1.0); sparseMatrix.add_element(2, 2, 4.0);
    sparseMatrix.finalize();

    Vector rhsVector(DIM);
    rhsVector[0] = 3.0;
    rhsVector[1] = 2.0;
    rhsVector[2] = 3.0;

    double tolerance = 1e-8;
    size_t maxSteps = 100;

    Vector solutionVec = CG(sparseMatrix, rhsVector, tolerance, maxSteps);

    
    EXPECT_NEAR(solutionVec[0], 1.0, 1e-6);
    EXPECT_NEAR(solutionVec[1], 1.0, 1e-6);
    EXPECT_NEAR(solutionVec[2], 1.0, 1e-6);

   
    Vector residual = sparseMatrix * solutionVec - rhsVector;
    EXPECT_LT(residual.norm(), 1e-6);
}


TEST(PoissonMatrixTest, GenerateAndPrint) {
    PoissonMatrix matrix(3, 3);
    matrix.Generate();
    
    
     matrix.Print();
    
    EXPECT_TRUE(true);
}
