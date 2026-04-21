#include <iostream>
#include <vector>
#include "utils.hpp"
#include "poisson_matrix.hpp"
#include "preconditioners.hpp"

int main() {
  

   
    
    std::vector<std::vector<double>> A_levels = {
        {5, 0, 1, 0, 0},
        {1, 5, 1, 0, 1},
        {0, 1, 5, 1, 1},
        {0, 0, 1, 5, 1},
        {1, 1, 0, 0, 5}
    };
    
    printDenseMatrix("Matrix A", A_levels);
    CSRMatrix csrA = convertToCSR(A_levels);
    auto levels = computeLevelMatrix(csrA);
    printLevelMatrix("Level Matrix", levels);
    
   
    
    std::vector<std::vector<double>> A_spd = {
        {4.0, 1.0, 0.0},
        {1.0, 3.0, 1.0},
        {0.0, 1.0, 2.0}
    };
    
    printDenseMatrix("SPD Matrix", A_spd);
    CSRMatrix csrSpd = convertToCSR(A_spd);
    
    try {
        CSRMatrix L = incompleteCholesky0(csrSpd);
        printCSRMatrix("Matrix L (Cholesky)", L);
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
    }
    
    
    PoissonMatrix P(2, 2);
    P.Generate();
    
    
    std::vector<std::vector<double>> poissonDense(P.N, std::vector<double>(P.N));
    for(size_t i=0; i<P.N; ++i)
        for(size_t j=0; j<P.N; ++j)
            poissonDense[i][j] = P.data[i * P.N + j];

    printDenseMatrix("Poisson Matrix", poissonDense);
    
  
    CSRMatrix csrPoisson = convertToCSR(poissonDense);
    
    try {
        CSRMatrix LP = incompleteCholesky0(csrPoisson);
        printCSRMatrix("Matrix L (Poisson)", LP);
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
    }

    return 0;
}
