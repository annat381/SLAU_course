#include <iostream>
#include <iomanip>
#include <cmath>
#include "jacobi_step.h"
#include "chebyshev.h"

int main() {
    CSRMatrix A(3,3);
    A.add_element(0,0,4); A.add_element(0,1,-1);
    A.add_element(1,0,-1); A.add_element(1,1,4); A.add_element(1,2,-1);
    A.add_element(2,1,-1); A.add_element(2,2,4);
    A.finalize();
    
    Vector b(3);
    b[0]=1;
    b[1]=2;
    b[2]=3;
    Vector x0(3);
    
    double tol = 1e-6;
    size_t max_iter = 1000;
    double rho = std::sqrt(2.0) / 4.0;
    
    std::cout << "Chebyshev acceleration test" << std::endl;
    std::cout << "rho = " << rho << std::endl << std::endl;
    
    Vector x_jac = x0;
    size_t it_jac = run_plain(jacobi_step, A, b, x_jac, tol, max_iter);
    
    Vector x_cheb = x0;
    size_t it_cheb = chebyshev_accelerate(jacobi_step, A, b, x_cheb, tol, rho, max_iter);
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Jacobi:    iter=" << it_jac << " x=[" << x_jac[0] << ", " << x_jac[1] << ", " << x_jac[2] << "]" << std::endl;
    std::cout << "Chebyshev: iter=" << it_cheb << " x=[" << x_cheb[0] << ", " << x_cheb[1] << ", " << x_cheb[2] << "]" << std::endl;
    std::cout << "Expected:         x=[0.464286, 0.857143, 0.964286]" << std::endl;
    std::cout << "Speedup: " << (double)it_jac/it_cheb << "x" << std::endl;
    
    return 0;
}
