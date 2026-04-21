#include <iomanip>
#include "csr_matrix.hpp"
#include "poisson_matrix.hpp"
#include "gmres.hpp"
int main() {
    size_t nx = 20, ny = 20;
    PoissonMatrix P(nx, ny);
    P.Generate();

    CSRMatrix A = poissonToCSR(P);
    size_t N = nx * ny;

    Vector b(N);
    for (int i=0; i<N; i++){

        b[i]=1.0;
    }
    Vector x(N);

    double tol = 1e-8;
    size_t max_iter = 50;

    bool converged = gmres(A, b, x, tol, max_iter);
    std::cout << "Converged: " << (converged ? "YES" : "NO") << "\n";


    Vector Ax = A * x;
    double res_norm = 0.0;
    for (size_t i = 0; i < N; ++i) res_norm += (Ax[i] - b[i]) * (Ax[i] - b[i]);
    std::cout << "Final residual norm: " << std::sqrt(res_norm) << "\n";


    std::cout << "Solution x (first 10 elements):\n";
    for (size_t i = 0; i < 10; ++i) {
        std::cout << x[i] << " ";
    }

    return 0;
}

