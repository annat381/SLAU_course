#include "iterative_methods.hpp"
#include <iostream>


Vector Jacobi(const CSRMatrix& A, const Vector& b, double epsilon, size_t max_iter) {

    size_t n = A.get_rows();


    Vector xprev(n);
    Vector xcurr(n);

    size_t iter = 0;
    double error = 0.0;


    while (iter < max_iter) {


        for (size_t i = 0; i < n; i++) {

            double diag = A.get(i, i);

            if (diag == 0.0) {
                throw std::runtime_error("Zero on diagonal");
            }

            double sum = 0.0;


            size_t start = A.row_ptr[i];
            size_t stop = A.row_ptr[i + 1];

            for (size_t k = start; k < stop; k++) {
                size_t j = A.col_indices[k];
                double val = A.values[k];


                if (j != i) {
                    sum += val * xprev[j];
                }
            }


            xcurr[i] = (b[i] - sum) / diag;
        }


        error = (xcurr - xprev).norm();
        if (error < epsilon) {

            return xcurr;
        }


        xprev = xcurr;
        iter++;
    }

    std::cerr << "Warning: method did not converge";
    return xcurr;
}



Vector SimpleIteration(const CSRMatrix& A, const Vector& b, double tau, double epsilon, size_t max_iter) {
    size_t n = A.get_rows();
    Vector x(n);

    for (size_t iter = 0; iter < max_iter; iter++) {

        Vector Ax = A * x;
        Vector residual = Ax - b;


        Vector x_new = x - residual * tau;


        double error = (x_new - x).norm();
        if (error < epsilon) {
            return x_new;
        }

        x = x_new;
    }

    std::cerr << " did not converge ";
    return x;
}


Vector GaussZeidel(const CSRMatrix& A, const Vector& b, double epsilon, size_t max_iter) {

    size_t n = A.get_rows();
    Vector x(n);

    size_t iter = 0;
    double error = 0.0;

    while (iter < max_iter) {


        Vector x_old = x;

        for (size_t i = 0; i < n; i++) {

            double diag = A.get(i, i);

            if (diag == 0.0) {
                throw std::runtime_error("Zero on diagonal");
            }

            double sum = 0.0;

            size_t start = A.row_ptr[i];
            size_t stop = A.row_ptr[i + 1];

            for (size_t k = start; k < stop; k++) {
                size_t j = A.col_indices[k];
                double val = A.values[k];

                if (j != i) {
                    sum += val * x[j];
                }
            }

            x[i] = (b[i] - sum) / diag;
        }


        error = (x - x_old).norm();

        if (error < epsilon) {
            std::cout << "Converged ";
            return x;
        }

        iter++;
    }

    std::cerr << "did not converge ";
    return x;
}



