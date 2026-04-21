#pragma once
#include <vector>
#include <cmath>
#include "csr_matrix.hpp"
bool gmres(const CSRMatrix& A, const Vector& b, Vector& x, double tol, size_t max_iter) {
    size_t n = A.get_rows();
    if (n != A.get_cols() || n != b.size() || n != x.size()) return false;


    Vector r = b - (A * x);
    double rho = r.norm();
    if (rho < tol) return true;


    std::vector<Vector> V(max_iter + 1, Vector(n));
    for (size_t i = 0; i < n; ++i) V[0][i] = r[i] / rho;


    std::vector<std::vector<double>> H(max_iter + 1, std::vector<double>(max_iter, 0.0));


    std::vector<double> g(max_iter + 1, 0.0);
    g[0] = rho;


    std::vector<double> cs(max_iter, 0.0), sn(max_iter, 0.0);

    size_t m = 0;
    for (; m < max_iter; ++m) {

        Vector w = A * V[m];


        for (size_t i = 0; i <= m; ++i) {
            H[i][m] = w * V[i];
            for (size_t k = 0; k < n; ++k) {
                w[k] -= H[i][m] * V[i][k];
            }
        }
        H[m + 1][m] = w.norm();


        if (H[m + 1][m] < 1e-14) {
            m++;
            break;
        }


        double inv_h = 1.0 / H[m + 1][m];
        for (size_t k = 0; k < n; ++k) V[m + 1][k] = w[k] * inv_h;


        for (size_t i = 0; i < m; ++i) {
            double temp = cs[i] * H[i][m] + sn[i] * H[i + 1][m];
            H[i + 1][m] = -sn[i] * H[i][m] + cs[i] * H[i + 1][m];
            H[i][m] = temp;
        }


        double a = H[m][m], b_val = H[m + 1][m];
        double r_val = std::hypot(a, b_val);

        if (r_val < 1e-14) { cs[m] = 1.0; sn[m] = 0.0; }
        else { cs[m] = a / r_val; sn[m] = b_val / r_val; }


        double temp_H = cs[m] * H[m][m] + sn[m] * H[m + 1][m];
        H[m + 1][m] = -sn[m] * H[m][m] + cs[m] * H[m + 1][m];
        H[m][m] = temp_H;

        double temp_g = cs[m] * g[m] + sn[m] * g[m + 1];
        g[m + 1] = -sn[m] * g[m] + cs[m] * g[m + 1];
        g[m] = temp_g;


        if (std::abs(g[m + 1]) < tol) {
            m++;
            break;
        }
    }


    Vector y(m);
    for (int i = m - 1; i >= 0; --i) {
        double sum = g[i];
        for (size_t j = i + 1; j < m; ++j) sum -= H[i][j] * y[j];
        y[i] = sum / H[i][i];
    }


    for (size_t j = 0; j < m; ++j) {
        double coef = y[j];
        for (size_t i = 0; i < n; ++i) x[i] += coef * V[j][i];
    }

    return std::abs(g[m]) < tol;
}



CSRMatrix poissonToCSR(const PoissonMatrix& P) {
    size_t N = P.Nx * P.Ny;
    CSRMatrix csr(N, N);

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            double val = P.data[i * N + j];
            if (std::abs(val) > 1e-14) {
                csr.add_element(i, j, val);
            }
        }
    }
    csr.finalize();
    return csr;
}



