#include "Matrix.hpp"
#include <cmath>
#include <algorithm>

std::pair<Matrix, Matrix> householderQR(Matrix A) {
    const size_t n = A.rows();
    const size_t m = A.cols();
    Matrix Q(n, n);

    for (size_t i = 0; i < n; ++i)
        Q(i, i) = 1.0;

    Matrix R = A;

    for (size_t k = 0; k < std::min(n, m); ++k) {

        size_t sub_size = n - k;
        Vector v(sub_size);
        for (size_t i = 0; i < sub_size; ++i)
            v[i] = R(i + k, k);


        double norm_v = v.norm();
        if (norm_v < 1e-12) continue;

        double sign = (v[0] >= 0) ? 1.0 : -1.0;
        double alpha = -sign * norm_v;


        v[0] -= alpha;
        double beta = v.norm();
        if (beta < 1e-12) continue;
        for (size_t i = 0; i < sub_size; i++)
            v[i] /= beta;

        //  R = (I - 2vv^T) * R
        for (size_t j = k; j < m; ++j) {
            double dot = 0.0;
            for (size_t i = 0; i < sub_size; i++)
                dot += v[i] * R(i + k, j);
            for (size_t i = 0; i < sub_size; i++)
                R(i + k, j) -= 2.0 * v[i] * dot;
        }

        // Q: Q = Q * (I - 2vv^T)
        for (size_t j = 0; j < n; ++j) {
            double dot = 0.0;
            for (size_t i = 0; i < sub_size; i++)
                dot += v[i] * Q(j, i + k);
            for (size_t i = 0; i < sub_size; i++)
                Q(j, i + k) -= 2.0 * dot * v[i];
        }
    }

    return {Q, R};
}

Vector backGauss(const Matrix& R, const Vector& b) {
    const size_t n = b.size();
    Vector x(n);
    const double EPS = 1e-12;

    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (size_t j = i + 1; j < n; j++)
            sum += R(i, j) * x[j];

        double diag = R(i, i);

        x[i] = (b[i] - sum) / diag;
    }
    return x;
}


Vector solveQR(const Matrix& A, const Vector& b) {
    auto [Q, R] = householderQR(A);
    Matrix Qt = Q.transpose();
    Vector d = Qt * b;
    return backGauss(R, d);
}




