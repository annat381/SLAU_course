#include "solver.hpp"
#include <cmath>
using namespace std;



vector<double> solve_tridiagonal(
    const TridiagonalMatrix& matrix,
    const vector<double>& d
){
    size_t n = matrix.n;


    if (matrix.n == 0) {
        throw invalid_argument("Empty!");
    }


    if (d.size() != matrix.n) {
        throw invalid_argument("Wrong d size");
    }


    if (n == 1) {
    return {d[0] / matrix.b[0]};
    }

    vector<double> p(n), q(n), x(n);

    p[0]=-matrix.c[0]/matrix.b[0];
    q[0]=d[0]/matrix.b[0];



    for (int i=1; i<n-1; i++){
        p[i]=-matrix.c[i]/(matrix.a[i-1]*p[i-1]+matrix.b[i]);
        q[i]=(d[i]-matrix.a[i-1]*q[i-1])/(matrix.a[i-1]*p[i-1]+matrix.b[i]);

    }
    q[n-1] = (d[n-1] - matrix.a[n-2] * q[n-2]) / (matrix.b[n-1] + matrix.a[n-2] * p[n-2]);

    x[n-1] = q[n-1];

    for (int i=n-2; i>=0; i--){
            x[i]=p[i]*x[i+1]+q[i];
    }


    return x;
}


bool is_close(const vector<double>& a, const vector<double>& b, double eps) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (abs(a[i] - b[i]) > eps) return false;
    }
    return true;
}

bool diag_dom(const TridiagonalMatrix& matrix) {
    size_t n = matrix.n;
    if (n == 1) return abs(matrix.b[0]) > 1e-12;

    if (abs(matrix.b[0]) + 1e-12 < abs(matrix.c[0]))
        return false;


    for (size_t i = 1; i < n - 1; i++) {
        if (abs(matrix.b[i]) + 1e-12 < abs(matrix.a[i-1]) + abs(matrix.c[i]))
            return false;
    }

    if (abs(matrix.b[n-1]) + 1e-12 < abs(matrix.a[n-2]))
        return false;

    return true;
}

