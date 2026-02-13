#pragma once

#include <vector>
#include <stdexcept>


using namespace std;


class TridiagonalMatrix {
    public:
    vector<double> a; //lower
    vector<double> b ;//main
    vector<double> c; //upper
    size_t n;

    TridiagonalMatrix(size_t n) : n(n), b(n), c(n-1), a(n-1) {}

};


vector<double> solve_tridiagonal(
    const TridiagonalMatrix& matrix,
    const vector<double>& d
);



bool is_close(const vector<double>& a, const vector<double>& b, double eps = 1e-9);

bool diag_dom(const TridiagonalMatrix& matrix);



