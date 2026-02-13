#include "solver.hpp"
#include <iostream>

using namespace std;


int main() {


    TridiagonalMatrix mat(3);
    mat.b = {2, 2, 2};
    mat.c = {1, 1};
    mat.a = {1, 1};

    vector<double> d = {3, 4, 3};
    vector<double> result = solve_tridiagonal(mat, d);
    vector<double> expected = {1, 1, 1};

    if (is_close(result, expected)) {
        cout << " Right solution! \n"; 
    } else {
        cout << " Fail\n ";
        for (double v : result) cout << v << " ";
        cout << "\n";
    }
    if (!diag_dom(mat)) {
    cout << "wrong matrix(no diagonal domination)\n";
    }

    if (diag_dom(mat)) {
    cout << "OK! matrix with diag. domination\n";
    }

    return 0;
}


