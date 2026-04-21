#pragma once
#include <iostream>
#include <vector>
#include <iomanip>
#include "csr_matrix.hpp"
#include "preconditioners.hpp"


inline void printDenseMatrix(const std::string& name, const std::vector<std::vector<double>>& m) {
    std::cout << name;
    for (auto& row : m) {
        for (double v : row) std::cout << std::setw(6) << std::fixed << std::setprecision(2) << v << " ";
        std::cout << "\n";
    }
    std::cout << "\n";
}


inline void printLevelMatrix(const std::string& name, const LevelMatrix& lvl) {
    std::cout << name;
    for (auto& row : lvl) {
        for (size_t v : row) std::cout << std::setw(3) << v << " ";
        std::cout << "\n";
    }
    std::cout << "\n";
}


inline void printCSRMatrix(const std::string& name, const CSRMatrix& csr) {
    size_t n = csr.get_rows();
    std::cout << name ;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double val = csr.get(i, j);
            std::cout << std::setw(8) << std::fixed << std::setprecision(4) 
                      << (std::abs(val) < 1e-9 ? 0.0 : val) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}


inline CSRMatrix convertToCSR(const std::vector<std::vector<double>>& dense) {
    size_t n = dense.size();
    CSRMatrix csr(n, n);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            if (std::abs(dense[i][j]) > 1e-14)
                csr.add_element(i, j, dense[i][j]);
    csr.finalize();
    return csr;
}
