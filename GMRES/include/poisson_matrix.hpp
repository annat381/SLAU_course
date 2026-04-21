#pragma once
#include <iostream>
#include <vector>
#include <iomanip>

class PoissonMatrix {
public:
    size_t Nx;
    size_t Ny;
    size_t N;
    std::vector<double> data;

    PoissonMatrix(size_t nx, size_t ny)
        : Nx(nx), Ny(ny), N(nx * ny) {
        data.resize(N * N, 0.0);
    }

    double& a(size_t row, size_t col) {
        return data[row * N + col];
    }

    void Generate() {
        for (size_t iy = 0; iy < Ny; ++iy) {
            for (size_t ix = 0; ix < Nx; ++ix) {
                size_t i = iy * Nx + ix;
                a(i, i) = 4.0;

                if (ix > 0) {
                    a(i, i - 1) = -1.0;
                }
                if (ix < Nx - 1) {
                    a(i, i + 1) = -1.0;
                }
                if (iy > 0) {
                    a(i, i - Nx) = -1.0;
                }
                if (iy < Ny - 1) {
                    a(i, i + Nx) = -1.0;
                }
            }
        }
    }

    void Print() {
        for (size_t r = 0; r < N; ++r) {
            for (size_t c = 0; c < N; ++c) {
                std::cout << std::setw(4) << a(r, c) << " ";
            }
            std::cout << "\n";
        }
    }
};
