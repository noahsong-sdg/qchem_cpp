#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "hf_fns.h"

void file_error_msg(std::fstream& inputFile) {
    if (!inputFile.is_open()) {
        throw std::runtime_error("Unable to open file: ");
    };
    return;
}

void cleanVector(longVector vector) {
    vector = (vector.array().abs() < 1e-8).select(0.0, vector);
    // Clean inf and nan values
    for(int i = 0; i < vector.size(); i++) {
        if(std::isinf(vector(i)) || 
           std::isnan(vector(i)) || 
           std::abs(vector(i)) < 1e-8 ||    // catch very small numbers
           std::abs(vector(i)) > 1e8 ||     // catch very large numbers
           std::abs(std::log10(std::abs(vector(i)))) > 10) {  // catch extreme exponents
            vector(i) = 0.0;
        }
    }
}

longMatrix initialize(int matrixSize) {
    return longMatrix::Zero(matrixSize, matrixSize);
}

std::pair<longMatrix, longVector> solve_eigen(longMatrix m) {
    matrixSolver solver(m);  
    longMatrix L = solver.eigenvectors();
    longVector lambda = solver.eigenvalues();
    return std::make_pair(L, lambda);
}

// MatrixFromFile implementations
longMatrix MatrixFromFile::get_data(std::fstream& inputFile, int row, int col, double value) {
    longMatrix zeroes = initialize(7);
    while (inputFile >> row >> col >> value) {
        zeroes(row-1, col-1) = value;
        if (row != col) {
            zeroes(col-1, row-1) = value;
        }
    }
    return zeroes;
}

longMatrix MatrixFromFile::run(const std::string& fileName) {
    std::fstream inputFile(fileName);
    file_error_msg(inputFile);
    int row, col; double value;
    longMatrix data = get_data(inputFile, row, col, value);
    inputFile.close();
    return data;
}

// double_integral implementations
uint64_t double_integral::calc_ij(uint64_t i, uint64_t j) {
    return i*(i+1)/2 + j;
}

void double_integral::resize_storage(int index, longVector& vector) {
    if (index >= vector.size()) {
                vector.conservativeResize(index + 1);  // Eigen's dynamic resize
            }
    return;
} 

longVector double_integral::get_data(std::fstream& inputFile, std::string line, longVector& storage) {
    uint64_t i, j, k, l;
    long double val;
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        if (iss >> i >> j >> l >> k >> val) {
            i--; j--; k--; l--;
            uint64_t index = compoundIndex(i, j, k, l);
            resize_storage(index, storage);
            storage[index] = val;
        }
    }
    return storage;
}

uint64_t double_integral::compoundIndex(uint64_t i, uint64_t j, uint64_t k, uint64_t l) {
    uint64_t ij = calc_ij(i, j);
    uint64_t kl = calc_ij(k, l);
    return ij*(ij+1)/2 + kl;
}

longVector double_integral::readDoubleFromFile(const std::string& fileName) {
    std::fstream inputFile(fileName);
    file_error_msg(inputFile);
    std::string line;
    longVector storage(1);
    longVector data = double_integral::get_data(inputFile, line, storage);
    inputFile.close();
    return data;
}

// ComputeNewFock implementations
long double ComputeNewFock::sumDoubles(longMatrix D, longVector doubles, uint64_t mu, uint64_t nu) {
    long double value {0};
    for(uint64_t lmda = 0; lmda < 7; lmda++) {
        for(uint64_t sig = 0; sig < 7; sig++) {
            value += D(lmda, sig) * 
                (2.0*doubles(double_integral::compoundIndex(mu, nu, lmda, sig))
                - doubles(double_integral::compoundIndex(mu, lmda, nu, sig)));
        }
    }
    return value;
}

longMatrix ComputeNewFock::computeFock(longMatrix H, longMatrix D, longVector doubles) {
    longMatrix fock_n = initialize(7);
    for(uint64_t mu = 0; mu < 7; mu++) {
        for(uint64_t nu = 0; nu < 7; nu++) {
            fock_n(mu, nu) = H(mu, nu) + sumDoubles(D, doubles, mu, nu);
        }
    }
    return fock_n;  
}
