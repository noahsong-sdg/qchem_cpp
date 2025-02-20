#include "hf_fns.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>


// Utility function implementations
namespace hf_utils {
    void cleanVector(longVector& vector) {
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

    void cleanMatrix(longMatrix& matrix) {
        matrix = (matrix.array().abs() < CLEAN_THRESHOLD).select(0.0, matrix);
    }
    
    void file_error_msg(std::fstream& inputFile) {
        if (!inputFile.is_open()) {
            throw std::runtime_error("Unable to open file");
        }
    }

    size_t countLines(const std::string& fileName) {
        std::ifstream file(fileName);
        size_t lineCount = 0;
        std::string unused;
        while (std::getline(file, unused)) {
            ++lineCount;
        }
        return lineCount;
    }

    longMatrix initialize(int matrixSize) {
        return longMatrix::Zero(matrixSize, matrixSize);  // Eigen's built-in Zero initialization
    }

    std::pair<longMatrix, longVector> solve_eigen(longMatrix m) {
        matrixSolver solver(m);
        return std::make_pair(solver.eigenvectors(), solver.eigenvalues());
    }

    double calculateRMSDensity(const longMatrix& D_new, const longMatrix& D_old) {
        double sum = 0.0;
        int size = D_new.rows();
        
        for(int mu = 0; mu < size; mu++) {
            for(int nu = 0; nu < size; nu++) {
                double diff = D_new(mu, nu) - D_old(mu, nu);
                sum += diff * diff;
            }
        }
        
        return std::sqrt(sum);
    }
    
    void releaseMemory(longVector& vec) {
        longVector().swap(vec);  // Force deallocation
    }
    
    longMatrix fastMatrixMultiply(const longMatrix& A, const longMatrix& B) {
        return A * B;  // Eigen optimizes this internally
    }
    
    // ...other utility function implementations...
}

// MatrixFromFile implementations
longMatrix MatrixFromFile::get_data(std::fstream& inputFile, int row, int col, double value) {
    longMatrix zeroes = hf_utils::initialize(7);
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
    hf_utils::file_error_msg(inputFile);
    int row, col; double value;
    longMatrix data = get_data(inputFile, row, col, value);
    inputFile.close();
    return data;
}

// DoubleIntegral implementations

uint64_t DoubleIntegral::compoundIndex(uint64_t i, uint64_t j, uint64_t k, uint64_t l) {
    // removed calc_ij usage; something in the calculation is not working
    uint64_t ij, kl;

    // Calculate ij index
    if(i > j) {
        ij = i * (i + 1)/2 + j;
    } else {
        ij = j * (j + 1)/2 + i;
    }

    // Calculate kl index
    if(k > l) {
        kl = k * (k + 1)/2 + l;
    } else {
        kl = l * (l + 1)/2 + k;
    }

    // Calculate final compound index
    if(ij > kl) {
        return ij * (ij + 1)/2 + kl;
    } else {
        return kl * (kl + 1)/2 + ij;
    }
}

uint64_t DoubleIntegral::calculate_array_size(uint64_t basis_size) {
    uint64_t n = basis_size;
    return (n * (n + 1) / 2) * ((n * (n + 1) / 2) + 1) / 2;
}

uint64_t DoubleIntegral::detect_basis_size(const std::string& fileName) {
    std::ifstream file(fileName);
    uint64_t max_index = 0;
    uint64_t i, j, k, l;
    long double value;
    std::string line;
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (iss >> i >> j >> k >> l >> value) {
            max_index = std::max({max_index, i, j, k, l});
        }
    }
    return max_index;
}

std::vector<std::array<uint64_t, 4>> DoubleIntegral::generate_permutations(uint64_t i, uint64_t j, uint64_t k, uint64_t l) {
    std::vector<std::array<uint64_t, 4>> perms;
    perms.reserve(8); // We know we'll have 8 permutations
    
    // Permutations of (i,j) and (k,l)
    perms.push_back({i,j,k,l});
    perms.push_back({j,i,k,l});
    perms.push_back({i,j,l,k});
    perms.push_back({j,i,l,k});
    perms.push_back({k,l,i,j});
    perms.push_back({l,k,i,j});
    perms.push_back({k,l,j,i});
    perms.push_back({l,k,j,i});
    
    return perms;
}

void DoubleIntegral::apply_symmetry_batch(longVector& storage, 
    const std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, long double>>& batch) {
    
    #pragma omp parallel for
    for (size_t b = 0; b < batch.size(); b++) {
        auto [i, j, k, l, value] = batch[b];
        // Don't subtract 1 here anymore, keep original 1-based indices
        auto perms = generate_permutations(i, j, k, l);
        
        for (const auto& p : perms) {
            // Subtract 1 from indices only when calculating final index
            uint64_t idx = compoundIndex(p[0]-1, p[1]-1, p[2]-1, p[3]-1);
            if (idx < storage.size()) {
                #pragma omp atomic
                storage[idx] = value;
            }
        }
    }
}

longVector DoubleIntegral::readDoubleFromFile(const std::string& fileName) {
    // Detect basis size automatically (no need to subtract 1 since we want actual size)
    uint64_t basis_size = detect_basis_size(fileName);
    uint64_t array_size = calculate_array_size(basis_size);
    
    // Pre-allocate storage
    longVector storage = longVector::Zero(array_size);
    
    try {
        std::ifstream file(fileName);
        if (!file.is_open()) {
            throw std::runtime_error("Unable to open file: " + fileName);
        }

        // Batch processing for better performance
        std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, long double>> batch;
        batch.reserve(1000); // Process 1000 integrals at a time
        
        std::string line;
        uint64_t i, j, k, l;
        long double value;

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            if (iss >> i >> j >> k >> l >> value) {
                // Store original 1-based indices in batch
                batch.emplace_back(i, j, k, l, value);
                
                if (batch.size() >= 1000) {
                    apply_symmetry_batch(storage, batch);
                    batch.clear();
                }
            }
        }
        
        // Process remaining integrals
        if (!batch.empty()) {
            apply_symmetry_batch(storage, batch);
        }

        file.close();
    } catch (const std::exception& e) {
        std::cerr << "Error reading integral file: " << e.what() << std::endl;
        throw;
    }

    return storage;
}

// Remove unused methods:
// void double_integral::resize_storage...
// void double_integral::apply_symmetries...
// longVector double_integral::get_data...

// FockBuilder implementations
long double FockBuilder::sumDoubles(const longMatrix& D, const longVector& doubles, uint64_t mu, uint64_t nu) {
    long double value {0};
    const int size = D.rows();
    
    hf_utils::parallel_for(0, size, [&](int lmda) {
        for(uint64_t sig = 0; sig < size; sig++) {
            value += D(lmda, sig) * 
                (2.0*doubles(DoubleIntegral::compoundIndex(mu, nu, lmda, sig))
                - doubles(DoubleIntegral::compoundIndex(mu, lmda, nu, sig)));
        }
    });
    
    return value;
}

longMatrix FockBuilder::build(const longMatrix& H, const longMatrix& D, const longVector& doubles) {
    longMatrix fock = hf_utils::initialize(H.rows());
    
    hf_utils::parallel_for(0, H.rows(), [&](int mu) {
        for(uint64_t nu = 0; nu < H.rows(); nu++) {
            fock(mu, nu) = H(mu, nu) + sumDoubles(D, doubles, mu, nu);
        }
    });
    
    return fock;
}


