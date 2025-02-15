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
    for(uint64_t lmda = 0; lmda < size; lmda++) {
        for(uint64_t sig = 0; sig < size; sig++) {
            value += D(lmda, sig) * 
                (2.0*doubles(DoubleIntegral::compoundIndex(mu, nu, lmda, sig))
                - doubles(DoubleIntegral::compoundIndex(mu, lmda, nu, sig)));
        }
    }
    return value;
}

longMatrix FockBuilder::build(const longMatrix& H, const longMatrix& D, const longVector& doubles) {
    longMatrix fock = hf_utils::initialize(H.rows());
    for(uint64_t mu = 0; mu < H.rows(); mu++) {
        for(uint64_t nu = 0; nu < H.rows(); nu++) {
            fock(mu, nu) = H(mu, nu) + sumDoubles(D, doubles, mu, nu);
        }
    }
    return fock;
}

SCF::SCF(const longMatrix& H_in, const longMatrix& S_minus_half_in, 
         const longVector& dintegrals_in, long double nuclear_repulsion_in)
    : H(H_in), S_minus_half(S_minus_half_in), dintegrals(dintegrals_in),
      nuclear_repulsion(nuclear_repulsion_in), converged(false), 
      iteration_count(0), E_electronic(0), E_total(0) {
    
    // Initialize density matrix
    D = longMatrix::Zero(H.rows(), H.rows());
}

longMatrix SCF::buildNewDensity(const longMatrix& fock) {
    // Transform Fock matrix
    
    longMatrix fock_prime = S_minus_half.transpose() * fock * S_minus_half;
    std::cout << "fockprime success" << std::endl;
    // Solve eigenvalue problem
    auto eigen_solution = hf_utils::solve_eigen(fock_prime);
    longMatrix C = S_minus_half * eigen_solution.first;
    
    // Build new density matrix
    longMatrix D_new = longMatrix::Zero(H.rows(), H.rows());
    for(int mu = 0; mu < H.rows(); mu++) {
        for(int nu = 0; nu < H.rows(); nu++) {
            for(int m = 0; m < num_occupied; m++) {
                D_new(mu, nu) += C(mu, m) * C(nu, m);
            }
        }
    }
    
    return D_new;
}

double SCF::calculateNewEnergy(const longMatrix& D_new, const longMatrix& fock) {
    double E_new = 0.0;
    for(int mu = 0; mu < H.rows(); mu++) {
        for(int nu = 0; nu < H.rows(); nu++) {
            E_new += D_new(mu, nu) * (H(mu, nu) + fock(mu, nu));
        }
    }
    return E_new;
}

void SCF::printIterationInfo(int iteration, double rmsd, double delta_E, double E_new) {
    std::cout << "Iteration " << iteration << ":" << std::endl
              << "  RMSD: " << rmsd << std::endl
              << "  ΔE: " << delta_E << std::endl
              << "  E_elec: " << E_new << std::endl
              << "  E_total: " << E_new + nuclear_repulsion << std::endl;
}

void SCF::runSCF() {
    double E_old = 0.0;
    iteration_count = 0;
    
    while(!converged && iteration_count < MAX_ITERATIONS) {
        // Build new Fock matrix
        longMatrix fock = FockBuilder::build(H, D, dintegrals);
        // Build new density matrix
        longMatrix D_new = buildNewDensity(fock);
        std::cout << "dnew success" << std::endl;

        // Calculate new energy
        double E_new = calculateNewEnergy(D_new, fock);
        
        // Check convergence
        double rmsd = hf_utils::calculateRMSDensity(D_new, D);
        double delta_E = std::abs(E_new - E_old);
        
        printIterationInfo(iteration_count, rmsd, delta_E, E_new);
        
        // Check convergence criteria
        if(rmsd < DENSITY_THRESHOLD && delta_E < ENERGY_THRESHOLD) {
            converged = true;
            E_electronic = E_new;
            E_total = E_new + nuclear_repulsion;
        }
        
        // Update for next iteration
        D = D_new;
        E_old = E_new;
        iteration_count++;
    }
    
    if(converged) {
        std::cout << "\nSCF Converged in " << iteration_count << " iterations!" << std::endl
                  << "Final Electronic Energy: " << E_electronic << std::endl
                  << "Final Total Energy: " << E_total << std::endl;
    } else {
        std::cout << "\nSCF Failed to converge in " << MAX_ITERATIONS << " iterations." << std::endl;
    }
}
