#ifndef HF_FNS_H
#define HF_FNS_H

#include "./Eigen/Dense"
#include "./Eigen/Sparse"
#include "Config.h"
#include "Logger.h"
#include <omp.h>
#include <string>
#include <vector>
#include <array>
#include <tuple>

// Type definitions
using longMatrix = Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>;
using longVector = Eigen::Matrix<long double, Eigen::Dynamic, 1>;
using matrixSolver = Eigen::SelfAdjointEigenSolver<longMatrix>;
using SparseMatrix = Eigen::SparseMatrix<long double>;

// Constants
constexpr double DENSITY_THRESHOLD = 1e-12;    // δ1
constexpr double ENERGY_THRESHOLD = 1e-12;     // δ2
constexpr int MAX_ITERATIONS = 50;
constexpr double CLEAN_THRESHOLD = 1e-10;

// Utility functions
namespace hf_utils {
    void cleanMatrix(longMatrix& matrix);
    void cleanVector(longVector& vector);
    size_t countLines(const std::string& fileName);
    void file_error_msg(std::fstream& inputFile);
    longMatrix initialize(int matrixSize);
    std::pair<longMatrix, longVector> solve_eigen(longMatrix m);
    double calculateRMSDensity(const longMatrix& D_new, const longMatrix& D_old);
    
    // Memory management
    void releaseMemory(longVector& vec);
    
    // Matrix operations
    longMatrix fastMatrixMultiply(const longMatrix& A, const longMatrix& B);
    
    template<typename Func>
    inline void parallel_for(int start, int end, Func f) {
        #pragma omp parallel for num_threads(HFConfig::NUM_THREADS)
        for(int i = start; i < end; i++) {
            f(i);
        }
    }
}

// Core functionality classes
class MatrixFromFile {
public:
    static longMatrix run(const std::string& fileName);
private:
    static longMatrix get_data(std::fstream& inputFile, int row, int col, double value);
};

class DoubleIntegral {
public:
    static longVector readDoubleFromFile(const std::string& fileName);
    static uint64_t compoundIndex(uint64_t i, uint64_t j, uint64_t k, uint64_t l);

private:
    static uint64_t calculate_array_size(uint64_t basis_size);
    static uint64_t detect_basis_size(const std::string& fileName);
    static void apply_symmetry_batch(longVector& storage, 
                                   const std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, long double>>& batch);
    static std::vector<std::array<uint64_t, 4>> generate_permutations(uint64_t i, uint64_t j, uint64_t k, uint64_t l);
};

class FockBuilder {
public:
    static longMatrix build(const longMatrix& H, const longMatrix& D, const longVector& doubles);
private:
    static long double sumDoubles(const longMatrix& D, const longVector& doubles, uint64_t mu, uint64_t nu);
};

#endif // HF_FNS_H
