#ifndef HF_FNS_H
#define HF_FNS_H

#include "./Eigen/Dense"
#include <string>
#include <vector>
#include <array>
#include <tuple>

// Type definitions
using longMatrix = Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>;
using longVector = Eigen::Matrix<long double, Eigen::Dynamic, 1>;
using matrixSolver = Eigen::SelfAdjointEigenSolver<longMatrix>;

// Constants
constexpr double DENSITY_THRESHOLD = 1e-12;    // δ1
constexpr double ENERGY_THRESHOLD = 1e-12;     // δ2
constexpr int MAX_ITERATIONS = 50;
constexpr double CLEAN_THRESHOLD = 1e-10;

// Utility functions
namespace hf_utils {
    void cleanMatrix(longMatrix& matrix); // Only declare it here
    
    void cleanVector(longVector& vector);
    size_t countLines(const std::string& fileName);
    void file_error_msg(std::fstream& inputFile);
    longMatrix initialize(int matrixSize);
    std::pair<longMatrix, longVector> solve_eigen(longMatrix m);
    double calculateRMSDensity(const longMatrix& D_new, const longMatrix& D_old);
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

class SCF { 
public:
    SCF(const longMatrix& H, const longMatrix& S_minus_half, 
        const longVector& dintegrals, long double nuclear_repulsion);
    
    void runSCF();
    double getFinalEnergy() const { return E_total; }
    int getIterationCount() const { return iteration_count; }
    bool hasConverged() const { return converged; }

private:
    // Core matrices
    longMatrix H;
    longMatrix S_minus_half;
    longVector dintegrals;
    longMatrix D;  // Current density matrix
    
    // State variables
    bool converged;
    int iteration_count;
    double E_electronic;
    double E_total;
    long double nuclear_repulsion;
    static constexpr int num_occupied = 5;  // Could make this configurable

    // Private methods
    longMatrix buildNewDensity(const longMatrix& fock);
    double calculateNewEnergy(const longMatrix& D_new, const longMatrix& fock);
    void printIterationInfo(int iteration, double rmsd, double delta_E, double E_new);
};

#endif // HF_FNS_H
