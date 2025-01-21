#ifndef HF_FNS_H
#define HF_FNS_H

#include "./Eigen/Dense"
#include <string>

// Define types for convenience
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> longMatrix;
typedef Eigen::Matrix<long double, Eigen::Dynamic, 1> longVector;
typedef Eigen::SelfAdjointEigenSolver<longMatrix> matrixSolver;

template<typename Derived>
void cleanMatrix(Eigen::MatrixBase<Derived>& matrix) {  // removed default value
    matrix = (matrix.array().abs() < 1e-10).select(0.0, matrix);
}


// Function declarations
void file_error_msg(std::fstream& inputFile);
longMatrix initialize(int matrixSize);
std::pair<longMatrix, longVector> solve_eigen(longMatrix m);
void cleanVector(longVector vector);
double calculateRMSDensity(const longMatrix& D_new, const longMatrix& D_old);
// Class declarations
class MatrixFromFile {
public:
    static longMatrix run(const std::string& fileName);
private:
    static longMatrix get_data(std::fstream& inputFile, int row, int col, double value);
};

class double_integral {
public:
    static uint64_t compoundIndex(uint64_t i, uint64_t j, uint64_t k, uint64_t l);
    static longVector readDoubleFromFile(const std::string& fileName);
// private:
    static uint64_t calc_ij(uint64_t i, uint64_t j);
    static void resize_storage(int index, longVector& vector);
    static void apply_symmetries(longVector& doubles, uint64_t i, uint64_t j, uint64_t k, uint64_t l, long double val);
    static longVector get_data(std::fstream& inputFile, std::string line, longVector& storage);
};

class ComputeNewFock {
public:
    static long double sumDoubles(longMatrix D, longVector doubles, uint64_t mu, uint64_t nu);
    static longMatrix computeFock(longMatrix H, longMatrix D, longVector doubles);
};

#endif // HF_FNS_H
