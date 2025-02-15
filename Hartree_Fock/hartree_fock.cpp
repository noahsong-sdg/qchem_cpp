#include "hartree_fock.h"
#include <fstream>
#include <iostream>

HartreeFockSolver::HartreeFockSolver(const std::string& data_directory) 
    : data_path(data_directory), converged(false), final_energy(0.0) {
    loadIntegrals();
}

void HartreeFockSolver::loadIntegrals() {
    // Load nuclear repulsion
    std::fstream enuc_file(data_path + "/enuc.dat");
    enuc_file >> enuc;
    enuc_file.close();
    
    // Load one-electron integrals
    longMatrix T = MatrixFromFile::run(data_path + "/t.dat");
    longMatrix V = MatrixFromFile::run(data_path + "/v.dat");
    S = MatrixFromFile::run(data_path + "/s.dat");
    H = T + V;
    // Load two-electron integrals
    dintegrals = DoubleIntegral::readDoubleFromFile(data_path + "/eri.dat");
}

double HartreeFockSolver::solve() {
    buildInitialGuess();
    runSCF();
    return final_energy;
}

void HartreeFockSolver::buildInitialGuess() {
    // Get eigenvectors and eigenvalues of overlap matrix
    auto s_ortho = hf_utils::solve_eigen(S);
    longMatrix L = s_ortho.first;
    longVector lambda = s_ortho.second;
    // Calculate S^(-1/2)
    // ISSUE IS HAPPENING HERE
    longVector lambda_minus_half = lambda.array().pow(-0.5);
    longMatrix lambda_minus_half_matrix = lambda_minus_half.asDiagonal();

    S_minus_half = L * lambda_minus_half_matrix * L.transpose();  // Remove local declaration
    hf_utils::cleanMatrix(S_minus_half);
    // Build initial Fock matrix and transform to orthogonal basis
    longMatrix fock_0 = S_minus_half.transpose() * H * S_minus_half;
    hf_utils::cleanMatrix(fock_0);

    // Solve for initial C coefficients
    auto f_ortho = hf_utils::solve_eigen(fock_0);
    longMatrix C_prime = f_ortho.first;
    // Transform coefficients back to original basis
    longMatrix C = S_minus_half * C_prime;
    hf_utils::cleanMatrix(C);
    
    // Build initial density matrix
    longMatrix D = longMatrix::Zero(H.rows(), H.rows());
    constexpr int num_occupied = 5;  // for H2O
    for(int mu = 0; mu < H.rows(); mu++) {
        for(int nu = 0; nu < H.rows(); nu++) {
            for(int m = 0; m < num_occupied; m++) {
                D(mu, nu) += C(mu, m) * C(nu, m);
            }
        }
    }
}

void HartreeFockSolver::runSCF() {
    std::cout << "initial build succeeded" << std::endl;
    SCF scf_solver(H, S_minus_half, dintegrals, enuc);
    scf_solver.runSCF();

    if(scf_solver.hasConverged()) {
        std::cout << "Final Total Energy: " << scf_solver.getFinalEnergy() << std::endl;
    }
}
