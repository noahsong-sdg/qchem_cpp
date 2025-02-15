#include <iostream>
#include <vector>
#include <fstream>
#include <cstdint>
#include <cmath>
#include <stdexcept>
#include <sstream>

#include "hf_fns.h"
#include "hf_debug_fns.h"

namespace hfu = hf_utils;

int main() 
{
    // Step One: Read Enuc
    std::string fileName = "data/enuc.dat";
    std::fstream inputFile(fileName);
    // Check if the file is open
    hfu::file_error_msg(inputFile);
    long enuc{};
    inputFile >> enuc;
    inputFile.close();
    // std::cout << enuc << std::endl;

    // Step Two: One-Electon Integrals
    // Overlap integral S, KE integral T, nucleus-nucleus attraction V:
    longMatrix s = MatrixFromFile::run("data/s.dat");
    longMatrix t = MatrixFromFile::run("data/t.dat");
    longMatrix v = MatrixFromFile::run("data/v.dat");
    // H = T + V
    longMatrix H = t + v;

    hfu::cleanMatrix(H);
    //std::cout << "Hamiltonian matrix:\n" << H << std::endl;



    // Step Three: Two-Electron Integrals
    // TODO: Speedup https://github.com/CrawfordGroup/ProgrammingProjects/blob/master/Project%2303/hints/hint3-2.md
    std::string fileName_d = "data/eri.dat";
    longVector dintegrals = DoubleIntegral::readDoubleFromFile(fileName_d); 
    //hf_debug::findZeroElements(dintegrals);
    //cleanVector(dintegrals);
    /*for(int i = 0; i < dintegrals.size(); i++) {
        if(std::isinf(dintegrals(i)) || std::isnan(dintegrals(i))) {
            dintegrals(i) = 100.0;
            }
        }*/
     //std::cout << dintegrals << std::endl;


    // Step Four: Build Orthogonalization Matrix (S * L = L * Lambda)
    // L is the matrix of eigenvectors
    // Lambda is the diagonal matrix of corresponding eigvals
    // Declare solver
    auto s_ortho = hf_utils::solve_eigen(s);
    longMatrix L = s_ortho.first;
    longVector lambda = s_ortho.second;
    // Calculate λ^(-1/2)
    longVector lambda_minus_half = lambda.array().pow(-0.5);
    // Construct diagonal matrix from λ^(-1/2)
    longMatrix lambda_minus_half_matrix = lambda_minus_half.asDiagonal();
    // Calculate S^(-1/2) = L λ^(-1/2) L^T
    longMatrix S_minus_half = L * lambda_minus_half_matrix * L.transpose();
    hfu::cleanMatrix(S_minus_half);
    //std::cout << "S^-1/2 matrix\n" << S_minus_half << std::endl;



    // Step Five: build initial density guess
    // Initial Fock matrix
    longMatrix fock_0 = S_minus_half.transpose() * H * S_minus_half;
    hfu::cleanMatrix(fock_0);
    // diagonalize
    auto f_ortho = hfu::solve_eigen(fock_0);
    longMatrix C_prime = f_ortho.first;
    longVector epsilon = f_ortho.second;
    //std::cout << "Initial Fock matrix:\n" << fock_0 << std::endl;

    // transform eigenvectors to original AO basis
    longMatrix C = S_minus_half * C_prime;
    hfu::cleanMatrix(C);
    // std::cout << "Initial MO Coeffs??\n" << C << std::endl;

    // Build density matrix using the occupied MOs:
    // Assuming number of occupied orbitals (for a closed-shell system, num_electrons/2)
    int num_occupied = 5;  // for H2O
    longMatrix D = hfu::initialize(7);  // hf_utils::initialize 
    // Build density matrix using occupied MOs
    for(int mu = 0; mu < 7; mu++) {
        for(int nu = 0; nu < 7; nu++) {
            for(int m = 0; m < num_occupied; m++) {
                D(mu, nu) += C(mu, m) * C(nu, m);  
            }
        }
    }
    // validated
    // std::cout << "Initial density matrix:\n" << D << std::endl; 



    // Step Six: Computing Initial SCF Energy
    double E_elec {0};
    for(int mu = 0; mu < 5; mu++) {
        for(int nu = 0; nu < 5; nu++) {
                E_elec += D(mu, nu) * (H(mu, nu) + fock_0(mu, nu));  
        }
    }
    double E_tot = E_elec + enuc; 
    // Result is -124.441 Hartrees - the guide's result was -125.842077437699 Hartrees. Maybe increased floating point precision?
    // std::cout << "E Total\n" << E_tot << std::endl;
    
    //Step Seven: Compute new Fock Matrix
    longMatrix fock_n = FockBuilder::build(H, D, dintegrals);

    //cleanMatrix(fock_n);
    //std::cout << "NewFock matrix:\n" << fock_n << std::endl;
    
    // Step Eight onwards: Use SCF class
    SCF scf_solver(H, S_minus_half, dintegrals, enuc);
    scf_solver.runSCF();

    if(scf_solver.hasConverged()) {
        std::cout << "Final Total Energy: " << scf_solver.getFinalEnergy() << std::endl;
    }

    return 0;
}
