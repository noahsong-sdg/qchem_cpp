#include <iostream>
#include <vector>
#include <fstream>
#include <cstdint>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include "hf_fns.h"

int main() 
{
    // Step One: Read Enuc
    std::string fileName = "data/enuc.dat";
    std::fstream inputFile(fileName);
    // Check if the file is open
    file_error_msg(inputFile);
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
    cleanMatrix<longMatrix>(H);
    //std::cout << "Hamiltonian matrix:\n" << H << std::endl;



    // Step Three: Two-Electron Integrals
    // TODO: Speedup https://github.com/CrawfordGroup/ProgrammingProjects/blob/master/Project%2303/hints/hint3-2.md
    std::string fileName_d = "data/eri.dat";
    longVector dintegrals = double_integral::readDoubleFromFile(fileName_d); 
    cleanVector(dintegrals);
    for(int i = 0; i < dintegrals.size(); i++) {
        if(std::isinf(dintegrals(i)) || std::isnan(dintegrals(i))) {
            dintegrals(i) = 0.0;
            }
        }
     std::cout << dintegrals << std::endl;


    // Step Four: Build Orthogonalization Matrix (S * L = L * Lambda)
    // L is the matrix of eigenvectors
    // Lambda is the diagonal matrix of corresponding eigvals
    // Declare solver
    auto s_ortho = solve_eigen(s);
    longMatrix L = s_ortho.first;
    longVector lambda = s_ortho.second;
    // Calculate λ^(-1/2)
    longVector lambda_minus_half = lambda.array().pow(-0.5);
    // Construct diagonal matrix from λ^(-1/2)
    longMatrix lambda_minus_half_matrix = lambda_minus_half.asDiagonal();
    // Calculate S^(-1/2) = L λ^(-1/2) L^T
    longMatrix S_minus_half = L * lambda_minus_half_matrix * L.transpose();
    cleanMatrix<longMatrix>(S_minus_half);
    //std::cout << "S^-1/2 matrix\n" << S_minus_half << std::endl;



    // Step Five: build initial density guess
    // Initial Fock matrix
    longMatrix fock_0 = S_minus_half.transpose() * H * S_minus_half;
    cleanMatrix<longMatrix>(fock_0);
    // diagonalize
    auto f_ortho = solve_eigen(fock_0);
    longMatrix C_prime = f_ortho.first;
    longVector epsilon = f_ortho.second;
    //std::cout << "Initial Fock matrix:\n" << fock_0 << std::endl;

    // transform eigenvectors to original AO basis
    longMatrix C = S_minus_half * C_prime;
    cleanMatrix<longMatrix>(C);
    // std::cout << "Initial MO Coeffs??\n" << C << std::endl;

    // Build density matrix using the occupied MOs:
    // Assuming number of occupied orbitals (for a closed-shell system, num_electrons/2)
    int num_occupied = 5;  // for H2O
    longMatrix D = initialize(7);  // Initialize 
    // Build density matrix using occupied MOs
    for(int mu = 0; mu < 7; mu++) {
        for(int nu = 0; nu < 7; nu++) {
            for(int m = 0; m < num_occupied; m++) {
                D(mu, nu) += C(mu, m) * C(nu, m);  
            }
        }
    }
    //std::cout << "Initial density matrix:\n" << D << std::endl; 



    // Step Six: Computing Initial SCF Energy
    double E_elec {0};
    for(int mu = 0; mu < 5; mu++) {
        for(int nu = 0; nu < 5; nu++) {
                E_elec += D(mu, nu) * (H(mu, nu) + fock_0(mu, nu));  
        }
    }
    double E_tot = E_elec + enuc; 
    // Result is -124.441 Hartrees - the guide's result was -125.842077437699 Hartrees. Maybe increased floating point precision?
    std::cout << "E Total\n" << E_tot << std::endl;
    


    //Step Seven: Compute new Fock Matrix
    longMatrix fock_n = ComputeNewFock::computeFock(H, D, dintegrals);
/*    longMatrix fock_n = initialize(7);
    for(uint16_t mu = 0; mu < 7; mu++) {
        for(uint16_t nu = 0; nu < 7; nu++) {
            fock_n(mu, nu) = H(mu, nu);
            long double value {0};
            for(uint16_t lmda = 0; lmda < 7; lmda++) {
                for(uint16_t sig = 0; sig < 7; sig++) {
                    value += D(lmda, sig) * 
                        (2.0*dintegrals(double_integral::compoundIndex(mu, nu, lmda, sig))
                        - dintegrals(double_integral::compoundIndex(mu, lmda, nu, sig)));
                        }
                    }
            fock_n(mu, nu) += value;
        }
    } */
    cleanMatrix(fock_n);
    std::cout << "NewFock matrix:\n" << fock_n << std::endl;



    // Step Eight: orthog F_new
    longMatrix fock_n_prime = S_minus_half.transpose() * fock_n * S_minus_half;
    auto f_prime_ortho = solve_eigen(fock_n_prime);
    longMatrix C_n_prime = f_prime_ortho.first;
    longMatrix C_n = S_minus_half * C_n_prime;
    longMatrix D_n = initialize(7);  // Initialize 
    // Build density matrix using occupied MOs
    for(int mu = 0; mu < 7; mu++) {
        for(int nu = 0; nu < 7; nu++) {
            for(int m = 0; m < num_occupied; m++) {
                D_n(mu, nu) += C_n(mu, m) * C_n(nu, m);  
            }
        }
    }
    // Step Nine: Compute new electronic and total energies:
    double E_elec_n {0};
    for(int mu = 0; mu < 5; mu++) {
        for(int nu = 0; nu < 5; nu++) {
                E_elec_n += D_n(mu, nu) * (H(mu, nu) + fock_n(mu, nu));  
        }
    }
    double E_tot = E_elec + enuc;

    // Step Ten: Test convergence
    // move these consts to header file
    const double density_threshold = 1e-12;    // δ1
    const double energy_threshold = 1e-12;     // δ2
    const int max_iterations = 50;
    bool converged = false;
    int iteration = 0;
    longMatrix D_old = D_n;
    double E_elec_old = E_elec_n;

    while(!converged && iteration < max_iterations) {
        // Build new Fock matrix
        longMatrix fock_n = ComputeNewFock::computeFock(H, D_old, dintegrals);
        
        // Transform and diagonalize Fock matrix
        longMatrix fock_n_prime = S_minus_half.transpose() * fock_n * S_minus_half;
        auto f_prime_ortho = solve_eigen(fock_n_prime);
        longMatrix C_n_prime = f_prime_ortho.first;
        longMatrix C_n = S_minus_half * C_n_prime;
        
        // Build new density matrix
        longMatrix D_new = initialize(7);
        for(int mu = 0; mu < 7; mu++) {
            for(int nu = 0; nu < 7; nu++) {
                for(int m = 0; m < num_occupied; m++) {
                    D_new(mu, nu) += C_n(mu, m) * C_n(nu, m);
                }
            }
        }
        
        // Calculate new electronic energy
        double E_elec_new = 0.0;
        for(int mu = 0; mu < 7; mu++) {
            for(int nu = 0; nu < 7; nu++) {
                E_elec_new += D_new(mu, nu) * (H(mu, nu) + fock_n(mu, nu));
            }
        }
        
        // Check both convergence criteria
        double rmsd = calculateRMSDensity(D_new, D_old);
        double delta_E = std::abs(E_elec_new - E_elec_old);
        
        // Print iteration info
        std::cout << "Iteration " << iteration << ":" << std::endl;
        std::cout << "  RMSD: " << rmsd << std::endl;
        std::cout << "  ΔE: " << delta_E << std::endl;
        std::cout << "  E_elec: " << E_elec_new << std::endl;
        std::cout << "  E_total: " << E_elec_new + enuc << std::endl;
        
        // Check if both convergence criteria are met
        if(rmsd < density_threshold && delta_E < energy_threshold) {
            converged = true;
        }
        
        // Update matrices and energy for next iteration
        D_old = D_new;
        E_elec_old = E_elec_new;
        iteration++;
}

if(converged) {
    std::cout << "\nSCF Converged in " << iteration << " iterations!" << std::endl;
    std::cout << "Final Electronic Energy: " << E_elec_old << std::endl;
    std::cout << "Final Total Energy: " << E_elec_old + enuc << std::endl;
} else {
    std::cout << "\nSCF Failed to converge in " << max_iterations << " iterations." << std::endl;
}
    // g++ main.cpp hf_fns.cpp -I. -o main

    return 0;
}
