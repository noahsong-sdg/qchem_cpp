#include "hartree_fock.h"
#include "post_hf.h"
#include <fstream>
#include <iostream>

// Private implementation classes
class HartreeFockSolver::IntegralLoader {
public:
    void loadIntegrals(const std::string& data_path, double& enuc, 
                      longMatrix& H, longMatrix& S, longVector& dintegrals) {
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
};

class HartreeFockSolver::SCFIterator {
public:
    void iterate(const longMatrix& H, const longMatrix& S, longVector& dintegrals,
                double enuc, longMatrix& C, longVector& epsilon, bool& converged,
                double& total_energy) {
        // Get eigenvectors and eigenvalues of overlap matrix
        auto s_ortho = hf_utils::solve_eigen(S);
        longMatrix L = s_ortho.first;
        longVector lambda = s_ortho.second;
        
        // Calculate S^(-1/2)
        longVector lambda_minus_half = lambda.array().pow(-0.5);
        longMatrix S_minus_half = L * lambda_minus_half.asDiagonal() * L.transpose();
        hf_utils::cleanMatrix(S_minus_half);
        
        // Initialize density matrix
        longMatrix D = longMatrix::Zero(H.rows(), H.rows());
        double E_old = 0.0;
        int iteration = 0;
        
        while(!converged && iteration < HFConfig::MAX_ITERATIONS) {
            // Build Fock matrix
            longMatrix fock = FockBuilder::build(H, D, dintegrals);
            
            // Transform Fock matrix
            longMatrix fock_prime = S_minus_half.transpose() * fock * S_minus_half;
            
            // Solve eigenvalue problem
            auto eigen_solution = hf_utils::solve_eigen(fock_prime);
            C = S_minus_half * eigen_solution.first;
            epsilon = eigen_solution.second;
            
            // Build new density matrix
            longMatrix D_new = longMatrix::Zero(H.rows(), H.rows());
            for(int mu = 0; mu < H.rows(); mu++) {
                for(int nu = 0; nu < H.rows(); nu++) {
                    for(int m = 0; m < HFConfig::NUM_OCCUPIED; m++) {
                        D_new(mu, nu) += C(mu, m) * C(nu, m);
                    }
                }
            }
            // Calculate energy and check convergence
            double E_new = 0.0;
            for(int mu = 0; mu < H.rows(); mu++) {
                for(int nu = 0; nu < H.rows(); nu++) {
                    E_new += D_new(mu, nu) * (H(mu, nu) + fock(mu, nu));
                }
            }
            
            double rmsd = hf_utils::calculateRMSDensity(D_new, D);
            double delta_E = std::abs(E_new - E_old);
            
            // Print iteration info
            Logger::log("Iteration " + std::to_string(iteration) + 
                       ": E = " + std::to_string(E_new + enuc) +
                       ", RMSD = " + std::to_string(rmsd));
            
            if(rmsd < HFConfig::DENSITY_THRESHOLD && delta_E < HFConfig::ENERGY_THRESHOLD) {
                converged = true;
                total_energy = E_new + enuc;
            }
            
            D = D_new;
            E_old = E_new;
            iteration++;
        }
    }
};

// Main class implementation
HartreeFockSolver::HartreeFockSolver(const std::string& data_directory)
    : data_path(data_directory), converged(false), total_energy(0.0),
      loader(std::make_unique<IntegralLoader>()),
      scf_iterator(std::make_unique<SCFIterator>()) {
    initializeSystem();
}

HartreeFockSolver::~HartreeFockSolver() = default;

void HartreeFockSolver::initializeSystem() {
    loader->loadIntegrals(data_path, enuc, H, S, dintegrals);
}

double HartreeFockSolver::solve() {
    Logger::init(data_path + "/hartree_fock.log");
    Logger::log("Starting Hartree-Fock calculation...");
    
    scf_iterator->iterate(H, S, dintegrals, enuc, C, epsilon, converged, total_energy);
    
    if(converged) {
        Logger::log("SCF Converged! Final Energy: " + std::to_string(total_energy));
        // Calculate MP2 correction
        mp2_energy = calculateMP2();
        total_energy += mp2_energy;
        Logger::log("MP2 Correction: " + std::to_string(mp2_energy));
        Logger::log("Total Energy with MP2: " + std::to_string(total_energy));
    } else {
        Logger::log("SCF Failed to converge.");
    }
    
    return total_energy;
}

double HartreeFockSolver::calculateMP2() {
    PostHF::MP2 mp2(getMOCoefficients(), 
                    getMOEnergies(),
                    getIntegrals(),
                    HFConfig::NUM_OCCUPIED,
                    C.rows() - HFConfig::NUM_OCCUPIED);
    return mp2.compute();
}


