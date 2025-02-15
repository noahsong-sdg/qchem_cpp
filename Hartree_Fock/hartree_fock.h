#ifndef HARTREE_FOCK_H
#define HARTREE_FOCK_H

#include "hf_fns.h"

class HartreeFockSolver {
public:
    HartreeFockSolver(const std::string& data_directory);
    double solve();
    bool hasConverged() const { return converged; }
    double getEnergy() const { return final_energy; }

private:
    std::string data_path;
    bool converged;
    double final_energy;
    
    // Core matrices
    longMatrix H;
    longMatrix S;
    longMatrix S_minus_half;  // Added this line
    longVector dintegrals;
    double enuc;
    
    // Helper methods
    void loadIntegrals();
    void buildInitialGuess();
    void runSCF();
};

#endif // HARTREE_FOCK_H
