#ifndef HARTREE_FOCK_H
#define HARTREE_FOCK_H

#include "hf_fns.h"
#include <memory>

class HartreeFockSolver {
public:
    HartreeFockSolver(const std::string& data_directory);
    ~HartreeFockSolver();
    
    // Main interface
    double solve();
    bool hasConverged() const { return converged; }
    double getEnergy() const { return total_energy; }
    
    // Additional getters for post-HF calculations
    const longMatrix& getMOCoefficients() const { return C; }
    const longVector& getMOEnergies() const { return epsilon; }
    const longVector& getIntegrals() const { return dintegrals; }
    
    // Add MP2 calculation method
    double calculateMP2();

private:
    // Forward declarations of implementation classes
    class IntegralLoader;
    class SCFIterator;
    
    // Private implementation classes (defined in cpp file)
    std::unique_ptr<IntegralLoader> loader;
    std::unique_ptr<SCFIterator> scf_iterator;
    
    // Core state
    bool converged;
    double total_energy;
    longMatrix C;          // MO coefficients
    longVector epsilon;    // Orbital energies
    
    // Core matrices (shared between implementations)
    longMatrix H;          // Core Hamiltonian
    longMatrix S;          // Overlap
    longMatrix D;          // Density
    longVector dintegrals; // Two-electron integrals
    double enuc;          // Nuclear repulsion
    
    // Configuration
    const std::string data_path;
    
    // Private methods
    void initializeSystem();
    void runSCF();
    
    // Add this to store MP2 correction
    double mp2_energy;
};

#endif
