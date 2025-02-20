#ifndef POST_HF_H
#define POST_HF_H

#include "hf_fns.h"

namespace PostHF {
    class MP2 {
    public:
        MP2(const longMatrix& mo_coeffs, 
            const longVector& mo_energies,
            const longVector& eri,
            int num_occupied,
            int num_virtuals);
        
        double compute();
        double getEnergy() const { return mp2_energy; }
        
    private:
        const longMatrix& C;     // MO coefficients
        const longVector& eps;   // Orbital energies
        const longVector& eri;   // Two-electron integrals
        const int nocc;         // Number of occupied orbitals
        const int nvirt;        // Number of virtual orbitals
        double mp2_energy;
        
        // Helper methods
        double transform_eri(int i, int j, int a, int b) const;
        void transform_batch(std::vector<double>& mp2_ints, 
                           const std::vector<size_t>& batch_indices) const;
    };
}

#endif
