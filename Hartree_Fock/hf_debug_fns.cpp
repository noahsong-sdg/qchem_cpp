#include "hf_debug_fns.h"
#include <iostream>

namespace hf_debug {

void findZeroElements(const longVector& dintegrals) {
    std::cout << "Indices where dintegrals = 0:" << std::endl;
    for(int i = 0; i < dintegrals.size(); i++) {
        if(dintegrals(i) == 0) {
            std::cout << "Index " << i << std::endl;
        }
    }
}

void debug_integral_mapping(uint64_t i, uint64_t j, uint64_t k, uint64_t l, long double value) {
    std::cout << "\nDebugging integral (" << i << "," << j << "," << k << "," << l << ") = " << value << std::endl;
    
    // Show all permutations and their compound indices
    auto perms = double_integral::generate_permutations(i, j, k, l);
    std::cout << "Permutations and their compound indices:" << std::endl;
    for (const auto& p : perms) {
        uint64_t idx = double_integral::compoundIndex(p[0]-1, p[1]-1, p[2]-1, p[3]-1);
        std::cout << "(" << p[0] << "," << p[1] << "," << p[2] << "," << p[3] << ") -> index " << idx << std::endl;
    }
    
    // Show intermediate calculations for the original indices
    uint64_t ij = (i-1) * ((i-1) + 1)/2 + (j-1);
    uint64_t kl = (k-1) * ((k-1) + 1)/2 + (l-1);
    std::cout << "\nDetailed calculation for original indices:" << std::endl;
    std::cout << "ij = " << i-1 << " * (" << i-1 << " + 1)/2 + " << j-1 << " = " << ij << std::endl;
    std::cout << "kl = " << k-1 << " * (" << k-1 << " + 1)/2 + " << l-1 << " = " << kl << std::endl;
    std::cout << "final index = (" << ij << " * (" << ij << " + 1)/2 + " << kl << ") = " 
              << double_integral::compoundIndex(i-1, j-1, k-1, l-1) << std::endl;
}

void verify_integral(uint64_t i, uint64_t j, uint64_t k, uint64_t l, 
                    long double expected_value, const longVector& storage) {
    uint64_t idx = double_integral::compoundIndex(i-1, j-1, k-1, l-1);
    if (idx >= storage.size()) {
        std::cout << "Index out of bounds for (" << i << "," << j << "," << k << "," << l << ")" << std::endl;
        return;
    }
    
    long double stored_value = storage[idx];
    std::cout << "\nVerifying (" << i << "," << j << "," << k << "," << l << "):" << std::endl;
    std::cout << "Expected: " << expected_value << std::endl;
    std::cout << "Stored at index " << idx << ": " << stored_value << std::endl;
    std::cout << "Difference: " << std::abs(expected_value - stored_value) << std::endl;
}

} // namespace hf_debug
