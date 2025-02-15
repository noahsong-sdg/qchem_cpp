#ifndef HF_DEBUG_FNS_H
#define HF_DEBUG_FNS_H

#include "./Eigen/Dense"
#include "hf_fns.h"

namespace hf_debug {
    void findZeroElements(const longVector& dintegrals);
    void debug_integral_mapping(uint64_t i, uint64_t j, uint64_t k, uint64_t l, long double value);
    void verify_integral(uint64_t i, uint64_t j, uint64_t k, uint64_t l, 
                        long double expected_value, const longVector& storage);
}

#endif // HF_DEBUG_FNS_H
