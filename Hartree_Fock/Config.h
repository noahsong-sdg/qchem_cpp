
#ifndef CONFIG_H
#define CONFIG_H

struct HFConfig {
    static constexpr double DENSITY_THRESHOLD = 1e-12;
    static constexpr double ENERGY_THRESHOLD = 1e-12;
    static constexpr int MAX_ITERATIONS = 50;
    static constexpr double CLEAN_THRESHOLD = 1e-10;
    static constexpr int NUM_OCCUPIED = 5;  // Make this configurable
    static constexpr int BATCH_SIZE = 1000;
    static constexpr bool USE_SPARSE = false;
    static constexpr int NUM_THREADS = 4;
};

#endif
