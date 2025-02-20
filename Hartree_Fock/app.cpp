#include "hartree_fock.h"
#include <iostream>
// build cmd: g++ -o app app.cpp Hartree_Fock.cpp hf_fns.cpp Logger.cpp post_hf.cpp -I. -std=c++17 -O3
int main(int argc, char* argv[]) {
    try {
        if (argc != 2) {
            std::cerr << "Usage: " << argv[0] << " <data_directory>\n";
            return 1;
        }

        std::cout << "Attempting to read data from: " << argv[1] << std::endl;

        HartreeFockSolver solver(argv[1]);
        double energy = solver.solve();
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
