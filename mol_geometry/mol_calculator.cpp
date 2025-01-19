#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#define M_PI 3.14159265358979323846

/* Functionality so far:
    - Reading info from .dat file, which includes atomic number, and cartesian coordinates of the atom
    - Calculation of bond lengths between atoms and storage into a distance matrix
    - Calculation of bond angles between atoms
    - Calculation of out of plane angles between subsets of four atoms
    - Calculation of Torsion angles
   Future additions:
    - Center of Mass
    - Moments of Inertia
    - Rotational Constants
*/

struct Atom {
    int atomicNumber;
    double x, y, z;
};

// calculates distances between atoms using euclidian norm
double calculate_distance(const Atom& a1, const Atom& a2)
{
    double dx {a1.x - a2.x};
    double dy {a1.y - a2.y};
    double dz {a1.z - a2.z};
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

// function for part two: calculating interatomic distances. creates a symmetric matrix with pairwise distances.
std::vector<std::vector<double>> calculate_distance_matrix(const std::vector<Atom>& coordinates) {
    size_t n_atoms = coordinates.size();
    // initializes n_atoms x n_atoms matrix: n_atoms rows, columns are vectors with n_atoms 0.0 entries
    std::vector<std::vector<double>> distances(n_atoms, std::vector<double>(n_atoms, 0.0));
    
    for (size_t i = 0; i < n_atoms; i++) {
        for (size_t j = 0; j < n_atoms; j++) {
            distances[i][j] = calculate_distance(coordinates[i], coordinates[j]);
            distances[j][i] = calculate_distance(coordinates[i], coordinates[j]);
        }
    }
    return distances;
}

// Calculate unit vector between two atoms
std::vector<double> calculate_unit_vector(const Atom& a1, const Atom& a2) {
    double dx = a2.x - a1.x;
    double dy = a2.y - a1.y;
    double dz = a2.z - a1.z;
    double R = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    return {dx/R, dy/R, dz/R};
}

// Calculate dot product of two vectors
double dot_product(const std::vector<double>& v1, const std::vector<double>& v2) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

std::vector<double> cross_product(const std::vector<double>& v1, const std::vector<double>& v2) {
    double xhat = v1[1]*v2[2] - v1[2]*v2[1];
    double yhat = v1[2]*v2[0]- v1[0]*v2[2];
    double zhat = v1[0]*v2[1] - v1[1]*v2[0];
    return {xhat, yhat, zhat};
}


// Calculate bond angle between three atoms (i-j-k where j is central atom)
double calculate_bond_angle(const Atom& i, const Atom& j, const Atom& k) {
    // Calculate unit vectors
    auto eji = calculate_unit_vector(j, i);
    auto ejk = calculate_unit_vector(j, k);
    
    // Calculate cos(angle) using dot product
    double cos_angle = dot_product(eji, ejk);
    
    // Convert to degrees
    return std::acos(cos_angle) * 180.0 / M_PI;
}


// Part Four; Out of Plane Angles
double out_of_plane(const Atom& i, const Atom& j, const Atom& k, const Atom& l) {
    auto ekj = calculate_unit_vector(k, j);
    auto ekl = calculate_unit_vector(k, l);
    auto eki = calculate_unit_vector(k, i);

    double numerator = dot_product((cross_product(ekj, ekl)), eki);
    double denominator = std::sin(calculate_bond_angle(j, k, l) * M_PI / 180.0);

    return std::asin(numerator / denominator) * 180.0 / M_PI;
}
// Part Five; Torsion/Dihedral Angles
double torsion(const Atom& i, const Atom& j, const Atom& k, const Atom& l) {
    auto eij = calculate_unit_vector(i, j);
    auto ejk = calculate_unit_vector(j, k);
    auto ekl = calculate_unit_vector(j, k);
    
    double bond_angle_one = std::sin(calculate_bond_angle(i, j, k) * M_PI / 180.0);
    double bond_angle_two = std::sin(calculate_bond_angle(j, k, l) * M_PI / 180.0);

    double numerator = dot_product(cross_product(eij, ejk), cross_product(ejk, ekl));

    return std::asin(numerator / (bond_angle_one*bond_angle_two)) * 180.0 / M_PI;
}


// Part Six: Center of Mass (need to obtain a header file consisting of the masses of common isotopes of each element)
/*
std::vector<double> com(std::vector<Atom> atoms, int numAtoms) {
    int xcom {};
    for (int i = 0; i < numAtoms; i++) {
        xcom += ()
    }

} */

int main() {
    // File name
    std::string fileName = "acetaldehyde.dat";
    // ifstream is a class that represents an input file stream
    std::ifstream inputFile(fileName);

    // Check if the file is open
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open the file!" << std::endl;
        return 1;
    }

    // Read the number of atoms
    int numAtoms{};
    inputFile >> numAtoms;

    // Store atoms using std::vector with type <Atom>; vector allows dynamic changing of size
    std::vector<Atom> atoms;

    // Read each atom's data
    for (int i = 0; i < numAtoms; ++i) {
        Atom atom;
        inputFile >> atom.atomicNumber >> atom.x >> atom.y >> atom.z;
        // push_back() is basically an append operation on std::vector
        atoms.push_back(atom);
    }
    // Close the file
    inputFile.close();

    // Output the data to verify
    std::cout << "Number of atoms: " << numAtoms << std::endl;
    // const denotes read only iteration
    // auto denotes the type of each atom will be based on elements in atoms
    // 
    for (const auto& atom : atoms) {
        std::cout << "Atomic Number: " << atom.atomicNumber
                  << ", Coordinates: (" << atom.x << ", " << atom.y << ", " << atom.z << ")"
                  << std::endl;
    }


    // Part Two: intermolecular distances
    auto distances = calculate_distance_matrix(atoms);

    std::cout << "Distance Matrix (Angstroms):\n";
    for (const auto& row : distances) {
        for (double dist : row)  {
            std::cout << dist << "\t";
        }
        std::cout << '\n';
    }

    //Part Three: Bond Distances
    // Calculate and print all possible bond angles
    std::cout << "\nBond Angles (degrees):\n";
    for (size_t i = 0; i < atoms.size(); i++) {
        for (size_t j = 0; j < atoms.size(); j++) {
            if (i == j) continue;
            for (size_t k = j + 1; k < atoms.size(); k++) {
                if (k == i) continue;
                double angle = calculate_bond_angle(atoms[i], atoms[j], atoms[k]);
                std::cout << "Angle " << i << "-" << j << "-" << k
                        << ": " << angle << " degrees\n";
            }
        }
    }

    // Part Four: Out of Plane Angles
    std::cout << "\nOut of Plane Angle Calculation:\n";
    std::vector<int> idxs(4);
        
    // Get indices from user with bounds checking
    do {
        std::cout << "Enter four indices (0-" << atoms.size()-1 << ") separated by spaces: ";
        for (int i = 0; i < 4; i++) {
            if (!(std::cin >> idxs[i]) || idxs[i] < 0 || idxs[i] >= atoms.size()) {
                std::cout << "Invalid input! Indices must be between 0 and " << atoms.size()-1 << "\n";
                std::cin.clear();
                std::cin.ignore(10000, '\n');
                continue;
            }
        }
        break;
    } while (true);

    // Calculate and display the out of plane angle
    double out_angle = out_of_plane(atoms[idxs[0]], atoms[idxs[1]], 
                            atoms[idxs[2]], atoms[idxs[3]]);
    std::cout << "Out of plane angle for atoms " 
            << idxs[0] << "-" << idxs[1] << "-" 
            << idxs[2] << "-" << idxs[3] << ": " 
            << out_angle << " degrees\n";

    return 0; 
}
