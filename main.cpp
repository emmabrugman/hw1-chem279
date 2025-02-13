#include <iostream>
#include <armadillo>
#include <fstream>
#include <sstream>
#include "lj_energy.h"
#include "lj_force.h"
#include "steepest_descent.h"

// Function to read coordinates
arma::mat read_coordinates_from_file(const std::string & file_path) 
{
    std::ifstream file(file_path);
    if (!file.is_open()) 
    {
        std::cerr << "Error: Unable to open file" << std::endl;
        exit(1);
    }

    size_t num_atoms;
    file >> num_atoms;
    // Initialize coordinate matrix
    arma::mat coords(num_atoms, 3); 
    for (size_t i = 0; i < num_atoms; i++) 
    {
        int atomic_number;
        double x, y, z;
        file >> atomic_number >> x >> y >> z;
        if (atomic_number != 79) 
        {  // Only processing gold atoms
            std::cerr << "Error: Only Au atoms are allowed." << std::endl;
            exit(1);
        }
        // Store the coordinates
        coords(i, 0) = x;
        coords(i, 1) = y;
        coords(i, 2) = z;
    }

    return coords;
}

int main() 
{
    double epsilon = 5.29;  // kcal/mol for Au
    double sigma = 2.951;   // angstroms for Au

    // Initialize object
    LennardJones lj(epsilon, sigma);

    // 1.1: THE LENNARD-JONES ENERGY OF A SET OF ATOMS

    // File processing
    std::string input_file_1_1 = "sample_input/Energy/3.txt"; 
    arma::mat coords_1_1 = read_coordinates_from_file(input_file_1_1);

    // Echo the input
    std::cout << "Coordinates- x, y, z:\n";
    for (size_t i = 0; i < coords_1_1.n_rows; i++) 
    {
        std::cout << "Coordinates for atom 79: " << coords_1_1(i, 0) << " " << coords_1_1(i, 1) << " " << coords_1_1(i, 2) << std::endl;
    }

    // Calculate the total energy of the cluster
    double energy = lj.calculate_total_energy(coords_1_1);
    std::cout << "Total energy of the cluster: " << energy << " kcal/mol\n";

    // 1.2: FINITE DIFFERENCE VERSUS ANALYTICAL FORCES

    // File processing
    std::string input_file_1_2 = "sample_input/Force/1.txt";
    arma::mat coords_1_2 = read_coordinates_from_file(input_file_1_2);

    // Initialize force object
    LennardJonesForce lj_force(epsilon, sigma);

    // Step sizes
    std::vector<double> step_sizes = {0.1, 0.01, 0.001, 0.0001};

    // Calculate the analytical forces
    arma::mat analytical_forces = lj_force.calculate_analytical_forces(coords_1_2);

    std::cout << "h Log_Error_Forward Log_Error_Central\n";

    for (double h : step_sizes) {
        // Calculate forward and central differences
        arma::mat forward_diff_forces = lj_force.calculate_forward_difference(coords_1_2, h);
        arma::mat central_diff_forces = lj_force.calculate_central_difference(coords_1_2, h);

        // Compute errors
        arma::mat forward_error = arma::abs(analytical_forces - forward_diff_forces);
        arma::mat central_error = arma::abs(analytical_forces - central_diff_forces);

        double mean_forward_error = arma::mean(arma::mean(forward_error));
        double mean_central_error = arma::mean(arma::mean(central_error));

        double epsilon = 1e-10; // to avoid inf values
        double log_forward_error = std::log(mean_forward_error + epsilon);
        double log_central_error = std::log(mean_central_error + epsilon);

        std::cout << h << " " << log_forward_error << " " << log_central_error << "\n";
    }

    // 1.3: OPTIMIZING THE GEOMETRY OF A SMALL GOLD CLUSTER

    // File processing
    std::string input_file_1_3 = "sample_input/SD_with_line_search/1.txt"; 
    arma::mat coords_1_3 = read_coordinates_from_file(input_file_1_3);
    
    // Perform steepest descent optimization with line search
    steepest_descent(epsilon, sigma, lj, coords_1_3, 1e-6, 100);

    // Print optimized coordinates
    std::cout << "\nOptimized coordinates for cluster:\n";
    coords_1_3.print();

    return 0;
}
