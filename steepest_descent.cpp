#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <sstream>
#include "lj_energy.h"
#include "lj_force.h"

// Function for line search
double line_search(LennardJones& lj, arma::mat& coords, const arma::mat& forces, double step_size) 
{

    // Define golden ratio
    const double golden_ratio = (1.0 + std::sqrt(5.0)) / 2.0;

    // Parameters
    double alpha = 0.0;
    double beta = step_size;

    double f_alpha = lj.calculate_total_energy(coords - alpha * forces);
    double f_beta = lj.calculate_total_energy(coords - beta * forces);

    // Loop through iterations
    for (int iter = 0; iter < 100; iter++) 
    {
        // Compute gamma
        double gamma = beta - (beta - alpha) / golden_ratio;
        // Calculate total energy
        double f_gamma = lj.calculate_total_energy(coords - gamma * forces);
        // Update interval based energy comparison
        if (f_gamma < f_alpha && f_gamma < f_beta) 
        {
            beta = gamma;
            f_beta = f_gamma;
        } else 
        {
            alpha = gamma;
            f_alpha = f_gamma;
        }

        // Break out of the loop
        if (std::abs(beta - alpha) < 1e-6) break;
    }

    return (alpha + beta) / 2.0;
}

// Function for steepest descent optimization with line search
void steepest_descent(double epsilon, double sigma, LennardJones& lj, arma::mat& coords, double tolerance = 1e-6, int max_iter = 100) 
{
    // LennardJonesForce object
    LennardJonesForce lj_force(epsilon, sigma);  
    // calculate energy
    double previous_energy = lj.calculate_total_energy(coords);
    double energy_change = 1.0;

    for (int iter = 0; iter < max_iter && energy_change > tolerance; iter++) 
    {
        // Calculate analytical forces
        arma::mat forces = lj_force.calculate_analytical_forces(coords);
        // Perform line search
        double step_size = line_search(lj, coords, forces, 0.1);
        // Update coordinates
        coords -= step_size * forces;
        // Calculate new energy
        double new_energy = lj.calculate_total_energy(coords);
        // Calculate the change in energy
        energy_change = std::abs(new_energy - previous_energy);
        previous_energy = new_energy;
        std::cout << "Iteration " << iter + 1 << ": Energy = " << new_energy << ", Energy change = " << energy_change << std::endl;
    }

    std::cout << "Final energy: " << previous_energy << std::endl;
}






