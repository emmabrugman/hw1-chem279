#include "lj_force.h"
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>

// Constructor for LennardJonesForce class
LennardJonesForce::LennardJonesForce(double eps, double sig) : LennardJones(eps, sig) {}

// Calculate analytical forces
arma::mat LennardJonesForce::calculate_analytical_forces(const arma::mat & coords) 
{
    size_t n_atoms = coords.n_rows;
    // Initialize a matrix to store forces 
    arma::mat forces(n_atoms, 3, arma::fill::zeros);

    // Loop through atom pairs
    for (size_t i = 0; i < n_atoms; i++) 
    {
        for (size_t j = 0; j < n_atoms; j++) 
        {
             if (i != j) 
             {
                // Calculate the vector between atoms i and j
                arma::rowvec r_ij = coords.row(i) - coords.row(j);  // Direction vector between atoms i and j
                // Calculate the distance between atoms i and j
                double distance = arma::norm(r_ij);
                // Calculate the force (Lennard-Jones formula)
                arma::rowvec force = 24 * epsilon * (2 * pow(sigma / distance, 12) - pow(sigma / distance, 6)) / pow(distance, 2) * r_ij;
                forces.row(i) += force;  // Add the force to atom i's force
            }
        }
    }
    return forces;
}

// Calculate forward difference force
arma::mat LennardJonesForce::calculate_forward_difference(const arma::mat & coords, double h) 
{
    size_t n_atoms = coords.n_rows;
    // Initialize a matrix to store forces
    arma::mat forces(n_atoms, 3, arma::fill::zeros);

    // Loop through atoms and coordinates(x,y,z)
    for (size_t i = 0; i < n_atoms; i++) 
    {
        for (size_t j = 0; j < 3; j++) 
        {
            // Create a copy of the coordinates
            arma::mat new_coords = coords;
            new_coords(i, j) += h;
            // Calculate energy for coordinates + h
            double energy_plus = calculate_total_energy(new_coords);
            // Calculate energy for original coordinates
            double energy_original = calculate_total_energy(coords);
            // Calculate force using forward difference approximation
            forces(i, j) = -(energy_plus - energy_original) / h;
        }
    }
    return forces;
}

// Calculate central difference force
arma::mat LennardJonesForce::calculate_central_difference(const arma::mat & coords, double h) 
{
    size_t n_atoms = coords.n_rows;
    // Initialize a matrix to store forces
    arma::mat forces(n_atoms, 3, arma::fill::zeros);

    // Loop through atoms and coordinates(x,y,z)
    for (size_t i = 0; i < n_atoms; i++) 
    {
        for (size_t j = 0; j < n_atoms; j++) 
        {
            if (i != j) 
            {
                // Calculate the vector between atoms i and j
                arma::rowvec r_ij = coords.row(i) - coords.row(j);
                // Distance between atoms i and j
                double distance = arma::norm(r_ij);

                // Calculate repulsive and attractive terms
                double repulsive_term = pow(sigma / distance, 12);
                double attractive_term = pow(sigma / distance, 6);

                // Calculate the total force between atoms i and j
                arma::rowvec force = 24 * epsilon * (2 * repulsive_term - attractive_term) / pow(distance, 2) * r_ij;

                forces.row(i) += force;
            } 
        }
    }
    return forces;
}
