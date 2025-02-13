#ifndef LJ_FORCE_H
#define LJ_FORCE_H

#include <armadillo>
#include "lj_energy.h" 

class LennardJonesForce : public LennardJones 
{
public:
    // Constructor
    LennardJonesForce(double eps, double sig);

    // Function to calculate the analytical forces
    arma::mat calculate_analytical_forces(const arma::mat & coords);

    // Function to calculate forward difference
    arma::mat calculate_forward_difference(const arma::mat & coords, double h);

    // Function to calculate central difference
    arma::mat calculate_central_difference(const arma::mat & coords, double h);

};

#endif
