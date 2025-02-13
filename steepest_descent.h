#ifndef STEEPEST_DESCENT_H
#define STEEPEST_DESCENT_H

#include <armadillo>
#include "lj_energy.h"

// Function for line search
double line_search(LennardJones & lj, arma::mat& coords, const arma::mat& forces, double step_size);

// Function for steepest descent optimization
void steepest_descent(double epsilon, double sigma, LennardJones & lj, arma::mat& coords, double tolerance = 1e-6, int max_iter = 100);

#endif // STEEPEST_DESCENT_H
