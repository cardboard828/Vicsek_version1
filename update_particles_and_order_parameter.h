// update_particles_and_order_parameter.h
#ifndef UPDATE_PARTICLES_AND_ORDER_PARAMETER_H
#define UPDATE_PARTICLES_AND_ORDER_PARAMETER_H

#include <vector>
#include <array>
#include "Particle.h"
#include <random>

void update_particles_and_order_parameter(
    std::vector<Particle>& particles, 
    std::vector<std::array<double, 2>>& new_orientation_2,
    std::vector<double>& order_parameter_x,
    std::vector<double>& order_parameter_y,
    std::uniform_real_distribution<>& white_noise,
    std::mt19937& engine,
    int N, 
    double v0, 
    double dt, 
    int L,
    double& order_parameter
);

#endif // UPDATE_PARTICLES_AND_ORDER_PARAMETER_H
