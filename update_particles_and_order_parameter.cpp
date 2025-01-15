// update_particles_and_order_parameter.cpp
#include "update_particles_and_order_parameter.h"
#include <cmath>
#include <omp.h>

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
)
{
    double parameter_x = 0;
    double parameter_y = 0;
    
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        
        #pragma omp for
        for (int i = 0; i < N; i++)
        {
            // Update orientation
            particles[i].orientation = atan2(new_orientation_2[i][1], new_orientation_2[i][0]);
            particles[i].orientation += white_noise(engine);  // Add white noise to the orientation

            // Update position
            particles[i].x += v0 * cos(particles[i].orientation) * dt;
            particles[i].y += v0 * sin(particles[i].orientation) * dt;

            // Periodic boundary conditions
            if (particles[i].x < 0) {
                particles[i].x += L;
            } else if (particles[i].x >= L) {
                particles[i].x -= L;
            }

            if (particles[i].y < 0) {
                particles[i].y += L;
            } else if (particles[i].y >= L) {
                particles[i].y -= L;
            }

            // Update order parameter contributions
            order_parameter_x[thread_id] += cos(particles[i].orientation);
            order_parameter_y[thread_id] += sin(particles[i].orientation);
        }
    }

    // Combine results from all threads to calculate the overall order parameter
    for (int i = 0; i < omp_get_num_threads(); i++) {
        parameter_x += order_parameter_x[i];
        parameter_y += order_parameter_y[i];
    }
    
    order_parameter = sqrt(parameter_x * parameter_x + parameter_y * parameter_y) / N;
}
