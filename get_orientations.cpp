// get_orientations.cpp
#include "get_orientations.h"
#include "getNeighboringCell.h"
#include "dist.h"
#include <cmath>
#include <omp.h>

void get_orientations(
    std::vector<Particle>& particles, 
    std::vector<std::vector<std::vector<Particle>>>& cells,
    std::vector<std::array<double, 2>>& new_orientation_2,
    int L
)
{
    // loop over all particles
    #pragma omp parallel for
    for (size_t i = 0; i < particles.size(); i++)
    {
        int cell_x = floor(particles[i].x);
        int cell_y = floor(particles[i].y);
        
        // Get neighboring cells
        std::vector<std::array<int, 2>> neighbors = getNeighboringCell(cell_x, cell_y, L);

        // Loop over all neighboring cells
        for (size_t j = 0; j < neighbors.size(); j++)
        {
            int neighbor_x = neighbors[j][0];
            int neighbor_y = neighbors[j][1];

            // Loop over all particles in the neighboring cell
            for (size_t k = 0; k < cells[neighbor_x][neighbor_y].size(); k++)
            {
                // Calculate distance between particles
                double d = dist(particles[i], cells[neighbor_x][neighbor_y][k]);
                if (d < 1)
                {
                    // Add angle to the sum
                    new_orientation_2[i][0] += cos(cells[neighbor_x][neighbor_y][k].orientation);
                    new_orientation_2[i][1] += sin(cells[neighbor_x][neighbor_y][k].orientation);
                }
            }
        }
    }
}
