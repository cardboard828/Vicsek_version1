// initialize_particles.cpp

#include "initialize_particles.h"
#include <cmath>
// #include <omp.h>

void initialize_particles(Particle* particles, int N, double L, std::mt19937& engine)
{
    // 乱数生成器の定義
    std::uniform_real_distribution<> init_pos(0, L);
    std::uniform_real_distribution<> dist_angle(0, 2 * M_PI);

    // #pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        particles[i].x = init_pos(engine);
        particles[i].y = init_pos(engine);
        particles[i].orientation = dist_angle(engine);
    }
}
