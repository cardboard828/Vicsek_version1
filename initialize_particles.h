// initialize_particles.h

#ifndef INITIALIZE_PARTICLES_H
#define INITIALIZE_PARTICLES_H

#include <random>
#include "Particle.h"  // Particleクラスの定義をインクルード

void initialize_particles(Particle* particles, int N, double L,std::mt19937& engine);

#endif
